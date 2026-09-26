function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
    newChannelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
    bsUnavailReasonVec, satUnavailReasonVec, policyKpis] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, prevChannelState)
% Για κάθε χρήστη: επιλέγει τον καλύτερο κόμβο (BS ή δορυφόρο) βάσει SINR και
% υπολογίζει χωρητικότητα/ενέργεια μετά την κατανομή εύρους ζώνης.
% Επιστρέφει και per-candidate διαγνωστικά (καλύτερο BS + δορυφόρος, ανεξάρτητα
% από την τελική επιλογή) για την εκπαίδευση του ML μοντέλου του Part 2.
%
% 23η έξοδος policyKpis: η ίδια πραγματοποίηση καναλιού αποτιμάται και με τις
% τρεις πολιτικές επιλογής κόμβου (η υπό εξέταση, αποκλειστικά επίγεια,
% αποκλειστικά δορυφορική), με την κατανομή πόρων να υπολογίζεται εκ νέου σε
% καθεμία. Τα δύο υποψήφια SINR έχουν ήδη υπολογιστεί για κάθε χρήστη ανεξάρτητα
% από το ποιος κερδίζει, οπότε οι τρεις πολιτικές διαφέρουν μόνο στην επιλογή
% και η σύγκριση είναι κατά ζεύγη πάνω στο ίδιο κανάλι.
%
% prevChannelState (προαιρετικό, 7ο όρισμα): αν δοθεί, LOS/NLOS και shadow
% fading κάθε ζεύξης συσχετίζονται χωρικά με την προηγούμενη κλήση αντί για
% i.i.d. δειγματοληψία. Το περνάει το runSimulation.m μέσα σε κάθε διέλευση
% (διαδοχικά βήματα 1 s με κινούμενους χρήστες) και το μηδενίζει στην αρχή
% κάθε νέας διέλευσης, ώστε οι διελεύσεις να είναι ανεξάρτητες μεταξύ τους.
% Τα scripts επαλήθευσης το παραλείπουν.
if nargin < 7
    prevChannelState = [];
end

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

%% ------------------ Bandwidth ------------------
BW_bs = simParameters.Carrier.NSizeGrid * 12 * ...
        simParameters.Carrier.SubcarrierSpacing * 1e3;   % Hz

%% ------------------ Noise power ------------------
kBoltz = physconst('Boltzmann');
NF = 10^(simParameters.RxNoiseFigure/10);
Teq = simParameters.RxAntTemperature + 290*(NF-1);
% kTB σε dBW
noisePowerBS_dBW  = 10*log10(kBoltz * Teq * BW_bs);
noisePowerSAT_dBW = 10*log10(kBoltz * Teq * satParameters.Bandwidth);

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SINR (κατάσταση outage) ------------------
% Κάτω από αυτό -> ο χρήστης θεωρείται outage αντί να ανατεθεί στον
% "λιγότερο κακό" κόμβο. = Shannon-ισοδύναμο SINR του MCS 0 (TS 38.214 Πίν. 5.1.3.1-1).
minSpectralEfficiency = 0.2344;                        % bits/s/Hz (MCS 0)
minUsableSinrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Κατώφλι σε επίπεδο υπηρεσίας ------------------
% 5ο εκατοστημόριο φασματικής απόδοσης χρήστη (απαίτηση ITU-R M.2410, όπως
% παρατίθεται στο TR 37.910). Η τιμή εξαρτάται από το περιβάλλον:
% 0.3 bit/s/Hz για Indoor Hotspot-eMBB DL (Πίν. 5.4.1.1.1-1) και
% 0.225 bit/s/Hz για Dense Urban-eMBB DL (Πίν. 5.4.1.2.1-1). Εφαρμόζεται η
% αυστηρότερη τιμή. Ορίζεται επί του ΣΥΝΟΛΙΚΟΥ εύρους καναλιού, άρα
% SE_ζεύξης >= targetNormalizedSe * L.
targetNormalizedSe = 0.3;   % bit/s/Hz

%% ------------------ Ισχύς εκπομπής ανά αλυσίδα πομποδέκτη ------------------
% Το μοντέλο EARTH (Auer et al. 2011, εξ. 1) ορίζει P_out ΑΝΑ αλυσίδα, με
% ανώτατο όριο P_max = 20 W για μακροκυψελικό σταθμό (Πίν. 2).
pOutTotalW    = 10^((simParameters.TxPower - 30)/10);
pOutPerChainW = pOutTotalW / simParameters.Power.NumTrx;
earthPmaxW    = 20;
if pOutPerChainW > earthPmaxW
    warning('simulateScenario:EarthOutOfRange', ...
        ['Ισχύς ανά αλυσίδα %.1f W > P_max = %.0f W του μοντέλου EARTH ' ...
         '(Auer et al. 2011, Πίν. 2). Αύξησε το simParameters.Power.NumTrx.'], ...
        pOutPerChainW, earthPmaxW);
end

%% ------------------ Διαλείψεις δορυφορικής ζεύξης ------------------
% Shadowed Rician (Abdi et al. 2003). Η σκίαση περιέχεται ήδη στο μοντέλο
% (τυχαίο πλάτος LOS κατά Nakagami-m), οπότε δεν προστίθεται χωριστός
% λογαριθμοκανονικός όρος. Παράμετροι (b0,m,Ω) από την ανύψωση, εξ. (19).
satFadeElevRangeDeg = [20 80];   % πεδίο ισχύος της προσαρμογής της εξ. (19)
terrKdBLos          = 9;         % Rician K επίγειο LOS, TR 38.901 Πίν. 7.5-6

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestBsNodeVec       = strings(numUsers,1);   % ποιος BS κέρδισε, για τη σύγκριση πολιτικών
bestDistanceVec     = nan(numUsers,1);
bestPathLossVec     = nan(numUsers,1);
bestSinrDbVec        = nan(numUsers,1);
capacityMbpsVec     = nan(numUsers,1);
nodePowerWattsVec   = nan(numUsers,1);
energyPerBitUJVec   = nan(numUsers,1);
bestElevationDegVec = nan(numUsers,1);
bestBsSinrDbVec      = nan(numUsers,1);
bestBsDistanceVec   = nan(numUsers,1);
bestBsPathLossVec   = nan(numUsers,1);
serviceStateVec     = strings(numUsers,1);
throughputMbpsVec   = zeros(numUsers,1);
bsUnavailReasonVec  = strings(numUsers,1);
satUnavailReasonVec = strings(numUsers,1);

% Διαγνωστικοί πίνακες
groundDistanceMat = nan(numUsers,numBs);
range3DMat        = nan(numUsers,numBs);
pathLossMat       = nan(numUsers,numBs);
snrDbMat          = nan(numUsers,numBs);   % χωρίς παρεμβολή, διαγνωστικό
sinrDbMat         = nan(numUsers,numBs);
pLosMat           = nan(numUsers,numBs);
losMat            = false(numUsers,numBs);
losLatentMat      = zeros(numUsers,numBs);
sfMat             = nan(numUsers,numBs);
satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSinrDbVec       = nan(numUsers,1);

hasPrevState = ~isempty(prevChannelState) && ...
    isfield(prevChannelState, 'LosLatent') && ...
    isequal(size(prevChannelState.IsLOS), [numUsers, numBs]);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SINR) ------------------
for u = 1:numUsers
    % Αρχικοποιήσεις
    userBestSinr      = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

    % Απόσταση μετακίνησης του χρήστη από την προηγούμενη κλήση· καθορίζει
    % πόσο συσχετίζεται το shadow fading με την προηγούμενη τιμή (βλ. correlatedLosState).
    if hasPrevState
        userMoveDistance = distance(prevChannelState.UserGeo(u,1), prevChannelState.UserGeo(u,2), ...
                                     user_geo(u,1), user_geo(u,2), wgs84);
    else
        userMoveDistance = Inf; % καμία προηγούμενη κατάσταση -> ανεξάρτητο δείγμα, όπως πριν
    end

    %% ===== Terrestrial BS candidates =====
    for b = 1:numBs
        % Οριζόντια απόσταση (γεωδαιτική) και φυσικά ύψη κεραιών. Η γεωμετρία
        % για τη nrPathLoss κατασκευάζεται απευθείας από αυτά: η κατακόρυφη
        % συντεταγμένη ENU αποκλίνει από το ύψος κεραίας λόγω καμπυλότητας.
        groundDistance = distance(bs_geo(b,1), bs_geo(b,2), ...
                                  user_geo(u,1), user_geo(u,2), wgs84);
        hBs = bs_geo(b,3);
        hUt = user_geo(u,3);
        d3d = hypot(groundDistance, hBs - hUt);

        txPosition = [0; 0; hBs];
        rxPosition = [groundDistance; 0; hUt];

        groundDistanceMat(u,b) = groundDistance;
        range3DMat(u,b)        = d3d;

        % Πεδίο ισχύος UMa/UMi (TR 38.901 Πίν. 7.4.1-1). Εκτός ορίων η ζεύξη
        % δεν υπολογίζεται: το SINR μένει NaN, ώστε να ξεχωρίζει από ζεύξη που
        % υπολογίστηκε και βρέθηκε ανεπαρκής.
        if ~isValidTerrestrialLink(groundDistance, hUt, simParameters.CarrierFrequency)
            sfMat(u,b) = 0;
            continue;
        end

        % LOS ανά ζεύξη βάσει πιθανότητας απόστασης (TR 38.901 §7.4.2).
        pLos = losProbability38901(groundDistance, user_geo(u,3), simParameters.PathLoss.Scenario);
        if hasPrevState
            prevIsLos  = prevChannelState.IsLOS(u,b);
            prevSF     = prevChannelState.ShadowFading_dB(u,b);
            prevLatent = prevChannelState.LosLatent(u,b);
        else
            prevIsLos  = false;
            prevSF     = 0;
            prevLatent = 0;
        end
        [isLos, losLatent] = spatiallyConsistentLos(pLos, userMoveDistance, ...
            hasPrevState, prevLatent);
        rho = shadowFadingCorrelation(userMoveDistance, ...
            simParameters.PathLoss.Scenario, isLos);
        useCorrelatedSF = hasPrevState && (isLos == prevIsLos);
        pLosMat(u,b)     = pLos;
        losMat(u,b)      = isLos;
        losLatentMat(u,b) = losLatent;

        [pathLoss, sigmaSF] = nrPathLoss(simParameters.PathLoss, ...
                              simParameters.CarrierFrequency, ...
                              isLos, ...
                              txPosition, rxPosition);

        % Shadow fading: log-normal δείγμα με τυπική απόκλιση sigmaSF (TR 38.901 §7.4.1),
        % AR(1)-συσχετισμένο με το προηγούμενο όταν useCorrelatedSF, αλλιώς ανεξάρτητο.
        if useCorrelatedSF
            sfSample = rho*prevSF + sqrt(1 - rho^2) * sigmaSF * randn();
        else
            sfSample = sigmaSF * randn();
        end
        sfMat(u,b) = sfSample;
        pathLoss = pathLoss + sfSample;

        % Small-scale fading: Rician (LOS, K=terrKdBLos) ή Rayleigh (NLOS).
        % Realization i.i.d. ανά κλήση (coherence time ~ ms << βήμα).
        pathLoss = pathLoss - smallScaleFadingDb(isLos, terrKdBLos);
        pathLossMat(u,b) = pathLoss;

        % Ο λόγος σήματος προς θόρυβο ΜΟΝΟ, χωρίς παρεμβολή, κρατείται ως
        % διαγνωστικό: είναι το άνω φράγμα του SINR και δείχνει πόσο κοστίζει
        % η παρεμβολή σε κάθε ζεύξη.
        snrDbMat(u,b) = (simParameters.EIRP - 30) - pathLoss - noisePowerBS_dBW;
    end

    %% ===== SINR ανά επίγειο υποψήφιο =====
    % Όλοι οι σταθμοί εκπέμπουν ταυτόχρονα στο ίδιο φάσμα: επαναχρησιμοποίηση
    % συχνοτήτων 1 και μοντέλο πλήρους απασχόλησης, όπως ορίζει ο ΠΙΝΑΚΑΣ 5 β)
    % του ITU-R M.2412-0 (§8.4, "Inter-site interference modeling: Explicitly
    % modelled", "Traffic model: Full buffer") και όπως προϋποθέτουν οι μετρικές
    % βαθμονόμησης του TR 38.901 §7.8 (Πίν. 7.8-1 "Geometry", Πίν. 7.8-2
    % "Wideband SIR", Πίν. 7.8-3 "Wideband SINR").
    %
    % Ορισμός: TR 38.821 Πίν. 6.1.1.2-1, ΣΗΜΕΙΩΣΗ:
    %   Geometry SINR = -10*log10(I/C + N/C)  <=>  SINR = C / (I + N)
    % με C, I, N μετρημένα στο ίδιο εύρος ζώνης. Επειδή σήμα, παρεμβολή και
    % θόρυβος κλιμακώνονται όλα με το εκχωρημένο εύρος, ο λόγος δεν εξαρτάται
    % από την κατανομή πόρων: υπολογίζεται μία φορά στο εύρος του κόμβου.
    %
    % Σταθμός εκτός του πεδίου ισχύος του μοντέλου δεν προσμετράται ως
    % παρεμβολέας, γιατί οι απώλειές του δεν είναι υπολογίσιμες με το UMa/UMi
    % χωρίς να παραβιαστεί ο έλεγχος εγκυρότητας. Η παρεμβολή είναι επομένως
    % κάτω φράγμα, όπως και λόγω του μικρού πλήθους σταθμών.
    noiseLinW   = 10^(noisePowerBS_dBW/10);
    rxPowLinW   = 10.^(((simParameters.EIRP - 30) - pathLossMat(u,:))/10);
    rxPowLinW(~isfinite(rxPowLinW)) = 0;
    totalRxLinW = sum(rxPowLinW);

    for b = 1:numBs
        if ~isfinite(pathLossMat(u,b))
            continue;
        end
        interfLinW     = totalRxLinW - rxPowLinW(b);
        sinr_db        = 10*log10(rxPowLinW(b) / (noiseLinW + interfLinW));
        sinrDbMat(u,b) = sinr_db;

        % Η επιλογή είναι μονότονη ως προς τη λαμβανόμενη ισχύ: με κοινό
        % άθροισμα ισχύων, το SINR αυξάνει με το C, οπότε ο καλύτερος σταθμός
        % είναι ο ίδιος με ή χωρίς παρεμβολή. Αλλάζει η τιμή, όχι η επιλογή.
        if sinr_db > userBestSinr
            userBestSinr      = sinr_db;
            userBestNode      = "BS" + string(b);
            userBestType      = "Terrestrial";
            userBestDistance  = range3DMat(u,b);
            userBestPathLoss  = pathLossMat(u,b);
            userBestElevation = NaN;
        end
    end

    % Στιγμιότυπο του καλύτερου υποψήφιου BS πριν τη σύγκριση με τον δορυφόρο (per-candidate διαγνωστικό).
    % NaN όταν καμία επίγεια ζεύξη δεν ήταν εντός του πεδίου ισχύος του μοντέλου.
    bestBsNodeVec(u) = userBestNode;
    if isfinite(userBestSinr)
        bestBsSinrDbVec(u) = userBestSinr;
        if userBestSinr < minUsableSinrDb
            bsUnavailReasonVec(u) = "BelowSinrFloor";
        end
    else
        % Καμία επίγεια ζεύξη δεν ήταν εντός του πεδίου ισχύος του μοντέλου.
        bsUnavailReasonVec(u) = "OutOfModelRange";
    end
    bestBsDistanceVec(u) = userBestDistance;
    bestBsPathLossVec(u) = userBestPathLoss;

    %% ===== Satellite candidate =====
    [~, elevSat, slantRangeSat] = geodetic2aer( ...
        sat_geo(1), sat_geo(2), sat_geo(3), ...
        user_geo(u,1), user_geo(u,2), user_geo(u,3), wgs84);

    satSlantRangeVec(u) = slantRangeSat;
    satElevationVec(u)  = elevSat;

    if elevSat >= satParameters.MinElevationDeg
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency;
        satPathLoss = fspl(slantRangeSat, lambdaSat);

        % Ατμοσφαιρική απόσβεση αερίων (οξυγόνο + υδρατμοί), TR 38.811 §6.6.4.
        % Βροχή/νέφωση & ιονοσφαιρική σπινθηρίδα δεν μοντελοποιούνται.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        % Shadowed Rician· i.i.d. ανά κλήση (ο δορυφόρος κινείται -> η γεωμετρία
        % σκίασης αποσυσχετίζεται γρήγορα). Η ανύψωση περιορίζεται στο πεδίο
        % ισχύος της προσαρμογής· εναλλακτικά σταθερή κατάσταση σκίασης μέσω
        % satParameters.ShadowingState (για ανάλυση ευαισθησίας).
        if isfield(satParameters, 'ShadowingState') && ~isempty(satParameters.ShadowingState)
            [b0, mNak, omega] = shadowedRicianStateParams(satParameters.ShadowingState);
        else
            elevClamped = min(max(elevSat, satFadeElevRangeDeg(1)), satFadeElevRangeDeg(2));
            [b0, mNak, omega] = shadowedRicianElevParams(elevClamped);
        end
        satPathLoss = satPathLoss - shadowedRicianFadingDb(b0, mNak, omega);

        % Ο όρος παρεμβολής είναι μηδενικός στο δορυφορικό σκέλος: 2.0 GHz
        % έναντι 3.5 GHz του επίγειου (TR 38.821 Πίν. 6.1.3.2-1 και TR 38.901),
        % άρα καμία ομοδιαυλική επικάλυψη, και ένας μόνο δορυφόρος, άρα καμία
        % παρεμβολή μεταξύ δεσμών. Το SINR ταυτίζεται εδώ με το SNR.
        satSinrDb = (satParameters.EIRP - 30) - satPathLoss - noisePowerSAT_dBW;
    else
        satPathLoss = inf;
        satSinrDb = -Inf;
        satUnavailReasonVec(u) = "NotVisible";
    end
    if satUnavailReasonVec(u) == "" && satSinrDb < minUsableSinrDb
        satUnavailReasonVec(u) = "BelowSinrFloor";
    end

    satPathLossVec(u) = satPathLoss;
    satSinrDbVec(u)    = satSinrDb;

    if satSinrDb > userBestSinr
        userBestSinr       = satSinrDb;
        userBestNode      = "SAT-1";
        userBestType      = "Satellite";
        userBestDistance  = slantRangeSat;
        userBestPathLoss  = satPathLoss;
        userBestElevation = elevSat;
    end

    % Σε outage κρατάμε τα διαγνωστικά του καλύτερου υποψηφίου, αλλά ο χρήστης
    % δεν ανατίθεται σε κόμβο.
    if userBestSinr < minUsableSinrDb
        bestNodeVec(u)     = "None";
        bestNodeTypeVec(u) = "Outage";
    else
        bestNodeVec(u)     = userBestNode;
        bestNodeTypeVec(u) = userBestType;
    end
    bestDistanceVec(u)     = userBestDistance;
    bestPathLossVec(u)     = userBestPathLoss;
    bestSinrDbVec(u)        = userBestSinr;
    bestElevationDegVec(u) = userBestElevation;
end

%% ------------------ Κατανομή πόρων, χωρητικότητα, ενέργεια ------------------
% Ο κανόνας ζει σε μία μόνο συνάρτηση (allocateAndTally), επειδή καλείται και
% για τις εναλλακτικές πολιτικές παρακάτω: αν ήταν αντιγραμμένος, οι δύο
% αντιγραφές θα μπορούσαν να αποκλίνουν χωρίς να φανεί.
alloc = struct('BW_bs', BW_bs, 'pOutPerChainW', pOutPerChainW, ...
               'targetNormalizedSe', targetNormalizedSe, 'numBs', numBs);

[capacityMbpsVec, throughputMbpsVec, serviceStateVec, nodePowerWattsVec, ...
 energyPerBitUJVec, networkEnergy] = allocateAndTally(bestNodeVec, ...
    bestNodeTypeVec, bestSinrDbVec, simParameters, satParameters, alloc);

%% ------------------ Σύγκριση πολιτικών στην ίδια πραγματοποίηση ------------------
% Τα δύο υποψήφια SINR είναι ήδη υπολογισμένα για κάθε χρήστη. Οι εναλλακτικές
% πολιτικές δεν ξαναδειγματοληπτούν τίποτα: κρατούν το ίδιο κανάλι και αλλάζουν
% μόνο ποιος υποψήφιος επιλέγεται, με την κατανομή πόρων να υπολογίζεται εκ νέου
% ώστε να αποτυπωθεί ο διαφορετικός φόρτος που προκύπτει.

% Αποκλειστικά επίγεια: κάθε χρήστης στον καλύτερο σταθμό του, αν υπάρχει και
% ξεπερνά το ελάχιστο χρησιμοποιήσιμο SINR.
terrNode = strings(numUsers,1);
terrType = strings(numUsers,1);
terrSnr  = bestBsSinrDbVec;
for u = 1:numUsers
    if isfinite(bestBsSinrDbVec(u)) && bestBsSinrDbVec(u) >= minUsableSinrDb
        terrNode(u) = bestBsNodeVec(u);
        terrType(u) = "Terrestrial";
    else
        terrNode(u) = "None";
        terrType(u) = "Outage";
    end
end

% Αποκλειστικά δορυφορική: όσοι χρήστες έχουν ορατό δορυφόρο πάνω από το
% κατώφλι μοιράζονται το εύρος ζώνης του, όπως ακριβώς και στην κανονική
% λειτουργία.
satNode = strings(numUsers,1);
satType = strings(numUsers,1);
satSnrPolicy = satSinrDbVec;
for u = 1:numUsers
    if isfinite(satSinrDbVec(u)) && satSinrDbVec(u) >= minUsableSinrDb
        satNode(u) = "SAT-1";
        satType(u) = "Satellite";
    else
        satNode(u) = "None";
        satType(u) = "Outage";
    end
end

policyKpis = struct();
policyKpis.Actual = policySummary(capacityMbpsVec, throughputMbpsVec, networkEnergy, numUsers);

[capT, thrT, ~, ~, ~, netT] = allocateAndTally(terrNode, terrType, terrSnr, ...
    simParameters, satParameters, alloc);
policyKpis.TerrestrialOnly = policySummary(capT, thrT, netT, numUsers);

[capS, thrS, ~, ~, ~, netS] = allocateAndTally(satNode, satType, satSnrPolicy, ...
    simParameters, satParameters, alloc);
policyKpis.SatelliteOnly = policySummary(capS, thrS, netS, numUsers);

%% ------------------ Κατάσταση καναλιού για την επόμενη κλήση ------------------
% Ό,τι χρειάζεται μια continuation κλήση για τη χωρική συσχέτιση του shadow fading.
newChannelState.UserGeo         = user_geo;
newChannelState.IsLOS           = losMat;
newChannelState.LosLatent       = losLatentMat;
newChannelState.ShadowFading_dB = sfMat;

end

function [capacityMbpsVec, throughputMbpsVec, serviceStateVec, ...
          nodePowerWattsVec, energyPerBitUJVec, networkEnergy] = ...
          allocateAndTally(nodeVec, typeVec, sinrDbVec, simParameters, satParameters, alloc)
% Κατανομή εύρους ζώνης, χωρητικότητα, κατάσταση υπηρεσίας, ισχύς και ενέργεια
% ανά bit, για ΔΕΔΟΜΕΝΗ ανάθεση χρηστών σε κόμβους. Απομονωμένη ώστε η ίδια
% λογική να εφαρμόζεται και στις εναλλακτικές πολιτικές.
numUsers = numel(nodeVec);
numBs    = alloc.numBs;

capacityMbpsVec   = nan(numUsers,1);
throughputMbpsVec = zeros(numUsers,1);
serviceStateVec   = strings(numUsers,1);
nodePowerWattsVec = nan(numUsers,1);
energyPerBitUJVec = inf(numUsers,1);

for u = 1:numUsers
    % Outage: μηδενική χωρητικότητα/ισχύς, ενέργεια/bit = Inf. Παραλείπονται
    % πριν το usersOnThisNode ώστε να μη μετρηθούν σαν να μοιράζονται κόμβο.
    if typeVec(u) == "Outage"
        % Καμία ενεργή ζεύξη: η χωρητικότητα δεν ορίζεται (δεν υπάρχει ζεύξη
        % να τη φέρει), η παραδοθείσα ρυθμαπόδοση είναι μηδενική, και η
        % ενέργεια ανά παραδοθέν bit απροσδιόριστη. Η κατανάλωση των κόμβων
        % συνεχίζεται και προσμετράται στο ισοζύγιο δικτύου παρακάτω.
        capacityMbpsVec(u)   = NaN;
        throughputMbpsVec(u) = 0;
        nodePowerWattsVec(u) = NaN;
        energyPerBitUJVec(u) = Inf;
        serviceStateVec(u)   = "Outage";
        continue;
    end

    servingNode = nodeVec(u);

    % Πόσοι χρήστες συνολικά εξυπηρετούνται από τον ΙΔΙΟ κόμβο
    usersOnThisNode = sum(nodeVec == servingNode);

    % Συνολικό bandwidth και κατανάλωση ισχύος του κόμβου: μοντέλο EARTH
    % (Auer et al. 2011) για BS, γραμμικό μοντέλο ενισχυτή ισχύος για δορυφόρο.
    if typeVec(u) == "Terrestrial"
        nodeBW = alloc.BW_bs;
        nodePowerW = simParameters.Power.NumTrx * ...
            (simParameters.Power.P0 + simParameters.Power.DeltaP * alloc.pOutPerChainW);
    else
        nodeBW = satParameters.Bandwidth;
        pOutW  = 10^((satParameters.TxPower - 30)/10);
        nodePowerW = satParameters.Power.Pfix + pOutW / satParameters.Power.EtaPA;
    end

    % Κατανομή πόρων (B_user = BW_grid / N_users)
    B_user = nodeBW / usersOnThisNode;

    % Χωρητικότητα Shannon, με clamp στη μέγιστη φασματική απόδοση του NR
    % (MCS 28 / 64QAM, TS 38.214 Πίν. 5.1.3.1-1) ώστε να μην υπερεκτιμάται σε υψηλό SINR.
    maxSpectralEfficiency = 5.5547;   % bits/s/Hz (MCS 28, 64QAM)
    sinr_lin = 10^(sinrDbVec(u)/10);
    spectralEfficiency = min(log2(1 + sinr_lin), maxSpectralEfficiency);
    capacity = B_user * spectralEfficiency;   % bits/s

    capacityMbpsVec(u)   = capacity * 1e-6;    % Mbps
    throughputMbpsVec(u) = capacity * 1e-6;    % ενεργή ζεύξη -> παραδίδεται

    % Κατάσταση υπηρεσίας: το κριτήριο ορίζεται επί του συνολικού εύρους
    % καναλιού, οπότε σφίγγει καθώς αυξάνεται ο φόρτος του κόμβου.
    if (capacity / nodeBW) >= alloc.targetNormalizedSe
        serviceStateVec(u) = "Served";
    else
        serviceStateVec(u) = "BelowTarget";
    end

    % Ενεργειακό proxy: ισομερής κατανομή ισχύος κόμβου ανά χρήστη (ίδια λογική
    % με το bandwidth split), διαιρεμένη με τον ρυθμό bit του χρήστη -> µJ/bit
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

%% ------------------ Ενεργειακή απόδοση σε επίπεδο δικτύου ------------------
% Η ενέργεια ανά bit της (eq. energy-per-bit) είναι ανεξάρτητη του πλήθους
% χρηστών: το L_n απλοποιείται. Ο δείκτης bit/J αθροίζει ισχύ ΑΝΑ ΚΟΜΒΟ (όχι
% ανά χρήστη) και συνολικό ρυθμό, οπότε αποτυπώνει τη σύνθεση των χρηστών.
% Δίνονται δύο εμβέλειες: πλήρης κόμβος και μόνο ενισχυτής (συμμετρική).
networkEnergy = struct();
totalRateBps      = sum(throughputMbpsVec) * 1e6;   % παραδοθέντα bits
pOutSatW          = 10^((satParameters.TxPower - 30)/10);
bsFullW           = simParameters.Power.NumTrx * (simParameters.Power.P0 + simParameters.Power.DeltaP*alloc.pOutPerChainW);
bsRfW             = simParameters.Power.NumTrx * simParameters.Power.DeltaP * alloc.pOutPerChainW;
bsSleepW          = simParameters.Power.NumTrx * simParameters.Power.Psleep;
satFullW          = satParameters.Power.Pfix + pOutSatW/satParameters.Power.EtaPA;
satRfW            = pOutSatW / satParameters.Power.EtaPA;

activeBs  = 0;
for b = 1:numBs
    if any(nodeVec == "BS" + string(b))
        activeBs = activeBs + 1;
    end
end
idleBs    = numBs - activeBs;
satActive = any(typeVec == "Satellite");

% Πλήρης εμβέλεια: ενεργοί BS κατά EARTH, αδρανείς σε Psleep. Ο δορυφόρος δεν
% έχει αντίστοιχο μέγεθος αδράνειας, οπότε προσμετράται μόνο όταν εξυπηρετεί.
networkEnergy.TotalPower_W    = activeBs*bsFullW + idleBs*bsSleepW + satActive*satFullW;
% Συμμετρική εμβέλεια: μόνο το τμήμα που εξαρτάται από τον ενισχυτή, και στα δύο σκέλη.
networkEnergy.TotalPowerRf_W  = activeBs*bsRfW + satActive*satRfW;
networkEnergy.TotalRate_bps   = totalRateBps;
networkEnergy.ActiveBs        = activeBs;
networkEnergy.IdleBs          = idleBs;
networkEnergy.SatActive       = satActive;
networkEnergy.BitPerJoule     = totalRateBps / networkEnergy.TotalPower_W;
networkEnergy.BitPerJouleRf   = totalRateBps / networkEnergy.TotalPowerRf_W;
networkEnergy.ServedUsers     = sum(serviceStateVec == "Served");
networkEnergy.BelowTargetUsers= sum(serviceStateVec == "BelowTarget");
networkEnergy.OutageUsers     = sum(serviceStateVec == "Outage");
networkEnergy.TargetNormSe    = alloc.targetNormalizedSe;
end

function s = policySummary(capacityMbpsVec, throughputMbpsVec, networkEnergy, numUsers)
% Δείκτες μιας πολιτικής για ένα βήμα, σε μορφή έτοιμη για συνάθροιση ανά
% διέλευση. Η μέση χωρητικότητα υπολογίζεται στους χρήστες που εξυπηρετούνται,
% ενώ η συνολική ρυθμαπόδοση και τα ποσοστά ορίζονται σε όλους - ώστε μια
% πολιτική να μην ωφελείται επειδή αφήνει χρήστες εκτός.
served = ~isnan(capacityMbpsVec);
s = struct();
if any(served)
    s.MeanCapacity_Mbps = mean(capacityMbpsVec(served));
else
    s.MeanCapacity_Mbps = NaN;
end
s.TotalRate_Mbps   = sum(throughputMbpsVec);
s.FracServed       = networkEnergy.ServedUsers / numUsers;
s.FracBelowTarget  = networkEnergy.BelowTargetUsers / numUsers;
s.FracOutage       = networkEnergy.OutageUsers / numUsers;
s.TotalPower_W     = networkEnergy.TotalPower_W;
s.TotalPowerRf_W   = networkEnergy.TotalPowerRf_W;
s.BitPerJoule      = networkEnergy.BitPerJoule;
s.BitPerJouleRf    = networkEnergy.BitPerJouleRf;
end

function [isLos, latent] = spatiallyConsistentLos(pLos, moveDistance, hasPrevState, prevLatent)
% Χωρικά συνεπής κατάσταση LOS/NLOS (TR 38.901 §7.6.3.3): η κατάσταση
% προκύπτει συγκρίνοντας μια χωρικά συσχετισμένη ομοιόμορφη μεταβλητή με την
% πιθανότητα LOS. Η υποκείμενη γκαουσιανή μεταβλητή εξελίσσεται ως AR(1) με
% εκθετική συσχέτιση (§7.4.4, εξ. 7.4-5) και απόσταση συσχέτισης 50 m για την
% κατάσταση LOS/NLOS (Πίν. 7.6.3.1-2, UMa & UMi).
dCorrLos = 50;
if ~hasPrevState
    latent = randn();
else
    r = exp(-moveDistance / dCorrLos);
    latent = r*prevLatent + sqrt(1 - r^2)*randn();
end
u = 0.5*erfc(-latent/sqrt(2));   % γκαουσιανή -> ομοιόμορφη στο (0,1)
isLos = u < pLos;
end

function rho = shadowFadingCorrelation(moveDistance, scenario, isLos)
% Συντελεστής αυτοσυσχέτισης σκίασης, TR 38.901 §7.4.4 εξ. (7.4-5):
% R(Δd) = exp(-Δd/d_corr). Αποστάσεις συσχέτισης από τον Πίν. 7.5-6:
% UMa LOS=37/NLOS=50 m, UMi LOS=10/NLOS=13 m.
switch scenario
    case 'UMi'
        dCorr = 10*isLos + 13*~isLos;
    otherwise
        dCorr = 37*isLos + 50*~isLos;
end
rho = exp(-moveDistance / dCorr);
end

function fadeDb = smallScaleFadingDb(isLos, KdB)
% Κέρδος small-scale fading σε dB, με E[|h|^2] = 1.
%   isLos=false -> Rayleigh: |h|^2 ~ Exp(1)
%   isLos=true  -> Rician με συντελεστή K (dB)· K->0 ανάγεται ομαλά σε Rayleigh
if ~isLos
    fadeDb = 10*log10(-log(rand()));
else
    Klin  = 10^(KdB/10);
    s     = sqrt(Klin/(Klin+1));       % πλάτος LOS συνιστώσας
    sigma = sqrt(1/(2*(Klin+1)));      % τυπ. απόκλιση ανά διάσταση scatter
    h     = (s + sigma*randn()) + 1i*(sigma*randn());
    fadeDb = 20*log10(abs(h));         % s^2 + 2*sigma^2 = 1
end
end

function isValid = isValidTerrestrialLink(d2D, hUT, fcHz)
% Πεδίο ισχύος των μοντέλων UMa/UMi (TR 38.901 Πίν. 7.4.1-1, στήλη
% "Applicability range"): 10m <= d2D <= 5km, 1.5m <= hUT <= 22.5m.
% Η συχνότητα ελέγχεται ως προς το εύρος του ίδιου του μοντέλου (§7.4.1).
isValid = d2D  >= 10   && d2D  <= 5000 && ...
          hUT  >= 1.5  && hUT  <= 22.5 && ...
          fcHz >= 0.5e9 && fcHz <= 100e9;
end

function [b0, m, omega] = shadowedRicianElevParams(elevDeg)
% Παράμετροι Shadowed Rician από τη γωνία ανύψωσης (Abdi et al. 2003, εξ. 19).
% Προσαρμογή πολυωνύμων σε πειραματικά δεδομένα, ισχύει για 20° < θ < 80°.
th = elevDeg;
b0    = -4.7943e-8*th^3 + 5.5784e-6*th^2 - 2.1344e-4*th + 3.2710e-2;
m     =  6.3739e-5*th^3 + 5.8533e-4*th^2 - 1.5973e-1*th + 3.5156;
omega =  1.4428e-5*th^3 - 2.3798e-3*th^2 + 1.2702e-1*th - 1.4864;
end

function [b0, m, omega] = shadowedRicianStateParams(state)
% Σταθερές καταστάσεις σκίασης (Abdi et al. 2003, Πίν. III) - ανάλυση ευαισθησίας.
switch lower(string(state))
    case "light"
        b0 = 0.158; m = 19.4;  omega = 1.29;
    case "average"
        b0 = 0.126; m = 10.1;  omega = 0.835;
    case "heavy"
        b0 = 0.063; m = 0.739; omega = 8.97e-4;
    otherwise
        error('shadowedRicianStateParams:UnknownState', ...
            'Άγνωστη κατάσταση σκίασης "%s" - δεκτές: "light", "average", "heavy".', state);
end
end

function fadeDb = shadowedRicianFadingDb(b0, m, omega)
% Κέρδος Shadowed Rician σε dB (Abdi et al. 2003, εξ. 1): σκεδαζόμενη
% συνιστώσα Rayleigh μέσης ισχύος 2*b0 συν συνιστώσα LOS με πλάτος
% κατά Nakagami-m μέσης ισχύος omega. E[|h|^2] = omega + 2*b0 < 1, δηλαδή
% η μέση εξασθένηση λόγω σκίασης περιέχεται στο ίδιο το μοντέλο.
losAmp  = sqrt(gammaRand(m, omega/m));            % |Z|, E[Z^2] = omega
scatter = sqrt(b0)*(randn() + 1i*randn());        % E[|A|^2] = 2*b0
fadeDb  = 20*log10(abs(losAmp + scatter));
end

function x = gammaRand(shape, scale)
% Δείγμα από κατανομή Gamma (Marsaglia & Tsang 2000). Υλοποιείται τοπικά
% ώστε να μη χρειάζεται το Statistics Toolbox· δέχεται και shape < 1.
if shape < 1
    x = gammaRand(shape + 1, scale) * rand()^(1/shape);
    return;
end
d = shape - 1/3;
c = 1/sqrt(9*d);
while true
    v = -1;
    while v <= 0
        z = randn();
        v = (1 + c*z)^3;
    end
    u = rand();
    if log(u) < 0.5*z^2 + d - d*v + d*log(v)
        x = d*v*scale;
        return;
    end
end
end

function pLos = losProbability38901(d2D, hUT, scenario)
% Πιθανότητα LOS ζεύξης BS-χρήστη (TR 38.901 Πίν. 7.4.2-1). d2D, hUT σε μέτρα.
switch scenario
    case 'UMi'
        if d2D <= 18
            pLos = 1;
        else
            pLos = 18/d2D + exp(-d2D/36) * (1 - 18/d2D);
        end
    case 'UMa'
        if d2D <= 18
            pLos = 1;
        else
            if hUT <= 13
                Cprime = 0;
            else
                g = 1.25e-6 * d2D^3 * exp(-d2D/150);
                Cprime = ((hUT - 13)/10)^1.5 * g;
            end
            pLos = (18/d2D + exp(-d2D/63) * (1 - 18/d2D)) * (1 + Cprime);
        end
    otherwise
        error('losProbability38901:UnsupportedScenario', ...
            'Άγνωστο PathLoss.Scenario "%s" - η πιθανότητα LOS (TR 38.901 §7.4.2) είναι ορισμένη μόνο για "UMa" και "UMi".', ...
            scenario);
end
end

function pla_dB = gasAttenuationSlantP676(freqHz, elevDeg)
% Απόσβεση ατμοσφαιρικών αερίων (O2 + υδρατμοί) σε ζεύξη δορυφόρου-χρήστη.
% TR 38.811 §6.6.4 εξ. (6.6-8): PLA = A_zenith / sin(ε).
% A_zenith = γ_o·h_o + γ_w·h_w· ισοδύναμα ύψη κατά ITU-R P.676-12 Annex 2·
% ειδικές αποσβέσεις γ_o/γ_w από gaspl (ITU-R P.676-13 Annex 1)·
% reference atmosphere κατά ITU-R P.835.
TcRef    = 15;        % °C (= 288.15 K)
TKRef    = 288.15;    % K
pPaRef   = 101325;    % Pa
pHpaRef  = 1013.25;   % hPa
rhoRef   = 7.5;       % g/m^3 (υδρατμοί)

freqGHz = freqHz / 1e9;

gammaDry = gaspl(1000, freqHz, TcRef, pPaRef, 0);       % dB/km, μόνο O2
gammaTot = gaspl(1000, freqHz, TcRef, pPaRef, rhoRef);  % dB/km, O2 + υδρατμοί
gammaWet = gammaTot - gammaDry;

[ho, hw] = equivalentHeightsP676(freqGHz, TKRef, pHpaRef, rhoRef);

Azenith = gammaDry*ho + gammaWet*hw;   % dB
pla_dB  = Azenith / sind(elevDeg);     % dB
end

function [ho, hw] = equivalentHeightsP676(freqGHz, T_K, p_hPa, rho)
% Ισοδύναμα ύψη O2/υδρατμών, ITU-R P.676-12 Annex 2 εξ. (30)-(38).
e_hPa = rho * T_K / 216.7;
rp = (p_hPa + e_hPa) / 1013.25;

% -- Οξυγόνο --
A_o = 0.7832 + 0.00709*(T_K - 273.15);

t1 = (5.1040 / (1 + 0.066*rp^(-2.3))) * ...
     exp(-((freqGHz - 59.7) / (2.87 + 12.4*exp(-7.9*rp)))^2);

fi3 = [118.750334 368.498246 424.763020 487.249273 715.392902 773.839490 834.145546];
ci3 = [0.1597 0.1066 0.1325 0.1242 0.0938 0.1448 0.1374];
t2 = sum( (ci3 * exp(2.12*rp)) ./ ((freqGHz - fi3).^2 + 0.025*exp(2.2*rp)) );

t3 = (0.0114*freqGHz * (15.02*freqGHz^2 - 1353*freqGHz + 5.333e4)) / ...
     ((1 + 0.14*rp^(-2.6)) * (freqGHz^3 - 151.3*freqGHz^2 + 9629*freqGHz - 6803));

ho = (6.1*A_o / (1 + 0.17*rp^(-1.1))) * (1 + t1 + t2 + t3);
if freqGHz >= 70
    ho = min(ho, 10.7*rp^0.3);
end

% -- Υδρατμοί --
A_w = 1.9298 - 0.04166*(T_K - 273.15) + 0.0517*e_hPa;
B_w = 1.1674 - 0.00622*(T_K - 273.15) + 0.0063*e_hPa;
sigma_w = 1 + 1.013 / (1 + exp(-8.6*(rp - 0.57)));

fi4 = [22.235080 183.310087 325.152888 380.197353 439.150807 448.001085 ...
       474.689092 488.490108 556.935985 620.700870 752.033113 916.171582 ...
       970.315022 987.926764];
ai4 = [1.52 7.62 1.56 4.15 0.20 1.63 0.76 0.26 7.81 1.25 16.2 1.47 1.36 1.60];
bi4 = [2.56 10.2 2.70 5.70 0.91 2.46 2.22 2.49 10.0 2.35 20.0 2.58 2.44 1.86];

hw = A_w + B_w * sum( (ai4*sigma_w) ./ ((freqGHz - fi4).^2 + bi4) );
end
