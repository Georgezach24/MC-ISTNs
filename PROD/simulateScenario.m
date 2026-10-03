function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
    newChannelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
    bsUnavailReasonVec, satUnavailReasonVec, policyKpis, newDecision] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, ...
                     prevChannelState, prevDecision)
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
% 8ο όρισμα prevDecision / 24η έξοδος newDecision: κατάσταση της απόφασης
% (εξυπηρετών κόμβος, υποψήφιος σε εκκρεμότητα, μετρητής επιβεβαίωσης) ανά
% πολιτική. Χρειάζεται επειδή η υστέρηση και ο χρόνος επιβεβαίωσης κάνουν την
% απόφαση εξαρτώμενη από το παρελθόν. Το runSimulation.m τη μεταφέρει μέσα σε
% κάθε διέλευση και τη μηδενίζει στην αρχή της επόμενης.
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
if nargin < 8
    prevDecision = [];
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

%% ------------------ Υστέρηση και χρόνος επιβεβαίωσης ------------------
% Μεταφορά του Event A3 (TS 38.331 §5.5.4.4) στον κανόνα της εργασίας. Το
% Hysteresis ορίζεται στην §6.3.2 ως INTEGER (0..30) με τιμή = πεδίο * 0.5 dB,
% δηλαδή 0 έως 15 dB· το TimeToTrigger ως απαριθμημένο σύνολο τιμών σε ms.
% Προεπιλογή 0 και 0: αναπαράγει τη συμπεριφορά χωρίς υστέρηση.
if isfield(simParameters, 'Mobility')
    hysDb    = simParameters.Mobility.HysteresisDb;
    tttMs    = simParameters.Mobility.TimeToTriggerMs;
    dtSecHys = simParameters.Mobility.DtSeconds;
else
    hysDb = 0; tttMs = 0; dtSecHys = 1;
end
% Ο χρόνος επιβεβαίωσης μετατρέπεται σε διαδοχικά βήματα. Η στρογγυλοποίηση
% είναι προς τα πάνω επειδή το πρότυπο απαιτεί η συνθήκη να ΙΣΧΥΕΙ για τη
% διάρκεια αυτή· με Δt = 1 s οι τιμές κάτω από 1024 ms δεν διακρίνονται.
nTtt = max(1, ceil(tttMs / (dtSecHys*1000)));

%% ------------------ Μέτρηση: μέση τιμή και φίλτρο L3 ------------------
% Το TS 38.331 §5.5.3.2 ορίζει ότι ο A3 δεν συγκρίνει τη στιγμιαία μέτρηση
% αλλά τη φιλτραρισμένη, F_n = (1-a)*F_{n-1} + a*M_n ("filter the measured
% result, BEFORE using for evaluation of reporting criteria"), και το ίδιο το
% M_n δεν είναι στιγμιαίο δείγμα: είναι αποτέλεσμα περιόδου μέτρησης
% τουλάχιστον 200 ms στο FR1 (TS 38.133 Πίν. 9.2.5.2-1). Αντίγραφα των δύο
% εδαφίων: Sources/TS38331_L3_filtering.md.
%
%   Meas.NAvg    πλήθος ανεξάρτητων πραγματώσεων γρήγορων διαλείψεων που
%                μεσοποιούνται μέσα στο βήμα. Προκύπτει από την αυτοσυσχέτιση
%                του Clarke και υπολογίζεται στο runSimulation.m.
%   Meas.L3Alpha ο συντελεστής a, ΗΔΗ προσαρμοσμένος στο βήμα της προσομοίωσης
%                κατά τη ρήτρα "adapt the filter such that the time
%                characteristics ... are preserved at different input rates".
%                Τιμή 1 = χωρίς φίλτρο (NOTE 1, k = 0).
%
% Προεπιλογές 1 και 1: αναπαράγουν ακριβώς τη συμπεριφορά πριν τη μέτρηση,
% ώστε τα scripts επαλήθευσης και η παλιά ρύθμιση να μένουν αναπαραγώγιμα.
measNAvgTerr = max(1, round(measParam(simParameters, 'NAvg',    1)));
measNAvgSat  = max(1, round(measParam(satParameters, 'NAvg',    1)));
l3Alpha      = min(max(measParam(simParameters, 'L3Alpha', 1), 0), 1);

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestBsNodeVec       = strings(numUsers,1);   % ποιος BS κέρδισε, για τη σύγκριση πολιτικών
policyNamesLocal    = {'Actual','TerrestrialOnly','SatelliteOnly'};
newDecision = struct();
for pol = policyNamesLocal
    newDecision.(pol{1}) = struct('Node', strings(numUsers,1), ...
        'Pending', strings(numUsers,1), 'Count', zeros(numUsers,1));
end
hasPrevDecision = ~isempty(prevDecision);
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

% Καταστάσεις εναλλακτικών πολιτικών (στήλες 2 και 3: επίγεια / δορυφορική)
polNodeMat = strings(numUsers,3);
polTypeMat = strings(numUsers,3);
polSinrMat = nan(numUsers,3);

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

% Κατάσταση του φίλτρου L3 ανά υποψήφιο. NaN σημαίνει "δεν υπάρχει
% προηγούμενη μέτρηση", οπότε F_0 = M_1 όπως ορίζει η §5.5.3.2. Ταξιδεύει
% μέσα στο prevChannelState, δηλαδή μηδενίζεται στην αρχή κάθε διέλευσης μαζί
% με το υπόλοιπο κανάλι. Η μέτρηση είναι ΚΟΙΝΗ για τις τρεις πολιτικές: δεν
% εξαρτάται από το ποιος κόμβος εξυπηρετεί, άρα η σύγκριση κατά ζεύγη μένει
% έγκυρη.
hasPrevFilt = hasPrevState && isfield(prevChannelState, 'FiltSinrBs_dB') && ...
    isequal(size(prevChannelState.FiltSinrBs_dB), [numUsers, numBs]);
if hasPrevFilt
    prevFiltBs  = prevChannelState.FiltSinrBs_dB;
    prevFiltSat = prevChannelState.FiltSinrSat_dB;
else
    prevFiltBs  = nan(numUsers, numBs);
    prevFiltSat = nan(numUsers, 1);
end
filtSinrBsMat  = nan(numUsers, numBs);
filtSinrSatVec = nan(numUsers, 1);

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
        % Ο χρόνος συνοχής (~43 ms στα 3.5 GHz με πεζό χρήστη) είναι πολύ
        % μικρότερος από το βήμα, οπότε μέσα σε ένα βήμα χωρούν πολλές
        % ανεξάρτητες πραγματώσεις. Η μέτρηση είναι η μέση τιμή τους σε
        % γραμμική ισχύ, όχι μία από αυτές (TS 38.133 Πίν. 9.2.5.2-1).
        pathLoss = pathLoss - smallScaleFadingDb(isLos, terrKdBLos, measNAvgTerr);
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
        filtSinrBsMat(u,b) = l3Filter(sinr_db, prevFiltBs(u,b), l3Alpha);

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
        satPathLoss = satPathLoss - shadowedRicianFadingDb(b0, mNak, omega, measNAvgSat);

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
    if isfinite(satSinrDb)
        filtSinrSatVec(u) = l3Filter(satSinrDb, prevFiltSat(u), l3Alpha);
    end

    if satSinrDb > userBestSinr
        userBestSinr       = satSinrDb;
        userBestNode      = "SAT-1";
        userBestType      = "Satellite";
        userBestDistance  = slantRangeSat;
        userBestPathLoss  = satPathLoss;
        userBestElevation = elevSat;
    end

    %% ===== Απόφαση με υστέρηση και χρόνο επιβεβαίωσης =====
    % Κατάλογος υποψηφίων: όλοι οι σταθμοί εντός πεδίου ισχύος συν ο δορυφόρος.
    % Η υστέρηση εφαρμόζεται σε ΟΛΟΥΣ, ώστε να καλύπτει και την εναλλαγή μεταξύ
    % επίγειων σταθμών, που αποτελεί το μεγαλύτερο μέρος του φαινομένου.
    nCand = 0;
    candNode = strings(numBs+1,1); candType = strings(numBs+1,1);
    candSinr = -inf(numBs+1,1);    candDist = nan(numBs+1,1);
    candSinrMeas = -inf(numBs+1,1);
    candPl   = nan(numBs+1,1);     candElev = nan(numBs+1,1);
    for b = 1:numBs
        if isfinite(sinrDbMat(u,b))
            nCand = nCand + 1;
            candNode(nCand) = "BS" + string(b);
            candType(nCand) = "Terrestrial";
            candSinr(nCand) = sinrDbMat(u,b);
            candSinrMeas(nCand) = filtSinrBsMat(u,b);
            candDist(nCand) = range3DMat(u,b);
            candPl(nCand)   = pathLossMat(u,b);
        end
    end
    if isfinite(satSinrDb)
        nCand = nCand + 1;
        candNode(nCand) = "SAT-1";
        candType(nCand) = "Satellite";
        candSinr(nCand) = satSinrDb;
        candSinrMeas(nCand) = filtSinrSatVec(u);
        candDist(nCand) = slantRangeSat;
        candPl(nCand)   = satPathLoss;
        candElev(nCand) = elevSat;
    end
    candNode = candNode(1:nCand); candType = candType(1:nCand);
    candSinr = candSinr(1:nCand); candDist = candDist(1:nCand);
    candSinrMeas = candSinrMeas(1:nCand);
    candPl   = candPl(1:nCand);   candElev = candElev(1:nCand);

    isTerr = (candType == "Terrestrial");
    isSat  = (candType == "Satellite");

    for pI = 1:numel(policyNamesLocal)
        polName = policyNamesLocal{pI};
        switch polName
            case 'Actual',          keepMask = true(nCand,1);
            case 'TerrestrialOnly', keepMask = isTerr;
            otherwise,              keepMask = isSat;
        end

        pN = ""; pP = ""; pC = 0;
        if hasPrevDecision
            pN = prevDecision.(polName).Node(u);
            pP = prevDecision.(polName).Pending(u);
            pC = prevDecision.(polName).Count(u);
        end

        [ci, pendN, pendC] = applyHysteresis(candNode, candSinr, candSinrMeas, ...
            keepMask, pN, pP, pC, minUsableSinrDb, hysDb, nTtt);

        if ci == 0
            selNode = "None"; selType = "Outage"; selSinr = NaN;
        else
            selNode = candNode(ci); selType = candType(ci); selSinr = candSinr(ci);
        end
        newDecision.(polName).Node(u)    = selNode;
        newDecision.(polName).Pending(u) = pendN;
        newDecision.(polName).Count(u)   = pendC;

        if strcmp(polName, 'Actual')
            if ci == 0
                % Σε outage κρατάμε τα διαγνωστικά του καλύτερου υποψηφίου,
                % αλλά ο χρήστης δεν ανατίθεται σε κόμβο.
                bestNodeVec(u)         = "None";
                bestNodeTypeVec(u)     = "Outage";
                bestDistanceVec(u)     = userBestDistance;
                bestPathLossVec(u)     = userBestPathLoss;
                bestSinrDbVec(u)       = userBestSinr;
                bestElevationDegVec(u) = userBestElevation;
            else
                % Με υστέρηση ο επιλεγμένος κόμβος δεν είναι κατ' ανάγκη ο
                % καλύτερος: η χωρητικότητα υπολογίζεται από το SINR του
                % κόμβου στον οποίο ο χρήστης είναι όντως συνδεδεμένος.
                bestNodeVec(u)         = candNode(ci);
                bestNodeTypeVec(u)     = candType(ci);
                bestDistanceVec(u)     = candDist(ci);
                bestPathLossVec(u)     = candPl(ci);
                bestSinrDbVec(u)       = candSinr(ci);
                bestElevationDegVec(u) = candElev(ci);
            end
        else
            polNodeMat(u,pI) = selNode;
            polTypeMat(u,pI) = selType;
            polSinrMat(u,pI) = selSinr;
        end
    end
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

% Οι δύο εναλλακτικές πολιτικές αποφασίστηκαν παραπάνω με τον ΙΔΙΟ κανόνα
% υστέρησης, καθεμία με τη δική της κατάσταση: μια πολιτική που αλλάζει κόμβο
% δεν πρέπει να συγκρίνεται με μία που δεν αλλάζει.
terrNode = polNodeMat(:,2);
terrType = polTypeMat(:,2);
terrSnr  = polSinrMat(:,2);

satNode = polNodeMat(:,3);
satType = polTypeMat(:,3);
satSnrPolicy = polSinrMat(:,3);

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
newChannelState.FiltSinrBs_dB   = filtSinrBsMat;
newChannelState.FiltSinrSat_dB  = filtSinrSatVec;

end

function [chosenIdx, pendNode, pendCount] = applyHysteresis(candNode, candSinr, ...
    candSinrMeas, keepMask, prevNode, prevPending, prevCount, minUsableSinrDb, hysDb, nTtt)
% Απόφαση επιλογής κόμβου με υστέρηση και χρόνο επιβεβαίωσης, κατά τη μεταφορά
% του Event A3 του TS 38.331 §5.5.4.4:
%   είσοδος (A3-1):  M_n - Hys > M_p     έξοδος (A3-2):  M_n + Hys < M_p
% με Off = Ofn = Ocn = Ofp = Ocp = 0, δηλαδή χωρίς μετατοπίσεις ανά κυψέλη ή
% ανά συχνότητα, και με τη συνθήκη εισόδου να πρέπει να ισχύει για nTtt
% διαδοχικά βήματα. Δεν πρόκειται για υλοποίηση της διαδικασίας RRC: δεν
% υπάρχει σηματοδοσία, αναφορά μέτρησης, ούτε έλεγχος επιτυχίας.
%
% Δύο διαφορετικά SINR, σκόπιμα:
%   candSinr     το πραγματικό SINR της ζεύξης. Καθορίζει ΜΟΝΟ αν η ζεύξη
%                μπορεί να σηκώσει δεδομένα, δηλαδή τον έλεγχο κατωφλίου.
%                Το κατώφλι είναι φυσικό όριο (MCS 0), όχι μέτρηση.
%   candSinrMeas η φιλτραρισμένη μέτρηση κατά §5.5.3.2. Καθορίζει ΠΟΙΟΝ
%                υποψήφιο βλέπει ο μηχανισμός απόφασης ως καλύτερο και αν
%                ικανοποιείται η συνθήκη A3. Αυτό βλέπει ένα πραγματικό
%                τερματικό· δεν βλέπει ποτέ στιγμιαίο δείγμα.
% Με l3Alpha = 1 τα δύο ταυτίζονται και η συμπεριφορά ανάγεται στην παλιά.
usable = keepMask & (candSinr >= minUsableSinrDb);
chosenIdx = 0; pendNode = ""; pendCount = 0;

if ~any(usable)
    return;   % κανένας επιλέξιμος υποψήφιος -> εκτός κάλυψης
end

masked = candSinrMeas; masked(~usable) = -Inf;
[~, iBest] = max(masked);

iPrev = 0;
if prevNode ~= "" && prevNode ~= "None"
    hit = find(candNode == prevNode & usable, 1);
    if ~isempty(hit)
        iPrev = hit;
    end
end

if iPrev == 0
    % Δεν υπάρχει εξυπηρετών κόμβος, είτε επειδή είναι το πρώτο βήμα είτε
    % επειδή η ζεύξη του έπεσε κάτω από το κατώφλι. Η υστέρηση καθυστερεί τη
    % μετάβαση προς καλύτερο κόμβο· δεν κρατά μια ζεύξη που έχει πάψει να
    % είναι αξιοποιήσιμη. Η επιλογή γίνεται άμεσα.
    chosenIdx = iBest;
    return;
end

if iBest == iPrev
    chosenIdx = iPrev;
    return;
end

if candSinrMeas(iBest) - hysDb > candSinrMeas(iPrev)
    % Η συνθήκη εισόδου ισχύει: ο μετρητής επιβεβαίωσης προχωρά μόνο αν ο
    % υποψήφιος είναι ο ίδιος με το προηγούμενο βήμα.
    if prevPending == candNode(iBest)
        pendCount = prevCount + 1;
    else
        pendCount = 1;
    end
    if pendCount >= nTtt
        chosenIdx = iBest;
        pendNode  = "";
        pendCount = 0;
    else
        chosenIdx = iPrev;
        pendNode  = candNode(iBest);
    end
else
    % Η συνθήκη έπαψε να ισχύει: ο μετρητής μηδενίζεται, όπως απαιτεί η
    % απαίτηση του προτύπου να ισχύει η συνθήκη συνεχώς για TimeToTrigger.
    chosenIdx = iPrev;
end
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

function fadeDb = smallScaleFadingDb(isLos, KdB, nAvg)
% Κέρδος small-scale fading σε dB, με E[|h|^2] = 1.
%   isLos=false -> Rayleigh: |h|^2 ~ Exp(1)
%   isLos=true  -> Rician με συντελεστή K (dB)· K->0 ανάγεται ομαλά σε Rayleigh
%
% nAvg: η τιμή που βγαίνει δεν είναι στιγμιαίο δείγμα αλλά μέση τιμή σε
% ΓΡΑΜΜΙΚΗ ισχύ πάνω σε nAvg ανεξάρτητες πραγματώσεις μέσα στο βήμα, επειδή
% η μέτρηση ορίζεται πάνω σε περίοδο μέτρησης και όχι σε στιγμή
% (TS 38.133 Πίν. 9.2.5.2-1). nAvg = 1 επιστρέφει το στιγμιαίο δείγμα.
if nargin < 3 || isempty(nAvg)
    nAvg = 1;
end
% Η περίπτωση nAvg = 1 κρατά κατά γράμμα την προηγούμενη έκφραση: το
% 10*log10(mean(x)) και το 20*log10(abs(h)) είναι μαθηματικά ταυτόσημα με
% τα αντίστοιχα, αλλά διαφέρουν κατά ένα ulp. Έτσι ο έλεγχος ουδετερότητας
% δίνει ταυτόσημο αρχείο και όχι «σχεδόν ταυτόσημο».
if ~isLos
    if nAvg == 1
        fadeDb = 10*log10(-log(rand()));
        return;
    end
    powLin = mean(-log(rand(1,nAvg)));
else
    Klin  = 10^(KdB/10);
    s     = sqrt(Klin/(Klin+1));       % πλάτος LOS συνιστώσας
    sigma = sqrt(1/(2*(Klin+1)));      % τυπ. απόκλιση ανά διάσταση scatter
    h     = (s + sigma*randn(1,nAvg)) + 1i*(sigma*randn(1,nAvg));
    if nAvg == 1
        fadeDb = 20*log10(abs(h));     % s^2 + 2*sigma^2 = 1
        return;
    end
    powLin = mean(abs(h).^2);
end
fadeDb = 10*log10(powLin);
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

function fadeDb = shadowedRicianFadingDb(b0, m, omega, nAvg)
% Κέρδος Shadowed Rician σε dB (Abdi et al. 2003, εξ. 1): σκεδαζόμενη
% συνιστώσα Rayleigh μέσης ισχύος 2*b0 συν συνιστώσα LOS με πλάτος
% κατά Nakagami-m μέσης ισχύος omega. E[|h|^2] = omega + 2*b0 < 1, δηλαδή
% η μέση εξασθένηση λόγω σκίασης περιέχεται στο ίδιο το μοντέλο.
%
% nAvg: μέση τιμή της μέτρησης μέσα στο βήμα, όπως και στο επίγειο σκέλος.
% Η συνιστώσα σκίασης κληρώνεται ΜΙΑ φορά και μένει σταθερή μέσα στο βήμα:
% η σκίαση δεν είναι γρήγορο φαινόμενο και δεν εξομαλύνεται από την περίοδο
% μέτρησης. Μέσος όρος παίρνεται μόνο στη σκεδαζόμενη συνιστώσα.
if nargin < 4 || isempty(nAvg)
    nAvg = 1;
end
losAmp  = sqrt(gammaRand(m, omega/m));            % |Z|, E[Z^2] = omega
scatter = sqrt(b0)*(randn(1,nAvg) + 1i*randn(1,nAvg));   % E[|A|^2] = 2*b0
if nAvg == 1
    fadeDb = 20*log10(abs(losAmp + scatter));   % ίδια έκφραση με πριν (βλ. smallScaleFadingDb)
else
    fadeDb = 10*log10(mean(abs(losAmp + scatter).^2));
end
end

function f = l3Filter(mNow, fPrev, alpha)
% Φίλτρο Layer 3, TS 38.331 §5.5.3.2:  F_n = (1 - a)*F_{n-1} + a*M_n
% Εφαρμόζεται σε dB, όπως απαιτεί η NOTE 2 ("logarithmic filtering for
% logarithmic measurements"). Χωρίς προηγούμενη μέτρηση ισχύει F_0 = M_1,
% ρητά στο ίδιο εδάφιο. Το alpha έρχεται ήδη προσαρμοσμένο στο βήμα της
% προσομοίωσης· alpha = 1 σημαίνει καθόλου φίλτρο (NOTE 1, k = 0).
if ~isfinite(fPrev)
    f = mNow;
else
    f = (1 - alpha)*fPrev + alpha*mNow;
end
end

function v = measParam(p, name, dflt)
% Ανάγνωση πεδίου από το προαιρετικό υπο-struct Meas, με προεπιλογή.
v = dflt;
if isfield(p, 'Meas') && isfield(p.Meas, name) && ~isempty(p.Meas.(name))
    v = p.Meas.(name);
end
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
