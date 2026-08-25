function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
    newChannelState] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, prevChannelState)
% Υπολογίζει, για κάθε χρήστη, τον καλύτερο κόμβο εξυπηρέτησης (BS ή δορυφόρο)
% βάσει SNR και την επιτευχθείσα χωρητικότητα Shannon μετά την κατανομή
% εύρους ζώνης. Εξάγει το βασικό μονοπάτι υπολογισμού από το test_simulation.m
% ώστε να μπορεί να κληθεί επανειλημμένα (π.χ. από έναν Monte-Carlo driver).
%
% Πέρα από τον τελικό (νικητή) κόμβο, επιστρέφει και τα per-candidate
% διαγνωστικά (καλύτερο BS ανεξάρτητα από το αν κέρδισε, και δορυφόρος)
% ώστε ένα μελλοντικό μοντέλο ML (Part 2) να μπορεί να εκπαιδευτεί στη
% σύγκριση των δύο υποψήφιων ζεύξεων αντί να διαβάζει απευθείας τον νικητή.
%
% prevChannelState (προαιρετικό, 7ο όρισμα): αν δοθεί, η γεωμετρική
% απόσταση μετακίνησης κάθε χρήστη από την προηγούμενη κλήση χρησιμοποιείται
% ώστε το shadow fading (και η κατάσταση LOS/NLOS) της κάθε ζεύξης BS-χρήστη
% να ΣΥΣΧΕΤΙΖΕΤΑΙ με την προηγούμενη τιμή αντί να επαναδειγματίζεται ανεξάρτητα
% (βλ. τοπική συνάρτηση correlatedLosState παρακάτω για τα βιβλιογραφικά
% θεμέλια). Αν παραλειφθεί, κάθε κλήση παράγει ανεξάρτητο δείγμα
% καναλιού (i.i.d.) όπως πριν - η συμπεριφορά αυτή είναι η σωστή για callers
% όπου κάθε κλήση αντιπροσωπεύει ένα νέο, ανεξάρτητο "drop" (π.χ.
% kpiRepeatedRuns.m, monteCarloDriver.m, test_simulation.m). Μόνο callers που
% προσομοιώνουν διαδοχικές μεταδόσεις της ΙΔΙΑΣ, ουσιαστικά ακίνητης τοπολογίας
% (π.χ. temporalPassSimulation.m) πρέπει να περνάνε το newChannelState της
% προηγούμενης κλήσης ως prevChannelState στην επόμενη.
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

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SNR (κατάσταση outage) ------------------
% Ο κανόνας επιλογής κόμβου (max-SNR) ανέθετε πάντα τον "λιγότερο κακό"
% υποψήφιο, ακόμα κι όταν το SNR του ήταν αδικαιολόγητα κακό (π.χ. -600dB
% σε NLOS ζεύξη εκτός εμβέλειας) - βλ. CLAUDE.md § Standards, gap "no
% outage state". Το ελάχιστο χρησιμοποιήσιμο SNR ορίζεται εδώ ως το
% Shannon-ισοδύναμο SNR για το πιο ανθεκτικό MCS που ορίζει το ίδιο
% πρότυπο ήδη χρησιμοποιούμενο για το capacity cap παραπάνω/παρακάτω
% (3GPP TS 38.214 v17.x, Πίνακας 5.1.3.1-2 "MCS Index Table 2 for PDSCH",
% MCS 0 -> QPSK, target code rate 120/1024 -> φασματική απόδοση
% 0.2344 bits/s/Hz): κάτω από αυτό το SNR, ούτε το πιο ανθεκτικό σχήμα
% διαμόρφωσης/κωδικοποίησης που ορίζει το NR δεν είναι θεωρητικά εφικτό.
minSpectralEfficiency = 0.2344;                        % bits/s/Hz, TS 38.214 §5.1.3.1, MCS 0
minUsableSnrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestDistanceVec     = nan(numUsers,1);
bestPathLossVec     = nan(numUsers,1);
bestSnrDbVec        = nan(numUsers,1);
capacityMbpsVec     = nan(numUsers,1);
nodePowerWattsVec   = nan(numUsers,1);
energyPerBitUJVec   = nan(numUsers,1);
bestElevationDegVec = nan(numUsers,1);
bestBsSnrDbVec      = nan(numUsers,1);
bestBsDistanceVec   = nan(numUsers,1);
bestBsPathLossVec   = nan(numUsers,1);

% Διαγνωστικοί πίνακες
groundDistanceMat = nan(numUsers,numBs);
range3DMat        = nan(numUsers,numBs);
pathLossMat       = nan(numUsers,numBs);
snrDbMat          = nan(numUsers,numBs);
pLosMat           = nan(numUsers,numBs);
losMat            = false(numUsers,numBs);
sfMat             = nan(numUsers,numBs);
satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSnrDbVec       = nan(numUsers,1);

hasPrevState = ~isempty(prevChannelState) && ...
    isequal(size(prevChannelState.IsLOS), [numUsers, numBs]);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SNR) ------------------
for u = 1:numUsers
    % Αρχικοποιήσεις
    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

    % Απόσταση μετακίνησης του χρήστη από την προηγούμενη κλήση (0 αν ο
    % χρήστης είναι ακίνητος μεταξύ διαδοχικών κλήσεων, όπως συμβαίνει σήμερα
    % στο temporalPassSimulation.m - μόνο ο δορυφόρος κινείται εκεί). Καθορίζει
    % πόσο "θυμάται" το shadow fading την προηγούμενη τιμή του (βλ. τοπική
    % συνάρτηση correlatedLosAndShadowFading).
    if hasPrevState
        userMoveDistance = distance(prevChannelState.UserGeo(u,1), prevChannelState.UserGeo(u,2), ...
                                     user_geo(u,1), user_geo(u,2), wgs84);
    else
        userMoveDistance = Inf; % καμία προηγούμενη κατάσταση -> ανεξάρτητο δείγμα, όπως πριν
    end

    %% ===== Terrestrial BS candidates =====
    for b = 1:numBs
        lat0 = bs_geo(b,1);
        lon0 = bs_geo(b,2);
        h0   = 0;

        [xBS, yBS, zBS] = geodetic2enu(bs_geo(b,1), bs_geo(b,2), bs_geo(b,3), ...
                                       lat0, lon0, h0, wgs84);
        [xUE, yUE, zUE] = geodetic2enu(user_geo(u,1), user_geo(u,2), user_geo(u,3), ...
                                       lat0, lon0, h0, wgs84);

        txPosition = [xBS; yBS; zBS];
        rxPosition = [xUE; yUE; zUE];

        groundDistance = distance(bs_geo(b,1), bs_geo(b,2), ...
                                  user_geo(u,1), user_geo(u,2), wgs84);
        d3d = norm(rxPosition - txPosition);
        groundDistanceMat(u,b) = groundDistance;
        range3DMat(u,b)        = d3d;

        % LOS ανά ζεύξη βάσει πιθανότητας απόστασης (3GPP TR 38.901 §7.4.2,
        % Πίνακας 7.4.2-1), αντί για μία σταθερή global τιμή LOS. Αν υπάρχει
        % προηγούμενη κατάσταση καναλιού (prevChannelState), η κατάσταση
        % LOS/NLOS και το shadow fading διατηρούν χωρική συσχέτιση με την
        % προηγούμενη κλήση αντί να επαναδειγματίζονται ανεξάρτητα - βλ.
        % correlatedLosState παρακάτω.
        pLos = losProbability38901(groundDistance, user_geo(u,3), simParameters.PathLoss.Scenario);
        if hasPrevState
            prevIsLos = prevChannelState.IsLOS(u,b);
            prevSF    = prevChannelState.ShadowFading_dB(u,b);
        else
            prevIsLos = false;
            prevSF    = 0;
        end
        [isLos, rho, useCorrelatedSF] = correlatedLosState(pLos, userMoveDistance, ...
            simParameters.PathLoss.Scenario, hasPrevState, prevIsLos);
        pLosMat(u,b) = pLos;
        losMat(u,b)  = isLos;

        [pathLoss, sigmaSF] = nrPathLoss(simParameters.PathLoss, ...
                              simParameters.CarrierFrequency, ...
                              isLos, ...
                              txPosition, rxPosition);

        % Shadow fading: log-normal δείγμα με τυπική απόκλιση sigmaSF (TR 38.901
        % §7.4.1). Αν η ζεύξη διατήρησε την ίδια κατάσταση LOS/NLOS από την
        % προηγούμενη κλήση, το δείγμα συσχετίζεται χωρικά με το προηγούμενο
        % κατά Gudmundson (1991, exponential autocorrelation, ρ όπως
        % υπολογίστηκε στο correlatedLosState) - αλλιώς είναι ανεξάρτητο
        % (νέο "drop" ή μετάβαση LOS<->NLOS, που ούτως ή άλλως ακυρώνει τη
        % στατιστική βάση της προηγούμενης τιμής).
        if useCorrelatedSF
            sfSample = rho*prevSF + sqrt(1 - rho^2) * sigmaSF * randn();
        else
            sfSample = sigmaSF * randn();
        end
        sfMat(u,b) = sfSample;
        pathLoss = pathLoss + sfSample;
        pathLossMat(u,b) = pathLoss;

        snr_db = (simParameters.EIRP - 30) - pathLoss - noisePowerBS_dBW;
        snrDbMat(u,b) = snr_db;

        if snr_db > userBestSNR
            userBestSNR       = snr_db;
            userBestNode      = "BS" + string(b);
            userBestType      = "Terrestrial";
            userBestDistance  = d3d;
            userBestPathLoss  = pathLoss;
            userBestElevation = NaN;
        end
    end

    % Στιγμιότυπο του καλύτερου υποψήφιου BS πριν συγκριθεί με τον δορυφόρο
    % (per-candidate διαγνωστικό, ανεξάρτητο από τον τελικό νικητή - βλ.
    % bestBsSnrDbVec/bestBsDistanceVec/bestBsPathLossVec στην έξοδο).
    bestBsSnrDbVec(u)   = userBestSNR;
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

        % Ατμοσφαιρική απόσβεση αερίων (οξυγόνο + υδρατμοί), TR 38.811 §6.6.4
        % εξ. (6.6-8) - βλ. τοπική συνάρτηση gasAttenuationSlantP676 παρακάτω
        % για τις πλήρεις βιβλιογραφικές αναφορές. Ισχύει για elevSat >= 5°
        % (όριο ισχύος ITU-R P.676 Annex 2) - πάντα αληθές εδώ αφού
        % MinElevationDeg=10°. Βροχή/νέφωση παραλείπεται σκόπιμα: το TR 38.811
        % §6.6.5 τη χαρακτηρίζει αμελητέα κάτω από 6GHz (εδώ 2.01GHz). Η
        % ιονοσφαιρική σπινθηρίδα (§6.6.6) παραμένει εκτός πεδίου - χρειάζεται
        % κλιματολογικό μοντέλο (γεωγρ. πλάτος/ώρα/εποχή/ηλιακή δραστηριότητα),
        % όχι κλειστή αναλυτική σχέση - future work, όχι σιωπηλή προσέγγιση.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        satSnrDb = (satParameters.EIRP - 30) - satPathLoss - noisePowerSAT_dBW;
    else
        satPathLoss = inf;
        satSnrDb = -Inf;
    end

    satPathLossVec(u) = satPathLoss;
    satSnrDbVec(u)    = satSnrDb;

    if satSnrDb > userBestSNR
        userBestSNR       = satSnrDb;
        userBestNode      = "SAT-1";
        userBestType      = "Satellite";
        userBestDistance  = slantRangeSat;
        userBestPathLoss  = satPathLoss;
        userBestElevation = elevSat;
    end

    % Αποθήκευση επιλογής κόμβου για τον χρήστη. Ακόμα κι όταν το καλύτερο
    % διαθέσιμο SNR δεν φτάνει το ελάχιστο χρησιμοποιήσιμο όριο, τα
    % διαγνωστικά (distance/pathloss/SNR/elevation) του "λιγότερο κακού"
    % υποψηφίου διατηρούνται (χρήσιμα για ανάλυση), απλά ο χρήστης δεν
    % ανατίθεται πλέον σε αυτόν ως εξυπηρετητή.
    if userBestSNR < minUsableSnrDb
        bestNodeVec(u)     = "None";
        bestNodeTypeVec(u) = "Outage";
    else
        bestNodeVec(u)     = userBestNode;
        bestNodeTypeVec(u) = userBestType;
    end
    bestDistanceVec(u)     = userBestDistance;
    bestPathLossVec(u)     = userBestPathLoss;
    bestSnrDbVec(u)        = userBestSNR;
    bestElevationDegVec(u) = userBestElevation;
end

%% ------------------ Υπολογισμός Χωρητικότητας & Ενέργειας (Κατανομή Πόρων) ------------------
for u = 1:numUsers
    % Οι χρήστες σε outage δεν εξυπηρετούνται από κανέναν κόμβο -
    % μηδενική χωρητικότητα, καμία κατανάλωση ισχύος να τους αποδοθεί, και
    % ενέργεια/bit μη ορισμένη (Inf). Παραλείπονται ΠΡΙΝ το usersOnThisNode
    % ώστε να μην μετρηθούν σαν να "μοιράζονται" τον ίδιο κόμβο μεταξύ
    % τους μέσω του κοινού sentinel ServingNode="None".
    if bestNodeTypeVec(u) == "Outage"
        capacityMbpsVec(u)   = 0;
        nodePowerWattsVec(u) = 0;
        energyPerBitUJVec(u) = Inf;
        continue;
    end

    servingNode = bestNodeVec(u);

    % Πόσοι χρήστες συνολικά εξυπηρετούνται από τον ΙΔΙΟ κόμβο
    usersOnThisNode = sum(bestNodeVec == servingNode);

    % Επιλέγουμε το συνολικό Bandwidth του κόμβου και την κατανάλωση ισχύος του
    % (μοντέλο EARTH για BS, γραμμικό μοντέλο ενισχυτή ισχύος για δορυφόρο -
    % βλ. CLAUDE.md § Standards & scientific grounding)
    if bestNodeTypeVec(u) == "Terrestrial"
        nodeBW = BW_bs;
        pOutW  = 10^((simParameters.TxPower - 30)/10);
        nodePowerW = simParameters.Power.NumTrx * ...
            (simParameters.Power.P0 + simParameters.Power.DeltaP * pOutW);
    else
        nodeBW = satParameters.Bandwidth;
        pOutW  = 10^((satParameters.TxPower - 30)/10);
        nodePowerW = satParameters.Power.Pfix + pOutW / satParameters.Power.EtaPA;
    end

    % Κατανομή πόρων (B_user = BW_grid / N_users)
    B_user = nodeBW / usersOnThisNode;

    % Υπολογισμός τελικής χωρητικότητας για το κομμάτι του B_user: η θεωρητική
    % χωρητικότητα Shannon περιορίζεται (clamp) στη μέγιστη φασματική απόδοση
    % που ορίζει το 3GPP για πραγματικό gNB/UE, ώστε να μην υπερεκτιμάται η
    % χωρητικότητα σε υψηλά SNR (3GPP TS 38.214 v17.x, Πίνακας 5.1.3.1-2
    % "MCS Index Table 2 for PDSCH", MCS 27 -> 256QAM, target code rate
    % 948/1024 -> φασματική απόδοση 5.5547 bits/s/Hz· ίδια ανώτατη τιμή στον
    % Πίνακα 5.2.2.1-4 "4-bit CQI Table", CQI index 15). Πάνω από αυτό το SNR
    % ένας πραγματικός δέκτης δεν μπορεί να εκμεταλλευτεί άλλα bits/s/Hz αφού
    % το 256QAM είναι ο υψηλότερος διαμορφωτής που ορίζει το NR.
    maxSpectralEfficiency = 5.5547;   % bits/s/Hz, TS 38.214 §5.1.3.1/§5.2.2.1
    snr_lin = 10^(bestSnrDbVec(u)/10);
    spectralEfficiency = min(log2(1 + snr_lin), maxSpectralEfficiency);
    capacity = B_user * spectralEfficiency;   % bits/s

    capacityMbpsVec(u) = capacity * 1e-6;    % Mbps

    % Ενεργειακό proxy: ισομερής κατανομή ισχύος κόμβου ανά χρήστη (ίδια λογική
    % με το bandwidth split), διαιρεμένη με τον ρυθμό bit του χρήστη -> µJ/bit
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

%% ------------------ Κατάσταση καναλιού για την επόμενη κλήση ------------------
% Ό,τι χρειάζεται η επόμενη κλήση (αν είναι continuation, π.χ. επόμενο
% χρονικό βήμα του temporalPassSimulation.m) ώστε να υπολογίσει τη χωρική
% συσχέτιση του shadow fading - βλ. correlatedLosState.
newChannelState.UserGeo         = user_geo;
newChannelState.IsLOS           = losMat;
newChannelState.ShadowFading_dB = sfMat;

end

function [isLos, rho, useCorrelatedSF] = correlatedLosState(pLos, moveDistance, scenario, hasPrevState, prevIsLos)
% Υπολογίζει τη συσχετισμένη κατάσταση LOS/NLOS μιας ζεύξης BS-χρήστη,
% αντί να την επαναδειγματίζει ανεξάρτητα σε κάθε κλήση.
%
% Χωρική συσχέτιση shadow fading κατά Gudmundson (1991, "Correlation model
% for shadow fading in mobile radio systems", Electronics Letters 27,
% 2145-2146): εκθετική αυτοσυσχέτιση ρ(Δd) = exp(-Δd/d_corr), όπου d_corr η
% "correlation distance" στο οριζόντιο επίπεδο. Οι τιμές του d_corr για το
% shadow fading (SF) λαμβάνονται από το 3GPP TR 38.901 v16.1.0, Πίνακας
% 7.5-6 Part-1 ("Correlation distance in the horizontal plane [m]", σειρά
% SF): UMa LOS=37m, UMa NLOS=50m, UMi-Street Canyon LOS=10m, NLOS=13m.
%
% Το TR 38.901 δεν ορίζει ξεχωριστή "correlation distance" για την ίδια την
% κατηγορική κατάσταση LOS/NLOS (μόνο για τις LSP παραμέτρους όπως SF/K/DS/
% κ.λπ. στον Πίνακα 7.5-6) - ως απλοποίηση, εδώ η ίδια απόσταση συσχέτισης
% (και το ίδιο ρ) χρησιμοποιείται και ως πιθανότητα διατήρησης της
% προηγούμενης κατάστασης LOS/NLOS (Bernoulli, με πιθανότητα ρ διατηρείται,
% με πιθανότητα 1-ρ επαναδειγματίζεται από το pLos). Αυτό είναι συνεπές με
% τη λογική "drop-based" παραγωγής LSP του TR 38.901 §7.5 (Βήματα 2 και 4:
% η κατάσταση LOS/NLOS και το SF παράγονται μαζί, ανά "drop"), και ανάγεται
% ορθά στα δύο ακραία σενάρια: χωρίς προηγούμενη κατάσταση (ρ=0) πάντα νέο
% δείγμα (i.i.d., όπως πριν)· με μηδενική μετακίνηση (ρ=1, π.χ. ακίνητος
% χρήστης στο temporalPassSimulation.m) πάντα διατήρηση της προηγούμενης
% κατάστασης.
if ~hasPrevState
    isLos = rand() < pLos;
    rho = 0;
    useCorrelatedSF = false;
    return;
end

if prevIsLos
    switch scenario
        case 'UMa'
            dCorr = 37;
        case 'UMi'
            dCorr = 10;
        otherwise
            dCorr = 37;
    end
else
    switch scenario
        case 'UMa'
            dCorr = 50;
        case 'UMi'
            dCorr = 13;
        otherwise
            dCorr = 50;
    end
end
rho = exp(-moveDistance / dCorr);

if rand() < rho
    isLos = prevIsLos;
else
    isLos = rand() < pLos;
end

% Το AR(1) δείγμα SF είναι έγκυρο μόνο αν η κατάσταση LOS/NLOS δεν άλλαξε -
% μια μετάβαση LOS<->NLOS αλλάζει το σ_SF (TR 38.901 Πίνακας 7.4.1-1) και
% ακυρώνει τη στατιστική βάση του προηγούμενου δείγματος.
useCorrelatedSF = (isLos == prevIsLos);
end

function pLos = losProbability38901(d2D, hUT, scenario)
% Πιθανότητα LOS για μία ζεύξη BS-χρήστη, βάσει 3GPP TR 38.901 v17.0.0,
% Πίνακας 7.4.2-1 (LOS probability). d2D σε μέτρα (οριζόντια απόσταση),
% hUT το ύψος του χρήστη σε μέτρα.
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
% Απόσβεση λόγω ατμοσφαιρικών αερίων (οξυγόνο + υδρατμοί) σε ζεύξη
% δορυφόρου-χρήστη, κατά 3GPP TR 38.811 v15.1.0 §6.6.4, εξ. (6.6-8):
%   PLA(ε,f) = A_zenith(f) / sin(ε),   ε >= 5° (όριο ισχύος της απλοποιημένης
%   μεθόδου Annex 2 της ITU-R P.676-12 - πάντα αληθές εδώ αφού το σενάριο
%   έχει MinElevationDeg=10°).
%
% Το A_zenith (ζενιθιακή απόσβεση) υπολογίζεται από τα ισοδύναμα ύψη
% οξυγόνου/υδρατμών (ITU-R P.676-12, Annex 2, εξ. 30-39):
%   A_zenith = γ_o·h_o + γ_w·h_w
% με "reference atmosphere" T=288.15K, p=1013.25hPa, ρ=7.5 g/m^3 (mean
% annual global reference atmosphere, ITU-R P.835) - το ίδιο baseline που
% ορίζει το TR 38.811 §6.6.4 για system-level προσομοιώσεις.
%
% Οι ειδικές αποσβέσεις γ_o (ξηρός αέρας) και γ_w (υδρατμοί) σε dB/km
% υπολογίζονται με την ενσωματωμένη συνάρτηση gaspl του MATLAB (Communications
% Toolbox, υλοποιεί το πλήρες line-by-line μοντέλο του Annex 1 της ITU-R
% P.676-13, πιο ακριβές από τη χειρωνακτική Annex 2 μέθοδο για το ίδιο
% βήμα)· η ξηρή/υγρή συνιστώσα διαχωρίζονται καλώντας τη με και χωρίς
% πυκνότητα υδρατμών.
%
% ΣΗΜΕΙΩΣΗ αξιοπιστίας πηγής: οι συντελεστές των Πινάκων 3/4 και οι
% εξισώσεις (30)-(37) επαληθεύτηκαν από το επίσημο κείμενο της ITU-R
% P.676-12. Ο διορθωτικός όρος σ_w της εξ. (38) ανασυντέθηκε από ένα
% τμηματικά κατεστραμμένο (OCR/εξαγωγή κειμένου) απόσπασμα του PDF - η
% συνεισφορά του στο h_w είναι όμως <1% (κυριαρχεί ο σταθερός όρος A_w),
% οπότε τυχόν μικρή ανακρίβεια εδώ δεν επηρεάζει ουσιωδώς το αποτέλεσμα.

TcRef    = 15;        % °C (= 288.15 K)
TKRef    = 288.15;    % K
pPaRef   = 101325;    % Pa
pHpaRef  = 1013.25;   % hPa
rhoRef   = 7.5;       % g/m^3 (υδρατμοί, mean annual global reference atmosphere)

freqGHz = freqHz / 1e9;

gammaDry = gaspl(1000, freqHz, TcRef, pPaRef, 0);       % dB/km, μόνο οξυγόνο
gammaTot = gaspl(1000, freqHz, TcRef, pPaRef, rhoRef);  % dB/km, οξυγόνο+υδρατμοί
gammaWet = gammaTot - gammaDry;

[ho, hw] = equivalentHeightsP676(freqGHz, TKRef, pHpaRef, rhoRef);

Azenith = gammaDry*ho + gammaWet*hw;   % dB, εξ. (39)

pla_dB = Azenith / sind(elevDeg);      % dB, εξ. (6.6-8)
end

function [ho, hw] = equivalentHeightsP676(freqGHz, T_K, p_hPa, rho)
% Ισοδύναμα ύψη οξυγόνου/υδρατμών, ITU-R P.676-12 Annex 2, εξ. (30)-(38).
e_hPa = rho * T_K / 216.7;
rp = (p_hPa + e_hPa) / 1013.25;

% -- Οξυγόνο (εξ. 30-35a, Πίνακας 3) --
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
    ho = min(ho, 10.7*rp^0.3);   % εξ. (35a)
end

% -- Υδρατμοί (εξ. 35b-38, Πίνακας 4) --
A_w = 1.9298 - 0.04166*(T_K - 273.15) + 0.0517*e_hPa;
B_w = 1.1674 - 0.00622*(T_K - 273.15) + 0.0063*e_hPa;
sigma_w = 1 + 1.013 / (1 + exp(-8.6*(rp - 0.57)));   % εξ. (38), βλ. σημείωση πηγής

fi4 = [22.235080 183.310087 325.152888 380.197353 439.150807 448.001085 ...
       474.689092 488.490108 556.935985 620.700870 752.033113 916.171582 ...
       970.315022 987.926764];
ai4 = [1.52 7.62 1.56 4.15 0.20 1.63 0.76 0.26 7.81 1.25 16.2 1.47 1.36 1.60];
bi4 = [2.56 10.2 2.70 5.70 0.91 2.46 2.22 2.49 10.0 2.35 20.0 2.58 2.44 1.86];

hw = A_w + B_w * sum( (ai4*sigma_w) ./ ((freqGHz - fi4).^2 + bi4) );
end
