function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
    newChannelState] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, prevChannelState)
% Υπολογίζει τον κόμβο εξυπηρέτησης (BS ή/και δορυφόρος) και τη
% χωρητικότητα κάθε χρήστη. Καλείται επανειλημμένα από test_simulation.m,
% monteCarloDriver.m κ.λπ. Επιστρέφει και per-candidate διαγνωστικά
% (καλύτερο BS, δορυφόρος) για το ML pipeline (Part 2).
%
% Πολυσυνδεσιμότητα (SS-SBS, Li & Shang, Κεφ.2): ο χρήστης συνδέεται σε
% κάθε κόμβο (καλύτερο BS, δορυφόρος) που ξεπερνά το minUsableSnrDb,
% ταυτόχρονα σε όσους το ξεπερνούν ("DualConnectivity"), αλλιώς "Outage".
% Χωρητικότητα/ενέργεια DualConnectivity = άθροισμα των ενεργών ζεύξεων
% (carrier aggregation).
%
% prevChannelState (προαιρετικό): αν δοθεί, το shadow fading/LOS κάθε
% ζεύξης συσχετίζεται χωρικά με την προηγούμενη κλήση (βλ. correlatedLosState)
% αντί να επαναδειγματίζεται ανεξάρτητα. Χρησιμοποιείται μόνο από callers
% διαδοχικών μεταδόσεων της ίδιας τοπολογίας (temporalPassSimulation.m).
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

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SNR (outage) ------------------
% Shannon-ισοδύναμο SNR για MCS 0 (3GPP TS 38.214 Πίνακας 5.1.3.1-2, QPSK,
% code rate 120/1024 -> 0.2344 bits/s/Hz). Κάτω από αυτό, outage.
minSpectralEfficiency = 0.2344;                        % bits/s/Hz, TS 38.214 §5.1.3.1, MCS 0
minUsableSnrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Υστέρηση (hysteresis) στην ενεργοποίηση/απενεργοποίηση ζεύξης ------------------
% Event A3-style (3GPP TS 38.331): μια ζεύξη ενεργοποιείται μόνο αφού
% SNR >= minUsableSnrDb+MarginDb για TimeToTriggerSteps διαδοχικές κλήσεις,
% απενεργοποιείται μόνο μετά από εξίσου επίμονη πτώση κάτω από
% minUsableSnrDb-MarginDb - νεκρή ζώνη γύρω από το κατώφλι, αποτρέπει
% ping-pong. Προαιρετικό: simParameters.Hysteresis μη ορισμένο ->
% MarginDb=0, TimeToTriggerSteps=0 (άμεση συμπεριφορά κατωφλίου).
if isfield(simParameters, 'Hysteresis') && isfield(simParameters.Hysteresis, 'MarginDb')
    hystMarginDb = simParameters.Hysteresis.MarginDb;
else
    hystMarginDb = 0;
end
if isfield(simParameters, 'Hysteresis') && isfield(simParameters.Hysteresis, 'TimeToTriggerSteps')
    hystTtt = simParameters.Hysteresis.TimeToTriggerSteps;
else
    hystTtt = 0;
end

hasPrevActivation = ~isempty(prevChannelState) && isfield(prevChannelState, 'ActiveBs') && ...
    isequal(size(prevChannelState.ActiveBs), [numUsers, 1]);
if hasPrevActivation
    prevActiveBs         = prevChannelState.ActiveBs;
    prevActiveSat         = prevChannelState.ActiveSat;
    prevBsPendingCounter  = prevChannelState.BsPendingCounter;
    prevSatPendingCounter = prevChannelState.SatPendingCounter;
else
    prevActiveBs          = false(numUsers,1);
    prevActiveSat          = false(numUsers,1);
    prevBsPendingCounter   = zeros(numUsers,1);
    prevSatPendingCounter  = zeros(numUsers,1);
end
newActiveBs         = false(numUsers,1);
newActiveSat         = false(numUsers,1);
newBsPendingCounter  = zeros(numUsers,1);
newSatPendingCounter = zeros(numUsers,1);

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

% Φορτίο ανά τύπο κόμβου: DualConnectivity = connected σε BS ΚΑΙ δορυφόρο.
connectedBsIdVec = strings(numUsers,1);   % "" αν όχι BS, αλλιώς "BSx"
connectedSatVec  = false(numUsers,1);     % true αν συνδεδεμένος στον δορυφόρο

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
    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;

    % Απόσταση μετακίνησης από την προηγούμενη κλήση (0 = ακίνητος χρήστης).
    if hasPrevState
        userMoveDistance = distance(prevChannelState.UserGeo(u,1), prevChannelState.UserGeo(u,2), ...
                                     user_geo(u,1), user_geo(u,2), wgs84);
    else
        userMoveDistance = Inf;
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

        % LOS ανά ζεύξη (TR 38.901 §7.4.2, Πίνακας 7.4.2-1), συσχετισμένο
        % χωρικά με prevChannelState αν υπάρχει - βλ. correlatedLosState.
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

        % Shadow fading: log-normal, TR 38.901 §7.4.1. Συσχετισμένο κατά
        % Gudmundson (1991) αν useCorrelatedSF, αλλιώς ανεξάρτητο δείγμα.
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
            userBestDistance  = d3d;
            userBestPathLoss  = pathLoss;
        end
    end

    % Καλύτερος υποψήφιος BS πριν συγκριθεί με τον δορυφόρο - per-candidate
    % διαγνωστικό (ML feature), ανεξάρτητο από τον τελικό νικητή.
    bestBsSnrDbVec(u)   = userBestSNR;
    bestBsDistanceVec(u) = userBestDistance;
    bestBsPathLossVec(u) = userBestPathLoss;
    bsWinnerNode          = userBestNode;

    %% ===== Satellite candidate =====
    [~, elevSat, slantRangeSat] = geodetic2aer( ...
        sat_geo(1), sat_geo(2), sat_geo(3), ...
        user_geo(u,1), user_geo(u,2), user_geo(u,3), wgs84);

    satSlantRangeVec(u) = slantRangeSat;
    satElevationVec(u)  = elevSat;

    if elevSat >= satParameters.MinElevationDeg
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency;
        satPathLoss = fspl(slantRangeSat, lambdaSat);

        % Ατμοσφαιρική απόσβεση αερίων, TR 38.811 §6.6.4 εξ. (6.6-8) - βλ.
        % gasAttenuationSlantP676. Βροχή/νέφωση παραλείπεται (§6.6.5,
        % αμελητέα <6GHz). Ιονοσφαιρική σπινθηρίδα (§6.6.6) εκτός πεδίου
        % (χρειάζεται κλιματολογικό μοντέλο) - future work.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        satSnrDb = (satParameters.EIRP - 30) - satPathLoss - noisePowerSAT_dBW;
    else
        satPathLoss = inf;
        satSnrDb = -Inf;
    end

    satPathLossVec(u) = satPathLoss;
    satSnrDbVec(u)    = satSnrDb;

    % ===== Απόφαση συνδεσιμότητας (SS-SBS, τέσσερις καταστάσεις) =====
    % Καλύτερος BS και δορυφόρος αξιολογούνται ανεξάρτητα έναντι του
    % minUsableSnrDb: DualConnectivity αν περνούν και οι δύο, Terrestrial/
    % Satellite αν μόνο ο ένας, Outage αν κανένας. Η ενεργοποίηση περνάει
    % από τη μηχανή υστέρησης (updateLinkActivation) εκτός από την πρώτη
    % κλήση μιας ακολουθίας (hasPrevActivation=false), όπου αποφασίζεται
    % άμεσα (η υστέρηση αφορά μεταβάσεις, όχι αρχική απόκτηση).
    if hasPrevActivation
        [bsUsable, newBsPendingCounter(u)] = updateLinkActivation( ...
            prevActiveBs(u), userBestSNR, minUsableSnrDb, hystMarginDb, hystTtt, prevBsPendingCounter(u));
        [satUsable, newSatPendingCounter(u)] = updateLinkActivation( ...
            prevActiveSat(u), satSnrDb, minUsableSnrDb, hystMarginDb, hystTtt, prevSatPendingCounter(u));
    else
        bsUsable  = userBestSNR >= minUsableSnrDb;
        satUsable = satSnrDb    >= minUsableSnrDb;
        newBsPendingCounter(u)  = 0;
        newSatPendingCounter(u) = 0;
    end
    newActiveBs(u)  = bsUsable;
    newActiveSat(u) = satUsable;

    if bsUsable && satUsable
        bestNodeTypeVec(u)  = "DualConnectivity";
        bestNodeVec(u)      = bsWinnerNode + "+SAT-1";
        connectedBsIdVec(u) = bsWinnerNode;
        connectedSatVec(u)  = true;
        % Κανένα ενιαίο SNR/distance/path loss νικητή σε DualConnectivity -
        % τα per-candidate διαγνωστικά παραμένουν διαθέσιμα ξεχωριστά.
        bestSnrDbVec(u)        = NaN;
        bestDistanceVec(u)     = NaN;
        bestPathLossVec(u)     = NaN;
        bestElevationDegVec(u) = elevSat;
    elseif bsUsable
        bestNodeTypeVec(u)  = "Terrestrial";
        bestNodeVec(u)      = bsWinnerNode;
        connectedBsIdVec(u) = bsWinnerNode;
        connectedSatVec(u)  = false;
        bestSnrDbVec(u)        = userBestSNR;
        bestDistanceVec(u)     = userBestDistance;
        bestPathLossVec(u)     = userBestPathLoss;
        bestElevationDegVec(u) = NaN;
    elseif satUsable
        bestNodeTypeVec(u)  = "Satellite";
        bestNodeVec(u)      = "SAT-1";
        connectedBsIdVec(u) = "";
        connectedSatVec(u)  = true;
        bestSnrDbVec(u)        = satSnrDb;
        bestDistanceVec(u)     = slantRangeSat;
        bestPathLossVec(u)     = satPathLoss;
        bestElevationDegVec(u) = elevSat;
    else
        % Outage - διαγνωστικά του λιγότερο κακού υποψηφίου διατηρούνται.
        bestNodeTypeVec(u)  = "Outage";
        bestNodeVec(u)      = "None";
        connectedBsIdVec(u) = "";
        connectedSatVec(u)  = false;
        if satSnrDb > userBestSNR
            bestSnrDbVec(u)        = satSnrDb;
            bestDistanceVec(u)     = slantRangeSat;
            bestPathLossVec(u)     = satPathLoss;
            bestElevationDegVec(u) = elevSat;
        else
            bestSnrDbVec(u)        = userBestSNR;
            bestDistanceVec(u)     = userBestDistance;
            bestPathLossVec(u)     = userBestPathLoss;
            bestElevationDegVec(u) = NaN;
        end
    end
end

% bsUsableMat: ποιοι BS είναι χρησιμοποιήσιμοι ανά χρήστη (όχι μόνο ο
% καλύτερος) - χρειάζεται στο joint load-balancing pass παρακάτω.
bsUsableMat = snrDbMat >= minUsableSnrDb;

% SNR της ζεύξης που ΠΡΑΓΜΑΤΙΚΑ εξυπηρετεί κάθε χρήστη. Ξεχωριστό από
% bestBsSnrDbVec (παραμένει αμετάβλητο ML diagnostic) γιατί το joint pass
% μπορεί να μετακινήσει έναν χρήστη σε διαφορετικό BS απ' ό,τι ήταν αρχικά
% ο καλύτερος.
servingBsSnrDbVec = bestBsSnrDbVec;

%% ------------------ Joint/fairness-aware εξισορρόπηση φορτίου (προαιρετική) ------------------
% Μέχρι εδώ κάθε χρήστης αποφάσισε ανεξάρτητα (greedy). Εδώ προστίθεται
% προαιρετικά μια load-aware φάση (simParameters.Fairness.MaxUsersPerBs,
% default Inf = off), σε δύο αλγορίθμους:
%
% - "PerBs" (simParameters.Fairness.Joint = false/μη ορισμένο, default):
%   διατρέχει τους BS με σταθερή σειρά· αν ένας υπερβαίνει το όριο, οι ήδη
%   DualConnectivity χρήστες του αποσυνδέονται προτεραιοποιημένα
%   (χαμηλότερο BS SNR πρώτα) μέχρι να επανέλθει εντός ορίου. Μοναδική
%   επιλογή: πτώση σε δορυφόρο-μόνο.
%
% - "Joint" (simParameters.Fairness.Joint = true): σε κάθε επανάληψη
%   βρίσκει τον πιο υπερφορτωμένο BS σε όλο το δίκτυο, και προτιμά lateral
%   handover σε άλλον χρησιμοποιήσιμο BS με ελεύθερη χωρητικότητα αντί για
%   πτώση σε δορυφόρο-μόνο· μπορεί έτσι να ωφελήσει και Terrestrial-only
%   χρήστες (όχι μόνο να τους προστατέψει, όπως το PerBs). Greedy heuristic
%   (worst-BS-first, worst-user-SNR-first), όχι πλήρης βελτιστοποίηση -
%   βλ. CLAUDE.md/Κεφ.6.
%
% Δεν αγγίζει την κατάσταση υστέρησης (newActiveBs κ.λπ.) - είναι
% διοικητική απόφαση αποδοχής πάνω σε ήδη SNR-επιλέξιμες ζεύξεις.
if isfield(simParameters, 'Fairness') && isfield(simParameters.Fairness, 'MaxUsersPerBs')
    maxUsersPerBs = simParameters.Fairness.MaxUsersPerBs;
else
    maxUsersPerBs = Inf;
end
if isfield(simParameters, 'Fairness') && isfield(simParameters.Fairness, 'Joint')
    jointFairness = simParameters.Fairness.Joint;
else
    jointFairness = false;
end

if isfinite(maxUsersPerBs) && ~jointFairness
    for b = 1:numBs
        bsName  = "BS" + string(b);
        bsUsers = find(connectedBsIdVec == bsName);
        overload = numel(bsUsers) - maxUsersPerBs;
        if overload <= 0
            continue;
        end

        eligible = bsUsers(bestNodeTypeVec(bsUsers) == "DualConnectivity");
        [~, order] = sort(bestBsSnrDbVec(eligible), 'ascend');
        eligible = eligible(order);
        numToOffload = min(overload, numel(eligible));

        for k = 1:numToOffload
            u = eligible(k);
            % Αποσύνδεση μόνο του BS - ο δορυφόρος παραμένει ενεργός.
            bestNodeTypeVec(u)     = "Satellite";
            bestNodeVec(u)         = "SAT-1";
            connectedBsIdVec(u)    = "";
            bestSnrDbVec(u)        = satSnrDbVec(u);
            bestDistanceVec(u)     = satSlantRangeVec(u);
            bestPathLossVec(u)     = satPathLossVec(u);
            bestElevationDegVec(u) = satElevationVec(u);
        end
    end
elseif isfinite(maxUsersPerBs) && jointFairness
    for iter = 1:numUsers   % κάθε επιτυχής επανάληψη μετακινεί έναν χρήστη
        bsLoadNow = zeros(numBs,1);
        for b = 1:numBs
            bsLoadNow(b) = sum(connectedBsIdVec == ("BS" + string(b)));
        end
        overloadNow = bsLoadNow - maxUsersPerBs;
        [worstOverload, bWorst] = max(overloadNow);
        if worstOverload <= 0
            break;
        end

        bsWorstUsers = find(connectedBsIdVec == ("BS" + string(bWorst)));
        bestAltBsForUser = zeros(numel(bsWorstUsers),1);   % 0 = καμία εναλλακτική BS
        for k = 1:numel(bsWorstUsers)
            u = bsWorstUsers(k);
            altBsIdx = find(bsUsableMat(u,:) & (1:numBs) ~= bWorst & bsLoadNow' < maxUsersPerBs);
            if ~isempty(altBsIdx)
                [~, pick] = min(bsLoadNow(altBsIdx));   % BS με τη μεγαλύτερη ελεύθερη χωρητικότητα
                bestAltBsForUser(k) = altBsIdx(pick);
            end
        end
        canDropToSat = connectedSatVec(bsWorstUsers);

        movable = find(bestAltBsForUser > 0 | canDropToSat);
        if isempty(movable)
            break;
        end

        [~, order] = sort(snrDbMat(bsWorstUsers(movable), bWorst), 'ascend');
        k = movable(order(1));
        u = bsWorstUsers(k);

        if bestAltBsForUser(k) > 0
            % Lateral handover σε λιγότερο φορτωμένο BS.
            bNew = bestAltBsForUser(k);
            connectedBsIdVec(u)  = "BS" + string(bNew);
            servingBsSnrDbVec(u) = snrDbMat(u, bNew);
            if connectedSatVec(u)
                bestNodeVec(u) = "BS" + string(bNew) + "+SAT-1";
            else
                bestNodeVec(u)     = "BS" + string(bNew);
                bestSnrDbVec(u)    = servingBsSnrDbVec(u);
                bestDistanceVec(u) = range3DMat(u, bNew);
                bestPathLossVec(u) = pathLossMat(u, bNew);
            end
        else
            % Καμία εναλλακτική BS -> πτώση σε δορυφόρο-μόνο.
            bestNodeTypeVec(u)     = "Satellite";
            bestNodeVec(u)         = "SAT-1";
            connectedBsIdVec(u)    = "";
            bestSnrDbVec(u)        = satSnrDbVec(u);
            bestDistanceVec(u)     = satSlantRangeVec(u);
            bestPathLossVec(u)     = satPathLossVec(u);
            bestElevationDegVec(u) = satElevationVec(u);
        end
    end
end

%% ------------------ Υπολογισμός Χωρητικότητας & Ενέργειας (Κατανομή Πόρων) ------------------
% Φορτίο ανά κόμβο - DualConnectivity φορτίζει BS ΚΑΙ δορυφόρο ταυτόχρονα.
bsLoadVec = zeros(numBs,1);
for b = 1:numBs
    bsLoadVec(b) = sum(connectedBsIdVec == ("BS" + string(b)));
end
satLoad = sum(connectedSatVec);

% Κατανάλωση ισχύος (EARTH model για BS, γραμμικό μοντέλο ενισχυτή για
% δορυφόρο - βλ. CLAUDE.md). Ίδια ανά χρήστη, υπολογίζεται μία φορά.
pOutW_bs  = 10^((simParameters.TxPower - 30)/10);
nodePowerW_bsActive = simParameters.Power.NumTrx * ...
    (simParameters.Power.P0 + simParameters.Power.DeltaP * pOutW_bs);

pOutW_sat = 10^((satParameters.TxPower - 30)/10);
nodePowerW_satActive = satParameters.Power.Pfix + pOutW_sat / satParameters.Power.EtaPA;

% Ανώτατη φασματική απόδοση (TS 38.214, MCS 27, 256QAM -> 5.5547 bits/s/Hz).
maxSpectralEfficiency = 5.5547;   % bits/s/Hz, TS 38.214 §5.1.3.1/§5.2.2.1

for u = 1:numUsers
    if bestNodeTypeVec(u) == "Outage"
        capacityMbpsVec(u)   = 0;
        nodePowerWattsVec(u) = 0;
        energyPerBitUJVec(u) = Inf;
        continue;
    end

    % Αθροιστική χωρητικότητα/ισχύς (carrier aggregation) στις ενεργές ζεύξεις.
    capacityBps = 0;
    powerW      = 0;

    if connectedBsIdVec(u) ~= ""
        bIdx   = str2double(extractAfter(connectedBsIdVec(u), "BS"));
        B_user = BW_bs / bsLoadVec(bIdx);
        snr_lin = 10^(servingBsSnrDbVec(u)/10);
        spectralEfficiency = min(log2(1 + snr_lin), maxSpectralEfficiency);
        capacityBps = capacityBps + B_user * spectralEfficiency;
        powerW      = powerW + nodePowerW_bsActive / bsLoadVec(bIdx);
    end

    if connectedSatVec(u)
        B_user = satParameters.Bandwidth / satLoad;
        snr_lin = 10^(satSnrDbVec(u)/10);
        spectralEfficiency = min(log2(1 + snr_lin), maxSpectralEfficiency);
        capacityBps = capacityBps + B_user * spectralEfficiency;
        powerW      = powerW + nodePowerW_satActive / satLoad;
    end

    capacityMbpsVec(u)   = capacityBps * 1e-6;   % Mbps
    nodePowerWattsVec(u) = powerW;
    energyPerBitUJVec(u) = powerW / capacityBps * 1e6;   % µJ/bit
end

%% ------------------ Κατάσταση καναλιού για την επόμενη κλήση ------------------
newChannelState.UserGeo         = user_geo;
newChannelState.IsLOS           = losMat;
newChannelState.ShadowFading_dB = sfMat;

newChannelState.ActiveBs          = newActiveBs;
newChannelState.ActiveSat          = newActiveSat;
newChannelState.BsPendingCounter   = newBsPendingCounter;
newChannelState.SatPendingCounter  = newSatPendingCounter;

end

function [isActive, pendingCounter] = updateLinkActivation(wasActive, snrDb, snrMinDb, marginDb, tttSteps, pendingCounter)
% Υστέρηση (hysteresis) + time-to-trigger (TTT) για ενεργοποίηση/
% απενεργοποίηση ζεύξης, κατά Event A3 (3GPP TS 38.331): ενεργοποίηση μόνο
% αφού SNR>=snrMinDb+marginDb για tttSteps+1 διαδοχικές κλήσεις,
% απενεργοποίηση μόνο μετά από εξίσου επίμονη πτώση κάτω από
% snrMinDb-marginDb. marginDb=0/tttSteps=0 -> άμεσο κατώφλι (backward-compatible).
if wasActive
    conditionForChange = snrDb < (snrMinDb - marginDb);
else
    conditionForChange = snrDb >= (snrMinDb + marginDb);
end

if conditionForChange
    pendingCounter = pendingCounter + 1;
else
    pendingCounter = 0;
end

if pendingCounter > tttSteps
    isActive = ~wasActive;
    pendingCounter = 0;
else
    isActive = wasActive;
end
end

function [isLos, rho, useCorrelatedSF] = correlatedLosState(pLos, moveDistance, scenario, hasPrevState, prevIsLos)
% Χωρικά συσχετισμένη κατάσταση LOS/NLOS, κατά Gudmundson (1991):
% εκθετική αυτοσυσχέτιση ρ(Δd)=exp(-Δd/d_corr), με d_corr από TR 38.901
% Πίνακα 7.5-6 (UMa LOS=37m/NLOS=50m, UMi LOS=10m/NLOS=13m). Το ίδιο ρ
% χρησιμοποιείται και ως πιθανότητα διατήρησης της προηγούμενης
% κατάστασης LOS/NLOS. ρ=0 (χωρίς προηγούμενη κατάσταση) -> ανεξάρτητο
% δείγμα· ρ=1 (μηδενική μετακίνηση) -> διατήρηση προηγούμενης κατάστασης.
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

% Το AR(1) δείγμα SF ισχύει μόνο αν η κατάσταση LOS/NLOS δεν άλλαξε.
useCorrelatedSF = (isLos == prevIsLos);
end

function pLos = losProbability38901(d2D, hUT, scenario)
% Πιθανότητα LOS, 3GPP TR 38.901 v17.0.0 Πίνακας 7.4.2-1. d2D σε μέτρα,
% hUT ύψος χρήστη σε μέτρα.
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
            'Άγνωστο PathLoss.Scenario "%s" - ορισμένο μόνο για "UMa"/"UMi".', ...
            scenario);
end
end

function pla_dB = gasAttenuationSlantP676(freqHz, elevDeg)
% Απόσβεση ατμοσφαιρικών αερίων (O2+υδρατμοί), 3GPP TR 38.811 v15.1.0
% §6.6.4 εξ. (6.6-8): PLA(ε,f) = A_zenith(f)/sin(ε), ε>=5°.
% A_zenith από ισοδύναμα ύψη οξυγόνου/υδρατμών (ITU-R P.676-12 Annex 2,
% εξ. 30-39), reference atmosphere T=288.15K/p=1013.25hPa/ρ=7.5g/m^3
% (ITU-R P.835), όπως ορίζει το TR 38.811 §6.6.4.
%
% Ειδικές αποσβέσεις γ_o/γ_w μέσω gaspl (Communications Toolbox, line-by-line
% μοντέλο ITU-R P.676-13 Annex 1).
%
% Σημείωση: ο διορθωτικός όρος σ_w της εξ. (38) ανασυντέθηκε από
% μερικώς κατεστραμμένο απόσπασμα πηγής· συνεισφορά <1% στο h_w, αμελητέο.

TcRef    = 15;        % °C (= 288.15 K)
TKRef    = 288.15;    % K
pPaRef   = 101325;    % Pa
pHpaRef  = 1013.25;   % hPa
rhoRef   = 7.5;       % g/m^3 (υδρατμοί)

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
sigma_w = 1 + 1.013 / (1 + exp(-8.6*(rp - 0.57)));   % εξ. (38)

fi4 = [22.235080 183.310087 325.152888 380.197353 439.150807 448.001085 ...
       474.689092 488.490108 556.935985 620.700870 752.033113 916.171582 ...
       970.315022 987.926764];
ai4 = [1.52 7.62 1.56 4.15 0.20 1.63 0.76 0.26 7.81 1.25 16.2 1.47 1.36 1.60];
bi4 = [2.56 10.2 2.70 5.70 0.91 2.46 2.22 2.49 10.0 2.35 20.0 2.58 2.44 1.86];

hw = A_w + B_w * sum( (ai4*sigma_w) ./ ((freqGHz - fi4).^2 + bi4) );
end
