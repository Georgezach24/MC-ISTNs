function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters)
% Υπολογίζει, για κάθε χρήστη, τον καλύτερο κόμβο εξυπηρέτησης (BS ή δορυφόρο)
% βάσει SNR και την επιτευχθείσα χωρητικότητα Shannon μετά την κατανομή
% εύρους ζώνης. Εξάγει το βασικό μονοπάτι υπολογισμού από το test_simulation.m
% ώστε να μπορεί να κληθεί επανειλημμένα (π.χ. από έναν Monte-Carlo driver).

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

% Διαγνωστικοί πίνακες
groundDistanceMat = nan(numUsers,numBs);
range3DMat        = nan(numUsers,numBs);
pathLossMat       = nan(numUsers,numBs);
snrDbMat          = nan(numUsers,numBs);
pLosMat           = nan(numUsers,numBs);
losMat            = false(numUsers,numBs);
satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSnrDbVec       = nan(numUsers,1);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SNR) ------------------
for u = 1:numUsers
    % Αρχικοποιήσεις
    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

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
        % Πίνακας 7.4.2-1), αντί για μία σταθερή global τιμή LOS.
        pLos  = losProbability38901(groundDistance, user_geo(u,3), simParameters.PathLoss.Scenario);
        isLos = rand() < pLos;
        pLosMat(u,b) = pLos;
        losMat(u,b)  = isLos;

        [pathLoss, sigmaSF] = nrPathLoss(simParameters.PathLoss, ...
                              simParameters.CarrierFrequency, ...
                              isLos, ...
                              txPosition, rxPosition);

        % Shadow fading: log-normal δείγμα με τυπική απόκλιση sigmaSF (TR 38.901 §7.4.1)
        pathLoss = pathLoss + sigmaSF * randn();
        pathLossMat(u,b) = pathLoss;

        snr_db = (simParameters.TxPower - 30) - pathLoss - noisePowerBS_dBW;
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

    %% ===== Satellite candidate =====
    [azSat, elevSat, slantRangeSat] = geodetic2aer( ...
        sat_geo(1), sat_geo(2), sat_geo(3), ...
        user_geo(u,1), user_geo(u,2), user_geo(u,3), wgs84);

    satSlantRangeVec(u) = slantRangeSat;
    satElevationVec(u)  = elevSat;

    if elevSat >= satParameters.MinElevationDeg
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency;
        satPathLoss = fspl(slantRangeSat, lambdaSat);
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

    % Αποθήκευση επιλογής κόμβου για τον χρήστη
    bestNodeVec(u)         = userBestNode;
    bestNodeTypeVec(u)     = userBestType;
    bestDistanceVec(u)     = userBestDistance;
    bestPathLossVec(u)     = userBestPathLoss;
    bestSnrDbVec(u)        = userBestSNR;
    bestElevationDegVec(u) = userBestElevation;
end

%% ------------------ Υπολογισμός Χωρητικότητας & Ενέργειας (Κατανομή Πόρων) ------------------
for u = 1:numUsers
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

    % Υπολογισμός τελικής χωρητικότητας βάσει Shannon για το κομμάτι του B_user
    snr_lin = 10^(bestSnrDbVec(u)/10);
    capacity = B_user * log2(1 + snr_lin);   % bits/s

    capacityMbpsVec(u) = capacity * 1e-6;    % Mbps

    % Ενεργειακό proxy: ισομερής κατανομή ισχύος κόμβου ανά χρήστη (ίδια λογική
    % με το bandwidth split), διαιρεμένη με τον ρυθμό bit του χρήστη -> µJ/bit
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

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
