function datasetTable = monteCarloDriver(numScenarios, outputCsvPath)
%MONTECARLODRIVER Παράγει labeled dataset τρέχοντας το simulateScenario
% πάνω σε τυχαιοποιημένες τοπολογίες (θέσεις/πλήθος BS/χρηστών, γεωμετρία
% δορυφόρου, UMa/UMi). Ραδιο-configuration σταθερό, ίδιο με test_simulation.m.
%
% Χρήση:
%   monteCarloDriver()                          % 200 σενάρια -> ../Dataset/dataset.csv
%   monteCarloDriver(500)                       % 500 σενάρια -> ../Dataset/dataset.csv
%   T = monteCarloDriver(500, 'C:\out\ds.csv')  % custom πλήθος + διαδρομή

if nargin < 1 || isempty(numScenarios)
    numScenarios = 200;
end
if nargin < 2 || isempty(outputCsvPath)
    outputCsvPath = fullfile(fileparts(mfilename('fullpath')), '..', 'Dataset', 'dataset.csv');
end

rng(42); % Σταθερό RNG seed για αναπαραγωγιμότητα ολόκληρου του batch

wgs84 = wgs84Ellipsoid;

%% ------------------ Κέντρο περιοχής & εύρη τυχαιοποίησης τοπολογίας ------------------
baseLat = 37.9838;   % περιοχή Αθήνας, ίδιο σημείο με test_simulation.m
baseLon = 23.7275;

numBsRange          = [1 4];      % πλήθος BS ανά σενάριο
numUsersRange       = [3 10];     % πλήθος χρηστών ανά σενάριο
bsClusterRadiusKm   = 5;          % ακτίνα τοποθέτησης BS γύρω από το κέντρο
% Ακτίνα "κοντινών" χρηστών γύρω από τυχαίο BS: εύρος (όχι σταθερή τιμή)
% ώστε να καλύπτεται και η οριακή ζώνη γύρω από το usable-SNR κατώφλι
% (~2.75km σε UMa NLOS), όχι μόνο σίγουρα εντός/εκτός εμβέλειας -
% διαφορετικά ένα μοντέλο ML με πρόσβαση μόνο σε γεωμετρία λύνει σχεδόν
% τέλεια το terrestrial σκέλος του ServingType. Δειγματοληψία ομοιόμορφη
% στην ίδια την ακτίνα (randRadiusOffsetDeg), όχι στο εμβαδόν δίσκου.
nearUserRadiusKmRange = [0.1 4.0];
farUserRadiusKmRange = [50 150];  % "μακρινοί" χρήστες, πρακτικά εκτός εμβέλειας BS
farUserProbability  = 0.3;
% Τυχαιοποίηση θέσης υποδορυφορικού σημείου. Footprint ακτίνας ορατότητας
% ~1.6-1.75 Mm στα 500-600km/10° μάσκα (γ=arccos(Re·cosε/(Re+h))-ε, d=Re·γ).
% ±20° (~1.75-2.2 Mm) ώστε ο δορυφόρος να ΜΗΝ είναι πάντα ορατός σε κάθε
% σενάριο (ένα LEO δεν είναι γεωστατικό) - παράγει γνήσιες Terrestrial-only
% περιπτώσεις αντί να είναι δομικά αδύνατες.
satLatJitterDeg     = 20;
satLonJitterDeg     = 20;
satAltitudeRangeM   = [500e3 600e3]; % τυπικό εύρος υψομέτρου LEO

%% ------------------ Σταθερό ραδιο-configuration (ίδιο με test_simulation.m) ------------------
simParametersBase.Carrier = nrCarrierConfig;
simParametersBase.Carrier.NSizeGrid = 51;
simParametersBase.Carrier.SubcarrierSpacing = 30;
simParametersBase.Carrier.CyclicPrefix = 'Normal';
simParametersBase.CarrierFrequency = 3.5e9;
simParametersBase.TxPower = 43;
simParametersBase.AntennaGain = 8;   % dBi, TR 38.901 §7.3, Πίνακας 7.3-1
simParametersBase.EIRP = simParametersBase.TxPower + simParametersBase.AntennaGain;
simParametersBase.RxNoiseFigure = 5;
simParametersBase.RxAntTemperature = 290;

simParametersBase.PathLossModel = '5G-NR';
simParametersBase.PathLoss = nrPathLossConfig;

simParametersBase.Power.NumTrx = 1;
simParametersBase.Power.P0     = 130;
simParametersBase.Power.DeltaP = 4.7;
simParametersBase.Power.Psleep = 75;

satParametersBase.CarrierFrequency = 2.01e9;
satParametersBase.TxPower = 34;
satParametersBase.AntennaGain = 30;
satParametersBase.EIRP = satParametersBase.TxPower + satParametersBase.AntennaGain;
satParametersBase.Bandwidth = 20e6;
satParametersBase.MinElevationDeg = 10;
satParametersBase.Power.Pfix  = 0;
satParametersBase.Power.EtaPA = 0.4;

%% ------------------ Κύριος βρόχος Monte-Carlo ------------------
allScenarioTables = cell(numScenarios,1);

for s = 1:numScenarios
    numBs    = randi(numBsRange);
    numUsers = randi(numUsersRange);

    if rand() < 0.5
        scenarioType = 'UMa';
        bs_height_m = 25;
    else
        scenarioType = 'UMi';
        bs_height_m = 10;
    end
    simParametersBase.PathLoss.Scenario = scenarioType;
    simParametersBase.PathLoss.EnvironmentHeight = 1;

    % -- Τοποθέτηση BS γύρω από το κέντρο --
    bs_geo = zeros(numBs,3);
    for b = 1:numBs
        [dLat, dLon] = randOffsetDeg(baseLat, bsClusterRadiusKm);
        bs_geo(b,:) = [baseLat + dLat, baseLon + dLon, bs_height_m];
    end

    % -- Τοποθέτηση χρηστών: κοντά σε τυχαίο BS, ή μακριά (μόνο δορυφόρος) --
    user_geo = zeros(numUsers,3);
    for u = 1:numUsers
        if rand() < farUserProbability
            radiusKm = farUserRadiusKmRange(1) + diff(farUserRadiusKmRange)*rand();
            [dLat, dLon] = randOffsetDeg(baseLat, radiusKm);
            user_geo(u,:) = [baseLat + dLat, baseLon + dLon, 1.5];
        else
            refB = randi(numBs);
            nearRadiusKm = nearUserRadiusKmRange(1) + diff(nearUserRadiusKmRange)*rand();
            [dLat, dLon] = randRadiusOffsetDeg(bs_geo(refB,1), nearRadiusKm);
            user_geo(u,:) = [bs_geo(refB,1) + dLat, bs_geo(refB,2) + dLon, 1.5];
        end
    end

    % -- Τυχαιοποίηση γεωμετρίας δορυφόρου --
    sat_geo = [baseLat + (2*rand()-1)*satLatJitterDeg, ...
               baseLon + (2*rand()-1)*satLonJitterDeg, ...
               satAltitudeRangeM(1) + diff(satAltitudeRangeM)*rand()];

    % -- Εκτέλεση σεναρίου --
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVecAll, satPathLossVec, satSnrDbVecAll] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParametersBase, satParametersBase);

    % Φορτίο κόμβου: πλήθος χρηστών του ίδιου σεναρίου στον ίδιο κόμβο.
    % DualConnectivity φορτίζει BS ΚΑΙ δορυφόρο ταυτόχρονα -> δύο στήλες.
    isBsConnected  = (bestNodeTypeVec == "Terrestrial") | (bestNodeTypeVec == "DualConnectivity");
    isSatConnected = (bestNodeTypeVec == "Satellite")   | (bestNodeTypeVec == "DualConnectivity");
    bsIdVec = extractBefore(bestNodeVec + "+", "+");   % "BSx" για Terrestrial/DualConnectivity
    bsLoadVec  = zeros(numUsers,1);
    satLoadVec = zeros(numUsers,1);
    for u = 1:numUsers
        if isBsConnected(u)
            bsLoadVec(u) = sum(isBsConnected & bsIdVec == bsIdVec(u));
        end
        if isSatConnected(u)
            satLoadVec(u) = sum(isSatConnected);
        end
    end

    scenarioID      = repmat(s, numUsers, 1);
    userID          = (1:numUsers)';
    scenarioTypeCol = repmat(string(scenarioType), numUsers, 1);
    numBsCol        = repmat(numBs, numUsers, 1);
    numUsersCol     = repmat(numUsers, numUsers, 1);
    userLat         = user_geo(:,1);
    userLon         = user_geo(:,2);

    allScenarioTables{s} = table(scenarioID, userID, scenarioTypeCol, numBsCol, numUsersCol, ...
        userLat, userLon, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, ...
        energyPerBitUJVec, bsLoadVec, satLoadVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSnrDbVecAll, satElevationVecAll, satSlantRangeVec, satPathLossVec, ...
        'VariableNames', {'ScenarioID','UserID','ScenarioType','NumBS','NumUsers', ...
        'UserLat','UserLon','ServingNode','ServingType','Distance_m','PathLoss_dB', ...
        'SNR_dB','Capacity_Mbps','SatElevation_deg','NodePower_W','EnergyPerBit_uJ', ...
        'BsLoad','SatLoad', ...
        'CandBS_SNR_dB','CandBS_Distance_m','CandBS_PathLoss_dB', ...
        'CandSat_SNR_dB','CandSat_Elevation_deg','CandSat_SlantRange_m','CandSat_PathLoss_dB'});
end

datasetTable = vertcat(allScenarioTables{:});

%% ------------------ Εγγραφή CSV ------------------
outDir = fileparts(outputCsvPath);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
writetable(datasetTable, outputCsvPath);

numTerrestrial = sum(datasetTable.ServingType == "Terrestrial");
numSatellite   = sum(datasetTable.ServingType == "Satellite");
numDual        = sum(datasetTable.ServingType == "DualConnectivity");
numOutage      = sum(datasetTable.ServingType == "Outage");
fprintf('Monte-Carlo dataset: %d σενάρια, %d γραμμές χρηστών (%d Terrestrial, %d Satellite, %d DualConnectivity, %d Outage) -> %s\n', ...
    numScenarios, height(datasetTable), numTerrestrial, numSatellite, numDual, numOutage, outputCsvPath);

end

function [dLat, dLon] = randOffsetDeg(refLat, radiusKm)
% Τυχαία μετατόπιση [dLat, dLon] σε μοίρες, ομοιόμορφη εντός δίσκου
% ακτίνας radiusKm γύρω από σημείο πλάτους refLat (επίπεδη προσέγγιση).
r = radiusKm * sqrt(rand());
theta = 2*pi*rand();
dNorthKm = r*cos(theta);
dEastKm  = r*sin(theta);
dLat = dNorthKm / 110.574;
dLon = dEastKm / (111.320*cosd(refLat));
end

function [dLat, dLon] = randRadiusOffsetDeg(refLat, radiusKm)
% Ίδιο με randOffsetDeg, αλλά σε ΑΚΡΙΒΩΣ απόσταση radiusKm (τυχαία γωνία
% μόνο) - ελέγχει ρητά την κατανομή αποστάσεων αντί να τη συγκεντρώνει
% προς το άνω άκρο όπως η ομοιόμορφη-σε-δίσκο δειγματοληψία.
theta = 2*pi*rand();
dNorthKm = radiusKm*cos(theta);
dEastKm  = radiusKm*sin(theta);
dLat = dNorthKm / 110.574;
dLon = dEastKm / (111.320*cosd(refLat));
end
