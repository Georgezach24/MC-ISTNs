function datasetTable = monteCarloDriver(numScenarios, outputCsvPath, label)
%MONTECARLODRIVER Παράγει ένα labeled dataset τρέχοντας το simulateScenario
% πολλές φορές πάνω σε τυχαιοποιημένες τοπολογίες: θέσεις BS/χρηστών, πλήθος
% BS/χρηστών, υποδορυφορικό σημείο (υψόμετρο σταθερό στο LEO-600), και
% σενάριο TR 38.901 (UMa/UMi). Το ραδιο-configuration παραμένει σταθερό,
% ίδιο με το test_simulation.m.
%
% Χρήση:
%   monteCarloDriver()                               % 200 σενάρια -> ../Dataset/dataset.csv
%   monteCarloDriver(500)                            % 500 σενάρια
%   T = monteCarloDriver(500, 'C:\out\ds.csv')       % custom πλήθος + διαδρομή
%   T = monteCarloDriver(500, [], 'tag')             % tag στο όνομα του versioned φακέλου

if nargin < 1 || isempty(numScenarios)
    numScenarios = 200;
end
if nargin < 2 || isempty(outputCsvPath)
    outputCsvPath = fullfile(fileparts(mfilename('fullpath')), '..', 'Dataset', 'dataset.csv');
end
if nargin < 3
    label = '';
end

rng(42); % Σταθερός σπόρος για αναπαραγωγιμότητα ολόκληρου του batch

wgs84 = wgs84Ellipsoid;

%% ------------------ Κέντρο περιοχής & εύρη τυχαιοποίησης τοπολογίας ------------------
% Ίδιο σημείο αναφοράς με το test_simulation.m (περιοχή Αθήνας)
baseLat = 37.9838;
baseLon = 23.7275;

numBsRange          = [1 4];      % πλήθος BS ανά σενάριο
numUsersRange       = [3 10];     % πλήθος χρηστών ανά σενάριο
bsClusterRadiusKm   = 5;          % ακτίνα τοποθέτησης BS γύρω από το κέντρο
userRadiusKmRange   = [0.01 5];   % απόσταση χρήστη από τον BS αναφοράς του
satLatJitterDeg     = 10;         % τυχαιοποίηση γεωγρ. πλάτους υποδορυφορικού σημείου
satLonJitterDeg     = 10;         % τυχαιοποίηση γεωγρ. μήκους υποδορυφορικού σημείου
satAltitudeM        = 600e3;      % m, LEO-600 (TR 38.821 Πίν. 6.1.1.1-1)

%% ------------------ Σταθερό ραδιο-configuration (ίδιο με test_simulation.m) ------------------
% Επίγειο: TR 38.901 §7.8 Πίν. 7.8-1. BS gain = element (8 dBi) + array
% 10·log10(N), N=10, ιδανική στόχευση -> composite 18 dBi. TxPower &
% EIRP ορίζονται ανά σενάριο μέσα στον βρόχο (UMa 49 / UMi 44 dBm).
simParametersBase.Carrier = nrCarrierConfig;
simParametersBase.Carrier.NSizeGrid = 51;
simParametersBase.Carrier.SubcarrierSpacing = 30;
simParametersBase.Carrier.CyclicPrefix = 'Normal';
simParametersBase.CarrierFrequency = 3.5e9;
simParametersBase.AntennaGain = 8;              % dBi, element gain (TR 38.901 §7.3 Πίν. 7.3-1)
simParametersBase.NumAntennaElements = 10;      % TR 38.901 Πίν. 7.8-1
simParametersBase.RxNoiseFigure = 9;           % dB, UE downlink NF (TR 38.901 Πίν. 7.8-1)
simParametersBase.RxAntTemperature = 290;

simParametersBase.PathLossModel = '5G-NR';
simParametersBase.PathLoss = nrPathLossConfig;

simParametersBase.Power.NumTrx = 4;   % αλυσίδες πομποδέκτη ανά τομέα· P_out/αλυσίδα <= 20 W (EARTH Πίν. 2)
simParametersBase.Power.P0     = 130;
simParametersBase.Power.DeltaP = 4.7;
simParametersBase.Power.Psleep = 75;

% Δορυφόρος: 3GPP TR 38.821 Set-1, LEO-600, S-band (Πίνακες 6.1.1.1-1 & 6.1.3.2-1).
% EIRP density (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
satParametersBase.CarrierFrequency = 2.0e9;
satParametersBase.Bandwidth = 30e6;
satParametersBase.AntennaGain = 30;
satParametersBase.EirpDensityDbwPerMHz = 34;
satParametersBase.EIRP = satParametersBase.EirpDensityDbwPerMHz + 10*log10(satParametersBase.Bandwidth/1e6) + 30;
satParametersBase.TxPower = satParametersBase.EIRP - satParametersBase.AntennaGain;
satParametersBase.MinElevationDeg = 20;
satParametersBase.Power.Pfix  = 0;   % W, εκτός ενισχυτή· Pfix=0 -> αισιόδοξη υπόθεση
satParametersBase.Power.EtaPA = 0.4;

%% ------------------ Κύριος βρόχος Monte-Carlo ------------------
allScenarioTables = cell(numScenarios,1);

for s = 1:numScenarios
    numBs    = randi(numBsRange);
    numUsers = randi(numUsersRange);

    if rand() < 0.5
        scenarioType = 'UMa';
        bs_height_m = 25;
        simParametersBase.TxPower = 49;     % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMa)
    else
        scenarioType = 'UMi';
        bs_height_m = 10;
        simParametersBase.TxPower = 44;     % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMi)
    end
    simParametersBase.PathLoss.Scenario = scenarioType;
    simParametersBase.PathLoss.EnvironmentHeight = 1;
    simParametersBase.EIRP = simParametersBase.TxPower + simParametersBase.AntennaGain + ...
                             10*log10(simParametersBase.NumAntennaElements);   % dBm, composite

    % -- Τοποθέτηση BS γύρω από το κέντρο --
    bs_geo = zeros(numBs,3);
    for b = 1:numBs
        [dLat, dLon] = randOffsetDeg(baseLat, bsClusterRadiusKm);
        bs_geo(b,:) = [baseLat + dLat, baseLon + dLon, bs_height_m];
    end

    % -- Τοποθέτηση χρηστών εντός της κυψέλης ενός BS αναφοράς --
    % Όλοι οι χρήστες βρίσκονται εντός του πεδίου ισχύος των UMa/UMi, ώστε η
    % επιλογή κόμβου να κρίνεται από την ποιότητα της ζεύξης και όχι από την
    % απόσταση. Το εύρος [10 m, 5 km] καλύπτει και τις δύο πλευρές του σημείου
    % όπου το επίγειο SNR συναντά το δορυφορικό.
    user_geo = zeros(numUsers,3);
    for u = 1:numUsers
        refB = randi(numBs);
        [dLat, dLon] = randOffsetDeg(bs_geo(refB,1), userRadiusKmRange);
        user_geo(u,:) = [bs_geo(refB,1) + dLat, bs_geo(refB,2) + dLon, 1.5];
    end

    % -- Υποδορυφορικό σημείο (υψόμετρο σταθερό στο LEO-600) --
    sat_geo = [baseLat + (2*rand()-1)*satLatJitterDeg, ...
               baseLon + (2*rand()-1)*satLonJitterDeg, ...
               satAltitudeM];

    % -- Εκτέλεση σεναρίου --
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVecAll, satPathLossVec, satSnrDbVecAll] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParametersBase, satParametersBase);

    % Φορτίο κόμβου = πλήθος χρηστών του σεναρίου στον ίδιο κόμβο (feature για
    % το Part 2). Outage -> NodeLoad=0.
    nodeLoadVec = nan(numUsers,1);
    for u = 1:numUsers
        if bestNodeTypeVec(u) == "Outage"
            nodeLoadVec(u) = 0;
        else
            nodeLoadVec(u) = sum(bestNodeVec == bestNodeVec(u));
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
        energyPerBitUJVec, nodeLoadVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSnrDbVecAll, satElevationVecAll, satSlantRangeVec, satPathLossVec, ...
        'VariableNames', {'ScenarioID','UserID','ScenarioType','NumBS','NumUsers', ...
        'UserLat','UserLon','ServingNode','ServingType','Distance_m','PathLoss_dB', ...
        'SNR_dB','Capacity_Mbps','SatElevation_deg','NodePower_W','EnergyPerBit_uJ','NodeLoad', ...
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
numOutage      = sum(datasetTable.ServingType == "Outage");
fprintf('Monte-Carlo dataset: %d σενάρια, %d γραμμές χρηστών (%d Terrestrial, %d Satellite, %d Outage) -> %s\n', ...
    numScenarios, height(datasetTable), numTerrestrial, numSatellite, numOutage, outputCsvPath);

%% ------------------ Versioning αποτελεσμάτων ------------------
runParams = struct();
runParams.numScenarios   = numScenarios;
runParams.numUserRows    = height(datasetTable);
runParams.rngSeed        = 42;
runParams.baseLat        = baseLat;
runParams.baseLon        = baseLon;
runParams.numBsRange     = numBsRange;
runParams.numUsersRange  = numUsersRange;
runParams.bsClusterRadiusKm    = bsClusterRadiusKm;
runParams.userRadiusKmRange    = userRadiusKmRange;
runParams.satLatJitterDeg      = satLatJitterDeg;
runParams.satLonJitterDeg      = satLonJitterDeg;
runParams.satAltitudeM         = satAltitudeM;
runParams.terrestrial   = struct('CarrierFrequency_Hz', simParametersBase.CarrierFrequency, ...
    'TxPower_dBm_UMa', 49, 'TxPower_dBm_UMi', 44, ...
    'AntennaGain_dBi', simParametersBase.AntennaGain, ...
    'NumAntennaElements', simParametersBase.NumAntennaElements, ...
    'RxNoiseFigure_dB', simParametersBase.RxNoiseFigure, ...
    'RxAntTemperature_K', simParametersBase.RxAntTemperature, ...
    'NSizeGrid', simParametersBase.Carrier.NSizeGrid, ...
    'SubcarrierSpacing_kHz', simParametersBase.Carrier.SubcarrierSpacing);
runParams.terrestrial.Power = simParametersBase.Power;
runParams.satellite     = satParametersBase;

saveRunVersion('monteCarloDriver', runParams, {outputCsvPath}, label);

end

function [dLat, dLon] = randOffsetDeg(refLat, radiusKm)
% Τυχαία μετατόπιση [dLat, dLon] σε μοίρες, ομοιόμορφα ως προς το εμβαδόν.
% Το radiusKm είναι είτε βαθμωτό (δίσκος [0, R]) είτε ζεύγος [rmin rmax]
% (δακτύλιος): r = sqrt(rmin^2 + U*(rmax^2 - rmin^2)), theta = 2*pi*V.
if isscalar(radiusKm)
    rMin = 0; rMax = radiusKm;
else
    rMin = radiusKm(1); rMax = radiusKm(2);
end
r = sqrt(rMin^2 + rand()*(rMax^2 - rMin^2));
theta = 2*pi*rand();
dNorthKm = r*cos(theta);
dEastKm  = r*sin(theta);
dLat = dNorthKm / 110.574;               % km ανά μοίρα γεωγρ. πλάτους
dLon = dEastKm / (111.320*cosd(refLat));  % km ανά μοίρα γεωγρ. μήκους στο refLat
end
