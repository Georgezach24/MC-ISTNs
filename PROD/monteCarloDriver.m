function datasetTable = monteCarloDriver(numScenarios, outputCsvPath)
%MONTECARLODRIVER Παράγει ένα labeled dataset τρέχοντας το simulateScenario
% πολλές φορές πάνω σε τυχαιοποιημένες τοπολογίες: θέσεις BS/χρηστών,
% πλήθος BS/χρηστών, γεωμετρία δορυφόρου (θέση/υψόμετρο), και σενάριο
% TR 38.901 (UMa/UMi). Το ραδιο-configuration (ισχύς, bandwidth, μοντέλα
% ισχύος κτλ.) παραμένει σταθερό - ίδιο με το test_simulation.m - ώστε το
% dataset να παραμένει συγκρίσιμο με το single-run σενάριο αναφοράς.
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

rng(42); % Σταθερός σπόρος για αναπαραγωγιμότητα ολόκληρου του batch

wgs84 = wgs84Ellipsoid;

%% ------------------ Κέντρο περιοχής & εύρη τυχαιοποίησης τοπολογίας ------------------
% Ίδιο σημείο αναφοράς με το test_simulation.m (περιοχή Αθήνας)
baseLat = 37.9838;
baseLon = 23.7275;

numBsRange          = [1 4];      % πλήθος BS ανά σενάριο
numUsersRange       = [3 10];     % πλήθος χρηστών ανά σενάριο
bsClusterRadiusKm   = 5;          % ακτίνα τοποθέτησης BS γύρω από το κέντρο
nearUserRadiusKm    = 1;          % ακτίνα "κοντινών" χρηστών γύρω από τυχαίο BS του σεναρίου
                                   % (τυπική κάλυψη UMi/UMa ISD - μεγαλύτερη ακτίνα οδηγεί
                                   % συστηματικά σε NLOS λόγω της TR 38.901 §7.4.2 LOS
                                   % probability, ευνοώντας τεχνητά τον δορυφόρο)
farUserRadiusKmRange = [50 150];  % εύρος απόστασης "μακρινών" χρηστών (πρακτικά εκτός εμβέλειας BS)
farUserProbability  = 0.3;        % πιθανότητα ένας χρήστης να τοποθετηθεί μακριά
satLatJitterDeg     = 3;          % τυχαιοποίηση γεωγρ. πλάτους υποδορυφορικού σημείου
satLonJitterDeg     = 3;          % τυχαιοποίηση γεωγρ. μήκους υποδορυφορικού σημείου
satAltitudeRangeM   = [500e3 600e3]; % τυπικό εύρος υψομέτρου LEO

%% ------------------ Σταθερό ραδιο-configuration (ίδιο με test_simulation.m) ------------------
simParametersBase.Carrier = nrCarrierConfig;
simParametersBase.Carrier.NSizeGrid = 51;
simParametersBase.Carrier.SubcarrierSpacing = 30;
simParametersBase.Carrier.CyclicPrefix = 'Normal';
simParametersBase.CarrierFrequency = 3.5e9;
simParametersBase.TxPower = 43;
simParametersBase.AntennaGain = 8;   % dBi, BS antenna element gain (TR 38.901 §7.3, Table 7.3-1, G_E,max)
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

    % -- Τοποθέτηση χρηστών: είτε κοντά σε τυχαίο BS, είτε μακριά (μόνο δορυφόρος) --
    user_geo = zeros(numUsers,3);
    for u = 1:numUsers
        if rand() < farUserProbability
            radiusKm = farUserRadiusKmRange(1) + diff(farUserRadiusKmRange)*rand();
            [dLat, dLon] = randOffsetDeg(baseLat, radiusKm);
            user_geo(u,:) = [baseLat + dLat, baseLon + dLon, 1.5];
        else
            refB = randi(numBs);
            [dLat, dLon] = randOffsetDeg(bs_geo(refB,1), nearUserRadiusKm);
            user_geo(u,:) = [bs_geo(refB,1) + dLat, bs_geo(refB,2) + dLon, 1.5];
        end
    end

    % -- Τυχαιοποίηση γεωμετρίας δορυφόρου (υποδορυφορικό σημείο + υψόμετρο) --
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

    % Φορτίο κόμβου: πλήθος χρηστών του ΙΔΙΟΥ σεναρίου που εξυπηρετούνται
    % από τον ίδιο ΣΥΓΚΕΚΡΙΜΕΝΟ κόμβο (χρήσιμο ως feature "node load" για
    % το Part 2). Ένας DualConnectivity χρήστης φορτίζει ΚΑΙ τον BS ΚΑΙ
    % τον δορυφόρο του ταυτόχρονα (simulateScenario.m), οπότε μία κοινή
    % στήλη "NodeLoad" δεν αρκεί πλέον - καταγράφονται δύο ξεχωριστές
    % στήλες, BsLoad και SatLoad, καθεμία 0 αν ο χρήστης δεν συνδέεται σε
    % εκείνον τον τύπο κόμβου (π.χ. outage, ή Satellite-only για BsLoad).
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
% Τυχαία μετατόπιση [dLat, dLon] σε μοίρες, ομοιόμορφα κατανεμημένη εντός
% δίσκου ακτίνας radiusKm γύρω από σημείο αναφοράς πλάτους refLat.
% Επίπεδη (Ευκλείδεια) προσέγγιση γύρω από το refLat - αρκετή ακρίβεια
% στην τοπική κλίμακα των km που χρησιμοποιείται εδώ.
r = radiusKm * sqrt(rand());
theta = 2*pi*rand();
dNorthKm = r*cos(theta);
dEastKm  = r*sin(theta);
dLat = dNorthKm / 110.574;               % km ανά μοίρα γεωγρ. πλάτους
dLon = dEastKm / (111.320*cosd(refLat));  % km ανά μοίρα γεωγρ. μήκους στο refLat
end
