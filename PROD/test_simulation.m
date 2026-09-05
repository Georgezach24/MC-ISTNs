clc;
clear;
rng(42);            % σταθερός σπόρος RNG για αναπαραγωγιμότητα
runLabel = '';      % προαιρετικό tag στο όνομα του versioned φακέλου αποτελεσμάτων
%% ------------------ Γεωγραφικές θέσεις [lat lon h(m)] ------------------
% Παράδειγμα συντεταγμένων κοντά στην Αθήνα
% BS: [latitude, longitude, height_m] (Τα ύψη θα ενημερωθούν αυτόματα από το σενάριο)
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

% Users: [latitude, longitude, height_m]
user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            37.0380 23.9550 1.5;
            38.0500 23.9500 1.5];

% SATs: [latitude, longitude, altitude_m]
sat_geo = [38.0200 23.8200 600e3];   % LEO-600 (3GPP TR 38.821 Πίνακας 6.1.1.1-1, στήλη LEO-600)

% User calculations
numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

% WGS84 spheroid
wgs84 = wgs84Ellipsoid;

%% ------------------ Parameters (Terrestrial NR) ------------------
% Reference σύνολο: TR 38.901 §7.8 Πίνακας 7.8-1 (large-scale calibration).
% BS antenna gain = element gain (8 dBi, §7.3 Πίν. 7.3-1) + array gain
% 10·log10(N) με N=10 στοιχεία ((M,N,P)=(10,1,1), single port) υπό παραδοχή
% ιδανικής στόχευσης δέσμης -> composite 18 dBi. UE gain = 0 dBi (isotropic).
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;     % Hz, n78
simParameters.AntennaGain = 8;              % dBi, element gain (TR 38.901 §7.3 Πίν. 7.3-1)
simParameters.NumAntennaElements = 10;      % TR 38.901 Πίν. 7.8-1, (M,N,P)=(10,1,1)
simParameters.RxNoiseFigure = 9;            % dB, UE downlink NF (TR 38.901 Πίν. 7.8-1)
simParameters.RxAntTemperature = 290;       % K

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;

% -- Επιλογή Σεναρίου βάσει TR 38.901 --
scenarioType = 'UMa';
simParameters.PathLoss.Scenario = scenarioType;

switch scenarioType
    case 'UMa'
        bs_height_m = 25;
        simParameters.TxPower = 49;         % dBm, conducted (TR 38.901 Πίν. 7.8-1)
        simParameters.PathLoss.EnvironmentHeight = 1;
    case 'UMi'
        bs_height_m = 10;
        simParameters.TxPower = 44;         % dBm, conducted (TR 38.901 Πίν. 7.8-1)
        simParameters.PathLoss.EnvironmentHeight = 1;
end

% Composite BS EIRP: conducted power + element gain + array gain (ιδανική στόχευση)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);   % dBm

% Ενημέρωση των υψομέτρων των BS στον πίνακα bs_geo αυτόματα βάσει σεναρίου
bs_geo(:, 3) = bs_height_m;

%% ------------------ Parameters (Satellite) ------------------
% Reference σύνολο: 3GPP TR 38.821, Set-1, LEO-600, S-band.
%   Πίνακας 6.1.1.1-1: EIRP density, Tx max gain, altitude.
%   Πίνακας 6.1.3.2-1: carrier frequency, system bandwidth (link budget).
% Δεδομένο του προτύπου είναι η πυκνότητα EIRP (dBW/MHz), όχι η ισχύς RF·
% το EIRP προκύπτει από το bandwidth και το TxPower ως EIRP - AntennaGain.
satParameters.CarrierFrequency = 2.0e9;      % Hz
satParameters.Bandwidth = 30e6;              % Hz, system bandwidth S-band
satParameters.AntennaGain = 30;              % dBi, Tx max gain LEO-600 S-band
satParameters.EirpDensityDbwPerMHz = 34;     % dBW/MHz
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + 10*log10(satParameters.Bandwidth/1e6) + 30;  % dBm
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;  % dBm, ισχύς RF στην είσοδο κεραίας
satParameters.MinElevationDeg = 10;          % visibility mask

%% ------------------ Parameters (Ενεργειακό μοντέλο) ------------------
% BS: μοντέλο EARTH (Auer et al. 2011), P = NumTrx*(P0 + DeltaP*Pout).
simParameters.Power.NumTrx = 1;      % TRX ανά BS
simParameters.Power.P0     = 130;    % W, σταθερή κατανάλωση σε ενεργή λειτουργία
simParameters.Power.DeltaP = 4.7;    % κλίση ως προς Pout
simParameters.Power.Psleep = 75;     % W, αδράνεια (δεν χρησιμοποιείται ακόμα)

% Δορυφόρος: γραμμικό μοντέλο PA, P = Pfix + Pout/EtaPA.
satParameters.Power.Pfix  = 0;       % W, κατανάλωση εκτός ενισχυτή· Pfix=0 -> αισιόδοξη υπόθεση
satParameters.Power.EtaPA = 0.4;     % απόδοση ενισχυτή (τυπικό εύρος 0.35-0.5)

%% ------------------ Εκτέλεση σεναρίου (επιλογή κόμβου + χωρητικότητα) ------------------
[bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

% Συνάρτηση για εμφάνιση του πίνακα (custom συνάρτηση χρήστη).
array(numUsers, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec)

% Call the visualization
visual(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec, bestNodeVec)

%% ------------------ Versioning αποτελεσμάτων ------------------
resultsTable = table((1:numUsers)', bestNodeVec, bestNodeTypeVec, bestDistanceVec, ...
    bestPathLossVec, bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    'VariableNames', {'User','ServingNode','ServingType','Distance_m','PathLoss_dB', ...
    'SNR_dB','Capacity_Mbps','SatElevation_deg','NodePower_W','EnergyPerBit_uJ'});
tmpDir = tempname; mkdir(tmpDir);
resultsCsv = fullfile(tmpDir, 'results.csv');
writetable(resultsTable, resultsCsv);
figPng = fullfile(tmpDir, 'network_3d.png');
saveas(gcf, figPng);

runParams = struct();
runParams.rngSeed          = 42;
runParams.scenarioType     = scenarioType;
runParams.bs_geo           = bs_geo;
runParams.user_geo         = user_geo;
runParams.sat_geo          = sat_geo;
runParams.terrestrial      = struct('CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'TxPower_dBm', simParameters.TxPower, 'AntennaGain_dBi', simParameters.AntennaGain, ...
    'NumAntennaElements', simParameters.NumAntennaElements, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'RxAntTemperature_K', simParameters.RxAntTemperature, ...
    'NSizeGrid', simParameters.Carrier.NSizeGrid, ...
    'SubcarrierSpacing_kHz', simParameters.Carrier.SubcarrierSpacing, ...
    'Scenario', simParameters.PathLoss.Scenario);
runParams.terrestrial.Power = simParameters.Power;
runParams.satellite        = satParameters;

saveRunVersion('test_simulation', runParams, {resultsCsv, figPng}, runLabel);
rmdir(tmpDir, 's');