clc;
clear;
rng(42); % Σταθερός σπόρος RNG για αναπαραγώγιμα αποτελέσματα (LOS draw + shadow fading είναι πλέον στοχαστικά)
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
sat_geo = [38.0200 23.8200 550e3];   % LEO (550 km)

% User calculations
numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

% WGS84 spheroid
wgs84 = wgs84Ellipsoid;

%% ------------------ Parameters (Terrestrial NR - 3GPP TR 38.901) ------------------
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;     % FR1
simParameters.TxPower = 43;                 % dBm ανά BS
simParameters.RxNoiseFigure = 5;            % dB
simParameters.RxAntTemperature = 290;       % K

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;

% -- Επιλογή Σεναρίου βάσει TR 38.901 --
scenarioType = 'UMa';
simParameters.PathLoss.Scenario = scenarioType;

switch scenarioType
    case 'UMa'
        bs_height_m = 25; 
        simParameters.PathLoss.EnvironmentHeight = 1; 
    case 'UMi'
        bs_height_m = 10;
        simParameters.PathLoss.EnvironmentHeight = 1;       
end

% Ενημέρωση των υψομέτρων των BS στον πίνακα bs_geo αυτόματα βάσει σεναρίου
bs_geo(:, 3) = bs_height_m;

%% ------------------ Parameters (Satellite) ------------------
satParameters.CarrierFrequency = 2.01e9;     % S-band
satParameters.TxPower = 34;                 % dBm (Ισχύς ενισχυτή)
satParameters.AntennaGain = 30;             % dBi (Κέρδος κατευθυντικής κεραίας LEO, TR 38.821)
satParameters.EIRP = satParameters.TxPower + satParameters.AntennaGain; 
satParameters.Bandwidth = 20e6;             % Hz
satParameters.MinElevationDeg = 10;         % visibility mask

%% ------------------ Parameters (Ενεργειακό μοντέλο) ------------------
% Γραμμικό μοντέλο κατανάλωσης ισχύος EARTH (Auer et al., "How much energy
% is needed to run a wireless network?", IEEE Wireless Commun., 2011) για
% τον σταθμό βάσης: P = NumTrx*(P0 + DeltaP*Pout) σε ενεργή λειτουργία,
% NumTrx*Psleep σε αδράνεια (τιμές αναφοράς macro cell, Pmax=20W <-> 43dBm
% ήδη ίδιο με το TxPower του σεναρίου).
simParameters.Power.NumTrx = 1;      % Αριθμός TRX ανά BS (μονο-sector μοντέλο)
simParameters.Power.P0     = 130;    % W, σταθερή κατανάλωση σε ενεργή λειτουργία
simParameters.Power.DeltaP = 4.7;    % κλίση κατανάλωσης ισχύος ως προς Pout
simParameters.Power.Psleep = 75;     % W, κατανάλωση σε αδράνεια (δεν χρησιμοποιείται ακόμα
                                      % στο per-user proxy - προορίζεται για μελλοντικό
                                      % network-wide accounting αδρανών κόμβων)

% Γραμμικό μοντέλο ενισχυτή ισχύος (PA) για τον δορυφόρο: P = Pfix + Pout/EtaPA
satParameters.Power.Pfix  = 0;       % W, σταθερή κατανάλωση εκτός ενισχυτή (μη τυποποιημένη
                                      % τιμή για payload - συντηρητική προσέγγιση 0)
satParameters.Power.EtaPA = 0.4;     % Απόδοση ενισχυτή ισχύος (τυπικό εύρος 0.35-0.5 SSPA/TWTA)

%% ------------------ Εκτέλεση σεναρίου (επιλογή κόμβου + χωρητικότητα) ------------------
[bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

% Συνάρτηση για εμφάνιση του πίνακα (custom συνάρτηση χρήστη).
array(numUsers, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec)

% Call the visualization
visual(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec, bestNodeVec)