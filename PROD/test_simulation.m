clc;
clear;
rng(42); % Σταθερό seed - LOS draw/shadow fading είναι στοχαστικά
%% ------------------ Γεωγραφικές θέσεις [lat lon h(m)] ------------------
% Περιοχή Αθήνας. Ύψη BS ενημερώνονται αυτόματα βάσει σεναρίου.
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

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
simParameters.AntennaGain = 8;              % dBi, στοιχείο κεραίας (TR 38.901 §7.3, Πίνακας 7.3-1)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain; % dBm
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
% EARTH model (Auer et al. 2011) για BS: P=NumTrx*(P0+DeltaP*Pout).
simParameters.Power.NumTrx = 1;
simParameters.Power.P0     = 130;    % W
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;     % W, αδράνεια (δεν χρησιμοποιείται ακόμα)

% Γραμμικό μοντέλο ενισχυτή ισχύος δορυφόρου: P = Pfix + Pout/EtaPA
satParameters.Power.Pfix  = 0;       % W
satParameters.Power.EtaPA = 0.4;     % απόδοση PA (τυπικό εύρος 0.35-0.5)

%% ------------------ Εκτέλεση σεναρίου (επιλογή κόμβου + χωρητικότητα) ------------------
[bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

array(numUsers, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec)

visual(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec, bestNodeVec)