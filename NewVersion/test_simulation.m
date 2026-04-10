clc;
clear;

%% ------------------ Γεωγραφικές θέσεις [lat lon h(m)] ------------------
% Παράδειγμα συντεταγμένων κοντά στην Αθήνα

% BS: [latitude, longitude, height_m]
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

% Users: [latitude, longitude, height_m]
user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            38.0500 23.9500 1.5];


%SATs: [latitude, longitude, altitude_m]
sat_geo = [38.0200 23.8200 550e3];   % LEO (550 km)

%User calculations
numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

% WGS84 spheroid
wgs84 = wgs84Ellipsoid;

%% ------------------ Parameters (Terrestrial NR) ------------------
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
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;

%% ------------------ Parameters (Satellite) ------------------
satParameters.CarrierFrequency = 2.01e9;     %  S-band
satParameters.TxPower = 34;                 % dBm
satParameters.Bandwidth = 20e6;             % Hz
satParameters.MinElevationDeg = 10;         % visibility mask

%% ------------------ Fading / LOS για terrestrial model ------------------
simParameters.DelayProfile = 'TDL-A';

if contains(simParameters.DelayProfile,'CDL','IgnoreCase',true)
    channel = nrCDLChannel;
    channel.DelayProfile = simParameters.DelayProfile;
    chInfo = info(channel);
    kFactor = chInfo.KFactorFirstCluster;
else
    channel = nrTDLChannel;
    channel.DelayProfile = simParameters.DelayProfile;
    chInfo = info(channel);
    kFactor = chInfo.KFactorFirstTap;
end

simParameters.LOS = kFactor > -Inf;

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
bestElevationDegVec = nan(numUsers,1);

% Διαγνωστικοί πίνακες
groundDistanceMat = nan(numUsers,numBs);
range3DMat        = nan(numUsers,numBs);
pathLossMat       = nan(numUsers,numBs);
snrDbMat          = nan(numUsers,numBs);

satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSnrDbVec       = nan(numUsers,1);


%% ------------------ Loop για κάθε χρήστη ------------------
for u = 1:numUsers

    %Αρικοποιήσεις.
    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestBW        = NaN;
    userBestElevation = NaN;

    %% ===== Terrestrial BS candidates =====
    for b = 1:numBs

        % latitude/longitude του BS και h0 = 0 m
        lat0 = bs_geo(b,1);
        lon0 = bs_geo(b,2);
        h0   = 0;

        % BS σε ENU
        [xBS, yBS, zBS] = geodetic2enu(bs_geo(b,1), bs_geo(b,2), bs_geo(b,3), ...
                                       lat0, lon0, h0, wgs84);

        % User σε ENU ως προς το ίδιο origin
        [xUE, yUE, zUE] = geodetic2enu(user_geo(u,1), user_geo(u,2), user_geo(u,3), ...
                                       lat0, lon0, h0, wgs84);

        txPosition = [xBS; yBS; zBS];
        rxPosition = [xUE; yUE; zUE];

        % 2D γεωδαιτική απόσταση πάνω στη γη
        groundDistance = distance(bs_geo(b,1), bs_geo(b,2), ...
                                  user_geo(u,1), user_geo(u,2), ...
                                  wgs84);

        % 3D απόσταση
        d3d = norm(rxPosition - txPosition);

        groundDistanceMat(u,b) = groundDistance;
        range3DMat(u,b)        = d3d;

        % Terrestrial NR path loss
        pathLoss = nrPathLoss(simParameters.PathLoss, ...
                              simParameters.CarrierFrequency, ...
                              simParameters.LOS, ...
                              txPosition, rxPosition);

        pathLossMat(u,b) = pathLoss;

        % SNR
        snr_db = (simParameters.TxPower - 30) ...
                 - pathLoss ...
                 - noisePowerBS_dBW;

        snrDbMat(u,b) = snr_db;

        %Σύγκριση καλύτερου SNR.
        if snr_db > userBestSNR
            userBestSNR       = snr_db;
            userBestNode      = "BS" + string(b);
            userBestType      = "Terrestrial";
            userBestDistance  = d3d;
            userBestPathLoss  = pathLoss;
            userBestBW        = BW_bs;
            userBestElevation = NaN;
        end
    end

    %% ===== Satellite candidate =====
    [azSat, elevSat, slantRangeSat] = geodetic2aer( ...
        sat_geo(1), sat_geo(2), sat_geo(3), ...
        user_geo(u,1), user_geo(u,2), user_geo(u,3), ...
        wgs84);

    satSlantRangeVec(u) = slantRangeSat;
    satElevationVec(u)  = elevSat;

    if elevSat >= satParameters.MinElevationDeg %Αν ο δορυφόρος είναι ορατός από τον χρήστη.
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency; %Υπολογισμός (λ)
        satPathLoss = fspl(slantRangeSat, lambdaSat); %Free space loss για τον δορυφόρο.
        satSnrDb = (satParameters.TxPower - 30) ...
                   - satPathLoss ...
                   - noisePowerSAT_dBW;
    else %Αλλιώς θέσε τα ως άπειρα για να μην επιλέγονται ποτέ.
        satPathLoss = inf;
        satSnrDb = -Inf;
    end

    satPathLossVec(u) = satPathLoss;
    satSnrDbVec(u)    = satSnrDb;

    %Σύγκριση των SNR για TN και NTN.
    if satSnrDb > userBestSNR
        userBestSNR       = satSnrDb;
        userBestNode      = "SAT-1";
        userBestType      = "Satellite";
        userBestDistance  = slantRangeSat;
        userBestPathLoss  = satPathLoss;
        userBestBW        = satParameters.Bandwidth;
        userBestElevation = elevSat;
    end

    %% ===== Capacity =====
    snr_lin = 10^(userBestSNR/10);
    capacity = userBestBW * log2(1 + snr_lin);   % bits/s

    bestNodeVec(u)         = userBestNode;
    bestNodeTypeVec(u)     = userBestType;
    bestDistanceVec(u)     = userBestDistance;
    bestPathLossVec(u)     = userBestPathLoss;
    bestSnrDbVec(u)        = userBestSNR;
    capacityMbpsVec(u)     = capacity * 1e-6;
    bestElevationDegVec(u) = userBestElevation;

end

%% ------------------ Πίνακας αποτελεσμάτων ------------------
userID = (1:numUsers)';

resultsTable = table(userID, bestNodeVec, bestNodeTypeVec, ...
    bestDistanceVec, bestPathLossVec, bestSnrDbVec, ...
    capacityMbpsVec, bestElevationDegVec, ...
    'VariableNames', {'User','ServingNode','ServingType', ...
    'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps','SatElevation_deg'});

disp(resultsTable)