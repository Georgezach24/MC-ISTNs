clc;
clear;

%% ------------------ Γεωγραφικές θέσεις [lat lon h(m)] ------------------
% Παράδειγμα συντεταγμένων κοντά στην Αθήνα
% Βάλε εδώ τις δικές σου πραγματικές συντεταγμένες

% BS: [latitude, longitude, height_m]
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

% Users: [latitude, longitude, height_m]
user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            38.0500 23.9500 1.5];


% [latitude, longitude, altitude_m]
sat_geo = [38.0200 23.8200 550e3];   % LEO (550 km)

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

% WGS84 spheroid
wgs84 = wgs84Ellipsoid;

%% ------------------ Parameters (Terrestrial NR) ------------------
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';

simParameters.CarrierFrequency = 3.5e9;
simParameters.TxPower = 40;                 % dBm ανά BS
simParameters.RxNoiseFigure = 6;            % dB
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

    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestBW        = NaN;
    userBestElevation = NaN;

    %% ===== Terrestrial BS candidates =====
    for b = 1:numBs

        % Χρησιμοποιούμε τοπικό ENU σύστημα με origin:
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

        fprintf("User %d - BS %d: GroundDist = %.2f m, 3D Range = %.2f m, PathLoss = %.2f dB, SNR = %.2f dB\n", ...
                u, b, groundDistance, d3d, pathLoss, snr_db);

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

    if elevSat >= satParameters.MinElevationDeg
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency;
        satPathLoss = fspl(slantRangeSat, lambdaSat);
        satSnrDb = (satParameters.TxPower - 30) ...
                   - satPathLoss ...
                   - noisePowerSAT_dBW;
    else
        satPathLoss = inf;
        satSnrDb = -Inf;
    end

    satPathLossVec(u) = satPathLoss;
    satSnrDbVec(u)    = satSnrDb;

    fprintf("User %d - SAT: Elevation = %.2f deg, SlantRange = %.2f m, PathLoss = %.2f dB, SNR = %.2f dB\n", ...
            u, elevSat, slantRangeSat, satPathLoss, satSnrDb);

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

    fprintf("=> User %d συνδέεται στο %s (%s) με SNR = %.2f dB\n\n", ...
            u, userBestNode, userBestType, userBestSNR);
end

%% ------------------ Πίνακας αποτελεσμάτων ------------------
userID = (1:numUsers)';

resultsTable = table(userID, bestNodeVec, bestNodeTypeVec, ...
    bestDistanceVec, bestPathLossVec, bestSnrDbVec, ...
    capacityMbpsVec, bestElevationDegVec, ...
    'VariableNames', {'User','ServingNode','ServingType', ...
    'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps','SatElevation_deg'});

disp(resultsTable)