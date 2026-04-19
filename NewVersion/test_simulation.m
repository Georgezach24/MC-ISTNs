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
    %(Για τώρα το σύστημα αφού υπολογίσει και επιλέξει τον καλύτερο κόμβο
    %βάση SNR απλά υπολογίζει το capacity αυτού μετά)
    % ToDo: να το βάλω να υπολογίζει για κάθε snr και το capacity και να
    % συγκρίνει μετά και το καλύτερο capacity.
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


%% ------------------ 3D Visualization ------------------
figure('Name', '3D Terrestrial & NTN Network', 'Color', 'w', 'Position', [100, 100, 900, 700]);
hold on; grid on;

% Ορίζουμε ως σημείο αναφοράς το 1ο BS
lat0 = bs_geo(1,1);
lon0 = bs_geo(1,2);
h0   = 0;

% Μετατροπή μόνο για X, Y (αγνοούμε το Z του ENU λόγω καμπυλότητας της γης κσι νσ μπορέσουμε ετσι να δούμε όλους τους χρήστες)
[xBS, yBS, ~] = geodetic2enu(bs_geo(:,1), bs_geo(:,2), bs_geo(:,3), lat0, lon0, h0, wgs84);
[xUE, yUE, ~] = geodetic2enu(user_geo(:,1), user_geo(:,2), user_geo(:,3), lat0, lon0, h0, wgs84);
[xSat, ySat, ~] = geodetic2enu(sat_geo(:,1), sat_geo(:,2), sat_geo(:,3), lat0, lon0, h0, wgs84);


zBS  = bs_geo(:,3);
zUE  = user_geo(:,3);
zSat = sat_geo(:,3);

% Σχεδίαση Base Stations (Κόκκινα τρίγωνα)
scatter3(xBS, yBS, zBS, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Base Stations (TN)');
for b = 1:numBs
    text(xBS(b), yBS(b), zBS(b)*1.3, " BS" + b, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r'); 
end

% Σχεδίαση Χρηστών (Μπλε κύκλοι)
scatter3(xUE, yUE, zUE, 80, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Users (UE)');
for u = 1:numUsers
    text(xUE(u), yUE(u), zUE(u)*1.5, " U" + u, 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
end

% Σχεδίαση Δορυφόρου (Κίτρινο αστέρι)
scatter3(xSat, ySat, zSat, 300, 'p', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'DisplayName', 'Satellite (Real Scale)');
text(xSat(1), ySat(1), zSat(1)*1.2, " LEO Sat", 'FontSize', 11, 'FontWeight', 'bold');

% Σχεδίαση Γραμμών Σύνδεσης (Με Interpolation για να έχουμε λογαριθμική κλίμακα και να χωράνε όλα στο γράφιμα)
for u = 1:numUsers
    if bestNodeTypeVec(u) == "Terrestrial"
        bs_idx = str2double(extractAfter(bestNodeVec(u), "BS"));
        
        num_pts = 50;
        xq = linspace(xUE(u), xBS(bs_idx), num_pts);
        yq = linspace(yUE(u), yBS(bs_idx), num_pts);
        zq = linspace(zUE(u), zBS(bs_idx), num_pts); 
        
        plot3(xq, yq, zq, 'g-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
    elseif bestNodeTypeVec(u) == "Satellite"
        num_pts = 100;
        xq = linspace(xUE(u), xSat, num_pts);
        yq = linspace(yUE(u), ySat, num_pts);
        
        % Logarithmic interpolation για τον άξονα Z
        zq = logspace(log10(zUE(u)), log10(zSat(1)), num_pts);
        
        plot3(xq, yq, zq, 'm-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

% Dummy plots για το Legend
plot3(nan, nan, nan, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Terrestrial Link (Green)');
plot3(nan, nan, nan, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Satellite Link (Magenta)');

% --- Μορφοποίηση Γραφήματος και Αξόνων ---
xlabel('East (meters)');
ylabel('North (meters)');
zlabel('Altitude (meters) - Log Scale');
title('3D Network Simulation: Real Altitudes with Logarithmic Z-Axis');

% Εφαρμογή Λογαριθμικής Κλίμακας MONO στον άξονα Z
set(gca, 'ZScale', 'log');

% Ορίζουμε το Z-axis να ξεκινάει κάτω από το 1.5m του χρήστη για να μη "κοπεί" το marker
zlim([0.5, 10^6]); 

legend('Location', 'northeastoutside');

% Έξυπνη γωνία θέασης
az = -25; 
el = 12;
view(az, el);
grid minor;