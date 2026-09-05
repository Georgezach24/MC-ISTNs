function T = temporalPassSimulation(dtSeconds, outputDir, label)
%TEMPORALPASSSIMULATION Τρέχει το στατικό σενάριο του test_simulation.m
% (2 BS, 6 χρήστες) επαναλαμβανόμενα σε διαδοχικά χρονικά βήματα,
% μετακινώντας το υποδορυφορικό σημείο του LEO κατά μήκος ενός
% απλοποιημένου ground track (σταθερό γεωγρ. πλάτος, ground-track speed
% κυκλικής Κεπλεριανής τροχιάς, χωρίς περιστροφή Γης - όχι πλήρης
% ορβιτογράφος). Η κατάσταση καναλιού περνάει από βήμα σε βήμα (χωρικά
% συσχετισμένο shadow fading, βλ. correlatedLosState).
%
% Χρήση:
%   T = temporalPassSimulation();              % dt = 5s -> ../Results
%   T = temporalPassSimulation(2);             % dt = 2s
%   T = temporalPassSimulation(5, 'C:\out')
%   T = temporalPassSimulation(5, [], 'tag')   % tag στο όνομα του versioned φακέλου

if nargin < 1 || isempty(dtSeconds)
    dtSeconds = 5;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if nargin < 3
    label = '';
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% ------------------ Ίδια στατική τοπολογία (BS + χρήστες) με το test_simulation.m ------------------
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            37.0380 23.9550 1.5;
            38.0500 23.9500 1.5];

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);
wgs84 = wgs84Ellipsoid;

%% ------------------ Ραδιο-configuration (ίδιο με test_simulation.m) ------------------
% Επίγειο: TR 38.901 §7.8 Πίν. 7.8-1. BS gain = element (8 dBi) + array
% 10·log10(N), N=10, ιδανική στόχευση -> composite 18 dBi.
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.AntennaGain = 8;              % dBi, element gain (TR 38.901 §7.3 Πίν. 7.3-1)
simParameters.NumAntennaElements = 10;      % TR 38.901 Πίν. 7.8-1
simParameters.RxNoiseFigure = 9;           % dB, UE downlink NF (TR 38.901 Πίν. 7.8-1)
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;                 % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMa)
bs_geo(:,3) = 25;

simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);   % dBm, composite

% Δορυφόρος: 3GPP TR 38.821 Set-1, LEO-600, S-band (Πίνακες 6.1.1.1-1 & 6.1.3.2-1).
% EIRP density (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
satAltitude = 600e3;                          % m, LEO-600
satParameters.CarrierFrequency = 2.0e9;
satParameters.Bandwidth = 30e6;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + 10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 10;

simParameters.Power.NumTrx = 1;
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.Power.Pfix  = 0;
satParameters.Power.EtaPA = 0.4;

%% ------------------ Ground track του LEO (απλοποιημένο μοντέλο διέλευσης) ------------------
muEarth = 3.986004418e14;  % m^3/s^2, βαρυτική παράμετρος Γης
Re      = 6371e3;          % m, μέση ακτίνα Γης
a       = Re + satAltitude;
orbitalPeriodS  = 2*pi*sqrt(a^3/muEarth);      % Κεπλεριανή περίοδος (s)
groundSpeedMps  = (2*pi/orbitalPeriodS) * Re;  % ταχύτητα ίχνους εδάφους (m/s)

centerLat = mean(bs_geo(:,1));
centerLon = mean(bs_geo(:,2));

startOffsetKm = -1500;  % km ανατολικά του κέντρου, αρχή της διέλευσης (δορυφόρος αόρατος)
maxSteps      = 2000;   % ασφαλιστικό όριο βημάτων

%% ------------------ Χρονικός βρόχος ------------------
allRows = cell(maxSteps,1);
wasVisible = false;
step = 0;
t = 0;
channelState = []; % καμία προηγούμενη κατάσταση πριν το πρώτο βήμα -> πρώτο δείγμα ανεξάρτητο (i.i.d.)

while step < maxSteps
    offsetKm = startOffsetKm + groundSpeedMps * t / 1000;
    subLon = centerLon + offsetKm / (111.320*cosd(centerLat));
    sat_geo = [centerLat, subLon, satAltitude];

    rng(step + 1); % RNG seed ανά χρονικό βήμα

    % channelState περνάει από βήμα σε βήμα -> χωρικά συσχετισμένο shadow fading.
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
        channelState] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, channelState);

    step = step + 1;
    userID = (1:numUsers)';
    tCol = repmat(t, numUsers, 1);
    stepCol = repmat(step, numUsers, 1);

    allRows{step} = table(stepCol, tCol, userID, bestNodeVec, bestNodeTypeVec, ...
        bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        satElevationVec, satSnrDbVec, satPathLossVec, ...
        bestBsSnrDbVec, ...
        'VariableNames', {'Step','Time_s','UserID','ServingNode','ServingType', ...
        'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps', ...
        'NodePower_W','EnergyPerBit_uJ', ...
        'SatElevation_deg','CandSat_SNR_dB','CandSat_PathLoss_dB', ...
        'CandBS_SNR_dB'});

    refElev = max(satElevationVec);
    if refElev >= satParameters.MinElevationDeg
        wasVisible = true;
    elseif wasVisible
        break; % τέλος διέλευσης (δορυφόρος αόρατος σε όλους μετά από ορατότητα)
    end

    t = t + dtSeconds;
end

T = vertcat(allRows{1:step});
writetable(T, fullfile(outputDir, 'temporal_pass_dataset.csv'));

fprintf('Temporal pass simulation: %d χρονικά βήματα (dt=%ds, ~%.0f s συνολικά), ταχύτητα ίχνους %.0f m/s (T_orbit=%.0f s)\n', ...
    step, dtSeconds, T.Time_s(end), groundSpeedMps, orbitalPeriodS);

%% ------------------ Handovers & outage events ανά χρήστη ------------------
% Μετάβαση προς/από ServingNode="None" μετράται ως outage event, όχι handover.
fprintf('\n--- Handovers & outage events ανά χρήστη ---\n');
handoverCounts    = zeros(numUsers,1);
outageEventCounts = zeros(numUsers,1);
for u = 1:numUsers
    userRows = T(T.UserID == u, :);
    userRows = sortrows(userRows, 'Step');
    prevNodes = userRows.ServingNode(1:end-1);
    currNodes = userRows.ServingNode(2:end);
    changed = currNodes ~= prevNodes;
    intoOutage = changed & (currNodes == "None");
    realHandover = changed & ~(prevNodes == "None" | currNodes == "None");

    handoverCounts(u)    = sum(realHandover);
    outageEventCounts(u) = sum(intoOutage);
    fprintf('  User %d: %d handovers, %d outage events (%s -> ... -> %s)\n', u, ...
        handoverCounts(u), outageEventCounts(u), ...
        userRows.ServingNode(1), userRows.ServingNode(end));
end

%% ------------------ Γραφήματα χρονοσειράς ανά χρήστη ------------------
kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB','SatElevation_deg'};
kpiLabel = {'Throughput (Mbps)','Energy per bit (\muJ/bit)','SNR (dB)','Sat. Elevation (deg)'};

for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    hold on;
    for u = 1:numUsers
        userRows = T(T.UserID == u, :);
        userRows = sortrows(userRows, 'Step');
        plot(userRows.Time_s, userRows.(kpiList{k}), 'DisplayName', sprintf('User %d', u));
    end
    hold off;
    xlabel('Χρόνος (s)');
    ylabel(kpiLabel{k});
    legend('Location','best');
    title(sprintf('%s κατά τη διάρκεια διέλευσης LEO', kpiList{k}), 'Interpreter','none');
    saveas(fig, fullfile(outputDir, ['temporal_' kpiList{k} '.png']));
    close(fig);
end

fprintf('\nTemporal pass simulation results -> %s\n', outputDir);

%% ------------------ Versioning αποτελεσμάτων ------------------
runParams = struct();
runParams.dtSeconds       = dtSeconds;
runParams.numSteps        = step;
runParams.rngScheme       = 'rng(step+1) ανά χρονικό βήμα';
runParams.startOffsetKm   = startOffsetKm;
runParams.satAltitude_m   = satAltitude;
runParams.orbitalPeriod_s = orbitalPeriodS;
runParams.groundSpeed_mps = groundSpeedMps;
runParams.scenarioType    = simParameters.PathLoss.Scenario;
runParams.bs_geo          = bs_geo;
runParams.user_geo        = user_geo;
runParams.terrestrial     = struct('CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'TxPower_dBm', simParameters.TxPower, 'AntennaGain_dBi', simParameters.AntennaGain, ...
    'NumAntennaElements', simParameters.NumAntennaElements, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'RxAntTemperature_K', simParameters.RxAntTemperature);
runParams.terrestrial.Power = simParameters.Power;
runParams.satellite       = satParameters;

outFiles = {fullfile(outputDir, 'temporal_pass_dataset.csv')};
for k = 1:numel(kpiList)
    outFiles{end+1} = fullfile(outputDir, ['temporal_' kpiList{k} '.png']); %#ok<AGROW>
end
saveRunVersion('temporalPassSimulation', runParams, outFiles, label);

end
