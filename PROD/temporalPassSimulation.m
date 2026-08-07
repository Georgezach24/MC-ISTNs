function T = temporalPassSimulation(dtSeconds, outputDir)
%TEMPORALPASSSIMULATION Τρέχει το ίδιο στατικό σενάριο (2 terrestrial BS,
% 6 χρήστες) του test_simulation.m αλλά επαναλαμβανόμενα σε διαδοχικά
% χρονικά βήματα (πολλαπλές "μεταδόσεις"), μετακινώντας το υποδορυφορικό
% σημείο του LEO κατά μήκος ενός απλοποιημένου ground track ώστε να
% αποτυπωθεί η κίνηση του δορυφόρου και η μεταβολή του καναλιού
% (elevation/path loss/SNR) προς κάθε χρήστη με τον χρόνο.
%
% Απλοποίηση γεωμετρίας πάσσου: το υποδορυφορικό σημείο κινείται σε
% σταθερό γεωγραφικό πλάτος (=πλάτος του κέντρου του BS cluster) προς
% ανατολάς, με σταθερή ταχύτητα ίση με την ταχύτητα εδάφους (ground-track
% speed) μιας κυκλικής τροχιάς στο υψόμετρο του δορυφόρου (Κεπλεριανή
% περίοδος, χωρίς αφαίρεση της περιστροφής της Γης). Δεν είναι πλήρης
% ορβιτογράφος (SGP4 κτλ.) - αρκεί όμως για να παραχθεί ρεαλιστική
% χρονική μεταβολή elevation/SNR κατά τη διάρκεια ενός πάσσου.
%
% Σε κάθε χρονικό βήμα ξαναδιαλέγεται LOS + shadow fading (νέο RNG seed),
% όπως θα συνέβαινε σε διαδοχικές πραγματικές μεταδόσεις.
%
% Χρήση:
%   T = temporalPassSimulation();          % dt = 5s -> ../Results
%   T = temporalPassSimulation(2);         % dt = 2s
%   T = temporalPassSimulation(5, 'C:\out')

if nargin < 1 || isempty(dtSeconds)
    dtSeconds = 5;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
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
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.TxPower = 43;
simParameters.RxNoiseFigure = 5;
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
bs_geo(:,3) = 25;

satAltitude = 550e3;
satParameters.CarrierFrequency = 2.01e9;
satParameters.TxPower = 34;
satParameters.AntennaGain = 30;
satParameters.EIRP = satParameters.TxPower + satParameters.AntennaGain;
satParameters.Bandwidth = 20e6;
satParameters.MinElevationDeg = 10;

simParameters.Power.NumTrx = 1;
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.Power.Pfix  = 0;
satParameters.Power.EtaPA = 0.4;

%% ------------------ Ground track του LEO (απλοποιημένο μοντέλο πάσσου) ------------------
muEarth = 3.986004418e14;  % m^3/s^2, βαρυτική παράμετρος Γης
Re      = 6371e3;          % m, μέση ακτίνα Γης
a       = Re + satAltitude;
orbitalPeriodS  = 2*pi*sqrt(a^3/muEarth);      % Κεπλεριανή περίοδος (s)
groundSpeedMps  = (2*pi/orbitalPeriodS) * Re;  % ταχύτητα ίχνους εδάφους (m/s)

centerLat = mean(bs_geo(:,1));
centerLon = mean(bs_geo(:,2));

startOffsetKm = -1500;  % km ανατολικά του κέντρου, αρχή του πάσσου (δορυφόρος αόρατος)
maxSteps      = 2000;   % ασφαλιστικό όριο βημάτων

%% ------------------ Χρονικός βρόχος ------------------
allRows = cell(maxSteps,1);
wasVisible = false;
step = 0;
t = 0;

while step < maxSteps
    offsetKm = startOffsetKm + groundSpeedMps * t / 1000;
    subLon = centerLon + offsetKm / (111.320*cosd(centerLat));
    sat_geo = [centerLat, subLon, satAltitude];

    rng(step + 1); % νέο LOS/shadow-fading draw ανά χρονικό βήμα (νέα μετάδοση)

    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

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
        break; % ο πάσσος έληξε (ο δορυφόρος έγινε αόρατος σε όλους μετά από ορατότητα)
    end

    t = t + dtSeconds;
end

T = vertcat(allRows{1:step});
writetable(T, fullfile(outputDir, 'temporal_pass_dataset.csv'));

fprintf('Temporal pass simulation: %d χρονικά βήματα (dt=%ds, ~%.0f s συνολικά), ταχύτητα ίχνους %.0f m/s (T_orbit=%.0f s)\n', ...
    step, dtSeconds, T.Time_s(end), groundSpeedMps, orbitalPeriodS);

%% ------------------ Handovers ανά χρήστη ------------------
fprintf('\n--- Handovers ανά χρήστη (αλλαγές ServingNode μεταξύ διαδοχικών βημάτων) ---\n');
handoverCounts = zeros(numUsers,1);
for u = 1:numUsers
    userRows = T(T.UserID == u, :);
    userRows = sortrows(userRows, 'Step');
    changes = sum(userRows.ServingNode(2:end) ~= userRows.ServingNode(1:end-1));
    handoverCounts(u) = changes;
    fprintf('  User %d: %d handovers (%s -> ... -> %s)\n', u, changes, ...
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
    title(sprintf('%s κατά τη διάρκεια πάσσου LEO', kpiList{k}), 'Interpreter','none');
    saveas(fig, fullfile(outputDir, ['temporal_' kpiList{k} '.png']));
    close(fig);
end

fprintf('\nTemporal pass simulation results -> %s\n', outputDir);

end
