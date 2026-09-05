function T = kpiRepeatedRuns(numRuns, outputDir, label)
%KPIREPEATEDRUNS Τρέχει την ίδια στατική τοπολογία του test_simulation.m
% (2 BS + 1 LEO, 6 χρήστες) πολλές φορές με διαφορετικό RNG seed, ώστε να
% αποτυπωθεί η διακύμανση των KPI (throughput, energy/bit, SNR) που
% οφείλεται αποκλειστικά στη στοχαστικότητα του καναλιού (LOS draw + shadow
% fading). Η γεωμετρία παραμένει σταθερή σε κάθε επανάληψη.
%
% Χρήση:
%   T = kpiRepeatedRuns();               % 500 επαναλήψεις -> ../Results
%   T = kpiRepeatedRuns(1000);           % 1000 επαναλήψεις
%   T = kpiRepeatedRuns(500, 'C:\out')
%   T = kpiRepeatedRuns(500, [], 'tag')  % tag στο όνομα του versioned φακέλου

if nargin < 1 || isempty(numRuns)
    numRuns = 500;
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

%% ------------------ Ίδια στατική τοπολογία με το test_simulation.m ------------------
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            37.0380 23.9550 1.5;
            38.0500 23.9500 1.5];

sat_geo = [38.0200 23.8200 600e3];   % LEO-600 (TR 38.821 Πίνακας 6.1.1.1-1), γεωμετρία σταθερή σε κάθε run

numUsers = size(user_geo,1);
wgs84 = wgs84Ellipsoid;

%% ------------------ Ραδιο-configuration (ίδιο με test_simulation.m) ------------------
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.TxPower = 43;
simParameters.AntennaGain = 8;   % dBi, BS antenna element gain (TR 38.901 §7.3, Table 7.3-1, G_E,max)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain;
simParameters.RxNoiseFigure = 5;
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;

scenarioType = 'UMa';
simParameters.PathLoss.Scenario = scenarioType;
bs_height_m = 25;
simParameters.PathLoss.EnvironmentHeight = 1;
bs_geo(:,3) = bs_height_m;

% Δορυφόρος: 3GPP TR 38.821 Set-1, LEO-600, S-band (Πίνακες 6.1.1.1-1 & 6.1.3.2-1).
% EIRP density (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
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

%% ------------------ Επαναλαμβανόμενα runs ------------------
allTables = cell(numRuns,1);

for r = 1:numRuns
    rng(r); % διαφορετικό seed ανά επανάληψη

    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

    runID  = repmat(r, numUsers, 1);
    userID = (1:numUsers)';

    allTables{r} = table(runID, userID, bestNodeVec, bestNodeTypeVec, ...
        bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, ...
        bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec, ...
        'VariableNames', {'RunID','UserID','ServingNode','ServingType', ...
        'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps','SatElevation_deg', ...
        'NodePower_W','EnergyPerBit_uJ'});
end

T = vertcat(allTables{:});
writetable(T, fullfile(outputDir, 'kpi_repeated_runs.csv'));

%% ------------------ Σύνοψη ανά τύπο εξυπηρέτησης ------------------
% Τυχόν γραμμές "Outage" εμφανίζονται αυτόματα ως ξεχωριστή ομάδα.
G = groupsummary(T, 'ServingType', {'mean','std','min','max'}, ...
    {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB'});
disp(G)
writetable(G, fullfile(outputDir, 'kpi_summary_by_type.csv'));

numOutage = sum(T.ServingType == "Outage");
fprintf('Outage: %d/%d γραμμές (%.2f%%)\n', numOutage, height(T), 100*numOutage/height(T));

%% ------------------ Γραφήματα (overlaid histograms ανά τύπο εξυπηρέτησης) ------------------
kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB'};
kpiLabel = {'Throughput (Mbps)','Energy per bit (\muJ/bit)','SNR (dB)'};

isTerr = T.ServingType == "Terrestrial";
isSat  = T.ServingType == "Satellite";

for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    hold on;
    vals = T.(kpiList{k});
    histogram(vals(isTerr), 'Normalization','probability', 'FaceAlpha',0.6, 'DisplayName','Terrestrial');
    histogram(vals(isSat),  'Normalization','probability', 'FaceAlpha',0.6, 'DisplayName','Satellite');
    hold off;
    xlabel(kpiLabel{k});
    ylabel('Σχετική συχνότητα');
    legend('Location','best');
    title(sprintf('%s ανά τύπο εξυπηρέτησης (%d runs, σταθερή τοπολογία 2 BS + 1 LEO)', ...
        kpiList{k}, numRuns), 'Interpreter','none');
    saveas(fig, fullfile(outputDir, ['kpi_hist_' kpiList{k} '.png']));
    close(fig);
end

fprintf('\nKPI repeated-run results (%d runs, σταθερή τοπολογία) -> %s\n', numRuns, outputDir);

%% ------------------ Versioning αποτελεσμάτων ------------------
runParams = struct();
runParams.numRuns      = numRuns;
runParams.rngScheme    = 'rng(r) ανά επανάληψη r';
runParams.scenarioType = simParameters.PathLoss.Scenario;
runParams.bs_geo       = bs_geo;
runParams.user_geo     = user_geo;
runParams.sat_geo      = sat_geo;
runParams.terrestrial  = struct('CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'TxPower_dBm', simParameters.TxPower, 'AntennaGain_dBi', simParameters.AntennaGain, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'RxAntTemperature_K', simParameters.RxAntTemperature, ...
    'NSizeGrid', simParameters.Carrier.NSizeGrid, ...
    'SubcarrierSpacing_kHz', simParameters.Carrier.SubcarrierSpacing);
runParams.terrestrial.Power = simParameters.Power;
runParams.satellite    = satParameters;

outFiles = {fullfile(outputDir, 'kpi_repeated_runs.csv'), ...
            fullfile(outputDir, 'kpi_summary_by_type.csv')};
for k = 1:numel(kpiList)
    outFiles{end+1} = fullfile(outputDir, ['kpi_hist_' kpiList{k} '.png']); %#ok<AGROW>
end
saveRunVersion('kpiRepeatedRuns', runParams, outFiles, label);

end
