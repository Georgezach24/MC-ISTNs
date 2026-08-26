function T = kpiRepeatedRuns(numRuns, outputDir)
%KPIREPEATEDRUNS Τρέχει την ίδια στατική τοπολογία του test_simulation.m
% πολλές φορές με διαφορετικό RNG seed, ώστε να αποτυπωθεί η διακύμανση
% των KPI από τη στοχαστικότητα του καναλιού (LOS/shadow fading) - η
% γεωμετρία παραμένει σταθερή.
%
% Χρήση:
%   T = kpiRepeatedRuns();          % 500 επαναλήψεις -> ../Results
%   T = kpiRepeatedRuns(1000);      % 1000 επαναλήψεις
%   T = kpiRepeatedRuns(500, 'C:\out')

if nargin < 1 || isempty(numRuns)
    numRuns = 500;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
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

sat_geo = [38.0200 23.8200 550e3];   % LEO 550km

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

%% ------------------ Επαναλαμβανόμενα runs ------------------
allTables = cell(numRuns,1);

for r = 1:numRuns
    rng(r); % διαφορετικό seed ανά run

    % Ζητείται και η πλήρης λίστα per-candidate εξόδων (CandBS_SNR_dB/
    % CandSat_SNR_dB), όχι μόνο ο νικητής - χρειάζονται παρακάτω.
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, ~, ~, ~, ~, ~, satSnrDbVec] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

    runID  = repmat(r, numUsers, 1);
    userID = (1:numUsers)';

    allTables{r} = table(runID, userID, bestNodeVec, bestNodeTypeVec, ...
        bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, ...
        bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, satSnrDbVec, ...
        'VariableNames', {'RunID','UserID','ServingNode','ServingType', ...
        'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps','SatElevation_deg', ...
        'NodePower_W','EnergyPerBit_uJ','CandBS_SNR_dB','CandSat_SNR_dB'});
end

T = vertcat(allTables{:});
writetable(T, fullfile(outputDir, 'kpi_repeated_runs.csv'));

%% ------------------ Σύνοψη ανά τύπο εξυπηρέτησης ------------------
G = groupsummary(T, 'ServingType', {'mean','std','min','max'}, ...
    {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB'});
disp(G)
writetable(G, fullfile(outputDir, 'kpi_summary_by_type.csv'));

numOutage = sum(T.ServingType == "Outage");
numDual   = sum(T.ServingType == "DualConnectivity");
fprintf('Outage: %d/%d γραμμές (%.2f%%)\n', numOutage, height(T), 100*numOutage/height(T));
fprintf('DualConnectivity: %d/%d γραμμές (%.2f%%)\n', numDual, height(T), 100*numDual/height(T));

%% ------------------ Διακύμανση SNR ανά ζεύξη (όχι ανά τύπο εξυπηρέτησης) ------------------
% Πηγή noise sigma για Model/train_model_noisy_snr.py. Χρησιμοποιεί το
% per-candidate CandBS_SNR_dB/CandSat_SNR_dB (πάντα υπολογισμένο,
% ανεξάρτητα από το ServingType), περιορισμένο σε χρησιμοποιήσιμους
% υποψηφίους (SNR >= minUsableSnrDb) ώστε μακρινοί, εκτός-εμβέλειας
% χρήστες να μη διογκώνουν τεχνητά τη std.
minSpectralEfficiency = 0.2344;
minUsableSnrDb = 10*log10(2^minSpectralEfficiency - 1);   % ≈ -7.53 dB, TS 38.214 MCS 0

bsUsableCand  = T.CandBS_SNR_dB  >= minUsableSnrDb;
satUsableCand = T.CandSat_SNR_dB >= minUsableSnrDb;
bsLinkStd  = std(T.CandBS_SNR_dB(bsUsableCand));
satLinkStd = std(T.CandSat_SNR_dB(satUsableCand));
linkSnrStd = table(["Terrestrial";"Satellite"], [bsLinkStd;satLinkStd], ...
    'VariableNames', {'LinkType','std_SNR_dB'});
disp(linkSnrStd)
writetable(linkSnrStd, fullfile(outputDir, 'kpi_link_snr_std.csv'));

%% ------------------ Γραφήματα (overlaid histograms ανά τύπο εξυπηρέτησης) ------------------
% boxplot() απαιτεί Statistics Toolbox (μη διαθέσιμο) - ιστογράμματα αντ' αυτού.
kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB'};
kpiLabel = {'Throughput (Mbps)','Energy per bit (\muJ/bit)','SNR (dB)'};

isTerr = T.ServingType == "Terrestrial";
isSat  = T.ServingType == "Satellite";
isDual = T.ServingType == "DualConnectivity";

for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    hold on;
    vals = T.(kpiList{k});
    histogram(vals(isTerr), 'Normalization','probability', 'FaceAlpha',0.6, 'DisplayName','Terrestrial');
    histogram(vals(isSat),  'Normalization','probability', 'FaceAlpha',0.6, 'DisplayName','Satellite');
    if any(isDual)
        histogram(vals(isDual), 'Normalization','probability', 'FaceAlpha',0.6, 'DisplayName','DualConnectivity');
    end
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

end
