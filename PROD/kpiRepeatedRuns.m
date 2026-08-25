function T = kpiRepeatedRuns(numRuns, outputDir)
%KPIREPEATEDRUNS Τρέχει την ΙΔΙΑ στατική τοπολογία του test_simulation.m
% (2 terrestrial BS + 1 LEO δορυφόρος, 6 χρήστες) πολλές φορές με
% διαφορετικό RNG seed κάθε φορά, ώστε να αποτυπωθεί η διακύμανση των KPI
% που οφείλεται αποκλειστικά στη στοχαστικότητα του καναλιού (LOS draw +
% shadow fading ανά ζεύξη, TR 38.901 §7.4.1/§7.4.2) - η γεωμετρία
% (θέσεις BS/χρηστών/δορυφόρου) παραμένει σταθερή σε κάθε run.
%
% Σκοπός: απάντηση στο αίτημα του επιβλέποντα να τρέξουμε το τρέχον setup
% (2 terrestrial BS + 1 LEO) και να δούμε τι αποτελέσματα βγάζει για
% 2-3 KPIs (throughput, energy per bit, SNR).
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

sat_geo = [38.0200 23.8200 550e3];   % LEO (550 km), γεωμετρία σταθερή σε κάθε run

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
    rng(r); % διαφορετικό seed ανά run -> διαφορετικό LOS/shadow-fading draw

    % Ζητείται η πλήρης λίστα εξόδων (όχι μόνο του νικητή) ώστε να είναι
    % διαθέσιμα τα per-candidate CandBS_SNR_dB/CandSat_SNR_dB - απαραίτητα
    % εδώ επειδή, μετά την προσθήκη DualConnectivity, δεν υπάρχει πλέον
    % εγγυημένα καμιά "Terrestrial"-only γραμμή στη σύνοψη ανά ServingType
    % (βλ. παρακάτω) από την οποία να αντληθεί η διακύμανση SNR του
    % επίγειου σκέλους.
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
% Το groupsummary ομαδοποιεί βάσει των ΠΡΑΓΜΑΤΙΚΩΝ τιμών ServingType, άρα
% τυχόν γραμμές "Outage" (simulateScenario.m - ελάχιστο χρησιμοποιήσιμο
% SNR) εμφανίζονται αυτόματα ως ξεχωριστή ομάδα, χωρίς να χρειάζεται
% ρητή διαχείριση εδώ.
G = groupsummary(T, 'ServingType', {'mean','std','min','max'}, ...
    {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB'});
disp(G)
writetable(G, fullfile(outputDir, 'kpi_summary_by_type.csv'));

numOutage = sum(T.ServingType == "Outage");
numDual   = sum(T.ServingType == "DualConnectivity");
fprintf('Outage: %d/%d γραμμές (%.2f%%)\n', numOutage, height(T), 100*numOutage/height(T));
fprintf('DualConnectivity: %d/%d γραμμές (%.2f%%)\n', numDual, height(T), 100*numDual/height(T));

%% ------------------ Διακύμανση SNR ανά ΖΕΥΞΗ (όχι ανά τύπο εξυπηρέτησης) ------------------
% Πηγή του noise sigma για το Model/train_model_noisy_snr.py. Πριν το
% DualConnectivity, η στήλη std_SNR_dB του kpi_summary_by_type.csv (groupby
% ServingType) αρκούσε, αφού κάθε νικητής κόμβος ήταν ακριβώς μία ζεύξη.
% Πλέον όμως δεν υπάρχει καμία εγγύηση ότι θα υπάρχουν "Terrestrial"-only
% γραμμές (στο τρέχον static topology ΔΕΝ υπάρχουν καθόλου - όλοι οι
% κοντινοί χρήστες παίρνουν DualConnectivity σε κάθε run, βλ.
% Ενότητα~\ref{sec:results-mc} στη διπλωματική) - το CandBS_SNR_dB/
% CandSat_SNR_dB όμως υπολογίζεται ΠΑΝΤΑ, ανεξάρτητα από την τελική
% κατάσταση σύνδεσης, οπότε η διακύμανση ανά ζεύξη αντλείται απευθείας
% από αυτά αντί από τον νικητή. Περιορίζεται στους ΧΡΗΣΙΜΟΠΟΙΗΣΙΜΟΥΣ
% υποψηφίους (SNR >= minUsableSnrDb, ίδιο κατώφλι με simulateScenario.m) -
% χωρίς αυτό τον περιορισμό, οι μακρινοί χρήστες 5-6 (πάντα εκτός εμβέλειας
% κάθε BS, SNR βαθιά αρνητικό λόγω γεωμετρίας/NLOS και όχι θορύβου
% μέτρησης) θα διόγκωναν τεχνητά τη std σε μη-ρεαλιστικά επίπεδα.
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
% boxplot() απαιτεί Statistics and Machine Learning Toolbox (μη διαθέσιμο
% εδώ) - χρησιμοποιούμε επικαλυπτόμενα ιστογράμματα (base MATLAB) αντ' αυτού.
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
