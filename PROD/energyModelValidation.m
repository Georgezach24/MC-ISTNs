function T = energyModelValidation(outputDir, label)
% Έλεγχος του ενεργειακού μοντέλου ως προς το πλήθος χρηστών.
%
% Δύο πειράματα, με μοναδική μεταβλητή το πλήθος χρηστών ανά κόμβο:
%
%   Α. Πανομοιότυποι χρήστες στην ίδια απόσταση. Επαληθεύεται ότι η
%      χωρητικότητα ανά χρήστη κλιμακώνεται ως 1/L ενώ η ενέργεια ανά bit
%      παραμένει σταθερή, δηλαδή ότι το L απλοποιείται αλγεβρικά:
%      E = (P/L) / ((B/L)*SE) = P/(B*SE).
%
%   Β. Χρήστες σε αυξανόμενες αποστάσεις. Δείχνει ότι ο δείκτης bit/J σε
%      επίπεδο δικτύου, σε αντίθεση με την ενέργεια ανά bit, αποτυπώνει τη
%      σύνθεση των εξυπηρετούμενων χρηστών.
%

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if nargin < 2
    label = '';
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

wgs84 = wgs84Ellipsoid;

%% ------------------ Παράμετροι (ίδιες με το runSimulation.m) ------------------
baseLat = 37.9838;
baseLon = 23.7275;
bs_geo  = [baseLat baseLon 25];

simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.AntennaGain = 8;
simParameters.NumAntennaElements = 10;
simParameters.RxNoiseFigure = 9;
simParameters.RxAntTemperature = 290;
simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);
simParameters.Power.NumTrx = 4;   % αλυσίδες πομποδέκτη ανά τομέα· P_out/αλυσίδα <= 20 W (EARTH Πίν. 2)
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.CarrierFrequency = 2.0e9;
satParameters.Bandwidth = 30e6;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + ...
                     10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 95;   % δορυφόρος μη προσβάσιμος: απομόνωση επίγειου
satParameters.Power.Pfix  = 0;
satParameters.Power.EtaPA = 0.4;

sat_geo = [baseLat baseLon 600e3];

userCounts = 1:8;
dRefM      = 800;    % απόσταση χρηστών στο πείραμα Α

%% ------------------ Πείραμα Α: πανομοιότυποι χρήστες ------------------
nA = numel(userCounts);
capA = nan(nA,1); eA = nan(nA,1); bitPerJouleA = nan(nA,1);

for k = 1:nA
    L = userCounts(k);
    rng(1);   % ίδιο κανάλι ανά εκτέλεση: μόνο το L μεταβάλλεται
    user_geo = zeros(L,3);
    for u = 1:L
        th = 2*pi*(u-1)/L;
        user_geo(u,:) = [baseLat + dRefM*cos(th)/111320, ...
                         baseLon + dRefM*sin(th)/(111320*cosd(baseLat)), 1.5];
    end
    [~,~,~,~,~, capMbps, ~, ~, ePerBit, ~,~,~,~,~,~,~, ~, netE] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);
    capA(k) = mean(capMbps);
    eA(k)   = mean(ePerBit);
    bitPerJouleA(k) = netE.BitPerJoule;
end

%% ------------------ Πείραμα Β: χρήστες σε αυξανόμενες αποστάσεις ------------------
capB = nan(nA,1); eB = nan(nA,1); bitPerJouleB = nan(nA,1);
distStepM = 450;

for k = 1:nA
    L = userCounts(k);
    rng(1);
    user_geo = zeros(L,3);
    for u = 1:L
        d = 300 + (u-1)*distStepM;
        th = 2*pi*(u-1)/max(L,1);
        user_geo(u,:) = [baseLat + d*cos(th)/111320, ...
                         baseLon + d*sin(th)/(111320*cosd(baseLat)), 1.5];
    end
    [~,~,~,~,~, capMbps, ~, ~, ePerBit, ~,~,~,~,~,~,~, ~, netE] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);
    capB(k) = mean(capMbps);
    eB(k)   = mean(ePerBit);
    bitPerJouleB(k) = netE.BitPerJoule;
end

T = table(userCounts(:), capA, eA, bitPerJouleA, capB, eB, bitPerJouleB, ...
    'VariableNames', {'NumUsers','CapA_Mbps','EnergyA_uJ','BitPerJouleA', ...
                      'CapB_Mbps','EnergyB_uJ','BitPerJouleB'});

%% ------------------ Έλεγχοι ------------------
% 1/L κλιμάκωση: C(L)*L σταθερό
scaled   = T.CapA_Mbps .* T.NumUsers;
scaleOk  = max(abs(scaled - scaled(1))) / scaled(1) < 0.02;
% ενέργεια ανά bit ανεξάρτητη του L
energyOk = max(abs(T.EnergyA_uJ - T.EnergyA_uJ(1))) / T.EnergyA_uJ(1) < 0.02;

fprintf('\n=== Ενεργειακό μοντέλο: εξάρτηση από το πλήθος χρηστών ===\n\n');
fprintf('Α. Πανομοιότυποι χρήστες σε %d m\n', dRefM);
fprintf('%8s %14s %14s %16s\n','L','C_u [Mbps]','E_u [uJ/bit]','δίκτυο [Mbit/J]');
for k = 1:nA
    fprintf('%8d %14.3f %14.4f %16.4f\n', T.NumUsers(k), T.CapA_Mbps(k), ...
        T.EnergyA_uJ(k), T.BitPerJouleA(k)/1e6);
end
fprintf('\nΒ. Χρήστες σε αυξανόμενες αποστάσεις (300 m + %d m ανά χρήστη)\n', distStepM);
fprintf('%8s %14s %14s %16s\n','L','C_u [Mbps]','E_u [uJ/bit]','δίκτυο [Mbit/J]');
for k = 1:nA
    fprintf('%8d %14.3f %14.4f %16.4f\n', T.NumUsers(k), T.CapB_Mbps(k), ...
        T.EnergyB_uJ(k), T.BitPerJouleB(k)/1e6);
end

fprintf('\n1. Χωρητικότητα ανά χρήστη ~ 1/L                 : %s\n', boolText(scaleOk));
fprintf('2. Ενέργεια ανά bit ανεξάρτητη του L             : %s\n', boolText(energyOk));
fprintf('   (το L απλοποιείται: E = P/(B*SE), βλ. κεφ. 3)\n');
fprintf('3. Μεταβολή bit/J δικτύου στο πείραμα Β          : %.1f%%\n\n', ...
    100*(T.BitPerJouleB(end)-T.BitPerJouleB(1))/T.BitPerJouleB(1));

if ~(scaleOk && energyOk)
    error('energyModelValidation:Failed', 'Ο έλεγχος του ενεργειακού μοντέλου απέτυχε.');
end

%% ------------------ Γράφημα ------------------
fig = figure('Visible','off','Position',[100 100 950 340]);

subplot(1,3,1);
plot(T.NumUsers, T.CapA_Mbps, '-o', 'LineWidth', 1.4); grid on;
xlabel('Πλήθος χρηστών L'); ylabel('C_u [Mbps]'); title('Χωρητικότητα ανά χρήστη');

subplot(1,3,2);
plot(T.NumUsers, T.EnergyA_uJ, '-o', 'LineWidth', 1.4); grid on;
xlabel('Πλήθος χρηστών L'); ylabel('E_u [\muJ/bit]'); title('Ενέργεια ανά bit');
ylim([0 max(T.EnergyA_uJ)*1.5]);

subplot(1,3,3);
plot(T.NumUsers, T.BitPerJouleA/1e6, '-o', 'LineWidth', 1.4); hold on;
plot(T.NumUsers, T.BitPerJouleB/1e6, '-s', 'LineWidth', 1.4); grid on;
xlabel('Πλήθος χρηστών L'); ylabel('Mbit/J'); title('Απόδοση δικτύου');
legend({'ίδια απόσταση','αυξανόμενη απόσταση'}, 'Location','best');

csvPath = fullfile(outputDir, 'energy_model_validation.csv');
pngPath = fullfile(outputDir, 'energy_model_validation.png');
writetable(T, csvPath);
saveas(fig, pngPath);
close(fig);

fprintf('Αποτελέσματα -> %s\n', csvPath);

runParams = struct('rngSeed', 1, 'scenario', 'UMa', 'userCounts', userCounts, ...
    'refDistance_m', dRefM, 'distStep_m', distStepM, ...
    'TxPower_dBm', simParameters.TxPower, 'EIRP_dBm', simParameters.EIRP, ...
    'MinElevationDeg', satParameters.MinElevationDeg);
runParams.terrestrialPower = simParameters.Power;
runParams.satellitePower   = satParameters.Power;
saveRunVersion('energyModelValidation', runParams, {csvPath, pngPath}, label);

end

function s = boolText(tf)
if tf
    s = 'ΟΚ';
else
    s = 'ΑΠΕΤΥΧΕ';
end
end
