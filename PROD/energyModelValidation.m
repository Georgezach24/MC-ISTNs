function T = energyModelValidation(outputDir, label)
%ENERGYMODELVALIDATION Ενεργειακό μοντέλο ως προς το πλήθος χρηστών L ανά κόμβο.
%   Α. Ίδια απόσταση: C_u = (B/L)*SE_u και E_u = P/(B*SE_u), ανεξάρτητη του L.
%   Β. Αυξανόμενες αποστάσεις: το bit/J του δικτύου αποτυπώνει τη σύνθεση των χρηστών.

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
simParameters.Power.NumTrx = 4;
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
capErrA = nan(nA,1); eErrA = nan(nA,1);
bwBsHz = simParameters.Carrier.NSizeGrid * 12 * simParameters.Carrier.SubcarrierSpacing * 1e3;
maxSe  = 5.5547;   % MCS 28, ίδιο όριο με το simulateScenario

for k = 1:nA
    L = userCounts(k);
    rng(1);   % ίδιο κανάλι ανά εκτέλεση: μόνο το L μεταβάλλεται
    user_geo = zeros(L,3);
    for u = 1:L
        th = 2*pi*(u-1)/L;
        user_geo(u,:) = [baseLat + dRefM*cos(th)/111320, ...
                         baseLon + dRefM*sin(th)/(111320*cosd(baseLat)), 1.5];
    end
    [~,~,~,~, sinrDb, capMbps, ~, nodePowW, ePerBit, ~,~,~,~,~,~,~, ~, netE] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);
    capA(k) = mean(capMbps, 'omitnan');
    eA(k)   = mean(ePerBit(isfinite(ePerBit)));
    bitPerJouleA(k) = netE.BitPerJoule;

    % Έλεγχος ανά χρήστη (κάθε χρήστης έχει δικό του κανάλι):
    %   C_u = (B/L_s)*SE_u,  E_u = P/(B*SE_u),  L_s = εξυπηρετούμενοι.
    served = ~isnan(capMbps);
    Ls  = sum(served);
    seU = min(log2(1 + 10.^(sinrDb(served)/10)), maxSe);
    capErrA(k) = max(abs(capMbps(served)*1e6 - bwBsHz/Ls*seU) ./ (bwBsHz/Ls*seU));
    eErrA(k)   = max(abs(ePerBit(served) - nodePowW(served)./(bwBsHz*seU)*1e6) ./ ...
                     (nodePowW(served)./(bwBsHz*seU)*1e6));
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
scaleOk  = all(capErrA < 1e-9);
energyOk = all(eErrA < 1e-9);

fprintf('\n=== Energy model: dependence on the number of users ===\n\n');
fprintf('A. Identical users at %d m\n', dRefM);
fprintf('%8s %14s %14s %16s\n','L','C_u [Mbps]','E_u [uJ/bit]','network [Mbit/J]');
for k = 1:nA
    fprintf('%8d %14.3f %14.4f %16.4f\n', T.NumUsers(k), T.CapA_Mbps(k), ...
        T.EnergyA_uJ(k), T.BitPerJouleA(k)/1e6);
end
fprintf('\nB. Users at increasing distances (300 m + %d m per user)\n', distStepM);
fprintf('%8s %14s %14s %16s\n','L','C_u [Mbps]','E_u [uJ/bit]','network [Mbit/J]');
for k = 1:nA
    fprintf('%8d %14.3f %14.4f %16.4f\n', T.NumUsers(k), T.CapB_Mbps(k), ...
        T.EnergyB_uJ(k), T.BitPerJouleB(k)/1e6);
end

fprintf('\n1. Capacity per user ~ 1/L                       : %s\n', boolText(scaleOk));
fprintf('2. Energy per bit independent of L               : %s\n', boolText(energyOk));
fprintf('   (L cancels out: E = P/(B*SE))\n');
fprintf('3. Change of network bit/J in experiment B       : %.1f%%\n\n', ...
    100*(T.BitPerJouleB(end)-T.BitPerJouleB(1))/T.BitPerJouleB(1));

if ~(scaleOk && energyOk)
    error('energyModelValidation:Failed', 'Energy model check failed.');
end

%% ------------------ Γράφημα ------------------
C = thesisColors();
fig = thesisFigure(17, 6.5);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
mk = {'MarkerSize', 4};

ax = nexttile(tl);
plot(ax, T.NumUsers, T.CapA_Mbps, '-o', 'Color', C.cat(1,:), 'MarkerFaceColor', C.cat(1,:), mk{:});
grid(ax, 'on'); xlabel(ax, 'Users per cell L'); ylabel(ax, 'C_u [Mbps]');
title(ax, '(a) Capacity per user (A)', 'FontWeight', 'normal');

ax = nexttile(tl);
plot(ax, T.NumUsers, T.EnergyA_uJ, '-o', 'Color', C.cat(1,:), 'MarkerFaceColor', C.cat(1,:), mk{:});
grid(ax, 'on'); xlabel(ax, 'Users per cell L'); ylabel(ax, 'E_u [\muJ/bit]');
title(ax, '(b) Energy per bit (A)', 'FontWeight', 'normal');
ylim(ax, [0 max(T.EnergyA_uJ)*1.5]);

ax = nexttile(tl);
plot(ax, T.NumUsers, T.BitPerJouleA/1e6, '-o', 'Color', C.cat(1,:), 'MarkerFaceColor', C.cat(1,:), mk{:});
hold(ax, 'on');
plot(ax, T.NumUsers, T.BitPerJouleB/1e6, '-s', 'Color', C.cat(2,:), 'MarkerFaceColor', C.cat(2,:), mk{:});
grid(ax, 'on'); xlabel(ax, 'Users per cell L'); ylabel(ax, 'Network efficiency [Mbit/J]');
title(ax, '(c) Network energy efficiency', 'FontWeight', 'normal');
legend(ax, {'A: same distance','B: increasing distance'}, 'Location', 'southwest');
title(tl, 'Energy model check (EARTH power model, Auer et al. 2011; equal resource sharing)', ...
    'FontWeight', 'bold', 'FontSize', 10);

pngPath = fullfile(outputDir, 'energy_model_validation.png');
saveThesisFigure(fig, pngPath);

fprintf('Figure -> %s\n', pngPath);

runParams = struct('rngSeed', 1, 'scenario', 'UMa', 'userCounts', userCounts, ...
    'refDistance_m', dRefM, 'distStep_m', distStepM, ...
    'TxPower_dBm', simParameters.TxPower, 'EIRP_dBm', simParameters.EIRP, ...
    'MinElevationDeg', satParameters.MinElevationDeg);
runParams.terrestrialPower = simParameters.Power;
runParams.satellitePower   = satParameters.Power;
runParams.results = T;
saveRunVersion('energyModelValidation', runParams, {pngPath}, label);

end

function s = boolText(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end
