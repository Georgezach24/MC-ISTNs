function T = geometryValidation(outputDir, label)
%GEOMETRYVALIDATION Έλεγχος της γεωμετρίας των μοντέλων UMa/UMi σε σάρωση d2D:
%   1. ύψος χρήστη σταθερό 1.5 m στη nrPathLoss,
%   2. d3D = sqrt(d2D^2 + (hBS - hUT)^2),
%   3. πεδίο ισχύος 10 m - 5 km (TR 38.901 Πίν. 7.4.1-1).
% Συγκρίνει με επίπεδη γεωμετρία ENU για το σφάλμα λόγω καμπυλότητας της Γης.

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if nargin < 2
    label = '';
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

rng(42);

%% ------------------ Παράμετροι ------------------
baseLat = 37.9838;
baseLon = 23.7275;
hBs     = 25;      % m, UMa
hUt     = 1.5;     % m
fcHz    = 3.5e9;
wgs84   = wgs84Ellipsoid;

d2dTargets = [10 50 100 250 500 1000 2000 3000 4000 5000 ...
              7500 10000 20000 50000 106881];

plCfg = nrPathLossConfig;
plCfg.Scenario = 'UMa';
plCfg.EnvironmentHeight = 1;

%% ------------------ Σάρωση ------------------
n = numel(d2dTargets);
d2dActual = nan(n,1);
hUtUsed   = nan(n,1);
hUtEnu    = nan(n,1);
d3dGeom   = nan(n,1);
d3dEnu    = nan(n,1);
plCorrect = nan(n,1);
plEnu     = nan(n,1);
inRange   = false(n,1);

for k = 1:n
    % Θέση χρήστη στην επιθυμητή οριζόντια απόσταση, βόρεια του σταθμού.
    userLat = baseLat + d2dTargets(k)/111320;
    d2dActual(k) = distance(baseLat, baseLon, userLat, baseLon, wgs84);

    % --- Σωστή γεωμετρία: φυσικά ύψη, οριζόντια απόσταση ---
    hUtUsed(k) = hUt;
    d3dGeom(k) = hypot(d2dActual(k), hBs - hUt);
    txCorrect  = [0; 0; hBs];
    rxCorrect  = [d2dActual(k); 0; hUt];

    % --- Γεωμετρία ENU, για σύγκριση ---
    [~, ~, zBs] = geodetic2enu(baseLat, baseLon, hBs, baseLat, baseLon, 0, wgs84);
    [xUe, yUe, zUe] = geodetic2enu(userLat, baseLon, hUt, baseLat, baseLon, 0, wgs84);
    hUtEnu(k) = zUe;
    txEnu     = [0; 0; zBs];
    rxEnu     = [xUe; yUe; zUe];
    d3dEnu(k) = norm(rxEnu - txEnu);

    inRange(k) = d2dActual(k) >= 10 && d2dActual(k) <= 5000 && ...
                 hUt >= 1.5 && hUt <= 22.5;

    plCorrect(k) = nrPathLoss(plCfg, fcHz, false, txCorrect, rxCorrect);
    plEnu(k)     = nrPathLoss(plCfg, fcHz, false, txEnu, rxEnu);
end

T = table(d2dActual, hUtUsed, hUtEnu, d3dGeom, d3dEnu, ...
    plCorrect, plEnu, plEnu - plCorrect, inRange, ...
    'VariableNames', {'d2D_m','hUT_used_m','hUT_enu_m','d3D_geom_m','d3D_enu_m', ...
                      'PathLoss_dB','PathLoss_enu_dB','Error_dB','WithinValidity'});

%% ------------------ Ελέγχοι ------------------
tolHeight = 1e-9;
tolRange  = 1e-6;

heightOk = all(abs(T.hUT_used_m - hUt) < tolHeight);
rangeOk  = all(abs(T.d3D_geom_m - sqrt(T.d2D_m.^2 + (hBs-hUt)^2)) < tolRange);
gateOk   = isequal(T.WithinValidity, T.d2D_m >= 10 & T.d2D_m <= 5000);

fprintf('\n=== Geometry check (UMa, hBS = %g m, hUT = %g m, fc = %.2f GHz) ===\n', ...
    hBs, hUt, fcHz/1e9);
disp(T);

fprintf('1. UE height constant at %.1f m at all distances       : %s\n', ...
    hUt, boolText(heightOk));
fprintf('2. d3D = sqrt(d2D^2 + (hBS-hUT)^2)                      : %s\n', boolText(rangeOk));
fprintf('3. Validity range 10 m <= d2D <= 5 km                   : %s\n', boolText(gateOk));
fprintf('   Max height deviation of the flat ENU geometry        : %.1f m\n', ...
    max(abs(T.hUT_enu_m - hUt)));
fprintf('   Max path-loss error of the flat ENU geometry         : %.1f dB\n\n', ...
    max(T.Error_dB));

if ~(heightOk && rangeOk && gateOk)
    error('geometryValidation:Failed', 'Geometry check failed.');
end

%% ------------------ Γράφημα ------------------
C = thesisColors();
fig = thesisFigure(17, 7);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
semilogx(ax, T.d2D_m, T.hUT_used_m, '-o', 'Color', C.cat(1,:), 'MarkerSize', 4, ...
    'MarkerFaceColor', C.cat(1,:)); hold(ax, 'on');
semilogx(ax, T.d2D_m, T.hUT_enu_m, '-s', 'Color', C.cat(2,:), 'MarkerSize', 4, ...
    'MarkerFaceColor', C.cat(2,:));
xline(ax, 5000, '--', 'validity limit 5 km', 'Color', C.text, 'FontSize', 7.5, ...
    'LabelVerticalAlignment', 'middle');
grid(ax, 'on'); xlabel(ax, 'Horizontal distance d_{2D} [m]'); ylabel(ax, 'UE height [m]');
legend(ax, {'geodetic (used)','flat ENU "up" axis'}, 'Location','southwest');
title(ax, '(a) UE antenna height', 'FontWeight', 'normal');

ax = nexttile(tl);
semilogx(ax, T.d2D_m, T.Error_dB, '-o', 'Color', C.cat(2,:), 'MarkerSize', 4, ...
    'MarkerFaceColor', C.cat(2,:));
xline(ax, 5000, '--', 'validity limit 5 km', 'Color', C.text, 'FontSize', 7.5, ...
    'LabelVerticalAlignment', 'middle');
grid(ax, 'on'); xlabel(ax, 'Horizontal distance d_{2D} [m]'); ylabel(ax, 'Path-loss error [dB]');
title(ax, '(b) Path-loss error of a flat ENU geometry', 'FontWeight', 'normal');
title(tl, sprintf('Geometry check (UMa, h_{BS} = %g m, h_{UT} = %g m, TR 38.901 Table 7.4.1-1)', hBs, hUt), ...
    'FontWeight', 'bold', 'FontSize', 10);

pngPath = fullfile(outputDir, 'geometry_validation.png');
saveThesisFigure(fig, pngPath);

fprintf('Figure -> %s\n', pngPath);

%% ------------------ Versioning ------------------
runParams = struct('rngSeed', 42, 'scenario', 'UMa', 'hBS_m', hBs, 'hUT_m', hUt, ...
    'CarrierFrequency_Hz', fcHz, 'd2D_targets_m', d2dTargets, ...
    'validity_d2D_min_m', 10, 'validity_d2D_max_m', 5000, ...
    'validity_hUT_min_m', 1.5, 'validity_hUT_max_m', 22.5);
runParams.results = T;
saveRunVersion('geometryValidation', runParams, {pngPath}, label);

end

function s = boolText(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end
