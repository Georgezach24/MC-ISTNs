function T = geometryValidation(outputDir, label)
% Έλεγχος της γεωμετρίας που τροφοδοτεί τα μοντέλα UMa/UMi.
%
% Σαρώνει την οριζόντια απόσταση κρατώντας σταθερό το ύψος του χρήστη και
% επαληθεύει τρία πράγματα:
%   1. Το ύψος που φτάνει στη nrPathLoss παραμένει 1.5 m σε κάθε απόσταση.
%   2. Η τρισδιάστατη απόσταση προκύπτει από sqrt(d2D^2 + (hBS-hUT)^2).
%   3. Ο έλεγχος πεδίου ισχύος αποκλείει τις ζεύξεις πέραν των 5 km.
%
% Συγκρίνει επίσης με τη γεωμετρία ENU που χρησιμοποιούνταν παλαιότερα, ώστε
% να ποσοτικοποιηθεί το σφάλμα ύψους λόγω καμπυλότητας της Γης.
%
% Επιστρέφει πίνακα με ένα σημείο ανά απόσταση και γράφει CSV/PNG.

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

    % --- Παλιά γεωμετρία ENU, για σύγκριση ---
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

fprintf('\n=== Έλεγχος γεωμετρίας (UMa, hBS=%g m, hUT=%g m, fc=%.2f GHz) ===\n', ...
    hBs, hUt, fcHz/1e9);
disp(T);

fprintf('1. Ύψος χρήστη σταθερό στα %.1f m σε όλες τις αποστάσεις : %s\n', ...
    hUt, boolText(heightOk));
fprintf('2. d3D = sqrt(d2D^2 + (hBS-hUT)^2)                        : %s\n', boolText(rangeOk));
fprintf('3. Πεδίο ισχύος 10 m <= d2D <= 5 km                        : %s\n', boolText(gateOk));
fprintf('   Μέγιστη απόκλιση ύψους της γεωμετρίας ENU              : %.1f m\n', ...
    max(abs(T.hUT_enu_m - hUt)));
fprintf('   Μέγιστο σφάλμα απωλειών από τη γεωμετρία ENU           : %.1f dB\n\n', ...
    max(T.Error_dB));

if ~(heightOk && rangeOk && gateOk)
    error('geometryValidation:Failed', 'Ο έλεγχος γεωμετρίας απέτυχε.');
end

%% ------------------ Γράφημα ------------------
fig = figure('Visible','off','Position',[100 100 900 380]);

subplot(1,2,1);
semilogx(T.d2D_m, T.hUT_used_m, '-o', 'LineWidth', 1.4); hold on;
semilogx(T.d2D_m, T.hUT_enu_m, '-s', 'LineWidth', 1.4);
xline(5000, '--', 'όριο 5 km');
grid on; xlabel('Οριζόντια απόσταση d_{2D} [m]'); ylabel('Ύψος χρήστη [m]');
legend({'φυσικό ύψος','κατακόρυφη ENU'}, 'Location','southwest');
title('Ύψος κεραίας χρήστη');

subplot(1,2,2);
semilogx(T.d2D_m, T.Error_dB, '-o', 'LineWidth', 1.4);
xline(5000, '--', 'όριο 5 km');
grid on; xlabel('Οριζόντια απόσταση d_{2D} [m]'); ylabel('Σφάλμα απωλειών [dB]');
title('Σφάλμα από τη γεωμετρία ENU');

csvPath = fullfile(outputDir, 'geometry_validation.csv');
pngPath = fullfile(outputDir, 'geometry_validation.png');
writetable(T, csvPath);
saveas(fig, pngPath);
close(fig);

fprintf('Αποτελέσματα -> %s\n', csvPath);

%% ------------------ Versioning ------------------
runParams = struct('rngSeed', 42, 'scenario', 'UMa', 'hBS_m', hBs, 'hUT_m', hUt, ...
    'CarrierFrequency_Hz', fcHz, 'd2D_targets_m', d2dTargets, ...
    'validity_d2D_min_m', 10, 'validity_d2D_max_m', 5000, ...
    'validity_hUT_min_m', 1.5, 'validity_hUT_max_m', 22.5);
saveRunVersion('geometryValidation', runParams, {csvPath, pngPath}, label);

end

function s = boolText(tf)
if tf
    s = 'ΟΚ';
else
    s = 'ΑΠΕΤΥΧΕ';
end
end
