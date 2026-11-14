function files = plotResults(resultsMat, outDir)
%PLOTRESULTS Όλα τα γραφήματα αποτελεσμάτων από το simulation_results.mat,
% χωρίς νέα εκτέλεση της προσομοίωσης.
%
%   plotResults()                                    % Results/simulation_results.mat
%   plotResults('Results/simulation_results.mat', 'Results')
%
% Γραφήματα (όλα με 95% CI στις 300 διελεύσεις, M.2412-0 §7.1):
%   policy_gain.png      υβριδικό TN-NTN έναντι μόνο TN / μόνο NTN (κατά ζεύγη)
%   mlb_gain.png         συνεισφορά του MLB (TS 38.300 §15.5.1)
%   ml_gain.png          συνεισφορά της απόφασης ML (TR 37.817 §5.3)
%   handover_rates.png   HO, ping-pong, RLF, αποτυχίες HO ανά πολιτική (TR 36.839)
%   convergence.png      σύγκλιση των 95% CI
%   reference_pass.png   χρονική εξέλιξη της διέλευσης αναφοράς
%   network_snapshot.png τοπολογία στη μέγιστη ανύψωση της διέλευσης αναφοράς
if nargin < 1 || isempty(resultsMat)
    resultsMat = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', 'simulation_results.mat');
end
if nargin < 2 || isempty(outDir)
    outDir = fileparts(resultsMat);
end
R = load(resultsMat);
if ~isfield(R, 'plotInfo')
    error('plotResults:NoPlotInfo', ...
        '%s has no plotInfo; it is written by runSimulation together with the results.', resultsMat);
end
P = R.plotInfo;
pol = policyArray(R.policyTable);

files = {};
files{end+1} = plotPolicyGain(R.pairedTable, P, fullfile(outDir, 'policy_gain.png'));
files{end+1} = plotMlbGain(R.pairedTable, P, fullfile(outDir, 'mlb_gain.png'));
if any(R.pairedTable.Policy == "ActualMLBML")
    files{end+1} = plotMlGain(R.pairedTable, P, fullfile(outDir, 'ml_gain.png'));
end
files{end+1} = plotHandoverRates(pol, P, fullfile(outDir, 'handover_rates.png'));
files{end+1} = plotConvergence(R.convTable, P, fullfile(outDir, 'convergence.png'));
if ~isempty(R.refPassTable)
    files{end+1} = plotReferencePass(R.refPassTable, P, fullfile(outDir, 'reference_pass.png'));
end
if isfield(P, 'refSnapshot') && ~isempty(fieldnames(P.refSnapshot))
    f = fullfile(outDir, 'network_snapshot.png');
    plotNetworkSnapshot(P.refSnapshot, P, f);
    files{end+1} = f;
end
for k = 1:numel(files)
    fprintf('  Figure -> %s\n', files{k});
end
end

% =====================================================================
% Γραφήματα
% =====================================================================

function f = plotPolicyGain(T, P, f)
% Υβριδικό (Actual) έναντι των δικτύων ενός στρώματος: διαφορά κατά ζεύγη.
C = thesisColors();
bases = ["TerrestrialOnly","SatelliteOnly"];
relRows = {'MeanCapacity_Mbps','Mean user rate',1; 'TotalRate_Mbps','Network throughput',1; ...
           'BitPerJouleRf','Energy efficiency (RF)',1; ...
           'SpectralEff_bpsHz','Spectral efficiency',1};
ppRows  = {'FracServed','Users meeting the SE target',1; ...
           'FracOutage','Users without a usable link',-1};
spec = zeros(1,2);
for j = 1:2
    [m, ~, b] = pairStat(T, "Actual", bases(j), "Spectrum_MHz");
    spec(j) = 100*m/abs(b);
end

fig = thesisFigure(17, 11);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
tag = 'abcd';
for j = 1:2
    ax = nexttile(tl, j);
    [V, H] = gainValues(T, "Actual", bases(j), relRows(:,1), 'rel');
    forestPanel(ax, relRows(:,2), cell2mat(relRows(:,3)), V, H, C.cat(j,:), 'rel', [], j == 1);
    title(ax, sprintf('(%s) vs %s: relative change', tag(j), policyLabel(bases(j))), ...
        'FontWeight', 'normal');
    xlabel(ax, 'Change [%]');
    ax = nexttile(tl, j + 2);
    [V, H] = gainValues(T, "Actual", bases(j), ppRows(:,1), 'pp');
    forestPanel(ax, ppRows(:,2), cell2mat(ppRows(:,3)), V, H, C.cat(j,:), 'pp', [], j == 1);
    title(ax, sprintf('(%s) vs %s: share of users', tag(j+2), policyLabel(bases(j))), ...
        'FontWeight', 'normal');
    xlabel(ax, 'Change [percentage points]');
end
title(tl, 'Hybrid TN–NTN network vs single-layer networks', 'FontWeight', 'bold', 'FontSize', 10);
subtitle(tl, {sprintf('Paired difference per pass (hybrid − baseline): mean and 95%% CI over %d passes', ...
    P.numPasses), sprintf(['Available spectrum of the hybrid network: %+.0f%% vs TN only, %+.0f%% vs NTN only; ' ...
    'spectral efficiency per available Hz'], spec(1), spec(2)), legendNote()}, ...
    'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

function f = plotMlbGain(T, P, f)
% Συνεισφορά του MLB: ίδια πολιτική με και χωρίς MLB, και υβριδικό έναντι TN με MLB.
C = thesisColors();
pairs = ["ActualMLB","Actual"; "TerrestrialOnlyMLB","TerrestrialOnly"; ...
         "ActualMLB","TerrestrialOnlyMLB"];
pairLbl = {'Hybrid: with vs without MLB', 'TN only: with vs without MLB', ...
           'With MLB: hybrid vs TN only'};
relRows = {'MeanCapacity_Mbps','Mean user rate',1; 'TotalRate_Mbps','Network throughput',1; ...
           'BitPerJouleRf','Energy efficiency (RF)',1};
mobRows = {'HoRate_perUserMin','Handover rate',-1; 'HoFailRate_perUserMin','HO failure rate',-1; ...
           'RlfRate_perUserMin','RLF rate',-1};
ppRow   = {'PingPongRate_perUserMin','Ping-pong rate',-1};
usrRows = {'FracServed','Users meeting the SE target',1; 'FracSatellite','Users served by NTN',0};

fig = thesisFigure(17, 14);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {relRows,'rel','(a) Throughput and energy','Change [%]'; ...
        mobRows,'rel','(b) Mobility events','Change [%]'; ...
        ppRow,'abs','(c) Ping-pong (TR 36.839 Def. 5)','Change [events per user per minute]'; ...
        usrRows,'pp','(d) User distribution','Change [percentage points]'};
for t = 1:4
    rows = spec{t,1};
    ax = nexttile(tl, t);
    [V, H, B] = gainMatrix(T, pairs, rows(:,1), spec{t,2});
    forestPanel(ax, rows(:,2), cell2mat(rows(:,3)), V, H, C.cat, spec{t,2}, B);
    title(ax, spec{t,3}, 'FontWeight', 'normal');
    xlabel(ax, spec{t,4});
end
addLegend(tl, C.cat, pairLbl);
title(tl, 'Mobility load balancing (MLB, TS 38.300 §15.5.1)', 'FontWeight', 'bold', 'FontSize', 10);
subtitle(tl, {sprintf(['A3 cell offset = %g dB × virtual-load difference (Nasri & Altman), ' ...
    'virtual load per Awada et al.'], P.mlbMaxOffsetDb), ...
    sprintf('Paired difference per pass: mean and 95%% CI over %d passes; (×r): ratio to the baseline', ...
    P.numPasses), legendNote()}, 'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

function f = plotMlGain(T, P, f)
% Συνεισφορά της απόφασης ML: ίδια πολιτική με και χωρίς το gate.
C = thesisColors();
pairs = ["ActualML","Actual"; "ActualMLBML","ActualMLB"];
pairLbl = {'Hybrid: with vs without ML gate', 'Hybrid + MLB: with vs without ML gate'};
mobRows = {'PingPongRate_perUserMin','Ping-pong rate',-1; 'HoRate_perUserMin','Handover rate',-1; ...
           'HoRateTerrTerr_perUserMin','  of which TN–TN',-1; 'RlfRate_perUserMin','RLF rate',-1; ...
           'HoFailRate_perUserMin','HO failure rate',-1};
relRows = {'MeanCapacity_Mbps','Mean user rate',1; 'TotalRate_Mbps','Network throughput',1; ...
           'BitPerJouleRf','Energy efficiency (RF)',1};
ppRows  = {'PingPongFrac','Ping-pong share of HOs',-1; 'HoFailFrac','HO failure ratio',-1};
heldRow = {'HoHeldRate_perUserMin','HO commands held by the gate',0};

fig = thesisFigure(17, 14);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {mobRows,'rel','(a) Mobility events','Change [%]'; ...
        relRows,'rel','(b) Throughput and energy','Change [%]'; ...
        ppRows,'pp','(c) Ratios (TR 36.839 Def. 4, 5)','Change [percentage points]'; ...
        heldRow,'abs','(d) Gate activity','Events per user per minute'};
for t = 1:4
    rows = spec{t,1};
    ax = nexttile(tl, t);
    [V, H, B] = gainMatrix(T, pairs, rows(:,1), spec{t,2});
    if t == 4
        B = [];   % ο ρυθμός αναμονών της βάσης είναι 0: χωρίς λόγο
    end
    forestPanel(ax, rows(:,2), cell2mat(rows(:,3)), V, H, C.cat(1:2,:), spec{t,2}, B);
    title(ax, spec{t,3}, 'FontWeight', 'normal');
    xlabel(ax, spec{t,4});
    ax.XAxis.Exponent = 0;
end
addLegend(tl, C.cat(1:2,:), pairLbl);
title(tl, 'ML handover gate at the gNB (TS 38.300 §9.2.3.2.1, TR 37.817 §5.3)', ...
    'FontWeight', 'bold', 'FontSize', 10);
gateTxt = '';
if isfield(P, 'gate') && ~isempty(P.gate)
    gateTxt = sprintf(['Random forest, %d trees: A3-triggered HO held if P(ping-pong) > %.3f, ' ...
        're-evaluated every %d ms (reportInterval)'], P.gate.NumTrees, P.gate.Threshold, ...
        P.mlReportIntervalMs);
end
subtitle(tl, {gateTxt, sprintf('Paired difference per pass: mean and 95%% CI over %d passes', P.numPasses), ...
    legendNote()}, 'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

function f = plotHandoverRates(pol, P, f)
% Απόλυτα επίπεδα ανά πολιτική (μέσος όρος διελεύσεων, 95% CI).
C = thesisColors();
names = pol.names;
nP = numel(names);
lbl = cellfun(@policyLabel, names, 'UniformOutput', false);
typeKey = {'TerrTerr','TerrSat','SatSat'};
typeLbl = {'TN–TN','TN–NTN','NTN–NTN'};
typeCol = [C.tn; C.ntn; C.third];

fig = thesisFigure(17, 12.5);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) επιτυχή HO ανά τύπο
ax = nexttile(tl, 1); hold(ax, 'on');
S = zeros(nP, 3);
for t = 1:3
    S(:,t) = mean(pol.get(['HoRate' typeKey{t} '_perUserMin']), 1, 'omitnan')';
end
keep = any(S > 0, 1);
b = barh(ax, 1:nP, S(:,keep), 'stacked', 'BarWidth', 0.62, 'EdgeColor', 'w', 'LineWidth', 0.8);
kc = typeCol(keep,:);
for k = 1:numel(b), b(k).FaceColor = kc(k,:); end
[m, hw] = meanCi(pol.get('HoRate_perUserMin'));
errorbar(ax, m, 1:nP, hw, 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
barLabels(ax, m, hw, '%.2f');
finishBarAxes(ax, lbl, true, m + hw);
xlim(ax, [0 1.75*max(m + hw)]);
title(ax, '(a) Successful handovers', 'FontWeight', 'normal');
xlabel(ax, 'Per user per minute');
legend(ax, b, typeLbl(keep), 'Location', 'northeast', 'FontSize', 7.5);

% (b)-(d) ping-pong, RLF, αποτυχίες HO
spec = {'PingPongRate_perUserMin', 1, '%.3f', ...
        sprintf('(b) Ping-pong (return within MTS = %g s)', P.pingPongMtsS), 'Per user per minute'; ...
        'RlfRate_perUserMin', 1, '%.3f', ...
        sprintf('(c) Radio link failures (Q_{out} %g dB, T310 %g s)', P.rlfQoutDb, P.rlfT310Ms/1000), ...
        'Per user per minute'; ...
        'HoFailFrac', 100, '%.1f', '(d) Handover failure ratio (TR 36.839 Def. 4)', ...
        'Share of HO attempts [%]'};
for t = 1:3
    ax = nexttile(tl, t + 1); hold(ax, 'on');
    [m, hw] = meanCi(spec{t,2} * pol.get(spec{t,1}));
    barh(ax, 1:nP, m, 'BarWidth', 0.62, 'FaceColor', C.neutral, 'EdgeColor', 'none');
    errorbar(ax, m, 1:nP, hw, 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    barLabels(ax, m, hw, spec{t,3});
    finishBarAxes(ax, lbl, t == 2, m + hw);
    title(ax, spec{t,4}, 'FontWeight', 'normal');
    xlabel(ax, spec{t,5});
end
title(tl, 'Mobility performance per policy', 'FontWeight', 'bold', 'FontSize', 10);
subtitle(tl, sprintf(['Mean over %d passes with 95%% CI; A3 with Hys = %g dB, TTT = %d ms; ' ...
    'definitions per TR 36.839'], P.numPasses, P.hysteresisDb, P.timeToTriggerMs), ...
    'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

function f = plotConvergence(conv, P, f)
% Ημιεύρος 95% CI του μέσου όρου ως συνάρτηση του πλήθους διελεύσεων.
C = thesisColors();
n = conv.PassID;
k = @(name) find(strcmp(P.kpiNames, name), 1);
fig = thesisFigure(17, 6.5);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on');
plotCi(ax, n, 100*conv.RelCI_MeanCapacity_Mbps, C.cat(1,:), 'Mean user rate');
plotCi(ax, n, 100*conv.RelCI_MeanBitPerJouleRf, C.cat(2,:), 'Energy efficiency');
yline(ax, 100*P.ciTolerance, '--', 'Color', C.text, 'LineWidth', 0.9, 'HandleVisibility', 'off');
text(ax, (n(1) + n(end))/2, 100*P.ciTolerance, sprintf('target %g%%', 100*P.ciTolerance), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7.5, 'Color', C.muted);
ylim(ax, [0 3*100*P.ciTolerance]);
convAxes(ax, '(a) Relative half-width [%]', n);

ax = nexttile(tl); hold(ax, 'on');
plotCi(ax, n, 100*conv.RelCI_FracServed, C.cat(1,:), 'SE target met');
plotCi(ax, n, 100*conv.RelCI_FracSatellite, C.cat(2,:), 'Served by NTN');
yline(ax, 100*P.ciToleranceFrac, '--', 'Color', C.text, 'LineWidth', 0.9, 'HandleVisibility', 'off');
text(ax, (n(1) + n(end))/2, 100*P.ciToleranceFrac, sprintf('target %g pp', 100*P.ciToleranceFrac), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 7.5, 'Color', C.muted);
ylim(ax, [0 3*100*P.ciToleranceFrac]);
convAxes(ax, '(b) Absolute half-width [pp]', n);

ax = nexttile(tl); hold(ax, 'on');
hwSinr = conv.RelCI_MeanSinr_dB .* abs(conv.Mean_MeanSinr_dB);
plotCi(ax, n, hwSinr, C.cat(1,:), 'Mean SINR');
ylim(ax, [0 max(1, 1.5*hwSinr(end))]);   % λίγες διελεύσεις: μεγαλύτερο ημιεύρος
convAxes(ax, '(c) Half-width [dB], no target', n);

title(tl, 'Convergence of the 95% confidence intervals (ITU-R M.2412-0 §7.1)', ...
    'FontWeight', 'bold', 'FontSize', 10);
subtitle(tl, sprintf('Hybrid TN–NTN policy; one sample per pass (the independent unit); %d passes', n(end)), ...
    'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

function f = plotReferencePass(T, P, f)
% Χρονική εξέλιξη της διέλευσης αναφοράς (διέλευση 1, πολιτική Actual).
C = thesisColors();
T = sortrows(T, {'Step','UserID'});
nU = T.NumUsers(1);
nS = height(T) / nU;
tS = T.Time_s(1:nU:end);
M = @(col) reshape(T.(col), nU, nS);
typ = reshape(T.ServingType, nU, nS);

fig = thesisFigure(17, 17);
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on');
plot(ax, tS, median(M('CandSat_Elevation_deg'), 1, 'omitnan'), 'Color', C.ntn);
yline(ax, P.satMinElevationDeg, '--', 'Color', C.text, 'LineWidth', 0.9);
text(ax, tS(1) + 0.01*(tS(end) - tS(1)), P.satMinElevationDeg, sprintf('elevation mask %g°', P.satMinElevationDeg), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 7.5, 'Color', C.muted);
ylabel(ax, 'Elevation [°]');
title(ax, '(a) Elevation of the best visible satellite', 'FontWeight', 'normal');
timeAxes(ax, tS, false);

ax = nexttile(tl); hold(ax, 'on');
cnt = [sum(typ == "Terrestrial", 1); sum(typ == "Satellite", 1); sum(typ == "Outage", 1)]';
a = area(ax, tS, cnt, 'EdgeColor', 'none');
a(1).FaceColor = C.tn; a(2).FaceColor = C.ntn; a(3).FaceColor = C.neutral;
ylim(ax, [0 1.3*nU]);   % κενό πάνω για το υπόμνημα
yticks(ax, 0:max(1, round(nU/4)):nU);
ylabel(ax, 'Users');
legend(ax, a([2 1 3]), {'NTN','TN','no usable link'}, 'Location', 'north', ...
    'Orientation', 'horizontal', 'FontSize', 7.5);
title(ax, '(b) Users per serving layer', 'FontWeight', 'normal');
timeAxes(ax, tS, false);

ax = nexttile(tl); hold(ax, 'on');
bandPlot(ax, tS, M('SINR_dB'), C);
ylabel(ax, 'SINR [dB]');
title(ax, '(c) SINR of the serving link', 'FontWeight', 'normal');
timeAxes(ax, tS, false);

ax = nexttile(tl); hold(ax, 'on');
bandPlot(ax, tS, M('Capacity_Mbps'), C);
ylabel(ax, 'Rate [Mbps]');
title(ax, '(d) User rate', 'FontWeight', 'normal');
timeAxes(ax, tS, true);

title(tl, sprintf('Reference pass (pass 1, %d users, hybrid TN–NTN policy)', nU), ...
    'FontWeight', 'bold', 'FontSize', 10);
subtitle(tl, {sprintf('Time step %g s; line: median over users; band: 10th–90th percentile', tS(2) - tS(1)), ...
    'Users without a usable link are excluded from (c) and (d)'}, 'FontSize', 7.5, 'Color', C.muted);
saveThesisFigure(fig, f);
end

% =====================================================================
% Βοηθητικές
% =====================================================================

function pol = policyArray(policyTable)
% Πίνακας (διέλευση x πολιτική) για κάθε δείκτη του policyTable.
names = unique(policyTable.Policy, 'stable');
nPass = max(policyTable.PassID);
pol = struct();
pol.names = cellstr(names(:))';
pol.get = @(metric) reshape(policyTable.(metric), numel(names), nPass)';
assert(isequal(policyTable.Policy, repmat(names(:), nPass, 1)), ...
    'plotResults:PolicyOrder', 'policyTable is not ordered by pass and policy.');
end

function [m, hw, base] = pairStat(T, policy, baseline, metric)
r = find(T.Policy == policy & T.Baseline == baseline & T.Metric == metric, 1);
if isempty(r)
    m = NaN; hw = NaN; base = NaN;
else
    m = T.MeanDiff(r); hw = T.CI95HalfWidth(r); base = T.BaselineMean(r);
end
end

function [V, H, B] = gainValues(T, policy, baseline, metrics, kind)
% Μία στήλη: μέση διαφορά και ημιεύρος, σε % (rel), ποσοστιαίες μονάδες (pp)
% ή στη μονάδα του δείκτη (abs). B = μέση τιμή της βάσης.
n = numel(metrics);
V = nan(n,1); H = V; B = V;
for i = 1:n
    [m, hw, b] = pairStat(T, policy, baseline, string(metrics{i}));
    B(i) = b;
    switch kind
        case 'rel'
            if abs(b) > eps
                V(i) = 100*m/abs(b); H(i) = 100*hw/abs(b);
            end
        case 'pp'
            V(i) = 100*m; H(i) = 100*hw;
        otherwise
            V(i) = m; H(i) = hw;
    end
end
end

function [V, H, B] = gainMatrix(T, pairs, metrics, kind)
n = numel(metrics);
V = nan(n, size(pairs,1)); H = V; B = V;
for j = 1:size(pairs,1)
    [V(:,j), H(:,j), B(:,j)] = gainValues(T, pairs(j,1), pairs(j,2), metrics, kind);
end
end

function forestPanel(ax, labels, dirs, V, H, cols, kind, B, showY)
% Γραμμές = δείκτες, ένα σημείο ανά σύγκριση με 95% CI. Γεμάτο σημείο: το
% CI δεν περιέχει το 0. Στα 'abs' με βάση B, σε παρένθεση ο λόγος προς τη βάση.
if nargin < 8
    B = [];
end
if nargin < 9
    showY = true;
end
hold(ax, 'on');
[nR, nC] = size(V);
off = 0;
if nC > 1
    off = linspace(-0.22, 0.22, nC);
end
lo = min([0; V(:) - H(:)], [], 'omitnan');
hi = max([0; V(:) + H(:)], [], 'omitnan');
span = max(hi - lo, eps);
pad = 0.02*span;
xline(ax, 0, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.7, 'HandleVisibility', 'off');
for r = 1:nR
    for c = 1:nC
        v = V(r,c); h = H(r,c);
        if ~isfinite(v)
            continue;
        end
        y = r + off(min(c, numel(off)));
        col = cols(min(c, size(cols,1)), :);
        sig = (v - h > 0) || (v + h < 0);
        plot(ax, [v - h, v + h], [y y], '-', 'Color', col, 'LineWidth', 1.3, 'HandleVisibility', 'off');
        face = 'w';
        if sig, face = col; end
        plot(ax, v, y, 'o', 'MarkerSize', 5, 'MarkerFaceColor', face, 'MarkerEdgeColor', col, ...
            'LineWidth', 1.1, 'HandleVisibility', 'off');
        txt = fmtVal(v, kind);
        if strcmp(kind, 'abs') && ~isempty(B) && isfinite(B(r,c)) && abs(B(r,c)) > eps
            txt = sprintf('%s (×%.2g)', txt, (B(r,c) + v)/B(r,c));
        end
        if ~sig
            txt = [txt ' ns']; %#ok<AGROW>
        end
        if v >= 0
            text(ax, v + h + pad, y, txt, 'HorizontalAlignment', 'left', 'FontSize', 7.5, ...
                'Color', [0.2 0.2 0.2]);
        else
            text(ax, v - h - pad, y, txt, 'HorizontalAlignment', 'right', 'FontSize', 7.5, ...
                'Color', [0.2 0.2 0.2]);
        end
    end
end
hold(ax, 'off');
lbl = labels(:);
for r = 1:nR
    switch dirs(r)
        case 1,  lbl{r} = [lbl{r} ' \uparrow'];
        case -1, lbl{r} = [lbl{r} ' \downarrow'];
    end
end
xlim(ax, [lo - (lo < 0)*0.6*span - pad, hi + 0.6*span + pad]);
ylim(ax, [0.4, nR + 0.6]);
set(ax, 'YTick', 1:nR, 'YTickLabel', lbl, 'YDir', 'reverse', 'TickLength', [0 0]);
if ~showY
    set(ax, 'YTickLabel', {});
end
ax.XGrid = 'on'; ax.YGrid = 'off';
ax.YAxis.Color = 'none';
ax.YAxis.TickLabelColor = [0.15 0.15 0.15];
end

function s = legendNote()
s = 'Filled marker: 95% CI excludes 0 (open: ns); \uparrow / \downarrow: direction of improvement';
end

function s = fmtVal(v, kind)
a = abs(v);
if a >= 10 || a == 0
    d = 1;
elseif a >= 1
    d = 1;
elseif a >= 0.1
    d = 2;
else
    d = 3;
end
switch kind
    case 'rel', s = sprintf(['%+.' num2str(d) 'f%%'], v);
    case 'pp',  s = sprintf(['%+.' num2str(d) 'f pp'], v);
    otherwise,  s = sprintf(['%+.' num2str(max(d,3)) 'f'], v);
end
end

function addLegend(tl, cols, labels)
% Κοινό υπόμνημα στο κάτω μέρος.
ax = nexttile(tl, 1);
hold(ax, 'on');
h = gobjects(numel(labels), 1);
for j = 1:numel(labels)
    h(j) = plot(ax, NaN, NaN, '-o', 'Color', cols(j,:), 'MarkerFaceColor', cols(j,:), ...
        'MarkerSize', 5, 'LineWidth', 1.3);
end
hold(ax, 'off');
lgd = legend(ax, h, labels, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
end

function [m, hw] = meanCi(X)
m  = mean(X, 1, 'omitnan')';
hw = (1.96 * std(X, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(X), 1)))';
end

function barLabels(ax, m, hw, fmt)
top = max(m + hw, [], 'omitnan');
for i = 1:numel(m)
    if isfinite(m(i))
        text(ax, m(i) + hw(i) + 0.02*top, i, sprintf(fmt, m(i)), 'FontSize', 7.5, ...
            'VerticalAlignment', 'middle', 'Color', [0.2 0.2 0.2]);
    end
end
end

function finishBarAxes(ax, lbl, showY, ends)
top = max(ends, [], 'omitnan');
if ~(top > 0), top = 1; end
xlim(ax, [0 1.25*top]);
ylim(ax, [0.4 numel(lbl) + 0.6]);
set(ax, 'YDir', 'reverse', 'YTick', 1:numel(lbl), 'TickLength', [0 0]);
if showY
    set(ax, 'YTickLabel', lbl);
else
    set(ax, 'YTickLabel', {});
end
ax.XGrid = 'on';
ax.YAxis.Color = 'none';
ax.YAxis.TickLabelColor = [0.15 0.15 0.15];
hold(ax, 'off');
end

function plotCi(ax, n, y, col, name)
plot(ax, n, y, '-', 'Color', col, 'LineWidth', 1.2, 'DisplayName', name);
end

function convAxes(ax, ttl, n)
grid(ax, 'on');
xlim(ax, [n(1) max(n(end), n(1) + 1)]);
xlabel(ax, 'Number of passes');
title(ax, ttl, 'FontWeight', 'normal');
legend(ax, 'Location', 'northeast', 'FontSize', 7.5);
hold(ax, 'off');
end

function bandPlot(ax, t, X, C)
q = prctile(X, [10 50 90], 1);
ok = all(isfinite(q), 1);
t = t(:)'; q = q(:, ok); t = t(ok);
fill(ax, [t fliplr(t)], [q(1,:) fliplr(q(3,:))], C.band, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.8);
plot(ax, t, q(2,:), '-', 'Color', C.tn, 'LineWidth', 0.8);
end

function timeAxes(ax, t, showX)
xlim(ax, [t(1) t(end)]);
grid(ax, 'on');
if showX
    xlabel(ax, 'Time [s]');
else
    set(ax, 'XTickLabel', {});
end
hold(ax, 'off');
end

function s = policyLabel(name)
switch char(name)
    case 'Actual',             s = 'Hybrid';
    case 'TerrestrialOnly',    s = 'TN only';
    case 'SatelliteOnly',      s = 'NTN only';
    case 'ActualMLB',          s = 'Hybrid + MLB';
    case 'TerrestrialOnlyMLB', s = 'TN only + MLB';
    case 'ActualML',           s = 'Hybrid + ML gate';
    case 'ActualMLBML',        s = 'Hybrid + MLB + ML gate';
    otherwise,                 s = char(name);
end
end
