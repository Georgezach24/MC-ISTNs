function plotNetworkSnapshot(snap, P, pngPath)
%PLOTNETWORKSNAPSHOT Στιγμιότυπο της διέλευσης αναφοράς στη μέγιστη ανύψωση:
% αριστερά ο ουρανός από το κέντρο του πλέγματος (αζιμούθιο/ανύψωση των
% δορυφόρων, μάσκα), δεξιά κάτοψη σε km με κυψέλες, τομείς και ζεύξεις.
%   snap  refSnapshot της runSimulation (user_geo, sat_geo, bestNodeVec, bestNodeTypeVec)
%   P     plotInfo: cellGeo (ένα κελί "BSc" ανά γραμμή), numCells, bsIsdM,
%         cellBoresightDeg ([] = χωρίς τομείς), satMinElevationDeg
C = thesisColors();
wgs84 = wgs84Ellipsoid;
cellGeo = P.cellGeo;
nU = size(snap.user_geo, 1);
typ = snap.bestNodeTypeVec(:);
node = snap.bestNodeVec(:);

% Τοπικό ENU με αρχή το κέντρο του πλέγματος (GCS: x ανατολικά, y βόρεια).
lat0 = mean(cellGeo(:,1));
lon0 = mean(cellGeo(:,2));
[xC, yC] = geodetic2enu(cellGeo(:,1), cellGeo(:,2), cellGeo(:,3), lat0, lon0, 0, wgs84);
[xU, yU] = geodetic2enu(snap.user_geo(:,1), snap.user_geo(:,2), snap.user_geo(:,3), lat0, lon0, 0, wgs84);
[azS, elS] = geodetic2aer(snap.sat_geo(:,1), snap.sat_geo(:,2), snap.sat_geo(:,3), lat0, lon0, 0, wgs84);
nSat = size(snap.sat_geo, 1);
km = 1e-3;

fig = thesisFigure(17, 9.5);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ---- (a) Ουρανός: r = 90° - ανύψωση, αζιμούθιο από βορρά δεξιόστροφα ----
pax = polaraxes(tl);
pax.Layout.Tile = 1;
hold(pax, 'on');
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
rMask = 90 - P.satMinElevationDeg;
th = linspace(0, 2*pi, 361);
polarplot(pax, th, rMask*ones(size(th)), '--', 'Color', C.text, 'LineWidth', 0.9);
used = false(nSat, 1);
for k = 1:nSat
    used(k) = any(node == "SAT-" + k);
end
for k = 1:nSat
    if elS(k) < 0
        continue;   % κάτω από τον ορίζοντα
    end
    face = 'w';
    if used(k), face = C.ntn; end
    polarplot(pax, deg2rad(azS(k)), 90 - elS(k), 'p', 'MarkerSize', 11, ...
        'MarkerFaceColor', face, 'MarkerEdgeColor', C.ntn, 'LineWidth', 1);
    text(pax, deg2rad(azS(k) + 9), 90 - elS(k), sprintf('%d', k), 'FontSize', 8, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', C.ntn);
end
pax.RLim = [0 90];
pax.RTick = [0 30 rMask];   % εξωτερικός κύκλος = ορίζοντας (0°)
pax.RTickLabel = {'90°', '60°', sprintf('%g°', P.satMinElevationDeg)};
pax.RAxisLocation = 135;
% Θέσεις των δορυφόρων σε κείμενο κάτω από τον ουρανό.
satTxt = strings(nSat, 1);
for k = 1:nSat
    if elS(k) < 0
        satTxt(k) = sprintf('SAT-%d: below the horizon', k);
    else
        satTxt(k) = sprintf('SAT-%d: az %.0f°, el %.0f°', k, mod(azS(k), 360), elS(k));
        if used(k)
            satTxt(k) = satTxt(k) + " (serving)";
        elseif elS(k) < P.satMinElevationDeg
            satTxt(k) = satTxt(k) + " (below mask)";
        end
    end
end
text(pax, deg2rad(180), 128, strjoin(satTxt, newline), 'FontSize', 7.5, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', C.muted);
pax.ThetaTick = 0:45:315;
pax.ThetaTickLabel = {'N','NE','E','SE','S','SW','W','NW'};
pax.FontSize = 8;
hold(pax, 'off');
title(pax, {'(a) Sky view from the network centre', ...
    sprintf('filled: serving satellite; dashed: %g° elevation mask', P.satMinElevationDeg)}, ...
    'FontWeight', 'normal', 'FontSize', 9);

%% ---- (b) Κάτοψη ----
ax = nexttile(tl, 2);
hold(ax, 'on');
[~, siteRow] = unique(round([xC yC]*1e3), 'rows', 'stable');   % ένα κελί ανά site
rHex = P.bsIsdM / sqrt(3) * km;
angHex = 30 + (0:6)' * 60;
for b = siteRow'
    plot(ax, xC(b)*km + rHex*cosd(angHex), yC(b)*km + rHex*sind(angHex), '-', ...
        'Color', [0.78 0.78 0.78], 'LineWidth', 0.6, 'HandleVisibility', 'off');
end
if ~isempty(P.cellBoresightDeg)
    lSec = P.bsIsdM / 4 * km;
    for c = 1:P.numCells
        plot(ax, xC(c)*km + [0 lSec*cosd(P.cellBoresightDeg(c))], ...
            yC(c)*km + [0 lSec*sind(P.cellBoresightDeg(c))], '-', ...
            'Color', [0.55 0.55 0.55], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
end
isTn  = typ == "Terrestrial";
isNtn = typ == "Satellite";
isOut = typ == "Outage";
for u = find(isTn)'
    c = str2double(extractAfter(node(u), "BS"));
    plot(ax, [xU(u) xC(c)]*km, [yU(u) yC(c)]*km, '-', 'Color', C.tn, 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end
hSite = plot(ax, xC(siteRow)*km, yC(siteRow)*km, '^', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none');
hTn  = plot(ax, xU(isTn)*km, yU(isTn)*km, 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', C.tn, 'MarkerEdgeColor', 'w');
hNtn = plot(ax, xU(isNtn)*km, yU(isNtn)*km, 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', C.ntn, 'MarkerEdgeColor', 'w');
for u = find(isNtn)'
    text(ax, xU(u)*km + 0.04, yU(u)*km, extractAfter(node(u), "SAT-"), 'FontSize', 7, ...
        'Color', C.ntn, 'VerticalAlignment', 'middle');
end
hOut = plot(ax, xU(isOut)*km, yU(isOut)*km, 'x', 'MarkerSize', 7, 'Color', C.text, 'LineWidth', 1.4);
hold(ax, 'off');
axis(ax, 'equal');
grid(ax, 'on');
xlabel(ax, 'East [km]');
ylabel(ax, 'North [km]');
nSites = numel(siteRow);
title(ax, {sprintf('(b) Ground view: %d sites, %d cells, %d users', nSites, P.numCells, nU), ...
    'number next to an NTN user: serving satellite'}, 'FontWeight', 'normal', 'FontSize', 9);
h = [hSite hTn hNtn];
lbl = {'gNB site (3 sectors)', sprintf('UE on TN (%d)', sum(isTn)), sprintf('UE on NTN (%d)', sum(isNtn))};
if any(isOut)
    h(end+1) = hOut;
    lbl{end+1} = sprintf('UE without usable link (%d)', sum(isOut));
end
legend(ax, h, lbl, 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);

title(tl, sprintf('Network snapshot at the peak elevation of the reference pass (ISD %g m)', P.bsIsdM), ...
    'FontWeight', 'bold', 'FontSize', 10);
saveThesisFigure(fig, pngPath);
end
