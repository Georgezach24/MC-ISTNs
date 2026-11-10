function plotNetworkSnapshot(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec, bestNodeVec, isdM, cellBoresightDeg)
%PLOTNETWORKSNAPSHOT Στιγμιότυπο του δικτύου: αριστερά 3D με λογαριθμικό
% άξονα ύψους (δορυφόροι σε πραγματική κλίμακα), δεξιά κάτοψη σε km.
%   bs_geo            μία γραμμή ανά κελί ("BSc")
%   isdM              απόσταση σταθμών (για τα εξάγωνα των κυψελών)
%   cellBoresightDeg  κατεύθυνση τομέα ανά κελί ([] = χωρίς τομείς)
fig = figure('Name', '3D Terrestrial & NTN Network', 'Color', 'w', 'Position', [100, 100, 1500, 700]);
try
    theme(fig, 'light');
catch
end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Τοπικό ENU με αρχή το κέντρο του πλέγματος.
lat0 = mean(bs_geo(:,1));
lon0 = mean(bs_geo(:,2));
h0   = 0;
[xBS, yBS, ~]   = geodetic2enu(bs_geo(:,1), bs_geo(:,2), bs_geo(:,3), lat0, lon0, h0, wgs84);
[xUE, yUE, ~]   = geodetic2enu(user_geo(:,1), user_geo(:,2), user_geo(:,3), lat0, lon0, h0, wgs84);
[xSat, ySat, ~] = geodetic2enu(sat_geo(:,1), sat_geo(:,2), sat_geo(:,3), lat0, lon0, h0, wgs84);
zBS  = bs_geo(:,3);
zUE  = user_geo(:,3);
zSat = sat_geo(:,3);

%% ---- 3D, λογαριθμικός άξονας ύψους ----
nexttile(tl);
hold on; grid on;
scatter3(xBS, yBS, zBS, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Base Stations (TN)');
text(xBS(1), yBS(1), zBS(1)*1.3, sprintf(' %d cells', numBs), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
scatter3(xUE, yUE, zUE, 80, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Users (UE)');
text(xUE(1), yUE(1), zUE(1)*1.5, sprintf(' %d UE', numUsers), 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
scatter3(xSat, ySat, zSat, 300, 'p', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'DisplayName', 'Satellites (Real Scale)');
for k = 1:numel(xSat)
    text(xSat(k), ySat(k), zSat(k)*1.2, " SAT-" + k, 'FontSize', 11, 'FontWeight', 'bold');
end

% Ζεύξεις· για τον δορυφόρο λογαριθμική παρεμβολή στο z.
for u = 1:numUsers
    if bestNodeTypeVec(u) == "Terrestrial"
        b = str2double(extractAfter(bestNodeVec(u), "BS"));
        plot3(linspace(xUE(u), xBS(b), 50), linspace(yUE(u), yBS(b), 50), ...
            linspace(zUE(u), zBS(b), 50), 'g-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    elseif bestNodeTypeVec(u) == "Satellite"
        k = str2double(extractAfter(bestNodeVec(u), "SAT-"));
        plot3(linspace(xUE(u), xSat(k), 100), linspace(yUE(u), ySat(k), 100), ...
            logspace(log10(zUE(u)), log10(zSat(k)), 100), 'm-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    elseif bestNodeTypeVec(u) == "Outage"
        scatter3(xUE(u), yUE(u), zUE(u), 200, 'x', 'MarkerEdgeColor', 'r', ...
            'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
end
plot3(nan, nan, nan, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Terrestrial Link (Green)');
plot3(nan, nan, nan, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Satellite Link (Magenta)');
scatter3(nan, nan, nan, 200, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 2.5, 'DisplayName', 'Outage (no usable SINR)');

xlabel('East (meters)');
ylabel('North (meters)');
zlabel('Altitude (meters) - Log Scale');
title({'3D: real altitudes, logarithmic z-axis', ...
    'the ground network (a few km) collapses to a point at this scale'});
set(gca, 'ZScale', 'log');
zlim([0.5, 10^6]);
legend('Location', 'southoutside', 'NumColumns', 3);
view(-25, 12);
grid minor;

%% ---- Κάτοψη σε km ----
nexttile(tl);
hold on; grid on; axis equal;
km = 1e-3;

% Κυψέλες: εξάγωνα ακτίνας ISD/sqrt(3), κορυφές σε 30° + 60°k.
[~, siteRow] = unique(round([xBS yBS]*1e3), 'rows', 'stable');   % μία γραμμή ανά site
rHex = isdM / sqrt(3) * km;
angHex = 30 + (0:6)' * 60;
for b = siteRow'
    plot(xBS(b)*km + rHex*cosd(angHex), yBS(b)*km + rHex*sind(angHex), ...
        '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
plot(nan, nan, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, ...
    'DisplayName', sprintf('cells (ISD %.0f m), where users start', isdM));

for u = 1:numUsers
    if bestNodeTypeVec(u) == "Terrestrial"
        b = str2double(extractAfter(bestNodeVec(u), "BS"));
        plot([xUE(u) xBS(b)]*km, [yUE(u) yBS(b)]*km, 'g-', ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

isTer = bestNodeTypeVec == "Terrestrial";
isSat = bestNodeTypeVec == "Satellite";
isOut = bestNodeTypeVec == "Outage";
scatter(xUE(isTer)*km, yUE(isTer)*km, 60, 'o', 'MarkerFaceColor', 'g', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'user on BS');
scatter(xUE(isSat)*km, yUE(isSat)*km, 60, 'o', 'MarkerFaceColor', 'm', ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'user on satellite');
for u = find(isSat(:))'
    text(xUE(u)*km, yUE(u)*km, "  " + extractAfter(bestNodeVec(u), "SAT-"), ...
        'FontSize', 8, 'Color', 'm');
end
scatter(xUE(isOut)*km, yUE(isOut)*km, 120, 'x', 'MarkerEdgeColor', 'r', ...
    'LineWidth', 2, 'DisplayName', 'outage');

scatter(xBS(siteRow)*km, yBS(siteRow)*km, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'base station site');
if isempty(cellBoresightDeg)
    for b = 1:numBs
        text(xBS(b)*km, yBS(b)*km - 0.12, "BS" + b, 'FontSize', 8, 'FontWeight', 'bold', ...
            'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end
else
    % Κατεύθυνση τομέα (μήκος ISD/3) με τον αριθμό του κελιού.
    lSec = isdM / 3 * km;
    for b = 1:numBs
        xe = xBS(b)*km + lSec*cosd(cellBoresightDeg(b));
        ye = yBS(b)*km + lSec*sind(cellBoresightDeg(b));
        plot([xBS(b)*km xe], [yBS(b)*km ye], 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
        text(xe, ye, "" + b, 'FontSize', 7, 'Color', 'r', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    plot(nan, nan, 'r-', 'DisplayName', 'sector boresight, cell number BSc');
end

xlabel('East (km)');
ylabel('North (km)');
title({sprintf('Ground view: %d sites, %d cells, %d users', numel(siteRow), numBs, numUsers), ...
    'number next to a magenta user = serving satellite'});
legend('Location', 'southoutside', 'NumColumns', 2);
end
