function visual(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec, bestNodeVec, userRingM)
    %% ------------------ 3D Visualization + κάτοψη εδάφους ------------------
    % Δύο πλαίσια: αριστερά η 3D εικόνα με τους δορυφόρους σε πραγματική
    % κλίμακα, όπου όμως το επίγειο δίκτυο (~1 km) χάνεται μέσα σε άξονες
    % χιλιάδων km· δεξιά η κάτοψη του εδάφους σε km, όπου φαίνονται οι θέσεις
    % των σταθμών, των χρηστών και οι επίγειες ζεύξεις.
    % userRingM = [rMin rMax]: όρια απόστασης χρήστη από τον πλησιέστερο
    % σταθμό, σχεδιάζονται στην κάτοψη (προαιρετικό).
    if nargin < 9
        userRingM = [];
    end
    fig = figure('Name', '3D Terrestrial & NTN Network', 'Color', 'w', 'Position', [100, 100, 1500, 700]);
    try
        theme(fig, 'light');
    catch
    end
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Σημείο αναφοράς: το κέντρο του πλέγματος των σταθμών
    lat0 = mean(bs_geo(:,1));
    lon0 = mean(bs_geo(:,2));
    h0   = 0;

    % Μετατροπή μόνο X, Y (αγνοούμε το Z του ENU λόγω καμπυλότητας της Γης)
    [xBS, yBS, ~] = geodetic2enu(bs_geo(:,1), bs_geo(:,2), bs_geo(:,3), lat0, lon0, h0, wgs84);
    [xUE, yUE, ~] = geodetic2enu(user_geo(:,1), user_geo(:,2), user_geo(:,3), lat0, lon0, h0, wgs84);
    [xSat, ySat, ~] = geodetic2enu(sat_geo(:,1), sat_geo(:,2), sat_geo(:,3), lat0, lon0, h0, wgs84);

    zBS  = bs_geo(:,3);
    zUE  = user_geo(:,3);
    zSat = sat_geo(:,3);

    %% ---- Αριστερά: 3D με λογαριθμικό άξονα ύψους ----
    nexttile(tl);
    hold on; grid on;

    % Σχεδίαση Base Stations (Κόκκινα τρίγωνα)· μία ετικέτα για όλο το πλέγμα,
    % αφού σε αυτή την κλίμακα οι σταθμοί συμπίπτουν
    scatter3(xBS, yBS, zBS, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Base Stations (TN)');
    text(xBS(1), yBS(1), zBS(1)*1.3, sprintf(' %d BS', numBs), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

    % Σχεδίαση Χρηστών (Μπλε κύκλοι)
    scatter3(xUE, yUE, zUE, 80, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Users (UE)');
    text(xUE(1), yUE(1), zUE(1)*1.5, sprintf(' %d UE', numUsers), 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');

    % Σχεδίαση Δορυφόρων (Κίτρινα αστέρια), ένα ανά γραμμή του sat_geo
    scatter3(xSat, ySat, zSat, 300, 'p', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'DisplayName', 'Satellites (Real Scale)');
    for k = 1:numel(xSat)
        text(xSat(k), ySat(k), zSat(k)*1.2, " SAT-" + k, 'FontSize', 11, 'FontWeight', 'bold');
    end

    % Σχεδίαση Γραμμών Σύνδεσης (Με Interpolation για να έχουμε λογαριθμική κλίμακα και να χωράνε όλα στο γράφιμα)
    for u = 1:numUsers
        if bestNodeTypeVec(u) == "Terrestrial"
            bs_idx = str2double(extractAfter(bestNodeVec(u), "BS"));

            num_pts = 50;
            xq = linspace(xUE(u), xBS(bs_idx), num_pts);
            yq = linspace(yUE(u), yBS(bs_idx), num_pts);
            zq = linspace(zUE(u), zBS(bs_idx), num_pts);

            plot3(xq, yq, zq, 'g-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        elseif bestNodeTypeVec(u) == "Satellite"
            sat_idx = str2double(extractAfter(bestNodeVec(u), "SAT-"));
            num_pts = 100;
            xq = linspace(xUE(u), xSat(sat_idx), num_pts);
            yq = linspace(yUE(u), ySat(sat_idx), num_pts);

            % Logarithmic interpolation για τον άξονα Z
            zq = logspace(log10(zUE(u)), log10(zSat(sat_idx)), num_pts);

            plot3(xq, yq, zq, 'm-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        elseif bestNodeTypeVec(u) == "Outage"
            % Outage: κόκκινο 'x' πάνω στον χρήστη, χωρίς γραμμή σύνδεσης.
            scatter3(xUE(u), yUE(u), zUE(u), 200, 'x', 'MarkerEdgeColor', 'r', ...
                'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
    end

    % Dummy plots για το Legend
    plot3(nan, nan, nan, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Terrestrial Link (Green)');
    plot3(nan, nan, nan, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Satellite Link (Magenta)');
    scatter3(nan, nan, nan, 200, 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 2.5, 'DisplayName', 'Outage (no usable SINR)');

    % --- Μορφοποίηση Γραφήματος και Αξόνων ---
    xlabel('East (meters)');
    ylabel('North (meters)');
    zlabel('Altitude (meters) - Log Scale');
    title({'3D: real altitudes, logarithmic z-axis', ...
        'the ground network (~1 km) collapses to a point at this scale'});

    % Εφαρμογή Λογαριθμικής Κλίμακας MONO στον άξονα Z
    set(gca, 'ZScale', 'log');

    % Ορίζουμε το Z-axis να ξεκινάει κάτω από το 1.5m του χρήστη για να μη "κοπεί" το marker
    zlim([0.5, 10^6]);

    legend('Location', 'southoutside', 'NumColumns', 3);

    % Έξυπνη γωνία θέασης
    az = -25;
    el = 12;
    view(az, el);
    grid minor;

    %% ---- Δεξιά: κάτοψη εδάφους σε km ----
    nexttile(tl);
    hold on; grid on; axis equal;
    km = 1e-3;

    % Όρια τοποθέτησης χρηστών: ισοϋψείς της απόστασης από τον πλησιέστερο
    % σταθμό (επίπεδη προσέγγιση, μόνο για την εικόνα)
    if numel(userRingM) == 2
        ext = userRingM(2) + 500;
        [gx, gy] = meshgrid(linspace(-ext, ext, 400));
        dNear = inf(size(gx));
        for b = 1:numBs
            dNear = min(dNear, hypot(gx - xBS(b), gy - yBS(b)));
        end
        contour(gx*km, gy*km, dNear, userRingM, 'k--', 'LineWidth', 1, ...
            'HandleVisibility', 'off');
        plot(nan, nan, 'k--', 'LineWidth', 1, ...
            'DisplayName', sprintf('user area: %.0f m – %.2f km from nearest BS', ...
            userRingM(1), userRingM(2)*km));
    end

    for u = 1:numUsers
        if bestNodeTypeVec(u) == "Terrestrial"
            bs_idx = str2double(extractAfter(bestNodeVec(u), "BS"));
            plot([xUE(u) xBS(bs_idx)]*km, [yUE(u) yBS(bs_idx)]*km, 'g-', ...
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

    scatter(xBS*km, yBS*km, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'base station');
    for b = 1:numBs
        text(xBS(b)*km, yBS(b)*km - 0.12, "BS" + b, 'FontSize', 8, 'FontWeight', 'bold', ...
            'Color', 'r', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    xlabel('East (km)');
    ylabel('North (km)');
    title({sprintf('Ground view: %d BS, %d users', numBs, numUsers), ...
        'number next to a magenta user = serving satellite'});
    legend('Location', 'southoutside', 'NumColumns', 2);
end
