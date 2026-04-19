function visual(bs_geo, user_geo, sat_geo, wgs84, numBs, numUsers, bestNodeTypeVec,bestNodeVec)
    %% ------------------ 3D Visualization ------------------
    figure('Name', '3D Terrestrial & NTN Network', 'Color', 'w', 'Position', [100, 100, 900, 700]);
    hold on; grid on;
    
    % Ορίζουμε ως σημείο αναφοράς το 1ο BS
    lat0 = bs_geo(1,1);
    lon0 = bs_geo(1,2);
    h0   = 0;
    
    % Μετατροπή μόνο για X, Y (αγνοούμε το Z του ENU λόγω καμπυλότητας της γης κσι νσ μπορέσουμε ετσι να δούμε όλους τους χρήστες)
    [xBS, yBS, ~] = geodetic2enu(bs_geo(:,1), bs_geo(:,2), bs_geo(:,3), lat0, lon0, h0, wgs84);
    [xUE, yUE, ~] = geodetic2enu(user_geo(:,1), user_geo(:,2), user_geo(:,3), lat0, lon0, h0, wgs84);
    [xSat, ySat, ~] = geodetic2enu(sat_geo(:,1), sat_geo(:,2), sat_geo(:,3), lat0, lon0, h0, wgs84);
    
    
    zBS  = bs_geo(:,3);
    zUE  = user_geo(:,3);
    zSat = sat_geo(:,3);
    
    % Σχεδίαση Base Stations (Κόκκινα τρίγωνα)
    scatter3(xBS, yBS, zBS, 150, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Base Stations (TN)');
    for b = 1:numBs
        text(xBS(b), yBS(b), zBS(b)*1.3, " BS" + b, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r'); 
    end
    
    % Σχεδίαση Χρηστών (Μπλε κύκλοι)
    scatter3(xUE, yUE, zUE, 80, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Users (UE)');
    for u = 1:numUsers
        text(xUE(u), yUE(u), zUE(u)*1.5, " U" + u, 'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
    end
    
    % Σχεδίαση Δορυφόρου (Κίτρινο αστέρι)
    scatter3(xSat, ySat, zSat, 300, 'p', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'DisplayName', 'Satellite (Real Scale)');
    text(xSat(1), ySat(1), zSat(1)*1.2, " LEO Sat", 'FontSize', 11, 'FontWeight', 'bold');
    
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
            num_pts = 100;
            xq = linspace(xUE(u), xSat, num_pts);
            yq = linspace(yUE(u), ySat, num_pts);
            
            % Logarithmic interpolation για τον άξονα Z
            zq = logspace(log10(zUE(u)), log10(zSat(1)), num_pts);
            
            plot3(xq, yq, zq, 'm-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % Dummy plots για το Legend
    plot3(nan, nan, nan, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Terrestrial Link (Green)');
    plot3(nan, nan, nan, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Satellite Link (Magenta)');
    
    % --- Μορφοποίηση Γραφήματος και Αξόνων ---
    xlabel('East (meters)');
    ylabel('North (meters)');
    zlabel('Altitude (meters) - Log Scale');
    title('3D Network Simulation: Real Altitudes with Logarithmic Z-Axis');
    
    % Εφαρμογή Λογαριθμικής Κλίμακας MONO στον άξονα Z
    set(gca, 'ZScale', 'log');
    
    % Ορίζουμε το Z-axis να ξεκινάει κάτω από το 1.5m του χρήστη για να μη "κοπεί" το marker
    zlim([0.5, 10^6]); 
    
    legend('Location', 'northeastoutside');
    
    % Έξυπνη γωνία θέασης
    az = -25; 
    el = 12;
    view(az, el);
    grid minor;
end