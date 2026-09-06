function T = temporalPassSimulation(dtSeconds, outputDir, label)
%TEMPORALPASSSIMULATION Τρέχει το στατικό σενάριο του test_simulation.m
% (2 BS, 6 χρήστες) επαναλαμβανόμενα σε διαδοχικά χρονικά βήματα,
% μετακινώντας τον δορυφόρο κατά μήκος κυκλικής Κεπλεριανής τροχιάς. Η
% θέση διαδίδεται σε αδρανειακό σύστημα και μετατρέπεται σε γεωδαιτικές
% συντεταγμένες λαμβάνοντας υπόψη την περιστροφή της Γης, ώστε ταχύτητα και
% γεωμετρία ίχνους να προέρχονται από το ίδιο μοντέλο. Η κατάσταση καναλιού
% περνάει από βήμα σε βήμα (χωρικά συσχετισμένο shadow fading).
%
% Καταγράφονται επίσης οι μεταπομπές ως ρητές μεταβάσεις κατάστασης, με
% χρόνο διακοπής 2*RTT κατά TR 38.821 §7.3.2.1.1.
%
% Χρήση:
%   T = temporalPassSimulation();              % dt = 5s -> ../Results
%   T = temporalPassSimulation(2);             % dt = 2s
%   T = temporalPassSimulation(5, 'C:\out')
%   T = temporalPassSimulation(5, [], 'tag')   % tag στο όνομα του versioned φακέλου

if nargin < 1 || isempty(dtSeconds)
    dtSeconds = 5;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if nargin < 3
    label = '';
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% ------------------ Ίδια στατική τοπολογία (BS + χρήστες) με το test_simulation.m ------------------
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];

user_geo = [37.9845 23.7288 1.5;
            37.9870 23.7325 1.5;
            37.9825 23.7268 1.5;
            37.9900 23.7450 1.5;
            37.9997 23.7476 1.5;
            37.9484 23.7111 1.5];

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);
wgs84 = wgs84Ellipsoid;

%% ------------------ Ραδιο-configuration (ίδιο με test_simulation.m) ------------------
% Επίγειο: TR 38.901 §7.8 Πίν. 7.8-1. BS gain = element (8 dBi) + array
% 10·log10(N), N=10, ιδανική στόχευση -> composite 18 dBi.
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.AntennaGain = 8;              % dBi, element gain (TR 38.901 §7.3 Πίν. 7.3-1)
simParameters.NumAntennaElements = 10;      % TR 38.901 Πίν. 7.8-1
simParameters.RxNoiseFigure = 9;           % dB, UE downlink NF (TR 38.901 Πίν. 7.8-1)
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;                 % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMa)
bs_geo(:,3) = 25;

simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);   % dBm, composite

% Δορυφόρος: 3GPP TR 38.821 Set-1, LEO-600, S-band (Πίνακες 6.1.1.1-1 & 6.1.3.2-1).
% EIRP density (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
satAltitude = 600e3;                          % m, LEO-600
satParameters.CarrierFrequency = 2.0e9;
satParameters.Bandwidth = 30e6;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + 10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 20;

simParameters.Power.NumTrx = 4;   % αλυσίδες πομποδέκτη ανά τομέα· P_out/αλυσίδα <= 20 W (EARTH Πίν. 2)
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.Power.Pfix  = 0;       % W, εκτός ενισχυτή· Pfix=0 -> αισιόδοξη υπόθεση
satParameters.Power.EtaPA = 0.4;

%% ------------------ Τροχιά LEO (κυκλική Κεπλεριανή) ------------------
% Στοιχεία εφημερίδας κατά TR 38.821 Πίν. 7.3.6.1-1, με εκκεντρότητα μηδέν.
% Η κλίση δεν ορίζεται από το πρότυπο για LEO-600 και επιλέγεται ώστε η
% διέλευση να φτάνει σε υψηλή ανύψωση πάνω από το σημείο αναφοράς.
muEarth   = 3.986004418e14;   % m^3/s^2, βαρυτική παράμετρος Γης
Re        = 6378137;          % m, ισημερινή ακτίνα WGS84
omegaEarth= 7.2921150e-5;     % rad/s, γωνιακή ταχύτητα περιστροφής Γης (WGS84)
inclDeg   = 53;               % μοίρες, κλίση τροχιάς
a         = Re + satAltitude;                  % m, μεγάλος ημιάξονας
orbitalPeriodS = 2*pi*sqrt(a^3/muEarth);       % s, Κεπλεριανή περίοδος
meanMotion     = 2*pi/orbitalPeriodS;          % rad/s
groundSpeedMps = meanMotion * Re;              % m/s, ταχύτητα ίχνους εδάφους

centerLat = mean(bs_geo(:,1));
centerLon = mean(bs_geo(:,2));

% Όρισμα πλάτους στο σημείο μέγιστης προσέγγισης: sin(lat) = sin(i)*sin(u).
uPeakRad  = asin(min(max(sind(centerLat)/sind(inclDeg), -1), 1));
% Η διέλευση ξεκινά πριν το σημείο αυτό, εκτός ορατότητας.
leadRad   = deg2rad(30);
u0Rad     = uPeakRad - leadRad;
tPeakS    = leadRad / meanMotion;
% Ορθή αναφορά ανερχόμενου δεσμού ώστε το ίχνος να περνά από το κέντρο.
lonPeakInertialRad = atan2(cosd(inclDeg)*sin(uPeakRad), cos(uPeakRad));
raanRad   = deg2rad(centerLon) + omegaEarth*tPeakS - lonPeakInertialRad;

maxSteps  = 4000;   % ασφαλιστικό όριο βημάτων

%% ------------------ Χρονικός βρόχος ------------------
allRows = cell(maxSteps,1);
wasVisible = false;
prevServingNode = strings(numUsers,1);   % κατάσταση μεταπομπής ανά χρήστη
cLight = physconst('LightSpeed');
step = 0;
t = 0;
channelState = []; % καμία προηγούμενη κατάσταση πριν το πρώτο βήμα -> πρώτο δείγμα ανεξάρτητο (i.i.d.)

while step < maxSteps
    sat_geo = orbitPositionLla(u0Rad + meanMotion*t, inclDeg, raanRad, a, omegaEarth, t);

    rng(step + 1); % RNG seed ανά χρονικό βήμα

    % channelState περνάει από βήμα σε βήμα -> χωρικά συσχετισμένο shadow fading.
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
        channelState, ~, serviceStateVec, throughputMbpsVec] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, channelState);

    % --- Κατάσταση ζεύξης και κόστος μεταπομπής ---
    % Χρόνος διακοπής = 2*RTT για την κατερχόμενη (TR 38.821 §7.3.2.1.1).
    % Δεν περιλαμβάνει καθυστέρηση επεξεργασίας RRC ούτε επανασυντονισμό,
    % όπως δηλώνει ρητά η ίδια η αναφορά.
    linkStateVec      = repmat("Stable", numUsers, 1);
    interruptionMsVec = zeros(numUsers,1);
    for u = 1:numUsers
        curr = bestNodeVec(u);
        prev = prevServingNode(u);
        if curr == "None"
            linkStateVec(u) = "Outage";
        elseif prev ~= "" && prev ~= "None" && curr ~= prev
            linkStateVec(u) = "InTransition";
            rttS = 2 * bestDistanceVec(u) / cLight;
            interruptionMsVec(u) = 2 * rttS * 1e3;   % 2*RTT σε ms
        end
    end
    prevServingNode = bestNodeVec;

    % Η διακοπή αφαιρείται από τον χρόνο του βήματος: τα bits που χάνονται
    % δεν παραδίδονται.
    lostFraction = min(interruptionMsVec/1e3/dtSeconds, 1);
    deliveredMbpsVec = throughputMbpsVec .* (1 - lostFraction);

    step = step + 1;
    userID = (1:numUsers)';
    tCol = repmat(t, numUsers, 1);
    stepCol = repmat(step, numUsers, 1);

    allRows{step} = table(stepCol, tCol, userID, bestNodeVec, bestNodeTypeVec, ...
        bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        satElevationVec, satSnrDbVec, satPathLossVec, ...
        bestBsSnrDbVec, serviceStateVec, deliveredMbpsVec, linkStateVec, interruptionMsVec, ...
        'VariableNames', {'Step','Time_s','UserID','ServingNode','ServingType', ...
        'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps', ...
        'NodePower_W','EnergyPerBit_uJ', ...
        'SatElevation_deg','CandSat_SNR_dB','CandSat_PathLoss_dB', ...
        'CandBS_SNR_dB','ServiceState','Throughput_Mbps','LinkState','Interruption_ms'});

    refElev = max(satElevationVec);
    if refElev >= satParameters.MinElevationDeg
        wasVisible = true;
    elseif wasVisible
        break; % τέλος διέλευσης (δορυφόρος αόρατος σε όλους μετά από ορατότητα)
    end

    t = t + dtSeconds;
end

T = vertcat(allRows{1:step});
writetable(T, fullfile(outputDir, 'temporal_pass_dataset.csv'));

fprintf('Temporal pass simulation: %d χρονικά βήματα (dt=%ds, ~%.0f s συνολικά), ταχύτητα ίχνους %.0f m/s (T_orbit=%.0f s)\n', ...
    step, dtSeconds, T.Time_s(end), groundSpeedMps, orbitalPeriodS);

%% ------------------ Handovers & outage events ανά χρήστη ------------------
% Μετάβαση προς/από ServingNode="None" μετράται ως outage event, όχι handover.
fprintf('\n--- Handovers & outage events ανά χρήστη ---\n');
handoverCounts    = zeros(numUsers,1);
outageEventCounts = zeros(numUsers,1);
for u = 1:numUsers
    userRows = T(T.UserID == u, :);
    userRows = sortrows(userRows, 'Step');
    prevNodes = userRows.ServingNode(1:end-1);
    currNodes = userRows.ServingNode(2:end);
    changed = currNodes ~= prevNodes;
    intoOutage = changed & (currNodes == "None");
    realHandover = changed & ~(prevNodes == "None" | currNodes == "None");

    handoverCounts(u)    = sum(realHandover);
    outageEventCounts(u) = sum(intoOutage);
    lostMs = sum(userRows.Interruption_ms);
    fprintf('  User %d: %d handovers, %d outage events, %.0f ms diakopis (%.4f%% tou xronou) (%s -> ... -> %s)\n', ...
        u, handoverCounts(u), outageEventCounts(u), lostMs, ...
        100*lostMs/1e3/T.Time_s(end), ...
        userRows.ServingNode(1), userRows.ServingNode(end));
end
nz = T.Interruption_ms(T.Interruption_ms > 0);
if ~isempty(nz)
    fprintf('\nDiakopi ana metapombi: %.1f - %.1f ms (2*RTT, TR 38.821 7.3.2.1.1)\n', min(nz), max(nz));
end

%% ------------------ Γραφήματα χρονοσειράς ανά χρήστη ------------------
kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SNR_dB','SatElevation_deg'};
kpiLabel = {'Throughput (Mbps)','Energy per bit (\muJ/bit)','SNR (dB)','Sat. Elevation (deg)'};

for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    hold on;
    for u = 1:numUsers
        userRows = T(T.UserID == u, :);
        userRows = sortrows(userRows, 'Step');
        plot(userRows.Time_s, userRows.(kpiList{k}), 'DisplayName', sprintf('User %d', u));
    end
    hold off;
    xlabel('Χρόνος (s)');
    ylabel(kpiLabel{k});
    legend('Location','best');
    title(sprintf('%s κατά τη διάρκεια διέλευσης LEO', kpiList{k}), 'Interpreter','none');
    saveas(fig, fullfile(outputDir, ['temporal_' kpiList{k} '.png']));
    close(fig);
end

fprintf('\nTemporal pass simulation results -> %s\n', outputDir);

%% ------------------ Versioning αποτελεσμάτων ------------------
runParams = struct();
runParams.dtSeconds       = dtSeconds;
runParams.numSteps        = step;
runParams.rngScheme       = 'rng(step+1) ανά χρονικό βήμα';
runParams.inclination_deg = inclDeg;
runParams.raan_deg        = rad2deg(raanRad);
runParams.u0_deg          = rad2deg(u0Rad);
runParams.omegaEarth_rads = omegaEarth;
runParams.semiMajorAxis_m = a;
runParams.satAltitude_m   = satAltitude;
runParams.orbitalPeriod_s = orbitalPeriodS;
runParams.groundSpeed_mps = groundSpeedMps;
runParams.scenarioType    = simParameters.PathLoss.Scenario;
runParams.bs_geo          = bs_geo;
runParams.user_geo        = user_geo;
runParams.terrestrial     = struct('CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'TxPower_dBm', simParameters.TxPower, 'AntennaGain_dBi', simParameters.AntennaGain, ...
    'NumAntennaElements', simParameters.NumAntennaElements, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'RxAntTemperature_K', simParameters.RxAntTemperature);
runParams.terrestrial.Power = simParameters.Power;
runParams.satellite       = satParameters;

outFiles = {fullfile(outputDir, 'temporal_pass_dataset.csv')};
for k = 1:numel(kpiList)
    outFiles{end+1} = fullfile(outputDir, ['temporal_' kpiList{k} '.png']); %#ok<AGROW>
end
saveRunVersion('temporalPassSimulation', runParams, outFiles, label);

end

function lla = orbitPositionLla(uRad, inclDeg, raanRad, aM, omegaEarth, tS)
% Θέση δορυφόρου σε κυκλική τροχιά -> γεωδαιτικές συντεταγμένες [lat lon alt].
% Διάδοση σε αδρανειακό σύστημα και στροφή σε γεωκεντρικό-σταθερό κατά
% omegaEarth*t, ώστε θέση και ταχύτητα να προκύπτουν από το ίδιο μοντέλο.
rPerifocal = aM * [cos(uRad); sin(uRad); 0];

% Στροφή κατά την κλίση (γύρω από τον άξονα των κόμβων) και κατά τη RAAN.
i = deg2rad(inclDeg);
Rx = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
Rz = [cos(raanRad) -sin(raanRad) 0; sin(raanRad) cos(raanRad) 0; 0 0 1];
rInertial = Rz * Rx * rPerifocal;

% Αδρανειακό -> γεωκεντρικό-σταθερό: στροφή κατά τη γωνία περιστροφής της Γης.
th = omegaEarth * tS;
Rg = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
rEcef = Rg * rInertial;

lla = ecef2lla(rEcef');   % [lat lon alt] σε deg/deg/m (ελλειψοειδές WGS84)
end
