function T = temporalPassSimulation(dtSeconds, outputDir)
%TEMPORALPASSSIMULATION Τρέχει το ίδιο στατικό σενάριο (2 terrestrial BS,
% 6 χρήστες) του test_simulation.m αλλά επαναλαμβανόμενα σε διαδοχικά
% χρονικά βήματα (πολλαπλές "μεταδόσεις"), μετακινώντας το υποδορυφορικό
% σημείο του LEO κατά μήκος ενός απλοποιημένου ground track ώστε να
% αποτυπωθεί η κίνηση του δορυφόρου και η μεταβολή του καναλιού
% (elevation/path loss/SNR) προς κάθε χρήστη με τον χρόνο.
%
% Απλοποίηση γεωμετρίας διέλευσης: το υποδορυφορικό σημείο κινείται σε
% σταθερό γεωγραφικό πλάτος (=πλάτος του κέντρου του BS cluster) προς
% ανατολάς, με σταθερή ταχύτητα ίση με την ταχύτητα εδάφους (ground-track
% speed) μιας κυκλικής τροχιάς στο υψόμετρο του δορυφόρου (Κεπλεριανή
% περίοδος, χωρίς αφαίρεση της περιστροφής της Γης). Δεν είναι πλήρης
% ορβιτογράφος (SGP4 κτλ.) - αρκεί όμως για να παραχθεί ρεαλιστική
% χρονική μεταβολή elevation/SNR κατά τη διάρκεια μίας διέλευσης.
%
% Η κατάσταση καναλιού (LOS/NLOS + shadow fading ανά ζεύξη) περνάει από
% βήμα σε βήμα και το νέο δείγμα σκίασης προκύπτει από το προηγούμενο μέσω
% χωρικής αυτοσυσχέτισης (Gudmundson 1991), αντί να ξαναδειγματίζεται
% ανεξάρτητα - βλ. correlatedLosState στο simulateScenario.m.
%
% Χρήση:
%   T = temporalPassSimulation();          % dt = 5s -> ../Results
%   T = temporalPassSimulation(2);         % dt = 2s
%   T = temporalPassSimulation(5, 'C:\out')

if nargin < 1 || isempty(dtSeconds)
    dtSeconds = 5;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
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
            37.0380 23.9550 1.5;
            38.0500 23.9500 1.5];

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);
wgs84 = wgs84Ellipsoid;

%% ------------------ Ραδιο-configuration (ίδιο με test_simulation.m) ------------------
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.TxPower = 43;
simParameters.AntennaGain = 8;   % dBi, BS antenna element gain (TR 38.901 §7.3, Table 7.3-1, G_E,max)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain;
simParameters.RxNoiseFigure = 5;
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
bs_geo(:,3) = 25;

satAltitude = 550e3;
satParameters.CarrierFrequency = 2.01e9;
satParameters.TxPower = 34;
satParameters.AntennaGain = 30;
satParameters.EIRP = satParameters.TxPower + satParameters.AntennaGain;
satParameters.Bandwidth = 20e6;
satParameters.MinElevationDeg = 10;

simParameters.Power.NumTrx = 1;
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.Power.Pfix  = 0;
satParameters.Power.EtaPA = 0.4;

% Υστέρηση (hysteresis) + time-to-trigger στην ενεργοποίηση/απενεργοποίηση
% ζεύξης (simulateScenario.m, κατά το πρότυπο Event A3 του 3GPP TS 38.331):
% ενεργό μόνο εδώ (temporal, sequential decisions) - δεν έχει νόημα σε
% ανεξάρτητα "drops" (test_simulation.m/monteCarloDriver.m/kpiRepeatedRuns.m).
% MarginDb=2dB είναι τυπική τιμή υστέρησης Event A3 (εύρος προτύπου 0-15dB,
% βήμα 0.5dB). TimeToTriggerSteps=1 (=1 βήμα * dt=5s = 5s επιμονή) είναι ήδη
% μεγαλύτερο από τις περισσότερες τυπικές τιμές TimeToTrigger του προτύπου
% (π.χ. 320ms) - το χρονικό βήμα dt=5s της προσομοίωσης είναι πολύ πιο αδρό
% από πραγματικά χρονικά διαστήματα σηματοδοσίας, οπότε ακόμα και 1 βήμα
% αντιπροσωπεύει ήδη μια σχετικά γενναιόδωρη περίοδο επιμονής.
simParameters.Hysteresis.MarginDb           = 2;
simParameters.Hysteresis.TimeToTriggerSteps = 1;

%% ------------------ Ground track του LEO (απλοποιημένο μοντέλο διέλευσης) ------------------
muEarth = 3.986004418e14;  % m^3/s^2, βαρυτική παράμετρος Γης
Re      = 6371e3;          % m, μέση ακτίνα Γης
a       = Re + satAltitude;
orbitalPeriodS  = 2*pi*sqrt(a^3/muEarth);      % Κεπλεριανή περίοδος (s)
groundSpeedMps  = (2*pi/orbitalPeriodS) * Re;  % ταχύτητα ίχνους εδάφους (m/s)

centerLat = mean(bs_geo(:,1));
centerLon = mean(bs_geo(:,2));

startOffsetKm = -1500;  % km ανατολικά του κέντρου, αρχή της διέλευσης (δορυφόρος αόρατος)
maxSteps      = 2000;   % ασφαλιστικό όριο βημάτων

%% ------------------ Χρονικός βρόχος ------------------
allRows = cell(maxSteps,1);
wasVisible = false;
step = 0;
t = 0;
channelState = []; % καμία προηγούμενη κατάσταση πριν το πρώτο βήμα -> πρώτο δείγμα ανεξάρτητο (i.i.d.)

while step < maxSteps
    offsetKm = startOffsetKm + groundSpeedMps * t / 1000;
    subLon = centerLon + offsetKm / (111.320*cosd(centerLat));
    sat_geo = [centerLat, subLon, satAltitude];

    rng(step + 1); % νέο RNG seed ανά χρονικό βήμα (νέα μετάδοση)

    % channelState περνάει από βήμα σε βήμα ώστε το shadow fading (και η
    % κατάσταση LOS/NLOS) κάθε ζεύξης BS-χρήστη να είναι χωρικά συσχετισμένο
    % με το προηγούμενο βήμα (Gudmundson 1991 + 3GPP TR 38.901 Πίνακας
    % 7.5-6), αντί να επαναδειγματίζεται ανεξάρτητα κάθε 5s παρόλο που οι
    % χρήστες/BS είναι ακίνητοι - βλ. correlatedLosState στο simulateScenario.m.
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
        channelState] = ...
        simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, channelState);

    step = step + 1;
    userID = (1:numUsers)';
    tCol = repmat(t, numUsers, 1);
    stepCol = repmat(step, numUsers, 1);

    allRows{step} = table(stepCol, tCol, userID, bestNodeVec, bestNodeTypeVec, ...
        bestDistanceVec, bestPathLossVec, bestSnrDbVec, capacityMbpsVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        satElevationVec, satSnrDbVec, satPathLossVec, ...
        bestBsSnrDbVec, ...
        'VariableNames', {'Step','Time_s','UserID','ServingNode','ServingType', ...
        'Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps', ...
        'NodePower_W','EnergyPerBit_uJ', ...
        'SatElevation_deg','CandSat_SNR_dB','CandSat_PathLoss_dB', ...
        'CandBS_SNR_dB'});

    refElev = max(satElevationVec);
    if refElev >= satParameters.MinElevationDeg
        wasVisible = true;
    elseif wasVisible
        break; % η διέλευση έληξε (ο δορυφόρος έγινε αόρατος σε όλους μετά από ορατότητα)
    end

    t = t + dtSeconds;
end

T = vertcat(allRows{1:step});
writetable(T, fullfile(outputDir, 'temporal_pass_dataset.csv'));

fprintf('Temporal pass simulation: %d χρονικά βήματα (dt=%ds, ~%.0f s συνολικά), ταχύτητα ίχνους %.0f m/s (T_orbit=%.0f s)\n', ...
    step, dtSeconds, T.Time_s(end), groundSpeedMps, orbitalPeriodS);

%% ------------------ Handovers, SN events & outage events ανά χρήστη ------------------
% Μετάβαση ΠΡΟΣ ServingNode="None" (outage - simulateScenario.m, ελάχιστο
% χρησιμοποιήσιμο SNR) δεν μετράει ως handover: δεν αποκτάται νέος κόμβος,
% χάνεται κάλυψη. Μετράται ξεχωριστά ως "outage event".
%
% Μεταξύ των υπόλοιπων μεταβάσεων γίνεται περαιτέρω διάκριση, με βάση την
% ορολογία master-node (MN, εδώ ο BS) / secondary-node (SN, εδώ ο
% δορυφόρος) που ήδη αναφέρεται στο Κεφ.2 (Majamaa, MR-DC):
%   - "SN event": ο δορυφόρος (SN) προστίθεται ή αφαιρείται (single- <->
%     dual-connectivity), ενώ ο BS (MN) παραμένει ο ΙΔΙΟΣ. Αντιστοιχεί στο
%     "προσθήκη/αφαίρεση δευτερεύοντος κόμβου" της βιβλιογραφίας NTN MC,
%     όχι σε πλήρες handover.
%   - "handover": αλλάζει η ταυτότητα του BS (MN) - είτε καθαρό BSx->BSy,
%     είτε πλήρης εναλλαγή Terrestrial<->Satellite, είτε ταυτόχρονη
%     αλλαγή BS ΚΑΙ SN κατάστασης.
bsPartOf  = @(s) regexprep(s, ["^None$" "^SAT-1$" "\+SAT-1$"], ["" "" ""]);
hasSatOf  = @(s) contains(s, "SAT-1");

fprintf('\n--- Handovers, SN events & outage events ανά χρήστη ---\n');
handoverCounts    = zeros(numUsers,1);
snEventCounts     = zeros(numUsers,1);
outageEventCounts = zeros(numUsers,1);
for u = 1:numUsers
    userRows = T(T.UserID == u, :);
    userRows = sortrows(userRows, 'Step');
    prevNodes = userRows.ServingNode(1:end-1);
    currNodes = userRows.ServingNode(2:end);
    changed = currNodes ~= prevNodes;
    intoOutage   = changed & (currNodes == "None");
    notOutageTxn = changed & ~(prevNodes == "None" | currNodes == "None");

    bsChanged  = notOutageTxn & (bsPartOf(prevNodes) ~= bsPartOf(currNodes));
    satChanged = notOutageTxn & (hasSatOf(prevNodes) ~= hasSatOf(currNodes));

    snEvent      = satChanged & ~bsChanged;
    realHandover = bsChanged;

    handoverCounts(u)    = sum(realHandover);
    snEventCounts(u)     = sum(snEvent);
    outageEventCounts(u) = sum(intoOutage);
    fprintf('  User %d: %d handovers, %d SN events, %d outage events (%s -> ... -> %s)\n', u, ...
        handoverCounts(u), snEventCounts(u), outageEventCounts(u), ...
        userRows.ServingNode(1), userRows.ServingNode(end));
end

%% ------------------ Γραφήματα χρονοσειράς ανά χρήστη ------------------
% Το SNR_dB (νικητή) εμφανίζει πλέον NaN σε κάθε βήμα όπου ο χρήστης είναι
% DualConnectivity (δεν υπάρχει ένα ενιαίο "SNR νικητή" όταν εξυπηρετείται
% ταυτόχρονα από δύο ζεύξεις - simulateScenario.m), κάτι που θα άφηνε το
% αντίστοιχο γράφημα σχεδόν κενό για τους περισσότερους χρήστες. Σχεδιάζεται
% αντ' αυτού το per-link CandBS_SNR_dB/CandSat_SNR_dB (πάντα πεπερασμένα -
% εκτός από CandSat_SNR_dB όταν ο δορυφόρος είναι εκτός ορατότητας),
% ανεξάρτητα από την τελική κατάσταση σύνδεσης.
kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','CandBS_SNR_dB','CandSat_SNR_dB','SatElevation_deg'};
kpiLabel = {'Throughput (Mbps)','Energy per bit (\muJ/bit)','Best-BS SNR (dB)','Satellite SNR (dB)','Sat. Elevation (deg)'};

for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    hold on;
    for u = 1:numUsers
        userRows = T(T.UserID == u, :);
        userRows = sortrows(userRows, 'Step');
        % -Inf (δορυφόρος εκτός ορατότητας, CandSat_SNR_dB) -> NaN για το
        % γράφημα, ώστε να αφήνει κενό στη γραμμή αντί να καταστρέφει την
        % κλίμακα του άξονα y.
        yVals = userRows.(kpiList{k});
        yVals(isinf(yVals)) = NaN;
        plot(userRows.Time_s, yVals, 'DisplayName', sprintf('User %d', u));
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

end
