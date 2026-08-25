function results = fairnessDemo(numUsers, maxUsersPerBs, outputDir)
%FAIRNESSDEMO Επιδεικνύει τη load-based εξισορρόπηση φορτίου
% (simParameters.Fairness.MaxUsersPerBs, βλ. simulateScenario.m) πάνω σε
% ένα σκόπιμα κατασκευασμένο σενάριο συμφόρησης: πολλοί χρήστες γύρω από
% έναν και μόνο σταθμό βάσης, με ΜΕΙΚΤΗ δορυφορική διαθεσιμότητα - κάποιοι
% έχουν χρησιμοποιήσιμο δορυφόρο (άρα εναλλακτική), κάποιοι όχι.
%
% Η μεικτή δορυφορική διαθεσιμότητα δεν είναι αυθαίρετη: το υποδορυφορικό
% σημείο τοποθετείται έτσι ώστε η γωνία ανύψωσης στην περιοχή του BS να
% βρίσκεται ΑΚΡΙΒΩΣ στο όριο ορατότητας (10°, satParameters.MinElevationDeg)
% - μια απόσταση/offset συντονισμένη εμπειρικά (via geodetic2aer πάνω στο
% ΙΔΙΟ μοντέλο, όχι εικασία) ώστε οι χρήστες σε δακτύλιο ακτίνας 2km γύρω
% από τον BS (ακόμα εντός εύλογης εμβέλειας UMa, ώστε να έχουν όλοι
% χρησιμοποιήσιμο BS SNR) να διασπώνται φυσικά σε ορατούς/μη-ορατούς προς
% τον δορυφόρο, ανάλογα με τη θέση τους στον δακτύλιο. Αυτό είναι
% γεωμετρικά γνήσιο (πραγματική μεταβολή γωνίας ανύψωσης ανά θέση), όχι
% επινοημένη ετικέτα "έχει/δεν έχει δορυφόρο".
%
% Τρέχει ΔΥΟ φορές πάνω στο ΙΔΙΟ σενάριο (ίδιο rng(42) πριν από κάθε
% κλήση, άρα ίδιες πραγματοποιήσεις καναλιού) - χωρίς εξισορρόπηση
% (MaxUsersPerBs=Inf) και με (MaxUsersPerBs=maxUsersPerBs) - ώστε η όποια
% διαφορά να οφείλεται αποκλειστικά στη νέα φάση εξισορρόπησης.
%
% Χρήση:
%   fairnessDemo();           % 16 χρήστες, όριο 6/BS -> ../Results
%   fairnessDemo(20, 5);

if nargin < 1 || isempty(numUsers)
    numUsers = 16;
end
if nargin < 2 || isempty(maxUsersPerBs)
    maxUsersPerBs = 6;
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

wgs84 = wgs84Ellipsoid;
bs_geo = [37.9838 23.7275 25];

% Χρήστες σε δακτύλιο ακτίνας 2km γύρω από τον BS (UMa, εντός εμβέλειας
% χρησιμοποιήσιμου SNR - βλ. hysteresisStressTest.m όπου το ίδιο μοντέλο
% δίνει μέσο SNR ≈ -7.5dB στα 2750m, άρα τα 2000m είναι άνετα εντός ορίου).
ringRadiusKm = 2;
angles = linspace(0, 360, numUsers+1);
angles(end) = [];
user_geo = zeros(numUsers,3);
for i = 1:numUsers
    dLat = (ringRadiusKm/110.574) * cosd(angles(i));
    dLon = (ringRadiusKm/(111.320*cosd(bs_geo(1)))) * sind(angles(i));
    user_geo(i,:) = [bs_geo(1)+dLat, bs_geo(2)+dLon, 1.5];
end

% Υποδορυφορικό σημείο συντονισμένο (βλ. σχόλιο κεφαλίδας) ώστε η γωνία
% ανύψωσης στον δακτύλιο των χρηστών να διασπάται γύρω από το όριο 10°.
sat_geo = [bs_geo(1) + 14.975, bs_geo(2), 550e3];

simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.TxPower = 43;
simParameters.AntennaGain = 8;
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain;
simParameters.RxNoiseFigure = 5;
simParameters.RxAntTemperature = 290;
simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.Power.NumTrx = 1;
simParameters.Power.P0 = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.CarrierFrequency = 2.01e9;
satParameters.TxPower = 34;
satParameters.AntennaGain = 30;
satParameters.EIRP = satParameters.TxPower + satParameters.AntennaGain;
satParameters.Bandwidth = 20e6;
satParameters.MinElevationDeg = 10;
satParameters.Power.Pfix = 0;
satParameters.Power.EtaPA = 0.4;

%% ------------------ Χωρίς εξισορρόπηση φορτίου ------------------
simParameters.Fairness.MaxUsersPerBs = Inf;
rng(42);
[bestNodeVec_off, bestNodeTypeVec_off, ~, ~, ~, capacityMbps_off, ~, ...
    nodePowerW_off, energyPerBit_off, bestBsSnrDbVec, ~, ~, ~, satElevationVec, ~, satSnrDbVec] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

%% ------------------ Με εξισορρόπηση φορτίου ------------------
simParameters.Fairness.MaxUsersPerBs = maxUsersPerBs;
rng(42);
[bestNodeVec_on, bestNodeTypeVec_on, ~, ~, ~, capacityMbps_on, ~, ...
    nodePowerW_on, energyPerBit_on] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters);

%% ------------------ Σύνοψη ------------------
satVisible = satElevationVec >= satParameters.MinElevationDeg;
hasAlternative = satVisible;   % "έχει δορυφορική εναλλακτική" - βλ. σχόλιο κεφαλίδας

T = table((1:numUsers)', hasAlternative, bestNodeTypeVec_off, capacityMbps_off, ...
    bestNodeTypeVec_on, capacityMbps_on, bestBsSnrDbVec, satSnrDbVec, satElevationVec, ...
    'VariableNames', {'User','HasSatAlternative','ServingType_off','Capacity_off_Mbps', ...
    'ServingType_on','Capacity_on_Mbps','CandBS_SNR_dB','CandSat_SNR_dB','SatElevation_deg'});
writetable(T, fullfile(outputDir, 'fairness_demo.csv'));
disp(T)

numBsUsers_off = sum(bestNodeVec_off == "BS1" | contains(bestNodeVec_off, "BS1+"));
numBsUsers_on  = sum(bestNodeVec_on  == "BS1" | contains(bestNodeVec_on,  "BS1+"));

numOffloaded = sum(bestNodeTypeVec_off == "DualConnectivity" & bestNodeTypeVec_on == "Satellite");
numProtected = sum(~hasAlternative & bestNodeTypeVec_off == "Terrestrial");
numProtectedStillOnBs = sum(~hasAlternative & bestNodeTypeVec_on == "Terrestrial");

% Δείκτης δικαιοσύνης Jain (Jain, Chiu, Hawe, 1984): J=1 -> τέλεια ίση
% κατανομή, J=1/n -> όλη η χωρητικότητα σε έναν μόνο χρήστη. Υπολογίζεται
% ΜΟΝΟ πάνω στους αρχικά BS1-connected χρήστες (τον πραγματικό "διεκδικούμενο
% πόρο" που ανακατανέμει η εξισορρόπηση) - όχι σε ολόκληρο τον πληθυσμό,
% όπου οι μόνιμα εκτός κάλυψης (Outage) χρήστες θα κυριαρχούσαν στον
% δείκτη ανεξάρτητα από το αν λειτουργεί η εξισορρόπηση ή όχι.
jainIndex = @(x) (sum(x)^2) / (numel(x) * sum(x.^2));
wasOnBs1 = (bestNodeTypeVec_off == "DualConnectivity") | (bestNodeTypeVec_off == "Terrestrial");

fprintf('\n--- Φορτίο BS1 ---\n');
fprintf('Χωρίς εξισορρόπηση: %d χρήστες στο BS1\n', numBsUsers_off);
fprintf('Με εξισορρόπηση (όριο %d):    %d χρήστες στο BS1\n', maxUsersPerBs, numBsUsers_on);

fprintf('\n--- Ποιοι μετακινήθηκαν ---\n');
fprintf('DualConnectivity χρήστες που αποσυνδέθηκαν από το BS (είχαν εναλλακτική): %d\n', numOffloaded);
fprintf('Terrestrial-only χρήστες (καμία εναλλακτική) που παρέμειναν στο BS: %d/%d (πριν: %d/%d)\n', ...
    numProtectedStillOnBs, sum(~hasAlternative), numProtected, sum(~hasAlternative));

fprintf('\n--- Δικαιοσύνη (Jain fairness index πάνω στους αρχικά BS1-connected χρήστες) ---\n');
fprintf('Χωρίς εξισορρόπηση: J=%.4f (μέσο throughput=%.2f Mbps, min=%.2f Mbps)\n', ...
    jainIndex(capacityMbps_off(wasOnBs1)), mean(capacityMbps_off(wasOnBs1)), min(capacityMbps_off(wasOnBs1)));
fprintf('Με εξισορρόπηση:    J=%.4f (μέσο throughput=%.2f Mbps, min=%.2f Mbps)\n', ...
    jainIndex(capacityMbps_on(wasOnBs1)), mean(capacityMbps_on(wasOnBs1)), min(capacityMbps_on(wasOnBs1)));

results.table = T;
results.numBsUsers_off = numBsUsers_off;
results.numBsUsers_on  = numBsUsers_on;
results.numOffloaded   = numOffloaded;
results.jain_off = jainIndex(capacityMbps_off(wasOnBs1));
results.jain_on  = jainIndex(capacityMbps_on(wasOnBs1));

%% ------------------ Γράφημα ------------------
fig = figure('Visible','off', 'Position', [100 100 800 500]);
bar([capacityMbps_off, capacityMbps_on]);
xlabel('Χρήστης'); ylabel('Χωρητικότητα (Mbps)');
legend('Χωρίς εξισορρόπηση', 'Με εξισορρόπηση', 'Location', 'best');
title(sprintf('Χωρητικότητα ανά χρήστη, BS1 όριο=%d χρήστες', maxUsersPerBs));
saveas(fig, fullfile(outputDir, 'fairness_demo_capacity.png'));
close(fig);

fprintf('\nfairnessDemo results -> %s\n', outputDir);

end
