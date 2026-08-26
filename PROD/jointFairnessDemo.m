function results = jointFairnessDemo(numOverlapUsers, maxUsersPerBs, outputDir)
%JOINTFAIRNESSDEMO Επιδεικνύει τον "Joint" αλγόριθμο εξισορρόπησης φορτίου
% έναντι του "PerBs", σε σενάριο ΔΥΟ σταθμών βάσης με επικαλυπτόμενη
% κάλυψη (BS1, BS2 ~1.3km απόσταση) - το fairnessDemo.m (ένας BS) δεν
% μπορεί να θέσει υπό δοκιμή το lateral handover του Joint.
%
% 1) DC-mix: δορυφόρος σχεδόν κατακόρυφα (πάντα ορατός), όλοι οι χρήστες
%    αρχικά DualConnectivity. Δείχνει πως το PerBs αφήνει τον BS2
%    αναξιοποίητο ενώ το Joint προτιμά lateral handover σε αυτόν.
% 2) Terrestrial-only: δορυφόρος σκόπιμα απρόσιτος (MinElevationDeg=95),
%    ίδια γεωμετρία. Το PerBs δεν έχει κανέναν επιλέξιμο χρήστη
%    (περιορίζεται σε DualConnectivity) - μηδενική ανακούφιση. Το Joint
%    παραμένει σε θέση να βοηθήσει.
%
% Χρήση:
%   jointFairnessDemo();          % 14 χρήστες στη ζώνη επικάλυψης, όριο 8/BS -> ../Results
%   jointFairnessDemo(20, 10);

if nargin < 1 || isempty(numOverlapUsers)
    numOverlapUsers = 14;
end
if nargin < 2 || isempty(maxUsersPerBs)
    maxUsersPerBs = 8;
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

wgs84 = wgs84Ellipsoid;

bs1_geo = [37.9838, 23.7275, 25];
bs2LonOffsetDeg = 1.3 / (111.320 * cosd(bs1_geo(1)));   % ~1.3km ανατολικά
bs2_geo = [37.9838, 23.7275 + bs2LonOffsetDeg, 25];
bs_geo = [bs1_geo; bs2_geo];

% Ζώνη επικάλυψης: κέντρο στο 35% BS1->BS2 (πιο κοντά στον BS1, ~455m/~845m
% αντίστοιχα, άνετα εντός usable-SNR εμβέλειας), διασπορά 200m ώστε οι
% χρήστες να μην έχουν όλοι ταυτόσημο SNR.
frac = 0.35;
centerLat = bs1_geo(1);
centerLon = bs1_geo(2) + frac * bs2LonOffsetDeg;
scatterRadiusKm = 0.2;
angles = linspace(0, 360, numOverlapUsers+1);
angles(end) = [];
user_geo = zeros(numOverlapUsers,3);
for i = 1:numOverlapUsers
    dLat = (scatterRadiusKm/110.574) * cosd(angles(i));
    dLon = (scatterRadiusKm/(111.320*cosd(centerLat))) * sind(angles(i));
    user_geo(i,:) = [centerLat+dLat, centerLon+dLon, 1.5];
end

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

satParametersBase.CarrierFrequency = 2.01e9;
satParametersBase.TxPower = 34;
satParametersBase.AntennaGain = 30;
satParametersBase.EIRP = satParametersBase.TxPower + satParametersBase.AntennaGain;
satParametersBase.Bandwidth = 20e6;
satParametersBase.Power.Pfix = 0;
satParametersBase.Power.EtaPA = 0.4;

%% ================== Υπο-σενάριο 1: DC-mix ==================
% Δορυφόρος σχεδόν κατακόρυφα - ορατός παντού εδώ, καμία λεπτή ρύθμιση.
satParameters1 = satParametersBase;
satParameters1.MinElevationDeg = 10;
sat_geo1 = [centerLat, centerLon, 550e3];

simParameters.Fairness.MaxUsersPerBs = Inf;
rng(42);
[bestNodeVec_off, bestNodeTypeVec_off, ~, ~, ~, capacityMbps_off] = ...
    simulateScenario(bs_geo, user_geo, sat_geo1, wgs84, simParameters, satParameters1);

simParameters.Fairness.MaxUsersPerBs = maxUsersPerBs;
simParameters.Fairness.Joint = false;
rng(42);
[bestNodeVec_perbs, bestNodeTypeVec_perbs, ~, ~, ~, capacityMbps_perbs] = ...
    simulateScenario(bs_geo, user_geo, sat_geo1, wgs84, simParameters, satParameters1);

simParameters.Fairness.Joint = true;
rng(42);
[bestNodeVec_joint, bestNodeTypeVec_joint, ~, ~, ~, capacityMbps_joint] = ...
    simulateScenario(bs_geo, user_geo, sat_geo1, wgs84, simParameters, satParameters1);

onBs = @(v,n) sum(v == ("BS"+string(n)) | contains(v, "BS"+string(n)+"+"));
numBs1_off   = onBs(bestNodeVec_off,   1); numBs2_off   = onBs(bestNodeVec_off,   2);
numBs1_perbs = onBs(bestNodeVec_perbs, 1); numBs2_perbs = onBs(bestNodeVec_perbs, 2);
numBs1_joint = onBs(bestNodeVec_joint, 1); numBs2_joint = onBs(bestNodeVec_joint, 2);
numSat_off   = sum(bestNodeTypeVec_off   == "Satellite");
numSat_perbs = sum(bestNodeTypeVec_perbs == "Satellite");
numSat_joint = sum(bestNodeTypeVec_joint == "Satellite");

fprintf('=== Υπο-σενάριο 1: DC-mix (%d χρήστες στη ζώνη επικάλυψης, όριο %d/BS) ===\n', ...
    numOverlapUsers, maxUsersPerBs);
fprintf('%-30s %6s %6s %10s\n', '', 'BS1', 'BS2', 'Satellite-only');
fprintf('%-30s %6d %6d %10d\n', 'Χωρίς εξισορρόπηση:', numBs1_off,   numBs2_off,   numSat_off);
fprintf('%-30s %6d %6d %10d\n', 'PerBs:',              numBs1_perbs, numBs2_perbs, numSat_perbs);
fprintf('%-30s %6d %6d %10d\n', 'Joint:',               numBs1_joint, numBs2_joint, numSat_joint);
fprintf('Μέση χωρητικότητα (Mbps) - Χωρίς: %.2f | PerBs: %.2f | Joint: %.2f\n', ...
    mean(capacityMbps_off), mean(capacityMbps_perbs), mean(capacityMbps_joint));

T1 = table((1:numOverlapUsers)', bestNodeVec_off, capacityMbps_off, ...
    bestNodeVec_perbs, capacityMbps_perbs, bestNodeVec_joint, capacityMbps_joint, ...
    'VariableNames', {'User','Node_off','Capacity_off_Mbps', ...
    'Node_perbs','Capacity_perbs_Mbps','Node_joint','Capacity_joint_Mbps'});
writetable(T1, fullfile(outputDir, 'joint_fairness_dcmix.csv'));

%% ================== Υπο-σενάριο 2: Terrestrial-only ==================
% Ίδια γεωμετρία, δορυφόρος σκόπιμα απρόσιτος (MinElevationDeg=95°) - το
% PerBs δεν έχει κανέναν DualConnectivity-επιλέξιμο χρήστη.
satParameters2 = satParametersBase;
satParameters2.MinElevationDeg = 95;
sat_geo2 = sat_geo1;

simParameters.Fairness.MaxUsersPerBs = Inf;
rng(42);
[bestNodeVec_off2, bestNodeTypeVec_off2, ~, ~, ~, capacityMbps_off2] = ...
    simulateScenario(bs_geo, user_geo, sat_geo2, wgs84, simParameters, satParameters2);

simParameters.Fairness.MaxUsersPerBs = maxUsersPerBs;
simParameters.Fairness.Joint = false;
rng(42);
[bestNodeVec_perbs2, bestNodeTypeVec_perbs2, ~, ~, ~, capacityMbps_perbs2] = ...
    simulateScenario(bs_geo, user_geo, sat_geo2, wgs84, simParameters, satParameters2);

simParameters.Fairness.Joint = true;
rng(42);
[bestNodeVec_joint2, bestNodeTypeVec_joint2, ~, ~, ~, capacityMbps_joint2] = ...
    simulateScenario(bs_geo, user_geo, sat_geo2, wgs84, simParameters, satParameters2);

numBs1_off2   = onBs(bestNodeVec_off2,   1); numBs2_off2   = onBs(bestNodeVec_off2,   2);
numBs1_perbs2 = onBs(bestNodeVec_perbs2, 1); numBs2_perbs2 = onBs(bestNodeVec_perbs2, 2);
numBs1_joint2 = onBs(bestNodeVec_joint2, 1); numBs2_joint2 = onBs(bestNodeVec_joint2, 2);
numOutage_perbs2 = sum(bestNodeTypeVec_perbs2 == "Outage");
numOutage_joint2 = sum(bestNodeTypeVec_joint2 == "Outage");

fprintf('\n=== Υπο-σενάριο 2: Terrestrial-only (χωρίς δορυφόρο, ίδια γεωμετρία) ===\n');
fprintf('%-30s %6s %6s %10s\n', '', 'BS1', 'BS2', 'Outage');
fprintf('%-30s %6d %6d %10d\n', 'Χωρίς εξισορρόπηση:', numBs1_off2,   numBs2_off2,   sum(bestNodeTypeVec_off2=="Outage"));
fprintf('%-30s %6d %6d %10d\n', 'PerBs:',              numBs1_perbs2, numBs2_perbs2, numOutage_perbs2);
fprintf('%-30s %6d %6d %10d\n', 'Joint:',               numBs1_joint2, numBs2_joint2, numOutage_joint2);
fprintf('Μέση χωρητικότητα (Mbps) - Χωρίς: %.2f | PerBs: %.2f | Joint: %.2f\n', ...
    mean(capacityMbps_off2), mean(capacityMbps_perbs2), mean(capacityMbps_joint2));
fprintf('PerBs μετακίνησε 0 χρηστών (καμία DualConnectivity επιλεξιμότητα) - BS1 παραμένει στο %d, %d πάνω από το όριο.\n', ...
    numBs1_perbs2, max(numBs1_perbs2 - maxUsersPerBs, 0));

T2 = table((1:numOverlapUsers)', bestNodeVec_off2, capacityMbps_off2, ...
    bestNodeVec_perbs2, capacityMbps_perbs2, bestNodeVec_joint2, capacityMbps_joint2, ...
    'VariableNames', {'User','Node_off','Capacity_off_Mbps', ...
    'Node_perbs','Capacity_perbs_Mbps','Node_joint','Capacity_joint_Mbps'});
writetable(T2, fullfile(outputDir, 'joint_fairness_terrestrial_only.csv'));

results.dcmix.numBs1 = [numBs1_off numBs1_perbs numBs1_joint];
results.dcmix.numBs2 = [numBs2_off numBs2_perbs numBs2_joint];
results.dcmix.numSat = [numSat_off numSat_perbs numSat_joint];
results.dcmix.meanCapacity = [mean(capacityMbps_off) mean(capacityMbps_perbs) mean(capacityMbps_joint)];
results.terrOnly.numBs1 = [numBs1_off2 numBs1_perbs2 numBs1_joint2];
results.terrOnly.numBs2 = [numBs2_off2 numBs2_perbs2 numBs2_joint2];
results.terrOnly.meanCapacity = [mean(capacityMbps_off2) mean(capacityMbps_perbs2) mean(capacityMbps_joint2)];

%% ------------------ Γραφήματα ------------------
fig = figure('Visible','off', 'Position', [100 100 1100 450]);

subplot(1,2,1);
bar([capacityMbps_off, capacityMbps_perbs, capacityMbps_joint]);
xlabel('Χρήστης'); ylabel('Χωρητικότητα (Mbps)');
legend('Χωρίς εξισορρόπηση', 'PerBs', 'Joint', 'Location', 'best');
title('DC-mix: BS1+BS2 επικάλυψη, δορυφόρος διαθέσιμος');

subplot(1,2,2);
bar([capacityMbps_off2, capacityMbps_perbs2, capacityMbps_joint2]);
xlabel('Χρήστης'); ylabel('Χωρητικότητα (Mbps)');
legend('Χωρίς εξισορρόπηση', 'PerBs', 'Joint', 'Location', 'best');
title('Terrestrial-only: ίδια γεωμετρία, χωρίς δορυφόρο');

sgtitle(sprintf('Joint vs PerBs load balancing (όριο %d χρήστες/BS)', maxUsersPerBs));
saveas(fig, fullfile(outputDir, 'joint_fairness_demo.png'));
close(fig);

fprintf('\njointFairnessDemo results -> %s\n', outputDir);

end
