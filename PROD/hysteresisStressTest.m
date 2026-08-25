function results = hysteresisStressTest(numSteps, outputDir)
%HYSTERESISSTRESSTEST Επικυρώνει τη μηχανή υστέρησης/TTT του
% simulateScenario.m (updateLinkActivation, βλ. σχόλια εκεί) πάνω σε ένα
% σκόπιμα κατασκευασμένο "boundary" σενάριο, όχι στο κύριο σενάριο
% αναφοράς (temporalPassSimulation.m): εκεί οι χρήστες είναι ακίνητοι, οπότε
% το χωρικά συσχετισμένο shadow fading (Gudmundson 1991) "παγώνει" στην
% τιμή του πρώτου βήματος (ρ=1, καμία νέα τυχαιότητα ανά βήμα) και δεν
% υπάρχει καθόλου ταλάντωση κοντά στο κατώφλι να καταστείλει η υστέρηση -
% ένα πριν/μετά πάνω σε εκείνο το σενάριο θα έδειχνε ταυτόσημους αριθμούς.
%
% Εδώ κατασκευάζεται αντ' αυτού ένας μοναδικός "boundary" χρήστης σε
% απόσταση ~2750m από έναν BS (UMa, NLOS) - απόσταση συντονισμένη
% εμπειρικά (μέσω nrPathLoss πάνω στο ΙΔΙΟ μοντέλο, όχι εικασία) ώστε το
% μέσο SNR να πέφτει ακριβώς πάνω στο minUsableSnrDb (≈-7.53dB) - και ο
% χρήστης μετακινείται ελαφρώς (~2m/βήμα, τυχαία κατεύθυνση, ίδιο
% μονοπάτι και στα δύο configs) ώστε το ρ<1 στο correlatedLosState να
% συνεχίζει να εγχέει γνήσια νέα τυχαιότητα shadow fading ανά βήμα -
% ίδιο μοντέλο καναλιού με την κύρια προσομοίωση (TR 38.901 §7.4.1,
% Gudmundson 1991), όχι επινοημένος θόρυβος μέτρησης.
%
% Ο δορυφόρος τίθεται σκόπιμα ΠΟΤΕ ορατός (MinElevationDeg=90), ώστε το
% τεστ να απομονώνει την απόφαση ενεργοποίησης ΜΟΝΟ στο επίγειο σκέλος.
% Πρόκειται ρητά για σενάριο ΕΠΙΚΥΡΩΣΗΣ (stress test) του μηχανισμού
% απόφασης, όχι για μια νέα claim ρεαλισμού πάνω στο ίδιο το φυσικό
% μοντέλο καναλιού (που παραμένει αμετάβλητο, standards-grounded).
%
% Τρέχει ΔΥΟ configs πάνω στο ΙΔΙΟ πρόγραμμα rng ανά βήμα (rng(step) πριν
% από κάθε κλήση simulateScenario, ίδιο και στα δύο) - άρα ΙΔΙΕΣ
% πραγματοποιήσεις καναλιού, με μόνη διαφορά το simParameters.Hysteresis:
%   'off': MarginDb=0, TimeToTriggerSteps=0  (ο ακατέργαστος κανόνας κατωφλίου)
%   'on' : MarginDb=2, TimeToTriggerSteps=1  (ίδιες τιμές με το κύριο
%          σενάριο αναφοράς, temporalPassSimulation.m)
% ώστε η όποια διαφορά στα αποτελέσματα να οφείλεται αποκλειστικά στον
% μηχανισμό απόφασης, όχι σε διαφορετική τυχαία πραγματοποίηση καναλιού.
%
% Χρήση:
%   hysteresisStressTest();      % 80 βήματα -> ../Results
%   hysteresisStressTest(120);

if nargin < 1 || isempty(numSteps)
    numSteps = 80;
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

wgs84 = wgs84Ellipsoid;
bs_geo = [37.9838 23.7275 25];

% Boundary user ~2750m από τον BS -> μέσο NLOS SNR ≈ minUsableSnrDb.
boundaryStartLat = bs_geo(1);
boundaryStartLon = bs_geo(2) + 2750/(111320*cosd(bs_geo(1)));
user_geo0 = [boundaryStartLat, boundaryStartLon, 1.5];

sat_geo = [0 0 550e3];   % αδιάφορη θέση - ο δορυφόρος είναι ούτως ή άλλως ποτέ ορατός

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
satParameters.MinElevationDeg = 90;   % σκόπιμα ποτέ ορατός - βλ. σχόλιο κεφαλίδας
satParameters.Power.Pfix = 0;
satParameters.Power.EtaPA = 0.4;

% Μονοπάτι κίνησης (~18m/βήμα, τυχαία κατεύθυνση κάθε βήμα) - υπολογισμένο
% ΜΙΑ φορά, κοινό και στα δύο configs, ανεξάρτητο από το rng(step) του
% καναλιού. Η correlation distance UMa NLOS είναι 50m (TR 38.901 Πίνακας
% 7.5-6) - μια μετατόπιση της τάξης των μερικών m/βήμα (αρχική δοκιμή:
% 2m) αποδείχθηκε πολύ μικρή για ουσιαστική αποσυσχέτιση (ρ=exp(-2/50)
% ≈0.96, σχεδόν παγωμένο ακόμα και σε 80 βήματα) - 18m/βήμα δίνει ρ≈0.69,
% αρκετή ανάμειξη ώστε να εμφανιστεί γνήσια ταλάντωση γύρω από το
% κατώφλι μέσα σε λίγες δεκάδες βήματα.
stepSizeM = 18;
rng(12345);
jitterAngles = 2*pi*rand(numSteps,1);

configNames    = ["off", "on"];
configMarginDb = [0, 2];
configTttSteps = [0, 1];
allTables = struct();
toggles   = struct();

for c = 1:numel(configNames)
    simParameters.Hysteresis.MarginDb           = configMarginDb(c);
    simParameters.Hysteresis.TimeToTriggerSteps = configTttSteps(c);

    user_geo = user_geo0;
    channelState = [];
    rows = cell(numSteps,1);
    for step = 1:numSteps
        if step > 1
            dNorthM = stepSizeM*cos(jitterAngles(step));
            dEastM  = stepSizeM*sin(jitterAngles(step));
            user_geo(1) = user_geo(1) + dNorthM/110574;
            user_geo(2) = user_geo(2) + dEastM/(111320*cosd(user_geo(1)));
        end

        rng(step);   % ίδιο πρόγραμμα rng και στα δύο configs -> ίδιο κανάλι
        [bestNodeVec, bestNodeTypeVec, ~, ~, bestSnrDbVec, ~, ~, ~, ~, ...
            bestBsSnrDbVec, ~, ~, ~, ~, ~, ~, channelState] = ...
            simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, channelState);

        rows{step} = table(step, bestNodeVec, bestNodeTypeVec, bestSnrDbVec, bestBsSnrDbVec, ...
            'VariableNames', {'Step','ServingNode','ServingType','SNR_dB','CandBS_SNR_dB'});
    end
    T = vertcat(rows{:});
    writetable(T, fullfile(outputDir, "hysteresis_stress_" + configNames(c) + ".csv"));

    changed = T.ServingType(2:end) ~= T.ServingType(1:end-1);
    toggles.(configNames(c))   = sum(changed);
    allTables.(configNames(c)) = T;

    fprintf('[%s] MarginDb=%g TimeToTriggerSteps=%g -> %d μεταβάσεις κατάστασης σε %d βήματα\n', ...
        configNames(c), configMarginDb(c), configTttSteps(c), toggles.(configNames(c)), numSteps);
end

reduction = 100*(1 - toggles.on/max(toggles.off,1));
fprintf('\nΧωρίς hysteresis: %d μεταβάσεις | Με hysteresis: %d μεταβάσεις (μείωση %.0f%%)\n', ...
    toggles.off, toggles.on, reduction);

results.toggles_off = toggles.off;
results.toggles_on  = toggles.on;
results.reduction_pct = reduction;
results.table_off = allTables.off;
results.table_on  = allTables.on;

%% ------------------ Γράφημα σύγκρισης ------------------
minUsableSnrDb = 10*log10(2^0.2344 - 1);
fig = figure('Visible','off', 'Position', [100 100 800 700]);

subplot(3,1,1);
plot(allTables.off.Step, allTables.off.CandBS_SNR_dB, 'k-');
yline(minUsableSnrDb, 'b--');
yline(minUsableSnrDb + configMarginDb(2), 'r:');
yline(minUsableSnrDb - configMarginDb(2), 'r:');
ylabel('BS SNR (dB)');
title('SNR ζεύξης boundary χρήστη (ίδιο και στα δύο configs)');
legend('SNR', 'SNR_{min}', '\pm margin', 'Location', 'best');

subplot(3,1,2);
stairs(allTables.off.Step, double(allTables.off.ServingType == "Terrestrial"), 'b-', 'LineWidth', 1.3);
ylim([-0.2 1.2]); yticks([0 1]); yticklabels({'Outage','Terrestrial'});
title(sprintf('Χωρίς hysteresis (%d μεταβάσεις)', toggles.off));

subplot(3,1,3);
stairs(allTables.on.Step, double(allTables.on.ServingType == "Terrestrial"), 'r-', 'LineWidth', 1.3);
ylim([-0.2 1.2]); yticks([0 1]); yticklabels({'Outage','Terrestrial'});
xlabel('Βήμα');
title(sprintf('Με hysteresis, MarginDb=2, TTT=1 (%d μεταβάσεις)', toggles.on));

saveas(fig, fullfile(outputDir, 'hysteresis_stress_comparison.png'));
close(fig);

fprintf('\nhysteresisStressTest results -> %s\n', outputDir);

end
