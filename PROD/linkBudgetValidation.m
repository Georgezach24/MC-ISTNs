function T = linkBudgetValidation(outputDir, label)
% Συστηματική επαλήθευση του ισοζυγίου ζεύξης και της κατανομής πόρων.
%
% Οι δύο υπάρχοντες έλεγχοι καλύπτουν τη γεωμετρία (geometryValidation.m) και
% την εξάρτηση της ενέργειας από το φορτίο (energyModelValidation.m). Εδώ
% ελέγχονται τα ενδιάμεσα μεγέθη, δηλαδή ό,τι μεσολαβεί ανάμεσα στη γεωμετρία
% και στο τελικό αποτέλεσμα:
%
%   Α. Μετατροπές μονάδων ισχύος (dBm, dBW, W) και τα δύο EIRP.
%   Β. Απώλειες ελεύθερου χώρου και η γνωστή κλίση των 6,02 dB ανά οκτάβα.
%   Γ. Ισχύς θορύβου kTB και η εξάρτησή της από το εύρος ζώνης.
%   Δ. Ταυτότητα SINR = EIRP - απώλειες - θόρυβος, πάνω στα μεγέθη που
%      επιστρέφει η ίδια η simulateScenario.
%   Ε. Γνωστή συμπεριφορά ισχύος: +-3 dB δίνουν +-3 dB χωρίς παρεμβολή και
%      αυστηρά λιγότερο με παρεμβολή.
%   ΣΤ. Διατήρηση εύρους ζώνης: το άθροισμα των μεριδίων ισούται με το
%      διαθέσιμο εύρος του κόμβου.
%   Ζ. Ισοζύγιο ισχύος δικτύου χωρίς διπλή καταμέτρηση.
%   Η. Οριακές περιπτώσεις με προκαθορισμένη συμπεριφορά: κανένας διαθέσιμος
%      κόμβος, ένας χρήστης, ισοπαλία υποψηφίων, ανύψωση ακριβώς στη μάσκα.
%
% Κάθε έλεγχος συγκρίνει μια αναμενόμενη τιμή, υπολογισμένη ανεξάρτητα από
% κλειστό τύπο, με την τιμή που παράγει ο κώδικας της προσομοίωσης. Τυπώνει
% γραμμή ΟΚ/ΑΠΕΤΥΧΕ ανά έλεγχο, γράφει CSV/PNG και τερματίζει με σφάλμα αν
% έστω ένας αποτύχει.

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
end
if nargin < 2
    label = '';
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

wgs84 = wgs84Ellipsoid;

%% ------------------ Παράμετροι (ίδιες με το runSimulation.m) ------------------
baseLat = 37.9838;
baseLon = 23.7275;

simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.AntennaGain = 8;
simParameters.NumAntennaElements = 10;
simParameters.RxNoiseFigure = 9;
simParameters.RxAntTemperature = 290;
simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);
simParameters.Power.NumTrx = 4;
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

satParameters.CarrierFrequency = 2.0e9;
satParameters.Bandwidth = 30e6;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + ...
                     10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 20;
satParameters.Power.Pfix  = 0;
satParameters.Power.EtaPA = 0.4;

% Μεγέθη που ξαναϋπολογίζονται εδώ, ανεξάρτητα από τη simulateScenario.
kBoltz = physconst('Boltzmann');
cLight = physconst('LightSpeed');
bwBsHz = simParameters.Carrier.NSizeGrid * 12 * ...
         simParameters.Carrier.SubcarrierSpacing * 1e3;
nfLin  = 10^(simParameters.RxNoiseFigure/10);
teqK   = simParameters.RxAntTemperature + 290*(nfLin - 1);
noiseBsDbw  = 10*log10(kBoltz * teqK * bwBsHz);
noiseSatDbw = 10*log10(kBoltz * teqK * satParameters.Bandwidth);
maxSe       = 5.5547;   % MCS 28, 64QAM, TS 38.214 Πίν. 5.1.3.1-1

% Συσσωρευτής ελέγχων
C = struct('name', {}, 'unit', {}, 'expected', {}, 'actual', {}, ...
           'tol', {}, 'pass', {});

%% ================== Α. Μονάδες ισχύος ==================
pOutTotalW    = 10^((simParameters.TxPower - 30)/10);
pOutPerChainW = pOutTotalW / simParameters.Power.NumTrx;

C = addCheck(C, 'A1 ισχυς εκπομπης σταθμου', 'W', 79.432823472428, pOutTotalW, 1e-9);
C = addCheck(C, 'A2 ισχυς ανα αλυσιδα πομποδεκτη', 'W', 19.858205868107, pOutPerChainW, 1e-9);
C = addCheck(C, 'A3 οριο EARTH ανα αλυσιδα (<= 20 W)', 'λογικο', 1, double(pOutPerChainW <= 20), 0);
C = addCheck(C, 'A4 συνθετο EIRP σταθμου', 'dBm', 67, simParameters.EIRP, 1e-12);
C = addCheck(C, 'A5 ιδιο EIRP σε γραμμικη ισχυ', 'W', 5011.872336272722, ...
             10^((simParameters.EIRP - 30)/10), 1e-6);
% Επιστροφή στην πυκνότητα EIRP: αν το 34 dBW/MHz είχε διαβαστεί ως ισχύς
% ενισχυτή, αυτός ο έλεγχος δεν θα έκλεινε (παρατήρηση #4 της αξιολόγησης).
densityBack = satParameters.EIRP - 30 - 10*log10(satParameters.Bandwidth/1e6);
C = addCheck(C, 'A6 πυκνοτητα EIRP δορυφορου (επιστροφη)', 'dBW/MHz', 34, densityBack, 1e-12);
C = addCheck(C, 'A7 EIRP δορυφορου', 'dBm', 78.771212547197, satParameters.EIRP, 1e-9);
C = addCheck(C, 'A8 ισχυς στην εισοδο κεραιας δορυφορου', 'dBm', 48.771212547197, ...
             satParameters.TxPower, 1e-9);

%% ================== Β. Απώλειες ελεύθερου χώρου ==================
lambdaSat  = cLight / satParameters.CarrierFrequency;
dFspl      = [600e3 1200e3 1754e3];
fsplClosed = 20*log10(4*pi*dFspl/lambdaSat);
fsplCode   = arrayfun(@(d) fspl(d, lambdaSat), dFspl);
C = addCheck(C, 'B1 fspl == 20log10(4*pi*d/lambda)', 'dB', 0, ...
             max(abs(fsplCode - fsplClosed)), 1e-9);
C = addCheck(C, 'B2 fspl στα 600 km, 2.0 GHz', 'dB', 154.031408142836, fsplCode(1), 1e-9);
octaveD    = 10.^(2:0.5:6);
fsplOctave = arrayfun(@(d) fspl(2*d, lambdaSat) - fspl(d, lambdaSat), octaveD);
C = addCheck(C, 'B3 διπλασιασμος αποστασης -> +6.02 dB', 'dB', 20*log10(2), ...
             max(fsplOctave), 1e-9);

% Το επίγειο μοντέλο δεν είναι ελεύθερος χώρος: ο εκθέτης απωλειών ξεπερνά το
% 2, οπότε η ίδια οκτάβα κοστίζει περισσότερο από 6,02 dB. Μετριέται χωρίς
% σκίαση και χωρίς διαλείψεις, απευθείας από τη nrPathLoss.
plCfgU = simParameters.PathLoss;
plU1 = nrPathLoss(plCfgU, simParameters.CarrierFrequency, false, [0;0;25], [500;0;1.5]);
plU2 = nrPathLoss(plCfgU, simParameters.CarrierFrequency, false, [0;0;25], [1000;0;1.5]);
umaOctaveDb = plU2 - plU1;
C = addCheck(C, 'B4 επιγεια οκταβα > 6.02 dB (εκθετης > 2)', 'λογικο', 1, ...
             double(umaOctaveDb > 20*log10(2)), 0);

%% ================== Γ. Ισχύς θορύβου ==================
C = addCheck(C, 'G1 ευρος ζωνης σταθμου', 'Hz', 18.36e6, bwBsHz, 1e-6);
C = addCheck(C, 'G2 ισοδυναμη θερμοκρασια θορυβου', 'K', ...
             290 + 290*(10^0.9 - 1), teqK, 1e-9);
% Δεύτερη, ανεξάρτητη διατύπωση του ίδιου μεγέθους: γινόμενο μέσα στον
% λογάριθμο έναντι αθροίσματος λογαρίθμων.
noiseBsSplit = 10*log10(kBoltz) + 10*log10(teqK) + 10*log10(bwBsHz);
C = addCheck(C, 'G3 kTB: γινομενο == αθροισμα σε dB', 'dB', 0, ...
             abs(noiseBsDbw - noiseBsSplit), 1e-9);
C = addCheck(C, 'G4 ισχυς θορυβου σταθμου', 'dBW', -122.336460425576, noiseBsDbw, 1e-9);
C = addCheck(C, 'G5 ισχυς θορυβου δορυφορου', 'dBW', -120.203974647031, noiseSatDbw, 1e-9);
noiseDouble = 10*log10(kBoltz * teqK * 2*bwBsHz);
C = addCheck(C, 'G6 διπλασιο ευρος ζωνης -> +3.01 dB θορυβου', 'dB', 10*log10(2), ...
             noiseDouble - noiseBsDbw, 1e-9);

%% ================== Δ. Ταυτότητα SINR ==================
% Σενάριο αναφοράς: δύο σταθμοί, πέντε χρήστες εντός πεδίου ισχύος, δορυφόρος
% στο ζενίθ. Η ταυτότητα EIRP - PL - N δίνει πλέον το SINR μόνο όταν δεν
% υπάρχει παρεμβολέας, οπότε ελέγχεται σε δύο βήματα: με έναν σταθμό πρέπει να
% ισχύει ακριβώς, με δύο πρέπει να παραβιάζεται προς τη σωστή κατεύθυνση.
bs2 = [baseLat baseLon 25; baseLat + 1200/111320 baseLon 25];
bs1 = bs2(1,:);
nU  = 5;
uD  = [300 700 1200 2500 4200];
uGeo = zeros(nU,3);
for u = 1:nU
    uGeo(u,:) = [baseLat + uD(u)*cosd(37*u)/111320, ...
                 baseLon + uD(u)*sind(37*u)/(111320*cosd(baseLat)), 1.5];
end
satOverhead = [baseLat baseLon 600e3];

rng(7);
[nodeVec, ~, ~, ~, sinrBest, capMbps, ~, ~, ~, ...
 bsSinr, ~, bsPl, ~, ~, satPl, satSinr, ...
 ~, netE] = ...
    simulateScenario(bs2, uGeo, satOverhead, wgs84, simParameters, satParameters);

% Ένας σταθμός: κανένας παρεμβολέας, άρα I = 0 και SINR = C/N ακριβώς.
rng(21);
[~,~,~,~,~,~,~,~,~, bsSinr1, ~, bsPl1] = ...
    simulateScenario(bs1, uGeo, satOverhead, wgs84, simParameters, satParameters);
sinrFromPl1 = (simParameters.EIRP - 30) - bsPl1 - noiseBsDbw;
C = addCheck(C, 'D1 ενας σταθμος: SINR == EIRP - PL - N', 'dB', 0, ...
             max(abs(bsSinr1 - sinrFromPl1)), 1e-9);

% Δορυφόρος: ο όρος παρεμβολής είναι μηδενικός (2.0 έναντι 3.5 GHz, ένας
% δορυφόρος), οπότε η ταυτότητα ισχύει ακριβώς και εκεί.
satSinrFromPl = (satParameters.EIRP - 30) - satPl - noiseSatDbw;
C = addCheck(C, 'D2 SINR δορυφορικο == EIRP - PL - N', 'dB', 0, ...
             max(abs(satSinr - satSinrFromPl)), 1e-9);
C = addCheck(C, 'D3 επιλεγμενο SINR == max(επιγειο, δορυφορικο)', 'dB', 0, ...
             max(abs(sinrBest - max(bsSinr, satSinr))), 1e-12);

% Δύο σταθμοί: το SINR πρέπει να είναι ΑΥΣΤΗΡΑ κάτω από την τιμή χωρίς
% παρεμβολή, για κάθε χρήστη. Αν η παρεμβολή δεν είχε μπει πραγματικά στον
% υπολογισμό, αυτός ο έλεγχος θα έδειχνε μηδενική διαφορά.
sinrNoInterf = (simParameters.EIRP - 30) - bsPl - noiseBsDbw;
interfDropDb = sinrNoInterf - bsSinr;
C = addCheck(C, 'D4 δυο σταθμοι: SINR < SNR για καθε χρηστη', 'λογικο', 1, ...
             double(all(interfDropDb > 0)), 0);
C = addCheck(C, 'D5 η πτωση αντιστοιχει σε θετικη ισχυ παρεμβολης', 'λογικο', 1, ...
             double(all(isfinite(interfDropDb)) && max(interfDropDb) < 60), 0);

%% ================== Ε. Μεταβολή ισχύος εκπομπής ==================
% Με έναν σταθμό η ζεύξη είναι περιορισμένη από τον θόρυβο: +-3 dB στην ισχύ
% μετατοπίζουν το SINR κατά ακριβώς +-3 dB. Με δύο σταθμούς η ίδια μεταβολή
% εφαρμόζεται ΚΑΙ στον παρεμβολέα, οπότε η ωφέλεια είναι μικρότερη: αυτό είναι
% το χαρακτηριστικό γνώρισμα ενός συστήματος περιορισμένου από παρεμβολή και
% ελέγχεται ρητά.
offsets    = [-3 0 3];
sinrShift1 = nan(numel(offsets),1);   % ένας σταθμός
sinrShift2 = nan(numel(offsets),1);   % δύο σταθμοί
for k = 1:numel(offsets)
    sp = simParameters;
    sp.TxPower = simParameters.TxPower + offsets(k);
    sp.EIRP = sp.TxPower + sp.AntennaGain + 10*log10(sp.NumAntennaElements);
    % Το πλήθος αλυσίδων δεν εισέρχεται στο SINR, μόνο στο ενεργειακό μοντέλο:
    % διπλασιάζεται στο +3 dB ώστε να μη σπάσει το όριο P_max των 20 W.
    if offsets(k) > 0
        sp.Power.NumTrx = 8;
    end
    rng(21);   % ίδιες κληρώσεις: μόνο η ισχύς μεταβάλλεται
    [~,~,~,~,~,~,~,~,~, bsSinrK1] = ...
        simulateScenario(bs1, uGeo, satOverhead, wgs84, sp, satParameters);
    sinrShift1(k) = mean(bsSinrK1 - bsSinr1);

    rng(7);
    [~,~,~,~,~,~,~,~,~, bsSinrK2] = ...
        simulateScenario(bs2, uGeo, satOverhead, wgs84, sp, satParameters);
    sinrShift2(k) = mean(bsSinrK2 - bsSinr);
end
C = addCheck(C, 'E1 ενας σταθμος: -3 dB ισχυος -> -3 dB SINR', 'dB', -3, sinrShift1(1), 1e-9);
C = addCheck(C, 'E2 ιδια ισχυς -> μηδενικη μετατοπιση', 'dB', 0, sinrShift1(2), 1e-12);
C = addCheck(C, 'E3 ενας σταθμος: +3 dB ισχυος -> +3 dB SINR', 'dB', 3, sinrShift1(3), 1e-9);
C = addCheck(C, 'E4 δυο σταθμοι: +3 dB ισχυος δινει ΛΙΓΟΤΕΡΟ απο +3 dB', 'λογικο', 1, ...
             double(sinrShift2(3) > 0 && sinrShift2(3) < 3 - 1e-6), 0);
C = addCheck(C, 'E5 δυο σταθμοι: -3 dB ισχυος κοστιζει ΛΙΓΟΤΕΡΟ απο 3 dB', 'λογικο', 1, ...
             double(sinrShift2(1) < 0 && sinrShift2(1) > -3 + 1e-6), 0);

%% ================== ΣΤ. Διατήρηση εύρους ζώνης ==================
% Αντίστροφος υπολογισμός του μεριδίου κάθε χρήστη από τη χωρητικότητα:
% C_u = B_u * SE_u  =>  B_u = C_u / SE_u. Το άθροισμα ανά κόμβο πρέπει να
% δίνει ακριβώς το εύρος ζώνης του κόμβου, ούτε λιγότερο ούτε περισσότερο.
seU = min(log2(1 + 10.^(sinrBest/10)), maxSe);
bwU = (capMbps*1e6) ./ seU;

activeNodes   = unique(nodeVec(nodeVec ~= "None"));
bwSumPerNode  = zeros(numel(activeNodes),1);
nodeBwPerNode = zeros(numel(activeNodes),1);
nodeLabels    = cell(numel(activeNodes),1);
for k = 1:numel(activeNodes)
    m = (nodeVec == activeNodes(k));
    bwSumPerNode(k) = sum(bwU(m));
    if startsWith(activeNodes(k), "SAT")
        nodeBwPerNode(k) = satParameters.Bandwidth;
    else
        nodeBwPerNode(k) = bwBsHz;
    end
    nodeLabels{k} = char(activeNodes(k));
end
C = addCheck(C, 'ST1 αθροισμα μεριδιων == ευρος ζωνης κομβου', 'Hz', 0, ...
             max(abs(bwSumPerNode - nodeBwPerNode)), 1e-3);
C = addCheck(C, 'ST2 κανενας κομβος δεν υπερδιαθετει', 'λογικο', 1, ...
             double(all(bwSumPerNode <= nodeBwPerNode + 1e-3)), 0);

%% ================== Ζ. Ισοζύγιο ισχύος δικτύου ==================
bsFullW  = simParameters.Power.NumTrx * ...
           (simParameters.Power.P0 + simParameters.Power.DeltaP*pOutPerChainW);
bsRfW    = simParameters.Power.NumTrx * simParameters.Power.DeltaP * pOutPerChainW;
bsSleepW = simParameters.Power.NumTrx * simParameters.Power.Psleep;
pOutSatW = 10^((satParameters.TxPower - 30)/10);
satFullW = satParameters.Power.Pfix + pOutSatW/satParameters.Power.EtaPA;

C = addCheck(C, 'Z1 ισχυς ενεργου σταθμου (EARTH)', 'W', 893.334270320412, bsFullW, 1e-9);
C = addCheck(C, 'Z2 τμημα σταθμου που εξαρταται απο τον ενισχυτη', 'W', 373.334270320412, bsRfW, 1e-9);
C = addCheck(C, 'Z3 ισχυς δορυφορου', 'W', 188.391482363219, satFullW, 1e-9);
expTotal = netE.ActiveBs*bsFullW + netE.IdleBs*bsSleepW + netE.SatActive*satFullW;
C = addCheck(C, 'Z4 συνολικη ισχυς == ενεργοι + αδρανεις + δορυφορος', 'W', 0, ...
             abs(netE.TotalPower_W - expTotal), 1e-9);
satRfW = pOutSatW / satParameters.Power.EtaPA;   % το τμήμα του ενισχυτή, χωρίς Pfix
expRf = netE.ActiveBs*bsRfW + netE.SatActive*satRfW;
C = addCheck(C, 'Z5 συμμετρικη ισχυς == μονο ενισχυτες', 'W', 0, ...
             abs(netE.TotalPowerRf_W - expRf), 1e-9);
C = addCheck(C, 'Z6 συμμετρικη <= πληρης εμβελεια', 'λογικο', 1, ...
             double(netE.TotalPowerRf_W <= netE.TotalPower_W), 0);
C = addCheck(C, 'Z7 καθε σταθμος μετριεται μια φορα', 'πληθος', size(bs2,1), ...
             netE.ActiveBs + netE.IdleBs, 0);

%% ================== Η. Οριακές περιπτώσεις ==================
% --- Η1/Η2: κανένας διαθέσιμος κόμβος, και ισοπαλία -Inf έναντι -Inf ---
farUser = [baseLat + 20000/111320, baseLon, 1.5];
satLow  = [baseLat + 2500e3/111320, baseLon, 600e3];   % ανύψωση ~1.5 deg, κάτω από τη μάσκα
rng(11);
[nodeN, typeN, ~, ~, sinrN, capN, ~, powN, eN, ...
 bsSinrN, ~, ~, ~, ~, ~, satSinrN, ~, netN, ~, thrN, bsRN, satRN] = ...
    simulateScenario([baseLat baseLon 25], farUser, satLow, wgs84, ...
                     simParameters, satParameters);

C = addCheck(C, 'H1 κανενας κομβος -> ServingNode "None"', 'λογικο', 1, ...
             double(nodeN == "None"), 0);
C = addCheck(C, 'H1 κανενας κομβος -> ServingType "Outage"', 'λογικο', 1, ...
             double(typeN == "Outage"), 0);
C = addCheck(C, 'H1 χωρητικοτητα απροσδιοριστη', 'λογικο', 1, double(isnan(capN)), 0);
C = addCheck(C, 'H1 ρυθμαποδοση μηδενικη', 'Mbps', 0, thrN, 0);
C = addCheck(C, 'H1 ισχυς κομβου απροσδιοριστη', 'λογικο', 1, double(isnan(powN)), 0);
C = addCheck(C, 'H1 ενεργεια ανα bit απειρη', 'λογικο', 1, double(isinf(eN) && eN > 0), 0);
C = addCheck(C, 'H1 αιτια επιγειας μη διαθεσιμοτητας', 'λογικο', 1, ...
             double(bsRN == "OutOfModelRange"), 0);
C = addCheck(C, 'H1 αιτια δορυφορικης μη διαθεσιμοτητας', 'λογικο', 1, ...
             double(satRN == "NotVisible"), 0);
C = addCheck(C, 'H1 η καταναλωση δικτυου δεν μηδενιζεται', 'W', bsSleepW, ...
             netN.TotalPower_W, 1e-9);
% Η σύγκριση των δύο υποψηφίων γίνεται με -Inf και στις δύο πλευρές. Ο
% τελεστής είναι αυστηρός (>), οπότε η ισοπαλία δεν αναθέτει κόμβο.
C = addCheck(C, 'H2 ισοπαλια -Inf: κανενας υποψηφιος δεν κερδιζει', 'λογικο', 1, ...
             double(isnan(bsSinrN) && satSinrN == -Inf && sinrN == -Inf && nodeN == "None"), 0);

% --- Η3: ένας χρήστης παίρνει ολόκληρο το εύρος ζώνης ---
oneUser = [baseLat + 400/111320, baseLon, 1.5];
satNo   = satParameters; satNo.MinElevationDeg = 95;   % δορυφόρος απρόσιτος
rng(13);
[~, type1, ~, ~, sinr1, cap1] = simulateScenario([baseLat baseLon 25], oneUser, ...
    satOverhead, wgs84, simParameters, satNo);
se1 = min(log2(1 + 10^(sinr1/10)), maxSe);
C = addCheck(C, 'H3 ενας χρηστης -> επιγεια εξυπηρετηση', 'λογικο', 1, ...
             double(type1 == "Terrestrial"), 0);
C = addCheck(C, 'H3 ενας χρηστης -> ολοκληρο το ευρος ζωνης', 'Hz', bwBsHz, ...
             cap1*1e6/se1, 1e-3);

% --- Η4: ανύψωση ακριβώς στη μάσκα ---
% Η ανύψωση υπολογίζεται με τον ίδιο τρόπο που τη βρίσκει η simulateScenario,
% ώστε η μάσκα να τεθεί ακριβώς πάνω της και να φανεί αν το όριο είναι κλειστό.
satEdge = [baseLat + 900e3/111320, baseLon, 600e3];
[~, elevEdge] = geodetic2aer(satEdge(1), satEdge(2), satEdge(3), ...
                             oneUser(1), oneUser(2), oneUser(3), wgs84);
satAt   = satParameters; satAt.MinElevationDeg   = elevEdge;
satJust = satParameters; satJust.MinElevationDeg = elevEdge + 1e-9;
rng(17);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, sinrAt, ~,~,~,~,~, reasonAt] = ...
    simulateScenario([baseLat baseLon 25], oneUser, satEdge, wgs84, simParameters, satAt);
rng(17);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, sinrJust, ~,~,~,~,~, reasonJust] = ...
    simulateScenario([baseLat baseLon 25], oneUser, satEdge, wgs84, simParameters, satJust);
C = addCheck(C, 'H4 ανυψωση ακριβως στη μασκα -> ορατος', 'λογικο', 1, ...
             double(isfinite(sinrAt) && reasonAt ~= "NotVisible"), 0);
C = addCheck(C, 'H4 ενα nanoβαθμο πιο πανω -> μη ορατος', 'λογικο', 1, ...
             double(sinrJust == -Inf && reasonJust == "NotVisible"), 0);

% --- Η5: πεπερασμένες ισοπαλίες είναι γεγονός μηδενικής πιθανότητας ---
% Ο κανόνας ορίζει ότι σε ισοπαλία κερδίζει ο επίγειος (αυστηρό >). Εδώ
% μετριέται πόσες φορές εμφανίζεται πεπερασμένη ισοπαλία σε μεγάλο δείγμα
% συγκρίσεων: αν εμφανιζόταν, ο κανόνας θα χρειαζόταν και αριθμητικό έλεγχο.
tieCount = 0; nCmp = 0;
for r = 1:40
    rng(100 + r);
    [~,~,~,~,~,~,~,~,~, bsS, ~,~,~,~,~, satS] = ...
        simulateScenario(bs2, uGeo, satOverhead, wgs84, simParameters, satParameters);
    fin = isfinite(bsS) & isfinite(satS);
    tieCount = tieCount + sum(bsS(fin) == satS(fin));
    nCmp = nCmp + sum(fin);
end
C = addCheck(C, 'H5 πεπερασμενες ισοπαλιες σε δειγμα συγκρισεων', 'πληθος', 0, tieCount, 0);

%% ------------------ Πίνακας και εκτύπωση ------------------
T = table(string({C.name}'), string({C.unit}'), [C.expected]', [C.actual]', ...
    abs([C.actual]' - [C.expected]'), [C.tol]', logical([C.pass]'), ...
    'VariableNames', {'Check','Unit','Expected','Actual','AbsError','Tolerance','Pass'});

fprintf('\n=== Συστηματική επαλήθευση ισοζυγίου ζεύξης ===\n');
fprintf('Σενάριο αναφοράς: UMa, fc = %.2f GHz, %d σταθμοί, %d χρήστες, δορυφόρος στα %g km\n\n', ...
    simParameters.CarrierFrequency/1e9, size(bs2,1), nU, satOverhead(3)/1e3);
fprintf('%-50s %14s %14s   %s\n', 'Έλεγχος', 'Αναμενόμενο', 'Μετρημένο', 'Αποτέλεσμα');
for k = 1:height(T)
    fprintf('%-50s %14.6g %14.6g   %s\n', T.Check(k), T.Expected(k), T.Actual(k), ...
        boolText(T.Pass(k)));
end

nPass = sum(T.Pass); nAll = height(T);
fprintf('\nΣύνολο: %d/%d έλεγχοι πέρασαν.\n', nPass, nAll);
fprintf('Συγκρίσεις υποψηφίων για ισοπαλία: %d, πεπερασμένες ισοπαλίες: %d\n', nCmp, tieCount);
fprintf('Επίγεια οκτάβα 500->1000 m: %.2f dB (ελεύθερος χώρος: %.2f dB)\n\n', ...
    umaOctaveDb, 20*log10(2));

%% ------------------ Γράφημα ------------------
fig = figure('Visible','off','Position',[100 100 1000 340]);

subplot(1,3,1);
dPlot = logspace(2,6,60);
semilogx(dPlot, arrayfun(@(d) fspl(d, lambdaSat), dPlot), 'LineWidth', 1.4); hold on;
semilogx(dPlot, fspl(dPlot(1), lambdaSat) + 20*log10(dPlot/dPlot(1)), '--', 'LineWidth', 1.2);
grid on; xlabel('Απόσταση [m]'); ylabel('Απώλειες [dB]');
legend({'fspl','20 dB ανά δεκάδα'}, 'Location','southeast');
title('Ελεύθερος χώρος: +6,02 dB ανά οκτάβα');

subplot(1,3,2);
plot(offsets, sinrShift1, '-o', 'LineWidth', 1.4); hold on;
plot(offsets, sinrShift2, '-s', 'LineWidth', 1.4);
plot(offsets, offsets, '--', 'LineWidth', 1.2);
grid on; xlabel('Μεταβολή ισχύος εκπομπής [dB]'); ylabel('Μεταβολή SINR [dB]');
legend({'ένας σταθμός','δύο σταθμοί','κλίση 1'}, 'Location','southeast');
title('Ισχύς εκπομπής προς SINR');

subplot(1,3,3);
bar([nodeBwPerNode(:)/1e6, bwSumPerNode(:)/1e6]);
grid on; ylabel('Εύρος ζώνης [MHz]');
set(gca, 'XTick', 1:numel(nodeLabels), 'XTickLabel', nodeLabels);
legend({'διαθέσιμο','άθροισμα μεριδίων'}, 'Location','northoutside', ...
    'Orientation','horizontal');
title('Διατήρηση εύρους ζώνης');

csvPath = fullfile(outputDir, 'link_budget_validation.csv');
pngPath = fullfile(outputDir, 'link_budget_validation.png');
writetable(T, csvPath);
saveas(fig, pngPath);
close(fig);

fprintf('Αποτελέσματα -> %s\n', csvPath);

%% ------------------ Versioning ------------------
runParams = struct('rngSeeds', [7 11 13 17], 'scenario', 'UMa', ...
    'CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'Bandwidth_Hz', bwBsHz, 'TxPower_dBm', simParameters.TxPower, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'NoiseTemperature_K', teqK, ...
    'NoisePowerBs_dBW', noiseBsDbw, 'NoisePowerSat_dBW', noiseSatDbw, ...
    'SatEIRP_dBm', satParameters.EIRP, 'SatBandwidth_Hz', satParameters.Bandwidth, ...
    'MinElevationDeg', satParameters.MinElevationDeg, ...
    'MaxSpectralEfficiency', maxSe, 'numChecks', nAll, 'numPassed', nPass, ...
    'tieComparisons', nCmp, 'finiteTies', tieCount);
runParams.terrestrialPower = simParameters.Power;
runParams.satellitePower   = satParameters.Power;
saveRunVersion('linkBudgetValidation', runParams, {csvPath, pngPath}, label);

if nPass < nAll
    error('linkBudgetValidation:Failed', ...
        'Απέτυχαν %d από %d έλεγχοι ισοζυγίου ζεύξης.', nAll - nPass, nAll);
end

end

%% ================== Τοπικές συναρτήσεις ==================
function C = addCheck(C, name, unit, expected, actual, tol)
k = numel(C) + 1;
C(k).name     = name;
C(k).unit     = unit;
C(k).expected = expected;
C(k).actual   = actual;
C(k).tol      = tol;
C(k).pass     = double(abs(actual - expected) <= tol);
end

function s = boolText(tf)
if tf
    s = 'ΟΚ';
else
    s = 'ΑΠΕΤΥΧΕ';
end
end
