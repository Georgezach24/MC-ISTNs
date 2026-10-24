function [convTable, summary] = runSimulation(opts)
%RUNSIMULATION Ενιαία προσομοίωση του σεναρίου ISTN: μία εκτέλεση που παράγει
% όλα τα δεδομένα  (αποτελέσματα προσομοίωσης + σύνολο εκπαίδευσης).
%
% Δομή του χρόνου, με την πηγή της κάθε επιλογής:
%
%   * Στοιχειώδης μονάδα = μία διέλευση του δορυφόρου. Το TR 38.821
%     Πίν. 4.2-3 NOTE 1 ορίζει ως φυσικό παράθυρο λειτουργίας τον χρόνο
%     ορατότητας του δορυφόρου ("a period of time corresponding to the
%     visibility time of the satellite"). Με numSats > 1 το παράθυρο
%     κεντράρεται στη διέλευση του μεσαίου δορυφόρου και περιέχει τις
%     διελεύσεις των διαδοχικών δορυφόρων του ίδιου τροχιακού επιπέδου.
%
%   * Βήμα Δt = 1 s. Δύο ανεξάρτητοι περιορισμοί δίνουν την ίδια τιμή:
%       - ITU-R M.2412-0, Παράρτημα 1 §5.3.2: η μετατόπιση του χρήστη ανά
%         βήμα πρέπει να μένει κάτω από 1 m. Στα 3 km/h -> 0.83 m ανά s.
%       - TR 38.821 §7.3.2.1.4, Πίν. 7.3.2.1.4-1: η ταχύτερη κλίμακα
%         κινητικότητας σε LEO NTN είναι 6.61 s (δέσμη 50 km).
%
%   * Συνολική διάρκεια: δεν ορίζεται από κανένα πρότυπο ως σταθερός αριθμός.
%     Το ITU-R M.2412-0 §7.1 ορίζει κριτήριο: "A sufficient number of drops is
%     simulated to ensure convergence in the UE and system performance metrics.
%     The proponent should provide information on the width of confidence
%     intervals". Το ITU-R M.2514-0 §8.2.4 επεκτείνει ρητά τη διαδικασία του
%     M.2412 στις δορυφορικές αξιολογήσεις. Η προσομοίωση τρέχει ΣΤΑΘΕΡΟ πλήθος
%     διελεύσεων (numPasses), ίδιο σε κάθε σενάριο, ώστε η διέλευση n να έχει
%     την ίδια τυχαιότητα παντού και τα σενάρια να συγκρίνονται κατά ζεύγη. Το
%     διάστημα εμπιστοσύνης δεν σταματά την εκτέλεση· αναφέρεται στο τέλος ως
%     η απαίτηση του §7.1, μαζί με την καμπύλη σύγκλισης που δείχνει ότι το
%     πλήθος αρκεί.
%
%   * Κίνηση χρηστών: ITU-R M.2412-0 §8.4, ΠΙΝΑΚΑΣ 5 b) και c): "Fixed and
%     identical speed |v| of all UEs of the same mobility class, randomly and
%     uniformly distributed direction", με 3 km/h για πεζό χρήστη. Κάθε
%     διέλευση ξεκινά με νέες θέσεις χρηστών (ένα "drop" κατά M.2412 §7.1)
%     και οι χρήστες περπατούν συνεχώς μέσα στη διέλευση.
%
% Χρήση:
%   runSimulation()                                   % προεπιλογές
%   runSimulation(struct('numPasses',10))             % σύντομο τρέξιμο δοκιμής
%   runSimulation(struct('geometryGrid',false))       % τυχαία γεωμετρία ανά διέλευση
%
% Η προσομοίωση προχωρά πάντα με βήμα dtSeconds και όλοι οι δείκτες
% υπολογίζονται σε αυτή την ανάλυση. Στο αρχείο όμως γράφεται μία γραμμή κάθε
% datasetStride βήματα: διαδοχικά δείγματα του ίδιου χρήστη απέχουν 0.83 m και
% είναι σχεδόν ταυτόσημα, οπότε η πλήρης ανάλυση θα παρήγαγε αρχείο εκατοντάδων
% MB χωρίς να προσθέτει πληροφορία. Η διέλευση αναφοράς αποθηκεύεται χωριστά σε
% πλήρη ανάλυση, όπου χρειάζεται για τις χρονοσειρές.
%
% Έξοδοι:
%   Dataset/dataset.csv        μία γραμμή ανά (διέλευση, δείγμα, χρήστη)· το
%                              ΜΟΝΟ csv της εκτέλεσης, είσοδος του ML μοντέλου
%   Results/policy_gain.png    πόσο καλύτερη είναι η υπό εξέταση πολιτική από
%                              καθεμία από τις δύο αποκλειστικές, με 95% CI
%   Results/*.png              χρονοσειρές διέλευσης αναφοράς + καμπύλη σύγκλισης
%   Results/simulation_results.mat  οι αριθμοί πίσω από τα γραφήματα (διέλευση
%                              αναφοράς, σύγκλιση, δείκτες πολιτικών, διαφορές)

%% ------------------ Επιλογές ------------------
if nargin < 1 || isempty(opts)
    opts = struct();
end
thisDir = fileparts(mfilename('fullpath'));
def = struct( ...
    'dtSeconds',    1, ...      % s, βήμα (M.2412 Παρ.1 §5.3.2 + TR 38.821 Πίν. 7.3.2.1.4-1)
    'ueSpeedKmh',   3, ...      % km/h, πεζός χρήστης (M.2412 ΠΙΝ. 5 b/c)
    'passWindowS',  900, ...    % s, παράθυρο ανά διέλευση γύρω από τη μέγιστη προσέγγιση
    'numPasses',    300, ...    % σταθερό πλήθος διελεύσεων (η εκτέλεση αναφοράς έπιασε το στόχο CI στις 227)
    'geometryGrid', true, ...   % στρωματοποιημένη εγκάρσια μετατόπιση ίχνους· false = rand() ανά διέλευση
    'ciTolerance',  0.05, ...   % στόχος αναφοράς: σχετικό ημιεύρος 95% CI (ρυθμός/ενέργεια)
    'ciToleranceFrac', 0.02, ...% στόχος αναφοράς: απόλυτο ημιεύρος 95% CI (ποσοστά)
    'datasetStride', 5, ...     % κάθε πόσα βήματα γράφεται γραμμή στο dataset
    'hysteresisDb',  3, ...      % TS 38.331 §6.3.2: 0..15 dB, βήμα 0.5 dB (τιμή πεδίου 6)
    'timeToTriggerMs', 2560, ... % TS 38.331 §6.3.2: απαριθμημένη τιμή, 3 βήματα στο Δt = 1 s
    'pingPongMtsS',  1, ...      % s, ελάχιστος χρόνος παραμονής (MTS) για τον ορισμό του ping-pong
    'l3FilterK',     4, ...      % TS 38.331 §6.3.2 filterCoefficient (DEFAULT fc4)· 0 = χωρίς φίλτρο
    'measPeriodMs',  200, ...    % TS 38.133 Πίν. 9.2.5.2-1, FR1 χωρίς DRX
    'measAveraging', true, ...   % μέση τιμή μέτρησης μέσα στο βήμα· false = στιγμιαίο δείγμα
    'numUsersRange', [10 30], ...% [min max] χρήστες ανά διέλευση, ομοιόμορφα
    'bsRings',       1, ...      % δακτύλιοι εξαγωνικού πλέγματος: 1 -> 7 σταθμοί
    'bsIsdM',        500, ...    % m, απόσταση σταθμών (TR 38.901 Πίν. 7.8-1, UMa)
    'numSats',       3, ...      % διαδοχικοί δορυφόροι στο ίδιο τροχιακό επίπεδο
    'satFreqReuse',  1, ...      % TR 38.821 Πίν. 6.1.1.1-5: Option 1 = 1, Option 2 = 3
    'mlbMaxOffsetDb', 6, ...     % dB ανά μονάδα διαφοράς φόρτου (Nasri & Altman: HMmax - f(0) = 12 - 6)
    'resume',       true, ...   % συνέχιση από σημείο ελέγχου, αν υπάρχει
    'rngSeed',      42, ...
    'datasetPath',  fullfile(thisDir, '..', 'Dataset', 'dataset.csv'), ...
    'handoverEventsPath', fullfile(thisDir, '..', 'Dataset', 'handover_events.mat'), ...
    'outputDir',    fullfile(thisDir, '..', 'Results'), ...
    'label',        '');
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = def.(fn{k});
    end
end

if ~isfolder(opts.outputDir)
    mkdir(opts.outputDir);
end
datasetDir = fileparts(opts.datasetPath);
if ~isempty(datasetDir) && ~isfolder(datasetDir)
    mkdir(datasetDir);
end
resultsMat     = fullfile(opts.outputDir, 'simulation_results.mat');
checkpointPath = fullfile(opts.outputDir, 'runSimulation_checkpoint.mat');

% Σημείο ελέγχου: μια μεγάλη εκτέλεση μπορεί να διακοπεί (π.χ. από έλλειψη
% μνήμης στο μηχάνημα). Αν υπάρχει συμβατό σημείο ελέγχου, η εκτέλεση
% συνεχίζει από την επόμενη διέλευση αντί να ξαναρχίσει από την αρχή.
resumeState = [];
if opts.resume && isfile(checkpointPath) && isfile(opts.datasetPath)
    S = load(checkpointPath);
    % Σημεία ελέγχου παλαιότερης μορφής δεν κρατούν τη διέλευση αναφοράς (τότε
    % ζούσε σε csv), οπότε δεν μπορούν να συνεχιστούν.
    if isfield(S,'ckpt') && isfield(S.ckpt,'refPassTable') && ...
            checkpointMatches(S.ckpt.opts, opts)
        % Το σημείο ελέγχου γράφεται αμέσως μετά την εγγραφή της διέλευσης στο
        % CSV. Αν η διακοπή έπεσε ακριβώς ανάμεσα στα δύο, το CSV έχει μία
        % διέλευση παραπάνω από όση ξέρει το σημείο ελέγχου· τότε η συνέχιση θα
        % παρήγαγε διπλές γραμμές, οπότε σταματάμε αντί να το αγνοήσουμε.
        lastInCsv = max(readmatrix(opts.datasetPath, 'Range', 'A:A', ...
            'NumHeaderLines', 1));
        if lastInCsv ~= S.ckpt.passIdx
            error('runSimulation:CheckpointMismatch', ...
                ['Το CSV φτάνει ως τη διέλευση %d ενώ το σημείο ελέγχου ως τη %d. ' ...
                 'Η προηγούμενη εκτέλεση διακόπηκε σε ακατάλληλη στιγμή. Σβήσε τα ' ...
                 '%s και %s και ξεκίνα από την αρχή.'], ...
                lastInCsv, S.ckpt.passIdx, opts.datasetPath, checkpointPath);
        end
        resumeState = S.ckpt;
        fprintf('Συνέχιση από σημείο ελέγχου: %d διελεύσεις ήδη ολοκληρωμένες.\n', ...
            resumeState.passIdx);
    else
        fprintf(['Βρέθηκε σημείο ελέγχου με διαφορετικές παραμέτρους - ' ...
                 'αγνοείται και η εκτέλεση ξεκινά από την αρχή.\n']);
    end
end

if isempty(resumeState)
    if isfile(opts.datasetPath),  delete(opts.datasetPath);  end
    if isfile(opts.handoverEventsPath), delete(opts.handoverEventsPath); end
    if isfile(resultsMat),        delete(resultsMat);        end
    if isfile(checkpointPath),    delete(checkpointPath);    end
end

wgs84 = wgs84Ellipsoid;
cLight = physconst('LightSpeed');

%% ------------------ Τοπολογία σταθμών βάσης (σταθερή: το σενάριο) ------------------
% Εξαγωνικό πλέγμα σταθμών με απόσταση ISD, όπως η διάταξη βαθμονόμησης του
% TR 38.901 §7.8 (Πίν. 7.8-1: UMa, ISD 500 m, ύψος BS 25 m). Με bsRings = 1:
% ένας σταθμός στο κέντρο και 6 στον πρώτο δακτύλιο. Το κέντρο είναι το ίδιο
% σημείο αναφοράς με το προηγούμενο σενάριο.
bsCenterGeo = [37.9838 23.7275];
bsHeightM   = 25;
hexOffsets  = hexSiteOffsetsM(opts.bsRings, opts.bsIsdM);   % [ανατολικά βόρεια] m
bs_geo = zeros(size(hexOffsets,1), 3);
for b = 1:size(hexOffsets,1)
    [bsLat, bsLon] = reckon(bsCenterGeo(1), bsCenterGeo(2), norm(hexOffsets(b,:)), ...
        atan2d(hexOffsets(b,1), hexOffsets(b,2)), wgs84);
    bs_geo(b,:) = [bsLat bsLon bsHeightM];
end
numBs = size(bs_geo,1);
bsClusterRadiusM = opts.bsRings * opts.bsIsdM;   % απόσταση του πιο μακρινού σταθμού από το κέντρο

numUsersRange = opts.numUsersRange;   % πλήθος χρηστών ανά διέλευση
ueHeightM     = 1.5;      % m, TR 38.901 §7.4.1 Πίν. 7.4.1-1 (υπαίθριος χρήστης)

%% ------------------ Ραδιο-configuration ------------------
% Επίγειο: TR 38.901 §7.8 Πίν. 7.8-1. BS gain = element (8 dBi) + array
% 10·log10(N), N=10, ιδανική στόχευση -> composite 18 dBi.
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;
simParameters.Carrier.SubcarrierSpacing = 30;
simParameters.Carrier.CyclicPrefix = 'Normal';
simParameters.CarrierFrequency = 3.5e9;
simParameters.AntennaGain = 8;              % dBi, element gain (TR 38.901 §7.3 Πίν. 7.3-1)
simParameters.NumAntennaElements = 10;      % TR 38.901 Πίν. 7.8-1
simParameters.RxNoiseFigure = 9;            % dB, UE downlink NF (TR 38.901 Πίν. 7.8-1)
simParameters.RxAntTemperature = 290;

simParameters.PathLossModel = '5G-NR';
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;                 % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMa)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);   % dBm, composite

simParameters.Power.NumTrx = 4;   % αλυσίδες πομποδέκτη ανά τομέα· P_out/αλυσίδα <= 20 W (EARTH Πίν. 2)
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

% Υστέρηση και χρόνος επιβεβαίωσης της απόφασης (TS 38.331 §5.5.4.4 Event A3,
% §6.3.2 Hysteresis / TimeToTrigger). Με μηδενικές τιμές ο κανόνας είναι ο
% στιγμιαίος argmax χωρίς περιθώριο.
simParameters.Mobility.HysteresisDb    = opts.hysteresisDb;
simParameters.Mobility.TimeToTriggerMs = opts.timeToTriggerMs;
simParameters.Mobility.DtSeconds       = opts.dtSeconds;
% Ping-pong: επιστροφή στον αμέσως προηγούμενο κόμβο μετά από παραμονή στον
% ενδιάμεσο το πολύ PingPongMtsS (ορισμός MTS του 3GPP TR 36.839). Με Δt = 1 s
% ο χρόνος παραμονής μετριέται σε ακέραια βήματα, οπότε το «< 1 s» του ορισμού
% γίνεται «<= 1 βήμα»: η μικρότερη παραμονή που μπορεί να παρατηρηθεί.
simParameters.Mobility.PingPongMtsS    = opts.pingPongMtsS;
% Mobility Load Balancing (TS 38.300 §15.5.1), μόνο στις πολιτικές *MLB: η
% μετατόπιση του A3 ως προς τον στόχο είναι MlbMaxOffsetDb επί τη διαφορά
% εικονικού φόρτου, περιορισμένη στο [-1, 1] (βλ. simulateScenario).
simParameters.Mobility.MlbMaxOffsetDb  = opts.mlbMaxOffsetDb;

% Δορυφόρος: TR 38.821 Set-1, LEO-600, S-band (Πίν. 6.1.1.1-1 & 6.1.3.2-1).
% Η πυκνότητα EIRP (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
satAltitude = 600e3;
satParameters.CarrierFrequency = 2.0e9;
% Το εύρος 30 MHz ορίζεται για το σύνολο των δεσμών και "must be adapted based
% on frequency reuse factor" (TR 38.821 Πίν. 6.1.1.1-5): κάθε δέσμη παίρνει
% 30/FRF MHz. Με σταθερή πυκνότητα EIRP η ισχύς ανά δέσμη ακολουθεί το εύρος.
satParameters.FreqReuse = opts.satFreqReuse;
satParameters.Bandwidth = 30e6 / satParameters.FreqReuse;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + 10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 20;

satParameters.Power.Pfix  = 0;       % W, εκτός ενισχυτή· Pfix=0 -> αισιόδοξη υπόθεση
satParameters.Power.EtaPA = 0.4;

%% ------------------ Πεδίο ισχύος επίγειων μοντέλων ------------------
% TR 38.901 §7.4.1, Πίν. 7.4.1-1. Ελέγχεται σε ΚΑΘΕ βήμα για κάθε χρήστη.
d2dMinM = 10;
d2dMaxM = 5000;

% Μέγιστη μετατόπιση χρήστη μέσα σε μία διέλευση. Οι αρχικές θέσεις
% δειγματοληπτούνται ώστε η απόσταση από τον ΠΛΗΣΙΕΣΤΕΡΟ σταθμό να μένει μέσα
% στα όρια σε όλη την τροχιά: ο χρήστης βρίσκεται στο όριο της επίγειας
% κάλυψης, όχι μέσα σε κυψέλη, όπως και στο σενάριο με τους δύο σταθμούς.
ueSpeedMps   = opts.ueSpeedKmh / 3.6;
maxWalkM     = ueSpeedMps * opts.passWindowS;
userRMinM    = d2dMinM + maxWalkM;
userRMaxM    = d2dMaxM - maxWalkM;
if userRMinM >= userRMaxM
    error('runSimulation:WindowTooLong', ...
        ['Με ταχύτητα %.2f m/s και παράθυρο %.0f s ο χρήστης διανύει %.0f m, ' ...
         'που δεν χωράει στο πεδίο ισχύος [%d, %d] m. Μείωσε το passWindowS.'], ...
        ueSpeedMps, opts.passWindowS, maxWalkM, d2dMinM, d2dMaxM);
end

%% ------------------ Μέτρηση πριν την απόφαση ------------------
% Το TS 38.331 §5.5.3.2 ορίζει ότι η σύγκριση A3 γίνεται πάνω σε
% φιλτραρισμένη μέτρηση, όχι στη στιγμιαία τιμή. Δύο σκέλη, και τα δύο
% τεκμηριωμένα (αντίγραφα των εδαφίων: Sources/TS38331_L3_filtering.md):
%
% (α) Το ίδιο το M_n δεν είναι στιγμιαίο δείγμα. Το TS 38.133 Πίν. 9.2.5.2-1
%     δίνει για FR1 χωρίς DRX περίοδο μέτρησης τουλάχιστον 200 ms, μέσα στην
%     οποία χωρούν πολλές ανεξάρτητες πραγματώσεις των γρήγορων διαλείψεων.
%     Το ισοδύναμο πλήθος τους υπολογίζεται απευθείας από την αυτοσυσχέτιση
%     ισχύος του μοντέλου Clarke και όχι από εμπειρικό χρόνο συνοχής
%     (βλ. effectiveSamplesClarke παρακάτω).
% (β) Πάνω σε αυτό εφαρμόζεται το φίλτρο L3. Ο συντελεστής k υποθέτει ρυθμό
%     δειγματοληψίας X ίσο με μία περίοδο μέτρησης L1, οπότε η ίδια η ρήτρα
%     απαιτεί προσαρμογή στο βήμα της προσομοίωσης:
%        a = 2^(-k/4) στα X ms   ->   a' = 1 - (1-a)^(Δt/X)
%     Με k = 0 δεν εφαρμόζεται φίλτρο (NOTE 1 του ίδιου εδαφίου).
cLight = physconst('LightSpeed');
if opts.measAveraging
    simParameters.Meas.NAvg = effectiveSamplesClarke( ...
        ueSpeedMps*simParameters.CarrierFrequency/cLight, opts.dtSeconds);
    satParameters.Meas.NAvg = effectiveSamplesClarke( ...
        ueSpeedMps*satParameters.CarrierFrequency/cLight, opts.dtSeconds);
else
    simParameters.Meas.NAvg = 1;
    satParameters.Meas.NAvg = 1;
end
if opts.l3FilterK <= 0
    l3Alpha = 1;                      % NOTE 1: k = 0 -> χωρίς φίλτρο
else
    aAtX    = 2^(-opts.l3FilterK/4);
    l3Alpha = 1 - (1 - aAtX)^(opts.dtSeconds*1000/opts.measPeriodMs);
end
simParameters.Meas.L3Alpha = l3Alpha;
satParameters.Meas.L3Alpha = l3Alpha;

%% ------------------ Τροχιά LEO (κυκλική Κεπλεριανή) ------------------
% Στοιχεία εφημερίδας κατά TR 38.821 Πίν. 7.3.6.1-1, με εκκεντρότητα μηδέν.
muEarth    = 3.986004418e14;   % m^3/s^2
Re         = 6378137;          % m, ισημερινή ακτίνα WGS84
omegaEarth = 7.2921150e-5;     % rad/s
inclDeg    = 53;               % μοίρες
a          = Re + satAltitude;
orbitalPeriodS = 2*pi*sqrt(a^3/muEarth);
meanMotion     = 2*pi/orbitalPeriodS;
groundSpeedMps = meanMotion * Re;

centerLat = mean(bs_geo(:,1));
centerLon = mean(bs_geo(:,2));

% Όρισμα πλάτους στο σημείο μέγιστης προσέγγισης: sin(lat) = sin(i)*sin(u).
uPeakRad = asin(min(max(sind(centerLat)/sind(inclDeg), -1), 1));
% Ορθή αναφορά ανερχόμενου δεσμού ώστε το ίχνος να περνά από το κέντρο.
lonPeakInertialRad = atan2(cosd(inclDeg)*sin(uPeakRad), cos(uPeakRad));
tPeakS   = uPeakRad / meanMotion;
raanBase = deg2rad(centerLon) + omegaEarth*tPeakS - lonPeakInertialRad;

% Κάθε διέλευση μετατοπίζει το ίχνος εγκάρσια, ώστε να καλυφθεί το φάσμα από
% οριακή έως ζενιθιακή διέλευση. Το μέγιστο χρήσιμο εύρος βρίσκεται αριθμητικά
% ως η μετατόπιση όπου η μέγιστη ανύψωση πέφτει στη μάσκα ορατότητας.
raanSpanRad = findMaxRaanOffset(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);

% Διαδοχικοί δορυφόροι στο ίδιο τροχιακό επίπεδο, σε ίσες αποστάσεις στο
% όρισμα πλάτους. Ο μεσαίος περνά από τη μέγιστη προσέγγιση στο μέσο του
% παραθύρου· ο SAT-1 προηγείται, άρα περνά πρώτος. Η απόσταση επιλέγεται ώστε
% οι numSats δορυφόροι να είναι ΑΚΡΙΒΩΣ όσοι μπορούν να φανούν μέσα στο
% παράθυρο: ο επόμενος του επιπέδου, που δεν προσομοιώνεται, ανατέλλει τη
% στιγμή που κλείνει το παράθυρο ακόμη και σε ζενιθιακή διέλευση, όπου η
% ορατότητα διαρκεί περισσότερο. Έτσι δεν λείπει από το σενάριο κανένας
% δορυφόρος που θα ήταν ορατός:
%     ((N+1)/2)*Δt - Tvis/2 = W/2   ->   Δt = (W/2 + Tvis/2) / ((N+1)/2)
% Όσο Δt < Tvis, διαδοχικοί δορυφόροι επικαλύπτονται στις ψηλές διελεύσεις
% (δορυφορικό handover χωρίς κενό), ενώ στις χαμηλές μένει κενό κάλυψης.
halfVisS = zenithHalfVisibilityS(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);
if opts.numSats > 1
    satSpacingS = (opts.passWindowS/2 + halfVisS) / ((opts.numSats+1)/2);
else
    satSpacingS = 0;
end
satUOffsetRad = ((opts.numSats+1)/2 - (1:opts.numSats)') * meanMotion*satSpacingS;

% Στρωματοποίηση της μετατόπισης: το [-1, 1] χωρίζεται σε numPasses ίσα
% διαστήματα και κάθε διέλευση παίρνει το μέσο ενός από αυτά. Η κατανομή-στόχος
% είναι η ίδια ομοιόμορφη με το rand(), αλλά με σταθερό πλήθος είναι εγγυημένο
% ότι κάθε περιοχή γεωμετριών, από οριακή ως ζενιθιακή, καλύπτεται ισόποσα.
% Η σειρά ανακατεύεται με δική της γεννήτρια, ώστε η καμπύλη σύγκλισης να μη
% σαρώνει τη γεωμετρία μονότονα και η τυχαιότητα κάθε διέλευσης να μην αλλάζει.
if opts.geometryGrid
    gridStream = RandStream('mt19937ar', 'Seed', opts.rngSeed);
    raanOffsetFrac = -1 + (2*randperm(gridStream, opts.numPasses) - 1) / opts.numPasses;
end

%% ------------------ Δείκτες σύγκλισης ------------------
% Παρακολουθούνται ανά διέλευση (η διέλευση είναι η ανεξάρτητη μονάδα: οι
% χρήστες μέσα στην ίδια διέλευση αλληλεξαρτώνται μέσω της κατανομής πόρων).
kpiNames = {'MeanCapacity_Mbps','MeanSinr_dB','FracServed','FracSatellite','MeanBitPerJouleRf'};
kpiIsFraction = [false false true true false];
% Ο στόχος ακρίβειας δεν περιλαμβάνει το μέσο SINR: είναι μέσος όρος
% λογαριθμικού μεγέθους πάνω σε ανάμεικτο πληθυσμό επίγειων και δορυφορικών
% χρηστών, με εύρος δεκάδων dB, και η εργασία δεν τον αναφέρει ως δείκτη (το
% SINR παρουσιάζεται ως κατανομή). Το διάστημα εμπιστοσύνης του καταγράφεται
% κανονικά στην καμπύλη σύγκλισης.
kpiHasTarget = [true false true true true];
numKpi = numel(kpiNames);
passKpi = nan(opts.numPasses, numKpi);

%% ------------------ Σύγκριση πολιτικών (ITU-R M.2412-0 §7.1) ------------------
% Κάθε βήμα αποτιμάται με όλες τις πολιτικές πάνω στην ΙΔΙΑ πραγματοποίηση
% καναλιού: τα υποψήφια SINR είναι κοινά, αλλάζει μόνο ποιος επιλέγεται και
% άρα ο φόρτος κάθε κόμβου. Οι πολιτικές *MLB είναι οι Actual/TerrestrialOnly
% με εξισορρόπηση φόρτου (TS 38.300 §15.5.1) πάνω στον ίδιο A3. Η σύγκριση είναι επομένως κατά ζεύγη και όχι μεταξύ
% ανεξάρτητων εκτελέσεων, οπότε η διακύμανση της τοπολογίας και του καναλιού
% απαλείφεται από τη διαφορά.
policyNames   = {'Actual','TerrestrialOnly','SatelliteOnly','ActualMLB','TerrestrialOnlyMLB'};
policyMetrics = {'MeanCapacity_Mbps','TotalRate_Mbps','FracServed', ...
                 'FracBelowTarget','FracOutage','BitPerJouleRf','Spectrum_MHz', ...
                 'HoRate_perUserMin','HoRateTerrTerr_perUserMin', ...
                 'HoRateTerrSat_perUserMin','HoRateSatSat_perUserMin', ...
                 'PingPongRate_perUserMin','HoRateForced_perUserMin', ...
                 'HoRateTerrTerrForced_perUserMin','HoRateTerrSatForced_perUserMin', ...
                 'HoRateSatSatForced_perUserMin','FracSatellite'};
numPolicy       = numel(policyNames);
numPolicyMetric = numel(policyMetrics);
policyPassKpi   = nan(opts.numPasses, numPolicy, numPolicyMetric);

convRows = cell(opts.numPasses,1);
checksRun = 0;
totalRows = 0;

fprintf('Ενιαία προσομοίωση ISTN\n');
fprintf('  Δt = %g s | ταχύτητα χρήστη = %g km/h (%.3f m/s) | παράθυρο διέλευσης = %g s\n', ...
    opts.dtSeconds, opts.ueSpeedKmh, ueSpeedMps, opts.passWindowS);
fprintf('  Μετατόπιση ανά βήμα = %.3f m (όριο M.2412 Παρ.1 §5.3.2: 1 m)\n', ueSpeedMps*opts.dtSeconds);
fprintf('  Τροχιακή περίοδος = %.0f s | ταχύτητα ίχνους = %.0f m/s\n', orbitalPeriodS, groundSpeedMps);
fprintf('  Εύρος εγκάρσιας μετατόπισης ίχνους = ±%.2f°\n', rad2deg(raanSpanRad));
fprintf('  Σταθμοί βάσης: %d (εξαγωνικό πλέγμα, ISD %g m) | χρήστες ανά διέλευση: %d-%d, %.0f-%.0f m από τον πλησιέστερο\n', ...
    numBs, opts.bsIsdM, numUsersRange(1), numUsersRange(2), userRMinM, userRMaxM);
fprintf('  Δορυφόροι: %d στο ίδιο επίπεδο, ανά %.1f s (%.2f°, ισοδύναμο %.1f ανά επίπεδο) | ζενιθιακή ορατότητα %.0f s | FRF = %d (%.0f MHz ανά δέσμη)\n', ...
    opts.numSats, satSpacingS, rad2deg(meanMotion*satSpacingS), ...
    ternary(satSpacingS > 0, orbitalPeriodS/max(satSpacingS,eps), NaN), 2*halfVisS, ...
    opts.satFreqReuse, satParameters.Bandwidth/1e6);
fprintf('  Διελεύσεις: %d (σταθερό πλήθος) | γεωμετρία: %s\n', opts.numPasses, ...
    ternary(opts.geometryGrid, 'στρωματοποιημένη', 'τυχαία'));
fprintf('  Στόχος ακρίβειας (ITU-R M.2412-0 §7.1, μόνο αναφορά): 95%% CI < %.0f%% σχετικό / %.0f ποσοστιαίες μονάδες\n', ...
    100*opts.ciTolerance, 100*opts.ciToleranceFrac);
fprintf('  Υστέρηση = %.1f dB | χρόνος επιβεβαίωσης = %d ms (%d βήματα)\n', ...
    opts.hysteresisDb, opts.timeToTriggerMs, ...
    max(1, ceil(opts.timeToTriggerMs/(opts.dtSeconds*1000))));
fprintf('  Μέτρηση: μέσος όρος %d (επίγ.) / %d (δορυφ.) δειγμάτων ανά βήμα | φίλτρο L3 k = %d -> a = %.3f\n', ...
    simParameters.Meas.NAvg, satParameters.Meas.NAvg, opts.l3FilterK, l3Alpha);
fprintf('  Δείκτες με στόχο ακρίβειας: %s\n\n', strjoin(kpiNames(kpiHasTarget), ', '));

%% ------------------ Βρόχος διελεύσεων ------------------
passIdx = 0;
refPassTable = table();
refSnapshot = struct();

% Καταγραφή handovers σε πλήρη ανάλυση 1 s (βλ. logHandoverEvents): ένα κελί
% ανά διέλευση, για τις πολιτικές που ενδιαφέρουν το ML.
hoLogPolicies = {'Actual','ActualMLB'};
hoHistW       = 11;          % βήματα ιστορικού μέτρησης: τρέχον + 10 s πίσω
hoEventCells  = cell(opts.numPasses,1);

if ~isempty(resumeState)
    passIdx     = resumeState.passIdx;
    passKpi(1:passIdx,:) = resumeState.passKpi;
    if isfield(resumeState, 'policyPassKpi') && ~isempty(resumeState.policyPassKpi)
        policyPassKpi(1:passIdx,:,:) = resumeState.policyPassKpi;
    end
    convRows(1:passIdx)  = resumeState.convRows;
    totalRows   = resumeState.totalRows;
    checksRun   = resumeState.checksRun;
    refSnapshot = resumeState.refSnapshot;
    if isfield(resumeState, 'refPassTable')
        refPassTable = resumeState.refPassTable;
    end
    if isfield(resumeState, 'hoEventCells')
        hoEventCells(1:passIdx) = resumeState.hoEventCells;
    end
end

while passIdx < opts.numPasses
    passIdx = passIdx + 1;

    % Ένας σπόρος ανά διέλευση: η τυχαιότητα μέσα στη διέλευση (τοπολογία,
    % κατευθύνσεις, κανάλι) εξελίσσεται συνεχόμενα, ώστε να μη σπάει η
    % χωρική συσχέτιση από επαναρχικοποίηση της γεννήτριας ανά βήμα.
    rng(opts.rngSeed + passIdx);

    % -- Θέσεις χρηστών στην αρχή της διέλευσης --
    % Ομοιόμορφα ως προς το εμβαδόν στην περιοχή όπου η απόσταση από τον
    % πλησιέστερο σταθμό είναι μέσα στο [userRMinM, userRMaxM]. Ισχύει και για
    % τη διέλευση αναφοράς, ώστε να μην αποτελεί ειδική περίπτωση στη στατιστική.
    nU = randi(numUsersRange);
    user_geo = zeros(nU,3);
    for u = 1:nU
        [uLat, uLon] = randPeripheryUser(bsCenterGeo, bs_geo, wgs84, ...
            bsClusterRadiusM + userRMaxM, userRMinM, userRMaxM);
        user_geo(u,:) = [uLat, uLon, ueHeightM];
    end
    numUsers = size(user_geo,1);

    % -- Κατεύθυνση κίνησης: σταθερή ταχύτητα, τυχαία ομοιόμορφη κατεύθυνση --
    % ITU-R M.2412-0 §8.4, ΠΙΝΑΚΑΣ 5 b)/c).
    walkAzimuthDeg = 360*rand(numUsers,1);

    % -- Γεωμετρία της συγκεκριμένης διέλευσης --
    if opts.geometryGrid
        raanRad = raanBase + raanOffsetFrac(passIdx)*raanSpanRad;
    else
        raanRad = raanBase + (2*rand()-1)*raanSpanRad;
    end
    tStart  = tPeakS - opts.passWindowS/2;
    numSteps = floor(opts.passWindowS/opts.dtSeconds) + 1;

    % Προδέσμευση σε απλούς πίνακες αντί για συσσώρευση αντικειμένων table:
    % ένα table ανά βήμα σήμαινε ~900 αντικείμενα ανά διέλευση και εξαντλούσε
    % τη μνήμη σε μεγάλες εκτελέσεις. Ο πίνακας φτιάχνεται μία φορά, στο τέλος.
    nRowsPass = numSteps * numUsers;
    numBuf = zeros(nRowsPass, 25);
    strBuf = strings(nRowsPass, 6);
    rowPtr = 0;

    prevServingNode = strings(numUsers,1);
    channelState = [];   % πρώτο βήμα της διέλευσης: ανεξάρτητο δείγμα
    decisionState = [];  % καμία προϊστορία απόφασης στην αρχή της διέλευσης
    stepBitPerJouleRf = nan(numSteps,1);
    policyStepBuf = nan(numSteps, numPolicy, numPolicyMetric);

    % Ιστορικό μετρήσεων και ακολουθία εξυπηρετούντων ανά πολιτική, για την
    % καταγραφή των handovers. Κυκλικός buffer hoHistW βημάτων.
    hoHist = struct('Bs', nan(hoHistW, numUsers, numBs), ...
        'Sat', nan(hoHistW, numUsers, opts.numSats), ...
        'Elev', nan(hoHistW, numUsers, opts.numSats), 'Dist', nan(numUsers, numBs));
    hoNodeSeq = struct();
    for pL = hoLogPolicies
        hoNodeSeq.(pL{1}) = strings(numSteps, numUsers);
    end
    hoPassNum = zeros(0, 27);
    hoPassStr = strings(0, 5);

    for s = 1:numSteps
        t = tStart + (s-1)*opts.dtSeconds;
        sat_geo = zeros(opts.numSats,3);
        for k = 1:opts.numSats
            sat_geo(k,:) = orbitPositionLla(uPeakRad + meanMotion*(t - tPeakS) + ...
                satUOffsetRad(k), inclDeg, raanRad, a, omegaEarth, t);
        end

        % --- Έλεγχος κατάστασης χρηστών πριν τον υπολογισμό ζεύξης ---
        checksRun = checksRun + 1;
        validateUserState(user_geo, bs_geo, wgs84, ueHeightM, d2dMinM, d2dMaxM, passIdx, s);

        prevDecisionState = decisionState;
        [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
            bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
            nodePowerWattsVec, energyPerBitUJVec, ...
            bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
            satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
            channelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
            bsReasonVec, satReasonVec, policyKpis, decisionState] = ...
            simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, ...
                             satParameters, channelState, prevDecisionState);

        % --- Καταγραφή handovers του βήματος ---
        slot = mod(s-1, hoHistW) + 1;
        hoHist.Bs(slot,:,:)   = channelState.FiltSinrBs_dB;
        hoHist.Sat(slot,:,:)  = channelState.FiltSinrSat_dB;
        hoHist.Elev(slot,:,:) = channelState.SatElev_deg;
        hoHist.Dist           = channelState.BsDist_m;
        for pL = hoLogPolicies
            hoNodeSeq.(pL{1})(s,:) = decisionState.(pL{1}).Node';
        end
        if ~isempty(prevDecisionState)
            [evNum, evStr] = logHandoverEvents(hoLogPolicies, prevDecisionState, ...
                decisionState, hoHist, hoHistW, s, opts.dtSeconds, ...
                passIdx, numUsers, hoPassNum, hoPassStr);
            hoPassNum = [hoPassNum; evNum]; %#ok<AGROW>
            hoPassStr = [hoPassStr; evStr]; %#ok<AGROW>
        end

        % --- Δείκτες των τριών πολιτικών για το βήμα αυτό ---
        for pI = 1:numPolicy
            ps = policyKpis.(policyNames{pI});
            for mI = 1:numPolicyMetric
                policyStepBuf(s, pI, mI) = ps.(policyMetrics{mI});
            end
        end

        % --- Φορτίο κόμβου ---
        nodeLoadVec = zeros(numUsers,1);
        for u = 1:numUsers
            if bestNodeTypeVec(u) == "Outage"
                nodeLoadVec(u) = 0;
            else
                nodeLoadVec(u) = sum(bestNodeVec == bestNodeVec(u));
            end
        end

        % --- Έλεγχος συνέπειας εξόδων ---
        validateOutputs(bestNodeTypeVec, capacityMbpsVec, throughputMbpsVec, ...
            nodeLoadVec, bestSinrDbVec, passIdx, s);

        % --- Κατάσταση ζεύξης και κόστος διακοπής (TR 38.821 §7.3.2.1.1) ---
        linkStateVec      = repmat("Stable", numUsers, 1);
        interruptionMsVec = zeros(numUsers,1);
        for u = 1:numUsers
            curr = bestNodeVec(u);
            prev = prevServingNode(u);
            if curr == "None"
                linkStateVec(u) = "Outage";
            elseif prev ~= "" && prev ~= "None" && curr ~= prev
                linkStateVec(u) = "InTransition";
                interruptionMsVec(u) = 2 * (2*bestDistanceVec(u)/cLight) * 1e3;
            end
        end
        prevServingNode = bestNodeVec;

        lostFraction = min(interruptionMsVec/1e3/opts.dtSeconds, 1);
        deliveredMbpsVec = throughputMbpsVec .* (1 - lostFraction);

        stepBitPerJouleRf(s) = networkEnergy.BitPerJouleRf;

        idx = rowPtr + (1:numUsers);
        numBuf(idx,:) = [ ...
            repmat(passIdx,numUsers,1), repmat(s,numUsers,1), repmat(t-tStart,numUsers,1), ...
            (1:numUsers)', repmat(numUsers,numUsers,1), ...
            user_geo(:,1), user_geo(:,2), walkAzimuthDeg, ...
            bestDistanceVec, bestPathLossVec, bestSinrDbVec, capacityMbpsVec, ...
            bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec, nodeLoadVec, ...
            deliveredMbpsVec, interruptionMsVec, ...
            bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
            satSinrDbVec, satElevationVec, satSlantRangeVec, satPathLossVec];
        strBuf(idx,:) = [bestNodeVec, bestNodeTypeVec, serviceStateVec, ...
            bsReasonVec, satReasonVec, linkStateVec];
        rowPtr = rowPtr + numUsers;

        % --- Στιγμιότυπο μέγιστης προσέγγισης της διέλευσης αναφοράς ---
        if passIdx == 1 && abs(t - tPeakS) <= opts.dtSeconds/2
            refSnapshot = struct('user_geo',user_geo,'sat_geo',sat_geo, ...
                'bestNodeVec',bestNodeVec,'bestNodeTypeVec',bestNodeTypeVec, ...
                'bestDistanceVec',bestDistanceVec,'bestPathLossVec',bestPathLossVec, ...
                'bestSinrDbVec',bestSinrDbVec,'capacityMbpsVec',capacityMbpsVec, ...
                'bestElevationDegVec',bestElevationDegVec, ...
                'nodePowerWattsVec',nodePowerWattsVec,'energyPerBitUJVec',energyPerBitUJVec);
        end

        % --- Κίνηση χρηστών προς το επόμενο βήμα ---
        if s < numSteps
            stepM = ueSpeedMps * opts.dtSeconds;
            [newLat, newLon] = reckon(user_geo(:,1), user_geo(:,2), ...
                repmat(stepM,numUsers,1), walkAzimuthDeg, wgs84);
            user_geo(:,1) = newLat;
            user_geo(:,2) = newLon;
            % Το ύψος δεν μεταβάλλεται ποτέ: παραμένει το φυσικό ύψος κεραίας.
            user_geo(:,3) = ueHeightM;
        end
    end

    % -- Δείκτες της διέλευσης: υπολογίζονται σε ΠΛΗΡΗ ανάλυση 1 s --
    % (η υποδειγματοληψία αφορά μόνο το τι γράφεται στο αρχείο)
    servedFull = strBuf(:,2) ~= "Outage";
    passKpi(passIdx,1) = mean(numBuf(servedFull,12), 'omitnan');   % Capacity_Mbps
    passKpi(passIdx,2) = mean(numBuf(servedFull,11), 'omitnan');   % SINR_dB
    passKpi(passIdx,3) = mean(strBuf(:,3) == "Served");
    passKpi(passIdx,4) = mean(strBuf(:,2) == "Satellite");
    passKpi(passIdx,5) = mean(stepBitPerJouleRf, 'omitnan');

    % Οι δείκτες κάθε πολιτικής συναθροίζονται πρώτα μέσα στη διέλευση και
    % μετά η διέλευση μπαίνει ως ΜΙΑ παρατήρηση στη στατιστική: τα βήματα της
    % ίδιας διέλευσης δεν είναι ανεξάρτητα μεταξύ τους.
    policyPassKpi(passIdx,:,:) = mean(policyStepBuf, 1, 'omitnan');

    % -- Handovers της διέλευσης: τι ακολούθησε το καθένα --
    hoEventCells{passIdx} = finishHandoverEvents(hoPassNum, hoPassStr, ...
        hoNodeSeq, opts.dtSeconds);

    % -- Εγγραφή: κρατούνται τα βήματα 1, 1+stride, ... --
    keep = mod(numBuf(:,2) - 1, opts.datasetStride) == 0;
    passTable = buildPassTable(numBuf(keep,:), strBuf(keep,:));
    if passIdx == 1
        writetable(passTable, opts.datasetPath);
        % Η διέλευση αναφοράς κρατείται και σε πλήρη ανάλυση, για τις
        % χρονοσειρές: εκεί χρειάζεται κάθε βήμα, όχι δείγμα. Ταξιδεύει μέσα
        % στο σημείο ελέγχου, ώστε να επιβιώνει μιας διακοπής.
        refPassTable = buildPassTable(numBuf, strBuf);
    else
        writetable(passTable, opts.datasetPath, 'WriteMode', 'append');
    end
    totalRows = totalRows + height(passTable);
    clear numBuf strBuf passTable;

    % -- Ακρίβεια: 95% CI του μέσου όρου πάνω στις διελεύσεις ως τώρα --
    % Καταγράφεται σε κάθε διέλευση για την καμπύλη σύγκλισης· δεν σταματά
    % την εκτέλεση.
    relHw = nan(1,numKpi);
    for k = 1:numKpi
        x = passKpi(1:passIdx,k);
        x = x(isfinite(x));
        if numel(x) < 2
            relHw(k) = Inf;
            continue;
        end
        hw = 1.96 * std(x) / sqrt(numel(x));
        if kpiIsFraction(k)
            relHw(k) = hw;                       % απόλυτο ημιεύρος για ποσοστά
        else
            relHw(k) = hw / max(abs(mean(x)), eps);
        end
    end

    convRows{passIdx} = [passIdx, totalRows, mean(passKpi(1:passIdx,:),1,'omitnan'), relHw];

    if mod(passIdx,10) == 0 || passIdx <= 3 || passIdx == opts.numPasses
        fprintf('  Διέλευση %3d/%d | γραμμές %8d | χωρητ. %6.2f Mbps | served %5.1f%% | sat %5.1f%% | max rel.CI %.4f | μνήμη %5.2f GB\n', ...
            passIdx, opts.numPasses, totalRows, mean(passKpi(1:passIdx,1),'omitnan'), ...
            100*mean(passKpi(1:passIdx,3),'omitnan'), ...
            100*mean(passKpi(1:passIdx,4),'omitnan'), max(relHw(kpiHasTarget)), ...
            matlabMemoryGb());
    end

    % Σημείο ελέγχου μετά από κάθε διέλευση: λίγα KB, ώστε μια διακοπή να
    % κοστίζει το πολύ μία διέλευση αντί για ολόκληρη την εκτέλεση.
    ckpt = struct('passIdx', passIdx, 'passKpi', passKpi(1:passIdx,:), ...
        'policyPassKpi', policyPassKpi(1:passIdx,:,:), ...
        'convRows', {convRows(1:passIdx)}, 'totalRows', totalRows, ...
        'checksRun', checksRun, 'refSnapshot', refSnapshot, ...
        'refPassTable', refPassTable, 'hoEventCells', {hoEventCells(1:passIdx)}, ...
        'opts', opts);
    save(checkpointPath, 'ckpt');
end

numPasses = passIdx;

%% ------------------ Σύγκριση πολιτικών κατά ζεύγη ------------------
% Κάθε διέλευση δίνει τρεις τιμές για τον ίδιο δείκτη, μία ανά πολιτική, πάνω
% στο ίδιο κανάλι και την ίδια τοπολογία. Η στατιστική γίνεται στη ΔΙΑΦΟΡΑ ανά
% διέλευση, όχι στους δύο μέσους όρους χωριστά.
pol = policyPassKpi(1:numPasses,:,:);

% Φασματική απόδοση: παράγωγο μέγεθος, υπολογισμένο σε επίπεδο ΔΙΕΛΕΥΣΗΣ ως
% λόγος των δύο μέσων και όχι ως μέσος των λόγων. Ο λόγος είναι ότι στα βήματα
% χωρίς ορατό δορυφόρο το φάσμα της αποκλειστικά δορυφορικής πολιτικής είναι
% μηδέν, οπότε ο ανά βήμα λόγος θα ήταν 0/0. Απομονώνει τη συνεισφορά της
% απόφασης από τη συνεισφορά του διαθέσιμου φάσματος (βλ. αξιολόγηση #3).
iRate = find(strcmp(policyMetrics, 'TotalRate_Mbps'), 1);
iSpec = find(strcmp(policyMetrics, 'Spectrum_MHz'), 1);
specEff = pol(:,:,iRate) ./ pol(:,:,iSpec);
specEff(~isfinite(specEff)) = NaN;
pol(:,:,end+1) = specEff;
policyMetrics{end+1} = 'SpectralEff_bpsHz';

% Ποσοστό ping-pong επί των handovers, με την ίδια λογική: λόγος των δύο ρυθμών
% της διέλευσης, αφού το πλήθος χρηστών είναι σταθερό μέσα σε αυτήν. Ανά βήμα
% θα ήταν 0/0 σε κάθε βήμα χωρίς handover. NaN σε διέλευση χωρίς handovers.
iHo = find(strcmp(policyMetrics, 'HoRate_perUserMin'), 1);
iPp = find(strcmp(policyMetrics, 'PingPongRate_perUserMin'), 1);
ppFrac = pol(:,:,iPp) ./ pol(:,:,iHo);
ppFrac(~isfinite(ppFrac)) = NaN;
pol(:,:,end+1) = ppFrac;
policyMetrics{end+1} = 'PingPongFrac';
numPolicyMetric = numel(policyMetrics);

polPassID = repelem((1:numPasses)', numPolicy, 1);
polName   = repmat(string(policyNames(:)), numPasses, 1);
polVals   = zeros(numPasses*numPolicy, numPolicyMetric);
r = 0;
for n = 1:numPasses
    for pI = 1:numPolicy
        r = r + 1;
        polVals(r,:) = reshape(pol(n,pI,:), 1, []);
    end
end
policyTable = array2table(polVals, 'VariableNames', policyMetrics);
policyTable = addvars(policyTable, polPassID, polName, 'Before', 1, ...
    'NewVariableNames', {'PassID','Policy'});

% Οι συγκρίσεις κατά ζεύγη (πολιτική, βάση): η υπό εξέταση πολιτική έναντι
% καθεμιάς από τις δύο αποκλειστικές, και τι προσθέτει το MLB στη μεικτή και
% στην επίγεια, καθώς και η μεικτή έναντι της επίγειας όταν έχουν και οι δύο MLB.
compNames = {'Actual','TerrestrialOnly'; 'Actual','SatelliteOnly'; ...
             'ActualMLB','Actual'; 'TerrestrialOnlyMLB','TerrestrialOnly'; ...
             'ActualMLB','TerrestrialOnlyMLB'};
compIdx = zeros(size(compNames));
for c = 1:numel(compNames)
    compIdx(c) = find(strcmp(policyNames, compNames{c}), 1);
end
nComp   = size(compIdx,1);
pairPolicy   = strings(nComp*numPolicyMetric,1);
pairBaseline = strings(nComp*numPolicyMetric,1);
pairMetric   = strings(nComp*numPolicyMetric,1);
pairN        = zeros(nComp*numPolicyMetric,1);
pairMean     = nan(nComp*numPolicyMetric,1);
pairHw       = nan(nComp*numPolicyMetric,1);
pairRelPct   = nan(nComp*numPolicyMetric,1);
pairBaseMean = nan(nComp*numPolicyMetric,1);
pairWinFrac  = nan(nComp*numPolicyMetric,1);
q = 0;
for c = 1:nComp
    for mI = 1:numPolicyMetric
        % Όχι a/b: το a είναι ο μεγάλος ημιάξονας της τροχιάς και καταγράφεται
        % παρακάτω στο runParams.
        valPol  = pol(:,compIdx(c,1),mI);
        valBase = pol(:,compIdx(c,2),mI);
        d = valPol - valBase;
        ok = isfinite(d);
        d = d(ok);
        q = q + 1;
        pairPolicy(q)   = string(policyNames{compIdx(c,1)});
        pairBaseline(q) = string(policyNames{compIdx(c,2)});
        pairMetric(q)   = string(policyMetrics{mI});
        pairN(q)        = numel(d);
        if numel(d) >= 2
            pairMean(q)    = mean(d);
            pairHw(q)      = 1.96*std(d)/sqrt(numel(d));
            baseMean       = mean(valBase(ok), 'omitnan');
            pairBaseMean(q) = baseMean;
            if abs(baseMean) > eps
                pairRelPct(q) = 100*mean(d)/abs(baseMean);
            end
            pairWinFrac(q) = mean(d > 0);
        end
    end
end
pairedTable = table(pairPolicy, pairBaseline, pairMetric, pairN, pairMean, pairHw, ...
    pairMean - pairHw, pairMean + pairHw, pairBaseMean, pairRelPct, pairWinFrac, ...
    'VariableNames', {'Policy','Baseline','Metric','NumPasses','MeanDiff','CI95HalfWidth', ...
                      'CI95Low','CI95High','BaselineMean','RelDiffPct','FracPassesPositive'});

fprintf('\n  Σύγκριση πολιτικών στην ίδια διέλευση (%d διελεύσεις)\n', numPasses);
iPpFrac = find(strcmp(policyMetrics, 'PingPongFrac'), 1);
iForced = find(strcmp(policyMetrics, 'HoRateForced_perUserMin'), 1);
iSatFrac = find(strcmp(policyMetrics, 'FracSatellite'), 1);
fprintf('    %-18s %12s %12s %12s %10s %14s %10s %12s\n', 'πολιτική', 'ρυθμ.[Mbps]', 'εξυπηρ.', 'εκτός', ...
    'δορυφ.', 'HO/χρήστη/min', 'ping-pong', 'αναγκαστικά');
for pI = 1:numPolicy
    hoM = mean(pol(:,pI,iHo),'omitnan');
    fprintf('    %-18s %12.2f %12.4f %12.4f %9.1f%% %14.3f %9.1f%% %11.1f%%\n', policyNames{pI}, ...
        mean(pol(:,pI,2),'omitnan'), mean(pol(:,pI,3),'omitnan'), ...
        mean(pol(:,pI,5),'omitnan'), 100*mean(pol(:,pI,iSatFrac),'omitnan'), hoM, ...
        100*mean(pol(:,pI,iPpFrac),'omitnan'), ...
        100*mean(pol(:,pI,iForced),'omitnan')/hoM);
end
fprintf('\n    Διαφορές ανά διέλευση (πολιτική μείον βάση), 95%% CI:\n');
for q = 1:height(pairedTable)
    signif = ternary(pairedTable.CI95Low(q) > 0 || pairedTable.CI95High(q) < 0, ...
        '', '   [το CI περιέχει το μηδέν]');
    fprintf('    %-18s - %-18s %-18s %+10.4f ± %.4f (%+6.1f%%, θετική στο %4.1f%% των διελεύσεων)%s\n', ...
        pairedTable.Policy(q), pairedTable.Baseline(q), pairedTable.Metric(q), pairedTable.MeanDiff(q), ...
        pairedTable.CI95HalfWidth(q), pairedTable.RelDiffPct(q), ...
        100*pairedTable.FracPassesPositive(q), signif);
end
fprintf('\n');

%% ------------------ Γράφημα σύγκρισης πολιτικών ------------------
polPng = fullfile(opts.outputDir, 'policy_gain.png');
plotPolicyGain(pairedTable, numPasses, polPng);
mlbPng = fullfile(opts.outputDir, 'mlb_gain.png');
plotMlbGain(pairedTable, numPasses, opts.mlbMaxOffsetDb, mlbPng);

hoPng = fullfile(opts.outputDir, 'handover_rates.png');
plotHandoverRates(pol, policyMetrics, policyNames, numPasses, opts.pingPongMtsS, hoPng);


%% ------------------ Καμπύλη σύγκλισης ------------------
convMat = vertcat(convRows{1:numPasses});
convVarNames = [{'PassID','CumulativeRows'}, ...
    strcat('Mean_', kpiNames), strcat('RelCI_', kpiNames)];
convTable = array2table(convMat, 'VariableNames', convVarNames);

%% ------------------ Ακρίβεια στο τέλος της εκτέλεσης ------------------
% Το πλήθος διελεύσεων είναι σταθερό· εδώ ελέγχεται αν ήταν αρκετό. Η
% διέλευση από την οποία ο στόχος ισχύει ΑΔΙΑΛΕΙΠΤΑ ως το τέλος δείχνει πόσο
% περιθώριο αφήνει το numPasses (NaN: ο στόχος δεν πιάστηκε).
tolVec = opts.ciTolerance * ones(1,numKpi);
tolVec(kpiIsFraction) = opts.ciToleranceFrac;
hwMat = convMat(:, end-numKpi+1:end);
withinTarget = all(hwMat(:,kpiHasTarget) < tolVec(kpiHasTarget), 2);
lastHw = hwMat(end,:);
ciTargetMet = withinTarget(end);
ciTargetFromPass = NaN;
if ciTargetMet
    ciTargetFromPass = find(~withinTarget, 1, 'last') + 1;
    if isempty(ciTargetFromPass), ciTargetFromPass = 1; end
end

%% ------------------ Σύνοψη ------------------
simulatedTimeS = numPasses * opts.passWindowS;
summary = struct();
summary.numPasses        = numPasses;
summary.numRows          = totalRows;
summary.simulatedTime_s  = simulatedTimeS;
summary.ciTargetMet      = ciTargetMet;
summary.ciTargetFromPass = ciTargetFromPass;
summary.checksRun        = checksRun;
for k = 1:numKpi
    summary.(kpiNames{k}) = mean(passKpi(1:numPasses,k), 'omitnan');
    summary.(['CI95_' kpiNames{k}]) = lastHw(k);
end

fprintf('\n--- Σύνοψη ---\n');
if ciTargetMet
    fprintf('  Στόχος ακρίβειας: ΕΠΙΤΕΥΧΘΗΚΕ με %d διελεύσεις (εντός στόχου αδιάλειπτα από τη διέλευση %d)\n', ...
        numPasses, ciTargetFromPass);
else
    fprintf('  Στόχος ακρίβειας: ΔΕΝ επιτεύχθηκε με %d διελεύσεις. Μέγιστο CI = %.4f. Αύξησε το numPasses.\n', ...
        numPasses, max(lastHw(kpiHasTarget)));
end
for k = 1:numKpi
    if kpiIsFraction(k)
        fprintf('    %-20s %10.4f  (95%% CI ±%.4f)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), lastHw(k), ternary(kpiHasTarget(k),'',' [χωρίς στόχο]'));
    else
        fprintf('    %-20s %10.4f  (95%% CI ±%.2f%%)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), 100*lastHw(k), ternary(kpiHasTarget(k),'',' [χωρίς στόχο]'));
    end
end
fprintf('  Διελεύσεις: %d | συνολικός χρόνος προσομοίωσης: %.0f s (%.2f h)\n', ...
    numPasses, simulatedTimeS, simulatedTimeS/3600);
fprintf('  Γραμμές συνόλου δεδομένων: %d -> %s\n', totalRows, opts.datasetPath);
fprintf('  Έλεγχοι κατάστασης χρηστών: %d βήματα, όλοι πέρασαν\n', checksRun);

%% ------------------ Πίνακας και γράφημα της διέλευσης αναφοράς ------------------
if ~isempty(fieldnames(refSnapshot))
    fprintf('\n--- Διέλευση αναφοράς, στιγμή μέγιστης ανύψωσης ---\n');
    array(size(refSnapshot.user_geo,1), refSnapshot.bestNodeVec, refSnapshot.bestNodeTypeVec, ...
        refSnapshot.bestDistanceVec, refSnapshot.bestPathLossVec, refSnapshot.bestSinrDbVec, ...
        refSnapshot.capacityMbpsVec, refSnapshot.bestElevationDegVec, ...
        refSnapshot.nodePowerWattsVec, refSnapshot.energyPerBitUJVec);
    visual(bs_geo, refSnapshot.user_geo, refSnapshot.sat_geo, wgs84, numBs, ...
        size(refSnapshot.user_geo,1), refSnapshot.bestNodeTypeVec, refSnapshot.bestNodeVec, ...
        [userRMinM userRMaxM]);
    saveas(gcf, fullfile(opts.outputDir, 'network_3d.png'));
    close(gcf);
end

% Το dataset ΔΕΝ αντιγράφεται στον versioned φάκελο: είναι δεκάδες MB και το
% Results/runs/ παρακολουθείται από το git. Αντ' αυτού καταγράφεται το άθροισμα
% ελέγχου SHA-256 του, που το συνδέει με τη συγκεκριμένη εκτέλεση και το commit.
% Το μόνο csv είναι το dataset· οι αριθμοί πίσω από τα γραφήματα κρατούνται σε
% ένα .mat, ώστε κάθε γράφημα να μπορεί να ξαναφτιαχτεί χωρίς νέα εκτέλεση.
save(resultsMat, 'refPassTable', 'convTable', 'policyTable', 'pairedTable', 'summary');
outFiles = {resultsMat, polPng, mlbPng, hoPng, fullfile(opts.outputDir,'network_3d.png')};

kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SINR_dB','SatElevation_deg'};
kpiLabel = {'Χωρητικότητα (Mbps)','Ενέργεια ανά bit (\muJ/bit)','SINR (dB)','Γωνία ανύψωσης (deg)'};
for k = 1:numel(kpiList)
    fig = figure('Visible','off');
    try
        theme(fig, 'light');
    catch
    end
    hold on;
    for u = 1:max(refPassTable.UserID)
        r = sortrows(refPassTable(refPassTable.UserID == u, :), 'Step');
        plot(r.Time_s, r.(kpiList{k}), 'DisplayName', sprintf('Χρήστης %d', u));
    end
    hold off; grid on;
    xlabel('Χρόνος (s)'); ylabel(kpiLabel{k});
    legend('Location','best');
    title(sprintf('%s κατά τη διέλευση αναφοράς', kpiList{k}), 'Interpreter','none');
    f = fullfile(opts.outputDir, ['temporal_' kpiList{k} '.png']);
    saveas(fig, f); close(fig);
    outFiles{end+1} = f; %#ok<AGROW>
end

% Καμπύλη σύγκλισης. Γραμμές χωρίς σημάδια (300 σημεία ανά καμπύλη) και άξονας
% περιορισμένος στην περιοχή των στόχων: τα πρώτα σημεία φτάνουν το ~0.9 και
% συμπίεζαν την περιοχή που ενδιαφέρει.
fig = figure('Visible','off','Position',[100 100 900 520]);
try
    theme(fig, 'light');
catch
end
plot(convTable.PassID, convTable.(['RelCI_' kpiNames{1}]), '-', 'LineWidth', 1.3, 'DisplayName', kpiNames{1});
hold on;
for k = 2:numKpi
    plot(convTable.PassID, convTable.(['RelCI_' kpiNames{k}]), '-', 'LineWidth', 1.3, 'DisplayName', kpiNames{k});
end
yline(opts.ciTolerance, '--k', 'LineWidth', 1.2, 'DisplayName', 'στόχος (σχετικό)');
yline(opts.ciToleranceFrac, ':k', 'LineWidth', 1.2, 'DisplayName', 'στόχος (ποσοστά)');
hold off; grid on;
ylim([0 max(0.25, 2*opts.ciTolerance)]);
xlabel('Πλήθος διελεύσεων');
ylabel('Ημιεύρος 95% CI (σχετικό· για τα ποσοστά απόλυτο)');
legend('Location','best','Interpreter','none');
title('Σύγκλιση δεικτών κατά ITU-R M.2412-0 §7.1');
f = fullfile(opts.outputDir, 'convergence.png');
saveas(fig, f); close(fig);
outFiles{end+1} = f;

%% ------------------ Versioning ------------------
runParams = struct();
runParams.dtSeconds        = opts.dtSeconds;
runParams.ueSpeedKmh       = opts.ueSpeedKmh;
runParams.ueSpeed_mps      = ueSpeedMps;
runParams.stepDisplacement_m = ueSpeedMps*opts.dtSeconds;
runParams.passWindow_s     = opts.passWindowS;
runParams.numPasses        = numPasses;
runParams.simulatedTime_s  = simulatedTimeS;
runParams.numRows          = totalRows;
runParams.datasetPath      = opts.datasetPath;
runParams.datasetBytes     = dir(opts.datasetPath).bytes;
runParams.datasetSha256    = fileSha256(opts.datasetPath);
saveHandoverEvents(hoEventCells(1:numPasses), opts.handoverEventsPath, struct( ...
    'dtSeconds', opts.dtSeconds, 'passWindowS', opts.passWindowS, ...
    'pingPongMtsS', opts.pingPongMtsS, 'hysteresisDb', opts.hysteresisDb, ...
    'timeToTriggerMs', opts.timeToTriggerMs, 'rngSeed', opts.rngSeed, ...
    'policies', {hoLogPolicies}));
runParams.handoverEventsPath   = opts.handoverEventsPath;
runParams.handoverEvents       = sum(cellfun(@(e) numel(e.PassID), hoEventCells(1:numPasses)));
runParams.handoverEventsSha256 = fileSha256(opts.handoverEventsPath);
runParams.datasetStride    = opts.datasetStride;
runParams.datasetSampling_s = opts.datasetStride*opts.dtSeconds;
runParams.durationRule     = 'σταθερό πλήθος διελεύσεων (numPasses)· CI μόνο αναφορά';
runParams.geometryGrid     = opts.geometryGrid;
runParams.ciTargetMet      = ciTargetMet;
runParams.ciTargetFromPass = ciTargetFromPass;
runParams.ciTolerance      = opts.ciTolerance;
runParams.ciToleranceFrac  = opts.ciToleranceFrac;
runParams.ciTargetKpis     = strjoin(kpiNames(kpiHasTarget), ', ');
runParams.ciHalfWidthFinal = lastHw;
runParams.hysteresis_dB    = opts.hysteresisDb;
runParams.timeToTrigger_ms = opts.timeToTriggerMs;
runParams.pingPongMts_s    = opts.pingPongMtsS;
runParams.timeToTrigger_steps = max(1, ceil(opts.timeToTriggerMs/(opts.dtSeconds*1000)));
runParams.measAveraging    = opts.measAveraging;
runParams.measNAvgTerrestrial = simParameters.Meas.NAvg;
runParams.measNAvgSatellite   = satParameters.Meas.NAvg;
runParams.l3FilterK        = opts.l3FilterK;
runParams.l3FilterAlphaAtX = ternary(opts.l3FilterK > 0, 2^(-opts.l3FilterK/4), 0);
runParams.l3FilterAlphaAtDt = l3Alpha;
runParams.measPeriod_ms    = opts.measPeriodMs;
runParams.policies         = strjoin(policyNames, ', ');
runParams.policyMetrics    = strjoin(policyMetrics, ', ');
runParams.rngSeed          = opts.rngSeed;
runParams.rngScheme        = 'rng(rngSeed+passIdx) στην αρχή κάθε διέλευσης';
runParams.checksRun        = checksRun;
runParams.userRadiusRange_m = [userRMinM userRMaxM];
runParams.maxWalkPerPass_m = maxWalkM;
runParams.d2dRange_m       = [d2dMinM d2dMaxM];
runParams.numUsersRange    = numUsersRange;
runParams.ueHeight_m       = ueHeightM;
runParams.inclination_deg  = inclDeg;
runParams.raanBase_deg     = rad2deg(raanBase);
runParams.raanSpan_deg     = rad2deg(raanSpanRad);
runParams.omegaEarth_rads  = omegaEarth;
runParams.semiMajorAxis_m  = a;
runParams.satAltitude_m    = satAltitude;
runParams.orbitalPeriod_s  = orbitalPeriodS;
runParams.groundSpeed_mps  = groundSpeedMps;
runParams.scenarioType     = simParameters.PathLoss.Scenario;
runParams.bs_geo           = bs_geo;
runParams.bsRings          = opts.bsRings;
runParams.bsIsd_m          = opts.bsIsdM;
runParams.numSats          = opts.numSats;
runParams.satSpacing_s     = satSpacingS;
runParams.satSpacing_deg   = rad2deg(meanMotion*satSpacingS);
runParams.zenithVisibility_s = 2*halfVisS;
runParams.satFreqReuse     = opts.satFreqReuse;
runParams.terrestrial      = struct('CarrierFrequency_Hz', simParameters.CarrierFrequency, ...
    'TxPower_dBm', simParameters.TxPower, 'AntennaGain_dBi', simParameters.AntennaGain, ...
    'NumAntennaElements', simParameters.NumAntennaElements, ...
    'EIRP_dBm', simParameters.EIRP, 'RxNoiseFigure_dB', simParameters.RxNoiseFigure, ...
    'RxAntTemperature_K', simParameters.RxAntTemperature, ...
    'NSizeGrid', simParameters.Carrier.NSizeGrid, ...
    'SubcarrierSpacing_kHz', simParameters.Carrier.SubcarrierSpacing);
runParams.terrestrial.Power = simParameters.Power;
runParams.satellite        = satParameters;
runParams.sources          = struct( ...
    'timeStep',   'ITU-R M.2412-0 Παρ.1 §5.3.2 (<1 m/βήμα) + TR 38.821 §7.3.2.1.4 Πίν. 7.3.2.1.4-1 (6.61 s)', ...
    'duration',   'ITU-R M.2412-0 §7.1 (επαρκές πλήθος drops + αναφορά εύρους διαστημάτων εμπιστοσύνης)· ITU-R M.2514-0 §8.2.4', ...
    'passWindow', 'TR 38.821 Πίν. 4.2-3 NOTE 1 (visibility time of the satellite)', ...
    'ueMobility', 'ITU-R M.2412-0 §8.4 ΠΙΝΑΚΑΣ 5 b)/c) (3 km/h, τυχαία ομοιόμορφη κατεύθυνση)', ...
    'validity',   'TR 38.901 §7.4.1 Πίν. 7.4.1-1 (10 m - 5 km, h_UT 1.5-22.5 m)', ...
    'hysteresis', 'TS 38.331 §5.5.4.4 (Event A3) + §6.3.2 (Hysteresis, TimeToTrigger)', ...
    'pingPong',   '3GPP TR 36.839 (ping-pong: επιστροφή στην αρχική κυψέλη με χρόνο παραμονής < MTS = 1 s)', ...
    'measFilter', 'TS 38.331 §5.5.3.2 (Layer 3 filtering, F_n = (1-a)F_{n-1} + a M_n· NOTE 1 k=0· NOTE 2 λογαριθμικό)', ...
    'measPeriod', 'TS 38.133 Πίν. 9.2.5.2-1 (FR1 χωρίς DRX: >= 200 ms)', ...
    'measAvg',    'Clarke 1968, αυτοσυσχέτιση ισχύος J0(2*pi*fD*tau)^2 -> ισοδύναμο πλήθος ανεξάρτητων δειγμάτων ανά βήμα', ...
    'bsLayout',   'TR 38.901 §7.8 Πίν. 7.8-1 (UMa: εξαγωνικό πλέγμα, ISD 500 m, h_BS 25 m)', ...
    'satReuse',   'TR 38.821 §6.1.1.1 Πίν. 6.1.1.1-5 (FRF Option 1 = 1 / Option 2 = 3, εύρος ανά δέσμη προσαρμοσμένο στο FRF) + Σχ. 6.1.1.1-2 (ομοδιαυλικές δέσμες = παρεμβολή)', ...
    'satSpacing', 'σχεδιαστική επιλογή: οι numSats δορυφόροι είναι ακριβώς όσοι του επιπέδου φαίνονται μέσα στο παράθυρο (βλ. σχόλιο satSpacingS)');

saveRunVersion('runSimulation', runParams, outFiles, opts.label);

% Η εκτέλεση ολοκληρώθηκε: το σημείο ελέγχου δεν χρειάζεται πια.
if isfile(checkpointPath)
    delete(checkpointPath);
end

end

% =====================================================================

function nEff = effectiveSamplesClarke(fDopplerHz, tStepS)
% Ισοδύναμο πλήθος ανεξάρτητων δειγμάτων διαλείψεων μέσα σε ένα βήμα.
%
% Η μέτρηση που βλέπει ο μηχανισμός απόφασης είναι μέση τιμή πάνω σε ένα
% χρονικό διάστημα, όχι στιγμιαίο δείγμα (TS 38.133 Πίν. 9.2.5.2-1). Πόσο
% μειώνεται η διασπορά από αυτή τη μεσοποίηση το καθορίζει η αυτοσυσχέτιση
% του καναλιού. Για διαλείψεις Rayleigh με ομοιόμορφη γωνιακή κατανομή
% σκεδαστών, το μοντέλο του Clarke (Bell Syst. Tech. J. 47(6), 1968) δίνει
% για την ΙΣΧΥ rho(tau) = J0(2*pi*fD*tau)^2. Η διασπορά της μέσης τιμής σε
% διάστημα T είναι τότε
%
%     Var[μ.ό.] / Var[δείγμα] = (2/T^2) * Int_0^T (T - tau) * rho(tau) dtau
%
% και το αντίστροφο αυτού είναι το ισοδύναμο πλήθος ανεξάρτητων δειγμάτων.
%
% Προτιμάται από τον συνήθη τύπο T_c ~ 0.423/fD, που είναι κατασκευή
% εγχειριδίου και όχι αποτέλεσμα του ίδιου του μοντέλου· δίνει επίσης
% μικρότερο πλήθος, δηλαδή πιο συντηρητική εκτίμηση της μεσοποίησης.
if ~isfinite(fDopplerHz) || fDopplerHz <= 0 || tStepS <= 0
    nEff = 1;
    return;
end
rho   = @(tau) besselj(0, 2*pi*fDopplerHz*tau).^2;
ratio = 2/tStepS^2 * integral(@(tau) (tStepS - tau).*rho(tau), 0, tStepS, ...
    'AbsTol', 1e-12, 'RelTol', 1e-9);
nEff = max(1, round(1/ratio));
end
% Τοπικές συναρτήσεις
% =====================================================================

function [evNum, evStr] = logHandoverEvents(pols, prevDec, newDec, H, W, s, dtS, ...
    passIdx, numUsers, pastNum, pastStr)
% Ένα δείγμα ανά handover, με την κατάσταση τη στιγμή της απόφασης. Ο ορισμός
% του handover είναι ο ίδιος με την καταμέτρηση της simulateScenario: αλλαγή
% από έναν εξυπηρετούντα κόμβο σε άλλον (όχι πρώτη σύνδεση, όχι επανασύνδεση
% μετά από outage). Καταγράφονται και τα αναγκαστικά, με σημαία, ώστε η ετικέτα
% να ορίζεται χωρίς κενά στην ακολουθία.
%
% Όλα τα χαρακτηριστικά είναι διαθέσιμα ΤΗ ΣΤΙΓΜΗ της απόφασης: οι μετρήσεις
% είναι οι φιλτραρισμένες (TS 38.331 §5.5.3.2) που βλέπει ο A3, του τρέχοντος
% βήματος και των προηγούμενων· ο φόρτος είναι αυτός που αναφέρθηκε στο
% προηγούμενο βήμα. Τίποτα από το μέλλον, που μπαίνει μόνο στις ετικέτες
% (finishHandoverEvents).
evNum = zeros(0, 27);
evStr = strings(0, 5);
tNow  = (s-1)*dtS;
l3    = round(3/dtS);
l10   = round(10/dtS);
for pI = 1:numel(pols)
    pol = pols{pI};
    pN  = prevDec.(pol).Node;
    sN  = newDec.(pol).Node;
    isHo = pN ~= "" & pN ~= "None" & sN ~= "None" & sN ~= pN;
    for u = find(isHo)'
        f = pN(u);
        g = sN(u);
        fromSat = startsWith(f, "SAT");
        toSat   = startsWith(g, "SAT");
        if ~fromSat && ~toSat
            hoType = "TerrTerr";
        elseif fromSat && toSat
            hoType = "SatSat";
        else
            hoType = "TerrSat";
        end
        prevNode = prevDec.(pol).Prev(u);
        mF = histMeas(H, W, s, 0, u, f);
        mT = histMeas(H, W, s, 0, u, g);
        eF = histElev(H, W, s, 0, u, f);
        eT = histElev(H, W, s, 0, u, g);
        past = pastStr(:,1) == string(pol) & pastNum(:,4) == u & pastNum(:,3) > tNow - 60;
        ldNode = prevDec.(pol).LoadNode;
        ldRho  = prevDec.(pol).LoadRho;
        row = [passIdx, s, tNow, u, numUsers, ...
            newDec.(pol).Forced(u), g == prevNode, ...
            (prevDec.(pol).Stay(u) + 1)*dtS, sum(past), ...
            mF, mT, mT - mF, ...
            mF - histMeas(H, W, s, l3, u, f),  mT - histMeas(H, W, s, l3, u, g), ...
            mF - histMeas(H, W, s, l10, u, f), mT - histMeas(H, W, s, l10, u, g), ...
            eF, eT, (eF - histElev(H, W, s, l3, u, f))/(l3*dtS), ...
            (eT - histElev(H, W, s, l3, u, g))/(l3*dtS), ...
            bsDist(H, u, f), bsDist(H, u, g), ...
            sum(ldRho(ldNode == f)), sum(ldRho(ldNode == g)), ...
            sum(pN == f), sum(sN == g), newDec.(pol).Offset(u)];
        evNum(end+1,:) = row;                                     %#ok<AGROW>
        evStr(end+1,:) = [string(pol), f, g, prevNode, hoType]; %#ok<AGROW>
    end
end
end

function ev = finishHandoverEvents(num, str, nodeSeq, dtS)
% Συμπληρώνει τι ακολούθησε κάθε handover, από την ακολουθία εξυπηρετούντων
% της διέλευσης σε πλήρη ανάλυση:
%   StayTo_s   βήματα στον νέο κόμβο από τη στιγμή του handover ως την αποχώρηση,
%              σε s (ίδια σύμβαση με την παραμονή του ορισμού ping-pong)
%   NextNode   πού πήγε μετά ("None" = outage, "" = τέλος διέλευσης πρώτα)
%   Censored   η διέλευση τελείωσε πριν φύγει από τον νέο κόμβο
% Ping-pong με παράθυρο T: NextNode == FromNode && StayTo_s <= T. Με T = 1 s
% ταυτίζεται με την καταμέτρηση της simulateScenario (MTS, TR 36.839).
n = size(num,1);
stayTo   = nan(n,1);
toEnd    = nan(n,1);
censored = false(n,1);
nextNode = strings(n,1);
for i = 1:n
    seq = nodeSeq.(char(str(i,1)))(:, num(i,4));
    s   = num(i,2);
    nSteps = numel(seq);
    toEnd(i) = (nSteps - s)*dtS;
    j = find(seq(s+1:end) ~= str(i,3), 1);
    if isempty(j)
        stayTo(i)   = (nSteps - s + 1)*dtS;
        censored(i) = true;
    else
        stayTo(i)   = j*dtS;
        nextNode(i) = seq(s+j);
    end
end
names = handoverEventNumNames();
ev = struct();
for c = 1:numel(names)
    ev.(names{c}) = num(:,c);
end
ev.Forced       = logical(ev.Forced);
ev.ReturnToPrev = logical(ev.ReturnToPrev);
ev.Policy   = str(:,1);
ev.FromNode = str(:,2);
ev.ToNode   = str(:,3);
ev.PrevNode = str(:,4);
ev.HoType   = str(:,5);
ev.TimeToPassEnd_s = toEnd;
ev.StayTo_s = stayTo;
ev.Censored = censored;
ev.NextNode = nextNode;
end

function names = handoverEventNumNames()
names = {'PassID','Step','Time_s','UserID','NumUsers','Forced','ReturnToPrev', ...
    'StayFrom_s','HoLast60s','MeasFrom_dB','MeasTo_dB','Margin_dB', ...
    'MeasFromDelta3s_dB','MeasToDelta3s_dB','MeasFromDelta10s_dB','MeasToDelta10s_dB', ...
    'ElevFrom_deg','ElevTo_deg','ElevRateFrom_degps','ElevRateTo_degps', ...
    'DistFrom_m','DistTo_m','LoadFrom','LoadTo','UsersOnFrom','UsersOnTo','MlbOffset_dB'};
end

function saveHandoverEvents(cells, path, meta)
% Ένα αρχείο .mat (-v7) με μία μεταβλητή ανά στήλη, ώστε να διαβάζεται
% απευθείας από την Python (scipy.io.loadmat). Τα κείμενα ως cellstr.
cells = cells(~cellfun(@isempty, cells));
out = struct();
if isempty(cells)
    return;
end
fn = fieldnames(cells{1});
for k = 1:numel(fn)
    col = cellfun(@(e) e.(fn{k}), cells, 'UniformOutput', false);
    col = vertcat(col{:});
    if isstring(col)
        col = cellstr(col);
    end
    out.(fn{k}) = col;
end
mf = fieldnames(meta);
for k = 1:numel(mf)
    out.(['Meta_' mf{k}]) = meta.(mf{k});
end
save(path, '-struct', 'out', '-v7');
end

function v = histMeas(H, W, s, lag, u, node)
% Φιλτραρισμένη μέτρηση του κόμβου lag βήματα πριν· NaN εκτός ιστορικού.
v = NaN;
if lag >= W || s - lag < 1
    return;
end
slot = mod(s - lag - 1, W) + 1;
[isSat, idx] = nodeIndex(node);
if isSat
    v = H.Sat(slot, u, idx);
else
    v = H.Bs(slot, u, idx);
end
end

function v = histElev(H, W, s, lag, u, node)
% Ανύψωση δορυφόρου lag βήματα πριν· NaN για σταθμό βάσης ή εκτός ιστορικού.
v = NaN;
[isSat, idx] = nodeIndex(node);
if ~isSat || lag >= W || s - lag < 1
    return;
end
v = H.Elev(mod(s - lag - 1, W) + 1, u, idx);
end

function v = bsDist(H, u, node)
% Απόσταση 3D από σταθμό βάσης στο τρέχον βήμα· NaN για δορυφόρο.
[isSat, idx] = nodeIndex(node);
if isSat
    v = NaN;
else
    v = H.Dist(u, idx);
end
end

function [isSat, idx] = nodeIndex(node)
isSat = startsWith(node, "SAT-");
if isSat
    idx = str2double(extractAfter(node, "SAT-"));
else
    idx = str2double(extractAfter(node, "BS"));
end
end

function tf = checkpointMatches(a, b)
% Το σημείο ελέγχου χρησιμοποιείται μόνο αν οι παράμετροι που καθορίζουν τι
% παράγει η εκτέλεση είναι ίδιες. Οτιδήποτε άλλο (π.χ. label) δεν πειράζει.
keys = {'dtSeconds','ueSpeedKmh','passWindowS','datasetStride','rngSeed', ...
        'numPasses','geometryGrid','datasetPath', ...
        'hysteresisDb','timeToTriggerMs','l3FilterK','measPeriodMs','measAveraging', ...
        'pingPongMtsS','numUsersRange','bsRings','bsIsdM','numSats','satFreqReuse', ...
        'mlbMaxOffsetDb','handoverEventsPath'};
tf = true;
for k = 1:numel(keys)
    if ~isfield(a,keys{k}) || ~isfield(b,keys{k}) || ~isequal(a.(keys{k}), b.(keys{k}))
        tf = false;
        return;
    end
end
end

function h = fileSha256(path)
% SHA-256 ενός αρχείου, ως συμβολοσειρά δεκαεξαδικών. Συνδέει το παραγόμενο
% σύνολο δεδομένων με τη συγκεκριμένη εκτέλεση και το commit που καταγράφει το
% saveRunVersion, χωρίς να χρειάζεται αντίγραφο του αρχείου στον φάκελο.
try
    md = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(path, 'r');
    if fid < 0
        h = '';
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    while true
        chunk = fread(fid, 1e7, '*uint8');
        if isempty(chunk)
            break;
        end
        md.update(chunk);
    end
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8')).', 1, []));
catch
    % MATLAB χωρίς JVM (π.χ. όταν η Java του συστήματος είναι headless): το
    % εργαλείο του λειτουργικού δίνει το ίδιο άθροισμα.
    h = systemSha256(path);
end
end

function h = systemSha256(path)
h = '';
if ispc
    cmd = sprintf('certutil -hashfile "%s" SHA256', path);
elseif ismac
    cmd = sprintf('shasum -a 256 "%s"', path);
else
    cmd = sprintf('sha256sum "%s"', path);
end
[st, out] = system(cmd);
if st ~= 0
    return;
end
% Το certutil σε παλαιότερα Windows χωρίζει τα bytes με κενά.
tok = regexp(strrep(out, ' ', ''), '[0-9a-fA-F]{64}', 'match', 'once');
if ~isempty(tok)
    h = lower(tok);
end
end

function g = matlabMemoryGb()
% Μνήμη που κρατά η MATLAB, σε GB. Χρήσιμο σε μεγάλες εκτελέσεις: αν ο
% αριθμός ανεβαίνει σταθερά ανά διέλευση, κάτι συσσωρεύεται και δεν
% ελευθερώνεται. Η συνάρτηση memory υπάρχει μόνο σε Windows.
try
    m = memory;
    g = m.MemUsedMATLAB / 2^30;
catch
    g = NaN;
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function plotPolicyGain(pairedTable, numPasses, pngPath)
% Πόσο καλύτερη είναι η υπό εξέταση πολιτική από καθεμία από τις δύο
% αποκλειστικές: μέση διαφορά κατά ζεύγη ανά διέλευση, με 95% CI. Θετική τιμή
% σημαίνει ΠΑΝΤΑ ότι η υπό εξέταση πολιτική είναι καλύτερη· για το ποσοστό
% εκτός κάλυψης το πρόσημο αντιστρέφεται ώστε να ισχύει το ίδιο.
% Ρυθμοί και ενέργεια ως σχετική μεταβολή ως προς τη βάση. Τα ποσοστά χρηστών
% σε ποσοστιαίες μονάδες: η σχετική μεταβολή ενός ποσοστού φουσκώνει όταν η
% βάση είναι μικρή, όπως το outage της αποκλειστικά επίγειας πολιτικής.
% Το φάσμα δεν είναι δείκτης απόδοσης αλλά διαθέσιμος πόρος, οπότε δεν
% σχεδιάζεται ως «βελτίωση»· αναγράφεται στον τίτλο ως πλαίσιο, γιατί χωρίς
% αυτό η αύξηση της ρυθμαπόδοσης διαβάζεται ως κέρδος της απόφασης.
baselines = ["TerrestrialOnly","SatelliteOnly"];
baseLabel = {'έναντι μόνο επίγειας','έναντι μόνο δορυφορικής'};
relMetric = ["MeanCapacity_Mbps","TotalRate_Mbps","BitPerJouleRf","SpectralEff_bpsHz"];
relLabel  = {'Μέση χωρητικότητα χρήστη','Συνολική ρυθμαπόδοση', ...
             'Ενεργειακή απόδοση (bit/J)','Φασματική απόδοση'};
ppMetric  = ["FracServed","FracOutage"];
ppSign    = [1 -1];
ppLabel   = {'Αύξηση εξυπηρετούμενων','Μείωση εκτός κάλυψης'};

T = pairedTable(pairedTable.Policy == "Actual", :);
[relVal, relHw] = gainMatrix(T, baselines, relMetric, ones(1,numel(relMetric)), true);
[ppVal,  ppHw]  = gainMatrix(T, baselines, ppMetric, ppSign, false);
specVal         = gainMatrix(T, baselines, "Spectrum_MHz", 1, true);

fig = figure('Visible','off','Position',[100 100 1400 520]);
try
    theme(fig, 'light');   % λευκό φόντο για το κείμενο της εργασίας (R2025a+)
catch
end
subplot(1,2,1);
gainPanel(relVal, relHw, relLabel, baseLabel, '%');
xlabel('Βελτίωση [%], σχετική ως προς τη βάση');
title('Ρυθμός και ενέργεια');

subplot(1,2,2);
gainPanel(ppVal, ppHw, ppLabel, baseLabel, ' μον.');
xlabel('Βελτίωση [ποσοστιαίες μονάδες]');
title('Κάλυψη χρηστών');

sgtitle(sprintf(['Υπό εξέταση πολιτική έναντι των αποκλειστικών: %d διελεύσεις, 95%% CI ' ...
    '(ns: το CI περιέχει το 0)\nΔιαθέσιμο φάσμα: %+.0f%% έναντι μόνο επίγειας, ' ...
    '%+.0f%% έναντι μόνο δορυφορικής'], numPasses, specVal(1), specVal(2)));
saveas(fig, pngPath);
close(fig);
end

function plotMlbGain(pairedTable, numPasses, mlbOffsetDb, pngPath)
% Τι προσθέτει η εξισορρόπηση φόρτου (MLB): η ίδια πολιτική με και χωρίς MLB,
% στο ίδιο κανάλι, και η μεικτή έναντι της επίγειας όταν έχουν και οι δύο MLB.
% Ίδια σύμβαση με το policy_gain: θετική τιμή = καλύτερη η πρώτη πολιτική.
% Στα handovers και στα ping-pong το πρόσημο αντιστρέφεται (λιγότερα = καλύτερα).
pairs = ["ActualMLB","Actual"; "TerrestrialOnlyMLB","TerrestrialOnly"; ...
         "ActualMLB","TerrestrialOnlyMLB"];
pairLabel = {'μεικτή: με MLB έναντι χωρίς','επίγεια: με MLB έναντι χωρίς', ...
             'με MLB: μεικτή έναντι επίγειας'};
relMetric = ["MeanCapacity_Mbps","TotalRate_Mbps","BitPerJouleRf","HoRate_perUserMin"];
relSign   = [1 1 1 -1];
relLabel  = {'Μέση χωρητικότητα χρήστη','Συνολική ρυθμαπόδοση', ...
             'Ενεργειακή απόδοση (bit/J)','Μείωση handovers'};
ppMetric  = ["FracServed","FracOutage","PingPongFrac"];
ppSign    = [1 -1 -1];
ppLabel   = {'Αύξηση εξυπηρετούμενων','Μείωση εκτός κάλυψης','Μείωση ping-pong (μερίδιο HO)'};

nP = size(pairs,1);
relVal = nan(numel(relMetric), nP); relHw = relVal;
ppVal  = nan(numel(ppMetric),  nP); ppHw  = ppVal;
for j = 1:nP
    T = pairedTable(pairedTable.Policy == pairs(j,1), :);
    [relVal(:,j), relHw(:,j)] = gainMatrix(T, pairs(j,2), relMetric, relSign, true);
    [ppVal(:,j),  ppHw(:,j)]  = gainMatrix(T, pairs(j,2), ppMetric,  ppSign,  false);
end

fig = figure('Visible','off','Position',[100 100 1400 560]);
try
    theme(fig, 'light');
catch
end
subplot(1,2,1);
gainPanel(relVal, relHw, relLabel, pairLabel, '%');
xlabel('Βελτίωση [%], σχετική ως προς τη βάση');
title('Ρυθμός, ενέργεια, handovers');

subplot(1,2,2);
gainPanel(ppVal, ppHw, ppLabel, pairLabel, ' μον.');
xlabel('Βελτίωση [ποσοστιαίες μονάδες]');
title('Κάλυψη χρηστών και ping-pong');

sgtitle(sprintf(['Εξισορρόπηση φόρτου (MLB, TS 38.300 §15.5.1): %d διελεύσεις, 95%% CI ' ...
    '(ns: το CI περιέχει το 0)\nμετατόπιση A3 = %g dB \\times διαφορά εικονικού φόρτου ' ...
    '(Nasri & Altman), φόρτος κατά Awada et al.'], numPasses, mlbOffsetDb));
saveas(fig, pngPath);
close(fig);
end

function plotHandoverRates(pol, policyMetrics, policyNames, numPasses, mtsS, pngPath)
% Ρυθμός handovers ανά πολιτική, σε απόλυτες τιμές και όχι ως «βελτίωση»: η
% αποκλειστικά δορυφορική πολιτική με έναν δορυφόρο δεν έχει κανένα handover,
% οπότε η σχετική διαφορά από αυτήν δεν ορίζεται. Στοιβαγμένα ανά τύπο και, μέσα
% σε κάθε τύπο, σε αποφάσεις του A3 (σκούρο) και αναγκαστικά (ανοιχτό, απώλεια
% του εξυπηρετούντος κόμβου). Το 95% CI αφορά το σύνολο, πάνω στις διελεύσεις.
% Τύποι που είναι μηδενικοί σε όλες τις πολιτικές παραλείπονται από το υπόμνημα
% (π.χ. δορυφορικό-δορυφορικό με έναν δορυφόρο).
typeKey   = {'TerrTerr','TerrSat','SatSat'};
typeLabel = {'επίγειο \leftrightarrow επίγειο','επίγειο \leftrightarrow δορυφορικό', ...
             'δορυφορικό \leftrightarrow δορυφορικό'};
polLabel  = {'υπό εξέταση','μόνο επίγεια','μόνο δορυφορική','υπό εξέταση + MLB','μόνο επίγεια + MLB'};
cols      = [0.20 0.45 0.75; 0.90 0.55 0.15; 0.45 0.65 0.30];
lightMix  = 0.55;   % ανάμειξη με λευκό για το αναγκαστικό τμήμα
numPolicy = numel(policyNames);

idx    = @(name) find(strcmp(policyMetrics, name), 1);
polMean = @(name) mean(pol(:,:,idx(name)), 1, 'omitnan')';
stackV  = zeros(numPolicy, 0);
segCol  = zeros(0, 3);
segLbl  = {};
for t = 1:numel(typeKey)
    allV    = polMean(['HoRate' typeKey{t} '_perUserMin']);
    forcedV = polMean(['HoRate' typeKey{t} 'Forced_perUserMin']);
    if ~any(allV > 0)
        continue;
    end
    stackV(:,end+1) = allV - forcedV; %#ok<AGROW>
    stackV(:,end+1) = forcedV;        %#ok<AGROW>
    segCol(end+1,:) = cols(t,:);                              %#ok<AGROW>
    segCol(end+1,:) = cols(t,:) + lightMix*(1 - cols(t,:));   %#ok<AGROW>
    segLbl{end+1} = [typeLabel{t} ' (A3)'];          %#ok<AGROW>
    segLbl{end+1} = [typeLabel{t} ' (αναγκαστικό)']; %#ok<AGROW>
end
totV    = pol(:,:,idx('HoRate_perUserMin'));
totM    = mean(totV, 1, 'omitnan');
totHw   = 1.96 * std(totV, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(totV), 1));
forcedM = mean(pol(:,:,idx('HoRateForced_perUserMin')), 1, 'omitnan');
ppM     = 100 * mean(pol(:,:,idx('PingPongFrac')), 1, 'omitnan');

fig = figure('Visible','off','Position',[100 100 1150 560]);
try
    theme(fig, 'light');
catch
end
if isempty(stackV)
    stackV = zeros(numPolicy, 1);   % κανένα handover σε καμία πολιτική
end
b = bar(stackV, 'stacked');
for k = 1:min(numel(b), size(segCol,1))
    b(k).FaceColor = segCol(k,:);
end
hold on;
errorbar(1:numPolicy, totM, totHw, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
top = max(totM + totHw, [], 'omitnan');
if ~(top > 0), top = 1; end
for pI = 1:numPolicy
    if totM(pI) > 0
        lbl = sprintf('%.3f\nαναγκαστικά %.0f%%\nping-pong %.1f%%', totM(pI), ...
            100*forcedM(pI)/totM(pI), ppM(pI));
    else
        lbl = sprintf('%.3f\nκανένα handover', totM(pI));
    end
    text(pI, totM(pI) + max(totHw(pI), 0) + 0.04*top, lbl, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
end
hold off;
ylim([0 1.45*top]);
set(gca, 'XTick', 1:numPolicy, 'XTickLabel', polLabel);
ylabel('Handovers ανά χρήστη ανά λεπτό');
if ~isempty(segLbl)
    legend(segLbl, 'Location', 'northoutside', 'NumColumns', 2);
end
grid on;
title(sprintf(['Ρυθμός handovers ανά πολιτική: %d διελεύσεις, 95%% CI\n' ...
    'αναγκαστικό = απώλεια του εξυπηρετούντος κόμβου· ping-pong = επιστροφή με παραμονή \\leq %g s'], ...
    numPasses, mtsS));
saveas(fig, pngPath);
close(fig);
end

function [val, hw] = gainMatrix(T, baselines, metrics, sgn, relative)
% Πίνακας (δείκτης x βάση) με τη μέση διαφορά και το ημιεύρος του CI της,
% είτε σχετικά (%) είτε σε ποσοστιαίες μονάδες.
val = nan(numel(metrics), numel(baselines));
hw  = nan(numel(metrics), numel(baselines));
for i = 1:numel(metrics)
    for j = 1:numel(baselines)
        r = find(T.Baseline == baselines(j) & T.Metric == metrics(i), 1);
        if isempty(r) || ~isfinite(T.MeanDiff(r))
            continue;
        end
        scale = 100;
        if relative
            base = abs(T.BaselineMean(r));
            if ~(base > eps)
                continue;   % χωρίς βάση η σχετική μεταβολή δεν ορίζεται
            end
            scale = 100 / base;
        end
        val(i,j) = sgn(i) * T.MeanDiff(r) * scale;
        hw(i,j)  = T.CI95HalfWidth(r) * scale;
    end
end
end

function gainPanel(val, hw, labels, seriesLabels, unit)
% Οριζόντιες ομαδοποιημένες μπάρες, μία ομάδα ανά δείκτη και μία μπάρα ανά βάση.
nM   = size(val,1);
cols = [0.20 0.45 0.75; 0.90 0.55 0.15; 0.45 0.65 0.30];
span = max([abs(val(:)) + hw(:); 1], [], 'omitnan');
pad  = 0.03 * span;

b = barh(val, 'grouped');
hold on;
for j = 1:numel(b)
    b(j).FaceColor = cols(j,:);
    y = b(j).XEndPoints;
    x = b(j).YEndPoints;
    errorbar(x, y, [], [], hw(:,j)', hw(:,j)', 'k', 'LineStyle', 'none', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    for i = 1:nM
        if ~isfinite(x(i))
            continue;
        end
        txt = sprintf('%+.1f%s', x(i), unit);
        if abs(x(i)) <= hw(i,j)
            txt = [txt ' (ns)']; %#ok<AGROW>
        end
        if x(i) >= 0
            text(x(i) + hw(i,j) + pad, y(i), txt, 'HorizontalAlignment', 'left', 'FontSize', 8);
        else
            text(x(i) - hw(i,j) - pad, y(i), txt, 'HorizontalAlignment', 'right', 'FontSize', 8);
        end
    end
end
xline(0, 'k-', 'HandleVisibility', 'off');
hold off;

lo = min([0; val(:) - hw(:)], [], 'omitnan');
hi = max([0; val(:) + hw(:)], [], 'omitnan');
xlim([lo - (lo < 0)*0.35*span - pad, hi + 0.35*span]);
set(gca, 'YTick', 1:nM, 'YTickLabel', labels, 'YDir', 'reverse');
grid on;
legend(seriesLabels, 'Location', 'southoutside', 'Orientation', 'horizontal');
end

function T = buildPassTable(numBuf, strBuf)
% Φτιάχνει τον πίνακα μιας διέλευσης από τους δύο προδεσμευμένους πίνακες.
% Οι στήλες του numBuf είναι, με τη σειρά: PassID, Step, Time_s, UserID,
% NumUsers, UserLat, UserLon, WalkAzimuth_deg, Distance_m, PathLoss_dB,
% SINR_dB, Capacity_Mbps, SatElevation_deg, NodePower_W, EnergyPerBit_uJ,
% NodeLoad, Throughput_Mbps, Interruption_ms, CandBS_SINR_dB,
% CandBS_Distance_m, CandBS_PathLoss_dB, CandSat_SINR_dB,
% CandSat_Elevation_deg, CandSat_SlantRange_m, CandSat_PathLoss_dB.
% Του strBuf: ServingNode, ServingType, ServiceState, BsUnavailReason,
% SatUnavailReason, LinkState.
T = table(numBuf(:,1), numBuf(:,2), numBuf(:,3), numBuf(:,4), numBuf(:,5), ...
    numBuf(:,6), numBuf(:,7), numBuf(:,8), ...
    strBuf(:,1), strBuf(:,2), ...
    numBuf(:,9), numBuf(:,10), numBuf(:,11), numBuf(:,12), numBuf(:,13), ...
    numBuf(:,14), numBuf(:,15), numBuf(:,16), numBuf(:,17), ...
    strBuf(:,3), strBuf(:,4), strBuf(:,5), strBuf(:,6), ...
    numBuf(:,18), ...
    numBuf(:,19), numBuf(:,20), numBuf(:,21), ...
    numBuf(:,22), numBuf(:,23), numBuf(:,24), numBuf(:,25), ...
    'VariableNames', {'PassID','Step','Time_s','UserID','NumUsers', ...
    'UserLat','UserLon','WalkAzimuth_deg', ...
    'ServingNode','ServingType','Distance_m','PathLoss_dB', ...
    'SINR_dB','Capacity_Mbps','SatElevation_deg', ...
    'NodePower_W','EnergyPerBit_uJ','NodeLoad', ...
    'Throughput_Mbps','ServiceState','BsUnavailReason','SatUnavailReason', ...
    'LinkState','Interruption_ms', ...
    'CandBS_SINR_dB','CandBS_Distance_m','CandBS_PathLoss_dB', ...
    'CandSat_SINR_dB','CandSat_Elevation_deg','CandSat_SlantRange_m','CandSat_PathLoss_dB'});
end

function validateUserState(user_geo, bs_geo, wgs84, ueHeightM, d2dMin, d2dMax, passIdx, stepIdx)
% Ελέγχει ότι σε κάθε βήμα οι μεταβλητές των χρηστών είναι οι αναμενόμενες:
% ύψος κεραίας αμετάβλητο, πεπερασμένες συντεταγμένες, και τουλάχιστον μία
% επίγεια ζεύξη εντός του πεδίου ισχύος των UMa/UMi (TR 38.901 Πίν. 7.4.1-1).
if any(user_geo(:,3) ~= ueHeightM)
    error('runSimulation:HeightDrift', ...
        'Διέλευση %d, βήμα %d: ύψος χρήστη %.4f m αντί για %.4f m.', ...
        passIdx, stepIdx, max(abs(user_geo(:,3))), ueHeightM);
end
if ~all(isfinite(user_geo(:)))
    error('runSimulation:NonFiniteGeo', ...
        'Διέλευση %d, βήμα %d: μη πεπερασμένη συντεταγμένη χρήστη.', passIdx, stepIdx);
end
if any(abs(user_geo(:,1)) > 90) || any(abs(user_geo(:,2)) > 180)
    error('runSimulation:GeoOutOfBounds', ...
        'Διέλευση %d, βήμα %d: συντεταγμένη χρήστη εκτός ορίων.', passIdx, stepIdx);
end

numUsers = size(user_geo,1);
for u = 1:numUsers
    d2d = distance(user_geo(u,1), user_geo(u,2), bs_geo(:,1), bs_geo(:,2), wgs84);
    if ~any(d2d >= d2dMin & d2d <= d2dMax)
        error('runSimulation:NoValidLink', ...
            ['Διέλευση %d, βήμα %d, χρήστης %d: καμία επίγεια ζεύξη εντός ' ...
             '[%d, %d] m (ελάχιστη απόσταση %.1f m, μέγιστη %.1f m).'], ...
            passIdx, stepIdx, u, d2dMin, d2dMax, min(d2d), max(d2d));
    end
end
end

function validateOutputs(typeVec, capacityVec, throughputVec, nodeLoadVec, sinrVec, passIdx, stepIdx)
% Ελέγχει τη συνέπεια των εξόδων: outage -> καμία χωρητικότητα και μηδενική
% ρυθμαπόδοση· εξυπηρετούμενος -> πεπερασμένη χωρητικότητα και φορτίο >= 1.
isOut = typeVec == "Outage";
if any(isOut & (~isnan(capacityVec) | throughputVec ~= 0 | nodeLoadVec ~= 0))
    error('runSimulation:OutageInconsistent', ...
        'Διέλευση %d, βήμα %d: χρήστης σε outage με μη μηδενική χωρητικότητα/ρυθμαπόδοση/φορτίο.', ...
        passIdx, stepIdx);
end
if any(~isOut & (~isfinite(capacityVec) | capacityVec <= 0))
    error('runSimulation:ServedInconsistent', ...
        'Διέλευση %d, βήμα %d: εξυπηρετούμενος χρήστης χωρίς πεπερασμένη θετική χωρητικότητα.', ...
        passIdx, stepIdx);
end
if any(~isOut & nodeLoadVec < 1)
    error('runSimulation:LoadInconsistent', ...
        'Διέλευση %d, βήμα %d: εξυπηρετούμενος χρήστης με φορτίο κόμβου < 1.', ...
        passIdx, stepIdx);
end
if any(~isOut & ~isfinite(sinrVec))
    error('runSimulation:SinrInconsistent', ...
        'Διέλευση %d, βήμα %d: εξυπηρετούμενος χρήστης χωρίς πεπερασμένο SINR.', ...
        passIdx, stepIdx);
end
end

function [lat, lon] = randPeripheryUser(centerGeo, bs_geo, wgs84, rDiskM, rMinM, rMaxM)
% Θέση χρήστη ομοιόμορφη ως προς το εμβαδόν, με απόσταση από τον πλησιέστερο
% σταθμό μέσα στο [rMinM, rMaxM]. Δειγματοληψία με απόρριψη από δίσκο ακτίνας
% rDiskM γύρω από το κέντρο του πλέγματος, που περιέχει όλη την επιτρεπτή
% περιοχή: η ομοιόμορφη κατανομή στον δίσκο, περιορισμένη στην περιοχή,
% είναι ομοιόμορφη στην περιοχή.
while true
    rM = rDiskM * sqrt(rand());
    az = 360 * rand();
    [lat, lon] = reckon(centerGeo(1), centerGeo(2), rM, az, wgs84);
    dNear = min(distance(lat, lon, bs_geo(:,1), bs_geo(:,2), wgs84));
    if dNear >= rMinM && dNear <= rMaxM
        return;
    end
end
end

function offsets = hexSiteOffsetsM(numRings, isdM)
% Θέσεις σταθμών εξαγωνικού πλέγματος [ανατολικά βόρεια] σε m γύρω από το
% κέντρο: 1 + 3*R*(R+1) σταθμοί για R δακτυλίους (R = 1 -> 7, R = 2 -> 19).
% Αξονικές συντεταγμένες (q, r) με διανύσματα βάσης απόστασης ίσης με isdM.
offsets = zeros(0,2);
for q = -numRings:numRings
    for r = max(-numRings, -q-numRings):min(numRings, -q+numRings)
        offsets(end+1,:) = isdM * [q + r/2, r*sqrt(3)/2]; %#ok<AGROW>
    end
end
% Ο κεντρικός σταθμός πρώτος (BS1), οι υπόλοιποι κατά απόσταση και αζιμούθιο.
distRel = round(hypot(offsets(:,1), offsets(:,2)) / isdM * 1e6) / 1e6;
azDeg   = mod(round(atan2d(offsets(:,1), offsets(:,2)) * 1e6) / 1e6, 360);
[~, order] = sortrows([distRel azDeg]);
offsets = offsets(order,:);
end

function spanRad = findMaxRaanOffset(raanBase, uPeakRad, meanMotion, inclDeg, aM, ...
    omegaEarth, tPeakS, windowS, centerLat, centerLon, wgs84, minElevDeg)
% Βρίσκει αριθμητικά τη μέγιστη εγκάρσια μετατόπιση του ίχνους για την οποία
% η διέλευση φτάνει ακόμη τη μάσκα ορατότητας. Διχοτόμηση στο [0, hi].
peak = @(dRaan) maxElevationOfPass(raanBase + dRaan, uPeakRad, meanMotion, ...
    inclDeg, aM, omegaEarth, tPeakS, windowS, centerLat, centerLon, wgs84);

lo = 0;
hi = deg2rad(25);
if peak(hi) > minElevDeg
    spanRad = hi;
    return;
end
for it = 1:30
    mid = 0.5*(lo+hi);
    if peak(mid) > minElevDeg
        lo = mid;
    else
        hi = mid;
    end
end
spanRad = lo;
end

function e = maxElevationOfPass(raanRad, uPeakRad, meanMotion, inclDeg, aM, ...
    omegaEarth, tPeakS, windowS, centerLat, centerLon, wgs84)
% Μέγιστη γωνία ανύψωσης του δορυφόρου, όπως τη βλέπει το κέντρο της περιοχής,
% μέσα στο παράθυρο της διέλευσης.
ts = linspace(tPeakS - windowS/2, tPeakS + windowS/2, 121);
e = -Inf;
for i = 1:numel(ts)
    lla = orbitPositionLla(uPeakRad + meanMotion*(ts(i)-tPeakS), inclDeg, raanRad, ...
        aM, omegaEarth, ts(i));
    [~, elev, ~] = geodetic2aer(lla(1), lla(2), lla(3), centerLat, centerLon, 0, wgs84);
    e = max(e, elev);
end
end

function halfS = zenithHalfVisibilityS(raanRad, uPeakRad, meanMotion, inclDeg, aM, ...
    omegaEarth, tPeakS, windowS, centerLat, centerLon, wgs84, minElevDeg)
% Μισή διάρκεια ορατότητας (πάνω από τη μάσκα) για τη διέλευση που περνά από
% το κέντρο της περιοχής, μετρημένη από τη μέγιστη προσέγγιση. Διχοτόμηση
% προς κάθε κατεύθυνση· κρατείται η μεγαλύτερη από τις δύο, που διαφέρουν
% ελάχιστα λόγω της περιστροφής της Γης.
elevAt = @(dt) elevationFromCenter(raanRad, uPeakRad + meanMotion*dt, inclDeg, aM, ...
    omegaEarth, tPeakS + dt, centerLat, centerLon, wgs84);
halfS = 0;
for sgn = [-1 1]
    lo = 0; hi = windowS;
    for it = 1:40
        mid = 0.5*(lo+hi);
        if elevAt(sgn*mid) >= minElevDeg
            lo = mid;
        else
            hi = mid;
        end
    end
    halfS = max(halfS, lo);
end
end

function e = elevationFromCenter(raanRad, uRad, inclDeg, aM, omegaEarth, tS, ...
    centerLat, centerLon, wgs84)
lla = orbitPositionLla(uRad, inclDeg, raanRad, aM, omegaEarth, tS);
[~, e, ~] = geodetic2aer(lla(1), lla(2), lla(3), centerLat, centerLon, 0, wgs84);
end

function lla = orbitPositionLla(uRad, inclDeg, raanRad, aM, omegaEarth, tS)
% Θέση δορυφόρου σε κυκλική τροχιά -> γεωδαιτικές συντεταγμένες [lat lon alt].
% Διάδοση σε αδρανειακό σύστημα και στροφή σε γεωκεντρικό-σταθερό κατά
% omegaEarth*t, ώστε θέση και ταχύτητα να προκύπτουν από το ίδιο μοντέλο.
rPerifocal = aM * [cos(uRad); sin(uRad); 0];

i = deg2rad(inclDeg);
Rx = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
Rz = [cos(raanRad) -sin(raanRad) 0; sin(raanRad) cos(raanRad) 0; 0 0 1];
rInertial = Rz * Rx * rPerifocal;

th = omegaEarth * tS;
Rg = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
rEcef = Rg * rInertial;

lla = ecef2lla(rEcef');
end
