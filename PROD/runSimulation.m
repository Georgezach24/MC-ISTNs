function [convTable, summary] = runSimulation(opts)
%RUNSIMULATION Ενιαία προσομοίωση του σεναρίου ISTN: μία εκτέλεση που παράγει
% όλα τα δεδομένα της εργασίας (αποτελέσματα προσομοίωσης + σύνολο εκπαίδευσης).
%
% Δομή του χρόνου, με την πηγή της κάθε επιλογής:
%
%   * Στοιχειώδης μονάδα = μία διέλευση του δορυφόρου. Το TR 38.821
%     Πίν. 4.2-3 NOTE 1 ορίζει ως φυσικό παράθυρο λειτουργίας τον χρόνο
%     ορατότητας του δορυφόρου ("a period of time corresponding to the
%     visibility time of the satellite").
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
%     M.2412 στις δορυφορικές αξιολογήσεις. Η προσομοίωση τρέχει επομένως
%     διέλευση-διέλευση μέχρι τα διαστήματα εμπιστοσύνης των δεικτών να
%     σταθεροποιηθούν, και η καμπύλη σύγκλισης αποθηκεύεται ως αποτέλεσμα.
%
%   * Κίνηση χρηστών: ITU-R M.2412-0 §8.4, ΠΙΝΑΚΑΣ 5 b) και c): "Fixed and
%     identical speed |v| of all UEs of the same mobility class, randomly and
%     uniformly distributed direction", με 3 km/h για πεζό χρήστη. Κάθε
%     διέλευση ξεκινά με νέες θέσεις χρηστών (ένα "drop" κατά M.2412 §7.1)
%     και οι χρήστες περπατούν συνεχώς μέσα στη διέλευση.
%
% Χρήση:
%   runSimulation()                                   % προεπιλογές
%   runSimulation(struct('maxPasses',10))             % σύντομο τρέξιμο δοκιμής
%   runSimulation(struct('ciTolerance',0.005))        % αυστηρότερη σύγκλιση
%
% Η προσομοίωση προχωρά πάντα με βήμα dtSeconds και όλοι οι δείκτες
% υπολογίζονται σε αυτή την ανάλυση. Στο αρχείο όμως γράφεται μία γραμμή κάθε
% datasetStride βήματα: διαδοχικά δείγματα του ίδιου χρήστη απέχουν 0.83 m και
% είναι σχεδόν ταυτόσημα, οπότε η πλήρης ανάλυση θα παρήγαγε αρχείο εκατοντάδων
% MB χωρίς να προσθέτει πληροφορία. Η διέλευση αναφοράς αποθηκεύεται χωριστά σε
% πλήρη ανάλυση, όπου χρειάζεται για τις χρονοσειρές.
%
% Έξοδοι:
%   Dataset/dataset.csv        μία γραμμή ανά (διέλευση, δείγμα, χρήστη)
%   Results/reference_pass.csv η διέλευση αναφοράς σε πλήρη ανάλυση 1 s
%   Results/convergence.csv    καμπύλη σύγκλισης ανά διέλευση
%   Results/policy_comparison.csv  δείκτες ανά διέλευση για τις τρεις πολιτικές
%   Results/policy_paired.csv      διαφορές κατά ζεύγη, με 95% CI
%   Results/*.png              χρονοσειρές διέλευσης αναφοράς + καμπύλη σύγκλισης

%% ------------------ Επιλογές ------------------
if nargin < 1 || isempty(opts)
    opts = struct();
end
thisDir = fileparts(mfilename('fullpath'));
def = struct( ...
    'dtSeconds',    1, ...      % s, βήμα (M.2412 Παρ.1 §5.3.2 + TR 38.821 Πίν. 7.3.2.1.4-1)
    'ueSpeedKmh',   3, ...      % km/h, πεζός χρήστης (M.2412 ΠΙΝ. 5 b/c)
    'passWindowS',  900, ...    % s, παράθυρο ανά διέλευση γύρω από τη μέγιστη προσέγγιση
    'minPasses',    30, ...     % ελάχιστες διελεύσεις πριν ελεγχθεί η σύγκλιση
    'maxPasses',    600, ...    % ασφαλιστικό άνω όριο
    'ciTolerance',  0.05, ...   % σχετικό ημιεύρος 95% CI (δείκτες ρυθμού/ενέργειας)
    'ciToleranceFrac', 0.02, ...% απόλυτο ημιεύρος 95% CI (δείκτες ποσοστού)
    'datasetStride', 5, ...     % κάθε πόσα βήματα γράφεται γραμμή στο dataset
    'resume',       true, ...   % συνέχιση από σημείο ελέγχου, αν υπάρχει
    'rngSeed',      42, ...
    'datasetPath',  fullfile(thisDir, '..', 'Dataset', 'dataset.csv'), ...
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
refPassCsv     = fullfile(opts.outputDir, 'reference_pass.csv');
checkpointPath = fullfile(opts.outputDir, 'runSimulation_checkpoint.mat');

% Σημείο ελέγχου: μια μεγάλη εκτέλεση μπορεί να διακοπεί (π.χ. από έλλειψη
% μνήμης στο μηχάνημα). Αν υπάρχει συμβατό σημείο ελέγχου, η εκτέλεση
% συνεχίζει από την επόμενη διέλευση αντί να ξαναρχίσει από την αρχή.
resumeState = [];
if opts.resume && isfile(checkpointPath) && isfile(opts.datasetPath)
    S = load(checkpointPath);
    if isfield(S,'ckpt') && checkpointMatches(S.ckpt.opts, opts)
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
    if isfile(refPassCsv),        delete(refPassCsv);        end
    if isfile(checkpointPath),    delete(checkpointPath);    end
end

wgs84 = wgs84Ellipsoid;
cLight = physconst('LightSpeed');

%% ------------------ Τοπολογία σταθμών βάσης (σταθερή: το σενάριο) ------------------
% Ίδιο σημείο αναφοράς και ίδια διάταξη BS με το προηγούμενο σενάριο αναφοράς.
bs_geo = [37.9838 23.7275 25;
          37.9865 23.7310 25];
numBs = size(bs_geo,1);

% Χρήστες της διέλευσης αναφοράς (διέλευση 1): οι ίδιες θέσεις με το
% προηγούμενο στατικό σενάριο, ώστε να παραμένουν συγκρίσιμα τα αποτελέσματα.
refUserGeo = [37.9845 23.7288 1.5;
              37.9870 23.7325 1.5;
              37.9825 23.7268 1.5;
              37.9900 23.7450 1.5;
              37.9997 23.7476 1.5;
              37.9484 23.7111 1.5];

numUsersRange = [3 10];   % πλήθος χρηστών ανά διέλευση (πλην της αναφοράς)
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

% Δορυφόρος: TR 38.821 Set-1, LEO-600, S-band (Πίν. 6.1.1.1-1 & 6.1.3.2-1).
% Η πυκνότητα EIRP (dBW/MHz) είναι το δεδομένο· EIRP και TxPower παράγωγα.
satAltitude = 600e3;
satParameters.CarrierFrequency = 2.0e9;
satParameters.Bandwidth = 30e6;
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
% δειγματοληπτούνται ώστε καμία τροχιά να μην μπορεί να βγει εκτός ορίων.
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

%% ------------------ Δείκτες σύγκλισης ------------------
% Παρακολουθούνται ανά διέλευση (η διέλευση είναι η ανεξάρτητη μονάδα: οι
% χρήστες μέσα στην ίδια διέλευση αλληλεξαρτώνται μέσω της κατανομής πόρων).
kpiNames = {'MeanCapacity_Mbps','MeanSinr_dB','FracServed','FracSatellite','MeanBitPerJouleRf'};
kpiIsFraction = [false false true true false];
% Ο κανόνας τερματισμού δεν περιλαμβάνει το μέσο SINR: είναι μέσος όρος
% λογαριθμικού μεγέθους πάνω σε ανάμεικτο πληθυσμό επίγειων και δορυφορικών
% χρηστών, με εύρος δεκάδων dB, και η εργασία δεν τον αναφέρει ως δείκτη (το
% SINR παρουσιάζεται ως κατανομή). Το διάστημα εμπιστοσύνης του καταγράφεται
% κανονικά στην καμπύλη σύγκλισης.
kpiInStopRule = [true false true true true];
numKpi = numel(kpiNames);
passKpi = nan(opts.maxPasses, numKpi);

%% ------------------ Σύγκριση πολιτικών (ITU-R M.2412-0 §7.1) ------------------
% Κάθε βήμα αποτιμάται και με τις τρεις πολιτικές πάνω στην ΙΔΙΑ πραγματοποίηση
% καναλιού: τα δύο υποψήφια SINR είναι κοινά, αλλάζει μόνο ποιος επιλέγεται και
% άρα ο φόρτος κάθε κόμβου. Η σύγκριση είναι επομένως κατά ζεύγη και όχι μεταξύ
% ανεξάρτητων εκτελέσεων, οπότε η διακύμανση της τοπολογίας και του καναλιού
% απαλείφεται από τη διαφορά.
policyNames   = {'Actual','TerrestrialOnly','SatelliteOnly'};
policyMetrics = {'MeanCapacity_Mbps','TotalRate_Mbps','FracServed', ...
                 'FracBelowTarget','FracOutage','BitPerJouleRf'};
numPolicy       = numel(policyNames);
numPolicyMetric = numel(policyMetrics);
policyPassKpi   = nan(opts.maxPasses, numPolicy, numPolicyMetric);

convRows = cell(opts.maxPasses,1);
checksRun = 0;
totalRows = 0;
converged = false;
convergedAtPass = NaN;

fprintf('Ενιαία προσομοίωση ISTN\n');
fprintf('  Δt = %g s | ταχύτητα χρήστη = %g km/h (%.3f m/s) | παράθυρο διέλευσης = %g s\n', ...
    opts.dtSeconds, opts.ueSpeedKmh, ueSpeedMps, opts.passWindowS);
fprintf('  Μετατόπιση ανά βήμα = %.3f m (όριο M.2412 Παρ.1 §5.3.2: 1 m)\n', ueSpeedMps*opts.dtSeconds);
fprintf('  Τροχιακή περίοδος = %.0f s | ταχύτητα ίχνους = %.0f m/s\n', orbitalPeriodS, groundSpeedMps);
fprintf('  Εύρος εγκάρσιας μετατόπισης ίχνους = ±%.2f°\n', rad2deg(raanSpanRad));
fprintf('  Σύγκλιση (ITU-R M.2412-0 §7.1): 95%% CI < %.0f%% σχετικό / %.0f ποσοστιαίες μονάδες\n', ...
    100*opts.ciTolerance, 100*opts.ciToleranceFrac);
fprintf('  Δείκτες κανόνα τερματισμού: %s | ελάχιστο %d διελεύσεις\n\n', ...
    strjoin(kpiNames(kpiInStopRule), ', '), opts.minPasses);

%% ------------------ Βρόχος διελεύσεων ------------------
passIdx = 0;
refPassTable = table();
refSnapshot = struct();

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
    if isfile(refPassCsv)
        refPassTable = readtable(refPassCsv, 'TextType', 'string');
    end
end

while passIdx < opts.maxPasses
    passIdx = passIdx + 1;

    % Ένας σπόρος ανά διέλευση: η τυχαιότητα μέσα στη διέλευση (τοπολογία,
    % κατευθύνσεις, κανάλι) εξελίσσεται συνεχόμενα, ώστε να μη σπάει η
    % χωρική συσχέτιση από επαναρχικοποίηση της γεννήτριας ανά βήμα.
    rng(opts.rngSeed + passIdx);

    % -- Θέσεις χρηστών στην αρχή της διέλευσης --
    if passIdx == 1
        user_geo = refUserGeo;
        user_geo = clampUsersIntoRange(user_geo, bs_geo, wgs84, userRMinM, userRMaxM);
    else
        nU = randi(numUsersRange);
        user_geo = zeros(nU,3);
        for u = 1:nU
            refB = randi(numBs);
            [dLat, dLon] = randAnnulusOffsetDeg(bs_geo(refB,1), userRMinM, userRMaxM);
            user_geo(u,:) = [bs_geo(refB,1)+dLat, bs_geo(refB,2)+dLon, ueHeightM];
        end
    end
    numUsers = size(user_geo,1);

    % -- Κατεύθυνση κίνησης: σταθερή ταχύτητα, τυχαία ομοιόμορφη κατεύθυνση --
    % ITU-R M.2412-0 §8.4, ΠΙΝΑΚΑΣ 5 b)/c).
    walkAzimuthDeg = 360*rand(numUsers,1);

    % -- Γεωμετρία της συγκεκριμένης διέλευσης --
    raanRad = raanBase + (2*rand()-1)*raanSpanRad;
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
    stepBitPerJouleRf = nan(numSteps,1);
    policyStepBuf = nan(numSteps, numPolicy, numPolicyMetric);

    for s = 1:numSteps
        t = tStart + (s-1)*opts.dtSeconds;
        sat_geo = orbitPositionLla(uPeakRad + meanMotion*(t - tPeakS), inclDeg, ...
            raanRad, a, omegaEarth, t);

        % --- Έλεγχος κατάστασης χρηστών πριν τον υπολογισμό ζεύξης ---
        checksRun = checksRun + 1;
        validateUserState(user_geo, bs_geo, wgs84, ueHeightM, d2dMinM, d2dMaxM, passIdx, s);

        [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
            bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
            nodePowerWattsVec, energyPerBitUJVec, ...
            bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
            satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
            channelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
            bsReasonVec, satReasonVec, policyKpis] = ...
            simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, ...
                             satParameters, channelState);

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

    % -- Εγγραφή: κρατούνται τα βήματα 1, 1+stride, ... --
    keep = mod(numBuf(:,2) - 1, opts.datasetStride) == 0;
    passTable = buildPassTable(numBuf(keep,:), strBuf(keep,:));
    if passIdx == 1
        writetable(passTable, opts.datasetPath);
        % Η διέλευση αναφοράς αποθηκεύεται και σε πλήρη ανάλυση, για τις
        % χρονοσειρές: εκεί χρειάζεται κάθε βήμα, όχι δείγμα.
        refPassTable = buildPassTable(numBuf, strBuf);
        writetable(refPassTable, refPassCsv);
    else
        writetable(passTable, opts.datasetPath, 'WriteMode', 'append');
    end
    totalRows = totalRows + height(passTable);
    clear numBuf strBuf passTable;

    % -- Σύγκλιση: 95% CI του μέσου όρου πάνω στις διελεύσεις --
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

    tolVec = opts.ciTolerance * ones(1,numKpi);
    tolVec(kpiIsFraction) = opts.ciToleranceFrac;
    if passIdx >= opts.minPasses && all(relHw(kpiInStopRule) < tolVec(kpiInStopRule))
        converged = true;
        convergedAtPass = passIdx;
    end

    if mod(passIdx,10) == 0 || passIdx <= 3 || converged
        fprintf('  Διέλευση %3d | γραμμές %8d | χωρητ. %6.2f Mbps | served %5.1f%% | sat %5.1f%% | max rel.CI %.4f | μνήμη %5.2f GB\n', ...
            passIdx, totalRows, mean(passKpi(1:passIdx,1),'omitnan'), ...
            100*mean(passKpi(1:passIdx,3),'omitnan'), ...
            100*mean(passKpi(1:passIdx,4),'omitnan'), max(relHw(kpiInStopRule)), ...
            matlabMemoryGb());
    end

    % Σημείο ελέγχου μετά από κάθε διέλευση: λίγα KB, ώστε μια διακοπή να
    % κοστίζει το πολύ μία διέλευση αντί για ολόκληρη την εκτέλεση.
    ckpt = struct('passIdx', passIdx, 'passKpi', passKpi(1:passIdx,:), ...
        'policyPassKpi', policyPassKpi(1:passIdx,:,:), ...
        'convRows', {convRows(1:passIdx)}, 'totalRows', totalRows, ...
        'checksRun', checksRun, 'refSnapshot', refSnapshot, 'opts', opts);
    save(checkpointPath, 'ckpt');

    if converged
        break;
    end
end

numPasses = passIdx;

%% ------------------ Σύγκριση πολιτικών κατά ζεύγη ------------------
% Κάθε διέλευση δίνει τρεις τιμές για τον ίδιο δείκτη, μία ανά πολιτική, πάνω
% στο ίδιο κανάλι και την ίδια τοπολογία. Η στατιστική γίνεται στη ΔΙΑΦΟΡΑ ανά
% διέλευση, όχι στους δύο μέσους όρους χωριστά.
pol = policyPassKpi(1:numPasses,:,:);

polPassID = repelem((1:numPasses)', numPolicy);
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
policyCsv = fullfile(opts.outputDir, 'policy_comparison.csv');
writetable(policyTable, policyCsv);

% Οι δύο συγκρίσεις: η υπό εξέταση πολιτική έναντι καθεμιάς από τις δύο
% αποκλειστικές.
compIdx = [1 2; 1 3];
nComp   = size(compIdx,1);
pairBaseline = strings(nComp*numPolicyMetric,1);
pairMetric   = strings(nComp*numPolicyMetric,1);
pairN        = zeros(nComp*numPolicyMetric,1);
pairMean     = nan(nComp*numPolicyMetric,1);
pairHw       = nan(nComp*numPolicyMetric,1);
pairRelPct   = nan(nComp*numPolicyMetric,1);
pairWinFrac  = nan(nComp*numPolicyMetric,1);
q = 0;
for c = 1:nComp
    for mI = 1:numPolicyMetric
        a = pol(:,compIdx(c,1),mI);
        b = pol(:,compIdx(c,2),mI);
        d = a - b;
        ok = isfinite(d);
        d = d(ok);
        q = q + 1;
        pairBaseline(q) = string(policyNames{compIdx(c,2)});
        pairMetric(q)   = string(policyMetrics{mI});
        pairN(q)        = numel(d);
        if numel(d) >= 2
            pairMean(q)    = mean(d);
            pairHw(q)      = 1.96*std(d)/sqrt(numel(d));
            baseMean       = mean(b(ok), 'omitnan');
            if abs(baseMean) > eps
                pairRelPct(q) = 100*mean(d)/abs(baseMean);
            end
            pairWinFrac(q) = mean(d > 0);
        end
    end
end
pairedTable = table(pairBaseline, pairMetric, pairN, pairMean, pairHw, ...
    pairMean - pairHw, pairMean + pairHw, pairRelPct, pairWinFrac, ...
    'VariableNames', {'Baseline','Metric','NumPasses','MeanDiff','CI95HalfWidth', ...
                      'CI95Low','CI95High','RelDiffPct','FracPassesPositive'});
pairedCsv = fullfile(opts.outputDir, 'policy_paired.csv');
writetable(pairedTable, pairedCsv);

fprintf('\n  Σύγκριση πολιτικών στην ίδια διέλευση (%d διελεύσεις)\n', numPasses);
fprintf('    %-18s %12s %12s %12s\n', 'πολιτική', 'ρυθμ.[Mbps]', 'εξυπηρ.', 'εκτός');
for pI = 1:numPolicy
    fprintf('    %-18s %12.2f %12.4f %12.4f\n', policyNames{pI}, ...
        mean(pol(:,pI,2),'omitnan'), mean(pol(:,pI,3),'omitnan'), ...
        mean(pol(:,pI,5),'omitnan'));
end
fprintf('\n    Διαφορές ανά διέλευση (η υπό εξέταση μείον τη βάση), 95%% CI:\n');
for q = 1:height(pairedTable)
    signif = ternary(pairedTable.CI95Low(q) > 0 || pairedTable.CI95High(q) < 0, ...
        '', '   [το CI περιέχει το μηδέν]');
    fprintf('    %-16s %-18s %+10.4f ± %.4f (%+6.1f%%, θετική στο %4.1f%% των διελεύσεων)%s\n', ...
        pairedTable.Baseline(q), pairedTable.Metric(q), pairedTable.MeanDiff(q), ...
        pairedTable.CI95HalfWidth(q), pairedTable.RelDiffPct(q), ...
        100*pairedTable.FracPassesPositive(q), signif);
end
fprintf('\n');

%% ------------------ Γράφημα σύγκρισης πολιτικών ------------------
figPol = figure('Visible','off','Position',[100 100 1050 330]);

subplot(1,3,1);
mv = arrayfun(@(pI) mean(pol(:,pI,2),'omitnan'), 1:numPolicy);
ev = arrayfun(@(pI) 1.96*std(pol(:,pI,2),'omitnan')/sqrt(numPasses), 1:numPolicy);
bar(mv); hold on; errorbar(1:numPolicy, mv, ev, 'k', 'LineStyle','none', 'LineWidth',1.2);
set(gca,'XTickLabel',{'υπό εξέταση','μόνο επίγεια','μόνο δορυφ.'});
ylabel('Συνολική ρυθμαπόδοση [Mbps]'); grid on;
title('Ρυθμαπόδοση ανά πολιτική');

subplot(1,3,2);
mv2 = arrayfun(@(pI) 100*mean(pol(:,pI,3),'omitnan'), 1:numPolicy);
ev2 = arrayfun(@(pI) 100*1.96*std(pol(:,pI,3),'omitnan')/sqrt(numPasses), 1:numPolicy);
bar(mv2); hold on; errorbar(1:numPolicy, mv2, ev2, 'k', 'LineStyle','none', 'LineWidth',1.2);
set(gca,'XTickLabel',{'υπό εξέταση','μόνο επίγεια','μόνο δορυφ.'});
ylabel('Εξυπηρετούμενοι [%]'); grid on;
title('Εξυπηρέτηση ανά πολιτική');

subplot(1,3,3);
dT = pol(:,1,2) - pol(:,2,2);
dS = pol(:,1,2) - pol(:,3,2);
histogram(dT, 'DisplayName','έναντι μόνο επίγειας'); hold on;
histogram(dS, 'DisplayName','έναντι μόνο δορυφορικής');
xline(0,'k--','HandleVisibility','off');
xlabel('Διαφορά ρυθμαπόδοσης ανά διέλευση [Mbps]'); ylabel('Διελεύσεις');
legend('Location','northoutside'); grid on;
title('Κατανομή των διαφορών');

polPng = fullfile(opts.outputDir, 'policy_comparison.png');
saveas(figPol, polPng);
close(figPol);


%% ------------------ Καμπύλη σύγκλισης ------------------
convMat = vertcat(convRows{1:numPasses});
convVarNames = [{'PassID','CumulativeRows'}, ...
    strcat('Mean_', kpiNames), strcat('RelCI_', kpiNames)];
convTable = array2table(convMat, 'VariableNames', convVarNames);
convCsv = fullfile(opts.outputDir, 'convergence.csv');
writetable(convTable, convCsv);

%% ------------------ Σύνοψη ------------------
simulatedTimeS = numPasses * opts.passWindowS;
summary = struct();
summary.numPasses       = numPasses;
summary.numRows         = totalRows;
summary.simulatedTime_s = simulatedTimeS;
summary.converged       = converged;
summary.convergedAtPass = convergedAtPass;
summary.checksRun       = checksRun;
for k = 1:numKpi
    summary.(kpiNames{k}) = mean(passKpi(1:numPasses,k), 'omitnan');
end

fprintf('\n--- Σύνοψη ---\n');
lastHw = convMat(end, end-numKpi+1:end);
if converged
    fprintf('  Σύγκλιση στη διέλευση %d (95%% CI εντός ορίου για όλους τους δείκτες του κανόνα)\n', ...
        convergedAtPass);
else
    fprintf('  ΔΕΝ συνέκλινε μέσα σε %d διελεύσεις. Μέγιστο CI κανόνα = %.4f\n', ...
        numPasses, max(lastHw(kpiInStopRule)));
end
for k = 1:numKpi
    if kpiIsFraction(k)
        fprintf('    %-20s %10.4f  (95%% CI ±%.4f)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), lastHw(k), ternary(kpiInStopRule(k),'',' [εκτός κανόνα]'));
    else
        fprintf('    %-20s %10.4f  (95%% CI ±%.2f%%)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), 100*lastHw(k), ternary(kpiInStopRule(k),'',' [εκτός κανόνα]'));
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
        size(refSnapshot.user_geo,1), refSnapshot.bestNodeTypeVec, refSnapshot.bestNodeVec);
    saveas(gcf, fullfile(opts.outputDir, 'network_3d.png'));
    close(gcf);
end

% Το dataset ΔΕΝ αντιγράφεται στον versioned φάκελο: είναι δεκάδες MB και το
% Results/runs/ παρακολουθείται από το git. Αντ' αυτού καταγράφεται το άθροισμα
% ελέγχου SHA-256 του, που το συνδέει με τη συγκεκριμένη εκτέλεση και το commit.
outFiles = {refPassCsv, convCsv, policyCsv, pairedCsv, polPng, ...
            fullfile(opts.outputDir,'network_3d.png')};

kpiList  = {'Capacity_Mbps','EnergyPerBit_uJ','SINR_dB','SatElevation_deg'};
kpiLabel = {'Χωρητικότητα (Mbps)','Ενέργεια ανά bit (\muJ/bit)','SINR (dB)','Γωνία ανύψωσης (deg)'};
for k = 1:numel(kpiList)
    fig = figure('Visible','off');
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

% Καμπύλη σύγκλισης
fig = figure('Visible','off');
plot(convTable.PassID, convTable.(['RelCI_' kpiNames{1}]), '-o', 'DisplayName', kpiNames{1});
hold on;
for k = 2:numKpi
    plot(convTable.PassID, convTable.(['RelCI_' kpiNames{k}]), '-o', 'DisplayName', kpiNames{k});
end
yline(opts.ciTolerance, '--k', 'DisplayName', 'κατώφλι (σχετικό)');
yline(opts.ciToleranceFrac, ':k', 'DisplayName', 'κατώφλι (ποσοστά)');
hold off; grid on;
xlabel('Πλήθος διελεύσεων'); ylabel('Ημιεύρος 95% CI (σχετικό)');
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
runParams.datasetStride    = opts.datasetStride;
runParams.datasetSampling_s = opts.datasetStride*opts.dtSeconds;
runParams.converged        = converged;
runParams.convergedAtPass  = convergedAtPass;
runParams.ciTolerance      = opts.ciTolerance;
runParams.ciToleranceFrac  = opts.ciToleranceFrac;
runParams.stopRuleKpis     = strjoin(kpiNames(kpiInStopRule), ', ');
runParams.policies         = strjoin(policyNames, ', ');
runParams.policyMetrics    = strjoin(policyMetrics, ', ');
runParams.minPasses        = opts.minPasses;
runParams.maxPasses        = opts.maxPasses;
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
runParams.refUserGeo       = refUserGeo;
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
    'duration',   'ITU-R M.2412-0 §7.1 (σύγκλιση + διαστήματα εμπιστοσύνης)· ITU-R M.2514-0 §8.2.4', ...
    'passWindow', 'TR 38.821 Πίν. 4.2-3 NOTE 1 (visibility time of the satellite)', ...
    'ueMobility', 'ITU-R M.2412-0 §8.4 ΠΙΝΑΚΑΣ 5 b)/c) (3 km/h, τυχαία ομοιόμορφη κατεύθυνση)', ...
    'validity',   'TR 38.901 §7.4.1 Πίν. 7.4.1-1 (10 m - 5 km, h_UT 1.5-22.5 m)');

saveRunVersion('runSimulation', runParams, outFiles, opts.label);

% Η εκτέλεση ολοκληρώθηκε: το σημείο ελέγχου δεν χρειάζεται πια.
if isfile(checkpointPath)
    delete(checkpointPath);
end

end

% =====================================================================
% Τοπικές συναρτήσεις
% =====================================================================

function tf = checkpointMatches(a, b)
% Το σημείο ελέγχου χρησιμοποιείται μόνο αν οι παράμετροι που καθορίζουν τι
% παράγει η εκτέλεση είναι ίδιες. Οτιδήποτε άλλο (π.χ. label) δεν πειράζει.
keys = {'dtSeconds','ueSpeedKmh','passWindowS','datasetStride','rngSeed', ...
        'ciTolerance','ciToleranceFrac','minPasses','maxPasses','datasetPath'};
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
    h = '';   % π.χ. MATLAB χωρίς JVM
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

function user_geo = clampUsersIntoRange(user_geo, bs_geo, wgs84, rMinM, rMaxM)
% Φέρνει τους χρήστες της διέλευσης αναφοράς μέσα στο επιτρεπτό δακτύλιο
% γύρω από τον πλησιέστερο σταθμό βάσης, διατηρώντας τον αζιμούθιό τους, ώστε
% καμία τροχιά να μη βγαίνει εκτός πεδίου ισχύος κατά τη διάρκεια της διέλευσης.
for u = 1:size(user_geo,1)
    d2d = distance(user_geo(u,1), user_geo(u,2), bs_geo(:,1), bs_geo(:,2), wgs84);
    [dMin, b] = min(d2d);
    if dMin >= rMinM && dMin <= rMaxM
        continue;
    end
    az = azimuth(bs_geo(b,1), bs_geo(b,2), user_geo(u,1), user_geo(u,2), wgs84);
    rNew = min(max(dMin, rMinM), rMaxM);
    [lat, lon] = reckon(bs_geo(b,1), bs_geo(b,2), rNew, az, wgs84);
    user_geo(u,1) = lat;
    user_geo(u,2) = lon;
end
end

function [dLat, dLon] = randAnnulusOffsetDeg(refLat, rMinM, rMaxM)
% Τυχαία μετατόπιση σε μοίρες, ομοιόμορφη ως προς το εμβαδόν δακτυλίου:
% r = sqrt(rmin^2 + U*(rmax^2 - rmin^2)), theta = 2*pi*V.
rM = sqrt(rMinM^2 + rand()*(rMaxM^2 - rMinM^2));
th = 2*pi*rand();
mPerDegLat = 111132.92;
dLat = (rM*cos(th)) / mPerDegLat;
dLon = (rM*sin(th)) / (mPerDegLat*cosd(refLat));
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
