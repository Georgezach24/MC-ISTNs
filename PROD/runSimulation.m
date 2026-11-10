function [convTable, summary] = runSimulation(opts)
%RUNSIMULATION Προσομοίωση του σεναρίου ISTN (επίγειο UMa + LEO).
%
%   * Μονάδα = μία διέλευση των δορυφόρων (TR 38.821 Πίν. 4.2-3 NOTE 1)· κάθε
%     διέλευση είναι ένα drop με νέους χρήστες (ITU-R M.2412-0 §7.1).
%   * Βήμα dtSeconds = 0.2 s = περίοδος μέτρησης L1->L3.
%   * Σταθερό πλήθος διελεύσεων (numPasses), ίδιο σε όλες τις πολιτικές, ώστε
%     οι συγκρίσεις να γίνονται κατά ζεύγη· τα 95% CI αναφέρονται στο τέλος
%     (M.2412-0 §7.1, M.2514-0 §8.2.4).
%   * Χρήστες: 3 km/h, τυχαία σταθερή κατεύθυνση (M.2412-0 §8.4, Πίν. 5).
%
% Χρήση:
%   runSimulation()                          % προεπιλογές
%   runSimulation(struct('numPasses',10))    % σύντομη δοκιμή
%
% Έξοδοι:
%   Dataset/dataset.csv            μία γραμμή ανά (διέλευση, δείγμα, χρήστη),
%                                  κάθε datasetStride βήματα
%   Dataset/handover_events.mat    ένα γεγονός ανά handover (εκπαίδευση ML)
%   Results/*.png                  διέλευση αναφοράς, σύγκλιση, κέρδη πολιτικών
%   Results/simulation_results.mat οι αριθμοί πίσω από τα γραφήματα

%% ------------------ Επιλογές ------------------
if nargin < 1 || isempty(opts)
    opts = struct();
end
thisDir = fileparts(mfilename('fullpath'));
def = struct( ...
    'dtSeconds',    0.2, ...    % s, βήμα = περίοδος L1->L3 (TR 36.839 Πίν. 5.3.2.1, TS 38.133 Πίν. 9.2.5.2-1)
    'ueSpeedKmh',   3, ...      % km/h, πεζός χρήστης (M.2412 ΠΙΝ. 5 b/c)
    'passWindowS',  900, ...    % s, παράθυρο ανά διέλευση γύρω από τη μέγιστη προσέγγιση
    'numPasses',    300, ...    % σταθερό πλήθος διελεύσεων
    'geometryGrid', true, ...   % στρωματοποιημένη εγκάρσια μετατόπιση ίχνους· false = rand() ανά διέλευση
    'ciTolerance',  0.05, ...   % στόχος: σχετικό ημιεύρος 95% CI (ρυθμός/ενέργεια)
    'ciToleranceFrac', 0.02, ...% στόχος: απόλυτο ημιεύρος 95% CI (ποσοστά)
    'datasetStride', 25, ...    % κάθε πόσα βήματα γράφεται γραμμή στο dataset (25 x 0.2 s = 5 s)
    'hysteresisDb',  3, ...      % TS 38.331 §6.3.2: 0..15 dB, βήμα 0.5 dB (τιμή πεδίου 6)
    'timeToTriggerMs', 480, ...  % TR 36.839 Πίν. 5.3.2.1 Set 1 (με offset 3 dB, K = 4)· 3 βήματα στο Δt = 0.2 s
    'pingPongMtsS',  1, ...      % s, ελάχιστος χρόνος παραμονής (MTS) για τον ορισμό του ping-pong
    'rlfQoutDb',    -8, ...      % dB, κατώφλι out-of-sync (TR 36.839 Πίν. 5.2.1.3.1)
    'rlfQinDb',     -6, ...      % dB, κατώφλι in-sync (TR 36.839 Πίν. 5.2.1.3.1)
    'rlfT310Ms',  1000, ...      % ms, T310 (TR 36.839 Πίν. 5.2.1.3.1, προεπιλογή του 36.331)
    'rlfN310',       1, ...      % διαδοχικά out-of-sync που ξεκινούν τον T310 (TR 36.839)
    'rlfN311',       1, ...      % διαδοχικά in-sync που τον σταματούν (TR 36.839)
    'softLosState',  true, ...   % soft LOS (TR 38.901 §7.6.3.3, εξ. 7.6-18/19)· false = σκληρή μετάβαση
    'l3FilterK',     4, ...      % TS 38.331 §6.3.2 filterCoefficient (DEFAULT fc4)· 0 = χωρίς φίλτρο
    'measPeriodMs',  200, ...    % TS 38.133 Πίν. 9.2.5.2-1, FR1 χωρίς DRX
    'measAveraging', true, ...   % μέση τιμή μέτρησης μέσα στο βήμα· false = στιγμιαίο δείγμα
    'numUsersRange', [10 30], ...% [min max] χρήστες ανά διέλευση, ομοιόμορφα
    'bsRings',       2, ...      % δακτύλιοι εξαγωνικού πλέγματος: 2 -> 19 σταθμοί (TR 38.901 Πίν. 7.2-1)
    'bsIsdM',        500, ...    % m, απόσταση σταθμών (TR 38.901 Πίν. 7.8-1, UMa)
    'numSats',       3, ...      % διαδοχικοί δορυφόροι στο ίδιο τροχιακό επίπεδο
    'satFreqReuse',  1, ...      % TR 38.821 Πίν. 6.1.1.1-5: Option 1 = 1, Option 2 = 3
    'mlbReportPeriodMs', 1000, ...% ms, Reporting Periodicity της αναφοράς φόρτου (TS 38.423 §9.1.3.18)
    'mlbMaxOffsetDb', 6, ...     % dB ανά μονάδα διαφοράς φόρτου (Nasri & Altman: HMmax - f(0) = 12 - 6)
    'mlGatePath', fullfile(thisDir, '..', 'Model', 'results_handover_gate', 'handover_gate.mat'), ... % μοντέλο ML· χωρίς αρχείο = χωρίς τις πολιτικές *ML
    'mlReportIntervalMs', 240, ...% ms, reportInterval (TS 38.331 ReportConfigNR)
    'wrapAround',    true, ...   % wrap-around του πλέγματος (TR 38.901)
    'bsSectors',     3, ...      % τομείς ανά site: 3 (TR 38.901 Πίν. 7.8-1) ή 1
    'numWorkers',    6, ...      % παράλληλες διελεύσεις (Parallel Computing Toolbox)· 0 = σειριακά
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

% Σημείο ελέγχου: μια διακοπείσα εκτέλεση συνεχίζει από την επόμενη διέλευση,
% αν οι παράμετροι και το μοντέλο ML είναι τα ίδια.
resumeState = [];
if opts.resume && isfile(checkpointPath) && isfile(opts.datasetPath)
    S = load(checkpointPath);
    gateSha = '';
    if ~isempty(opts.mlGatePath) && isfile(opts.mlGatePath)
        gateSha = fileSha256(opts.mlGatePath);
    end
    if isfield(S,'ckpt') && isfield(S.ckpt,'mlGateSha256') && ...
            strcmp(S.ckpt.mlGateSha256, gateSha) && checkpointMatches(S.ckpt.opts, opts)
        % Διακοπή ανάμεσα στην εγγραφή του CSV και του σημείου ελέγχου θα
        % έδινε διπλές γραμμές: σταματάμε.
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

%% ------------------ Τοπολογία σταθμών βάσης ------------------
% Εξαγωνικό πλέγμα (TR 38.901 Πίν. 7.2-1, 7.8-1: UMa, 19 sites, ISD 500 m, 25 m).
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

% Θέσεις του επόμενου δακτυλίου (χωρίς σταθμό): ορίζουν το όριο της περιοχής·
% ένα σημείο ανήκει στο δίκτυο όταν ο πλησιέστερος κόμβος είναι σταθμός.
ringOffsets = hexSiteOffsetsM(opts.bsRings + 1, opts.bsIsdM);
ringOffsets = ringOffsets(numBs+1:end,:);
ghost_geo = zeros(size(ringOffsets,1), 2);
for b = 1:size(ringOffsets,1)
    [ghost_geo(b,1), ghost_geo(b,2)] = reckon(bsCenterGeo(1), bsCenterGeo(2), ...
        norm(ringOffsets(b,:)), atan2d(ringOffsets(b,1), ringOffsets(b,2)), wgs84);
end
cellRadiusM = opts.bsIsdM / sqrt(3);   % ακτίνα εξαγώνου (κορυφή)

% Wrap-around: 6 αντίγραφα του πλέγματος (αξονικές μετατοπίσεις ±(2R+1, -R),
% ±(R, R+1), ±(-R-1, 2R+1))· κάθε ζεύξη μετριέται από το πλησιέστερο αντίγραφο,
% οπότε όλες οι κυψέλες βλέπουν τον ίδιο περίγυρο παρεμβολέων.
R = opts.bsRings;
wrapQR = [0 0; 2*R+1 -R; R R+1; -R-1 2*R+1; -(2*R+1) R; -R -(R+1); R+1 -(2*R+1)];
wrapShiftM = opts.bsIsdM * [wrapQR(:,1) + wrapQR(:,2)/2, wrapQR(:,2)*sqrt(3)/2];
if ~opts.wrapAround
    wrapShiftM = wrapShiftM(1,:);
end
numImg = size(wrapShiftM,1);
bsImageGeo = zeros(numBs, numImg, 2);
for b = 1:numBs
    for k = 1:numImg
        o = hexOffsets(b,:) + wrapShiftM(k,:);
        [bsImageGeo(b,k,1), bsImageGeo(b,k,2)] = reckon(bsCenterGeo(1), bsCenterGeo(2), ...
            norm(o), atan2d(o(1), o(2)), wgs84);
    end
end
% Έλεγχος πλακόστρωσης: ο επόμενος δακτύλιος = αντίγραφα σταθμών, χωρίς
% επικαλύψεις, και ίδιος περίγυρος για κάθε σταθμό.
if opts.wrapAround
    allImg = reshape(reshape(hexOffsets, numBs, 1, 2) + reshape(wrapShiftM, 1, numImg, 2), [], 2);
    for g = 1:size(ringOffsets,1)
        if min(hypot(allImg(:,1) - ringOffsets(g,1), allImg(:,2) - ringOffsets(g,2))) > 1e-6
            error('runSimulation:WrapAround', 'Η θέση %d του επόμενου δακτυλίου δεν είναι αντίγραφο σταθμού.', g);
        end
    end
    dPair = hypot(allImg(:,1) - allImg(:,1)', allImg(:,2) - allImg(:,2)');
    dPair(1:size(dPair,1)+1:end) = Inf;
    if min(dPair(:)) < opts.bsIsdM - 1e-6
        error('runSimulation:WrapAround', 'Επικαλυπτόμενα αντίγραφα σταθμών.');
    end
    p0 = [123.4 56.7];
    dRef = [];
    for b = 1:numBs
        q = hexOffsets(b,:) + p0;
        dW = min(reshape(hypot(q(1) - allImg(:,1), q(2) - allImg(:,2)), numBs, numImg), [], 2);
        if isempty(dRef)
            dRef = sort(dW);
        elseif max(abs(sort(dW) - dRef)) > 1e-6
            error('runSimulation:WrapAround', 'Ο περίγυρος του σταθμού %d διαφέρει από του κεντρικού.', b);
        end
    end
end

numUsersRange = opts.numUsersRange;
ueHeightM     = 1.5;      % m, TR 38.901 Πίν. 7.4.1-1 (υπαίθριος χρήστης)

%% ------------------ Ραδιο-παράμετροι ------------------
% Επίγειο σκέλος: TR 38.901 Πίν. 7.8-1 (στοιχείο 8 dBi + στοιχειοκεραία 10 -> 18 dBi).
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
simParameters.BsImageGeo = [];
if opts.wrapAround
    simParameters.BsImageGeo = bsImageGeo;
end
% 3 τομείς ανά site στις 30°, 150°, 270° (TR 38.901 Πίν. 7.8-1· ITU-R M.2412-0
% §8.3). Κελί c = 3(site-1) + τομέας, κόμβος "BSc". Κεραία: TR 38.901 Πίν.
% 7.3-1 (στοιχείο) και 7.8-1 ((M,N,P) = (10,1,1), d_V = 0.5λ, κλίση 102°).
switch opts.bsSectors
    case 3
        sectorBoresightDeg = [30; 150; 270];
    case 1
        sectorBoresightDeg = [];
    otherwise
        error('runSimulation:Sectors', 'bsSectors: 1 ή 3.');
end
numSec   = opts.bsSectors;
numCells = numBs * numSec;
cell_geo = repelem(bs_geo, numSec, 1);
simParameters.CellSite         = repelem((1:numBs)', numSec, 1);
simParameters.CellBoresightDeg = repmat(sectorBoresightDeg, numBs, 1);
simParameters.BsAntenna = struct('ElementGainDbi', simParameters.AntennaGain, ...
    'Theta3dBDeg', 65, 'Phi3dBDeg', 65, 'SlaVDb', 30, 'AmDb', 30, ...
    'NumVertElements', simParameters.NumAntennaElements, 'VertSpacingLambda', 0.5, ...
    'ElecTiltDeg', 102);
simParameters.PathLoss.EnvironmentHeight = 1;
simParameters.TxPower = 49;                 % dBm, conducted (TR 38.901 Πίν. 7.8-1, UMa)
simParameters.EIRP = simParameters.TxPower + simParameters.AntennaGain + ...
                     10*log10(simParameters.NumAntennaElements);   % dBm, composite

% Μοντέλο ισχύος EARTH (Auer et al. 2011, Πίν. 2).
simParameters.Power.NumTrx = 4;   % αλυσίδες ανά τομέα
simParameters.Power.P0     = 130;
simParameters.Power.DeltaP = 4.7;
simParameters.Power.Psleep = 75;

% Κινητικότητα: Event A3 (TS 38.331 §5.5.4.4), ping-pong με MTS (TR 36.839).
simParameters.Mobility.HysteresisDb    = opts.hysteresisDb;
simParameters.Mobility.TimeToTriggerMs = opts.timeToTriggerMs;
simParameters.Mobility.DtSeconds       = opts.dtSeconds;
simParameters.Mobility.PingPongMtsS    = opts.pingPongMtsS;
% Soft LOS (TR 38.901 §7.6.3.3): αποφεύγει τα άλματα LOS <-> NLOS με κινητικότητα.
simParameters.SoftLos = opts.softLosState;
% RLF (TS 38.331 §5.3.10, τιμές TR 36.839 Πίν. 5.2.1.3.1).
simParameters.Mobility.Rlf = struct('QoutDb', opts.rlfQoutDb, 'QinDb', opts.rlfQinDb, ...
    'T310Ms', opts.rlfT310Ms, 'N310', opts.rlfN310, 'N311', opts.rlfN311);
% MLB (TS 38.300 §15.5.1), πολιτικές *MLB: Ocn = MlbMaxOffsetDb x διαφορά φόρτου.
simParameters.Mobility.MlbMaxOffsetDb  = opts.mlbMaxOffsetDb;
simParameters.Mobility.LoadReportPeriodMs = opts.mlbReportPeriodMs;
% Απόφαση handover με ML στον σταθμό (TS 38.300 §9.2.3.2.1, TR 37.817 §5.3),
% πολιτικές ActualML και ActualMLBML· χωρίς αρχείο μοντέλου δεν τρέχουν.
mlGateInfo = struct();
simParameters.Mobility.MlGate = [];
if ~isempty(opts.mlGatePath) && isfile(opts.mlGatePath)
    [simParameters.Mobility.MlGate, mlGateInfo] = loadHandoverGate(opts.mlGatePath, ...
        max(1, ceil(opts.mlReportIntervalMs/(opts.dtSeconds*1000))));
end
hasMlGate = ~isempty(simParameters.Mobility.MlGate);

% Δορυφόρος: TR 38.821 Set-1, LEO-600, S-band (Πίν. 6.1.1.1-1, 6.1.3.2-1).
% Δεδομένο η πυκνότητα EIRP· εύρος ανά δέσμη 30/FRF MHz (Πίν. 6.1.1.1-5).
satAltitude = 600e3;
satParameters.CarrierFrequency = 2.0e9;
satParameters.FreqReuse = opts.satFreqReuse;
satParameters.Bandwidth = 30e6 / satParameters.FreqReuse;
satParameters.AntennaGain = 30;
satParameters.EirpDensityDbwPerMHz = 34;
satParameters.EIRP = satParameters.EirpDensityDbwPerMHz + 10*log10(satParameters.Bandwidth/1e6) + 30;
satParameters.TxPower = satParameters.EIRP - satParameters.AntennaGain;
satParameters.MinElevationDeg = 20;

satParameters.Power.Pfix  = 0;       % W, ισχύς εκτός ενισχυτή (αισιόδοξη υπόθεση)
satParameters.Power.EtaPA = 0.4;

%% ------------------ Πεδίο ισχύος επίγειων μοντέλων ------------------
% TR 38.901 Πίν. 7.4.1-1, ελέγχεται σε κάθε βήμα.
d2dMinM = 10;
d2dMaxM = 5000;

% Χρήστες ομοιόμορφα στις κυψέλες (TR 38.901 Πίν. 7.2-1), με διαδρομή >= 35 m
% από κάθε σταθμό· όλη η διαδρομή πρέπει να μένει στο πεδίο ισχύος.
ueSpeedMps   = opts.ueSpeedKmh / 3.6;
maxWalkM     = ueSpeedMps * opts.passWindowS;
userRMinM    = 35;
userRMaxM    = cellRadiusM + maxWalkM;
if userRMinM >= userRMaxM || 2*bsClusterRadiusM + cellRadiusM + maxWalkM > d2dMaxM
    error('runSimulation:WindowTooLong', ...
        ['Με ταχύτητα %.2f m/s και παράθυρο %.0f s ο χρήστης διανύει %.0f m, ' ...
         'οπότε σε πλέγμα %d δακτυλίων με ISD %g m κάποια ζεύξη θα έβγαινε από ' ...
         'το πεδίο ισχύος [%d, %d] m. Μείωσε το passWindowS ή το bsRings.'], ...
        ueSpeedMps, opts.passWindowS, maxWalkM, opts.bsRings, opts.bsIsdM, d2dMinM, d2dMaxM);
end

%% ------------------ Μέτρηση ------------------
% (α) M_n = μέσος όρος στην περίοδο μέτρησης 200 ms (TS 38.133 Πίν. 9.2.5.2-1),
%     με ισοδύναμο πλήθος ανεξάρτητων δειγμάτων από το μοντέλο Clarke.
% (β) Φίλτρο L3 (TS 38.331 §5.5.3.2): a = 2^(-k/4) ανά περίοδο X, προσαρμοσμένο
%     στο βήμα: a' = 1 - (1-a)^(Δt/X)· k = 0 -> χωρίς φίλτρο.
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
    l3Alpha = 1;
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

% Εγκάρσια μετατόπιση ίχνους ανά διέλευση, από ζενιθιακή έως οριακή (μέγιστη
% ανύψωση = μάσκα).
raanSpanRad = findMaxRaanOffset(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);

% numSats διαδοχικοί δορυφόροι στο ίδιο επίπεδο· ο μεσαίος στη μέγιστη
% προσέγγιση. Απόσταση ώστε ο επόμενος (μη προσομοιωμένος) να ανατέλλει όταν
% κλείνει το παράθυρο: Δt = (W/2 + Tvis/2) / ((N+1)/2).
halfVisS = zenithHalfVisibilityS(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);
if opts.numSats > 1
    satSpacingS = (opts.passWindowS/2 + halfVisS) / ((opts.numSats+1)/2);
else
    satSpacingS = 0;
end
satUOffsetRad = ((opts.numSats+1)/2 - (1:opts.numSats)') * meanMotion*satSpacingS;

% Στρωματοποιημένη μετατόπιση: τα μέσα numPasses ίσων διαστημάτων του [-1, 1],
% σε τυχαία σειρά (δική της γεννήτρια).
if opts.geometryGrid
    gridStream = RandStream('mt19937ar', 'Seed', opts.rngSeed);
    raanOffsetFrac = -1 + (2*randperm(gridStream, opts.numPasses) - 1) / opts.numPasses;
end

%% ------------------ Δείκτες σύγκλισης ------------------
% Ανά διέλευση (η ανεξάρτητη μονάδα). Το μέσο SINR καταγράφεται χωρίς στόχο.
kpiNames = {'MeanCapacity_Mbps','MeanSinr_dB','FracServed','FracSatellite','MeanBitPerJouleRf'};
kpiIsFraction = [false false true true false];
kpiHasTarget = [true false true true true];
numKpi = numel(kpiNames);
passKpi = nan(opts.numPasses, numKpi);

%% ------------------ Πολιτικές ------------------
% Όλες οι πολιτικές στην ίδια πραγματοποίηση καναλιού: οι διαφορές τους είναι
% κατά ζεύγη. *MLB = με MLB, *ML = με την απόφαση ML (ίδια σειρά με τη simulateScenario).
policyNames   = {'Actual','TerrestrialOnly','SatelliteOnly','ActualMLB','TerrestrialOnlyMLB'};
if hasMlGate
    policyNames = [policyNames, {'ActualML','ActualMLBML'}];
end
policyMetrics = {'MeanCapacity_Mbps','TotalRate_Mbps','FracServed', ...
                 'FracBelowTarget','FracOutage','BitPerJouleRf','Spectrum_MHz', ...
                 'HoRate_perUserMin','HoRateTerrTerr_perUserMin', ...
                 'HoRateTerrSat_perUserMin','HoRateSatSat_perUserMin', ...
                 'PingPongRate_perUserMin','RlfRate_perUserMin', ...
                 'RlfRateTerr_perUserMin','RlfRateSat_perUserMin', ...
                 'HoFailRate_perUserMin','FracSatellite'};
if hasMlGate
    policyMetrics{end+1} = 'HoHeldRate_perUserMin';
end
numPolicy       = numel(policyNames);
numPolicyMetric = numel(policyMetrics);
policyPassKpi   = nan(opts.numPasses, numPolicy, numPolicyMetric);

convRows = cell(opts.numPasses,1);
checksRun = 0;
gateChecksRun = 0;   % handovers *ML με ελεγμένα χαρακτηριστικά
totalRows = 0;

fprintf('Ενιαία προσομοίωση ISTN\n');
fprintf('  Δt = %g s | ταχύτητα χρήστη = %g km/h (%.3f m/s) | παράθυρο διέλευσης = %g s\n', ...
    opts.dtSeconds, opts.ueSpeedKmh, ueSpeedMps, opts.passWindowS);
fprintf('  Μετατόπιση ανά βήμα = %.3f m (όριο M.2412 Παρ.1 §5.3.2: 1 m)\n', ueSpeedMps*opts.dtSeconds);
fprintf('  Τροχιακή περίοδος = %.0f s | ταχύτητα ίχνους = %.0f m/s\n', orbitalPeriodS, groundSpeedMps);
fprintf('  Εύρος εγκάρσιας μετατόπισης ίχνους = ±%.2f°\n', rad2deg(raanSpanRad));
fprintf('  Σταθμοί βάσης: %d sites x %d τομείς = %d κελιά (εξαγωνικό πλέγμα, ISD %g m) | χρήστες ανά διέλευση: %d-%d, ομοιόμορφα στις κυψέλες, διαδρομή >= %.0f m από κάθε σταθμό\n', ...
    numBs, numSec, numCells, opts.bsIsdM, numUsersRange(1), numUsersRange(2), userRMinM);
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

% Καταγραφή handovers (logHandoverEvents). Εκπαίδευση μόνο από Actual/ActualMLB·
% οι *ML καταγράφονται για τον έλεγχο ότι το μοντέλο βλέπει τα ίδια χαρακτηριστικά.
hoLogPolicies = {'Actual','ActualMLB'};
if hasMlGate
    hoLogPolicies = [hoLogPolicies, {'ActualML','ActualMLBML'}];
end
% Ιστορικό μέτρησης 2 s (δείγματα t, t-1 s, t-2 s· Aggarwal et al.).
hoHistW       = 2*round(1/opts.dtSeconds) + 1;
hoEventCells  = cell(opts.numPasses,1);

if ~isempty(resumeState)
    passIdx     = resumeState.passIdx;
    passKpi(1:passIdx,:) = resumeState.passKpi;
    policyPassKpi(1:passIdx,:,:) = resumeState.policyPassKpi;
    convRows(1:passIdx)  = resumeState.convRows;
    totalRows     = resumeState.totalRows;
    checksRun     = resumeState.checksRun;
    gateChecksRun = resumeState.gateChecksRun;
    refSnapshot   = resumeState.refSnapshot;
    refPassTable  = resumeState.refPassTable;
    hoEventCells(1:passIdx) = resumeState.hoEventCells;
end

%% ------------------ Παράλληλη εκτέλεση διελεύσεων ------------------
% Οι διελεύσεις είναι ανεξάρτητες (δικός τους σπόρος): τρέχουν σε ομάδες με
% parfor και γράφονται με τη σειρά τους, οπότε τα αρχεία είναι ίδια με τη
% σειριακή εκτέλεση.
if ~opts.geometryGrid
    raanOffsetFrac = [];
end
passCtx = struct();
passCtx.a                = a;
passCtx.bsCenterGeo      = bsCenterGeo;
passCtx.bsClusterRadiusM = bsClusterRadiusM;
passCtx.bsImageGeo       = bsImageGeo;
passCtx.bs_geo           = bs_geo;
passCtx.cLight           = cLight;
passCtx.cellRadiusM      = cellRadiusM;
passCtx.cell_geo         = cell_geo;
passCtx.d2dMaxM          = d2dMaxM;
passCtx.d2dMinM          = d2dMinM;
passCtx.ghost_geo        = ghost_geo;
passCtx.hoHistW          = hoHistW;
passCtx.hoLogPolicies    = hoLogPolicies;
passCtx.inclDeg          = inclDeg;
passCtx.maxWalkM         = maxWalkM;
passCtx.meanMotion       = meanMotion;
passCtx.numCells         = numCells;
passCtx.numKpi           = numKpi;
passCtx.numPolicy        = numPolicy;
passCtx.numPolicyMetric  = numPolicyMetric;
passCtx.numUsersRange    = numUsersRange;
passCtx.omegaEarth       = omegaEarth;
passCtx.opts             = opts;
passCtx.policyMetrics    = policyMetrics;
passCtx.policyNames      = policyNames;
passCtx.raanBase         = raanBase;
passCtx.raanOffsetFrac   = raanOffsetFrac;
passCtx.raanSpanRad      = raanSpanRad;
passCtx.satParameters    = satParameters;
passCtx.satUOffsetRad    = satUOffsetRad;
passCtx.simParameters    = simParameters;
passCtx.tPeakS           = tPeakS;
passCtx.uPeakRad         = uPeakRad;
passCtx.ueHeightM        = ueHeightM;
passCtx.ueSpeedMps       = ueSpeedMps;
passCtx.userRMinM        = userRMinM;
passCtx.wgs84            = wgs84;
numWorkers = 0;
if opts.numWorkers > 0 && ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox')
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', opts.numWorkers);
    end
    numWorkers = pool.NumWorkers;
end
fprintf('  Εκτέλεση διελεύσεων: %s\n\n', ternary(numWorkers > 0, ...
    sprintf('παράλληλα, %d workers', numWorkers), 'σειριακά'));

mlGateSha256 = '';
if hasMlGate
    mlGateSha256 = mlGateInfo.Sha256;
end

while passIdx < opts.numPasses
    batch = passIdx+1 : min(passIdx + max(numWorkers,1), opts.numPasses);
    res = cell(numel(batch), 1);
    if numWorkers > 0
        parfor jB = 1:numel(batch)
            res{jB} = runOnePass(batch(jB), passCtx);
        end
    else
        for jB = 1:numel(batch)
            res{jB} = runOnePass(batch(jB), passCtx);
        end
    end

    for jB = 1:numel(batch)
        passIdx = batch(jB);
        R = res{jB};
        res{jB} = [];
        passKpi(passIdx,:)          = R.passKpi;
        policyPassKpi(passIdx,:,:)  = R.policyPassKpi;
        hoEventCells{passIdx}       = R.hoEvents;
        checksRun                   = checksRun + R.checks;
        gateChecksRun               = gateChecksRun + R.gateChecks;
        if passIdx == 1
            refSnapshot = R.refSnapshot;
        end

        % -- Εγγραφή --
        passTable = R.passTable;
        if passIdx == 1
            writetable(passTable, opts.datasetPath);
            refPassTable = R.refPassTable;
        else
            writetable(passTable, opts.datasetPath, 'WriteMode', 'append');
        end
        totalRows = totalRows + height(passTable);
        clear passTable R;

        % -- 95% CI του μέσου όρου ως τώρα (καμπύλη σύγκλισης) --
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

        % Σημείο ελέγχου μετά από κάθε διέλευση.
        ckpt = struct('passIdx', passIdx, 'passKpi', passKpi(1:passIdx,:), ...
            'policyPassKpi', policyPassKpi(1:passIdx,:,:), ...
            'convRows', {convRows(1:passIdx)}, 'totalRows', totalRows, ...
            'checksRun', checksRun, 'gateChecksRun', gateChecksRun, 'refSnapshot', refSnapshot, ...
            'refPassTable', refPassTable, 'hoEventCells', {hoEventCells(1:passIdx)}, ...
            'opts', opts, 'mlGateSha256', mlGateSha256);
        save(checkpointPath, 'ckpt');
    end
end

numPasses = passIdx;

%% ------------------ Σύγκριση πολιτικών κατά ζεύγη ------------------
% Η στατιστική γίνεται στη διαφορά ανά διέλευση (ίδιο κανάλι, ίδια τοπολογία).
pol = policyPassKpi(1:numPasses,:,:);

% Παράγωγοι λόγοι ανά διέλευση (λόγος μέσων, όχι μέσος λόγων, για να
% αποφεύγεται το 0/0 ανά βήμα): φασματική απόδοση, ποσοστό ping-pong,
% ποσοστό αποτυχίας HO (TR 36.839 Ορισμός 4).
iRate = find(strcmp(policyMetrics, 'TotalRate_Mbps'), 1);
iSpec = find(strcmp(policyMetrics, 'Spectrum_MHz'), 1);
specEff = pol(:,:,iRate) ./ pol(:,:,iSpec);
specEff(~isfinite(specEff)) = NaN;
pol(:,:,end+1) = specEff;
policyMetrics{end+1} = 'SpectralEff_bpsHz';

iHo = find(strcmp(policyMetrics, 'HoRate_perUserMin'), 1);
iPp = find(strcmp(policyMetrics, 'PingPongRate_perUserMin'), 1);
ppFrac = pol(:,:,iPp) ./ pol(:,:,iHo);
ppFrac(~isfinite(ppFrac)) = NaN;
pol(:,:,end+1) = ppFrac;
policyMetrics{end+1} = 'PingPongFrac';

iHof = find(strcmp(policyMetrics, 'HoFailRate_perUserMin'), 1);
hofFrac = pol(:,:,iHof) ./ (pol(:,:,iHof) + pol(:,:,iHo));
hofFrac(~isfinite(hofFrac)) = NaN;
pol(:,:,end+1) = hofFrac;
policyMetrics{end+1} = 'HoFailFrac';
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

% Ζεύγη (πολιτική, βάση): μεικτή έναντι αποκλειστικών, συνεισφορά MLB,
% συνεισφορά ML.
compNames = {'Actual','TerrestrialOnly'; 'Actual','SatelliteOnly'; ...
             'ActualMLB','Actual'; 'TerrestrialOnlyMLB','TerrestrialOnly'; ...
             'ActualMLB','TerrestrialOnlyMLB'};
if hasMlGate
    compNames = [compNames; {'ActualML','Actual'; 'ActualMLBML','ActualMLB'; ...
                             'ActualMLBML','Actual'}];
end
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
iRlf = find(strcmp(policyMetrics, 'RlfRate_perUserMin'), 1);
iHofFrac = find(strcmp(policyMetrics, 'HoFailFrac'), 1);
iSatFrac = find(strcmp(policyMetrics, 'FracSatellite'), 1);
fprintf('    %-18s %12s %12s %12s %10s %14s %10s %14s %10s\n', 'πολιτική', 'ρυθμ.[Mbps]', 'εξυπηρ.', 'εκτός', ...
    'δορυφ.', 'HO/χρήστη/min', 'ping-pong', 'RLF/χρήστη/min', 'αποτ. HO');
for pI = 1:numPolicy
    hoM = mean(pol(:,pI,iHo),'omitnan');
    fprintf('    %-18s %12.2f %12.4f %12.4f %9.1f%% %14.3f %9.1f%% %14.3f %9.1f%%\n', policyNames{pI}, ...
        mean(pol(:,pI,2),'omitnan'), mean(pol(:,pI,3),'omitnan'), ...
        mean(pol(:,pI,5),'omitnan'), 100*mean(pol(:,pI,iSatFrac),'omitnan'), hoM, ...
        100*mean(pol(:,pI,iPpFrac),'omitnan'), mean(pol(:,pI,iRlf),'omitnan'), ...
        100*mean(pol(:,pI,iHofFrac),'omitnan'));
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
mlPng = '';
if hasMlGate
    mlPng = fullfile(opts.outputDir, 'ml_gain.png');
    plotMlGain(pairedTable, numPasses, mlGateInfo, opts.mlReportIntervalMs, mlPng);
    fprintf('  Απόφαση handover με ML: %d handovers με χαρακτηριστικά ίδια με την καταγραφή\n', ...
        gateChecksRun);
end

%% ------------------ Καμπύλη σύγκλισης ------------------
convMat = vertcat(convRows{1:numPasses});
convVarNames = [{'PassID','CumulativeRows'}, ...
    strcat('Mean_', kpiNames), strcat('RelCI_', kpiNames)];
convTable = array2table(convMat, 'VariableNames', convVarNames);

%% ------------------ Ακρίβεια στο τέλος της εκτέλεσης ------------------
% Αν ο στόχος CI πιάστηκε, και από ποια διέλευση ισχύει αδιάλειπτα.
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
    printUserTable(size(refSnapshot.user_geo,1), refSnapshot.bestNodeVec, refSnapshot.bestNodeTypeVec, ...
        refSnapshot.bestDistanceVec, refSnapshot.bestPathLossVec, refSnapshot.bestSinrDbVec, ...
        refSnapshot.capacityMbpsVec, refSnapshot.bestElevationDegVec, ...
        refSnapshot.nodePowerWattsVec, refSnapshot.energyPerBitUJVec);
    plotNetworkSnapshot(cell_geo, refSnapshot.user_geo, refSnapshot.sat_geo, wgs84, numCells, ...
        size(refSnapshot.user_geo,1), refSnapshot.bestNodeTypeVec, refSnapshot.bestNodeVec, ...
        opts.bsIsdM, simParameters.CellBoresightDeg);
    saveas(gcf, fullfile(opts.outputDir, 'network_3d.png'));
    close(gcf);
end

% Οι αριθμοί πίσω από τα γραφήματα· το dataset καταγράφεται μόνο με SHA-256.
save(resultsMat, 'refPassTable', 'convTable', 'policyTable', 'pairedTable', 'summary');
outFiles = {resultsMat, polPng, mlbPng, hoPng, fullfile(opts.outputDir,'network_3d.png')};
if ~isempty(mlPng)
    outFiles{end+1} = mlPng;
end

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

% Καμπύλη σύγκλισης.
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
    'pingPongMtsS', opts.pingPongMtsS, ...
    'pingPongMaxStay_s', max(1, ceil(opts.pingPongMtsS/opts.dtSeconds - 1e-9) - 1)*opts.dtSeconds, ...
    'hysteresisDb', opts.hysteresisDb, ...
    'timeToTriggerMs', opts.timeToTriggerMs, 'rngSeed', opts.rngSeed, ...
    'rlfQoutDb', opts.rlfQoutDb, 'rlfQinDb', opts.rlfQinDb, 'rlfT310Ms', opts.rlfT310Ms, ...
    'rlfN310', opts.rlfN310, 'rlfN311', opts.rlfN311, ...
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
runParams.pingPongMaxStay_s = max(1, ceil(opts.pingPongMtsS/opts.dtSeconds - 1e-9) - 1)*opts.dtSeconds;
runParams.softLosState = opts.softLosState;
runParams.mlbReportPeriod_ms = opts.mlbReportPeriodMs;
runParams.rlf = struct('Qout_dB', opts.rlfQoutDb, 'Qin_dB', opts.rlfQinDb, ...
    'T310_ms', opts.rlfT310Ms, 'N310', opts.rlfN310, 'N311', opts.rlfN311, ...
    'T310_steps', max(1, ceil(opts.rlfT310Ms/(opts.dtSeconds*1000))), ...
    'source', 'TR 36.839 V11.1.0 Πίν. 5.2.1.3.1· διαδικασία TS 38.331 §5.3.10');
runParams.timeToTrigger_steps = max(1, ceil(opts.timeToTriggerMs/(opts.dtSeconds*1000)));
runParams.measAveraging    = opts.measAveraging;
runParams.measNAvgTerrestrial = simParameters.Meas.NAvg;
runParams.measNAvgSatellite   = satParameters.Meas.NAvg;
runParams.l3FilterK        = opts.l3FilterK;
runParams.l3FilterAlphaAtX = ternary(opts.l3FilterK > 0, 2^(-opts.l3FilterK/4), 0);
runParams.l3FilterAlphaAtDt = l3Alpha;
runParams.measPeriod_ms    = opts.measPeriodMs;
runParams.mlGate           = mlGateInfo;
runParams.mlReportInterval_ms = opts.mlReportIntervalMs;
runParams.mlGateChecks     = gateChecksRun;
runParams.policies         = strjoin(policyNames, ', ');
runParams.policyMetrics    = strjoin(policyMetrics, ', ');
runParams.rngSeed          = opts.rngSeed;
runParams.rngScheme        = 'rng(rngSeed+passIdx, twister) στην αρχή κάθε διέλευσης· διελεύσεις παράλληλα σε ομάδες, εγγραφή με τη σειρά';
runParams.numWorkers       = numWorkers;
runParams.checksRun        = checksRun;
runParams.userArea         = 'ομοιόμορφα στις κυψέλες των σταθμών (TR 38.901 Πίν. 7.2-1)';
runParams.userMinDistBs_m  = userRMinM;
runParams.userMaxDistNearestBs_m = userRMaxM;
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
runParams.bsSectors        = opts.bsSectors;
runParams.cellSite         = simParameters.CellSite;
runParams.cellBoresight_deg = simParameters.CellBoresightDeg;
runParams.bsAntenna        = simParameters.BsAntenna;
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
    'timeStep',   'περίοδος L1->L3 200 ms (TR 36.839 Πίν. 5.3.2.1, TS 38.133 Πίν. 9.2.5.2-1)', ...
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
    'satSpacing', 'σχεδιαστική επιλογή: οι numSats δορυφόροι είναι ακριβώς όσοι του επιπέδου φαίνονται μέσα στο παράθυρο', ...
    'mlGate',     'TS 38.300 §9.2.3.2.1 βήμα 3 (απόφαση HO στον σταθμό) + TR 37.817 §5.3 (AI/ML βελτιστοποίηση κινητικότητας) + Aggarwal et al. ICSigSys 2026 (πρόβλεψη ασταθών HO)· ετικέτα ping-pong TR 36.839 Ορισμός 5', ...
    'mlReport',   'TS 38.331 ReportConfigNR reportInterval (ms120, ms240, ...): επανάληψη της αναφοράς όσο ισχύει η συνθήκη εισόδου');

saveRunVersion('runSimulation', runParams, outFiles, opts.label);

if isfile(checkpointPath)
    delete(checkpointPath);
end

end

% =====================================================================
% Τοπικές συναρτήσεις
% =====================================================================

function nEff = effectiveSamplesClarke(fDopplerHz, tStepS)
% Ισοδύναμο πλήθος ανεξάρτητων δειγμάτων διαλείψεων σε διάστημα T (Clarke,
% Bell Syst. Tech. J. 1968): ρ(τ) = J0(2π fD τ)^2 για την ισχύ και
%     1/nEff = (2/T^2) ∫_0^T (T - τ) ρ(τ) dτ.
if ~isfinite(fDopplerHz) || fDopplerHz <= 0 || tStepS <= 0
    nEff = 1;
    return;
end
rho   = @(tau) besselj(0, 2*pi*fDopplerHz*tau).^2;
ratio = 2/tStepS^2 * integral(@(tau) (tStepS - tau).*rho(tau), 0, tStepS, ...
    'AbsTol', 1e-12, 'RelTol', 1e-9);
nEff = max(1, round(1/ratio));
end

function [evNum, evStr, gateChecks] = logHandoverEvents(pols, prevDec, newDec, H, W, s, dtS, ...
    passIdx, numUsers)
% Ένα δείγμα ανά εντολή handover (όχι μόνο τις επιτυχείς: η επιτυχία είναι
% αποτέλεσμα), με χαρακτηριστικά διαθέσιμα μόνο τη στιγμή της απόφασης
% (handoverFeatures). Το αποτέλεσμα στη στήλη Outcome (TR 36.839 Ορισμός 3):
%   "HO"         επιτυχές handover
%   "HOF_T310"   η εντολή έφτασε ενώ έτρεχε ο T310
%   "HOF_TARGET" ο στόχος δεν ήταν προσβάσιμος
% Για τις *ML ελέγχεται ότι το μοντέλο είδε ακριβώς τα ίδια χαρακτηριστικά.
evNum = zeros(0, numel(handoverEventNumNames()));
evStr = strings(0, 6);
gateChecks = 0;
tNow  = (s-1)*dtS;
for pI = 1:numel(pols)
    pol = pols{pI};
    pN  = prevDec.(pol).Node;
    isAtt = newDec.(pol).Attempt ~= "";
    for u = find(isAtt)'
        f = pN(u);
        g = newDec.(pol).Attempt(u);
        switch newDec.(pol).Cause(u)
            case "HO",  outcome = "HO";
            case "HOF", outcome = "HOF_TARGET";
            otherwise,  outcome = "HOF_T310";
        end
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
        % Μετρήσεις στα t, t-1 s, t-2 s (Aggarwal et al.).
        lagStep = (W - 1)/2;
        mF = nan(1, 3);
        mT = nan(1, 3);
        for k = 1:3
            mF(k) = histMeas(H, W, s, (k-1)*lagStep, u, f);
            mT(k) = histMeas(H, W, s, (k-1)*lagStep, u, g);
        end
        ldNode = prevDec.(pol).LoadNode;
        ldRho  = prevDec.(pol).LoadRho;
        feat = handoverFeatures(struct('SinrFrom', mF, 'SinrTo', mT, ...
            'LagS', lagStep*dtS, 'NumUsers', numUsers, ...
            'ReturnToPrev', g == prevNode, 'StayFrom_s', (prevDec.(pol).Stay(u) + 1)*dtS, ...
            'ElevFrom', histElev(H, W, s, 0, u, f), 'ElevTo', histElev(H, W, s, 0, u, g), ...
            'DistFrom', bsDist(H, u, f), 'DistTo', bsDist(H, u, g), ...
            'LoadFrom', sum(ldRho(ldNode == f)), 'LoadTo', sum(ldRho(ldNode == g)), ...
            'UsersOnFrom', sum(pN == f), 'UsersOnTo', sum(pN == g), ...
            'MlbOffset', newDec.(pol).AttemptOffset(u)));
        % *ML: τα χαρακτηριστικά της απόφασης πρέπει να είναι ίδια bit προς bit.
        if isfield(newDec.(pol), 'GateFeat') && any(isfinite(newDec.(pol).GateFeat(u,:)))
            gf = newDec.(pol).GateFeat(u,:);
            if ~isequaln(gf, feat)
                error('runSimulation:GateFeatureMismatch', ...
                    ['Διέλευση %d, βήμα %d, χρήστης %d (%s): τα χαρακτηριστικά της ' ...
                     'απόφασης ML διαφέρουν από τα καταγεγραμμένα (μέγιστη διαφορά %g).'], ...
                    passIdx, s, u, pol, max(abs(gf - feat), [], 'omitnan'));
            end
            gateChecks = gateChecks + 1;
        end
        row = [passIdx, s, tNow, u, feat];
        evNum(end+1,:) = row;                                     %#ok<AGROW>
        evStr(end+1,:) = [string(pol), f, g, prevNode, hoType, outcome]; %#ok<AGROW>
    end
end
end

function ev = finishHandoverEvents(num, str, nodeSeq, causeSeq, dtS)
% Ετικέτες από το τι ακολούθησε κάθε handover (μόνο για την εκπαίδευση):
%   StayTo_s   παραμονή στον νέο κόμβο (s)
%   NextNode   επόμενος κόμβος ("None" = outage, "" = τέλος διέλευσης)
%   NextCause  "HO", "RLF" ή "HOF" ("" = τέλος διέλευσης)
%   Censored   η διέλευση τελείωσε πριν την αποχώρηση
% Ping-pong (TR 36.839 Ορισμός 5): NextNode == FromNode, NextCause == "HO",
% StayTo_s <= pingPongMaxStay_s. Για αποτυχημένη εντολή: NaN / "".
n = size(num,1);
stayTo   = nan(n,1);
toEnd    = nan(n,1);
censored = false(n,1);
nextNode = strings(n,1);
nextCause = strings(n,1);
for i = 1:n
    seq = nodeSeq.(char(str(i,1)))(:, num(i,4));
    cau = causeSeq.(char(str(i,1)))(:, num(i,4));
    s   = num(i,2);
    nSteps = numel(seq);
    toEnd(i) = (nSteps - s)*dtS;
    if str(i,6) ~= "HO"
        continue;
    end
    j = find(seq(s+1:end) ~= str(i,3) | cau(s+1:end) == "RLF" | cau(s+1:end) == "HOF", 1);
    if isempty(j)
        stayTo(i)   = (nSteps - s + 1)*dtS;
        censored(i) = true;
    else
        stayTo(i)    = j*dtS;
        nextNode(i)  = seq(s+j);
        nextCause(i) = cau(s+j);
    end
end
names = handoverEventNumNames();
ev = struct();
for c = 1:numel(names)
    ev.(names{c}) = num(:,c);
end
ev.ReturnToPrev = logical(ev.ReturnToPrev);
ev.Policy   = str(:,1);
ev.FromNode = str(:,2);
ev.ToNode   = str(:,3);
ev.PrevNode = str(:,4);
ev.HoType   = str(:,5);
ev.Outcome  = str(:,6);
ev.TimeToPassEnd_s = toEnd;
ev.StayTo_s = stayTo;
ev.Censored = censored;
ev.NextNode = nextNode;
ev.NextCause = nextCause;
end

function names = handoverEventNumNames()
names = [{'PassID','Step','Time_s','UserID'}, handoverFeatures()];
end

function saveHandoverEvents(cells, path, meta)
% .mat (-v7), μία μεταβλητή ανά στήλη, για scipy.io.loadmat.
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
% Ίδιες οι παράμετροι που καθορίζουν τα αποτελέσματα (όχι π.χ. το label).
keys = {'dtSeconds','ueSpeedKmh','passWindowS','datasetStride','rngSeed', ...
        'numPasses','geometryGrid','datasetPath', ...
        'hysteresisDb','timeToTriggerMs','l3FilterK','measPeriodMs','measAveraging', ...
        'pingPongMtsS','numUsersRange','bsRings','bsIsdM','numSats','satFreqReuse', ...
        'mlbMaxOffsetDb','handoverEventsPath','wrapAround','bsSectors', ...
        'rlfQoutDb','rlfQinDb','rlfT310Ms','rlfN310','rlfN311','softLosState','mlbReportPeriodMs', ...
        'mlGatePath','mlReportIntervalMs'};
tf = true;
for k = 1:numel(keys)
    if ~isfield(a,keys{k}) || ~isfield(b,keys{k}) || ~isequal(a.(keys{k}), b.(keys{k}))
        tf = false;
        return;
    end
end
end

function h = fileSha256(path)
% SHA-256 ενός αρχείου (δεκαεξαδικό).
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
    h = systemSha256(path);   % MATLAB χωρίς JVM
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
tok = regexp(strrep(out, ' ', ''), '[0-9a-fA-F]{64}', 'match', 'once');
if ~isempty(tok)
    h = lower(tok);
end
end

function g = matlabMemoryGb()
% Μνήμη της MATLAB σε GB (μόνο Windows).
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
% Μεικτή πολιτική έναντι των αποκλειστικών: μέση διαφορά κατά ζεύγη με 95% CI.
% Θετική τιμή = καλύτερη η μεικτή. Ρυθμοί/ενέργεια σε %, ποσοστά χρηστών σε
% ποσοστιαίες μονάδες· το διαθέσιμο φάσμα αναγράφεται στον τίτλο.
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
    theme(fig, 'light');
catch
end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose');
nexttile(tl);
b = gainPanel(relVal, relHw, relLabel, '%');
xlabel('Βελτίωση [%], σχετική ως προς τη βάση');
title('Ρυθμός και ενέργεια');

nexttile(tl);
gainPanel(ppVal, ppHw, ppLabel, ' μον.');
xlabel('Βελτίωση [ποσοστιαίες μονάδες]');
title('Κάλυψη χρηστών');
gainLegend(b, baseLabel);

title(tl, sprintf(['Υπό εξέταση πολιτική έναντι των αποκλειστικών: %d διελεύσεις, 95%% CI ' ...
    '(ns: το CI περιέχει το 0)\nΔιαθέσιμο φάσμα: %+.0f%% έναντι μόνο επίγειας, ' ...
    '%+.0f%% έναντι μόνο δορυφορικής'], numPasses, specVal(1), specVal(2)));
saveas(fig, pngPath);
close(fig);
end

function plotMlbGain(pairedTable, numPasses, mlbOffsetDb, pngPath)
% Συνεισφορά του MLB. Θετική τιμή = καλύτερη η πρώτη πολιτική (για handovers
% και ping-pong: λιγότερα).
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
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose');
nexttile(tl);
b = gainPanel(relVal, relHw, relLabel, '%');
xlabel('Βελτίωση [%], σχετική ως προς τη βάση');
title('Ρυθμός, ενέργεια, handovers');

nexttile(tl);
gainPanel(ppVal, ppHw, ppLabel, ' μον.');
xlabel('Βελτίωση [ποσοστιαίες μονάδες]');
title('Κάλυψη χρηστών και ping-pong');
gainLegend(b, pairLabel);

title(tl, sprintf(['Εξισορρόπηση φόρτου (MLB, TS 38.300 §15.5.1): %d διελεύσεις, 95%% CI ' ...
    '(ns: το CI περιέχει το 0)\nμετατόπιση A3 = %g dB \\times διαφορά εικονικού φόρτου ' ...
    '(Nasri & Altman), φόρτος κατά Awada et al.'], numPasses, mlbOffsetDb));
saveas(fig, pngPath);
close(fig);
end

function plotMlGain(pairedTable, numPasses, info, reportMs, pngPath)
% Συνεισφορά της απόφασης ML. Αριστερά οι δείκτες κινητικότητας (TR 36.839),
% δεξιά handovers, RLF και ρυθμός. Θετική τιμή = καλύτερη η πρώτη πολιτική.
pairs = ["ActualML","Actual"; "ActualMLBML","ActualMLB"; "ActualMLBML","Actual"];
pairLabel = {'A3: με ML έναντι χωρίς','A3 + MLB: με ML έναντι χωρίς', ...
             'A3 + MLB + ML έναντι σκέτου A3'};
mobMetric = ["PingPongFrac","HoFailFrac","FracOutage"];
mobSign   = [-1 -1 -1];
mobLabel  = {'Μείωση ping-pong (μερίδιο HO)','Μείωση αποτυχιών HO (μερίδιο)', ...
             'Μείωση εκτός κάλυψης'};
relMetric = ["HoRate_perUserMin","RlfRate_perUserMin","MeanCapacity_Mbps","TotalRate_Mbps"];
relSign   = [-1 -1 1 1];
relLabel  = {'Μείωση handovers','Μείωση RLF','Μέση χωρητικότητα χρήστη', ...
             'Συνολική ρυθμαπόδοση'};

nP = size(pairs,1);
mobVal = nan(numel(mobMetric), nP); mobHw = mobVal;
relVal = nan(numel(relMetric), nP); relHw = relVal;
for j = 1:nP
    T = pairedTable(pairedTable.Policy == pairs(j,1), :);
    [mobVal(:,j), mobHw(:,j)] = gainMatrix(T, pairs(j,2), mobMetric, mobSign, false);
    [relVal(:,j), relHw(:,j)] = gainMatrix(T, pairs(j,2), relMetric, relSign, true);
end

fig = figure('Visible','off','Position',[100 100 1400 560]);
try
    theme(fig, 'light');
catch
end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'loose');
nexttile(tl);
b = gainPanel(mobVal, mobHw, mobLabel, ' μον.');
xlabel('Βελτίωση [ποσοστιαίες μονάδες]');
title('Ping-pong, αποτυχίες, κάλυψη');

nexttile(tl);
gainPanel(relVal, relHw, relLabel, '%');
xlabel('Βελτίωση [%], σχετική ως προς τη βάση');
title('Handovers, RLF, ρυθμός');
gainLegend(b, pairLabel);

title(tl, sprintf(['Απόφαση handover στον σταθμό με ML (TS 38.300 §9.2.3.2.1, TR 37.817 §5.3): ' ...
    '%d διελεύσεις, 95%% CI (ns: το CI περιέχει το 0)\n%d δέντρα, αναμονή αν P(ping-pong) > %.2f, ' ...
    'νέα αναφορά μετά από %d ms (reportInterval)'], numPasses, info.NumTrees, ...
    info.Threshold, reportMs));
saveas(fig, pngPath);
close(fig);
end

function plotHandoverRates(pol, policyMetrics, policyNames, numPasses, mtsS, pngPath)
% Ρυθμός επιτυχών handovers ανά πολιτική (TR 36.839 Ορισμός 6), στοιβαγμένος
% ανά τύπο, με 95% CI στο σύνολο· ping-pong, αποτυχίες και RLF στις ετικέτες.
typeKey   = {'TerrTerr','TerrSat','SatSat'};
typeLabel = {'επίγειο \leftrightarrow επίγειο','επίγειο \leftrightarrow δορυφορικό', ...
             'δορυφορικό \leftrightarrow δορυφορικό'};
labelOf   = containers.Map( ...
    {'Actual','TerrestrialOnly','SatelliteOnly','ActualMLB','TerrestrialOnlyMLB', ...
     'ActualML','ActualMLBML'}, ...
    {'υπό εξέταση','μόνο επίγεια','μόνο δορυφορική','υπό εξέταση + MLB','μόνο επίγεια + MLB', ...
     'υπό εξέταση + ML','υπό εξέταση + MLB + ML'});
polLabel  = cellfun(@(n) labelOf(n), policyNames, 'UniformOutput', false);
cols      = [0.20 0.45 0.75; 0.90 0.55 0.15; 0.45 0.65 0.30];
numPolicy = numel(policyNames);

idx    = @(name) find(strcmp(policyMetrics, name), 1);
polMean = @(name) mean(pol(:,:,idx(name)), 1, 'omitnan')';
stackV  = zeros(numPolicy, 0);
segCol  = zeros(0, 3);
segLbl  = {};
for t = 1:numel(typeKey)
    allV = polMean(['HoRate' typeKey{t} '_perUserMin']);
    if ~any(allV > 0)
        continue;
    end
    stackV(:,end+1) = allV;         %#ok<AGROW>
    segCol(end+1,:) = cols(t,:);    %#ok<AGROW>
    segLbl{end+1}   = typeLabel{t}; %#ok<AGROW>
end
totV    = pol(:,:,idx('HoRate_perUserMin'));
totM    = mean(totV, 1, 'omitnan');
totHw   = 1.96 * std(totV, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(totV), 1));
rlfM    = mean(pol(:,:,idx('RlfRate_perUserMin')), 1, 'omitnan');
hofM    = 100 * mean(pol(:,:,idx('HoFailFrac')), 1, 'omitnan');
ppM     = 100 * mean(pol(:,:,idx('PingPongFrac')), 1, 'omitnan');

fig = figure('Visible','off','Position',[100 100 max(1150, 215*numPolicy) 560]);
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
        lbl = sprintf('%.3f\nping-pong %.1f%%\nαποτυχίες %.1f%%\nRLF %.3f/min', totM(pI), ...
            ppM(pI), hofM(pI), rlfM(pI));
    else
        lbl = sprintf('%.3f\nκανένα handover\nRLF %.3f/min', totM(pI), rlfM(pI));
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
    'επιτυχή handovers του A3· ping-pong = επιστροφή με παραμονή \\leq %g s· ' ...
    'αποτυχίες και RLF κατά TR 36.839 (Qout -8 dB, Qin -6 dB, T310 1 s)'], ...
    numPasses, mtsS));
saveas(fig, pngPath);
close(fig);
end

function [val, hw] = gainMatrix(T, baselines, metrics, sgn, relative)
% (δείκτης x βάση): μέση διαφορά και ημιεύρος CI, σε % ή ποσοστιαίες μονάδες.
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
                continue;
            end
            scale = 100 / base;
        end
        val(i,j) = sgn(i) * T.MeanDiff(r) * scale;
        hw(i,j)  = T.CI95HalfWidth(r) * scale;
    end
end
end

function b = gainPanel(val, hw, labels, unit)
% Οριζόντιες μπάρες: μία ομάδα ανά δείκτη, μία μπάρα ανά βάση.
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
end

function gainLegend(b, seriesLabels)
% Ένα κοινό υπόμνημα κάτω από τα δύο panels.
lgd = legend(b, seriesLabels, 'Orientation', 'horizontal');
lgd.Layout.Tile = 'south';
end

function [G, info] = loadHandoverGate(path, reportSteps)
% Φορτώνει το μοντέλο του Model/train_handover_gate.py και ελέγχει ότι τα
% χαρακτηριστικά υπάρχουν και ότι η handoverGatePredict δίνει τις ίδιες
% πιθανότητες με το scikit-learn στα δείγματα ελέγχου.
S = load(path);
allNames = [handoverFeatures(), {'FromSat','ToSat'}];
featNames = cellstr(S.FeatureNames);
[ok, featIdx] = ismember(featNames, allNames);
if ~all(ok)
    error('runSimulation:GateFeatures', ...
        'Το μοντέλο %s ζητά στήλες που δεν υπάρχουν στην προσομοίωση: %s', ...
        path, strjoin(featNames(~ok), ', '));
end
M = struct('Roots', double(S.Roots(:)), 'Feature', double(S.Feature(:)), ...
    'Threshold', double(S.Threshold(:)), 'Left', double(S.Left(:)), ...
    'Right', double(S.Right(:)), 'Value', double(S.Value(:)), ...
    'MaxDepth', double(S.MaxDepth), 'Sentinel', double(S.Sentinel));
testX = double(S.TestX);
testP = double(S.TestP(:));
pM = zeros(size(testP));
for i = 1:numel(testP)
    pM(i) = handoverGatePredict(M, testX(i,:));
end
maxErr = max(abs(pM - testP));
if ~(maxErr <= 1e-12)
    error('runSimulation:GatePredict', ...
        ['Η handoverGatePredict διαφέρει από το scikit-learn στα δείγματα ' ...
         'ελέγχου του %s (μέγιστη διαφορά %g).'], path, maxErr);
end
G = struct('Model', M, 'Threshold', double(S.Threshold_p), ...
    'FeatIdx', featIdx(:)', 'ReportSteps', reportSteps);
info = struct('Path', path, 'Sha256', fileSha256(path), ...
    'Threshold', G.Threshold, 'NumTrees', numel(M.Roots), ...
    'Features', strjoin(featNames, ', '), 'ReportSteps', reportSteps, ...
    'SelfTestSamples', numel(testP), 'SelfTestMaxErr', maxErr);
if isfield(S, 'Meta')
    info.Meta = S.Meta;
end
fprintf('  Μοντέλο απόφασης handover: %s\n', path);
fprintf('    %d δέντρα, %d χαρακτηριστικά, κατώφλι %.3f, reportInterval %d βήματα\n', ...
    numel(M.Roots), numel(featIdx), G.Threshold, reportSteps);
fprintf('    ίδιες πιθανότητες με το scikit-learn σε %d δείγματα (μέγιστη διαφορά %.1e)\n', ...
    numel(testP), maxErr);
end

function R = runOnePass(passIdx, C)
% Μία διέλευση. Εξαρτάται μόνο από το passIdx και το πλαίσιο C, οπότε δίνει
% τα ίδια αποτελέσματα σειριακά ή σε worker.
a              = C.a;
bsCenterGeo      = C.bsCenterGeo;
bsClusterRadiusM = C.bsClusterRadiusM;
bsImageGeo       = C.bsImageGeo;
bs_geo           = C.bs_geo;
cLight           = C.cLight;
cellRadiusM      = C.cellRadiusM;
cell_geo         = C.cell_geo;
d2dMaxM          = C.d2dMaxM;
d2dMinM          = C.d2dMinM;
ghost_geo        = C.ghost_geo;
hoHistW          = C.hoHistW;
hoLogPolicies    = C.hoLogPolicies;
inclDeg          = C.inclDeg;
maxWalkM         = C.maxWalkM;
meanMotion       = C.meanMotion;
numCells         = C.numCells;
numKpi           = C.numKpi;
numPolicy        = C.numPolicy;
numPolicyMetric  = C.numPolicyMetric;
numUsersRange    = C.numUsersRange;
omegaEarth       = C.omegaEarth;
opts             = C.opts;
policyMetrics    = C.policyMetrics;
policyNames      = C.policyNames;
raanBase         = C.raanBase;
raanOffsetFrac   = C.raanOffsetFrac;
raanSpanRad      = C.raanSpanRad;
satParameters    = C.satParameters;
satUOffsetRad    = C.satUOffsetRad;
simParameters    = C.simParameters;
tPeakS           = C.tPeakS;
uPeakRad         = C.uPeakRad;
ueHeightM        = C.ueHeightM;
ueSpeedMps       = C.ueSpeedMps;
userRMinM        = C.userRMinM;
wgs84            = C.wgs84;
checksRun   = 0;
gateChecks  = 0;
refSnapshot = struct();
passKpi     = nan(1, numKpi);

% Ένας σπόρος ανά διέλευση (ρητά twister: οι workers έχουν άλλη προεπιλογή).
rng(opts.rngSeed + passIdx, 'twister');

% -- Χρήστες: θέση ομοιόμορφη στις κυψέλες, τυχαία κατεύθυνση (M.2412-0 §8.4) --
nU = randi(numUsersRange);
user_geo = zeros(nU,3);
walkAzimuthDeg = zeros(nU,1);
for u = 1:nU
    [uLat, uLon, walkAzimuthDeg(u)] = randUserAndHeading(bsCenterGeo, bs_geo, ghost_geo, ...
        reshape(bsImageGeo, [], 2), wgs84, bsClusterRadiusM + cellRadiusM, userRMinM, maxWalkM);
    user_geo(u,:) = [uLat, uLon, ueHeightM];
end
numUsers = size(user_geo,1);

% -- Γεωμετρία της συγκεκριμένης διέλευσης --
if opts.geometryGrid
    raanRad = raanBase + raanOffsetFrac(passIdx)*raanSpanRad;
else
    raanRad = raanBase + (2*rand()-1)*raanSpanRad;
end
tStart  = tPeakS - opts.passWindowS/2;
numSteps = floor(opts.passWindowS/opts.dtSeconds) + 1;

% Προδεσμευμένοι πίνακες· ο table φτιάχνεται μία φορά στο τέλος.
nRowsPass = numSteps * numUsers;
numBuf = zeros(nRowsPass, 25);
strBuf = strings(nRowsPass, 6);
rowPtr = 0;

prevServingNode = strings(numUsers,1);
channelState = [];
decisionState = [];
stepBitPerJouleRf = nan(numSteps,1);
policyStepBuf = nan(numSteps, numPolicy, numPolicyMetric);

% Κυκλικό ιστορικό μετρήσεων και ακολουθία εξυπηρετούντων, για τα handovers.
hoHist = struct('Bs', nan(hoHistW, numUsers, numCells), ...
    'Sat', nan(hoHistW, numUsers, opts.numSats), ...
    'Elev', nan(hoHistW, numUsers, opts.numSats), 'Dist', nan(numUsers, numCells));
hoNodeSeq = struct();
hoCauseSeq = struct();
for pL = hoLogPolicies
    hoNodeSeq.(pL{1})  = strings(numSteps, numUsers);
    hoCauseSeq.(pL{1}) = strings(numSteps, numUsers);
end
hoPassNum = zeros(0, numel(handoverEventNumNames()));
hoPassStr = strings(0, 6);

for s = 1:numSteps
    t = tStart + (s-1)*opts.dtSeconds;
    sat_geo = zeros(opts.numSats,3);
    for k = 1:opts.numSats
        sat_geo(k,:) = orbitPositionLla(uPeakRad + meanMotion*(t - tPeakS) + ...
            satUOffsetRad(k), inclDeg, raanRad, a, omegaEarth, t);
    end

    % --- Έλεγχος κατάστασης χρηστών ---
    checksRun = checksRun + 1;
    validateUserState(user_geo, bsImageGeo, wgs84, ueHeightM, d2dMinM, d2dMaxM, ...
        userRMinM, passIdx, s);

    prevDecisionState = decisionState;
    [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
        bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
        nodePowerWattsVec, energyPerBitUJVec, ...
        bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
        satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
        channelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
        bsReasonVec, satReasonVec, policyKpis, decisionState] = ...
        simulateScenario(cell_geo, user_geo, sat_geo, wgs84, simParameters, ...
                         satParameters, channelState, prevDecisionState);

    % --- Καταγραφή handovers του βήματος ---
    slot = mod(s-1, hoHistW) + 1;
    hoHist.Bs(slot,:,:)   = channelState.FiltSinrBs_dB;
    hoHist.Sat(slot,:,:)  = channelState.FiltSinrSat_dB;
    hoHist.Elev(slot,:,:) = channelState.SatElev_deg;
    hoHist.Dist           = channelState.BsDist_m;
    for pL = hoLogPolicies
        hoNodeSeq.(pL{1})(s,:)  = decisionState.(pL{1}).Node';
        hoCauseSeq.(pL{1})(s,:) = decisionState.(pL{1}).Cause';
    end
    if ~isempty(prevDecisionState)
        [evNum, evStr, nGate] = logHandoverEvents(hoLogPolicies, prevDecisionState, ...
            decisionState, hoHist, hoHistW, s, opts.dtSeconds, passIdx, numUsers);
        gateChecks = gateChecks + nGate;
        hoPassNum = [hoPassNum; evNum]; %#ok<AGROW>
        hoPassStr = [hoPassStr; evStr]; %#ok<AGROW>
    end

    % --- Δείκτες πολιτικών του βήματος ---
    for pI = 1:numPolicy
        ps = policyKpis.(policyNames{pI});
        for mI = 1:numPolicyMetric
            policyStepBuf(s, pI, mI) = ps.(policyMetrics{mI});
        end
    end

    % --- Φορτίο κόμβου ---
    [~, ~, nodeGrp] = unique(bestNodeVec);
    nodeCount   = accumarray(nodeGrp(:), 1);
    nodeLoadVec = nodeCount(nodeGrp(:));
    nodeLoadVec(bestNodeTypeVec == "Outage") = 0;

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
        user_geo(:,3) = ueHeightM;
    end
end

% -- Δείκτες της διέλευσης, σε πλήρη χρονική ανάλυση --
servedFull = strBuf(:,2) ~= "Outage";
passKpi(1,1) = mean(numBuf(servedFull,12), 'omitnan');   % Capacity_Mbps
passKpi(1,2) = mean(numBuf(servedFull,11), 'omitnan');   % SINR_dB
passKpi(1,3) = mean(strBuf(:,3) == "Served");
passKpi(1,4) = mean(strBuf(:,2) == "Satellite");
passKpi(1,5) = mean(stepBitPerJouleRf, 'omitnan');

% Μία παρατήρηση ανά διέλευση (τα βήματα δεν είναι ανεξάρτητα).
policyPassRow = mean(policyStepBuf, 1, 'omitnan');

% -- Handovers της διέλευσης: τι ακολούθησε το καθένα --
hoEvents = finishHandoverEvents(hoPassNum, hoPassStr, ...
    hoNodeSeq, hoCauseSeq, opts.dtSeconds);

% -- Γραμμές του αρχείου: βήματα 1, 1+stride, ...· η διέλευση 1 και πλήρης --
keep = mod(numBuf(:,2) - 1, opts.datasetStride) == 0;
R = struct();
R.passTable    = buildPassTable(numBuf(keep,:), strBuf(keep,:));
R.refPassTable = table();
if passIdx == 1
    R.refPassTable = buildPassTable(numBuf, strBuf);
end
R.passKpi       = passKpi;
R.policyPassKpi = policyPassRow;
R.hoEvents      = hoEvents;
R.checks        = checksRun;
R.gateChecks    = gateChecks;
R.refSnapshot   = refSnapshot;
end

function T = buildPassTable(numBuf, strBuf)
% Πίνακας μιας διέλευσης από τους αριθμητικούς (numBuf) και τους κειμενικούς
% (strBuf) πίνακες.
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

function validateUserState(user_geo, bsImageGeo, wgs84, ueHeightM, d2dMin, d2dMax, minDistBs, passIdx, stepIdx)
% Έλεγχος χρηστών: σταθερό ύψος, έγκυρες συντεταγμένες, ζεύξεις στο πεδίο
% ισχύος UMa (TR 38.901 Πίν. 7.4.1-1).
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
% Απόσταση από το πλησιέστερο αντίγραφο κάθε σταθμού, όπως στο κανάλι.
imgLatRow = reshape(bsImageGeo(:,:,1), 1, []);
imgLonRow = reshape(bsImageGeo(:,:,2), 1, []);
nPts = numel(imgLatRow);
dAll = distance(repmat(user_geo(:,1), 1, nPts), repmat(user_geo(:,2), 1, nPts), ...
    repmat(imgLatRow, numUsers, 1), repmat(imgLonRow, numUsers, 1), wgs84);
d2dAll = min(reshape(dAll, numUsers, size(bsImageGeo,1), []), [], 3);
for u = 1:numUsers
    d2d = d2dAll(u,:)';
    if min(d2d) < minDistBs - 1e-6
        error('runSimulation:TooCloseToBs', ...
            'Διέλευση %d, βήμα %d, χρήστης %d: %.2f m από σταθμό (< %g m).', ...
            passIdx, stepIdx, u, min(d2d), minDistBs);
    end
    if ~any(d2d >= d2dMin & d2d <= d2dMax)
        error('runSimulation:NoValidLink', ...
            ['Διέλευση %d, βήμα %d, χρήστης %d: καμία επίγεια ζεύξη εντός ' ...
             '[%d, %d] m (ελάχιστη απόσταση %.1f m, μέγιστη %.1f m).'], ...
            passIdx, stepIdx, u, d2dMin, d2dMax, min(d2d), max(d2d));
    end
end
end

function validateOutputs(typeVec, capacityVec, throughputVec, nodeLoadVec, sinrVec, passIdx, stepIdx)
% Συνέπεια εξόδων: outage -> χωρίς χωρητικότητα/ρυθμό/φορτίο· συνδεδεμένος ->
% πεπερασμένη χωρητικότητα (> 0 πάνω από το MCS 0) και φορτίο >= 1.
isOut = typeVec == "Outage";
if any(isOut & (~isnan(capacityVec) | throughputVec ~= 0 | nodeLoadVec ~= 0))
    error('runSimulation:OutageInconsistent', ...
        'Διέλευση %d, βήμα %d: χρήστης σε outage με μη μηδενική χωρητικότητα/ρυθμαπόδοση/φορτίο.', ...
        passIdx, stepIdx);
end
if any(~isOut & (~isfinite(capacityVec) | capacityVec < 0))
    error('runSimulation:ServedInconsistent', ...
        'Διέλευση %d, βήμα %d: συνδεδεμένος χρήστης χωρίς πεπερασμένη μη αρνητική χωρητικότητα.', ...
        passIdx, stepIdx);
end
if any(~isOut & sinrVec >= 10*log10(2^0.2344 - 1) & capacityVec <= 0)
    error('runSimulation:ServedInconsistent', ...
        'Διέλευση %d, βήμα %d: χρήστης πάνω από το MCS 0 με μηδενική χωρητικότητα.', ...
        passIdx, stepIdx);
end
if any(~isOut & nodeLoadVec < 1)
    error('runSimulation:LoadInconsistent', ...
        'Διέλευση %d, βήμα %d: εξυπηρετούμενος χρήστης με φορτίο κόμβου < 1.', ...
        passIdx, stepIdx);
end
if any(~isOut & isnan(sinrVec) & capacityVec ~= 0)
    error('runSimulation:SinrInconsistent', ...
        'Διέλευση %d, βήμα %d: χρήστης χωρίς ζεύξη στο μοντέλο με μη μηδενική χωρητικότητα.', ...
        passIdx, stepIdx);
end
end

function [lat, lon, walkAz] = randUserAndHeading(centerGeo, bs_geo, ghost_geo, img_geo, wgs84, rDiskM, dMinM, walkM)
% Θέση ομοιόμορφη στις κυψέλες και ομοιόμορφη κατεύθυνση, με απόρριψη: δεκτή
% αν ο πλησιέστερος κόμβος είναι σταθμός (όχι του κενού δακτυλίου) και η
% διαδρομή walkM δεν περνά πιο κοντά από dMinM σε σταθμό ή αντίγραφό του.
while true
    rM = rDiskM * sqrt(rand());
    az = 360 * rand();
    [lat, lon] = reckon(centerGeo(1), centerGeo(2), rM, az, wgs84);
    if min(distance(lat, lon, bs_geo(:,1), bs_geo(:,2), wgs84)) >= ...
            min(distance(lat, lon, ghost_geo(:,1), ghost_geo(:,2), wgs84))
        continue;
    end
    [dBs, azBs] = distance(lat, lon, img_geo(:,1), img_geo(:,2), wgs84);
    walkAz = 360 * rand();
    dir = [sind(walkAz) cosd(walkAz)];
    bsXY = dBs .* [sind(azBs) cosd(azBs)];
    along = min(max(bsXY * dir', 0), walkM);
    dPath = hypot(bsXY(:,1) - along*dir(1), bsXY(:,2) - along*dir(2));
    if min(dPath) >= dMinM
        return;
    end
end
end

function offsets = hexSiteOffsetsM(numRings, isdM)
% Θέσεις [ανατολικά βόρεια] (m) εξαγωνικού πλέγματος R δακτυλίων: 1 + 3R(R+1) σταθμοί.
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
% Μέγιστη εγκάρσια μετατόπιση ίχνους που φτάνει ακόμη τη μάσκα (διχοτόμηση).
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
% Μέγιστη ανύψωση από το κέντρο της περιοχής μέσα στο παράθυρο.
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
% Μισή διάρκεια ορατότητας της ζενιθιακής διέλευσης (διχοτόμηση, μέγιστο των δύο πλευρών).
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
% Θέση σε κυκλική τροχιά -> [lat lon alt]: αδρανειακό σύστημα, στροφή κατά
% omegaEarth*t σε ECEF.
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

function printUserTable(numUsers, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec)
% Πίνακας χρηστών στη στιγμή μέγιστης ανύψωσης της διέλευσης αναφοράς.
userID = (1:numUsers)';
disp(table(userID, bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, nodePowerWattsVec, energyPerBitUJVec, ...
    'VariableNames', {'User','ServingNode','ServingType', ...
    'Distance_m','PathLoss_dB','SINR_dB','Capacity_Mbps','SatElevation_deg', ...
    'NodePower_W','EnergyPerBit_uJ'}));
end
