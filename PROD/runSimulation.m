function [convTable, summary] = runSimulation(opts)
%RUNSIMULATION Προσομοίωση ISTN (επίγειο UMa + LEO).
%   Μονάδα: μία διέλευση = ένα drop με νέους χρήστες (TR 38.821 Πίν. 4.2-3·
%   M.2412-0 §7.1). Βήμα 0.2 s. Ίδιες διελεύσεις σε όλες τις πολιτικές
%   (κατά ζεύγη, 95% CI: M.2412-0 §7.1, M.2514-0 §8.2.4).
%   Χρήστες 3 km/h, τυχαία κατεύθυνση (M.2412-0 §8.4).
%
%   runSimulation()                          % προεπιλογές
%   runSimulation(struct('numPasses',10))    % σύντομη δοκιμή
%
% Έξοδοι: Dataset/dataset.csv, Dataset/handover_events.mat,
%         Results/simulation_results.mat, Results/*.png

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

% Σημείο ελέγχου: συνέχιση αν δεν άλλαξαν παράμετροι και μοντέλο.
resumeState = [];
if opts.resume && isfile(checkpointPath) && isfile(opts.datasetPath)
    S = load(checkpointPath);
    gateSha = '';
    if ~isempty(opts.mlGatePath) && isfile(opts.mlGatePath)
        gateSha = fileSha256(opts.mlGatePath);
    end
    if isfield(S,'ckpt') && isfield(S.ckpt,'mlGateSha256') && ...
            strcmp(S.ckpt.mlGateSha256, gateSha) && checkpointMatches(S.ckpt.opts, opts)
        % Διακοπή μεταξύ CSV και σημείου ελέγχου: διπλές γραμμές, άρα σταματάμε.
        lastInCsv = max(readmatrix(opts.datasetPath, 'Range', 'A:A', ...
            'NumHeaderLines', 1));
        if lastInCsv ~= S.ckpt.passIdx
            error('runSimulation:CheckpointMismatch', ...
                ['The dataset reaches pass %d but the checkpoint pass %d: the previous ' ...
                 'run was interrupted between the two writes. Delete %s and %s and ' ...
                 'start again.'], ...
                lastInCsv, S.ckpt.passIdx, opts.datasetPath, checkpointPath);
        end
        resumeState = S.ckpt;
        fprintf('Resuming from checkpoint: %d passes already completed.\n', ...
            resumeState.passIdx);
    else
        fprintf(['Checkpoint found with different parameters: ignored, ' ...
                 'the run starts from pass 1.\n']);
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

% Επόμενος δακτύλιος (χωρίς σταθμό): ορίζει το όριο της περιοχής.
ringOffsets = hexSiteOffsetsM(opts.bsRings + 1, opts.bsIsdM);
ringOffsets = ringOffsets(numBs+1:end,:);
ghost_geo = zeros(size(ringOffsets,1), 2);
for b = 1:size(ringOffsets,1)
    [ghost_geo(b,1), ghost_geo(b,2)] = reckon(bsCenterGeo(1), bsCenterGeo(2), ...
        norm(ringOffsets(b,:)), atan2d(ringOffsets(b,1), ringOffsets(b,2)), wgs84);
end
cellRadiusM = opts.bsIsdM / sqrt(3);   % ακτίνα εξαγώνου (κορυφή)

% Wrap-around: 6 αντίγραφα του πλέγματος· κάθε ζεύξη από το πλησιέστερο.
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
% Έλεγχος πλακόστρωσης: χωρίς επικαλύψεις, ίδιος περίγυρος.
if opts.wrapAround
    allImg = reshape(reshape(hexOffsets, numBs, 1, 2) + reshape(wrapShiftM, 1, numImg, 2), [], 2);
    for g = 1:size(ringOffsets,1)
        if min(hypot(allImg(:,1) - ringOffsets(g,1), allImg(:,2) - ringOffsets(g,2))) > 1e-6
            error('runSimulation:WrapAround', 'Wrap-around: position %d of the next ring is not an image of a site.', g);
        end
    end
    dPair = hypot(allImg(:,1) - allImg(:,1)', allImg(:,2) - allImg(:,2)');
    dPair(1:size(dPair,1)+1:end) = Inf;
    if min(dPair(:)) < opts.bsIsdM - 1e-6
        error('runSimulation:WrapAround', 'Wrap-around: overlapping site images.');
    end
    p0 = [123.4 56.7];
    dRef = [];
    for b = 1:numBs
        q = hexOffsets(b,:) + p0;
        dW = min(reshape(hypot(q(1) - allImg(:,1), q(2) - allImg(:,2)), numBs, numImg), [], 2);
        if isempty(dRef)
            dRef = sort(dW);
        elseif max(abs(sort(dW) - dRef)) > 1e-6
            error('runSimulation:WrapAround', 'Wrap-around: the surroundings of site %d differ from those of the central site.', b);
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
% 3 τομείς ανά site (30°, 150°, 270°)· κελί c = 3(site-1) + τομέας.
% TR 38.901 Πίν. 7.3-1, 7.8-1· ITU-R M.2412-0 §8.3.
switch opts.bsSectors
    case 3
        sectorBoresightDeg = [30; 150; 270];
    case 1
        sectorBoresightDeg = [];
    otherwise
        error('runSimulation:Sectors', 'bsSectors must be 1 or 3.');
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
% Απόφαση HO με ML (TS 38.300 §9.2.3.2.1, TR 37.817 §5.3)· θέλει αρχείο μοντέλου.
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

% Χρήστες ομοιόμορφα στις κυψέλες (TR 38.901 Πίν. 7.2-1), >= 35 m από σταθμό.
ueSpeedMps   = opts.ueSpeedKmh / 3.6;
maxWalkM     = ueSpeedMps * opts.passWindowS;
userRMinM    = 35;
userRMaxM    = cellRadiusM + maxWalkM;
if userRMinM >= userRMaxM || 2*bsClusterRadiusM + cellRadiusM + maxWalkM > d2dMaxM
    error('runSimulation:WindowTooLong', ...
        ['At %.2f m/s over a %.0f s window a UE moves %.0f m, so with %d rings ' ...
         'and ISD %g m some link would leave the UMa validity range [%d, %d] m ' ...
         '(TR 38.901 Table 7.4.1-1). Reduce passWindowS or bsRings.'], ...
        ueSpeedMps, opts.passWindowS, maxWalkM, opts.bsRings, opts.bsIsdM, d2dMinM, d2dMaxM);
end

%% ------------------ Μέτρηση ------------------
% M_n: μέσος όρος 200 ms (TS 38.133 Πίν. 9.2.5.2-1, δείγματα κατά Clarke).
% Φίλτρο L3 (TS 38.331 §5.5.3.2) προσαρμοσμένο στο βήμα· k = 0: χωρίς φίλτρο.
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

% Εγκάρσια μετατόπιση ίχνους ανά διέλευση (ζενιθιακή έως οριακή).
raanSpanRad = findMaxRaanOffset(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);

% numSats δορυφόροι στο ίδιο επίπεδο, ο μεσαίος στη μέγιστη προσέγγιση·
% ο επόμενος ανατέλλει όταν κλείνει το παράθυρο.
halfVisS = zenithHalfVisibilityS(raanBase, uPeakRad, meanMotion, inclDeg, a, ...
    omegaEarth, tPeakS, opts.passWindowS, centerLat, centerLon, wgs84, ...
    satParameters.MinElevationDeg);
if opts.numSats > 1
    satSpacingS = (opts.passWindowS/2 + halfVisS) / ((opts.numSats+1)/2);
else
    satSpacingS = 0;
end
satUOffsetRad = ((opts.numSats+1)/2 - (1:opts.numSats)') * meanMotion*satSpacingS;

% Στρωματοποιημένη μετατόπιση, σε τυχαία σειρά.
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
% Όλες στο ίδιο κανάλι (κατά ζεύγη). *MLB: με MLB, *ML: με ML.
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

fprintf('ISTN system-level simulation: 5G NR UMa (TN) + LEO-600 (NTN)\n');
fprintf('  Time step = %g s | UE speed = %g km/h (%.3f m/s) | pass window = %g s\n', ...
    opts.dtSeconds, opts.ueSpeedKmh, ueSpeedMps, opts.passWindowS);
fprintf('  UE displacement per step = %.3f m (ITU-R M.2412-0 Annex 1 §5.3.2 limit: 1 m)\n', ueSpeedMps*opts.dtSeconds);
fprintf('  Orbital period = %.0f s | ground-track speed = %.0f m/s\n', orbitalPeriodS, groundSpeedMps);
fprintf('  Cross-track offset range of the ground track = ±%.2f°\n', rad2deg(raanSpanRad));
fprintf('  TN: %d sites x %d sectors = %d cells (hexagonal layout, ISD %g m) | UEs per pass: %d-%d, uniform over the cells, path >= %.0f m from every site\n', ...
    numBs, numSec, numCells, opts.bsIsdM, numUsersRange(1), numUsersRange(2), userRMinM);
fprintf('  NTN: %d satellites in one orbital plane, spaced %.1f s (%.2f°, i.e. %.1f per plane) | zenith-pass visibility %.0f s | FRF = %d (%.0f MHz per beam)\n', ...
    opts.numSats, satSpacingS, rad2deg(meanMotion*satSpacingS), ...
    ternary(satSpacingS > 0, orbitalPeriodS/max(satSpacingS,eps), NaN), 2*halfVisS, ...
    opts.satFreqReuse, satParameters.Bandwidth/1e6);
fprintf('  Passes: %d (fixed number) | ground-track offsets: %s\n', opts.numPasses, ...
    ternary(opts.geometryGrid, 'stratified', 'random'));
fprintf('  Precision target (ITU-R M.2412-0 §7.1, reported only): 95%% CI half-width < %.0f%% relative / %.0f percentage points\n', ...
    100*opts.ciTolerance, 100*opts.ciToleranceFrac);
fprintf('  Event A3: hysteresis = %.1f dB | time-to-trigger = %d ms (%d steps)\n', ...
    opts.hysteresisDb, opts.timeToTriggerMs, ...
    max(1, ceil(opts.timeToTriggerMs/(opts.dtSeconds*1000))));
fprintf('  Measurement: mean of %d (TN) / %d (NTN) independent samples per step | L3 filter k = %d -> a = %.3f\n', ...
    simParameters.Meas.NAvg, satParameters.Meas.NAvg, opts.l3FilterK, l3Alpha);
fprintf('  KPIs with a precision target: %s\n\n', strjoin(kpiNames(kpiHasTarget), ', '));

%% ------------------ Βρόχος διελεύσεων ------------------
passIdx = 0;
refPassTable = table();
refSnapshot = struct();

% Καταγραφή handovers· εκπαίδευση μόνο από Actual/ActualMLB.
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
% Ανεξάρτητες διελεύσεις: parfor σε ομάδες, εγγραφή με τη σειρά.
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
fprintf('  Pass execution: %s\n\n', ternary(numWorkers > 0, ...
    sprintf('parallel, %d workers', numWorkers), 'serial'));

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
            fprintf('  Pass %3d/%d | rows %8d | mean UE rate %6.2f Mbps | SE target met %5.1f%% | on NTN %5.1f%% | max CI half-width %.4f | memory %5.2f GB\n', ...
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

% Λόγοι ανά διέλευση ως λόγος μέσων (αποφεύγει το 0/0)· HOF: TR 36.839 Ορ. 4.
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

% Ζεύγη (πολιτική, βάση).
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

fprintf('\n  Policies on the same passes (%d passes, mean over passes)\n', numPasses);
iPpFrac = find(strcmp(policyMetrics, 'PingPongFrac'), 1);
iRlf = find(strcmp(policyMetrics, 'RlfRate_perUserMin'), 1);
iHofFrac = find(strcmp(policyMetrics, 'HoFailFrac'), 1);
iSatFrac = find(strcmp(policyMetrics, 'FracSatellite'), 1);
fprintf('    %-18s %12s %12s %12s %10s %14s %10s %14s %10s\n', 'Policy', 'Thr. [Mbps]', 'SE target', 'No link', ...
    'On NTN', 'HO/UE/min', 'Ping-pong', 'RLF/UE/min', 'HOF ratio');
for pI = 1:numPolicy
    hoM = mean(pol(:,pI,iHo),'omitnan');
    fprintf('    %-18s %12.2f %11.1f%% %11.2f%% %9.1f%% %14.3f %9.1f%% %14.3f %9.1f%%\n', policyNames{pI}, ...
        mean(pol(:,pI,2),'omitnan'), 100*mean(pol(:,pI,3),'omitnan'), ...
        100*mean(pol(:,pI,5),'omitnan'), 100*mean(pol(:,pI,iSatFrac),'omitnan'), hoM, ...
        100*mean(pol(:,pI,iPpFrac),'omitnan'), mean(pol(:,pI,iRlf),'omitnan'), ...
        100*mean(pol(:,pI,iHofFrac),'omitnan'));
end
fprintf('\n    Paired differences per pass (policy minus baseline), mean ± 95%% CI:\n');
for q = 1:height(pairedTable)
    signif = ternary(pairedTable.CI95Low(q) > 0 || pairedTable.CI95High(q) < 0, ...
        '', '   [CI includes 0]');
    fprintf('    %-18s - %-18s %-26s %+10.4f ± %.4f (%+6.1f%%, positive in %4.1f%% of passes)%s\n', ...
        pairedTable.Policy(q), pairedTable.Baseline(q), pairedTable.Metric(q), pairedTable.MeanDiff(q), ...
        pairedTable.CI95HalfWidth(q), pairedTable.RelDiffPct(q), ...
        100*pairedTable.FracPassesPositive(q), signif);
end
fprintf('\n');

if hasMlGate
    fprintf('  ML handover gate: features identical to the logged ones in %d handover decisions\n', ...
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

fprintf('\n--- Summary ---\n');
if ciTargetMet
    fprintf('  Precision target: MET with %d passes (continuously within target from pass %d on)\n', ...
        numPasses, ciTargetFromPass);
else
    fprintf('  Precision target: NOT met with %d passes (largest CI half-width %.4f). Increase numPasses.\n', ...
        numPasses, max(lastHw(kpiHasTarget)));
end
for k = 1:numKpi
    if kpiIsFraction(k)
        fprintf('    %-20s %10.4f  (95%% CI ±%.4f)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), lastHw(k), ternary(kpiHasTarget(k),'',' [no target]'));
    else
        fprintf('    %-20s %10.4f  (95%% CI ±%.2f%%)%s\n', kpiNames{k}, ...
            mean(passKpi(1:numPasses,k),'omitnan'), 100*lastHw(k), ternary(kpiHasTarget(k),'',' [no target]'));
    end
end
fprintf('  Passes: %d | simulated time: %.0f s (%.2f h)\n', ...
    numPasses, simulatedTimeS, simulatedTimeS/3600);
fprintf('  Dataset rows: %d -> %s\n', totalRows, opts.datasetPath);
fprintf('  UE state checks: %d steps, all passed\n', checksRun);

%% ------------------ Πίνακας της διέλευσης αναφοράς ------------------
if ~isempty(fieldnames(refSnapshot))
    fprintf('\n--- Reference pass, instant of peak elevation ---\n');
    printUserTable(size(refSnapshot.user_geo,1), refSnapshot.bestNodeVec, refSnapshot.bestNodeTypeVec, ...
        refSnapshot.bestDistanceVec, refSnapshot.bestPathLossVec, refSnapshot.bestSinrDbVec, ...
        refSnapshot.capacityMbpsVec, refSnapshot.bestElevationDegVec, ...
        refSnapshot.nodePowerWattsVec, refSnapshot.energyPerBitUJVec);
end

%% ------------------ Αποθήκευση και γραφήματα ------------------
% Αποτελέσματα + plotInfo για την plotResults· το dataset μόνο με SHA-256.
plotInfo = struct('numPasses', numPasses, 'ciTolerance', opts.ciTolerance, ...
    'ciToleranceFrac', opts.ciToleranceFrac, 'kpiNames', {kpiNames}, ...
    'kpiIsFraction', kpiIsFraction, 'kpiHasTarget', kpiHasTarget, ...
    'mlbMaxOffsetDb', opts.mlbMaxOffsetDb, 'pingPongMtsS', opts.pingPongMtsS, ...
    'rlfQoutDb', opts.rlfQoutDb, 'rlfQinDb', opts.rlfQinDb, 'rlfT310Ms', opts.rlfT310Ms, ...
    'hysteresisDb', opts.hysteresisDb, 'timeToTriggerMs', opts.timeToTriggerMs, ...
    'gate', [], 'mlReportIntervalMs', opts.mlReportIntervalMs, ...
    'satMinElevationDeg', satParameters.MinElevationDeg, 'refSnapshot', refSnapshot, ...
    'cellGeo', cell_geo, 'numCells', numCells, 'bsIsdM', opts.bsIsdM, ...
    'cellBoresightDeg', simParameters.CellBoresightDeg);
if hasMlGate
    plotInfo.gate = struct('NumTrees', mlGateInfo.NumTrees, 'Threshold', mlGateInfo.Threshold);
end
save(resultsMat, 'refPassTable', 'convTable', 'policyTable', 'pairedTable', 'summary', 'plotInfo');
fprintf('\n');
outFiles = [{resultsMat}, plotResults(resultsMat, opts.outputDir)];

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
runParams.durationRule     = 'fixed number of passes (numPasses); CI reported only';
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
    'source', 'TR 36.839 V11.1.0 Table 5.2.1.3.1; procedure TS 38.331 §5.3.10');
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
runParams.rngScheme        = 'rng(rngSeed+passIdx, twister) at the start of each pass; passes run in parallel batches and are written in order';
runParams.numWorkers       = numWorkers;
runParams.checksRun        = checksRun;
runParams.userArea         = 'uniform over the cells of the sites (TR 38.901 Table 7.2-1)';
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
    'timeStep',   'L1->L3 measurement period 200 ms (TR 36.839 Table 5.3.2.1, TS 38.133 Table 9.2.5.2-1)', ...
    'duration',   'ITU-R M.2412-0 §7.1 (sufficient number of drops, confidence intervals reported); ITU-R M.2514-0 §8.2.4', ...
    'passWindow', 'TR 38.821 Table 4.2-3 NOTE 1 (visibility time of the satellite)', ...
    'ueMobility', 'ITU-R M.2412-0 §8.4 Table 5 b)/c) (3 km/h, uniformly random direction)', ...
    'validity',   'TR 38.901 §7.4.1 Table 7.4.1-1 (10 m - 5 km, h_UT 1.5-22.5 m)', ...
    'hysteresis', 'TS 38.331 §5.5.4.4 (Event A3) + §6.3.2 (Hysteresis, TimeToTrigger)', ...
    'pingPong',   'TR 36.839 Definition 5 (return to the source cell with time-of-stay < MTS = 1 s)', ...
    'measFilter', 'TS 38.331 §5.5.3.2 (Layer 3 filtering, F_n = (1-a)F_{n-1} + a M_n; NOTE 1 k=0; NOTE 2 logarithmic domain)', ...
    'measPeriod', 'TS 38.133 Table 9.2.5.2-1 (FR1 without DRX: >= 200 ms)', ...
    'measAvg',    'Clarke 1968, power autocorrelation J0(2*pi*fD*tau)^2 -> equivalent number of independent samples per step', ...
    'bsLayout',   'TR 38.901 §7.8 Table 7.8-1 (UMa: hexagonal layout, ISD 500 m, h_BS 25 m)', ...
    'satReuse',   'TR 38.821 §6.1.1.1 Table 6.1.1.1-5 (FRF Option 1 = 1 / Option 2 = 3, bandwidth per beam scaled by FRF) + Figure 6.1.1.1-2 (co-channel beams = interference)', ...
    'satSpacing', 'design choice: the numSats satellites are exactly those of the plane visible within the window', ...
    'mlGate',     'TS 38.300 §9.2.3.2.1 step 3 (HO decision at the source gNB) + TR 37.817 §5.3 (AI/ML mobility optimisation) + Aggarwal et al. ICSigSys 2026 (prediction of unstable HOs); ping-pong label TR 36.839 Definition 5', ...
    'mlReport',   'TS 38.331 ReportConfigNR reportInterval (ms120, ms240, ...): the report repeats while the entering condition holds');

saveRunVersion('runSimulation', runParams, outFiles, opts.label);

if isfile(checkpointPath)
    delete(checkpointPath);
end

end

% =====================================================================
% Τοπικές συναρτήσεις
% =====================================================================

function nEff = effectiveSamplesClarke(fDopplerHz, tStepS)
% Ισοδύναμο πλήθος ανεξάρτητων δειγμάτων διαλείψεων σε T (Clarke 1968).
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
% Ένα δείγμα ανά εντολή HO, με χαρακτηριστικά τη στιγμή της απόφασης.
% Outcome (TR 36.839 Ορ. 3): "HO", "HOF_T310", "HOF_TARGET".
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
                    ['Pass %d, step %d, UE %d (%s): the features seen by the ML gate ' ...
                     'differ from the logged ones (max difference %g).'], ...
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
% Τι ακολούθησε κάθε HO (μόνο για την ετικέτα): StayTo_s, NextNode ("None" =
% outage, "" = τέλος), NextCause, Censored. Ping-pong: TR 36.839 Ορ. 5.
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

function [G, info] = loadHandoverGate(path, reportSteps)
% Φόρτωση μοντέλου· έλεγχος ίδιων πιθανοτήτων με το scikit-learn.
S = load(path);
allNames = [handoverFeatures(), {'FromSat','ToSat'}];
featNames = cellstr(S.FeatureNames);
[ok, featIdx] = ismember(featNames, allNames);
if ~all(ok)
    error('runSimulation:GateFeatures', ...
        'Model %s requires features that the simulation does not provide: %s', ...
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
        ['handoverGatePredict differs from scikit-learn on the self-test ' ...
         'samples of %s (max difference %g).'], path, maxErr);
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
fprintf('  ML handover gate: %s\n', path);
fprintf('    %d trees, %d features, threshold %.3f, reportInterval %d steps\n', ...
    numel(M.Roots), numel(featIdx), G.Threshold, reportSteps);
fprintf('    same probabilities as scikit-learn on %d self-test samples (max difference %.1e)\n', ...
    numel(testP), maxErr);
end

function R = runOnePass(passIdx, C)
% Μία διέλευση· ίδιο αποτέλεσμα σειριακά ή σε worker.
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
% Πίνακας μιας διέλευσης από numBuf/strBuf.
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
% Έλεγχος χρηστών: ύψος, συντεταγμένες, πεδίο ισχύος (TR 38.901 Πίν. 7.4.1-1).
if any(user_geo(:,3) ~= ueHeightM)
    error('runSimulation:HeightDrift', ...
        'Pass %d, step %d: UE height %.4f m instead of %.4f m.', ...
        passIdx, stepIdx, max(abs(user_geo(:,3))), ueHeightM);
end
if ~all(isfinite(user_geo(:)))
    error('runSimulation:NonFiniteGeo', ...
        'Pass %d, step %d: non-finite UE coordinate.', passIdx, stepIdx);
end
if any(abs(user_geo(:,1)) > 90) || any(abs(user_geo(:,2)) > 180)
    error('runSimulation:GeoOutOfBounds', ...
        'Pass %d, step %d: UE coordinate out of bounds.', passIdx, stepIdx);
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
            'Pass %d, step %d, UE %d: %.2f m from a site (< %g m).', ...
            passIdx, stepIdx, u, min(d2d), minDistBs);
    end
    if ~any(d2d >= d2dMin & d2d <= d2dMax)
        error('runSimulation:NoValidLink', ...
            ['Pass %d, step %d, UE %d: no terrestrial link within ' ...
             '[%d, %d] m (min distance %.1f m, max %.1f m).'], ...
            passIdx, stepIdx, u, d2dMin, d2dMax, min(d2d), max(d2d));
    end
end
end

function validateOutputs(typeVec, capacityVec, throughputVec, nodeLoadVec, sinrVec, passIdx, stepIdx)
% Συνέπεια εξόδων (outage / συνδεδεμένος).
isOut = typeVec == "Outage";
if any(isOut & (~isnan(capacityVec) | throughputVec ~= 0 | nodeLoadVec ~= 0))
    error('runSimulation:OutageInconsistent', ...
        'Pass %d, step %d: UE without a usable link has non-zero capacity, throughput or load.', ...
        passIdx, stepIdx);
end
if any(~isOut & (~isfinite(capacityVec) | capacityVec < 0))
    error('runSimulation:ServedInconsistent', ...
        'Pass %d, step %d: connected UE without a finite non-negative capacity.', ...
        passIdx, stepIdx);
end
if any(~isOut & sinrVec >= 10*log10(2^0.2344 - 1) & capacityVec <= 0)
    error('runSimulation:ServedInconsistent', ...
        'Pass %d, step %d: UE above MCS 0 with zero capacity.', ...
        passIdx, stepIdx);
end
if any(~isOut & nodeLoadVec < 1)
    error('runSimulation:LoadInconsistent', ...
        'Pass %d, step %d: connected UE on a node with load < 1.', ...
        passIdx, stepIdx);
end
if any(~isOut & isnan(sinrVec) & capacityVec ~= 0)
    error('runSimulation:SinrInconsistent', ...
        'Pass %d, step %d: UE without a modelled link has non-zero capacity.', ...
        passIdx, stepIdx);
end
end

function [lat, lon, walkAz] = randUserAndHeading(centerGeo, bs_geo, ghost_geo, img_geo, wgs84, rDiskM, dMinM, walkM)
% Ομοιόμορφη θέση και κατεύθυνση· απόρριψη αν η διαδρομή πλησιάζει σταθμό < dMinM.
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
% Θέση σε κυκλική τροχιά -> [lat lon alt] (ECEF).
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
