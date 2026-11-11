function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
    newChannelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
    bsUnavailReasonVec, satUnavailReasonVec, policyKpis, newDecision] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, ...
                     prevChannelState, prevDecision)
% Ένα βήμα της προσομοίωσης. Για κάθε χρήστη υπολογίζει το SINR κάθε
% υποψήφιου κόμβου (κελιά BS και δορυφόροι), εφαρμόζει τον κανόνα επιλογής
% (A3, RLF) ανά πολιτική και τη χωρητικότητα/ενέργεια μετά την κατανομή
% εύρους ζώνης. sat_geo: μία γραμμή [lat lon alt] ανά δορυφόρο.
%
% Όλες οι πολιτικές αποτιμώνται στην ίδια πραγματοποίηση καναλιού, οπότε η
% σύγκρισή τους (policyKpis) είναι κατά ζεύγη. Τα prevChannelState και
% prevDecision μεταφέρουν το κανάλι και την απόφαση από το προηγούμενο βήμα
% της διέλευσης· είναι κενά στην αρχή κάθε διέλευσης και στα scripts
% επαλήθευσης.
if nargin < 7
    prevChannelState = [];
end
if nargin < 8
    prevDecision = [];
end

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);   % κελιά (TRxP): κόμβοι "BS1".."BSnumBs"
% Τομείς (TR 38.901 Πίν. 7.8-1): CellSite = site κάθε κελιού, CellBoresightDeg
% = κατεύθυνση τομέα (από ανατολή, αριστερόστροφα). Χωρίς αυτά, ένα κελί ανά
% σταθμό με σταθερό κέρδος (simParameters.EIRP).
if isfield(simParameters, 'CellSite') && ~isempty(simParameters.CellSite)
    cellSite         = simParameters.CellSite(:);
    cellBoresightDeg = simParameters.CellBoresightDeg(:);
else
    cellSite         = (1:numBs)';
    cellBoresightDeg = [];
end
sectorized = ~isempty(cellBoresightDeg);
numSites   = max(cellSite);
siteGeo    = zeros(numSites, 3);
for k = 1:numSites
    siteGeo(k,:) = bs_geo(find(cellSite == k, 1), :);
end
numSats  = size(sat_geo,1);   % κόμβοι "SAT-k"

% Επαναχρησιμοποίηση συχνότητας (TR 38.821 Πίν. 6.1.1.1-5): ο δορυφόρος k
% εκπέμπει στη ζώνη mod(k-1, FRF)· το satParameters.Bandwidth είναι ήδη ανά δέσμη.
if isfield(satParameters, 'FreqReuse') && ~isempty(satParameters.FreqReuse)
    satFreqReuse = satParameters.FreqReuse;
else
    satFreqReuse = 1;
end
satSubBand = mod((0:numSats-1)', satFreqReuse);

%% ------------------ Bandwidth ------------------
BW_bs = simParameters.Carrier.NSizeGrid * 12 * ...
        simParameters.Carrier.SubcarrierSpacing * 1e3;   % Hz

%% ------------------ Noise power ------------------
kBoltz = physconst('Boltzmann');
NF = 10^(simParameters.RxNoiseFigure/10);
Teq = simParameters.RxAntTemperature + 290*(NF-1);
% kTB σε dBW
noisePowerBS_dBW  = 10*log10(kBoltz * Teq * BW_bs);
noisePowerSAT_dBW = 10*log10(kBoltz * Teq * satParameters.Bandwidth);

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SINR ------------------
% Shannon-ισοδύναμο του MCS 0 (TS 38.214 Πίν. 5.1.3.1-1)· κάτω από αυτό ο
% χρήστης είναι εκτός κάλυψης.
minSpectralEfficiency = 0.2344;                        % bit/s/Hz (MCS 0)
minUsableSinrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Στόχος υπηρεσίας ------------------
% 5ο εκατοστημόριο φασματικής απόδοσης χρήστη (ITU-R M.2410 / TR 37.910
% Πίν. 5.4.1.1.1-1), επί του συνολικού εύρους: SE_ζεύξης >= target * L.
targetNormalizedSe = 0.3;   % bit/s/Hz

%% ------------------ Ισχύς ανά αλυσίδα πομποδέκτη ------------------
% EARTH (Auer et al. 2011, εξ. 1, Πίν. 2): P_out ανά αλυσίδα, P_max = 20 W.
pOutTotalW    = 10^((simParameters.TxPower - 30)/10);
pOutPerChainW = pOutTotalW / simParameters.Power.NumTrx;
earthPmaxW    = 20;
if pOutPerChainW > earthPmaxW
    warning('simulateScenario:EarthOutOfRange', ...
        ['Ισχύς ανά αλυσίδα %.1f W > P_max = %.0f W του μοντέλου EARTH ' ...
         '(Auer et al. 2011, Πίν. 2). Αύξησε το simParameters.Power.NumTrx.'], ...
        pOutPerChainW, earthPmaxW);
end

%% ------------------ Διαλείψεις ------------------
% Δορυφόρος: Shadowed Rician (Abdi et al. 2003), παράμετροι από την ανύψωση
% (εξ. 19)· η σκίαση περιέχεται στο μοντέλο. Συσχέτιση σκίασης 37 m (TR 38.811
% Πίν. 6.7.2-3a, Urban, S-band, LOS).
satFadeElevRangeDeg = [20 80];   % πεδίο ισχύος της εξ. (19)
satShadowDcorrM     = 37;
terrKdBLos          = 9;         % Rician K επίγειο LOS, TR 38.901 Πίν. 7.5-6

%% ------------------ A3: υστέρηση και χρόνος επιβεβαίωσης ------------------
% TS 38.331 §5.5.4.4 / §6.3.2. Χωρίς simParameters.Mobility: 0 dB, 0 ms.
if isfield(simParameters, 'Mobility')
    hysDb    = simParameters.Mobility.HysteresisDb;
    tttMs    = simParameters.Mobility.TimeToTriggerMs;
    dtSecHys = simParameters.Mobility.DtSeconds;
else
    hysDb = 0; tttMs = 0; dtSecHys = 1;
end
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'PingPongMtsS')
    pingPongMtsS = simParameters.Mobility.PingPongMtsS;
else
    pingPongMtsS = 1;
end
% Ping-pong: παραμονή μικρότερη από το MTS (TR 36.839 Ορισμός 5), σε βήματα
% (Δt = 0.2 s: έως 4 βήματα = 0.8 s).
ppMaxStaySteps = max(1, ceil(pingPongMtsS/dtSecHys - 1e-9) - 1);
% TTT σε βήματα, προς τα πάνω: η συνθήκη πρέπει να ισχύει όλη τη διάρκεια.
nTtt = max(1, ceil(tttMs / (dtSecHys*1000)));

%% ------------------ Παρακολούθηση ζεύξης και RLF ------------------
% TS 38.331 §5.3.10 με τις τιμές του TR 36.839 Πίν. 5.2.1.3.1: N310 ενδείξεις
% SINR < Qout ξεκινούν τον T310, N311 ενδείξεις SINR > Qin τον σταματούν· αν
% λήξει, RLF και επανεγκατάσταση (§5.3.7). Μία ένδειξη ανά βήμα.
rlfP = struct('QoutDb', -8, 'QinDb', -6, 'T310Ms', 1000, 'N310', 1, 'N311', 1);
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'Rlf')
    rlfP = simParameters.Mobility.Rlf;
end
mob = struct('MinUsableDb', minUsableSinrDb, 'HysDb', hysDb, 'NTtt', nTtt, ...
    'QoutDb', rlfP.QoutDb, 'QinDb', rlfP.QinDb, 'N310', rlfP.N310, 'N311', rlfP.N311, ...
    'NT310', max(1, ceil(rlfP.T310Ms / (dtSecHys*1000))), 'NReport', 0);

%% ------------------ Εξισορρόπηση φόρτου (MLB) ------------------
% TS 38.300 §15.5.1: μετατόπιση του A3 (Ocn) ανάλογη της διαφοράς φόρτου,
% κατά Nasri & Altman (arXiv:1307.1212, εξ. 4-5):
%     offset = MlbMaxOffsetDb * (χ_εξυπηρετούντος - χ_στόχου),  διαφορά ∈ [-1, 1]
% στρογγυλεμένη στις τιμές του Q-OffsetRange (TS 38.331).
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'MlbMaxOffsetDb')
    mlbMaxOffsetDb = simParameters.Mobility.MlbMaxOffsetDb;
else
    mlbMaxOffsetDb = 6;   % = HMmax - f(0) του Nasri & Altman
end

%% ------------------ Απόφαση handover με ML ------------------
% Ο σταθμός αποφασίζει το handover μετά την αναφορά μέτρησης (TS 38.300
% §9.2.3.2.1, βήμα 3· TR 37.817 §5.3). Το μοντέλο εκτιμά P(ping-pong) (TR 36.839
% Ορισμός 5)· αν ξεπερνά το κατώφλι, ο σταθμός περιμένει την επόμενη αναφορά
% (reportInterval, TS 38.331 ReportConfigNR) και ξανακρίνει.
%
% simParameters.Mobility.MlGate (κενό = χωρίς τις πολιτικές *ML):
%   Model        δέντρα του μοντέλου (handoverGatePredict)
%   Threshold    κατώφλι πιθανότητας
%   FeatIdx      στήλες του [handoverFeatures(), FromSat, ToSat] που δέχεται
%   ReportSteps  reportInterval σε βήματα
mlGate = [];
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'MlGate')
    mlGate = simParameters.Mobility.MlGate;
end
hasMl = ~isempty(mlGate);
% Ιστορικό της φιλτραρισμένης μέτρησης (t, t-1 s, t-2 s) για τα χαρακτηριστικά.
histW   = 2*round(1/dtSecHys) + 1;
histLag = (histW - 1)/2;
if hasMl
    mob.NReport = mlGate.ReportSteps;
end

%% ------------------ Μέτρηση: μέση τιμή και φίλτρο L3 ------------------
% Ο A3 συγκρίνει τη φιλτραρισμένη μέτρηση F_n = (1-a) F_{n-1} + a M_n (TS 38.331
% §5.5.3.2)· το M_n είναι μέση τιμή περιόδου 200 ms (TS 38.133 Πίν. 9.2.5.2-1).
%   Meas.NAvg    ανεξάρτητες πραγματώσεις γρήγορων διαλείψεων ανά βήμα (Clarke)
%   Meas.L3Alpha συντελεστής a προσαρμοσμένος στο βήμα· 1 = χωρίς φίλτρο
measNAvgTerr = max(1, round(measParam(simParameters, 'NAvg',    1)));
measNAvgSat  = max(1, round(measParam(satParameters, 'NAvg',    1)));
l3Alpha      = min(max(measParam(simParameters, 'L3Alpha', 1), 0), 1);

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestBsNodeVec       = strings(numUsers,1);   % καλύτερος BS, για τη σύγκριση πολιτικών
% Πολιτικές: οι τρεις βασικές και οι δύο εκδοχές με MLB (ίδιος A3 + μετατόπιση).
policyNamesLocal    = {'Actual','TerrestrialOnly','SatelliteOnly', ...
                       'ActualMLB','TerrestrialOnlyMLB'};
policyBaseLocal     = {'Actual','TerrestrialOnly','SatelliteOnly', ...
                       'Actual','TerrestrialOnly'};
policyMlbLocal      = [false false false true true];
% Με μοντέλο: η μεικτή χωρίς και με MLB, με την απόφαση του σταθμού από το ML.
if hasMl
    policyNamesLocal = [policyNamesLocal, {'ActualML','ActualMLBML'}];
    policyBaseLocal  = [policyBaseLocal,  {'Actual','Actual'}];
    policyMlbLocal   = [policyMlbLocal, false, true];
end
% Επιπλέον αποφάσεις του σταθμού για σύγκριση, όλες με MLB
% (simParameters.Mobility.ExtraGates: Name, Kind, Gate):
%   Kind "history"  κανόνας UE History (TS 38.300 §15.5.2.4): αναμονή αν ο στόχος
%                   είναι ο προηγούμενος κόμβος και η παραμονή < MTS.
%   Kind "ml"       άλλο μοντέλο, π.χ. από άλλον σπόρο εκπαίδευσης.
extraGates = [];
if hasMl && isfield(simParameters.Mobility, 'ExtraGates')
    extraGates = simParameters.Mobility.ExtraGates;
end
policyGateLocal = [cell(1, 5), repmat({mlGate}, 1, numel(policyNamesLocal) - 5)];
policyKindLocal = [repmat("", 1, 5), repmat("ml", 1, numel(policyNamesLocal) - 5)];
for k = 1:numel(extraGates)
    policyNamesLocal{end+1} = extraGates(k).Name; %#ok<AGROW>
    policyBaseLocal{end+1}  = 'Actual'; %#ok<AGROW>
    policyMlbLocal(end+1)   = true; %#ok<AGROW>
    policyGateLocal{end+1}  = extraGates(k).Gate; %#ok<AGROW>
    policyKindLocal(end+1)  = string(extraGates(k).Kind); %#ok<AGROW>
end
% Dual connectivity (TS 37.340, μεταφορά στο ISTN): MCG = επίγειο κελί με την
% κινητικότητα της TerrestrialOnly, SCG = δορυφόρος (scgStep).
%   DcSplit  split bearer: ρυθμός R_MCG + R_SCG
%   DcDup    PDCP duplication (TS 38.323 §5.11): ρυθμός max(R_MCG, R_SCG)
% simParameters.Mobility.Dc: Enabled, SnAddThrDb (κατώφλι του A4).
hasDc = isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'Dc') && ...
    simParameters.Mobility.Dc.Enabled;
if hasDc
    dcAddThrDb = simParameters.Mobility.Dc.SnAddThrDb;
    policyNamesLocal = [policyNamesLocal, {'DcSplit','DcDup'}];
    policyBaseLocal  = [policyBaseLocal,  {'TerrestrialOnly','TerrestrialOnly'}];
    policyMlbLocal   = [policyMlbLocal, false, false];
    policyGateLocal  = [policyGateLocal, {[], []}];
    policyKindLocal  = [policyKindLocal, "dc", "dc"];
end
numPolLocal         = numel(policyNamesLocal);
policyMlLocal       = policyKindLocal == "ml" | policyKindLocal == "history";
policyDcLocal       = policyKindLocal == "dc";
nGateFeat           = numel(handoverFeatures());
% Κατάσταση απόφασης ανά πολιτική και χρήστη:
%   Node/Pending/Count   εξυπηρετών, υποψήφιος του A3, βήματα επιβεβαίωσης
%   Wait                 βήματα ως την επόμενη αναφορά μετά από αναμονή (*ML)
%   Prev/Stay            κόμβος πριν το τελευταίο HO και βήματα παραμονής (ping-pong)
%   Cause/Offset         πώς άλλαξε ο κόμβος ("HO", "RLF", "HOF", "") και η μετατόπιση
%   Attempt/AttemptOffset στόχος και μετατόπιση της εντολής HO του βήματος
%   Oos/Is/T310          μετρητές N310/N311 και βήματα του T310 (-1 = σταματημένος)
%   LoadNode/LoadRho     εικονικός φόρτος ανά κόμβο (αναφορά φόρτου, TS 38.300 §15.5.1)
%   GateFeat             χαρακτηριστικά που είδε το μοντέλο (*ML· NaN = δεν ρωτήθηκε)
newDecision = struct();
for pI = 1:numPolLocal
    newDecision.(policyNamesLocal{pI}) = struct('Node', strings(numUsers,1), ...
        'Pending', strings(numUsers,1), 'Count', zeros(numUsers,1), 'Wait', zeros(numUsers,1), ...
        'Prev', strings(numUsers,1), 'Stay', zeros(numUsers,1), ...
        'Cause', strings(numUsers,1), 'Offset', zeros(numUsers,1), ...
        'Attempt', strings(numUsers,1), 'AttemptOffset', zeros(numUsers,1), ...
        'Oos', zeros(numUsers,1), 'Is', zeros(numUsers,1), 'T310', -ones(numUsers,1), ...
        'LoadNode', strings(0,1), 'LoadRho', zeros(0,1), 'LoadAge', 0);
    if policyMlLocal(pI)
        newDecision.(policyNamesLocal{pI}).GateFeat = nan(numUsers, nGateFeat);
    end
    if policyDcLocal(pI)
        % SCG: PSCell, υποψήφιος και βήματα του A4/A3, RLM του PSCell
        newDecision.(policyNamesLocal{pI}).SnNode    = repmat("None", numUsers, 1);
        newDecision.(policyNamesLocal{pI}).SnPending = strings(numUsers,1);
        newDecision.(policyNamesLocal{pI}).SnCount   = zeros(numUsers,1);
        newDecision.(policyNamesLocal{pI}).SnOos     = zeros(numUsers,1);
        newDecision.(policyNamesLocal{pI}).SnIs      = zeros(numUsers,1);
        newDecision.(policyNamesLocal{pI}).SnT310    = -ones(numUsers,1);
    end
end
hasPrevDecision = ~isempty(prevDecision);

% Μετρητές του βήματος ανά πολιτική: επιτυχή HO ανά τύπο, ping-pong, RLF ανά
% τύπο κόμβου (TR 36.839 Ορισμός 2), αποτυχίες HO (Ορισμός 3), αναμονές ML.
hoCount = struct();
for pol = policyNamesLocal
    hoCount.(pol{1}) = struct('TerrTerr', 0, 'TerrSat', 0, 'SatSat', 0, 'PingPong', 0, ...
        'RlfTerr', 0, 'RlfSat', 0, 'HoFail', 0, 'Held', 0, ...
        'SnAdd', 0, 'SnChange', 0, 'ScgFail', 0, 'McgRlfCovered', 0);
end
bestDistanceVec     = nan(numUsers,1);
bestPathLossVec     = nan(numUsers,1);
bestSinrDbVec        = nan(numUsers,1);
capacityMbpsVec     = nan(numUsers,1);
nodePowerWattsVec   = nan(numUsers,1);
energyPerBitUJVec   = nan(numUsers,1);
bestElevationDegVec = nan(numUsers,1);
bestBsSinrDbVec      = nan(numUsers,1);
bestBsDistanceVec   = nan(numUsers,1);
bestBsPathLossVec   = nan(numUsers,1);
serviceStateVec     = strings(numUsers,1);
throughputMbpsVec   = zeros(numUsers,1);
bsUnavailReasonVec  = strings(numUsers,1);
satUnavailReasonVec = strings(numUsers,1);

% Αναθέσεις ανά πολιτική, μία στήλη ανά policyNamesLocal.
polNodeMat = strings(numUsers,numPolLocal);
polTypeMat = strings(numUsers,numPolLocal);
polSinrMat = nan(numUsers,numPolLocal);
% PSCell και SINR του ανά πολιτική DC ("None" = χωρίς SCG).
snNodeMat  = repmat("None", numUsers, numPolLocal);
snSinrMat  = nan(numUsers,numPolLocal);

% Η διάδοση ορίζεται ανά site· τα κελιά ενός site διαφέρουν μόνο στο κέρδος κεραίας.
pathLossMat       = nan(numUsers,numSites);
pLosMat           = nan(numUsers,numSites);
losMat            = false(numUsers,numSites);
losLatentMat      = zeros(numUsers,numSites);
sfMat             = nan(numUsers,numSites);
% Soft LOS (TR 38.901 §7.6.3.3): μία διεργασία σκίασης για κάθε κατάσταση.
sfLosMat          = zeros(numUsers,numSites);
sfNlosMat         = zeros(numUsers,numSites);
losSoftMat        = nan(numUsers,numSites);
useSoftLos = isfield(simParameters, 'SoftLos') && simParameters.SoftLos;
% sqrt(d_LOS/λ) της εξ. 7.6-18, με d_LOS = 50 m (απόσταση συσχέτισης LOS, Πίν. 7.6.3.1-2).
softLosK = sqrt(50 / (physconst('LightSpeed')/simParameters.CarrierFrequency));
snrDbMat          = nan(numUsers,numBs);   % χωρίς παρεμβολή
sinrDbMat         = nan(numUsers,numBs);
cellPathLossMat   = nan(numUsers,numBs);   % = απώλειες του site του κελιού
cellRange3DMat    = nan(numUsers,numBs);
% Καλύτερος δορυφόρος ανά χρήστη και SINR κάθε δορυφόρου.
satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSinrDbVec       = nan(numUsers,1);
satSinrDbMat      = -inf(numUsers,numSats);
satElevMat        = nan(numUsers,numSats);

hasPrevState = ~isempty(prevChannelState) && ...
    isfield(prevChannelState, 'LosLatent') && ...
    isequal(size(prevChannelState.IsLOS), [numUsers, numSites]) && ...
    (~useSoftLos || isfield(prevChannelState, 'SfLos_dB'));
% Κρυφή γκαουσιανή μεταβλητή της δορυφορικής σκίασης ανά (χρήστη, δορυφόρο).
hasPrevSatShadow = hasPrevState && isfield(prevChannelState, 'SatShadowLatent') && ...
    isequal(size(prevChannelState.SatShadowLatent), [numUsers, numSats]);
satShadowLatentMat = zeros(numUsers,numSats);

% Κατάσταση του φίλτρου L3 ανά υποψήφιο (NaN: καμία προηγούμενη μέτρηση,
% F_0 = M_1, §5.5.3.2). Κοινή για όλες τις πολιτικές.
hasPrevFilt = hasPrevState && isfield(prevChannelState, 'FiltSinrBs_dB') && ...
    isequal(size(prevChannelState.FiltSinrBs_dB), [numUsers, numBs]) && ...
    isequal(size(prevChannelState.FiltSinrSat_dB), [numUsers, numSats]);
if hasPrevFilt
    prevFiltBs  = prevChannelState.FiltSinrBs_dB;
    prevFiltSat = prevChannelState.FiltSinrSat_dB;
else
    prevFiltBs  = nan(numUsers, numBs);
    prevFiltSat = nan(numUsers, numSats);
end
filtSinrBsMat  = nan(numUsers, numBs);
filtSinrSatMat = nan(numUsers, numSats);

% Κυκλικό ιστορικό histW βημάτων της φιλτραρισμένης μέτρησης (μόνο *ML)·
% histStep = αύξων αριθμός του βήματος στη διέλευση.
if hasMl
    if hasPrevFilt && isfield(prevChannelState, 'HistStep')
        histBs   = prevChannelState.HistFiltBs;
        histSat  = prevChannelState.HistFiltSat;
        histStep = prevChannelState.HistStep + 1;
    else
        histBs   = nan(histW, numUsers, numBs);
        histSat  = nan(histW, numUsers, numSats);
        histStep = 1;
    end
end

lambdaSat    = physconst('LightSpeed') / satParameters.CarrierFrequency;
noiseSatLinW = 10^(noisePowerSAT_dBW/10);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SINR) ------------------
useWrapAround = isfield(simParameters, 'BsImageGeo') && ~isempty(simParameters.BsImageGeo);
if useWrapAround
    bsImageGeo = simParameters.BsImageGeo;
end

% Ντετερμινιστική γεωμετρία και απώλειες όλων των ζεύξεων, διανυσματικά· στον
% βρόχο των χρηστών μένουν μόνο οι κληρώσεις.
uLat  = user_geo(:,1);
uLon  = user_geo(:,2);
hUtV  = user_geo(:,3);
% Wrap-around: BsImageGeo (numSites x K x 2) = [lat lon] του site και των
% αντιγράφων του· κάθε ζεύξη μετριέται από το πλησιέστερο αντίγραφο.
if useWrapAround
    imgLatRow = reshape(bsImageGeo(:,:,1), 1, []);
    imgLonRow = reshape(bsImageGeo(:,:,2), 1, []);
    nImgPts   = numel(imgLatRow);
    dImgAll = distance(repmat(uLat, 1, nImgPts), repmat(uLon, 1, nImgPts), ...
        repmat(imgLatRow, numUsers, 1), repmat(imgLonRow, numUsers, 1), wgs84);
    [gdAll, nearImgAll] = min(reshape(dImgAll, numUsers, numSites, []), [], 3);
else
    gdAll = distance(repmat(siteGeo(:,1)', numUsers, 1), repmat(siteGeo(:,2)', numUsers, 1), ...
        repmat(uLat, 1, numSites), repmat(uLon, 1, numSites), wgs84);
end
% Γεωδαιτική οριζόντια απόσταση και φυσικά ύψη κεραιών (όχι ENU, που αποκλίνει
% λόγω καμπυλότητας).
groundDistanceMat = gdAll;
range3DMat        = hypot(gdAll, siteGeo(:,3)' - hUtV);
% Πεδίο ισχύος UMa/UMi (TR 38.901 Πίν. 7.4.1-1)· εκτός ορίων το SINR μένει NaN.
validLinkMat = isValidTerrestrialLink(gdAll, repmat(hUtV, 1, numSites), simParameters.CarrierFrequency);
% Πιθανότητα LOS κάθε ζεύξης (TR 38.901 §7.4.2).
pLosAll = losProbability38901(gdAll, repmat(hUtV, 1, numSites), simParameters.PathLoss.Scenario);
% Απώλειες LOS και NLOS κάθε έγκυρης ζεύξης.
plLosMat  = nan(numUsers, numSites);  sfSigmaLosMat  = nan(numUsers, numSites);
plNlosMat = nan(numUsers, numSites);  sfSigmaNlosMat = nan(numUsers, numSites);
for b = 1:numSites
    iv = find(validLinkMat(:,b));
    if isempty(iv)
        continue;
    end
    txPosition = [0; 0; siteGeo(b,3)];
    rxPosition = [gdAll(iv,b)'; zeros(1, numel(iv)); hUtV(iv)'];
    [pl, sg] = nrPathLoss(simParameters.PathLoss, simParameters.CarrierFrequency, ...
        true(1, numel(iv)), txPosition, rxPosition);
    plLosMat(iv,b) = pl';  sfSigmaLosMat(iv,b) = sg';
    [pl, sg] = nrPathLoss(simParameters.PathLoss, simParameters.CarrierFrequency, ...
        false(1, numel(iv)), txPosition, rxPosition);
    plNlosMat(iv,b) = pl'; sfSigmaNlosMat(iv,b) = sg';
end
% Αζιμούθιο site -> χρήστη (από βορρά, δεξιόστροφα) για την κεραία τομέα.
if sectorized
    if useWrapAround
        siteIdxM = repmat(1:numSites, numUsers, 1);
        imgLatM = bsImageGeo(sub2ind(size(bsImageGeo), siteIdxM, nearImgAll, ones(numUsers, numSites)));
        imgLonM = bsImageGeo(sub2ind(size(bsImageGeo), siteIdxM, nearImgAll, 2*ones(numUsers, numSites)));
    else
        imgLatM = repmat(siteGeo(:,1)', numUsers, 1);
        imgLonM = repmat(siteGeo(:,2)', numUsers, 1);
    end
    [~, azNorthAll] = distance(imgLatM, imgLonM, repmat(uLat, 1, numSites), ...
        repmat(uLon, 1, numSites), wgs84);
end
% Μετακίνηση κάθε χρήστη από το προηγούμενο βήμα (συσχέτιση της σκίασης).
if hasPrevState
    userMoveAll = distance(prevChannelState.UserGeo(:,1), prevChannelState.UserGeo(:,2), ...
        uLat, uLon, wgs84);
else
    userMoveAll = inf(numUsers, 1);   % καμία προηγούμενη κατάσταση -> ανεξάρτητο δείγμα
end
% Ανύψωση και απόσταση κάθε δορυφόρου από κάθε χρήστη.
satElevAll  = nan(numUsers, numSats);
satRangeAll = nan(numUsers, numSats);
for k = 1:numSats
    [~, satElevAll(:,k), satRangeAll(:,k)] = geodetic2aer( ...
        repmat(sat_geo(k,1), numUsers, 1), repmat(sat_geo(k,2), numUsers, 1), ...
        repmat(sat_geo(k,3), numUsers, 1), uLat, uLon, hUtV, wgs84);
end
cellNames = "BS" + string((1:numBs)');
satNames  = "SAT-" + string((1:numSats)');
% MLB: φόρτος κάθε κόμβου από την αναφορά του προηγούμενου βήματος (κόμβος
% χωρίς χρήστες: 0) και ο δείκτης του εξυπηρετούντος κάθε χρήστη.
nodeLoadAll = struct();
prevNodeIdx = struct();
if hasPrevDecision
    allNodeNames = [cellNames; satNames];
    for pI = find(policyMlbLocal)
        polName = policyNamesLocal{pI};
        ldAll = zeros(numBs + numSats, 1);
        [tf, loc] = ismember(prevDecision.(polName).LoadNode, allNodeNames);
        ldAll(loc(tf)) = prevDecision.(polName).LoadRho(tf);
        nodeLoadAll.(polName) = ldAll;
        [~, prevNodeIdx.(polName)] = ismember(prevDecision.(polName).Node, allNodeNames);
    end
end

for u = 1:numUsers
    userBestSinr      = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

    userMoveDistance = userMoveAll(u);

    %% ===== Επίγειοι υποψήφιοι =====
    % Κληρώσεις ανά ζεύξη (LOS, σκίαση, γρήγορες διαλείψεις), σειρά χρήστης -> site.
    for b = 1:numSites
        if ~validLinkMat(u,b)
            sfMat(u,b) = 0;
            continue;
        end

        % Κατάσταση LOS με χωρική συνέπεια (TR 38.901 §7.4.2, §7.6.3.3).
        pLos = pLosAll(u,b);
        if hasPrevState
            prevIsLos  = prevChannelState.IsLOS(u,b);
            prevSF     = prevChannelState.ShadowFading_dB(u,b);
            prevLatent = prevChannelState.LosLatent(u,b);
        else
            prevIsLos  = false;
            prevSF     = 0;
            prevLatent = 0;
        end
        [isLos, losLatent] = spatiallyConsistentLos(pLos, userMoveDistance, ...
            hasPrevState, prevLatent);
        pLosMat(u,b)      = pLos;
        losMat(u,b)       = isLos;
        losLatentMat(u,b) = losLatent;

        if useSoftLos
            % Soft LOS (TR 38.901 §7.6.3.3, εξ. 7.6-18/19):
            %   LOS_soft = 1/2 + (1/π) arctan(sqrt(d_LOS/λ) (G + F(d))),
            %   F(d) = sqrt(2) erf^-1(2 P_LOS(d) - 1),  G = -latent,
            %   P = LOS_soft P_LOS + (1 - LOS_soft) P_NLOS  (ισχύς, ανεξάρτητες φάσεις).
            losSoft = 0.5 + atan(softLosK*(-losLatent + sqrt(2)*erfinv(2*pLos - 1)))/pi;
            losSoftMat(u,b) = losSoft;
            if hasPrevState
                rhoL = shadowFadingCorrelation(userMoveDistance, simParameters.PathLoss.Scenario, true);
                rhoN = shadowFadingCorrelation(userMoveDistance, simParameters.PathLoss.Scenario, false);
                sfL = rhoL*prevChannelState.SfLos_dB(u,b)  + sqrt(1 - rhoL^2)*sfSigmaLosMat(u,b)*randn();
                sfN = rhoN*prevChannelState.SfNlos_dB(u,b) + sqrt(1 - rhoN^2)*sfSigmaNlosMat(u,b)*randn();
            else
                sfL = sfSigmaLosMat(u,b)*randn();
                sfN = sfSigmaNlosMat(u,b)*randn();
            end
            sfLosMat(u,b) = sfL; sfNlosMat(u,b) = sfN;
            sfMat(u,b) = isLos*sfL + (~isLos)*sfN;   % διαγνωστικό
            plL = plLosMat(u,b)  + sfL - smallScaleFadingDb(true,  terrKdBLos, measNAvgTerr);
            plN = plNlosMat(u,b) + sfN - smallScaleFadingDb(false, terrKdBLos, measNAvgTerr);
            pathLossMat(u,b) = -10*log10(losSoft*10^(-plL/10) + (1 - losSoft)*10^(-plN/10));
            continue;
        end

        rho = shadowFadingCorrelation(userMoveDistance, ...
            simParameters.PathLoss.Scenario, isLos);
        useCorrelatedSF = hasPrevState && (isLos == prevIsLos);

        if isLos
            pathLoss = plLosMat(u,b);  sigmaSF = sfSigmaLosMat(u,b);
        else
            pathLoss = plNlosMat(u,b); sigmaSF = sfSigmaNlosMat(u,b);
        end

        % Σκίαση log-normal (TR 38.901 §7.4.1), AR(1) με το προηγούμενο βήμα (§7.4.4).
        if useCorrelatedSF
            sfSample = rho*prevSF + sqrt(1 - rho^2) * sigmaSF * randn();
        else
            sfSample = sigmaSF * randn();
        end
        sfMat(u,b) = sfSample;
        pathLoss = pathLoss + sfSample;

        % Γρήγορες διαλείψεις: Rician (LOS) ή Rayleigh (NLOS), μέση τιμή NAvg
        % πραγματώσεων μέσα στο βήμα (TS 38.133 Πίν. 9.2.5.2-1).
        pathLoss = pathLoss - smallScaleFadingDb(isLos, terrKdBLos, measNAvgTerr);
        pathLossMat(u,b) = pathLoss;
    end

    %% ===== Λαμβανόμενη ισχύς ανά κελί (dBW) =====
    % Με τομείς: κέρδος κεραίας στη γεωμετρική κατεύθυνση του χρήστη
    % (bsSectorGainDb, TR 38.901 Πίν. 7.8-1). Χωρίς τομείς: σταθερό EIRP.
    cellPathLossMat(u,:) = pathLossMat(u, cellSite);
    cellRange3DMat(u,:)  = range3DMat(u, cellSite);
    if sectorized
        phiSiteDeg   = 90 - azNorthAll(u,:)';   % GCS: από ανατολή, αριστερόστροφα
        % Ζενιθιακή γωνία στο site (90° = οριζόντια, > 90° προς τα κάτω).
        thetaSiteDeg = 90 + atand((siteGeo(:,3) - user_geo(u,3)) ./ groundDistanceMat(u,:)');
        phiCell = mod(phiSiteDeg(cellSite) - cellBoresightDeg + 180, 360) - 180;
        gainDb  = bsSectorGainDb(thetaSiteDeg(cellSite), phiCell, simParameters.BsAntenna);
        rxDbW   = (simParameters.TxPower - 30) + gainDb' - cellPathLossMat(u,:);
    else
        rxDbW   = (simParameters.EIRP - 30) - cellPathLossMat(u,:);
    end
    snrDbMat(u,:) = rxDbW - noisePowerBS_dBW;   % χωρίς παρεμβολή

    %% ===== SINR ανά επίγειο υποψήφιο =====
    % SINR = C / (I + N) (TR 38.821 Πίν. 6.1.1.2-1). Όλοι οι σταθμοί εκπέμπουν
    % συνεχώς στο ίδιο φάσμα: reuse 1, full buffer (ITU-R M.2412-0 §8.4).
    % Σταθμοί εκτός πεδίου ισχύος του μοντέλου δεν μετρούν ως παρεμβολή.
    noiseLinW   = 10^(noisePowerBS_dBW/10);
    rxPowLinW   = 10.^(rxDbW/10);
    rxPowLinW(~isfinite(rxPowLinW)) = 0;
    totalRxLinW = sum(rxPowLinW);

    vb = isfinite(cellPathLossMat(u,:));
    sinrRow = 10*log10(rxPowLinW(vb) ./ (noiseLinW + (totalRxLinW - rxPowLinW(vb))));
    sinrDbMat(u,vb) = sinrRow;
    % Φίλτρο L3.
    fPrev = prevFiltBs(u,vb);
    fNew  = sinrRow;
    hasF  = isfinite(fPrev);
    fNew(hasF) = (1 - l3Alpha)*fPrev(hasF) + l3Alpha*sinrRow(hasF);
    filtSinrBsMat(u,vb) = fNew;

    % Καλύτερο κελί: ο πρώτος με το μέγιστο SINR.
    if any(vb)
        vbIdx = find(vb);
        [mxSinr, iMx] = max(sinrRow);
        if mxSinr > userBestSinr
            b = vbIdx(iMx);
            userBestSinr      = mxSinr;
            userBestNode      = cellNames(b);
            userBestType      = "Terrestrial";
            userBestDistance  = cellRange3DMat(u,b);
            userBestPathLoss  = cellPathLossMat(u,b);
            userBestElevation = NaN;
        end
    end

    % Καλύτερος BS, πριν τη σύγκριση με τους δορυφόρους.
    bestBsNodeVec(u) = userBestNode;
    if isfinite(userBestSinr)
        bestBsSinrDbVec(u) = userBestSinr;
        if userBestSinr < minUsableSinrDb
            bsUnavailReasonVec(u) = "BelowSinrFloor";
        end
    else
        bsUnavailReasonVec(u) = "OutOfModelRange";
    end
    bestBsDistanceVec(u) = userBestDistance;
    bestBsPathLossVec(u) = userBestPathLoss;

    %% ===== Δορυφορικοί υποψήφιοι =====
    % Πρώτα οι απώλειες όλων: το SINR του καθενός εξαρτάται από τους υπόλοιπους.
    satPlRow    = inf(1,numSats);
    satElevRow  = nan(1,numSats);
    satRangeRow = nan(1,numSats);
    for k = 1:numSats
        % Σκίαση AR(1) ανά δορυφόρο, ανεξάρτητη μεταξύ δορυφόρων· εξελίσσεται
        % σε κάθε βήμα, και κάτω από τη μάσκα.
        if hasPrevSatShadow
            rhoSat = exp(-userMoveDistance / satShadowDcorrM);
            satShadowLatentMat(u,k) = rhoSat*prevChannelState.SatShadowLatent(u,k) + ...
                sqrt(1 - rhoSat^2)*randn();
        else
            satShadowLatentMat(u,k) = randn();
        end

        elevSat       = satElevAll(u,k);
        slantRangeSat = satRangeAll(u,k);
        satElevRow(k)  = elevSat;
        satRangeRow(k) = slantRangeSat;
        satElevMat(u,k) = elevSat;

        if elevSat < satParameters.MinElevationDeg
            continue;
        end
        satPathLoss = fspl(slantRangeSat, lambdaSat);

        % Απόσβεση αερίων (TR 38.811 §6.6.4)· βροχή και σπινθηρισμός δεν μοντελοποιούνται.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        % Shadowed Rician: η κρυφή γκαουσιανή της σκίασης μετατρέπεται στην
        % κατανομή Gamma του Abdi μέσω της αθροιστικής της, οπότε η κατανομή
        % μένει ακριβής και προστίθεται μόνο η χρονική συσχέτιση.
        if isfield(satParameters, 'ShadowingState') && ~isempty(satParameters.ShadowingState)
            [b0, mNak, omega] = shadowedRicianStateParams(satParameters.ShadowingState);
        else
            elevClamped = min(max(elevSat, satFadeElevRangeDeg(1)), satFadeElevRangeDeg(2));
            [b0, mNak, omega] = shadowedRicianElevParams(elevClamped);
        end
        uShadow = 0.5*erfc(-satShadowLatentMat(u,k)/sqrt(2));   % γκαουσιανή -> ομοιόμορφη
        satPlRow(k) = satPathLoss - shadowedRicianFadingDb(b0, mNak, omega, measNAvgSat, uShadow);
    end

    % SINR ανά δορυφόρο: παρεμβολή από τις άλλες ορατές δέσμες της ίδιας ζώνης
    % (TR 38.821 §6.1.1.1)· όχι από πλευρικούς λοβούς (κάτω φράγμα), όχι από
    % το επίγειο σκέλος (2.0 έναντι 3.5 GHz).
    satRxLinW = 10.^(((satParameters.EIRP - 30) - satPlRow)/10);   % PL = Inf -> 0
    for k = 1:numSats
        if ~isfinite(satPlRow(k))
            continue;
        end
        coChannel = (satSubBand == satSubBand(k))';
        coChannel(k) = false;
        satSinrDbMat(u,k) = 10*log10(satRxLinW(k) / (noiseSatLinW + sum(satRxLinW(coChannel))));
        filtSinrSatMat(u,k) = l3Filter(satSinrDbMat(u,k), prevFiltSat(u,k), l3Alpha);
    end

    % Καλύτερος δορυφόρος: ο ορατός με το μέγιστο SINR, αλλιώς με τη μέγιστη ανύψωση.
    if any(isfinite(satPlRow))
        [~, kBest] = max(satSinrDbMat(u,:));
    else
        [~, kBest] = max(satElevRow);
        satUnavailReasonVec(u) = "NotVisible";
    end
    satSlantRangeVec(u) = satRangeRow(kBest);
    satElevationVec(u)  = satElevRow(kBest);
    satPathLossVec(u)   = satPlRow(kBest);
    satSinrDbVec(u)     = satSinrDbMat(u,kBest);
    if satUnavailReasonVec(u) == "" && satSinrDbVec(u) < minUsableSinrDb
        satUnavailReasonVec(u) = "BelowSinrFloor";
    end

    if satSinrDbVec(u) > userBestSinr
        userBestSinr       = satSinrDbVec(u);
        userBestNode      = satNames(kBest);
        userBestType      = "Satellite";
        userBestDistance  = satRangeRow(kBest);
        userBestPathLoss  = satPlRow(kBest);
        userBestElevation = satElevRow(kBest);
    end

    %% ===== Απόφαση ανά πολιτική (A3, RLF) =====
    % Υποψήφιοι: τα κελιά εντός πεδίου ισχύος και οι ορατοί δορυφόροι.
    nCand = 0;
    nMax = numBs + numSats;
    candNode = strings(nMax,1); candType = strings(nMax,1);
    candSinr = -inf(nMax,1);    candDist = nan(nMax,1);
    candSinrMeas = -inf(nMax,1);
    candPl   = nan(nMax,1);     candElev = nan(nMax,1);
    for b = 1:numBs
        if isfinite(sinrDbMat(u,b))
            nCand = nCand + 1;
            candNode(nCand) = cellNames(b);
            candType(nCand) = "Terrestrial";
            candSinr(nCand) = sinrDbMat(u,b);
            candSinrMeas(nCand) = filtSinrBsMat(u,b);
            candDist(nCand) = cellRange3DMat(u,b);
            candPl(nCand)   = cellPathLossMat(u,b);
        end
    end
    for k = 1:numSats
        if isfinite(satSinrDbMat(u,k))
            nCand = nCand + 1;
            candNode(nCand) = satNames(k);
            candType(nCand) = "Satellite";
            candSinr(nCand) = satSinrDbMat(u,k);
            candSinrMeas(nCand) = filtSinrSatMat(u,k);
            candDist(nCand) = satRangeRow(k);
            candPl(nCand)   = satPlRow(k);
            candElev(nCand) = satElevRow(k);
        end
    end
    candIdx  = [find(isfinite(sinrDbMat(u,:)))'; numBs + find(isfinite(satSinrDbMat(u,:)))'];
    candNode = candNode(1:nCand); candType = candType(1:nCand);
    candSinr = candSinr(1:nCand); candDist = candDist(1:nCand);
    candSinrMeas = candSinrMeas(1:nCand);
    candPl   = candPl(1:nCand);   candElev = candElev(1:nCand);

    isTerr = (candType == "Terrestrial");
    isSat  = (candType == "Satellite");

    for pI = 1:numel(policyNamesLocal)
        polName = policyNamesLocal{pI};
        switch policyBaseLocal{pI}
            case 'Actual',          keepMask = true(nCand,1);
            case 'TerrestrialOnly', keepMask = isTerr;
            otherwise,              keepMask = isSat;
        end

        pN = ""; pPrev = ""; pStay = 0;
        st = struct('Node', "", 'Pending', "", 'Count', 0, 'Wait', 0, 'Oos', 0, 'Is', 0, 'T310', -1);
        if hasPrevDecision
            pd = prevDecision.(polName);
            pN = pd.Node(u);
            pPrev = pd.Prev(u);
            pStay = pd.Stay(u);
            st = struct('Node', pN, 'Pending', pd.Pending(u), 'Count', pd.Count(u), ...
                'Wait', pd.Wait(u), 'Oos', pd.Oos(u), 'Is', pd.Is(u), 'T310', pd.T310(u));
        end

        % Μετατόπιση MLB κάθε υποψηφίου (0 για τον εξυπηρετούντα και χωρίς MLB).
        candOffsetDb = zeros(nCand,1);
        if policyMlbLocal(pI) && hasPrevDecision && pN ~= "" && pN ~= "None"
            ldAll   = nodeLoadAll.(polName);
            rhoServ = ldAll(prevNodeIdx.(polName)(u));
            x = min(max(rhoServ - ldAll(candIdx), -1), 1);
            candOffsetDb = quantizeQOffsetDb(mlbMaxOffsetDb * x);
            candOffsetDb(candNode == pN) = 0;
        end

        % Πολιτική με ML: δεδομένα για τα χαρακτηριστικά της απόφασης (mlGateDecision).
        gateCtx = [];
        if policyMlLocal(pI) && hasPrevDecision
            gateCtx = struct('Gate', policyGateLocal{pI}, 'Kind', policyKindLocal(pI), ...
                'PpMaxStay', ppMaxStaySteps, 'User', u, 'Step', histStep, ...
                'W', histW, 'Lag', histLag, 'DtS', dtSecHys, 'NumUsers', numUsers, ...
                'HistBs', histBs, 'HistSat', histSat, ...
                'FiltBs', filtSinrBsMat(u,:), 'FiltSat', filtSinrSatMat(u,:), ...
                'Elev', satElevMat(u,:), 'Dist', cellRange3DMat(u,:), ...
                'Prev', pPrev, 'Stay', pStay, 'PrevNodes', pd.Node, ...
                'LoadNode', pd.LoadNode, 'LoadRho', pd.LoadRho);
        end

        [ci, st, ev] = mobilityStep(candNode, candSinr, candSinrMeas, keepMask, ...
            st, mob, candOffsetDb, gateCtx);

        % ci > 0: υποψήφιος του καταλόγου. ci = 0: εκτός κάλυψης ("None"), ή
        % εξυπηρετών χωρίς ζεύξη στο μοντέλο ενώ τρέχει ο T310 (χωρίς SINR).
        selNode = st.Node;
        if ci > 0
            selType = candType(ci); selSinr = candSinr(ci);
        elseif selNode == "None"
            selType = "Outage"; selSinr = NaN;
        else
            selSinr = NaN;
            if startsWith(selNode, "SAT")
                selType = "Satellite";
            else
                selType = "Terrestrial";
            end
        end
        newDecision.(polName).Node(u)    = selNode;
        newDecision.(polName).Pending(u) = st.Pending;
        newDecision.(polName).Count(u)   = st.Count;
        newDecision.(polName).Wait(u)    = st.Wait;
        newDecision.(polName).Oos(u)     = st.Oos;
        newDecision.(polName).Is(u)      = st.Is;
        newDecision.(polName).T310(u)    = st.T310;
        newDecision.(polName).Cause(u)   = ev.Cause;
        if ci > 0
            newDecision.(polName).Offset(u) = candOffsetDb(ci);
        end
        if ev.Target ~= ""
            newDecision.(polName).Attempt(u)       = ev.Target;
            newDecision.(polName).AttemptOffset(u) = candOffsetDb(candNode == ev.Target);
        end
        if ~isempty(ev.GateFeat)
            newDecision.(polName).GateFeat(u,:) = ev.GateFeat;
        end
        hoCount.(polName).Held = hoCount.(polName).Held + ev.Held;
        if ev.Rlf
            if startsWith(pN, "SAT")
                hoCount.(polName).RlfSat = hoCount.(polName).RlfSat + 1;
            else
                hoCount.(polName).RlfTerr = hoCount.(polName).RlfTerr + 1;
            end
        end
        hoCount.(polName).HoFail = hoCount.(polName).HoFail + ev.HoFail;

        % DC: βήμα του SCG μετά το MCG του ίδιου βήματος.
        if policyDcLocal(pI)
            sn = struct('Node', "None", 'Pending', "", 'Count', 0, 'Wait', 0, ...
                'Oos', 0, 'Is', 0, 'T310', -1);
            if hasPrevDecision
                sn = struct('Node', pd.SnNode(u), 'Pending', pd.SnPending(u), ...
                    'Count', pd.SnCount(u), 'Wait', 0, 'Oos', pd.SnOos(u), ...
                    'Is', pd.SnIs(u), 'T310', pd.SnT310(u));
            end
            snBefore = sn.Node ~= "None";
            [si, sn, sev] = scgStep(candNode, candSinr, candSinrMeas, isSat, sn, mob, ...
                selNode ~= "None", dcAddThrDb);
            newDecision.(polName).SnNode(u)    = sn.Node;
            newDecision.(polName).SnPending(u) = sn.Pending;
            newDecision.(polName).SnCount(u)   = sn.Count;
            newDecision.(polName).SnOos(u)     = sn.Oos;
            newDecision.(polName).SnIs(u)      = sn.Is;
            newDecision.(polName).SnT310(u)    = sn.T310;
            hoCount.(polName).SnAdd    = hoCount.(polName).SnAdd    + sev.Add;
            hoCount.(polName).SnChange = hoCount.(polName).SnChange + sev.Change;
            hoCount.(polName).ScgFail  = hoCount.(polName).ScgFail  + sev.Fail;
            % Fast MCG link recovery (TS 38.331 §5.7.3b): RLF του MCG με ενεργό SCG,
            % ο χρήστης μένει συνδεδεμένος μέσω του SCG.
            if ev.Rlf && snBefore && sn.Node ~= "None"
                hoCount.(polName).McgRlfCovered = hoCount.(polName).McgRlfCovered + 1;
            end
            snNodeMat(u,pI) = sn.Node;
            if si > 0
                snSinrMat(u,pI) = candSinr(si);
            end
        end

        % Handover = επιτυχής εκτέλεση του A3 (TR 36.839 Ορισμός 6)· όχι η
        % επανεγκατάσταση μετά από RLF/αποτυχία ή η πρώτη σύνδεση.
        if ev.Cause == "HO"
            fromSat = startsWith(pN, "SAT");
            toSat   = startsWith(selNode, "SAT");
            if ~fromSat && ~toSat
                hoType = 'TerrTerr';
            elseif fromSat && toSat
                hoType = 'SatSat';
            else
                hoType = 'TerrSat';
            end
            hoCount.(polName).(hoType) = hoCount.(polName).(hoType) + 1;
            % Ping-pong (TR 36.839 Ορισμός 5): επιστροφή στον προηγούμενο κόμβο
            % με παραμονή μικρότερη από το MTS.
            if selNode == pPrev && (pStay + 1) <= ppMaxStaySteps
                hoCount.(polName).PingPong = hoCount.(polName).PingPong + 1;
            end
            newDecision.(polName).Prev(u) = pN;
            newDecision.(polName).Stay(u) = 0;
        elseif ev.Cause == "" && selNode ~= "None" && selNode == pN
            newDecision.(polName).Prev(u) = pPrev;
            newDecision.(polName).Stay(u) = pStay + 1;
        end
        % Αλλιώς (RLF, αποτυχία, outage, πρώτη σύνδεση) η αλυσίδα παραμονής
        % ξεκινά από την αρχή (TR 36.839 §5.2.2).

        if strcmp(polName, 'Actual')
            if selNode == "None"
                % Εκτός κάλυψης: διαγνωστικά του καλύτερου υποψηφίου.
                bestNodeVec(u)         = "None";
                bestNodeTypeVec(u)     = "Outage";
                bestDistanceVec(u)     = userBestDistance;
                bestPathLossVec(u)     = userBestPathLoss;
                bestSinrDbVec(u)       = userBestSinr;
                bestElevationDegVec(u) = userBestElevation;
            elseif ci > 0
                % Ο κόμβος που εξυπηρετεί (όχι κατ' ανάγκη ο καλύτερος).
                bestNodeVec(u)         = candNode(ci);
                bestNodeTypeVec(u)     = candType(ci);
                bestDistanceVec(u)     = candDist(ci);
                bestPathLossVec(u)     = candPl(ci);
                bestSinrDbVec(u)       = candSinr(ci);
                bestElevationDegVec(u) = candElev(ci);
            else
                bestNodeVec(u)         = selNode;
                bestNodeTypeVec(u)     = selType;
                bestDistanceVec(u)     = NaN;
                bestPathLossVec(u)     = NaN;
                bestSinrDbVec(u)       = NaN;
                bestElevationDegVec(u) = NaN;
            end
        end
        polNodeMat(u,pI) = selNode;
        polTypeMat(u,pI) = selType;
        polSinrMat(u,pI) = selSinr;
    end
end

%% ------------------ Κατανομή πόρων, χωρητικότητα, ενέργεια ------------------
% maxSe: μέγιστη φασματική απόδοση NR (MCS 28, TS 38.214 Πίν. 5.1.3.1-1).
alloc = struct('BW_bs', BW_bs, 'pOutPerChainW', pOutPerChainW, ...
               'targetNormalizedSe', targetNormalizedSe, 'numBs', numBs, ...
               'maxSe', 5.5547, 'minUsableDb', minUsableSinrDb);

[capacityMbpsVec, throughputMbpsVec, serviceStateVec, nodePowerWattsVec, ...
 energyPerBitUJVec, networkEnergy] = allocateAndTally(bestNodeVec, ...
    bestNodeTypeVec, bestSinrDbVec, simParameters, satParameters, alloc);

%% ------------------ KPIs ανά πολιτική (ίδιο κανάλι) ------------------
% Η κατανομή πόρων υπολογίζεται ξανά για κάθε πολιτική. Διαθέσιμο φάσμα (MHz):
% η ζώνη των BS και κάθε δορυφορική ζώνη με τουλάχιστον έναν ορατό δορυφόρο.
spectrumTerrMHz = BW_bs/1e6;
satVisible = any(isfinite(satSinrDbMat), 1)';
spectrumSatMHz = numel(unique(satSubBand(satVisible))) * satParameters.Bandwidth/1e6;

policyKpis = struct();
policyKpis.Actual = policySummary(capacityMbpsVec, throughputMbpsVec, networkEnergy, numUsers, ...
    spectrumTerrMHz + spectrumSatMHz, hoCount.Actual, dtSecHys);

for pI = 2:numPolLocal
    polName = policyNamesLocal{pI};
    switch policyBaseLocal{pI}
        case 'Actual',          specMHz = spectrumTerrMHz + spectrumSatMHz;
        case 'TerrestrialOnly', specMHz = spectrumTerrMHz;
        otherwise,              specMHz = spectrumSatMHz;
    end
    [capP, thrP, ssP, ~, ~, netP] = allocateAndTally(polNodeMat(:,pI), polTypeMat(:,pI), ...
        polSinrMat(:,pI), simParameters, satParameters, alloc);
    if policyDcLocal(pI)
        % Ανεξάρτητος χρονοπρογραμματισμός σε MN και SN (TR 36.842 Παράρτ. F):
        % ίσο μερίδιο σε κάθε κόμβο, χωριστά για MCG και SCG.
        snType = repmat("Satellite", numUsers, 1);
        snType(snNodeMat(:,pI) == "None") = "Outage";
        [capS, thrS, ssS, ~, ~, netS] = allocateAndTally(snNodeMat(:,pI), snType, ...
            snSinrMat(:,pI), simParameters, satParameters, alloc);
        policyKpis.(polName) = dcSummary(strcmp(polName, 'DcSplit'),capP, thrP, ssP, netP, ...
            capS, thrS, ssS, netS, simParameters, numUsers, spectrumTerrMHz + spectrumSatMHz, ...
            hoCount.(polName), dtSecHys);
        continue;
    end
    policyKpis.(polName) = policySummary(capP, thrP, netP, numUsers, specMHz, ...
        hoCount.(polName), dtSecHys);
end

% Αναφορά φόρτου για το επόμενο βήμα, περιοδική (Reporting Periodicity, TS 38.423
% §9.1.3.18)· ενδιάμεσα ισχύει η τελευταία. LoadAge = βήματα από την αναφορά.
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'LoadReportPeriodMs')
    nLoadRep = max(1, round(simParameters.Mobility.LoadReportPeriodMs / (dtSecHys*1000)));
else
    nLoadRep = 1;
end
for pI = 1:numPolLocal
    polName = policyNamesLocal{pI};
    if hasPrevDecision && isfield(prevDecision.(polName), 'LoadAge') && ...
            prevDecision.(polName).LoadAge + 1 < nLoadRep
        newDecision.(polName).LoadNode = prevDecision.(polName).LoadNode;
        newDecision.(polName).LoadRho  = prevDecision.(polName).LoadRho;
        newDecision.(polName).LoadAge  = prevDecision.(polName).LoadAge + 1;
    else
        [ldNode, ldRho] = virtualLoad(polNodeMat(:,pI), polSinrMat(:,pI), ...
            targetNormalizedSe, alloc.maxSe, minUsableSinrDb);
        newDecision.(polName).LoadNode = ldNode;
        newDecision.(polName).LoadRho  = ldRho;
        newDecision.(polName).LoadAge  = 0;
    end
end

%% ------------------ Κατάσταση καναλιού για το επόμενο βήμα ------------------
newChannelState.UserGeo         = user_geo;
newChannelState.IsLOS           = losMat;
newChannelState.LosLatent       = losLatentMat;
newChannelState.SatShadowLatent = satShadowLatentMat;
newChannelState.ShadowFading_dB = sfMat;
newChannelState.SfLos_dB        = sfLosMat;
newChannelState.SfNlos_dB       = sfNlosMat;
newChannelState.LosSoft         = losSoftMat;
newChannelState.FiltSinrBs_dB   = filtSinrBsMat;
newChannelState.FiltSinrSat_dB  = filtSinrSatMat;
if hasMl
    slot = mod(histStep - 1, histW) + 1;
    histBs(slot,:,:)  = filtSinrBsMat;
    histSat(slot,:,:) = filtSinrSatMat;
    newChannelState.HistFiltBs  = histBs;
    newChannelState.HistFiltSat = histSat;
    newChannelState.HistStep    = histStep;
end
% Γεωμετρία των υποψηφίων, για την καταγραφή των handovers.
newChannelState.SatElev_deg     = satElevMat;
newChannelState.BsDist_m        = cellRange3DMat;

end

function [chosenIdx, st, ev] = mobilityStep(candNode, candSinr, candSinrMeas, keepMask, st, P, candOffsetDb, gateCtx)
% Ένα βήμα κινητικότητας για έναν χρήστη σε μία πολιτική.
%   (1) RLM (TS 38.331 §5.3.10): N310/N311/T310 -> RLF.
%   (2) Event A3 (TS 38.331 §5.5.4.4): M_n + Ocn - Hys > M_p για NTtt βήματα,
%       στη φιλτραρισμένη μέτρηση.
%   (3) Αποτυχία HO (TR 36.839 Ορισμός 3): εντολή ενώ τρέχει ο T310, RLF ενώ
%       ισχύει η συνθήκη εισόδου, ή στόχος κάτω από το Qout.
%   (4) Μετά από RLF/αποτυχία: επανεγκατάσταση (§5.3.7) με τον κανόνα της
%       πρώτης σύνδεσης (cellSelection).
%   (5) Μόνο με gateCtx (*ML): το μοντέλο αποφασίζει αν θα σταλεί η εντολή·
%       αλλιώς νέα απόφαση μετά από NReport βήματα (st.Wait).
% candSinr: πραγματικό SINR (RLM, στόχος, επιλογή κελιού)· candSinrMeas:
% φιλτραρισμένη μέτρηση που βλέπει ο A3.
% chosenIdx: δείκτης υποψηφίου (0 = εκτός κάλυψης ή εξυπηρετών χωρίς ζεύξη).
% ev: Cause ("", "HO", "RLF", "HOF"), Rlf, HoFail, Held, GateFeat, Target.
ev = struct('Cause', "", 'Rlf', false, 'HoFail', false, 'Held', false, 'GateFeat', [], ...
    'Target', "");
valid = keepMask & isfinite(candSinrMeas);
pN = st.Node;

if pN == "" || pN == "None"
    [chosenIdx, st] = cellSelection(candNode, candSinr, candSinrMeas, valid, P, st);
    return;
end

iPrev = find(candNode == pN & valid, 1);
if isempty(iPrev)
    iPrev = 0;
    sServ = -Inf;   % η ζεύξη δεν μετριέται: out-of-sync
else
    sServ = candSinr(iPrev);
end

% ---- (1) RLM ----
if st.T310 < 0
    if sServ < P.QoutDb
        st.Oos = st.Oos + 1;
    else
        st.Oos = 0;
    end
    if st.Oos >= P.N310
        st.T310 = 0; st.Oos = 0; st.Is = 0;
    end
else
    st.T310 = st.T310 + 1;
    if sServ > P.QinDb
        st.Is = st.Is + 1;
    else
        st.Is = 0;
    end
    if st.Is >= P.N311
        st.T310 = -1; st.Is = 0;
    elseif st.T310 >= P.NT310
        ev.Rlf = true;
        ev.HoFail = st.Count > 0;    % RLF στην κατάσταση 2
        ev.Cause = "RLF";
        [chosenIdx, st] = cellSelection(candNode, candSinr, candSinrMeas, valid, P, st);
        return;
    end
end

% ---- (2) A3 ----
chosenIdx = iPrev;
masked = candSinrMeas + candOffsetDb;
masked(~valid) = -Inf;
if iPrev > 0
    masked(iPrev) = candSinrMeas(iPrev);   % η μετατόπιση του ίδιου είναι 0
    mServ = candSinrMeas(iPrev);
else
    mServ = -Inf;
end
[mBest, iBest] = max(masked);
if ~isfinite(mBest) || iBest == iPrev
    st.Pending = ""; st.Count = 0; st.Wait = 0;
    return;
end
if mBest - P.HysDb > mServ
    if st.Pending == candNode(iBest)
        st.Count = st.Count + 1;
        st.Wait  = max(st.Wait - 1, 0);
    else
        st.Pending = candNode(iBest);
        st.Count = 1;
        st.Wait  = 0;
    end
else
    st.Pending = ""; st.Count = 0; st.Wait = 0;   % η συνθήκη έπαψε να ισχύει
    return;
end
if st.Count < P.NTtt || st.Wait > 0
    return;
end

% ---- (5) Απόφαση του σταθμού ----
if ~isempty(gateCtx)
    [holdHo, ev.GateFeat] = mlGateDecision(gateCtx, candNode, iBest, candOffsetDb, pN);
    if holdHo
        ev.Held = true;
        st.Wait = P.NReport;
        return;
    end
end

% ---- (3) Εντολή handover ----
ev.Target = candNode(iBest);
st.Pending = ""; st.Count = 0; st.Wait = 0;
if st.T310 >= 0
    ev.HoFail = true;   % κριτήριο 1: ο T310 τρέχει, ο χρήστης μένει
    return;
end
if candSinr(iBest) < P.QoutDb
    ev.HoFail = true;   % κατάσταση 3: ο στόχος δεν είναι προσβάσιμος
    ev.Cause = "HOF";
    [chosenIdx, st] = cellSelection(candNode, candSinr, candSinrMeas, valid, P, st);
    return;
end
ev.Cause = "HO";
chosenIdx = iBest;
st.Node = candNode(iBest);
st.Oos = 0; st.Is = 0; st.T310 = -1;   % νέα ζεύξη, νέα παρακολούθηση
end

function [chosenIdx, st] = cellSelection(candNode, candSinr, candSinrMeas, valid, P, st)
% Επιλογή κελιού χωρίς εξυπηρετούντα (πρώτη σύνδεση, μετά από outage/RLF/
% αποτυχία): ο καλύτερος στη μέτρηση με SINR >= MCS 0· κανένας: outage.
usable = valid & (candSinr >= P.MinUsableDb);
st.Pending = ""; st.Count = 0; st.Wait = 0; st.Oos = 0; st.Is = 0; st.T310 = -1;
if ~any(usable)
    chosenIdx = 0;
    st.Node = "None";
    return;
end
m = candSinrMeas;
m(~usable) = -Inf;
[~, chosenIdx] = max(m);
st.Node = candNode(chosenIdx);
end

function [holdHo, feat] = mlGateDecision(g, candNode, iBest, candOffsetDb, servNode)
% Απόφαση του σταθμού για αναφορά A3 servNode -> candNode(iBest): ίδια
% χαρακτηριστικά με την καταγραφή (handoverFeatures), P(ping-pong) από το
% μοντέλο, αναμονή αν ξεπερνά το κατώφλι. Κανόνας "history": αναμονή αν η
% εκτέλεση θα ήταν ping-pong κατά το ιστορικό του χρήστη (TS 38.300 §15.5.2.4).
tgtNode = candNode(iBest);
if g.Kind == "history"
    holdHo = tgtNode == g.Prev && (g.Stay + 1) <= g.PpMaxStay;
    feat = [];
    return;
end
[satF, idxF] = gateNodeIndex(servNode);
[satT, idxT] = gateNodeIndex(tgtNode);
mF = [gateMeas(g, satF, idxF, 0), gateMeas(g, satF, idxF, g.Lag), gateMeas(g, satF, idxF, 2*g.Lag)];
mT = [gateMeas(g, satT, idxT, 0), gateMeas(g, satT, idxT, g.Lag), gateMeas(g, satT, idxT, 2*g.Lag)];
feat = handoverFeatures(struct('SinrFrom', mF, 'SinrTo', mT, ...
    'LagS', g.Lag*g.DtS, 'NumUsers', g.NumUsers, ...
    'ReturnToPrev', tgtNode == g.Prev, 'StayFrom_s', (g.Stay + 1)*g.DtS, ...
    'ElevFrom', gateGeom(g.Elev, satF, idxF), 'ElevTo', gateGeom(g.Elev, satT, idxT), ...
    'DistFrom', gateGeom(g.Dist, ~satF, idxF), 'DistTo', gateGeom(g.Dist, ~satT, idxT), ...
    'LoadFrom', sum(g.LoadRho(g.LoadNode == servNode)), ...
    'LoadTo', sum(g.LoadRho(g.LoadNode == tgtNode)), ...
    'UsersOnFrom', sum(g.PrevNodes == servNode), 'UsersOnTo', sum(g.PrevNodes == tgtNode), ...
    'MlbOffset', candOffsetDb(iBest)));
x = [feat, satF, satT];
p = handoverGatePredict(g.Gate.Model, x(g.Gate.FeatIdx));
holdHo = p > g.Gate.Threshold;
end

function v = gateMeas(g, isSat, idx, lag)
% Φιλτραρισμένη μέτρηση lag βήματα πριν· NaN εκτός ιστορικού.
v = NaN;
if lag >= g.W || g.Step - lag < 1
    return;
end
if lag == 0
    if isSat, v = g.FiltSat(idx); else, v = g.FiltBs(idx); end
    return;
end
slot = mod(g.Step - lag - 1, g.W) + 1;
if isSat
    v = g.HistSat(slot, g.User, idx);
else
    v = g.HistBs(slot, g.User, idx);
end
end

function v = gateGeom(row, applies, idx)
% Ανύψωση (δορυφόρος) ή απόσταση (σταθμός) του κόμβου· NaN για τον άλλο τύπο.
if applies
    v = row(idx);
else
    v = NaN;
end
end

function [isSat, idx] = gateNodeIndex(node)
isSat = startsWith(node, "SAT-");
if isSat
    idx = str2double(extractAfter(node, "SAT-"));
else
    idx = str2double(extractAfter(node, "BS"));
end
end

function [capacityMbpsVec, throughputMbpsVec, serviceStateVec, ...
          nodePowerWattsVec, energyPerBitUJVec, networkEnergy] = ...
          allocateAndTally(nodeVec, typeVec, sinrDbVec, simParameters, satParameters, alloc)
% Κατανομή εύρους ζώνης, χωρητικότητα, κατάσταση υπηρεσίας, ισχύς και ενέργεια
% ανά bit για δεδομένη ανάθεση χρηστών σε κόμβους.
numUsers = numel(nodeVec);
numBs    = alloc.numBs;

capacityMbpsVec   = nan(numUsers,1);
throughputMbpsVec = zeros(numUsers,1);
serviceStateVec   = strings(numUsers,1);
nodePowerWattsVec = nan(numUsers,1);
energyPerBitUJVec = inf(numUsers,1);
% Πλήθος χρηστών ανά κόμβο.
[~, ~, nodeGrp] = unique(nodeVec);
nodeCount = accumarray(nodeGrp(:), 1);

for u = 1:numUsers
    % Εκτός κάλυψης: χωρητικότητα NaN, ρυθμός 0, ενέργεια/bit Inf.
    if typeVec(u) == "Outage"
        capacityMbpsVec(u)   = NaN;
        throughputMbpsVec(u) = 0;
        nodePowerWattsVec(u) = NaN;
        energyPerBitUJVec(u) = Inf;
        serviceStateVec(u)   = "Outage";
        continue;
    end

    usersOnThisNode = nodeCount(nodeGrp(u));

    % Εύρος και ισχύς κόμβου: EARTH (Auer et al. 2011) για BS, γραμμικός
    % ενισχυτής για δορυφόρο.
    if typeVec(u) == "Terrestrial"
        nodeBW = alloc.BW_bs;
        nodePowerW = simParameters.Power.NumTrx * ...
            (simParameters.Power.P0 + simParameters.Power.DeltaP * alloc.pOutPerChainW);
    else
        nodeBW = satParameters.Bandwidth;
        pOutW  = 10^((satParameters.TxPower - 30)/10);
        nodePowerW = satParameters.Power.Pfix + pOutW / satParameters.Power.EtaPA;
    end

    % Ίσος διαμοιρασμός εύρους.
    B_user = nodeBW / usersOnThisNode;

    % Shannon με όριο τη μέγιστη SE του NR (MCS 28)· κάτω από το MCS 0 (π.χ.
    % όσο τρέχει ο T310) ρυθμός 0.
    sinr_lin = 10^(sinrDbVec(u)/10);
    if sinrDbVec(u) >= alloc.minUsableDb
        spectralEfficiency = min(log2(1 + sinr_lin), alloc.maxSe);
    else
        spectralEfficiency = 0;
    end
    capacity = B_user * spectralEfficiency;   % bits/s

    capacityMbpsVec(u)   = capacity * 1e-6;    % Mbps
    throughputMbpsVec(u) = capacity * 1e-6;

    % Κατάσταση υπηρεσίας, επί του συνολικού εύρους του κόμβου.
    if (capacity / nodeBW) >= alloc.targetNormalizedSe
        serviceStateVec(u) = "Served";
    else
        serviceStateVec(u) = "BelowTarget";
    end

    % Ενέργεια ανά bit: ισομερές μερίδιο της ισχύος του κόμβου / ρυθμός (µJ/bit).
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

%% ------------------ Ενεργειακή απόδοση δικτύου ------------------
% bit/J = συνολικός ρυθμός / ισχύς ανά κόμβο, σε δύο εμβέλειες: πλήρης κόμβος
% και μόνο ενισχυτής (RF, συμμετρική για BS και δορυφόρο).
networkEnergy = struct();
totalRateBps      = sum(throughputMbpsVec) * 1e6;   % παραδοθέντα bits
pOutSatW          = 10^((satParameters.TxPower - 30)/10);
bsFullW           = simParameters.Power.NumTrx * (simParameters.Power.P0 + simParameters.Power.DeltaP*alloc.pOutPerChainW);
bsRfW             = simParameters.Power.NumTrx * simParameters.Power.DeltaP * alloc.pOutPerChainW;
bsSleepW          = simParameters.Power.NumTrx * simParameters.Power.Psleep;
satFullW          = satParameters.Power.Pfix + pOutSatW/satParameters.Power.EtaPA;
satRfW            = pOutSatW / satParameters.Power.EtaPA;

activeBs  = sum(ismember("BS" + string((1:numBs)'), nodeVec));
idleBs    = numBs - activeBs;
satActive = numel(unique(nodeVec(typeVec == "Satellite")));   % δορυφόροι με χρήστες

% Πλήρης: ενεργοί BS κατά EARTH, αδρανείς σε Psleep, δορυφόροι μόνο όταν εξυπηρετούν.
networkEnergy.TotalPower_W    = activeBs*bsFullW + idleBs*bsSleepW + satActive*satFullW;
networkEnergy.TotalPowerRf_W  = activeBs*bsRfW + satActive*satRfW;
networkEnergy.TotalRate_bps   = totalRateBps;
networkEnergy.ActiveBs        = activeBs;
networkEnergy.IdleBs          = idleBs;
networkEnergy.SatActive       = satActive;
networkEnergy.BitPerJoule     = totalRateBps / networkEnergy.TotalPower_W;
networkEnergy.BitPerJouleRf   = totalRateBps / networkEnergy.TotalPowerRf_W;
networkEnergy.ServedUsers     = sum(serviceStateVec == "Served");
networkEnergy.BelowTargetUsers= sum(serviceStateVec == "BelowTarget");
networkEnergy.OutageUsers     = sum(serviceStateVec == "Outage");
networkEnergy.SatelliteUsers  = sum(typeVec == "Satellite");
networkEnergy.TargetNormSe    = alloc.targetNormalizedSe;
end

function q = quantizeQOffsetDb(offsetDb)
% Πλησιέστερη τιμή του Q-OffsetRange (TS 38.331): -24, -22, ..., -8, -6..6, 8, ..., 24 dB.
% Στρογγύλευση στο μέτρο, ώστε q(-x) = -q(x).
levels = [0:6, 8:2:24];
[~, k] = min(abs(levels - abs(offsetDb(:))), [], 2);
q = reshape(sign(offsetDb(:)) .* levels(k)', size(offsetDb));
end

function [nodes, rho] = virtualLoad(nodeVec, sinrDbVec, targetNormalizedSe, maxSe, minUsableDb)
% Εικονικός φόρτος (Awada et al., PIMRC 2010, εξ. 1-3):
%     κ_u = target / SE_u,   ρ_c = Σ_{u στον c} κ_u   (ρ_c <= 1: εξυπηρετεί όλους).
% Με full buffer η κατάληψη PRB δεν ξεχωρίζει γεμάτο από άδειο κόμβο· ο ρ_c ναι.
% Χρήστες κάτω από το MCS 0 δεν μετρούν.
served = nodeVec ~= "None" & nodeVec ~= "" & isfinite(sinrDbVec) & sinrDbVec >= minUsableDb;
nodes  = unique(nodeVec(served));
rho    = zeros(numel(nodes),1);
se     = min(log2(1 + 10.^(sinrDbVec/10)), maxSe);
for i = 1:numel(nodes)
    onNode = served & nodeVec == nodes(i);
    rho(i) = sum(targetNormalizedSe ./ se(onNode));
end
end

function s = policySummary(capacityMbpsVec, throughputMbpsVec, networkEnergy, numUsers, spectrumMHz, ho, dtS)
% KPIs μιας πολιτικής για ένα βήμα. Μέση χωρητικότητα στους εξυπηρετούμενους·
% συνολικός ρυθμός και ποσοστά επί όλων των χρηστών.
served = ~isnan(capacityMbpsVec);
s = struct();
if any(served)
    s.MeanCapacity_Mbps = mean(capacityMbpsVec(served));
else
    s.MeanCapacity_Mbps = NaN;
end
s.TotalRate_Mbps   = sum(throughputMbpsVec);
s.FracServed       = networkEnergy.ServedUsers / numUsers;
s.FracBelowTarget  = networkEnergy.BelowTargetUsers / numUsers;
s.FracOutage       = networkEnergy.OutageUsers / numUsers;
s.FracSatellite    = networkEnergy.SatelliteUsers / numUsers;
s.TotalPower_W     = networkEnergy.TotalPower_W;
s.TotalPowerRf_W   = networkEnergy.TotalPowerRf_W;
s.BitPerJoule      = networkEnergy.BitPerJoule;
s.BitPerJouleRf    = networkEnergy.BitPerJouleRf;
s.Spectrum_MHz     = spectrumMHz;   % διαθέσιμο (όχι χρησιμοποιούμενο) φάσμα
% Ρυθμοί γεγονότων ανά χρήστη ανά λεπτό (TR 36.839 Ορισμοί 2, 3, 5).
perUserMin = 60 / (numUsers * dtS);
s.HoRate_perUserMin         = (ho.TerrTerr + ho.TerrSat + ho.SatSat) * perUserMin;
s.HoRateTerrTerr_perUserMin = ho.TerrTerr * perUserMin;
s.HoRateTerrSat_perUserMin  = ho.TerrSat  * perUserMin;
s.HoRateSatSat_perUserMin   = ho.SatSat   * perUserMin;
s.PingPongRate_perUserMin   = ho.PingPong * perUserMin;
s.RlfRate_perUserMin        = (ho.RlfTerr + ho.RlfSat) * perUserMin;
s.RlfRateTerr_perUserMin    = ho.RlfTerr * perUserMin;
s.RlfRateSat_perUserMin     = ho.RlfSat  * perUserMin;
s.HoFailRate_perUserMin     = ho.HoFail  * perUserMin;
s.HoHeldRate_perUserMin     = ho.Held    * perUserMin;   % αναμονές του μοντέλου (*ML)
% Διακοπή υπηρεσίας: χρήστες με μηδενικό ρυθμό (εκτός κάλυψης ή κάτω από το MCS 0,
% π.χ. όσο τρέχει ο T310).
s.FracZeroRate              = sum(throughputMbpsVec <= 0) / numUsers;
% SCG (μόνο DC): προσθήκες, αλλαγές PSCell, αποτυχίες SCG (TS 38.331 §5.7.3),
% RLF του MCG που καλύφθηκαν από το SCG (§5.7.3b).
s.SnAddRate_perUserMin         = ho.SnAdd         * perUserMin;
s.SnChangeRate_perUserMin      = ho.SnChange      * perUserMin;
s.ScgFailRate_perUserMin       = ho.ScgFail       * perUserMin;
s.McgRlfCoveredRate_perUserMin = ho.McgRlfCovered * perUserMin;
end

function s = dcSummary(isSplit, capM, thrM, ssM, netM, capS, thrS, ssS, netS, ...
    simParameters, numUsers, spectrumMHz, ho, dtS)
% KPIs πολιτικής DC από τις δύο κατανομές (MCG, SCG) του βήματος.
%   split bearer: R = R_MCG + R_SCG (TS 37.340 §4.2.2)
%   duplication:  R = max(R_MCG, R_SCG): υπόθεση, η ταχύτερη ζεύξη παραδίδει
%                 πρώτη και τα αντίγραφα απορρίπτονται (TS 38.323 §5.11)
% Εξυπηρετούμενος: τουλάχιστον μία ζεύξη πάνω από τον στόχο· εκτός κάλυψης:
% καμία ζεύξη. Ισχύς: οι σταθμοί από το MCG, οι δορυφόροι από το SCG.
legM = ~isnan(capM);
legS = ~isnan(capS);
if isSplit
    thr = thrM + thrS;
else
    thr = max(thrM, thrS);
end
cap = thr;
cap(~legM & ~legS) = NaN;
bsSleepW = simParameters.Power.NumTrx * simParameters.Power.Psleep;
net = struct();
% Το netS μετρά όλους τους σταθμούς ως αδρανείς· αφαιρούνται.
net.TotalPower_W     = netM.TotalPower_W + netS.TotalPower_W - netS.IdleBs*bsSleepW;
net.TotalPowerRf_W   = netM.TotalPowerRf_W + netS.TotalPowerRf_W;
net.BitPerJoule      = sum(thr)*1e6 / net.TotalPower_W;
net.BitPerJouleRf    = sum(thr)*1e6 / net.TotalPowerRf_W;
net.ServedUsers      = sum(ssM == "Served" | ssS == "Served");
net.OutageUsers      = sum(~legM & ~legS);
net.BelowTargetUsers = numUsers - net.ServedUsers - net.OutageUsers;
net.SatelliteUsers   = sum(legS);   % χρήστες με SCG
s = policySummary(cap, thr, net, numUsers, spectrumMHz, ho, dtS);
end

function [chosenIdx, st, ev] = scgStep(candNode, candSinr, candSinrMeas, isSat, st, P, mcgUp, addThrDb)
% Ένα βήμα του SCG (δορυφόρος) ενός χρήστη σε πολιτική DC:
%   χωρίς MCG:  αποδέσμευση (χωρίς σύνδεση RRC δεν υπάρχει SCG)
%   χωρίς SN:   προσθήκη με A4 (TS 38.331 §5.5.4.5): M_n - Hys > Thresh για NTtt
%               βήματα· στόχος κάτω από το Qout στην εκτέλεση = αποτυχία SCG
%               (αποτυχία reconfigurationWithSync, §5.7.3.2)
%   με SN:      RLM του PSCell και A3 μεταξύ δορυφόρων (mobilityStep)· RLF ή
%               αποτυχία αλλαγής = αποτυχία SCG (§5.3.10.3, §5.7.3) και αποδέσμευση
% chosenIdx: δείκτης του PSCell στους υποψηφίους (0 = χωρίς SCG ή χωρίς ζεύξη).
ev = struct('Add', false, 'Change', false, 'Fail', false);
chosenIdx = 0;
if ~mcgUp
    st = resetScg(st);
    return;
end
valid = isSat & isfinite(candSinrMeas);
if st.Node == "None"
    m = candSinrMeas;
    m(~valid) = -Inf;
    [mBest, iBest] = max(m);
    if ~isfinite(mBest) || mBest - P.HysDb <= addThrDb
        st = resetScg(st);
        return;
    end
    if st.Pending == candNode(iBest)
        st.Count = st.Count + 1;
    else
        st.Pending = candNode(iBest);
        st.Count = 1;
    end
    if st.Count < P.NTtt
        return;
    end
    if candSinr(iBest) < P.QoutDb
        ev.Fail = true;
        st = resetScg(st);
        return;
    end
    ev.Add = true;
    chosenIdx = iBest;
    st = resetScg(st);
    st.Node = candNode(iBest);
    return;
end
[chosenIdx, st, mev] = mobilityStep(candNode, candSinr, candSinrMeas, valid, st, P, ...
    zeros(numel(candNode), 1), []);
if mev.Cause == "RLF" || mev.Cause == "HOF"
    ev.Fail = true;
    chosenIdx = 0;
    st = resetScg(st);
elseif mev.Cause == "HO"
    ev.Change = true;
end
end

function st = resetScg(st)
st.Node = "None"; st.Pending = ""; st.Count = 0; st.Wait = 0;
st.Oos = 0; st.Is = 0; st.T310 = -1;
end

function [isLos, latent] = spatiallyConsistentLos(pLos, moveDistance, hasPrevState, prevLatent)
% Χωρικά συνεπής κατάσταση LOS (TR 38.901 §7.6.3.3): γκαουσιανή AR(1) με
% απόσταση συσχέτισης 50 m (Πίν. 7.6.3.1-2), σε ομοιόμορφη, σύγκριση με pLos.
dCorrLos = 50;
if ~hasPrevState
    latent = randn();
else
    r = exp(-moveDistance / dCorrLos);
    latent = r*prevLatent + sqrt(1 - r^2)*randn();
end
u = 0.5*erfc(-latent/sqrt(2));   % γκαουσιανή -> ομοιόμορφη στο (0,1)
isLos = u < pLos;
end

function rho = shadowFadingCorrelation(moveDistance, scenario, isLos)
% Αυτοσυσχέτιση σκίασης exp(-Δd/d_corr) (TR 38.901 εξ. 7.4-5, Πίν. 7.5-6).
switch scenario
    case 'UMi'
        dCorr = 10*isLos + 13*~isLos;
    otherwise
        dCorr = 37*isLos + 50*~isLos;
end
rho = exp(-moveDistance / dCorr);
end

function fadeDb = smallScaleFadingDb(isLos, KdB, nAvg)
% Small-scale fading σε dB, E[|h|^2] = 1: Rayleigh (NLOS) ή Rician με K dB (LOS).
% nAvg: μέσος όρος γραμμικής ισχύος σε nAvg δείγματα της περιόδου μέτρησης
% (TS 38.133 Πίν. 9.2.5.2-1).
if nargin < 3 || isempty(nAvg)
    nAvg = 1;
end
if ~isLos
    if nAvg == 1
        fadeDb = 10*log10(-log(rand()));
        return;
    end
    powLin = mean(-log(rand(1,nAvg)));
else
    Klin  = 10^(KdB/10);
    s     = sqrt(Klin/(Klin+1));       % πλάτος LOS συνιστώσας
    sigma = sqrt(1/(2*(Klin+1)));      % τυπ. απόκλιση ανά διάσταση scatter
    h     = (s + sigma*randn(1,nAvg)) + 1i*(sigma*randn(1,nAvg));
    if nAvg == 1
        fadeDb = 20*log10(abs(h));     % s^2 + 2*sigma^2 = 1
        return;
    end
    powLin = mean(abs(h).^2);
end
fadeDb = 10*log10(powLin);
end

function isValid = isValidTerrestrialLink(d2D, hUT, fcHz)
% Πεδίο ισχύος UMa/UMi (TR 38.901 Πίν. 7.4.1-1): 10 m <= d2D <= 5 km,
% 1.5 m <= hUT <= 22.5 m, 0.5-100 GHz.
isValid = d2D  >= 10   & d2D  <= 5000 & ...
          hUT  >= 1.5  & hUT  <= 22.5 & ...
          fcHz >= 0.5e9 & fcHz <= 100e9;
end

function [b0, m, omega] = shadowedRicianElevParams(elevDeg)
% Παράμετροι Shadowed Rician από την ανύψωση (Abdi et al. 2003, εξ. 19), 20°-80°.
th = elevDeg;
b0    = -4.7943e-8*th^3 + 5.5784e-6*th^2 - 2.1344e-4*th + 3.2710e-2;
m     =  6.3739e-5*th^3 + 5.8533e-4*th^2 - 1.5973e-1*th + 3.5156;
omega =  1.4428e-5*th^3 - 2.3798e-3*th^2 + 1.2702e-1*th - 1.4864;
end

function [b0, m, omega] = shadowedRicianStateParams(state)
% Σταθερές καταστάσεις σκίασης (Abdi et al. 2003, Πίν. III) - ανάλυση ευαισθησίας.
switch lower(string(state))
    case "light"
        b0 = 0.158; m = 19.4;  omega = 1.29;
    case "average"
        b0 = 0.126; m = 10.1;  omega = 0.835;
    case "heavy"
        b0 = 0.063; m = 0.739; omega = 8.97e-4;
    otherwise
        error('shadowedRicianStateParams:UnknownState', ...
            'Άγνωστη κατάσταση σκίασης "%s" - δεκτές: "light", "average", "heavy".', state);
end
end

function fadeDb = shadowedRicianFadingDb(b0, m, omega, nAvg, uShadow)
% Shadowed Rician σε dB (Abdi et al. 2003, εξ. 1): Rayleigh ισχύος 2*b0 +
% LOS με πλάτος Nakagami-m ισχύος omega. Η σκίαση μένει σταθερή μέσα στο
% βήμα· ο μέσος όρος nAvg αφορά μόνο τη σκεδαζόμενη συνιστώσα.
if nargin < 4 || isempty(nAvg)
    nAvg = 1;
end
% Z^2 ~ Gamma(m, omega/m): από το εκατοστημόριο uShadow, αλλιώς ανεξάρτητο δείγμα.
if nargin >= 5 && ~isempty(uShadow)
    uShadow = min(max(uShadow, eps), 1 - eps);
    losAmp  = sqrt((omega/m) * gammaincinv(uShadow, m));   % |Z|, E[Z^2] = omega
else
    losAmp  = sqrt(gammaRand(m, omega/m));
end
scatter = sqrt(b0)*(randn(1,nAvg) + 1i*randn(1,nAvg));   % E[|A|^2] = 2*b0
if nAvg == 1
    fadeDb = 20*log10(abs(losAmp + scatter));
else
    fadeDb = 10*log10(mean(abs(losAmp + scatter).^2));
end
end

function f = l3Filter(mNow, fPrev, alpha)
% Φίλτρο L3 (TS 38.331 §5.5.3.2): F_n = (1 - a)*F_{n-1} + a*M_n, σε dB· F_0 = M_1.
if ~isfinite(fPrev)
    f = mNow;
else
    f = (1 - alpha)*fPrev + alpha*mNow;
end
end

function v = measParam(p, name, dflt)
% Ανάγνωση πεδίου από το προαιρετικό υπο-struct Meas, με προεπιλογή.
v = dflt;
if isfield(p, 'Meas') && isfield(p.Meas, name) && ~isempty(p.Meas.(name))
    v = p.Meas.(name);
end
end

function x = gammaRand(shape, scale)
% Δείγμα Gamma (Marsaglia & Tsang 2000), χωρίς Statistics Toolbox.
if shape < 1
    x = gammaRand(shape + 1, scale) * rand()^(1/shape);
    return;
end
d = shape - 1/3;
c = 1/sqrt(9*d);
while true
    v = -1;
    while v <= 0
        z = randn();
        v = (1 + c*z)^3;
    end
    u = rand();
    if log(u) < 0.5*z^2 + d - d*v + d*log(v)
        x = d*v*scale;
        return;
    end
end
end

function pLos = losProbability38901(d2D, hUT, scenario)
% Πιθανότητα LOS ζεύξης BS-χρήστη (TR 38.901 Πίν. 7.4.2-1).
pLos = ones(size(d2D));
far  = d2D > 18;
d    = d2D(far);
switch scenario
    case 'UMi'
        pLos(far) = 18./d + exp(-d/36) .* (1 - 18./d);
    case 'UMa'
        h = hUT(far);
        Cprime = zeros(size(d));
        hi = h > 13;
        g = 1.25e-6 * d(hi).^3 .* exp(-d(hi)/150);
        Cprime(hi) = ((h(hi) - 13)/10).^1.5 .* g;
        pLos(far) = (18./d + exp(-d/63) .* (1 - 18./d)) .* (1 + Cprime);
    otherwise
        error('losProbability38901:UnsupportedScenario', ...
            'Άγνωστο PathLoss.Scenario "%s" - η πιθανότητα LOS (TR 38.901 §7.4.2) είναι ορισμένη μόνο για "UMa" και "UMi".', ...
            scenario);
end
end

function pla_dB = gasAttenuationSlantP676(freqHz, elevDeg)
% Απόσβεση αερίων (TR 38.811 εξ. 6.6-8): PLA = A_zenith / sin(ε), με
% A_zenith = γ_o·h_o + γ_w·h_w (γ από gaspl, ITU-R P.676 Annex 1· ύψη από
% Annex 2· ατμόσφαιρα αναφοράς ITU-R P.835). Το A_zenith κρατείται ανά συχνότητα.
persistent cacheFreqHz cacheAzenith
if ~isempty(cacheFreqHz) && cacheFreqHz == freqHz
    pla_dB = cacheAzenith / sind(elevDeg);
    return;
end
TcRef    = 15;        % °C (= 288.15 K)
TKRef    = 288.15;    % K
pPaRef   = 101325;    % Pa
pHpaRef  = 1013.25;   % hPa
rhoRef   = 7.5;       % g/m^3 (υδρατμοί)

freqGHz = freqHz / 1e9;

gammaDry = gaspl(1000, freqHz, TcRef, pPaRef, 0);       % dB/km, μόνο O2
gammaTot = gaspl(1000, freqHz, TcRef, pPaRef, rhoRef);  % dB/km, O2 + υδρατμοί
gammaWet = gammaTot - gammaDry;

[ho, hw] = equivalentHeightsP676(freqGHz, TKRef, pHpaRef, rhoRef);

Azenith = gammaDry*ho + gammaWet*hw;   % dB
cacheFreqHz  = freqHz;
cacheAzenith = Azenith;
pla_dB  = Azenith / sind(elevDeg);     % dB
end

function [ho, hw] = equivalentHeightsP676(freqGHz, T_K, p_hPa, rho)
% Ισοδύναμα ύψη O2/υδρατμών, ITU-R P.676-12 Annex 2 εξ. (30)-(38).
e_hPa = rho * T_K / 216.7;
rp = (p_hPa + e_hPa) / 1013.25;

% -- Οξυγόνο --
A_o = 0.7832 + 0.00709*(T_K - 273.15);

t1 = (5.1040 / (1 + 0.066*rp^(-2.3))) * ...
     exp(-((freqGHz - 59.7) / (2.87 + 12.4*exp(-7.9*rp)))^2);

fi3 = [118.750334 368.498246 424.763020 487.249273 715.392902 773.839490 834.145546];
ci3 = [0.1597 0.1066 0.1325 0.1242 0.0938 0.1448 0.1374];
t2 = sum( (ci3 * exp(2.12*rp)) ./ ((freqGHz - fi3).^2 + 0.025*exp(2.2*rp)) );

t3 = (0.0114*freqGHz * (15.02*freqGHz^2 - 1353*freqGHz + 5.333e4)) / ...
     ((1 + 0.14*rp^(-2.6)) * (freqGHz^3 - 151.3*freqGHz^2 + 9629*freqGHz - 6803));

ho = (6.1*A_o / (1 + 0.17*rp^(-1.1))) * (1 + t1 + t2 + t3);
if freqGHz >= 70
    ho = min(ho, 10.7*rp^0.3);
end

% -- Υδρατμοί --
A_w = 1.9298 - 0.04166*(T_K - 273.15) + 0.0517*e_hPa;
B_w = 1.1674 - 0.00622*(T_K - 273.15) + 0.0063*e_hPa;
sigma_w = 1 + 1.013 / (1 + exp(-8.6*(rp - 0.57)));

fi4 = [22.235080 183.310087 325.152888 380.197353 439.150807 448.001085 ...
       474.689092 488.490108 556.935985 620.700870 752.033113 916.171582 ...
       970.315022 987.926764];
ai4 = [1.52 7.62 1.56 4.15 0.20 1.63 0.76 0.26 7.81 1.25 16.2 1.47 1.36 1.60];
bi4 = [2.56 10.2 2.70 5.70 0.91 2.46 2.22 2.49 10.0 2.35 20.0 2.58 2.44 1.86];

hw = A_w + B_w * sum( (ai4*sigma_w) ./ ((freqGHz - fi4).^2 + bi4) );
end
