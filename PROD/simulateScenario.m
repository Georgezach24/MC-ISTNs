function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSinrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSinrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSinrDbVec, ...
    newChannelState, networkEnergy, serviceStateVec, throughputMbpsVec, ...
    bsUnavailReasonVec, satUnavailReasonVec, policyKpis, newDecision] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, ...
                     prevChannelState, prevDecision)
% Για κάθε χρήστη: επιλέγει τον καλύτερο κόμβο (BS ή δορυφόρο) βάσει SINR και
% υπολογίζει χωρητικότητα/ενέργεια μετά την κατανομή εύρους ζώνης.
% Επιστρέφει και per-candidate διαγνωστικά (καλύτερος BS + καλύτερος δορυφόρος,
% ανεξάρτητα από την τελική επιλογή) για την εκπαίδευση του ML μοντέλου του Part 2.
% Το sat_geo έχει μία γραμμή [lat lon alt] ανά δορυφόρο.
%
% 23η έξοδος policyKpis: η ίδια πραγματοποίηση καναλιού αποτιμάται με όλες τις
% πολιτικές επιλογής κόμβου (η υπό εξέταση, αποκλειστικά επίγεια, αποκλειστικά
% δορυφορική, και οι δύο πρώτες με εξισορρόπηση φόρτου, *MLB), με την κατανομή
% πόρων να υπολογίζεται εκ νέου σε καθεμία. Τα υποψήφια SINR έχουν ήδη
% υπολογιστεί για κάθε χρήστη ανεξάρτητα από το ποιος κερδίζει, οπότε οι
% πολιτικές διαφέρουν μόνο στην επιλογή και η σύγκριση είναι κατά ζεύγη πάνω
% στο ίδιο κανάλι.
%
% 8ο όρισμα prevDecision / 24η έξοδος newDecision: κατάσταση της απόφασης
% (εξυπηρετών κόμβος, υποψήφιος σε εκκρεμότητα, μετρητής επιβεβαίωσης) ανά
% πολιτική. Χρειάζεται επειδή η υστέρηση και ο χρόνος επιβεβαίωσης κάνουν την
% απόφαση εξαρτώμενη από το παρελθόν. Το runSimulation.m τη μεταφέρει μέσα σε
% κάθε διέλευση και τη μηδενίζει στην αρχή της επόμενης.
%
% prevChannelState (προαιρετικό, 7ο όρισμα): αν δοθεί, LOS/NLOS και shadow
% fading κάθε ζεύξης συσχετίζονται χωρικά με την προηγούμενη κλήση αντί για
% i.i.d. δειγματοληψία. Το περνάει το runSimulation.m μέσα σε κάθε διέλευση
% (διαδοχικά βήματα 1 s με κινούμενους χρήστες) και το μηδενίζει στην αρχή
% κάθε νέας διέλευσης, ώστε οι διελεύσεις να είναι ανεξάρτητες μεταξύ τους.
% Τα scripts επαλήθευσης το παραλείπουν.
if nargin < 7
    prevChannelState = [];
end
if nargin < 8
    prevDecision = [];
end

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);   % κελιά (TRxP): κόμβοι "BS1".."BSnumBs"
% Τομείς: με simParameters.CellSite (site κάθε κελιού) και
% simParameters.CellBoresightDeg (κατεύθυνση τομέα, αζιμούθιο GCS από τον
% άξονα x = ανατολή, αριστερόστροφα) κάθε site έχει πολλά κελιά στην ίδια θέση,
% με κατευθυντική κεραία (TR 38.901 Πίν. 7.8-1: 3 τομείς ανά site, 30°, 150°,
% 270°). Χωρίς αυτά, ένα κελί ανά σταθμό με σταθερό κέρδος (simParameters.EIRP).
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
% Μία γραμμή του sat_geo ανά δορυφόρο, κόμβος "SAT-k".
numSats  = size(sat_geo,1);

% Επαναχρησιμοποίηση συχνότητας μεταξύ δεσμών (TR 38.821 Πίν. 6.1.1.1-5:
% Option 1 = 1, Option 2 = 3). Κάθε δορυφόρος φωτίζει την περιοχή με μία
% δέσμη· ο δορυφόρος k εκπέμπει στη ζώνη mod(k-1, FRF), και το εύρος ανά δέσμη
% (satParameters.Bandwidth) έχει ήδη προσαρμοστεί στο FRF από τον καλούντα.
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

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SINR (κατάσταση outage) ------------------
% Κάτω από αυτό -> ο χρήστης θεωρείται outage αντί να ανατεθεί στον
% "λιγότερο κακό" κόμβο. = Shannon-ισοδύναμο SINR του MCS 0 (TS 38.214 Πίν. 5.1.3.1-1).
minSpectralEfficiency = 0.2344;                        % bits/s/Hz (MCS 0)
minUsableSinrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Κατώφλι σε επίπεδο υπηρεσίας ------------------
% 5ο εκατοστημόριο φασματικής απόδοσης χρήστη (απαίτηση ITU-R M.2410, όπως
% παρατίθεται στο TR 37.910). Η τιμή εξαρτάται από το περιβάλλον:
% 0.3 bit/s/Hz για Indoor Hotspot-eMBB DL (Πίν. 5.4.1.1.1-1) και
% 0.225 bit/s/Hz για Dense Urban-eMBB DL (Πίν. 5.4.1.2.1-1). Εφαρμόζεται η
% αυστηρότερη τιμή. Ορίζεται επί του ΣΥΝΟΛΙΚΟΥ εύρους καναλιού, άρα
% SE_ζεύξης >= targetNormalizedSe * L.
targetNormalizedSe = 0.3;   % bit/s/Hz

%% ------------------ Ισχύς εκπομπής ανά αλυσίδα πομποδέκτη ------------------
% Το μοντέλο EARTH (Auer et al. 2011, εξ. 1) ορίζει P_out ΑΝΑ αλυσίδα, με
% ανώτατο όριο P_max = 20 W για μακροκυψελικό σταθμό (Πίν. 2).
pOutTotalW    = 10^((simParameters.TxPower - 30)/10);
pOutPerChainW = pOutTotalW / simParameters.Power.NumTrx;
earthPmaxW    = 20;
if pOutPerChainW > earthPmaxW
    warning('simulateScenario:EarthOutOfRange', ...
        ['Ισχύς ανά αλυσίδα %.1f W > P_max = %.0f W του μοντέλου EARTH ' ...
         '(Auer et al. 2011, Πίν. 2). Αύξησε το simParameters.Power.NumTrx.'], ...
        pOutPerChainW, earthPmaxW);
end

%% ------------------ Διαλείψεις δορυφορικής ζεύξης ------------------
% Shadowed Rician (Abdi et al. 2003). Η σκίαση περιέχεται ήδη στο μοντέλο
% (τυχαίο πλάτος LOS κατά Nakagami-m), οπότε δεν προστίθεται χωριστός
% λογαριθμοκανονικός όρος. Παράμετροι (b0,m,Ω) από την ανύψωση, εξ. (19).
satFadeElevRangeDeg = [20 80];   % πεδίο ισχύος της προσαρμογής της εξ. (19)
% Απόσταση συσχέτισης της σκίασης στο οριζόντιο επίπεδο: TR 38.811 Πίν.
% 6.7.2-3a (Urban, S-band), SF 37 m για LOS. Το Shadowed Rician έχει πάντα
% συνιστώσα LOS (σκιασμένη κατά Nakagami-m), οπότε ισχύει η τιμή LOS.
satShadowDcorrM     = 37;
terrKdBLos          = 9;         % Rician K επίγειο LOS, TR 38.901 Πίν. 7.5-6

%% ------------------ Υστέρηση και χρόνος επιβεβαίωσης ------------------
% Μεταφορά του Event A3 (TS 38.331 §5.5.4.4) στον κανόνα της εργασίας. Το
% Hysteresis ορίζεται στην §6.3.2 ως INTEGER (0..30) με τιμή = πεδίο * 0.5 dB,
% δηλαδή 0 έως 15 dB· το TimeToTrigger ως απαριθμημένο σύνολο τιμών σε ms.
% Προεπιλογή 0 και 0: αναπαράγει τη συμπεριφορά χωρίς υστέρηση.
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
% Ping-pong: παραμονή ΜΙΚΡΟΤΕΡΗ από το MTS (TR 36.839 Ορισμός 5), σε ακέραια
% βήματα. Με Δt = 1 s η μικρότερη παρατηρήσιμη παραμονή είναι 1 βήμα = MTS,
% οπότε κρατείται το 1 βήμα· με Δt = 0.2 s: έως 4 βήματα = 0.8 s.
ppMaxStaySteps = max(1, ceil(pingPongMtsS/dtSecHys - 1e-9) - 1);
% Ο χρόνος επιβεβαίωσης μετατρέπεται σε διαδοχικά βήματα. Η στρογγυλοποίηση
% είναι προς τα πάνω επειδή το πρότυπο απαιτεί η συνθήκη να ΙΣΧΥΕΙ για τη
% διάρκεια αυτή· με Δt = 1 s οι τιμές κάτω από 1024 ms δεν διακρίνονται.
nTtt = max(1, ceil(tttMs / (dtSecHys*1000)));

%% ------------------ Παρακολούθηση ζεύξης και RLF ------------------
% Ο εξυπηρετών κόμβος δεν εγκαταλείπεται με την πρώτη βύθιση του SINR. Η
% απώλεια της ζεύξης κρίνεται με τη διαδικασία του TS 38.331 §5.3.10:
% N310 διαδοχικές ενδείξεις out-of-sync (SINR < Qout) ξεκινούν τον T310·
% N311 διαδοχικές in-sync (SINR > Qin) όσο τρέχει τον σταματούν· αν λήξει,
% δηλώνεται RLF και ο χρήστης κάνει επανεγκατάσταση (§5.3.7) σε κελί που
% επιλέγει εκ νέου. Τιμές: TR 36.839 Πίν. 5.2.1.3.1 (Qout -8 dB, Qin -6 dB,
% T310 1 s = προεπιλογή του 36.331, N310 = N311 = 1). Μία ένδειξη ανά βήμα:
% με N310 = N311 = 1 η καταμέτρηση είναι η ίδια· ο T310 γίνεται βήματα με
% στρογγυλοποίηση προς τα πάνω, όπως ο χρόνος επιβεβαίωσης.
rlfP = struct('QoutDb', -8, 'QinDb', -6, 'T310Ms', 1000, 'N310', 1, 'N311', 1);
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'Rlf')
    rlfP = simParameters.Mobility.Rlf;
end
mob = struct('MinUsableDb', minUsableSinrDb, 'HysDb', hysDb, 'NTtt', nTtt, ...
    'QoutDb', rlfP.QoutDb, 'QinDb', rlfP.QinDb, 'N310', rlfP.N310, 'N311', rlfP.N311, ...
    'NT310', max(1, ceil(rlfP.T310Ms / (dtSecHys*1000))));

%% ------------------ Εξισορρόπηση φόρτου (Mobility Load Balancing) ------------------
% TS 38.300 §15.5.1: το MLB μοιράζει τον φόρτο αλλάζοντας τις παραμέτρους του
% handover, με την αλλαγή εκφρασμένη ως διαφορά (delta) στο κατώφλι ενεργοποίησης
% (§15.5.1.4) και πάντα μέσα στα όρια που επιτρέπει το OAM. Στον τύπο του A3
% (TS 38.331 §5.5.4.4) αυτό είναι η μετατόπιση ανά κυψέλη (Ocn/Ocp). Το πρότυπο
% δεν ορίζει ΠΟΣΟ αλλάζει η μετατόπιση· ο κανόνας είναι του Nasri & Altman
% (arXiv:1307.1212, εξ. 4-5): το περιθώριο από την κυψέλη e προς την k είναι
%     HM(e,k) = f(x),  f(x) = f(0) + (f(0) - HMmax)*x,  x = χe - χk ∈ [-1, 1]
% με f(0) = 6 dB και HMmax = 12 dB, δηλαδή το περιθώριο μετατοπίζεται γύρω από
% το προγραμματισμένο κατά -6 dB ανά μονάδα διαφοράς φόρτου. Η αντισυμμετρία
% f(x) + f(-x) = 2 f(0) είναι ο όρος του θεωρήματός τους που αποτρέπει το
% ping-pong: η γεμάτη κυψέλη σπρώχνει χρήστες και η άδεια καθυστερεί την
% επιστροφή τους. Εδώ το προγραμματισμένο περιθώριο είναι η υστέρηση του A3,
% οπότε μεταφέρεται μόνο ο όρος του φόρτου: μετατόπιση υπέρ του στόχου ίση με
% MlbMaxOffsetDb * (χ_εξυπηρετούντος - χ_στόχου). Η μετατόπιση είναι το Ocn
% που ο εξυπηρετών κόμβος δίνει στον γείτονα, οπότε στρογγυλεύεται στις
% τιμές που μπορεί να πάρει το πεδίο (Q-OffsetRange, TS 38.331).
if isfield(simParameters, 'Mobility') && isfield(simParameters.Mobility, 'MlbMaxOffsetDb')
    mlbMaxOffsetDb = simParameters.Mobility.MlbMaxOffsetDb;
else
    mlbMaxOffsetDb = 6;   % = HMmax - f(0) του Nasri & Altman
end

%% ------------------ Μέτρηση: μέση τιμή και φίλτρο L3 ------------------
% Το TS 38.331 §5.5.3.2 ορίζει ότι ο A3 δεν συγκρίνει τη στιγμιαία μέτρηση
% αλλά τη φιλτραρισμένη, F_n = (1-a)*F_{n-1} + a*M_n ("filter the measured
% result, BEFORE using for evaluation of reporting criteria"), και το ίδιο το
% M_n δεν είναι στιγμιαίο δείγμα: είναι αποτέλεσμα περιόδου μέτρησης
% τουλάχιστον 200 ms στο FR1 (TS 38.133 Πίν. 9.2.5.2-1). Αντίγραφα των δύο
% εδαφίων: Sources/TS38331_L3_filtering.md.
%
%   Meas.NAvg    πλήθος ανεξάρτητων πραγματώσεων γρήγορων διαλείψεων που
%                μεσοποιούνται μέσα στο βήμα. Προκύπτει από την αυτοσυσχέτιση
%                του Clarke και υπολογίζεται στο runSimulation.m.
%   Meas.L3Alpha ο συντελεστής a, ΗΔΗ προσαρμοσμένος στο βήμα της προσομοίωσης
%                κατά τη ρήτρα "adapt the filter such that the time
%                characteristics ... are preserved at different input rates".
%                Τιμή 1 = χωρίς φίλτρο (NOTE 1, k = 0).
%
% Προεπιλογές 1 και 1: αναπαράγουν ακριβώς τη συμπεριφορά πριν τη μέτρηση,
% ώστε τα scripts επαλήθευσης και η παλιά ρύθμιση να μένουν αναπαραγώγιμα.
measNAvgTerr = max(1, round(measParam(simParameters, 'NAvg',    1)));
measNAvgSat  = max(1, round(measParam(satParameters, 'NAvg',    1)));
l3Alpha      = min(max(measParam(simParameters, 'L3Alpha', 1), 0), 1);

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestBsNodeVec       = strings(numUsers,1);   % ποιος BS κέρδισε, για τη σύγκριση πολιτικών
% Οι τρεις βασικές πολιτικές και οι δύο εκδοχές τους με MLB: ίδιοι υποψήφιοι
% και ίδιος A3, με επιπλέον τη μετατόπιση φόρτου.
policyNamesLocal    = {'Actual','TerrestrialOnly','SatelliteOnly', ...
                       'ActualMLB','TerrestrialOnlyMLB'};
policyBaseLocal     = {'Actual','TerrestrialOnly','SatelliteOnly', ...
                       'Actual','TerrestrialOnly'};
policyMlbLocal      = [false false false true true];
numPolLocal         = numel(policyNamesLocal);
% Prev/Stay: ο κόμβος πριν το τελευταίο handover και τα βήματα παραμονής στον
% τρέχοντα από τότε. Χρειάζονται μόνο για την καταμέτρηση των ping-pong· δεν
% επηρεάζουν την απόφαση.
newDecision = struct();
% LoadNode/LoadRho: εικονικός φόρτος κάθε κόμβου στο τέλος του βήματος, όπως
% τον αναφέρουν οι σταθμοί ο ένας στον άλλον (TS 38.300 §15.5.1, load
% reporting). Η απόφαση του επόμενου βήματος βλέπει τον φόρτο αυτού, δηλαδή με
% καθυστέρηση ενός βήματος, και όχι τον φόρτο που θα προκύψει από την ίδια.
% Υπολογίζεται σε όλες τις πολιτικές, αλλά μόνο οι *MLB τον χρησιμοποιούν.
% Cause/Offset: πώς άλλαξε ο κόμβος στο βήμα ("HO" επιτυχές handover του A3,
% "RLF" επανεγκατάσταση μετά από RLF, "HOF" επανεγκατάσταση μετά από αποτυχία
% στον στόχο, "" καμία αλλαγή) και ποια μετατόπιση είχε ο επιλεγμένος κόμβος.
% Oos/Is/T310: μετρητές N310/N311 και βήματα από την εκκίνηση του T310 (-1 =
% δεν τρέχει), TS 38.331 §5.3.10.1.
for pol = policyNamesLocal
    newDecision.(pol{1}) = struct('Node', strings(numUsers,1), ...
        'Pending', strings(numUsers,1), 'Count', zeros(numUsers,1), ...
        'Prev', strings(numUsers,1), 'Stay', zeros(numUsers,1), ...
        'Cause', strings(numUsers,1), 'Offset', zeros(numUsers,1), ...
        'Oos', zeros(numUsers,1), 'Is', zeros(numUsers,1), 'T310', -ones(numUsers,1), ...
        'LoadNode', strings(0,1), 'LoadRho', zeros(0,1), 'LoadAge', 0);
end
hasPrevDecision = ~isempty(prevDecision);

% Επιτυχή handovers του βήματος ανά πολιτική, ανά τύπο: επίγειο-επίγειο,
% επίγειο-δορυφορικό (και προς τις δύο κατευθύνσεις), δορυφορικό-δορυφορικό.
% RLF ανά τύπο εξυπηρετούντος κόμβου (TR 36.839 Ορισμός 2) και αποτυχίες
% handover (Ορισμός 3).
hoCount = struct();
for pol = policyNamesLocal
    hoCount.(pol{1}) = struct('TerrTerr', 0, 'TerrSat', 0, 'SatSat', 0, 'PingPong', 0, ...
        'RlfTerr', 0, 'RlfSat', 0, 'HoFail', 0);
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

% Αναθέσεις ανά πολιτική, μία στήλη ανά policyNamesLocal. Η στήλη 1 (Actual)
% συμπληρώνεται κι αυτή, για τον υπολογισμό του φόρτου της.
polNodeMat = strings(numUsers,numPolLocal);
polTypeMat = strings(numUsers,numPolLocal);
polSinrMat = nan(numUsers,numPolLocal);

% Διαγνωστικοί πίνακες. Η διάδοση (απόσταση, LOS, σκίαση, απώλειες) ορίζεται
% ανά site, επειδή εξαρτάται μόνο από τις θέσεις και τα κελιά ενός site είναι
% στο ίδιο σημείο· ανά κελί διαφέρει μόνο το κέρδος της κεραίας.
pathLossMat       = nan(numUsers,numSites);
pLosMat           = nan(numUsers,numSites);
losMat            = false(numUsers,numSites);
losLatentMat      = zeros(numUsers,numSites);
sfMat             = nan(numUsers,numSites);
% Soft LOS (TR 38.901 §7.6.3.3, εξ. 7.6-18/19): δύο διεργασίες σκίασης, μία
% για κάθε κατάσταση, ώστε το κανάλι να υπολογίζεται και ως LOS και ως NLOS.
sfLosMat          = zeros(numUsers,numSites);
sfNlosMat         = zeros(numUsers,numSites);
losSoftMat        = nan(numUsers,numSites);
useSoftLos = isfield(simParameters, 'SoftLos') && simParameters.SoftLos;
% Συντελεστής της εξ. 7.6-18: sqrt(d_LOS/λ). Το d_LOS δεν ορίζεται ρητά στο
% εδάφιο· λαμβάνεται η απόσταση συσχέτισης της κατάστασης LOS/NLOS του Πίν.
% 7.6.3.1-2 (UMa/UMi: 50 m), με την οποία παράγεται και η G.
softLosK = sqrt(50 / (physconst('LightSpeed')/simParameters.CarrierFrequency));
snrDbMat          = nan(numUsers,numBs);   % χωρίς παρεμβολή, διαγνωστικό
sinrDbMat         = nan(numUsers,numBs);
cellPathLossMat   = nan(numUsers,numBs);   % = απώλειες του site του κελιού
cellRange3DMat    = nan(numUsers,numBs);
% Διαγνωστικά του καλύτερου δορυφόρου ανά χρήστη (ένας ανά γραμμή του
% dataset, όπως και για τους BS), και το SINR κάθε δορυφόρου χωριστά.
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
% Κρυφή γκαουσιανή μεταβλητή της δορυφορικής σκίασης ανά (χρήστη, δορυφόρο)
% (βλ. «Satellite candidates»). Ταξιδεύει μέσα στο prevChannelState και
% μηδενίζεται στην αρχή κάθε διέλευσης μαζί με το υπόλοιπο κανάλι.
hasPrevSatShadow = hasPrevState && isfield(prevChannelState, 'SatShadowLatent') && ...
    isequal(size(prevChannelState.SatShadowLatent), [numUsers, numSats]);
satShadowLatentMat = zeros(numUsers,numSats);

% Κατάσταση του φίλτρου L3 ανά υποψήφιο. NaN σημαίνει "δεν υπάρχει
% προηγούμενη μέτρηση", οπότε F_0 = M_1 όπως ορίζει η §5.5.3.2. Ταξιδεύει
% μέσα στο prevChannelState, δηλαδή μηδενίζεται στην αρχή κάθε διέλευσης μαζί
% με το υπόλοιπο κανάλι. Η μέτρηση είναι ΚΟΙΝΗ για τις τρεις πολιτικές: δεν
% εξαρτάται από το ποιος κόμβος εξυπηρετεί, άρα η σύγκριση κατά ζεύγη μένει
% έγκυρη.
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

lambdaSat    = physconst('LightSpeed') / satParameters.CarrierFrequency;
noiseSatLinW = 10^(noisePowerSAT_dBW/10);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SINR) ------------------
useWrapAround = isfield(simParameters, 'BsImageGeo') && ~isempty(simParameters.BsImageGeo);
if useWrapAround
    bsImageGeo = simParameters.BsImageGeo;
end

% Γεωμετρία και απώλειες διάδοσης όλων των ζεύξεων, πριν τον βρόχο των
% χρηστών. Είναι ντετερμινιστικές, οπότε υπολογίζονται με μία διανυσματική
% κλήση αντί για μία ανά ζεύξη· οι συναρτήσεις (distance, nrPathLoss,
% geodetic2aer) είναι στοιχείο-προς-στοιχείο και δίνουν τα ίδια ακριβώς bits.
% Μέσα στον βρόχο μένουν μόνο οι κληρώσεις, με την ίδια σειρά όπως πριν.
uLat  = user_geo(:,1);
uLon  = user_geo(:,2);
hUtV  = user_geo(:,3);
% Wrap-around: με simParameters.BsImageGeo (numSites x K x 2, [lat lon] του
% site και των K-1 μετατοπισμένων αντιγράφων του) η απόσταση κάθε ζεύξης
% μετριέται από το πλησιέστερο αντίγραφο. Όλο το κανάλι (απώλειες, LOS,
% σκίαση, παρεμβολή) εξαρτάται από το site μόνο μέσω αυτής.
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
% Οριζόντια απόσταση (γεωδαιτική) και φυσικά ύψη κεραιών. Η γεωμετρία για τη
% nrPathLoss κατασκευάζεται απευθείας από αυτά: η κατακόρυφη συντεταγμένη ENU
% αποκλίνει από το ύψος κεραίας λόγω καμπυλότητας.
groundDistanceMat = gdAll;
range3DMat        = hypot(gdAll, siteGeo(:,3)' - hUtV);
% Πεδίο ισχύος UMa/UMi (TR 38.901 Πίν. 7.4.1-1). Εκτός ορίων η ζεύξη δεν
% υπολογίζεται: το SINR μένει NaN, ώστε να ξεχωρίζει από ζεύξη που
% υπολογίστηκε και βρέθηκε ανεπαρκής.
validLinkMat = isValidTerrestrialLink(gdAll, repmat(hUtV, 1, numSites), simParameters.CarrierFrequency);
% Πιθανότητα LOS κάθε ζεύξης (TR 38.901 §7.4.2).
pLosAll = losProbability38901(gdAll, repmat(hUtV, 1, numSites), simParameters.PathLoss.Scenario);
% Απώλειες LOS και NLOS κάθε έγκυρης ζεύξης· η κατάσταση LOS κληρώνεται στον
% βρόχο και διαλέγει μία από τις δύο.
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
% Αζιμούθιο από το site (ή το πλησιέστερο αντίγραφό του) προς τον χρήστη, από
% βορρά δεξιόστροφα, για την κεραία τομέα.
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
% Απόσταση μετακίνησης κάθε χρήστη από την προηγούμενη κλήση· καθορίζει πόσο
% συσχετίζεται η σκίαση με την προηγούμενη τιμή.
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
% MLB: φόρτος κάθε κόμβου (κελιά, μετά δορυφόροι) από την αναφορά του
% προηγούμενου βήματος, και ο δείκτης του εξυπηρετούντος κάθε χρήστη. Κόμβος
% που δεν εμφανίζεται στην αναφορά δεν εξυπηρετούσε κανέναν: φόρτος 0.
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
    % Αρχικοποιήσεις
    userBestSinr      = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

    userMoveDistance = userMoveAll(u);

    %% ===== Terrestrial BS candidates =====
    % Κληρώσεις ανά ζεύξη (LOS, σκίαση, γρήγορες διαλείψεις) με τη σειρά
    % χρήστης -> site, όπως πάντα.
    for b = 1:numSites
        if ~validLinkMat(u,b)
            sfMat(u,b) = 0;
            continue;
        end

        % LOS ανά ζεύξη βάσει πιθανότητας απόστασης (TR 38.901 §7.4.2).
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
            % Soft LOS, TR 38.901 §7.6.3.3: LOS_soft = 1/2 + (1/π) arctan(
            % sqrt(d_LOS/λ) (G + F(d))), F(d) = sqrt(2) erf^-1(2 P_LOS(d) - 1).
            % G = -latent: ίδια κατανομή, και το πρόσημο του G + F δίνει την
            % ίδια σκληρή απόφαση με το spatiallyConsistentLos (όριο K -> ∞).
            % Τα βήματα του καναλιού γίνονται δύο φορές, μία ως LOS και μία ως
            % NLOS, και H = sqrt(LOS_soft) H_LOS + sqrt(1-LOS_soft) H_NLOS
            % (εξ. 7.6-19). Σε ισχύ, με ανεξάρτητες φάσεις των δύο σκελών:
            % P = LOS_soft P_LOS + (1 - LOS_soft) P_NLOS.
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

        % Shadow fading: log-normal δείγμα με τυπική απόκλιση sigmaSF (TR 38.901 §7.4.1),
        % AR(1)-συσχετισμένο με το προηγούμενο όταν useCorrelatedSF, αλλιώς ανεξάρτητο.
        if useCorrelatedSF
            sfSample = rho*prevSF + sqrt(1 - rho^2) * sigmaSF * randn();
        else
            sfSample = sigmaSF * randn();
        end
        sfMat(u,b) = sfSample;
        pathLoss = pathLoss + sfSample;

        % Small-scale fading: Rician (LOS, K=terrKdBLos) ή Rayleigh (NLOS).
        % Ο χρόνος συνοχής (~43 ms στα 3.5 GHz με πεζό χρήστη) είναι πολύ
        % μικρότερος από το βήμα, οπότε μέσα σε ένα βήμα χωρούν πολλές
        % ανεξάρτητες πραγματώσεις. Η μέτρηση είναι η μέση τιμή τους σε
        % γραμμική ισχύ, όχι μία από αυτές (TS 38.133 Πίν. 9.2.5.2-1).
        pathLoss = pathLoss - smallScaleFadingDb(isLos, terrKdBLos, measNAvgTerr);
        pathLossMat(u,b) = pathLoss;
    end

    %% ===== Λαμβανόμενη ισχύς ανά κελί =====
    % Χωρίς τομείς: σταθερό σύνθετο κέρδος, EIRP = P + G_E,max + 10log10(M).
    % Με τομείς: κέρδος της κεραίας του κελιού στην κατεύθυνση του χρήστη,
    % G(θ,φ) = G_E,max + A''(θ,φ) + κέρδος συστοιχίας (βλ. bsSectorGainDb).
    % Η κατεύθυνση είναι η γεωμετρική (LOS) από το site, ή από το πλησιέστερο
    % αντίγραφό του με wrap-around, όπως στη βαθμονόμηση μεγάλης κλίμακας του
    % TR 38.901 Πίν. 7.8-1 ("UT attachment: based on pathloss considering the
    % LOS angle") και του TR 36.873 Πίν. 8.2-1. Τα κελιά του ίδιου site έχουν
    % τις ίδιες απώλειες διάδοσης· ισχύς ανά κελί σε dBW.
    cellPathLossMat(u,:) = pathLossMat(u, cellSite);
    cellRange3DMat(u,:)  = range3DMat(u, cellSite);
    if sectorized
        % Αζιμούθιο από βορρά, δεξιόστροφα -> αζιμούθιο GCS φ (από ανατολή,
        % αριστερόστροφα).
        phiSiteDeg   = 90 - azNorthAll(u,:)';
        % Ζενιθιακή γωνία στο site (90° = οριζόντια, > 90° προς τα κάτω).
        thetaSiteDeg = 90 + atand((siteGeo(:,3) - user_geo(u,3)) ./ groundDistanceMat(u,:)');
        phiCell = mod(phiSiteDeg(cellSite) - cellBoresightDeg + 180, 360) - 180;
        gainDb  = bsSectorGainDb(thetaSiteDeg(cellSite), phiCell, simParameters.BsAntenna);
        rxDbW   = (simParameters.TxPower - 30) + gainDb' - cellPathLossMat(u,:);
    else
        rxDbW   = (simParameters.EIRP - 30) - cellPathLossMat(u,:);
    end
    % Ο λόγος σήματος προς θόρυβο ΜΟΝΟ, χωρίς παρεμβολή, κρατείται ως
    % διαγνωστικό: είναι το άνω φράγμα του SINR και δείχνει πόσο κοστίζει
    % η παρεμβολή σε κάθε ζεύξη.
    snrDbMat(u,:) = rxDbW - noisePowerBS_dBW;

    %% ===== SINR ανά επίγειο υποψήφιο =====
    % Όλοι οι σταθμοί εκπέμπουν ταυτόχρονα στο ίδιο φάσμα: επαναχρησιμοποίηση
    % συχνοτήτων 1 και μοντέλο πλήρους απασχόλησης, όπως ορίζει ο ΠΙΝΑΚΑΣ 5 β)
    % του ITU-R M.2412-0 (§8.4, "Inter-site interference modeling: Explicitly
    % modelled", "Traffic model: Full buffer") και όπως προϋποθέτουν οι μετρικές
    % βαθμονόμησης του TR 38.901 §7.8 (Πίν. 7.8-1 "Geometry", Πίν. 7.8-2
    % "Wideband SIR", Πίν. 7.8-3 "Wideband SINR").
    %
    % Ορισμός: TR 38.821 Πίν. 6.1.1.2-1, ΣΗΜΕΙΩΣΗ:
    %   Geometry SINR = -10*log10(I/C + N/C)  <=>  SINR = C / (I + N)
    % με C, I, N μετρημένα στο ίδιο εύρος ζώνης. Επειδή σήμα, παρεμβολή και
    % θόρυβος κλιμακώνονται όλα με το εκχωρημένο εύρος, ο λόγος δεν εξαρτάται
    % από την κατανομή πόρων: υπολογίζεται μία φορά στο εύρος του κόμβου.
    %
    % Σταθμός εκτός του πεδίου ισχύος του μοντέλου δεν προσμετράται ως
    % παρεμβολέας, γιατί οι απώλειές του δεν είναι υπολογίσιμες με το UMa/UMi
    % χωρίς να παραβιαστεί ο έλεγχος εγκυρότητας. Η παρεμβολή είναι επομένως
    % κάτω φράγμα, όπως και λόγω του μικρού πλήθους σταθμών.
    noiseLinW   = 10^(noisePowerBS_dBW/10);
    rxPowLinW   = 10.^(rxDbW/10);
    rxPowLinW(~isfinite(rxPowLinW)) = 0;
    totalRxLinW = sum(rxPowLinW);

    vb = isfinite(cellPathLossMat(u,:));
    sinrRow = 10*log10(rxPowLinW(vb) ./ (noiseLinW + (totalRxLinW - rxPowLinW(vb))));
    sinrDbMat(u,vb) = sinrRow;
    % Φίλτρο L3 (βλ. l3Filter), στοιχείο προς στοιχείο.
    fPrev = prevFiltBs(u,vb);
    fNew  = sinrRow;
    hasF  = isfinite(fPrev);
    fNew(hasF) = (1 - l3Alpha)*fPrev(hasF) + l3Alpha*sinrRow(hasF);
    filtSinrBsMat(u,vb) = fNew;

    % Η επιλογή είναι μονότονη ως προς τη λαμβανόμενη ισχύ: με κοινό
    % άθροισμα ισχύων, το SINR αυξάνει με το C, οπότε ο καλύτερος σταθμός
    % είναι ο ίδιος με ή χωρίς παρεμβολή. Αλλάζει η τιμή, όχι η επιλογή.
    % Ο πρώτος με το μέγιστο SINR, όπως η σειριακή σύγκριση με «>».
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

    % Στιγμιότυπο του καλύτερου υποψήφιου BS πριν τη σύγκριση με τον δορυφόρο (per-candidate διαγνωστικό).
    % NaN όταν καμία επίγεια ζεύξη δεν ήταν εντός του πεδίου ισχύος του μοντέλου.
    bestBsNodeVec(u) = userBestNode;
    if isfinite(userBestSinr)
        bestBsSinrDbVec(u) = userBestSinr;
        if userBestSinr < minUsableSinrDb
            bsUnavailReasonVec(u) = "BelowSinrFloor";
        end
    else
        % Καμία επίγεια ζεύξη δεν ήταν εντός του πεδίου ισχύος του μοντέλου.
        bsUnavailReasonVec(u) = "OutOfModelRange";
    end
    bestBsDistanceVec(u) = userBestDistance;
    bestBsPathLossVec(u) = userBestPathLoss;

    %% ===== Satellite candidates =====
    % Κάθε δορυφόρος είναι χωριστός υποψήφιος. Πρώτα υπολογίζονται οι απώλειες
    % όλων, επειδή το SINR του καθενός εξαρτάται από την ισχύ των υπολοίπων.
    satPlRow    = inf(1,numSats);
    satElevRow  = nan(1,numSats);
    satRangeRow = nan(1,numSats);
    for k = 1:numSats
        % Σκίαση με χρονική συσχέτιση: η μεταβλητή εξελίσσεται σε ΚΑΘΕ βήμα,
        % ακόμη και όταν ο δορυφόρος είναι κάτω από τη μάσκα, ώστε όταν
        % ξαναφανεί να συνεχίζει από τη σωστή κατάσταση και όχι από νέα
        % κλήρωση. Μία μεταβλητή ανά δορυφόρο: κάθε δορυφόρος φαίνεται από
        % άλλη κατεύθυνση, άρα τα εμπόδια γύρω από τον χρήστη που τον σκιάζουν
        % είναι άλλα. Το TR 38.811 δεν δίνει συσχέτιση μεταξύ δορυφόρων, οπότε
        % οι μεταβλητές είναι ανεξάρτητες μεταξύ τους.
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

        % Ατμοσφαιρική απόσβεση αερίων (οξυγόνο + υδρατμοί), TR 38.811 §6.6.4.
        % Βροχή/νέφωση & ιονοσφαιρική σπινθηρίδα δεν μοντελοποιούνται.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        % Shadowed Rician. Η σκίαση (πλάτος LOS κατά Nakagami-m) συσχετίζεται
        % με τη μετατόπιση του χρήστη, με εκθετική αυτοσυσχέτιση και απόσταση
        % satShadowDcorrM (TR 38.811 Πίν. 6.7.2-3a), όπως και η επίγεια
        % σκίαση. Η κρυφή γκαουσιανή μεταβλητή μετατρέπεται στην κατανομή του
        % μοντέλου μέσω της αθροιστικής της (copula), οπότε η κατανομή μένει
        % ακριβώς η Gamma του Abdi και αλλάζει ΜΟΝΟ η συσχέτιση στον χρόνο.
        % Όταν αλλάζει η ανύψωση αλλάζουν οι παράμετροι, αλλά το εκατοστημόριο
        % της σκίασης (πόσο «σκιασμένος» είναι ο χρήστης) συνεχίζει ομαλά.
        % Δεν μοντελοποιείται η αποσυσχέτιση από την κίνηση του δορυφόρου· το
        % πρότυπο δίνει μόνο την οριζόντια απόσταση, και η NOTE 2 του ίδιου
        % εδαφίου σημειώνει ότι στην πράξη οι αποστάσεις μπορεί να είναι
        % μικρότερες.
        % Η ανύψωση περιορίζεται στο πεδίο ισχύος της προσαρμογής· εναλλακτικά
        % σταθερή κατάσταση σκίασης μέσω satParameters.ShadowingState (για
        % ανάλυση ευαισθησίας).
        if isfield(satParameters, 'ShadowingState') && ~isempty(satParameters.ShadowingState)
            [b0, mNak, omega] = shadowedRicianStateParams(satParameters.ShadowingState);
        else
            elevClamped = min(max(elevSat, satFadeElevRangeDeg(1)), satFadeElevRangeDeg(2));
            [b0, mNak, omega] = shadowedRicianElevParams(elevClamped);
        end
        uShadow = 0.5*erfc(-satShadowLatentMat(u,k)/sqrt(2));   % γκαουσιανή -> ομοιόμορφη
        satPlRow(k) = satPathLoss - shadowedRicianFadingDb(b0, mNak, omega, measNAvgSat, uShadow);
    end

    % SINR ανά δορυφόρο. Καμία ομοδιαυλική επικάλυψη με το επίγειο σκέλος
    % (2.0 GHz έναντι 3.5 GHz, TR 38.821 Πίν. 6.1.3.2-1 και TR 38.901). Μεταξύ
    % δορυφόρων ισχύει ο κανόνας του TR 38.821 §6.1.1.1 (Σχ. 6.1.1.1-2): κάθε
    % άλλη δέσμη στην ίδια ζώνη είναι παρεμβολή, σε μοντέλο πλήρους
    % απασχόλησης όπως και στους BS. Παρεμβάλλουν μόνο οι δορυφόροι πάνω από τη
    % μάσκα, δηλαδή όσοι φωτίζουν την περιοχή· η ισχύς των πλευρικών λοβών
    % δεσμών που στοχεύουν αλλού δεν μοντελοποιείται, οπότε η παρεμβολή είναι
    % κάτω φράγμα. Με έναν ορατό δορυφόρο ανά ζώνη το SINR ταυτίζεται με το SNR.
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

    % Διαγνωστικά του καλύτερου δορυφόρου: ο ορατός με το μεγαλύτερο SINR, ή,
    % αν κανείς δεν είναι ορατός, αυτός με τη μεγαλύτερη ανύψωση.
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

    %% ===== Απόφαση με υστέρηση και χρόνο επιβεβαίωσης =====
    % Κατάλογος υποψηφίων: όλοι οι σταθμοί εντός πεδίου ισχύος συν οι ορατοί
    % δορυφόροι. Η υστέρηση εφαρμόζεται σε ΟΛΟΥΣ, ώστε να καλύπτει και την
    % εναλλαγή μεταξύ επίγειων σταθμών και μεταξύ δορυφόρων.
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
        st = struct('Node', "", 'Pending', "", 'Count', 0, 'Oos', 0, 'Is', 0, 'T310', -1);
        if hasPrevDecision
            pd = prevDecision.(polName);
            pN = pd.Node(u);
            pPrev = pd.Prev(u);
            pStay = pd.Stay(u);
            st = struct('Node', pN, 'Pending', pd.Pending(u), 'Count', pd.Count(u), ...
                'Oos', pd.Oos(u), 'Is', pd.Is(u), 'T310', pd.T310(u));
        end

        % Μετατόπιση φόρτου κάθε υποψηφίου ως προς τον εξυπηρετούντα (0 για
        % τον ίδιο και για τις πολιτικές χωρίς MLB). Κόμβος που δεν εμφανίζεται
        % στην αναφορά φόρτου δεν εξυπηρετούσε κανέναν: φόρτος 0.
        candOffsetDb = zeros(nCand,1);
        if policyMlbLocal(pI) && hasPrevDecision && pN ~= "" && pN ~= "None"
            % Φόρτος ανά κόμβο από την αναφορά (nodeLoadAll, βλ. πριν τον
            % βρόχο)· υποψήφιοι: τα έγκυρα κελιά με τη σειρά τους, μετά οι
            % ορατοί δορυφόροι.
            ldAll   = nodeLoadAll.(polName);
            rhoServ = ldAll(prevNodeIdx.(polName)(u));
            x = min(max(rhoServ - ldAll(candIdx), -1), 1);
            candOffsetDb = quantizeQOffsetDb(mlbMaxOffsetDb * x);
            candOffsetDb(candNode == pN) = 0;
        end

        [ci, st, ev] = mobilityStep(candNode, candSinr, candSinrMeas, keepMask, ...
            st, mob, candOffsetDb);

        % ci > 0: υποψήφιος του καταλόγου· ci = 0 με st.Node = κόμβος: ο
        % χρήστης μένει στον εξυπηρετούντα ενώ τρέχει ο T310, αλλά η ζεύξη
        % του δεν υπάρχει πια στο μοντέλο (π.χ. δορυφόρος κάτω από τη γωνία
        % αποκοπής), οπότε δεν έχει SINR· ci = 0 με "None": εκτός κάλυψης.
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
        newDecision.(polName).Oos(u)     = st.Oos;
        newDecision.(polName).Is(u)      = st.Is;
        newDecision.(polName).T310(u)    = st.T310;
        newDecision.(polName).Cause(u)   = ev.Cause;
        if ci > 0
            newDecision.(polName).Offset(u) = candOffsetDb(ci);
        end
        if ev.Rlf
            if startsWith(pN, "SAT")
                hoCount.(polName).RlfSat = hoCount.(polName).RlfSat + 1;
            else
                hoCount.(polName).RlfTerr = hoCount.(polName).RlfTerr + 1;
            end
        end
        hoCount.(polName).HoFail = hoCount.(polName).HoFail + ev.HoFail;

        % Handover = επιτυχής εκτέλεση του A3 από έναν εξυπηρετούντα κόμβο σε
        % άλλον. Η επανεγκατάσταση μετά από RLF ή αποτυχία δεν είναι handover
        % (TR 36.839 Ορισμός 6: «successful handovers excl. handover
        % failures»), ούτε η πρώτη σύνδεση ή η επανασύνδεση μετά από outage.
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
            % Ping-pong (TR 36.839 Ορισμός 5): επιστροφή στον κόμβο πριν το
            % προηγούμενο handover, με παραμονή στον pN μικρότερη από το MTS.
            if selNode == pPrev && (pStay + 1) <= ppMaxStaySteps
                hoCount.(polName).PingPong = hoCount.(polName).PingPong + 1;
            end
            newDecision.(polName).Prev(u) = pN;
            newDecision.(polName).Stay(u) = 0;
        elseif ev.Cause == "" && selNode ~= "None" && selNode == pN
            newDecision.(polName).Prev(u) = pPrev;
            newDecision.(polName).Stay(u) = pStay + 1;
        end
        % Αλλιώς (RLF, αποτυχία, outage ή πρώτη σύνδεση) Prev = "" και Stay = 0
        % από την αρχικοποίηση: η αλυσίδα ξεκινά από την αρχή, αφού «whenever
        % there is a handover failure, the time of stay should not be logged»
        % (TR 36.839 §5.2.2).

        if strcmp(polName, 'Actual')
            if selNode == "None"
                % Σε outage κρατάμε τα διαγνωστικά του καλύτερου υποψηφίου,
                % αλλά ο χρήστης δεν ανατίθεται σε κόμβο.
                bestNodeVec(u)         = "None";
                bestNodeTypeVec(u)     = "Outage";
                bestDistanceVec(u)     = userBestDistance;
                bestPathLossVec(u)     = userBestPathLoss;
                bestSinrDbVec(u)       = userBestSinr;
                bestElevationDegVec(u) = userBestElevation;
            elseif ci > 0
                % Ο επιλεγμένος κόμβος δεν είναι κατ' ανάγκη ο καλύτερος: η
                % χωρητικότητα υπολογίζεται από το SINR του κόμβου στον οποίο
                % ο χρήστης είναι όντως συνδεδεμένος.
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
% Ο κανόνας ζει σε μία μόνο συνάρτηση (allocateAndTally), επειδή καλείται και
% για τις εναλλακτικές πολιτικές παρακάτω: αν ήταν αντιγραμμένος, οι δύο
% αντιγραφές θα μπορούσαν να αποκλίνουν χωρίς να φανεί.
% Μέγιστη φασματική απόδοση του NR (MCS 28 / 64QAM, TS 38.214 Πίν.
% 5.1.3.1-1): κοινή για τη χωρητικότητα και για τον φόρτο του MLB.
alloc = struct('BW_bs', BW_bs, 'pOutPerChainW', pOutPerChainW, ...
               'targetNormalizedSe', targetNormalizedSe, 'numBs', numBs, ...
               'maxSe', 5.5547, 'minUsableDb', minUsableSinrDb);

[capacityMbpsVec, throughputMbpsVec, serviceStateVec, nodePowerWattsVec, ...
 energyPerBitUJVec, networkEnergy] = allocateAndTally(bestNodeVec, ...
    bestNodeTypeVec, bestSinrDbVec, simParameters, satParameters, alloc);

%% ------------------ Σύγκριση πολιτικών στην ίδια πραγματοποίηση ------------------
% Τα δύο υποψήφια SINR είναι ήδη υπολογισμένα για κάθε χρήστη. Οι εναλλακτικές
% πολιτικές δεν ξαναδειγματοληπτούν τίποτα: κρατούν το ίδιο κανάλι και αλλάζουν
% μόνο ποιος υποψήφιος επιλέγεται, με την κατανομή πόρων να υπολογίζεται εκ νέου
% ώστε να αποτυπωθεί ο διαφορετικός φόρτος που προκύπτει.

% Οι εναλλακτικές πολιτικές αποφασίστηκαν παραπάνω με τον ΙΔΙΟ κανόνα
% υστέρησης, καθεμία με τη δική της κατάσταση: μια πολιτική που αλλάζει κόμβο
% δεν πρέπει να συγκρίνεται με μία που δεν αλλάζει.
% Διαθέσιμο φάσμα ανά σκέλος σε αυτό το βήμα (MHz). Στο δορυφορικό σκέλος
% μετράει κάθε διακριτή ζώνη που έχει τουλάχιστον έναν ορατό δορυφόρο: με
% FRF = 1 όλοι οι δορυφόροι μοιράζονται την ίδια ζώνη, όπως οι BS τη δική τους.
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
    [capP, thrP, ~, ~, ~, netP] = allocateAndTally(polNodeMat(:,pI), polTypeMat(:,pI), ...
        polSinrMat(:,pI), simParameters, satParameters, alloc);
    policyKpis.(polName) = policySummary(capP, thrP, netP, numUsers, specMHz, ...
        hoCount.(polName), dtSecHys);
end

% Αναφορά φόρτου για το επόμενο βήμα. Στέλνεται περιοδικά, με Reporting
% Periodicity του Resource Status Request (TS 38.423 §9.1.3.18: 500, 1000,
% 2000, 5000, 10000 ms)· ανάμεσα στις αναφορές ισχύει η τελευταία. LoadAge =
% βήματα από την τελευταία αναφορά.
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

%% ------------------ Κατάσταση καναλιού για την επόμενη κλήση ------------------
% Ό,τι χρειάζεται μια continuation κλήση για τη χωρική συσχέτιση του shadow fading.
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
% Γεωμετρία όλων των υποψηφίων, μόνο για την καταγραφή των handovers.
newChannelState.SatElev_deg     = satElevMat;
newChannelState.BsDist_m        = cellRange3DMat;

end

function [chosenIdx, st, ev] = mobilityStep(candNode, candSinr, candSinrMeas, keepMask, st, P, candOffsetDb)
% Ένα βήμα κινητικότητας για έναν χρήστη σε μία πολιτική: παρακολούθηση της
% ζεύξης (RLM) του εξυπηρετούντος κόμβου και απόφαση handover κατά Event A3.
%
% (1) RLM, TS 38.331 §5.3.10.1 και §5.3.10.3, τιμές TR 36.839 Πίν. 5.2.1.3.1:
%     out-of-sync όταν SINR < Qout, in-sync όταν SINR > Qin. P.N310
%     διαδοχικά out-of-sync ξεκινούν τον T310· P.N311 διαδοχικά in-sync όσο
%     τρέχει τον σταματούν· στη λήξη του (P.NT310 βήματα) δηλώνεται RLF. Η
%     ένδειξη in-sync ελέγχεται πριν τη λήξη: η μέτρηση του βήματος της λήξης
%     ανήκει στο διάστημα που έτρεχε ο χρονομετρητής.
% (2) Event A3, TS 38.331 §5.5.4.4: είσοδος M_n + Ocn - Hys > M_p για P.NTtt
%     διαδοχικά βήματα (Ocn = μετατόπιση MLB). Συγκρίνεται η φιλτραρισμένη
%     μέτρηση (§5.5.3.2), όχι το στιγμιαίο SINR.
% (3) Αποτυχία handover, TR 36.839 Ορισμός 3:
%     - η εντολή handover φτάνει ενώ τρέχει ο T310 (κριτήριο 1 της
%       κατάστασης 2): ο χρήστης μένει στον εξυπηρετούντα με τον T310 να
%       συνεχίζει (§5.4.1)·
%     - RLF ενώ ισχύει η συνθήκη εισόδου του A3 (κριτήριο 2 της κατάστασης 2)·
%     - SINR του στόχου < Qout στο τέλος της εκτέλεσης (κατάσταση 3): ο
%       χρήστης δεν συνδέεται στον στόχο και κάνει επανεγκατάσταση.
% (4) Μετά από RLF ή αποτυχία στον στόχο ο χρήστης κάνει επανεγκατάσταση
%     (TS 38.331 §5.3.7) σε κελί που επιλέγει εκ νέου. Το TR 36.839 αφήνει
%     το μοντέλο ανάκαμψης ελεύθερο («companies should have the flexibility to
%     choose a realistic RLF recovery model», §5.4.1)· εδώ είναι ο ίδιος
%     κανόνας με την πρώτη σύνδεση: ο καλύτερος υποψήφιος στη μέτρηση με SINR
%     τουλάχιστον όσο του MCS 0, αλλιώς outage. Η διάρκεια της
%     επανεγκατάστασης δεν μοντελοποιείται (μικρότερη του βήματος).
%
% Δύο SINR, σκόπιμα: candSinr είναι το πραγματικό SINR της ζεύξης (RLM,
% έλεγχος στόχου, επιλογή κελιού), candSinrMeas η φιλτραρισμένη μέτρηση του
% L3 που βλέπει ο A3.
%
% Επιστρέφει: chosenIdx δείκτη στον κατάλογο υποψηφίων (0 = κανένας: outage
% αν st.Node == "None", αλλιώς ο εξυπηρετών χωρίς ζεύξη στο μοντέλο), τη νέα
% κατάσταση st και τα γεγονότα ev (Cause "", "HO", "RLF", "HOF"· Rlf· HoFail).
ev = struct('Cause', "", 'Rlf', false, 'HoFail', false);
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
    st.Pending = ""; st.Count = 0;
    return;
end
if mBest - P.HysDb > mServ
    if st.Pending == candNode(iBest)
        st.Count = st.Count + 1;
    else
        st.Pending = candNode(iBest);
        st.Count = 1;
    end
else
    % Η συνθήκη έπαψε να ισχύει: ο μετρητής μηδενίζεται.
    st.Pending = ""; st.Count = 0;
    return;
end
if st.Count < P.NTtt
    return;
end

% ---- (3) Εντολή handover ----
st.Pending = ""; st.Count = 0;
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
% Επιλογή κελιού χωρίς εξυπηρετούντα (πρώτη σύνδεση, μετά από outage ή μετά
% από RLF/αποτυχία): ο καλύτερος υποψήφιος στη μέτρηση, μεταξύ όσων έχουν
% SINR τουλάχιστον όσο του MCS 0. Κανένας: outage.
usable = valid & (candSinr >= P.MinUsableDb);
st.Pending = ""; st.Count = 0; st.Oos = 0; st.Is = 0; st.T310 = -1;
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

function [capacityMbpsVec, throughputMbpsVec, serviceStateVec, ...
          nodePowerWattsVec, energyPerBitUJVec, networkEnergy] = ...
          allocateAndTally(nodeVec, typeVec, sinrDbVec, simParameters, satParameters, alloc)
% Κατανομή εύρους ζώνης, χωρητικότητα, κατάσταση υπηρεσίας, ισχύς και ενέργεια
% ανά bit, για ΔΕΔΟΜΕΝΗ ανάθεση χρηστών σε κόμβους. Απομονωμένη ώστε η ίδια
% λογική να εφαρμόζεται και στις εναλλακτικές πολιτικές.
numUsers = numel(nodeVec);
numBs    = alloc.numBs;

capacityMbpsVec   = nan(numUsers,1);
throughputMbpsVec = zeros(numUsers,1);
serviceStateVec   = strings(numUsers,1);
nodePowerWattsVec = nan(numUsers,1);
energyPerBitUJVec = inf(numUsers,1);
% Πλήθος χρηστών ανά κόμβο, μία φορά για όλους.
[~, ~, nodeGrp] = unique(nodeVec);
nodeCount = accumarray(nodeGrp(:), 1);

for u = 1:numUsers
    % Outage: μηδενική χωρητικότητα/ισχύς, ενέργεια/bit = Inf. Παραλείπονται
    % πριν το usersOnThisNode ώστε να μη μετρηθούν σαν να μοιράζονται κόμβο.
    if typeVec(u) == "Outage"
        % Καμία ενεργή ζεύξη: η χωρητικότητα δεν ορίζεται (δεν υπάρχει ζεύξη
        % να τη φέρει), η παραδοθείσα ρυθμαπόδοση είναι μηδενική, και η
        % ενέργεια ανά παραδοθέν bit απροσδιόριστη. Η κατανάλωση των κόμβων
        % συνεχίζεται και προσμετράται στο ισοζύγιο δικτύου παρακάτω.
        capacityMbpsVec(u)   = NaN;
        throughputMbpsVec(u) = 0;
        nodePowerWattsVec(u) = NaN;
        energyPerBitUJVec(u) = Inf;
        serviceStateVec(u)   = "Outage";
        continue;
    end

    % Πόσοι χρήστες συνολικά εξυπηρετούνται από τον ΙΔΙΟ κόμβο
    usersOnThisNode = nodeCount(nodeGrp(u));

    % Συνολικό bandwidth και κατανάλωση ισχύος του κόμβου: μοντέλο EARTH
    % (Auer et al. 2011) για BS, γραμμικό μοντέλο ενισχυτή ισχύος για δορυφόρο.
    if typeVec(u) == "Terrestrial"
        nodeBW = alloc.BW_bs;
        nodePowerW = simParameters.Power.NumTrx * ...
            (simParameters.Power.P0 + simParameters.Power.DeltaP * alloc.pOutPerChainW);
    else
        nodeBW = satParameters.Bandwidth;
        pOutW  = 10^((satParameters.TxPower - 30)/10);
        nodePowerW = satParameters.Power.Pfix + pOutW / satParameters.Power.EtaPA;
    end

    % Κατανομή πόρων (B_user = BW_grid / N_users)
    B_user = nodeBW / usersOnThisNode;

    % Χωρητικότητα Shannon, με clamp στη μέγιστη φασματική απόδοση του NR
    % (MCS 28 / 64QAM, TS 38.214 Πίν. 5.1.3.1-1) ώστε να μην υπερεκτιμάται σε υψηλό SINR.
    % Κάτω από το SINR του MCS 0 (ή χωρίς ζεύξη στο μοντέλο, NaN) δεν
    % αποκωδικοποιείται τίποτα: ο χρήστης είναι συνδεδεμένος, όσο τρέχει ο T310,
    % αλλά χωρίς ρυθμό. Πριν το RLF αυτό δεν συνέβαινε ποτέ, γιατί κάθε
    % επιλεγμένος κόμβος είχε SINR τουλάχιστον όσο του MCS 0.
    sinr_lin = 10^(sinrDbVec(u)/10);
    if sinrDbVec(u) >= alloc.minUsableDb
        spectralEfficiency = min(log2(1 + sinr_lin), alloc.maxSe);
    else
        spectralEfficiency = 0;
    end
    capacity = B_user * spectralEfficiency;   % bits/s

    capacityMbpsVec(u)   = capacity * 1e-6;    % Mbps
    throughputMbpsVec(u) = capacity * 1e-6;    % ενεργή ζεύξη -> παραδίδεται

    % Κατάσταση υπηρεσίας: το κριτήριο ορίζεται επί του συνολικού εύρους
    % καναλιού, οπότε σφίγγει καθώς αυξάνεται ο φόρτος του κόμβου.
    if (capacity / nodeBW) >= alloc.targetNormalizedSe
        serviceStateVec(u) = "Served";
    else
        serviceStateVec(u) = "BelowTarget";
    end

    % Ενεργειακό proxy: ισομερής κατανομή ισχύος κόμβου ανά χρήστη (ίδια λογική
    % με το bandwidth split), διαιρεμένη με τον ρυθμό bit του χρήστη -> µJ/bit
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

%% ------------------ Ενεργειακή απόδοση σε επίπεδο δικτύου ------------------
% Η ενέργεια ανά bit της (eq. energy-per-bit) είναι ανεξάρτητη του πλήθους
% χρηστών: το L_n απλοποιείται. Ο δείκτης bit/J αθροίζει ισχύ ΑΝΑ ΚΟΜΒΟ (όχι
% ανά χρήστη) και συνολικό ρυθμό, οπότε αποτυπώνει τη σύνθεση των χρηστών.
% Δίνονται δύο εμβέλειες: πλήρης κόμβος και μόνο ενισχυτής (συμμετρική).
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
% Πλήθος δορυφόρων που εξυπηρετούν τουλάχιστον έναν χρήστη.
satActive = numel(unique(nodeVec(typeVec == "Satellite")));

% Πλήρης εμβέλεια: ενεργοί BS κατά EARTH, αδρανείς σε Psleep. Ο δορυφόρος δεν
% έχει αντίστοιχο μέγεθος αδράνειας, οπότε προσμετράται μόνο όταν εξυπηρετεί,
% ένας όρος ανά ενεργό δορυφόρο.
networkEnergy.TotalPower_W    = activeBs*bsFullW + idleBs*bsSleepW + satActive*satFullW;
% Συμμετρική εμβέλεια: μόνο το τμήμα που εξαρτάται από τον ενισχυτή, και στα δύο σκέλη.
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
% Στρογγυλεύει μια μετατόπιση στην πλησιέστερη τιμή του Q-OffsetRange
% (TS 38.331: "a cell, beam or measurement object specific offset to be
% applied ... when evaluating triggering conditions for measurement reporting"):
%     -24, -22, ..., -8, -6, -5, ..., 5, 6, 8, ..., 24 dB.
% Η στρογγύλευση γίνεται στην απόλυτη τιμή και επανέρχεται το πρόσημο, ώστε
% να διατηρείται η αντισυμμετρία του κανόνα (q(-x) = -q(x)) και στις ισοπαλίες,
% όπου κρατείται η μικρότερη κατά μέτρο τιμή.
% Δέχεται και διάνυσμα στήλη (στοιχείο προς στοιχείο).
levels = [0:6, 8:2:24];
[~, k] = min(abs(levels - abs(offsetDb(:))), [], 2);
q = reshape(sign(offsetDb(:)) .* levels(k)', size(offsetDb));
end

function [nodes, rho] = virtualLoad(nodeVec, sinrDbVec, targetNormalizedSe, maxSe, minUsableDb)
% Εικονικός φόρτος κάθε κόμβου κατά Awada et al. (PIMRC 2010, εξ. 1-3): το
% μερίδιο των πόρων του κόμβου που χρειάζεται ο χρήστης u για να ικανοποιηθεί,
%     κ_u = D_u / (R(SINR_u) * N_tot),     ρ_c = Σ_{u στον c} κ_u,
% που μπορεί να ξεπεράσει το 1 ("virtual" metric). Με μοντέλο πλήρους
% απασχόλησης όλοι οι πόροι είναι πάντα κατειλημμένοι, οπότε η κατάληψη PRB
% δεν ξεχωρίζει έναν γεμάτο κόμβο από έναν άδειο· ο εικονικός φόρτος ναι. Ως
% απαίτηση D_u χρησιμοποιείται το ίδιο κριτήριο εξυπηρέτησης με την κατανομή
% πόρων (targetNormalizedSe επί του εύρους του κόμβου), άρα κ_u = target/SE_u
% και ρ_c <= 1 σημαίνει ότι ο κόμβος μπορεί να εξυπηρετήσει όλους τους χρήστες του.
% Χρήστης κάτω από το MCS 0 (όσο τρέχει ο T310) δεν δρομολογείται και δεν
% καταναλώνει πόρους· με R = 0 το κ_u δεν ορίζεται.
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
% Δείκτες μιας πολιτικής για ένα βήμα, σε μορφή έτοιμη για συνάθροιση ανά
% διέλευση. Η μέση χωρητικότητα υπολογίζεται στους χρήστες που εξυπηρετούνται,
% ενώ η συνολική ρυθμαπόδοση και τα ποσοστά ορίζονται σε όλους - ώστε μια
% πολιτική να μην ωφελείται επειδή αφήνει χρήστες εκτός.
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
s.FracSatellite    = networkEnergy.SatelliteUsers / numUsers;   % μερίδιο χρηστών στο δορυφορικό σκέλος
s.TotalPower_W     = networkEnergy.TotalPower_W;
s.TotalPowerRf_W   = networkEnergy.TotalPowerRf_W;
s.BitPerJoule      = networkEnergy.BitPerJoule;
s.BitPerJouleRf    = networkEnergy.BitPerJouleRf;
% Φάσμα που έχει στη διάθεσή της η πολιτική σε αυτό το βήμα. Ορίζεται από τη
% ΔΙΑΘΕΣΙΜΟΤΗΤΑ και όχι από τις επιλογές της, ώστε να μένει συγκρίσιμο: οι
% σταθμοί βάσης εκπέμπουν πάντα (όλοι στο ίδιο φάσμα, επαναχρησιμοποίηση 1),
% ενώ μια δορυφορική ζώνη μετράει μόνο όταν κάποιος δορυφόρος της είναι πάνω
% από τη γωνία αποκοπής.
% Χρησιμεύει για να χωριστεί το πλεονέκτημα της μεικτής πολιτικής σε μέρος που
% οφείλεται στην απόφαση και σε μέρος που οφείλεται στο διαθέσιμο φάσμα.
s.Spectrum_MHz     = spectrumMHz;
% Handovers ανά χρήστη ανά λεπτό. Ο ρυθμός (και όχι το πλήθος) κάνει τα βήματα
% και τις διελεύσεις με διαφορετικό πλήθος χρηστών συγκρίσιμα, και ο μέσος
% όρος του πάνω στα βήματα της διέλευσης είναι ο ρυθμός της διέλευσης.
perUserMin = 60 / (numUsers * dtS);
s.HoRate_perUserMin         = (ho.TerrTerr + ho.TerrSat + ho.SatSat) * perUserMin;
s.HoRateTerrTerr_perUserMin = ho.TerrTerr * perUserMin;
s.HoRateTerrSat_perUserMin  = ho.TerrSat  * perUserMin;
s.HoRateSatSat_perUserMin   = ho.SatSat   * perUserMin;
s.PingPongRate_perUserMin   = ho.PingPong * perUserMin;
% RLF ανά χρήστη ανά λεπτό (TR 36.839 Ορισμός 2 επί 60 s), συνολικά και ανά
% τύπο εξυπηρετούντος κόμβου, και αποτυχίες handover (Ορισμός 3). Το ποσοστό
% αποτυχίας (Ορισμός 4) υπολογίζεται ανά διέλευση στο runSimulation.m.
s.RlfRate_perUserMin        = (ho.RlfTerr + ho.RlfSat) * perUserMin;
s.RlfRateTerr_perUserMin    = ho.RlfTerr * perUserMin;
s.RlfRateSat_perUserMin     = ho.RlfSat  * perUserMin;
s.HoFailRate_perUserMin     = ho.HoFail  * perUserMin;
end

function [isLos, latent] = spatiallyConsistentLos(pLos, moveDistance, hasPrevState, prevLatent)
% Χωρικά συνεπής κατάσταση LOS/NLOS (TR 38.901 §7.6.3.3): η κατάσταση
% προκύπτει συγκρίνοντας μια χωρικά συσχετισμένη ομοιόμορφη μεταβλητή με την
% πιθανότητα LOS. Η υποκείμενη γκαουσιανή μεταβλητή εξελίσσεται ως AR(1) με
% εκθετική συσχέτιση (§7.4.4, εξ. 7.4-5) και απόσταση συσχέτισης 50 m για την
% κατάσταση LOS/NLOS (Πίν. 7.6.3.1-2, UMa & UMi).
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
% Συντελεστής αυτοσυσχέτισης σκίασης, TR 38.901 §7.4.4 εξ. (7.4-5):
% R(Δd) = exp(-Δd/d_corr). Αποστάσεις συσχέτισης από τον Πίν. 7.5-6:
% UMa LOS=37/NLOS=50 m, UMi LOS=10/NLOS=13 m.
switch scenario
    case 'UMi'
        dCorr = 10*isLos + 13*~isLos;
    otherwise
        dCorr = 37*isLos + 50*~isLos;
end
rho = exp(-moveDistance / dCorr);
end

function fadeDb = smallScaleFadingDb(isLos, KdB, nAvg)
% Κέρδος small-scale fading σε dB, με E[|h|^2] = 1.
%   isLos=false -> Rayleigh: |h|^2 ~ Exp(1)
%   isLos=true  -> Rician με συντελεστή K (dB)· K->0 ανάγεται ομαλά σε Rayleigh
%
% nAvg: η τιμή που βγαίνει δεν είναι στιγμιαίο δείγμα αλλά μέση τιμή σε
% ΓΡΑΜΜΙΚΗ ισχύ πάνω σε nAvg ανεξάρτητες πραγματώσεις μέσα στο βήμα, επειδή
% η μέτρηση ορίζεται πάνω σε περίοδο μέτρησης και όχι σε στιγμή
% (TS 38.133 Πίν. 9.2.5.2-1). nAvg = 1 επιστρέφει το στιγμιαίο δείγμα.
if nargin < 3 || isempty(nAvg)
    nAvg = 1;
end
% Η περίπτωση nAvg = 1 κρατά κατά γράμμα την προηγούμενη έκφραση: το
% 10*log10(mean(x)) και το 20*log10(abs(h)) είναι μαθηματικά ταυτόσημα με
% τα αντίστοιχα, αλλά διαφέρουν κατά ένα ulp. Έτσι ο έλεγχος ουδετερότητας
% δίνει ταυτόσημο αρχείο και όχι «σχεδόν ταυτόσημο».
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
% Πεδίο ισχύος των μοντέλων UMa/UMi (TR 38.901 Πίν. 7.4.1-1, στήλη
% "Applicability range"): 10m <= d2D <= 5km, 1.5m <= hUT <= 22.5m.
% Η συχνότητα ελέγχεται ως προς το εύρος του ίδιου του μοντέλου (§7.4.1).
% Δέχεται και πίνακες (στοιχείο προς στοιχείο).
isValid = d2D  >= 10   & d2D  <= 5000 & ...
          hUT  >= 1.5  & hUT  <= 22.5 & ...
          fcHz >= 0.5e9 & fcHz <= 100e9;
end

function [b0, m, omega] = shadowedRicianElevParams(elevDeg)
% Παράμετροι Shadowed Rician από τη γωνία ανύψωσης (Abdi et al. 2003, εξ. 19).
% Προσαρμογή πολυωνύμων σε πειραματικά δεδομένα, ισχύει για 20° < θ < 80°.
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
% Κέρδος Shadowed Rician σε dB (Abdi et al. 2003, εξ. 1): σκεδαζόμενη
% συνιστώσα Rayleigh μέσης ισχύος 2*b0 συν συνιστώσα LOS με πλάτος
% κατά Nakagami-m μέσης ισχύος omega. E[|h|^2] = omega + 2*b0 < 1, δηλαδή
% η μέση εξασθένηση λόγω σκίασης περιέχεται στο ίδιο το μοντέλο.
%
% nAvg: μέση τιμή της μέτρησης μέσα στο βήμα, όπως και στο επίγειο σκέλος.
% Η συνιστώσα σκίασης κληρώνεται ΜΙΑ φορά και μένει σταθερή μέσα στο βήμα:
% η σκίαση δεν είναι γρήγορο φαινόμενο και δεν εξομαλύνεται από την περίοδο
% μέτρησης. Μέσος όρος παίρνεται μόνο στη σκεδαζόμενη συνιστώσα.
if nargin < 4 || isempty(nAvg)
    nAvg = 1;
end
% Z^2 ~ Gamma(m, omega/m). Με uShadow (εκατοστημόριο της σκίασης από τη
% συσχετισμένη μεταβλητή) η τιμή παίρνεται από την αντίστροφη αθροιστική· η
% gammaincinv ανήκει στη βασική MATLAB. Χωρίς uShadow, ανεξάρτητο δείγμα.
if nargin >= 5 && ~isempty(uShadow)
    uShadow = min(max(uShadow, eps), 1 - eps);
    losAmp  = sqrt((omega/m) * gammaincinv(uShadow, m));   % |Z|, E[Z^2] = omega
else
    losAmp  = sqrt(gammaRand(m, omega/m));
end
scatter = sqrt(b0)*(randn(1,nAvg) + 1i*randn(1,nAvg));   % E[|A|^2] = 2*b0
if nAvg == 1
    fadeDb = 20*log10(abs(losAmp + scatter));   % ίδια έκφραση με πριν (βλ. smallScaleFadingDb)
else
    fadeDb = 10*log10(mean(abs(losAmp + scatter).^2));
end
end

function f = l3Filter(mNow, fPrev, alpha)
% Φίλτρο Layer 3, TS 38.331 §5.5.3.2:  F_n = (1 - a)*F_{n-1} + a*M_n
% Εφαρμόζεται σε dB, όπως απαιτεί η NOTE 2 ("logarithmic filtering for
% logarithmic measurements"). Χωρίς προηγούμενη μέτρηση ισχύει F_0 = M_1,
% ρητά στο ίδιο εδάφιο. Το alpha έρχεται ήδη προσαρμοσμένο στο βήμα της
% προσομοίωσης· alpha = 1 σημαίνει καθόλου φίλτρο (NOTE 1, k = 0).
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
% Δείγμα από κατανομή Gamma (Marsaglia & Tsang 2000). Υλοποιείται τοπικά
% ώστε να μη χρειάζεται το Statistics Toolbox· δέχεται και shape < 1.
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
% Πιθανότητα LOS ζεύξης BS-χρήστη (TR 38.901 Πίν. 7.4.2-1). d2D, hUT σε μέτρα,
% βαθμωτά ή πίνακες ίδιου μεγέθους (στοιχείο προς στοιχείο).
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
% Απόσβεση ατμοσφαιρικών αερίων (O2 + υδρατμοί) σε ζεύξη δορυφόρου-χρήστη.
% TR 38.811 §6.6.4 εξ. (6.6-8): PLA = A_zenith / sin(ε).
% A_zenith = γ_o·h_o + γ_w·h_w· ισοδύναμα ύψη κατά ITU-R P.676-12 Annex 2·
% ειδικές αποσβέσεις γ_o/γ_w από gaspl (ITU-R P.676-13 Annex 1)·
% reference atmosphere κατά ITU-R P.835.
% Το A_zenith δεν εξαρτάται από την ανύψωση: υπολογίζεται μία φορά ανά
% συχνότητα και κρατείται, και ανά ζεύξη μένει μόνο η διαίρεση με sin(ε).
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
