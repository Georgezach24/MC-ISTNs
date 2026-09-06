function [bestNodeVec, bestNodeTypeVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, bestElevationDegVec, ...
    nodePowerWattsVec, energyPerBitUJVec, ...
    bestBsSnrDbVec, bestBsDistanceVec, bestBsPathLossVec, ...
    satSlantRangeVec, satElevationVec, satPathLossVec, satSnrDbVec, ...
    newChannelState] = ...
    simulateScenario(bs_geo, user_geo, sat_geo, wgs84, simParameters, satParameters, prevChannelState)
% Για κάθε χρήστη: επιλέγει τον καλύτερο κόμβο (BS ή δορυφόρο) βάσει SNR και
% υπολογίζει χωρητικότητα/ενέργεια μετά την κατανομή εύρους ζώνης.
% Επιστρέφει και per-candidate διαγνωστικά (καλύτερο BS + δορυφόρος, ανεξάρτητα
% από την τελική επιλογή) για την εκπαίδευση του ML μοντέλου του Part 2.
%
% prevChannelState (προαιρετικό, 7ο όρισμα): αν δοθεί, LOS/NLOS και shadow
% fading κάθε ζεύξης συσχετίζονται χωρικά με την προηγούμενη κλήση αντί για
% i.i.d. δειγματοληψία. Το περνάνε μόνο callers που προσομοιώνουν διαδοχικές
% μεταδόσεις της ίδιας τοπολογίας (temporalPassSimulation.m)· οι υπόλοιποι
% (test_simulation.m, monteCarloDriver.m, kpiRepeatedRuns.m) το παραλείπουν.
if nargin < 7
    prevChannelState = [];
end

numUsers = size(user_geo,1);
numBs    = size(bs_geo,1);

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

%% ------------------ Ελάχιστο χρησιμοποιήσιμο SNR (κατάσταση outage) ------------------
% Κάτω από αυτό -> ο χρήστης θεωρείται outage αντί να ανατεθεί στον
% "λιγότερο κακό" κόμβο. = Shannon-ισοδύναμο SNR του MCS 0 (TS 38.214 Πίν. 5.1.3.1-2).
minSpectralEfficiency = 0.2344;                        % bits/s/Hz (MCS 0)
minUsableSnrDb = 10*log10(2^minSpectralEfficiency - 1); % ≈ -7.53 dB

%% ------------------ Διαλείψεις δορυφορικής ζεύξης ------------------
% Shadowed Rician (Abdi et al. 2003). Η σκίαση περιέχεται ήδη στο μοντέλο
% (τυχαίο πλάτος LOS κατά Nakagami-m), οπότε δεν προστίθεται χωριστός
% λογαριθμοκανονικός όρος. Παράμετροι (b0,m,Ω) από την ανύψωση, εξ. (19).
satFadeElevRangeDeg = [20 80];   % πεδίο ισχύος της προσαρμογής της εξ. (19)
terrKdBLos          = 9;         % Rician K επίγειο LOS, TR 38.901 Πίν. 7.5-6

%% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestNodeVec         = strings(numUsers,1);
bestNodeTypeVec     = strings(numUsers,1);
bestDistanceVec     = nan(numUsers,1);
bestPathLossVec     = nan(numUsers,1);
bestSnrDbVec        = nan(numUsers,1);
capacityMbpsVec     = nan(numUsers,1);
nodePowerWattsVec   = nan(numUsers,1);
energyPerBitUJVec   = nan(numUsers,1);
bestElevationDegVec = nan(numUsers,1);
bestBsSnrDbVec      = nan(numUsers,1);
bestBsDistanceVec   = nan(numUsers,1);
bestBsPathLossVec   = nan(numUsers,1);

% Διαγνωστικοί πίνακες
groundDistanceMat = nan(numUsers,numBs);
range3DMat        = nan(numUsers,numBs);
pathLossMat       = nan(numUsers,numBs);
snrDbMat          = nan(numUsers,numBs);
pLosMat           = nan(numUsers,numBs);
losMat            = false(numUsers,numBs);
sfMat             = nan(numUsers,numBs);
satSlantRangeVec  = nan(numUsers,1);
satElevationVec   = nan(numUsers,1);
satPathLossVec    = nan(numUsers,1);
satSnrDbVec       = nan(numUsers,1);

hasPrevState = ~isempty(prevChannelState) && ...
    isequal(size(prevChannelState.IsLOS), [numUsers, numBs]);

%% ------------------ Επιλογή Καλύτερου Κόμβου (βάσει SNR) ------------------
for u = 1:numUsers
    % Αρχικοποιήσεις
    userBestSNR       = -Inf;
    userBestNode      = "";
    userBestType      = "";
    userBestDistance  = NaN;
    userBestPathLoss  = NaN;
    userBestElevation = NaN;

    % Απόσταση μετακίνησης του χρήστη από την προηγούμενη κλήση· καθορίζει
    % πόσο συσχετίζεται το shadow fading με την προηγούμενη τιμή (βλ. correlatedLosState).
    if hasPrevState
        userMoveDistance = distance(prevChannelState.UserGeo(u,1), prevChannelState.UserGeo(u,2), ...
                                     user_geo(u,1), user_geo(u,2), wgs84);
    else
        userMoveDistance = Inf; % καμία προηγούμενη κατάσταση -> ανεξάρτητο δείγμα, όπως πριν
    end

    %% ===== Terrestrial BS candidates =====
    for b = 1:numBs
        lat0 = bs_geo(b,1);
        lon0 = bs_geo(b,2);
        h0   = 0;

        [xBS, yBS, zBS] = geodetic2enu(bs_geo(b,1), bs_geo(b,2), bs_geo(b,3), ...
                                       lat0, lon0, h0, wgs84);
        [xUE, yUE, zUE] = geodetic2enu(user_geo(u,1), user_geo(u,2), user_geo(u,3), ...
                                       lat0, lon0, h0, wgs84);

        txPosition = [xBS; yBS; zBS];
        rxPosition = [xUE; yUE; zUE];

        groundDistance = distance(bs_geo(b,1), bs_geo(b,2), ...
                                  user_geo(u,1), user_geo(u,2), wgs84);
        d3d = norm(rxPosition - txPosition);
        groundDistanceMat(u,b) = groundDistance;
        range3DMat(u,b)        = d3d;

        % LOS ανά ζεύξη βάσει πιθανότητας απόστασης (TR 38.901 §7.4.2).
        pLos = losProbability38901(groundDistance, user_geo(u,3), simParameters.PathLoss.Scenario);
        if hasPrevState
            prevIsLos = prevChannelState.IsLOS(u,b);
            prevSF    = prevChannelState.ShadowFading_dB(u,b);
        else
            prevIsLos = false;
            prevSF    = 0;
        end
        [isLos, rho, useCorrelatedSF] = correlatedLosState(pLos, userMoveDistance, ...
            simParameters.PathLoss.Scenario, hasPrevState, prevIsLos);
        pLosMat(u,b) = pLos;
        losMat(u,b)  = isLos;

        [pathLoss, sigmaSF] = nrPathLoss(simParameters.PathLoss, ...
                              simParameters.CarrierFrequency, ...
                              isLos, ...
                              txPosition, rxPosition);

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
        % Realization i.i.d. ανά κλήση (coherence time ~ ms << βήμα).
        pathLoss = pathLoss - smallScaleFadingDb(isLos, terrKdBLos);
        pathLossMat(u,b) = pathLoss;

        snr_db = (simParameters.EIRP - 30) - pathLoss - noisePowerBS_dBW;
        snrDbMat(u,b) = snr_db;

        if snr_db > userBestSNR
            userBestSNR       = snr_db;
            userBestNode      = "BS" + string(b);
            userBestType      = "Terrestrial";
            userBestDistance  = d3d;
            userBestPathLoss  = pathLoss;
            userBestElevation = NaN;
        end
    end

    % Στιγμιότυπο του καλύτερου υποψήφιου BS πριν τη σύγκριση με τον δορυφόρο (per-candidate διαγνωστικό).
    bestBsSnrDbVec(u)   = userBestSNR;
    bestBsDistanceVec(u) = userBestDistance;
    bestBsPathLossVec(u) = userBestPathLoss;

    %% ===== Satellite candidate =====
    [~, elevSat, slantRangeSat] = geodetic2aer( ...
        sat_geo(1), sat_geo(2), sat_geo(3), ...
        user_geo(u,1), user_geo(u,2), user_geo(u,3), wgs84);

    satSlantRangeVec(u) = slantRangeSat;
    satElevationVec(u)  = elevSat;

    if elevSat >= satParameters.MinElevationDeg
        lambdaSat = physconst('LightSpeed') / satParameters.CarrierFrequency;
        satPathLoss = fspl(slantRangeSat, lambdaSat);

        % Ατμοσφαιρική απόσβεση αερίων (οξυγόνο + υδρατμοί), TR 38.811 §6.6.4.
        % Βροχή/νέφωση & ιονοσφαιρική σπινθηρίδα δεν μοντελοποιούνται.
        gasAttenuationDb = gasAttenuationSlantP676(satParameters.CarrierFrequency, elevSat);
        satPathLoss = satPathLoss + gasAttenuationDb;

        % Shadowed Rician· i.i.d. ανά κλήση (ο δορυφόρος κινείται -> η γεωμετρία
        % σκίασης αποσυσχετίζεται γρήγορα). Η ανύψωση περιορίζεται στο πεδίο
        % ισχύος της προσαρμογής· εναλλακτικά σταθερή κατάσταση σκίασης μέσω
        % satParameters.ShadowingState (για ανάλυση ευαισθησίας).
        if isfield(satParameters, 'ShadowingState') && ~isempty(satParameters.ShadowingState)
            [b0, mNak, omega] = shadowedRicianStateParams(satParameters.ShadowingState);
        else
            elevClamped = min(max(elevSat, satFadeElevRangeDeg(1)), satFadeElevRangeDeg(2));
            [b0, mNak, omega] = shadowedRicianElevParams(elevClamped);
        end
        satPathLoss = satPathLoss - shadowedRicianFadingDb(b0, mNak, omega);

        satSnrDb = (satParameters.EIRP - 30) - satPathLoss - noisePowerSAT_dBW;
    else
        satPathLoss = inf;
        satSnrDb = -Inf;
    end

    satPathLossVec(u) = satPathLoss;
    satSnrDbVec(u)    = satSnrDb;

    if satSnrDb > userBestSNR
        userBestSNR       = satSnrDb;
        userBestNode      = "SAT-1";
        userBestType      = "Satellite";
        userBestDistance  = slantRangeSat;
        userBestPathLoss  = satPathLoss;
        userBestElevation = elevSat;
    end

    % Σε outage κρατάμε τα διαγνωστικά του καλύτερου υποψηφίου, αλλά ο χρήστης
    % δεν ανατίθεται σε κόμβο.
    if userBestSNR < minUsableSnrDb
        bestNodeVec(u)     = "None";
        bestNodeTypeVec(u) = "Outage";
    else
        bestNodeVec(u)     = userBestNode;
        bestNodeTypeVec(u) = userBestType;
    end
    bestDistanceVec(u)     = userBestDistance;
    bestPathLossVec(u)     = userBestPathLoss;
    bestSnrDbVec(u)        = userBestSNR;
    bestElevationDegVec(u) = userBestElevation;
end

%% ------------------ Υπολογισμός Χωρητικότητας & Ενέργειας (Κατανομή Πόρων) ------------------
for u = 1:numUsers
    % Outage: μηδενική χωρητικότητα/ισχύς, ενέργεια/bit = Inf. Παραλείπονται
    % πριν το usersOnThisNode ώστε να μη μετρηθούν σαν να μοιράζονται κόμβο.
    if bestNodeTypeVec(u) == "Outage"
        capacityMbpsVec(u)   = 0;
        nodePowerWattsVec(u) = 0;
        energyPerBitUJVec(u) = Inf;
        continue;
    end

    servingNode = bestNodeVec(u);

    % Πόσοι χρήστες συνολικά εξυπηρετούνται από τον ΙΔΙΟ κόμβο
    usersOnThisNode = sum(bestNodeVec == servingNode);

    % Συνολικό bandwidth και κατανάλωση ισχύος του κόμβου: μοντέλο EARTH
    % (Auer et al. 2011) για BS, γραμμικό μοντέλο ενισχυτή ισχύος για δορυφόρο.
    if bestNodeTypeVec(u) == "Terrestrial"
        nodeBW = BW_bs;
        pOutW  = 10^((simParameters.TxPower - 30)/10);
        nodePowerW = simParameters.Power.NumTrx * ...
            (simParameters.Power.P0 + simParameters.Power.DeltaP * pOutW);
    else
        nodeBW = satParameters.Bandwidth;
        pOutW  = 10^((satParameters.TxPower - 30)/10);
        nodePowerW = satParameters.Power.Pfix + pOutW / satParameters.Power.EtaPA;
    end

    % Κατανομή πόρων (B_user = BW_grid / N_users)
    B_user = nodeBW / usersOnThisNode;

    % Χωρητικότητα Shannon, με clamp στη μέγιστη φασματική απόδοση του NR
    % (MCS 27 / 256QAM, TS 38.214 Πίν. 5.1.3.1-2) ώστε να μην υπερεκτιμάται σε υψηλό SNR.
    maxSpectralEfficiency = 5.5547;   % bits/s/Hz (MCS 27)
    snr_lin = 10^(bestSnrDbVec(u)/10);
    spectralEfficiency = min(log2(1 + snr_lin), maxSpectralEfficiency);
    capacity = B_user * spectralEfficiency;   % bits/s

    capacityMbpsVec(u) = capacity * 1e-6;    % Mbps

    % Ενεργειακό proxy: ισομερής κατανομή ισχύος κόμβου ανά χρήστη (ίδια λογική
    % με το bandwidth split), διαιρεμένη με τον ρυθμό bit του χρήστη -> µJ/bit
    nodePowerWattsVec(u) = nodePowerW;
    energyPerBitUJVec(u) = (nodePowerW / usersOnThisNode) / capacity * 1e6;
end

%% ------------------ Κατάσταση καναλιού για την επόμενη κλήση ------------------
% Ό,τι χρειάζεται μια continuation κλήση για τη χωρική συσχέτιση του shadow fading.
newChannelState.UserGeo         = user_geo;
newChannelState.IsLOS           = losMat;
newChannelState.ShadowFading_dB = sfMat;

end

function [isLos, rho, useCorrelatedSF] = correlatedLosState(pLos, moveDistance, scenario, hasPrevState, prevIsLos)
% Συσχετισμένη κατάσταση LOS/NLOS ζεύξης BS-χρήστη μεταξύ διαδοχικών κλήσεων.
% ρ(Δd) = exp(-Δd/d_corr), Gudmundson (1991)· d_corr από TR 38.901 v16.1.0
% Πίν. 7.5-6 (SF correlation distance): UMa LOS=37/NLOS=50, UMi LOS=10/NLOS=13 m.
% Το ίδιο ρ χρησιμοποιείται και ως πιθανότητα διατήρησης της προηγούμενης
% κατάστασης LOS/NLOS (το TR 38.901 δεν ορίζει ξεχωριστό d_corr γι' αυτήν).
if ~hasPrevState
    isLos = rand() < pLos;
    rho = 0;
    useCorrelatedSF = false;
    return;
end

if prevIsLos
    switch scenario
        case 'UMa'
            dCorr = 37;
        case 'UMi'
            dCorr = 10;
        otherwise
            dCorr = 37;
    end
else
    switch scenario
        case 'UMa'
            dCorr = 50;
        case 'UMi'
            dCorr = 13;
        otherwise
            dCorr = 50;
    end
end
rho = exp(-moveDistance / dCorr);

if rand() < rho
    isLos = prevIsLos;
else
    isLos = rand() < pLos;
end

% Το AR(1) δείγμα SF ισχύει μόνο αν δεν άλλαξε η κατάσταση LOS/NLOS (αλλάζει το σ_SF).
useCorrelatedSF = (isLos == prevIsLos);
end

function fadeDb = smallScaleFadingDb(isLos, KdB)
% Κέρδος small-scale fading σε dB, με E[|h|^2] = 1.
%   isLos=false -> Rayleigh: |h|^2 ~ Exp(1)
%   isLos=true  -> Rician με συντελεστή K (dB)· K->0 ανάγεται ομαλά σε Rayleigh
if ~isLos
    fadeDb = 10*log10(-log(rand()));
else
    Klin  = 10^(KdB/10);
    s     = sqrt(Klin/(Klin+1));       % πλάτος LOS συνιστώσας
    sigma = sqrt(1/(2*(Klin+1)));      % τυπ. απόκλιση ανά διάσταση scatter
    h     = (s + sigma*randn()) + 1i*(sigma*randn());
    fadeDb = 20*log10(abs(h));         % s^2 + 2*sigma^2 = 1
end
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

function fadeDb = shadowedRicianFadingDb(b0, m, omega)
% Κέρδος Shadowed Rician σε dB (Abdi et al. 2003, εξ. 1): σκεδαζόμενη
% συνιστώσα Rayleigh μέσης ισχύος 2*b0 συν συνιστώσα LOS με πλάτος
% κατά Nakagami-m μέσης ισχύος omega. E[|h|^2] = omega + 2*b0 < 1, δηλαδή
% η μέση εξασθένηση λόγω σκίασης περιέχεται στο ίδιο το μοντέλο.
losAmp  = sqrt(gammaRand(m, omega/m));            % |Z|, E[Z^2] = omega
scatter = sqrt(b0)*(randn() + 1i*randn());        % E[|A|^2] = 2*b0
fadeDb  = 20*log10(abs(losAmp + scatter));
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
% Πιθανότητα LOS ζεύξης BS-χρήστη (TR 38.901 Πίν. 7.4.2-1). d2D, hUT σε μέτρα.
switch scenario
    case 'UMi'
        if d2D <= 18
            pLos = 1;
        else
            pLos = 18/d2D + exp(-d2D/36) * (1 - 18/d2D);
        end
    case 'UMa'
        if d2D <= 18
            pLos = 1;
        else
            if hUT <= 13
                Cprime = 0;
            else
                g = 1.25e-6 * d2D^3 * exp(-d2D/150);
                Cprime = ((hUT - 13)/10)^1.5 * g;
            end
            pLos = (18/d2D + exp(-d2D/63) * (1 - 18/d2D)) * (1 + Cprime);
        end
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
