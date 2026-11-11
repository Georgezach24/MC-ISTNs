function T = linkBudgetValidation(outputDir, label)
%LINKBUDGETVALIDATION Επαλήθευση ισοζυγίου ζεύξης και κατανομής πόρων.
%   Α. Μονάδες ισχύος και EIRP          ΣΤ. Διατήρηση εύρους ζώνης
%   Β. Απώλειες ελεύθερου χώρου         Ζ.  Ισοζύγιο ισχύος δικτύου
%   Γ. Θόρυβος kTB                      Η.  Οριακές περιπτώσεις
%   Δ. SINR = EIRP - PL - N, παρεμβολή  Θ.  Κεραία τομέα (TR 38.901 Πίν. 7.3-1, 7.8-1)
%   Ε. Μεταβολή ισχύος εκπομπής         Ι.  Soft LOS και RLF
% Κάθε έλεγχος συγκρίνει τιμή κλειστού τύπου με την έξοδο της simulateScenario·
% σφάλμα αν αποτύχει έστω ένας.

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

C = addCheck(C, 'A1 gNB transmit power', 'W', 79.432823472428, pOutTotalW, 1e-9);
C = addCheck(C, 'A2 power per transceiver chain', 'W', 19.858205868107, pOutPerChainW, 1e-9);
C = addCheck(C, 'A3 EARTH limit per chain (<= 20 W)', 'logical', 1, double(pOutPerChainW <= 20), 0);
C = addCheck(C, 'A4 composite gNB EIRP', 'dBm', 67, simParameters.EIRP, 1e-12);
C = addCheck(C, 'A5 same EIRP in linear power', 'W', 5011.872336272722, ...
             10^((simParameters.EIRP - 30)/10), 1e-6);
densityBack = satParameters.EIRP - 30 - 10*log10(satParameters.Bandwidth/1e6);
C = addCheck(C, 'A6 satellite EIRP density (round trip)', 'dBW/MHz', 34, densityBack, 1e-12);
C = addCheck(C, 'A7 satellite EIRP', 'dBm', 78.771212547197, satParameters.EIRP, 1e-9);
C = addCheck(C, 'A8 power at the satellite antenna input', 'dBm', 48.771212547197, ...
             satParameters.TxPower, 1e-9);

%% ================== Β. Απώλειες ελεύθερου χώρου ==================
lambdaSat  = cLight / satParameters.CarrierFrequency;
dFspl      = [600e3 1200e3 1754e3];
fsplClosed = 20*log10(4*pi*dFspl/lambdaSat);
fsplCode   = arrayfun(@(d) fspl(d, lambdaSat), dFspl);
C = addCheck(C, 'B1 fspl == 20log10(4*pi*d/lambda)', 'dB', 0, ...
             max(abs(fsplCode - fsplClosed)), 1e-9);
C = addCheck(C, 'B2 FSPL at 600 km, 2.0 GHz', 'dB', 154.031408142836, fsplCode(1), 1e-9);
octaveD    = 10.^(2:0.5:6);
fsplOctave = arrayfun(@(d) fspl(2*d, lambdaSat) - fspl(d, lambdaSat), octaveD);
C = addCheck(C, 'B3 doubled distance -> +6.02 dB', 'dB', 20*log10(2), ...
             max(fsplOctave), 1e-9);

% UMa: εκθέτης απωλειών > 2, άρα οκτάβα > 6.02 dB (nrPathLoss χωρίς σκίαση).
plCfgU = simParameters.PathLoss;
plU1 = nrPathLoss(plCfgU, simParameters.CarrierFrequency, false, [0;0;25], [500;0;1.5]);
plU2 = nrPathLoss(plCfgU, simParameters.CarrierFrequency, false, [0;0;25], [1000;0;1.5]);
umaOctaveDb = plU2 - plU1;
C = addCheck(C, 'B4 terrestrial octave > 6.02 dB (exponent > 2)', 'logical', 1, ...
             double(umaOctaveDb > 20*log10(2)), 0);

%% ================== Γ. Ισχύς θορύβου ==================
C = addCheck(C, 'G1 gNB bandwidth', 'Hz', 18.36e6, bwBsHz, 1e-6);
C = addCheck(C, 'G2 equivalent noise temperature', 'K', ...
             290 + 290*(10^0.9 - 1), teqK, 1e-9);
noiseBsSplit = 10*log10(kBoltz) + 10*log10(teqK) + 10*log10(bwBsHz);
C = addCheck(C, 'G3 kTB: product == sum in dB', 'dB', 0, ...
             abs(noiseBsDbw - noiseBsSplit), 1e-9);
C = addCheck(C, 'G4 noise power (gNB link)', 'dBW', -122.336460425576, noiseBsDbw, 1e-9);
C = addCheck(C, 'G5 noise power (satellite link)', 'dBW', -120.203974647031, noiseSatDbw, 1e-9);
noiseDouble = 10*log10(kBoltz * teqK * 2*bwBsHz);
C = addCheck(C, 'G6 doubled bandwidth -> +3.01 dB noise', 'dB', 10*log10(2), ...
             noiseDouble - noiseBsDbw, 1e-9);

%% ================== Δ. Ταυτότητα SINR ==================
% Δύο σταθμοί, πέντε χρήστες, δορυφόρος στο ζενίθ. Με έναν σταθμό SINR = SNR
% ακριβώς· με δύο SINR < SNR (παρεμβολή).
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

rng(21);
[~,~,~,~,~,~,~,~,~, bsSinr1, ~, bsPl1] = ...
    simulateScenario(bs1, uGeo, satOverhead, wgs84, simParameters, satParameters);
sinrFromPl1 = (simParameters.EIRP - 30) - bsPl1 - noiseBsDbw;
C = addCheck(C, 'D1 single gNB: SINR == EIRP - PL - N', 'dB', 0, ...
             max(abs(bsSinr1 - sinrFromPl1)), 1e-9);

% Ένας δορυφόρος, άλλη ζώνη από τους σταθμούς: SINR = SNR.
satSinrFromPl = (satParameters.EIRP - 30) - satPl - noiseSatDbw;
C = addCheck(C, 'D2 satellite SINR == EIRP - PL - N', 'dB', 0, ...
             max(abs(satSinr - satSinrFromPl)), 1e-9);
C = addCheck(C, 'D3 selected SINR == max(TN, NTN)', 'dB', 0, ...
             max(abs(sinrBest - max(bsSinr, satSinr))), 1e-12);

% Δύο σταθμοί: SINR < SNR όπου ισχύουν και οι δύο ζεύξεις, αλλιώς SINR = SNR.
sinrNoInterf = (simParameters.EIRP - 30) - bsPl - noiseBsDbw;
interfDropDb = sinrNoInterf - bsSinr;
d2dBoth = [distance(bs2(1,1), bs2(1,2), uGeo(:,1), uGeo(:,2), wgs84), ...
           distance(bs2(2,1), bs2(2,2), uGeo(:,1), uGeo(:,2), wgs84)];
twoValid = all(d2dBoth >= 10 & d2dBoth <= 5000, 2);
C = addCheck(C, 'D4 two valid gNBs: SINR < SNR', 'logical', 1, ...
             double(any(twoValid) && all(interfDropDb(twoValid) > 0)), 0);
C = addCheck(C, 'D4 one valid gNB: SINR == SNR', 'dB', 0, ...
             max([0; abs(interfDropDb(~twoValid))]), 1e-9);
C = addCheck(C, 'D5 the drop equals a positive interference power', 'logical', 1, ...
             double(all(isfinite(interfDropDb)) && max(interfDropDb) < 60), 0);

%% ================== Ε. Μεταβολή ισχύος εκπομπής ==================
% Ένας σταθμός (θόρυβος): ±3 dB ισχύος -> ±3 dB SINR. Δύο σταθμοί
% (παρεμβολή): μικρότερη μεταβολή.
offsets    = [-3 0 3];
sinrShift1 = nan(numel(offsets),1);   % ένας σταθμός
sinrShift2 = nan(numel(offsets),1);   % δύο σταθμοί
for k = 1:numel(offsets)
    sp = simParameters;
    sp.TxPower = simParameters.TxPower + offsets(k);
    sp.EIRP = sp.TxPower + sp.AntennaGain + 10*log10(sp.NumAntennaElements);
    % Στο +3 dB διπλές αλυσίδες, για το όριο 20 W/αλυσίδα (δεν επηρεάζει το SINR).
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
C = addCheck(C, 'E1 single gNB: -3 dB power -> -3 dB SINR', 'dB', -3, sinrShift1(1), 1e-9);
C = addCheck(C, 'E2 same power -> zero shift', 'dB', 0, sinrShift1(2), 1e-12);
C = addCheck(C, 'E3 single gNB: +3 dB power -> +3 dB SINR', 'dB', 3, sinrShift1(3), 1e-9);
C = addCheck(C, 'E4 two gNBs: +3 dB power gives LESS than +3 dB', 'logical', 1, ...
             double(sinrShift2(3) > 0 && sinrShift2(3) < 3 - 1e-6), 0);
C = addCheck(C, 'E5 two gNBs: -3 dB power costs LESS than 3 dB', 'logical', 1, ...
             double(sinrShift2(1) < 0 && sinrShift2(1) > -3 + 1e-6), 0);

%% ================== ΣΤ. Διατήρηση εύρους ζώνης ==================
% B_u = C_u / SE_u· το άθροισμα ανά κόμβο = εύρος του κόμβου.
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
C = addCheck(C, 'ST1 sum of shares == node bandwidth', 'Hz', 0, ...
             max(abs(bwSumPerNode - nodeBwPerNode)), 1e-3);
C = addCheck(C, 'ST2 no node over-allocates', 'logical', 1, ...
             double(all(bwSumPerNode <= nodeBwPerNode + 1e-3)), 0);

%% ================== Ζ. Ισοζύγιο ισχύος δικτύου ==================
bsFullW  = simParameters.Power.NumTrx * ...
           (simParameters.Power.P0 + simParameters.Power.DeltaP*pOutPerChainW);
bsRfW    = simParameters.Power.NumTrx * simParameters.Power.DeltaP * pOutPerChainW;
bsSleepW = simParameters.Power.NumTrx * simParameters.Power.Psleep;
pOutSatW = 10^((satParameters.TxPower - 30)/10);
satFullW = satParameters.Power.Pfix + pOutSatW/satParameters.Power.EtaPA;

C = addCheck(C, 'Z1 active gNB power (EARTH)', 'W', 893.334270320412, bsFullW, 1e-9);
C = addCheck(C, 'Z2 PA-dependent part of the gNB power', 'W', 373.334270320412, bsRfW, 1e-9);
C = addCheck(C, 'Z3 satellite power', 'W', 188.391482363219, satFullW, 1e-9);
expTotal = netE.ActiveBs*bsFullW + netE.IdleBs*bsSleepW + netE.SatActive*satFullW;
C = addCheck(C, 'Z4 total power == active + idle + satellite', 'W', 0, ...
             abs(netE.TotalPower_W - expTotal), 1e-9);
satRfW = pOutSatW / satParameters.Power.EtaPA;
expRf = netE.ActiveBs*bsRfW + netE.SatActive*satRfW;
C = addCheck(C, 'Z5 symmetric (RF) power == PAs only', 'W', 0, ...
             abs(netE.TotalPowerRf_W - expRf), 1e-9);
C = addCheck(C, 'Z6 RF scope <= full scope', 'logical', 1, ...
             double(netE.TotalPowerRf_W <= netE.TotalPower_W), 0);
C = addCheck(C, 'Z7 every gNB counted once', 'count', size(bs2,1), ...
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

C = addCheck(C, 'H1 no node -> ServingNode "None"', 'logical', 1, ...
             double(nodeN == "None"), 0);
C = addCheck(C, 'H1 no node -> ServingType "Outage"', 'logical', 1, ...
             double(typeN == "Outage"), 0);
C = addCheck(C, 'H1 capacity undefined', 'logical', 1, double(isnan(capN)), 0);
C = addCheck(C, 'H1 throughput zero', 'Mbps', 0, thrN, 0);
C = addCheck(C, 'H1 node power undefined', 'logical', 1, double(isnan(powN)), 0);
C = addCheck(C, 'H1 energy per bit infinite', 'logical', 1, double(isinf(eN) && eN > 0), 0);
C = addCheck(C, 'H1 reason for TN unavailability', 'logical', 1, ...
             double(bsRN == "OutOfModelRange"), 0);
C = addCheck(C, 'H1 reason for NTN unavailability', 'logical', 1, ...
             double(satRN == "NotVisible"), 0);
C = addCheck(C, 'H1 network consumption does not vanish', 'W', bsSleepW, ...
             netN.TotalPower_W, 1e-9);
C = addCheck(C, 'H2 -Inf tie: no candidate wins', 'logical', 1, ...
             double(isnan(bsSinrN) && satSinrN == -Inf && sinrN == -Inf && nodeN == "None"), 0);

% --- Η3: ένας χρήστης παίρνει ολόκληρο το εύρος ζώνης ---
oneUser = [baseLat + 400/111320, baseLon, 1.5];
satNo   = satParameters; satNo.MinElevationDeg = 95;   % δορυφόρος απρόσιτος
rng(13);
[~, type1, ~, ~, sinr1, cap1] = simulateScenario([baseLat baseLon 25], oneUser, ...
    satOverhead, wgs84, simParameters, satNo);
se1 = min(log2(1 + 10^(sinr1/10)), maxSe);
C = addCheck(C, 'H3 one UE -> served by TN', 'logical', 1, ...
             double(type1 == "Terrestrial"), 0);
C = addCheck(C, 'H3 one UE -> full bandwidth', 'Hz', bwBsHz, ...
             cap1*1e6/se1, 1e-3);

% --- Η4: ανύψωση ακριβώς στη μάσκα (κλειστό όριο) ---
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
C = addCheck(C, 'H4 elevation exactly at the mask -> visible', 'logical', 1, ...
             double(isfinite(sinrAt) && reasonAt ~= "NotVisible"), 0);
C = addCheck(C, 'H4 mask raised by 1e-9 deg -> not visible', 'logical', 1, ...
             double(sinrJust == -Inf && reasonJust == "NotVisible"), 0);

% --- Η5: πεπερασμένες ισοπαλίες SINR δεν εμφανίζονται ---
tieCount = 0; nCmp = 0;
for r = 1:40
    rng(100 + r);
    [~,~,~,~,~,~,~,~,~, bsS, ~,~,~,~,~, satS] = ...
        simulateScenario(bs2, uGeo, satOverhead, wgs84, simParameters, satParameters);
    fin = isfinite(bsS) & isfinite(satS);
    tieCount = tieCount + sum(bsS(fin) == satS(fin));
    nCmp = nCmp + sum(fin);
end
C = addCheck(C, 'H5 finite ties in a sample of comparisons', 'count', 0, tieCount, 0);

% --- D6-D8: παρεμβολή μεταξύ δορυφόρων ---
% FRF 1: SINR < SNR· FRF 3: SINR = SNR· δορυφόρος κάτω από τη μάσκα δεν παρεμβάλλει.
satPair  = [satOverhead; baseLat + 200e3/111320, baseLon, 600e3];
satHidden = [satOverhead; baseLat + 3000e3/111320, baseLon, 600e3];
rng(31);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~, satPlPair, satSinrPair] = ...
    simulateScenario(bs1, uGeo, satPair, wgs84, simParameters, satParameters);
C = addCheck(C, 'D6 two visible satellites, FRF 1: SINR < SNR', 'logical', 1, ...
             double(all(satSinrPair < (satParameters.EIRP - 30) - satPlPair - noiseSatDbw)), 0);

sat3 = satParameters;
sat3.FreqReuse = 3;
sat3.Bandwidth = satParameters.Bandwidth / 3;
sat3.EIRP      = sat3.EirpDensityDbwPerMHz + 10*log10(sat3.Bandwidth/1e6) + 30;
sat3.TxPower   = sat3.EIRP - sat3.AntennaGain;
noiseSat3Dbw   = 10*log10(kBoltz * teqK * sat3.Bandwidth);
rng(31);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~, satPl3, satSinr3] = ...
    simulateScenario(bs1, uGeo, satPair, wgs84, simParameters, sat3);
C = addCheck(C, 'D7 two visible satellites, FRF 3: SINR == SNR', 'dB', 0, ...
             max(abs(satSinr3 - ((sat3.EIRP - 30) - satPl3 - noiseSat3Dbw))), 1e-9);

rng(31);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~, satPlHid, satSinrHid] = ...
    simulateScenario(bs1, uGeo, satHidden, wgs84, simParameters, satParameters);
C = addCheck(C, 'D8 satellite below the mask: SINR == SNR', 'dB', 0, ...
             max(abs(satSinrHid - ((satParameters.EIRP - 30) - satPlHid - noiseSatDbw))), 1e-9);

%% ================== Θ. Κεραία τομέα ==================
% Αναμενόμενες τιμές από τους τύπους των πινάκων, ανεξάρτητα από τη bsSectorGainDb.
ant = struct('ElementGainDbi', 8, 'Theta3dBDeg', 65, 'Phi3dBDeg', 65, 'SlaVDb', 30, ...
    'AmDb', 30, 'NumVertElements', 10, 'VertSpacingLambda', 0.5, 'ElecTiltDeg', 102);
aVref = @(th) -min(12*((th - 90)/65).^2, 30);
aHref = @(ph) -min(12*(ph/65).^2, 30);
afRef = @(th) afFromWeights(th, 10, 0.5, 102);
gRef  = @(th, ph) 8 + (-min(-(aVref(th) + aHref(ph)), 30)) + 10*log10(afRef(th));

C = addCheck(C, 'TH1 array maximum at theta_etilt = M (10 dB)', 'dB', 10, ...
             10*log10(afRef(102)), 1e-9);
C = addCheck(C, 'TH2 gain at (102°, 0°) = 8 + A_V(102°) + 10', 'dBi', ...
             18 - 12*(12/65)^2, bsSectorGainDb(102, 0, ant), 1e-9);
C = addCheck(C, 'TH3 horizontal half-width: phi = 32.5° -> -3 dB', 'dB', -3, ...
             bsSectorGainDb(102, 32.5, ant) - bsSectorGainDb(102, 0, ant), 1e-9);
C = addCheck(C, 'TH4 back lobe: A = -A_m = -30 dB', 'dBi', 8 - 30 + 10*log10(afRef(102)), ...
             bsSectorGainDb(102, 180, ant), 1e-9);
thGrid = 90:0.5:150;  phGrid = -180:5:180;
[TH, PH] = ndgrid(thGrid, phGrid);
C = addCheck(C, 'TH5 (theta, phi) grid: same as the table formulas', 'dB', 0, ...
             max(abs(bsSectorGainDb(TH, PH, ant) - gRef(TH, PH)), [], 'all'), 1e-9);
C = addCheck(C, 'TH6 symmetry phi -> -phi', 'dB', 0, ...
             max(abs(bsSectorGainDb(TH, PH, ant) - bsSectorGainDb(TH, -PH, ant)), [], 'all'), 1e-12);

% Ένα site, 3 τομείς: SINR = max_c S_c / (N + Σ_{c'≠c} S_c'),  S_c = P + G_c - PL.
spS = simParameters;
spS.CellSite = [1; 1; 1];
spS.CellBoresightDeg = [30; 150; 270];
spS.BsAntenna = ant;
siteGeo = [baseLat baseLon 25];
uAz = [30 75 150 200 270 330 10 110];   % από ανατολή, αριστερόστροφα
uDs = [120 300 450 200 800 60 950 350];
uSec = zeros(numel(uAz),3);
for u = 1:numel(uAz)
    [la, lo] = reckon(baseLat, baseLon, uDs(u), 90 - uAz(u), wgs84);
    uSec(u,:) = [la lo 1.5];
end
rng(11);
[~,~,~,~,~,~,~,~,~, secSinr, ~, secPl] = simulateScenario(repmat(siteGeo,3,1), uSec, ...
    satOverhead, wgs84, spS, satParameters);
[d2dS, azS] = distance(baseLat, baseLon, uSec(:,1), uSec(:,2), wgs84);
thS  = 90 + atand((25 - 1.5) ./ d2dS);
phS  = 90 - azS;
sExp = zeros(numel(uAz),3);
for c = 1:3
    phC = mod(phS - spS.CellBoresightDeg(c) + 180, 360) - 180;
    sExp(:,c) = 10.^(((simParameters.TxPower - 30) + gRef(thS, phC) - secPl)/10);
end
sinrExp = 10*log10(max(sExp, [], 2) ./ (10^(noiseBsDbw/10) + sum(sExp,2) - max(sExp, [], 2)));
C = addCheck(C, 'TH7 3-sector site: SINR == max S_c/(N + other sectors)', 'dB', 0, ...
             max(abs(secSinr - sinrExp)), 1e-9);
[~, bestSec] = max(sExp, [], 2);
C = addCheck(C, 'TH8 UE on a sector boresight -> that sector', 'logical', 1, ...
             double(isequal(bestSec([1 3 5])', [1 2 3])), 0);

%% ================== Ι. Soft LOS και RLF ==================
% Soft LOS (TR 38.901 εξ. 7.6-18): P_LOS = 1 -> LOS_soft = 1· LOS_soft > 1/2 <=> LOS.
spSoft = simParameters; spSoft.SoftLos = true;
satNoS = satParameters; satNoS.MinElevationDeg = 95;
rng(17);
uNear = [baseLat + 15/111320, baseLon, 1.5];   % d2D = 15 m < 18 m: P_LOS = 1 (UMa)
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, chNear] = simulateScenario([baseLat baseLon 25], ...
    uNear, satOverhead, wgs84, spSoft, satNoS);
C = addCheck(C, 'I1 P_LOS = 1 -> LOS_soft = 1', '-', 1, chNear.LosSoft(1), 0);
nU = 200;
uMany = zeros(nU,3);
for u = 1:nU
    [la, lo] = reckon(baseLat, baseLon, 20 + 1500*rand(), 360*rand(), wgs84);
    uMany(u,:) = [la lo 1.5];
end
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, chMany] = simulateScenario([baseLat baseLon 25], ...
    uMany, satOverhead, wgs84, spSoft, satNoS);
okS = ~isnan(chMany.LosSoft(:,1));
C = addCheck(C, 'I2 LOS_soft in [0, 1]', 'logical', 1, ...
             double(all(chMany.LosSoft(okS,1) >= 0 & chMany.LosSoft(okS,1) <= 1)), 0);
C = addCheck(C, 'I2 LOS_soft > 1/2 <=> hard LOS decision', 'logical', 1, ...
             double(isequal(chMany.LosSoft(okS,1) > 0.5, chMany.IsLOS(okS,1))), 0);

% RLF (TS 38.331 §5.3.10) με Δt = 1 s (T310 = 1 βήμα)· πτώση ισχύος 60 dB
% φέρνει το SINR κάτω από το Qout.
spR = simParameters;
spR.Mobility = struct('HysteresisDb', 3, 'TimeToTriggerMs', 480, 'DtSeconds', 1, ...
    'PingPongMtsS', 1, 'MlbMaxOffsetDb', 6, ...
    'Rlf', struct('QoutDb', -8, 'QinDb', -6, 'T310Ms', 1000, 'N310', 1, 'N311', 1));
spLow = spR; spLow.EIRP = spR.EIRP - 60;
uR = [baseLat + 300/111320, baseLon, 1.5];
bsR = [baseLat baseLon 25];
rng(19);
[n1,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~, ch1,~,~,~,~,~,~, d1] = simulateScenario(bsR, uR, ...
    satOverhead, wgs84, spR, satNoS);
[n2,t2,~,~,s2,c2,~,~,~,~,~,~,~,~,~,~, ch2,~,~,thr2,~,~,k2, d2] = simulateScenario(bsR, uR, ...
    satOverhead, wgs84, spLow, satNoS, ch1, d1);
C = addCheck(C, 'I3 SINR < Qout: the UE stays on the gNB', 'logical', 1, ...
             double(n1 == "BS1" && n2 == "BS1" && t2 == "Terrestrial" && s2 < -8), 0);
C = addCheck(C, 'I3 T310 starts (0 steps)', 'steps', 0, d2.Actual.T310, 0);
C = addCheck(C, 'I3 below MCS 0: zero rate', 'Mbps', 0, thr2 + c2, 0);
C = addCheck(C, 'I3 no RLF yet', '1/min', 0, k2.Actual.RlfRate_perUserMin, 0);
[n3,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,k3] = simulateScenario(bsR, uR, ...
    satOverhead, wgs84, spLow, satNoS, ch2, d2);
C = addCheck(C, 'I4 T310 expiry without in-sync -> RLF, no link', 'logical', 1, ...
             double(n3 == "None"), 0);
C = addCheck(C, 'I4 RLF per UE per minute (1 in 1 s)', '1/min', 60, ...
             k3.Actual.RlfRate_perUserMin, 1e-12);
C = addCheck(C, 'I4 an RLF is not a handover', '1/min', 0, k3.Actual.HoRate_perUserMin, 0);
[n3b,~,~,~,s3b,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,k3b, d3b] = simulateScenario(bsR, uR, ...
    satOverhead, wgs84, spR, satNoS, ch2, d2);
C = addCheck(C, 'I5 SINR > Qin before expiry -> stays, T310 stops', 'logical', 1, ...
             double(n3b == "BS1" && s3b > -6 && d3b.Actual.T310 == -1 && ...
                    k3b.Actual.RlfRate_perUserMin == 0), 0);

%% ------------------ Πίνακας και εκτύπωση ------------------
T = table(string({C.name}'), string({C.unit}'), [C.expected]', [C.actual]', ...
    abs([C.actual]' - [C.expected]'), [C.tol]', logical([C.pass]'), ...
    'VariableNames', {'Check','Unit','Expected','Actual','AbsError','Tolerance','Pass'});

fprintf('\n=== Systematic link-budget verification ===\n');
fprintf('Reference scenario: UMa, fc = %.2f GHz, %d gNBs, %d UEs, satellite at %g km\n\n', ...
    simParameters.CarrierFrequency/1e9, size(bs2,1), nU, satOverhead(3)/1e3);
fprintf('%-52s %14s %14s   %s\n', 'Check', 'Expected', 'Measured', 'Result');
for k = 1:height(T)
    fprintf('%-52s %14.6g %14.6g   %s\n', T.Check(k), T.Expected(k), T.Actual(k), ...
        boolText(T.Pass(k)));
end

nPass = sum(T.Pass); nAll = height(T);
fprintf('\nTotal: %d/%d checks passed.\n', nPass, nAll);
fprintf('Candidate comparisons checked for ties: %d, finite ties: %d\n', nCmp, tieCount);
fprintf('Terrestrial octave 500->1000 m: %.2f dB (free space: %.2f dB)\n\n', ...
    umaOctaveDb, 20*log10(2));

%% ------------------ Γράφημα ------------------
Cc = thesisColors();
fig = thesisFigure(17, 6.8);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);
dPlot = logspace(2,6,60);
semilogx(ax, dPlot, arrayfun(@(d) fspl(d, lambdaSat), dPlot), '-', 'Color', Cc.cat(1,:)); hold(ax, 'on');
semilogx(ax, dPlot, fspl(dPlot(1), lambdaSat) + 20*log10(dPlot/dPlot(1)), '--', 'Color', Cc.cat(2,:));
grid(ax, 'on'); xlabel(ax, 'Distance [m]'); ylabel(ax, 'Path loss [dB]');
legend(ax, {'FSPL (model)','20 dB per decade'}, 'Location','southeast');
title(ax, '(a) Free space: +6.02 dB per octave', 'FontWeight', 'normal');

ax = nexttile(tl);
plot(ax, offsets, sinrShift1, '-o', 'Color', Cc.cat(1,:), 'MarkerFaceColor', Cc.cat(1,:), 'MarkerSize', 4); hold(ax, 'on');
plot(ax, offsets, sinrShift2, '-s', 'Color', Cc.cat(2,:), 'MarkerFaceColor', Cc.cat(2,:), 'MarkerSize', 4);
plot(ax, offsets, offsets, '--', 'Color', Cc.text, 'LineWidth', 0.9);
grid(ax, 'on'); xlabel(ax, 'Change of transmit power [dB]'); ylabel(ax, 'Change of SINR [dB]');
legend(ax, {'single gNB (noise-limited)','two gNBs (interference)','slope 1'}, 'Location','northwest');
title(ax, '(b) Transmit power vs SINR', 'FontWeight', 'normal');

ax = nexttile(tl);
b = bar(ax, [nodeBwPerNode(:)/1e6, bwSumPerNode(:)/1e6], 'EdgeColor', 'w');
b(1).FaceColor = Cc.cat(1,:); b(2).FaceColor = Cc.cat(2,:);
grid(ax, 'on'); ylabel(ax, 'Bandwidth [MHz]');
set(ax, 'XTick', 1:numel(nodeLabels), 'XTickLabel', nodeLabels);
ylim(ax, [0 1.3*max(nodeBwPerNode(:)/1e6)]);
legend(ax, {'available','sum of user shares'}, 'Location','north', 'Orientation','horizontal');
title(ax, '(c) Bandwidth conservation', 'FontWeight', 'normal');
title(tl, sprintf('Link-budget verification: %d/%d checks passed', nPass, nAll), ...
    'FontWeight', 'bold', 'FontSize', 10);

pngPath = fullfile(outputDir, 'link_budget_validation.png');
saveThesisFigure(fig, pngPath);

fprintf('Figure -> %s\n', pngPath);

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
runParams.results = T;
saveRunVersion('linkBudgetValidation', runParams, {pngPath}, label);

if nPass < nAll
    error('linkBudgetValidation:Failed', ...
        '%d of %d link-budget checks failed.', nAll - nPass, nAll);
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

function af = afFromWeights(thetaDeg, K, dVLambda, tiltDeg)
% Παράγοντας συστοιχίας από τα βάρη w_m (TR 36.873 Πίν. 7.1-1), με βρόχο.
af = zeros(size(thetaDeg));
for i = 1:numel(thetaDeg)
    acc = 0;
    for m = 1:K
        w   = exp(-1j*2*pi*(m-1)*dVLambda*cosd(tiltDeg)) / sqrt(K);
        acc = acc + w * exp(1j*2*pi*(m-1)*dVLambda*cosd(thetaDeg(i)));
    end
    af(i) = abs(acc)^2;
end
end

function s = boolText(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end
