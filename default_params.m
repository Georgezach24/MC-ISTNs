
%Συνάρτηση που θέτει τις παραμέτρους της προσομοίωσης.
function P = default_params(scenario)

% --------- Προσομοίωση -----------------------------------------------
P.dt = 1e-3; % 1ms (Βήμα προσομοίωσης).
simTime = 20; % 20s (Συνολικός χρόνος προσομοίωσης).
P.T = round(simTime / P.dt); % Συνολικός αριθμός time slots.
P.areaLen = 20000; % 20km (Μήκος περιοχής προσομοίωσης). 
P.c = 3e8; % Ταχύτητα φωτός.
P.W = 20e6; % 20ΜΗz (Bandwidth).
P.N0 = 1e-20;% Φασματική πυκνότητα θορύβου.
%----------------------------------------------------------------------

% --------- UE - User Equipment ---------------------------------------
P.UE.v = 20; % 20m/s (Ταχύτητα Χρήστη).
P.UE.x0 = P.areaLen/2; % Αρχική θέση χρήστη (Στο κέντρο περιοχής).
%----------------------------------------------------------------------

% --------- BS - Terrestrial Base Stations ----------------------------
%Ορισμός 2 Base stations στα άκρα της περιοχής προσομοίωσης.
P.BS(1).x = 0;
P.BS(2).x = P.areaLen;
%----------------------------------------------------------------------

% --------- Satellites ------------------------------------------------
%Αρχικές οριζόντιες θέσεις των δορυφόρων.
P.SAT(1).x0 = 0.25*P.areaLen;
P.SAT(2).x0 = 0.75*P.areaLen;
%----------------------------------------------------------------------
%Για κάθε δορυφόρο.
for i=1:2
    P.SAT(i).v = 50; % Ταχύτητα (50m/s).
    P.SAT(i).alt = 500e3; % Υψόμετρο (500km).
end
%----------------------------------------------------------------------

% --------- Channel ---------------------------------------------------
% Terestrial Network:
P.TN.alpha = 3.5; % Pathloss exponent. [Urban τιμές 3-4]
P.TN.Pdl = 1; % Εκπεμπόμενη ισχύς DL.
P.TN.Pul = 1; % Εκπεμπόμενη ισχύς UL.
P.TN.fading = 'rayleigh'; % Small-scale fading.[Rayleigh:NLOS environment.]

% Non Terrestrial Network:
P.NTN.alpha = 2.2; % Pathloss exponent.[Μικρότερος λόγο LoS]
P.NTN.Pdl = 1; % Εκπεμπόμενη ισχύς DL.
P.NTN.Pul = 1; % Εκπεμπόμενη ισχύς UL.
P.NTN.fading = 'nakagami'; % nakagami fading. [Καλύτερο για satellites]
P.NTN.nak_m = 2; % Παράμετρος nakagami. [>1 :λιγότερο fading]
P.NTN.extraLoss_dB = 3; %Επιπλέον απώλειες. [Ατμοσφερικές, shadowing, polarization]

P.SINRmin = 10^(-5/10); % Κατώφλι SINR (-5dB).
%----------------------------------------------------------------------

% --------- Handover --------------------------------------------------
P.HO.hyst_dB = 2; % Hysteresis margin (2 dB).
P.HO.hyst = 10^(P.HO.hyst_dB/10); % Μετατροπή σε γραμμική κλίμακα.
P.HO.TTT_TN_s = 0.02; % Time-To-Trigger TN (20 ms).
P.HO.TTT_NTN_s = 0.08; % Time-To-Trigger NTN (80 ms).
P.HO.exec_TN_s = 0.005; % Execution interruption TN (5 ms).
P.HO.exec_NTN_s = 0.01; % Execution interruption NTN (10 ms).
P.HO.alpha = 0.95; % Φίλτρο εξομάλυνσης RSRP.
%----------------------------------------------------------------------

% --------- URLLC -----------------------------------------------------
P.URLLC.pktSize = 64*8; % Μέγεθος πακέτου (512bits).
P.URLLC.lambda = 100; % Arrival rate (100 packets/sec).

% Deadline:
if scenario == 'A'
    P.URLLC.deadline = 0.02;
else
    P.URLLC.deadline = 0.01;
end
%----------------------------------------------------------------------

end