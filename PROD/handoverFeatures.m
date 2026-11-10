function x = handoverFeatures(c)
%HANDOVERFEATURES Χαρακτηριστικά ενός handover τη στιγμή της απόφασης.
%
%   names = handoverFeatures()    ονόματα των στηλών, με τη σειρά τους
%   x     = handoverFeatures(c)   διάνυσμα γραμμή, μία τιμή ανά στήλη
%
% Κοινή για την καταγραφή (runSimulation/logHandoverEvents, εκπαίδευση) και την
% απόφαση ML (simulateScenario/mlGateDecision), ώστε το μοντέλο να βλέπει ό,τι
% είδε στην εκπαίδευση.
%
% Πεδία του c (From = εξυπηρετών, To = στόχος) και πηγές:
%   SinrFrom, SinrTo  [1x3] μέτρηση L3 (TS 38.331 §5.5.3.2) στα t, t-1 s, t-2 s
%                     -> SINR, gap, μεταβολή 1 s, μέσος όρος/κλίση 3 s
%                     (Aggarwal et al., IEEE ICSigSys 2026)
%   LagS              απόσταση δειγμάτων (s)
%   ReturnToPrev, StayFrom_s   UE History Information (TS 38.413 §9.3.1.97)
%   LoadFrom/To, UsersOnFrom/To  Resource Status (TS 38.423 §9.1.3.21)
%   ElevFrom/To, DistFrom/To     γεωμετρία (TR 37.817 §5.3.2.4)
%   NumUsers, MlbOffset          πλήθος χρηστών, μετατόπιση MLB του A3
names = {'NumUsers','ReturnToPrev','StayFrom_s', ...
    'SinrServ_dB','SinrCand_dB','Gap_dB','DeltaSinr1s_dB','DeltaGap1s_dB', ...
    'AvgSinr3s_dB','AvgGap3s_dB','SlopeSinr3s_dBps','SlopeGap3s_dBps', ...
    'ElevFrom_deg','ElevTo_deg','DistFrom_m','DistTo_m', ...
    'LoadFrom','LoadTo','UsersOnFrom','UsersOnTo','MlbOffset_dB'};
if nargin == 0
    x = names;
    return;
end
sinrK = c.SinrFrom;
gapK  = c.SinrTo - c.SinrFrom;
sinrK(~isfinite(sinrK)) = NaN;
gapK(~isfinite(gapK))   = NaN;
x = [c.NumUsers, c.ReturnToPrev, c.StayFrom_s, ...
    sinrK(1), sinrK(1) + gapK(1), gapK(1), ...
    sinrK(1) - sinrK(2), gapK(1) - gapK(2), ...
    mean(sinrK), mean(gapK), ...
    (sinrK(1) - sinrK(3))/(2*c.LagS), (gapK(1) - gapK(3))/(2*c.LagS), ...
    c.ElevFrom, c.ElevTo, c.DistFrom, c.DistTo, ...
    c.LoadFrom, c.LoadTo, c.UsersOnFrom, c.UsersOnTo, c.MlbOffset];
end
