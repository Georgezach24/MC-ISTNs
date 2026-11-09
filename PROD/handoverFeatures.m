function x = handoverFeatures(c)
%HANDOVERFEATURES Χαρακτηριστικά ενός handover τη στιγμή της απόφασης.
%
%   names = handoverFeatures()    ονόματα των στηλών, με τη σειρά τους
%   x     = handoverFeatures(c)   διάνυσμα γραμμή, μία τιμή ανά στήλη
%
% Μία συνάρτηση για δύο χρήστες: την καταγραφή των handovers στο
% Dataset/handover_events.mat (runSimulation, logHandoverEvents), από την
% οποία εκπαιδεύεται το μοντέλο, και την απόφαση με ML μέσα στην προσομοίωση
% (simulateScenario, mlGateDecision). Έτσι το μοντέλο βλέπει στην εφαρμογή
% ακριβώς ό,τι είδε στην εκπαίδευση· το runSimulation το επαληθεύει σε κάθε
% handover των πολιτικών *ML.
%
% Πεδία του c (f = κόμβος που αφήνει ο χρήστης, g = στόχος):
%   SinrFrom, SinrTo  [1x3] φιλτραρισμένη μέτρηση L3 (TS 38.331 §5.5.3.2) των
%                     f και g στα δείγματα t, t-1 s, t-2 s· NaN όπου δεν
%                     υπάρχει (αρχή διέλευσης, κόμβος εκτός μέτρησης)
%   LagS              s, απόσταση των δειγμάτων (1 s)
%   NumUsers, ReturnToPrev, StayFrom_s, ElevFrom, ElevTo, DistFrom, DistTo,
%   LoadFrom, LoadTo, UsersOnFrom, UsersOnTo, MlbOffset
% Οι πηγές κάθε στήλης είναι στο σχόλιο της logHandoverEvents (runSimulation).
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
