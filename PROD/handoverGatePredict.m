function p = handoverGatePredict(M, x)
%HANDOVERGATEPREDICT Πιθανότητα ping-pong ενός handover από το δάσος δέντρων
% που εκπαίδευσε το Model/train_handover_gate.py (scikit-learn
% RandomForestClassifier), για ένα διάνυσμα χαρακτηριστικών x.
%
% Ίδιος κανόνας με το scikit-learn, ώστε η πιθανότητα να είναι η ίδια ως το
% σφάλμα στρογγύλευσης του μέσου όρου (το runSimulation το ελέγχει με τα
% δείγματα ελέγχου του αρχείου του μοντέλου πριν από κάθε εκτέλεση):
%   - οι τιμές NaN (δεν ορίζονται: π.χ. ανύψωση για σταθμό βάσης, ιστορικό
%     πριν την αρχή της διέλευσης) αντικαθίστανται με M.Sentinel, όπως στην
%     εκπαίδευση·
%   - το scikit-learn συγκρίνει σε float32, οπότε και εδώ x -> single·
%   - σε κάθε κόμβο αριστερά αν x(f) <= κατώφλι·
%   - η πιθανότητα του δάσους είναι ο μέσος όρος των φύλλων των δέντρων.
%
% M (από loadHandoverGate στο runSimulation), όλα διανύσματα στήλη:
%   Roots      κόμβος-ρίζα κάθε δέντρου (δείκτης στους πίνακες κόμβων)
%   Feature    στήλη του x που ελέγχει ο κόμβος (0 = φύλλο)
%   Threshold  κατώφλι του κόμβου
%   Left/Right παιδιά του κόμβου (0 = φύλλο)
%   Value      πιθανότητα ping-pong στο φύλλο
%   MaxDepth   μέγιστο βάθος των δέντρων
x = x(:);
x(isnan(x)) = M.Sentinel;
x = double(single(x));
node = M.Roots;
for d = 1:M.MaxDepth
    f = M.Feature(node);
    internal = f > 0;
    if ~any(internal)
        break;
    end
    n = node(internal);
    goLeft = x(f(internal)) <= M.Threshold(n);
    node(internal) = M.Left(n).*goLeft + M.Right(n).*~goLeft;
end
p = mean(M.Value(node));
end
