function p = handoverGatePredict(M, x)
%HANDOVERGATEPREDICT Πιθανότητα ping-pong από το Random Forest του
% Model/train_handover_gate.py, για διάνυσμα χαρακτηριστικών x.
%
% Ίδιος κανόνας με το scikit-learn: NaN -> M.Sentinel, σύγκριση σε float32,
% αριστερά αν x(f) <= κατώφλι, πιθανότητα = μέσος όρος των φύλλων.
%
% M (από loadHandoverGate), διανύσματα στήλη:
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
