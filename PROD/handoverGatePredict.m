function p = handoverGatePredict(M, x)
%HANDOVERGATEPREDICT P(ping-pong) από το Random Forest του train_handover_gate.py.
% Ίδιο με το scikit-learn: NaN -> Sentinel, float32, αριστερά αν x <= κατώφλι,
% μέσος όρος φύλλων. M από loadHandoverGate (Feature/Left/Right = 0 στα φύλλα).
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
