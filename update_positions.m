function S = update_positions(S, P, k)
%UPDATE_POSITIONS Update UE and satellites

t = (k-1) * P.dt;

% UE moves right, bounces at ends
S.UE.x = S.UE.x + S.UE.v * P.dt;
if S.UE.x < 0
    S.UE.x = -S.UE.x;
    S.UE.v = -S.UE.v;
elseif S.UE.x > P.areaLen
    S.UE.x = 2*P.areaLen - S.UE.x;
    S.UE.v = -S.UE.v;
end

% Satellite ground-track proxy motion 
for i = 1:numel(P.SAT)
    S.SATx(i) = mod(P.SAT(i).x0 + P.SAT(i).v * t, P.areaLen);
end

end
