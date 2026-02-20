function S = traffic_step(S, P, k)
%TRAFFIC_STEP Generate URLLC packets (DL+UL) and keep eMBB backlogged

t = (k-1) * P.dt;

% URLLC periodic arrivals (DL & UL)
if t >= S.URLLC.nextGenTime
    % DL packet
    S.qDL_URLLC = S.qDL_URLLC + P.URLLC.pktBits;
    S.URLLC.inflight(end+1) = struct('genTime', t, 'bits', P.URLLC.pktBits, ...
                                     'delivered', false, 'delay', NaN);
    % UL packet
    S.qUL_URLLC = S.qUL_URLLC + P.URLLC.pktBits;

    S.URLLC.nextGenTime = S.URLLC.nextGenTime + P.URLLC.period;
end

% eMBB backlogged queues (DL & UL)
if P.eMBB.backlogged
    S.qDL_eMBB = 1e12; % effectively infinite bits
    S.qUL_eMBB = 1e12;
end

end
