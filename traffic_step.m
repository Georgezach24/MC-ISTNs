function S = traffic_step(S, P, k)
%TRAFFIC_STEP Generate URLLC packets (1 per 10ms) into preallocated queue

t = (k-1) * P.dt;

period = 10e-3;        % 10 ms -> 100 pkt/s
pktBits = 64*8;

S.URLLC.pktBits = pktBits;

if t >= S.URLLC.nextGenTime
    S.URLLC.genCount = S.URLLC.genCount + 1;

    S.URLLC.tail = S.URLLC.tail + 1;
    j = S.URLLC.tail;

    % Safety: in case maxPkts was underestimated
    if j > S.URLLC.maxPkts
        error('URLLC.maxPkts exceeded. Increase margin in init_scenario.');
    end

    S.URLLC.queue(j).genTime = t;
    S.URLLC.queue(j).delivered = false;
    S.URLLC.queue(j).expired = false;
    S.URLLC.queue(j).delay = NaN;

    S.URLLC.nextGenTime = S.URLLC.nextGenTime + period;
end

end