function [S, stepKPI] = apply_action_and_serve(S, P, L, action)
%APPLY_ACTION_AND_SERVE Serve URLLC + eMBB under HO interruptions.
% URLLC modes: 0 TN, 1 NTN, 2 SPLIT, 3 DUP
% eMBB modes:  0 TN, 1 NTN, 2 SPLIT (aggregate)

stepKPI = struct('URLLC_delivered',0,'URLLC_delay',NaN,'URLLC_success',0,...
                 'eMBB_bitsDL',0,'eMBB_bitsUL',0);

% Per-slot capacities (serving links)
capTNdl  = L.R_TN_DL  * P.dt;
capTNul  = L.R_TN_UL  * P.dt;
capNTNdl = L.R_NTN_DL * P.dt;
capNTNul = L.R_NTN_UL * P.dt;

% HO interruption: if active, that link has 0 capacity
if S.HO_TN.active
    capTNdl = 0; capTNul = 0;
end
if S.HO_NTN.active
    capNTNdl = 0; capNTNul = 0;
end

% -----------------------
% URLLC (DL, packet-level)
% -----------------------
% Efficient: serve only the head-of-line packet (FIFO)

% Safety: if tail doesn't exist, fall back (but you should have it)
if ~isfield(S.URLLC,'tail')
    S.URLLC.tail = numel(S.URLLC.queue);
end
if ~isfield(S.URLLC,'head')
    S.URLLC.head = 1;
end

h = S.URLLC.head;

% Move head forward if already delivered/expired
while (h <= S.URLLC.tail) && (S.URLLC.queue(h).delivered || S.URLLC.queue(h).expired)
    h = h + 1;
end
S.URLLC.head = h;

idx = [];
if h <= S.URLLC.tail
    idx = h;
end

% Serve at most 1 packet per slot
if ~isempty(idx)
    pktBits = S.URLLC.pktBits;

    switch action.modeU
        case 0 % TN
            if capTNdl >= pktBits
                delay = (S.t - S.URLLC.queue(idx).genTime) + L.Dprop_TN;
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        case 1 % NTN
            if capNTNdl >= pktBits
                delay = (S.t - S.URLLC.queue(idx).genTime) + L.Dprop_NTN;
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        case 2 % SPLIT (conservative: require both halves)
            if (capTNdl/2 >= pktBits/2) && (capNTNdl/2 >= pktBits/2)
                delay = (S.t - S.URLLC.queue(idx).genTime) + max(L.Dprop_TN, L.Dprop_NTN);
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end

        case 3 % DUP (either link can deliver)
            if (capTNdl >= pktBits) || (capNTNdl >= pktBits)
                delay = (S.t - S.URLLC.queue(idx).genTime) + min(L.Dprop_TN, L.Dprop_NTN);
                [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI);
            end
    end
end

% -----------------------
% Deadline timeout check (efficient FIFO expiry)
% -----------------------
while S.URLLC.head <= S.URLLC.tail
    if S.URLLC.queue(S.URLLC.head).delivered || S.URLLC.queue(S.URLLC.head).expired
        S.URLLC.head = S.URLLC.head + 1;
        continue;
    end

    if (S.t - S.URLLC.queue(S.URLLC.head).genTime) > P.URLLC.deadline
        S.URLLC.queue(S.URLLC.head).expired = true;
        S.URLLC.failCount = S.URLLC.failCount + 1;
        S.URLLC.head = S.URLLC.head + 1;
    else
        break;
    end
end

% -----------------------
% eMBB (DL/UL throughput)
% -----------------------
switch action.modeE
    case 0 % TN
        stepKPI.eMBB_bitsDL = capTNdl;
        stepKPI.eMBB_bitsUL = capTNul;

    case 1 % NTN
        stepKPI.eMBB_bitsDL = capNTNdl;
        stepKPI.eMBB_bitsUL = capNTNul;

    case 2 % SPLIT aggregate
        stepKPI.eMBB_bitsDL = capTNdl + capNTNdl;
        stepKPI.eMBB_bitsUL = capTNul + capNTNul;
end

end

function [S, stepKPI] = deliver_urlcc(S, P, idx, delay, stepKPI)

S.URLLC.queue(idx).delivered = true;
S.URLLC.queue(idx).delay = delay;

ok = (delay <= P.URLLC.deadline);
if ok
    S.URLLC.succCount = S.URLLC.succCount + 1;
else
    S.URLLC.failCount = S.URLLC.failCount + 1;
end

stepKPI.URLLC_delivered = 1;
stepKPI.URLLC_delay = delay;
stepKPI.URLLC_success = ok;

end