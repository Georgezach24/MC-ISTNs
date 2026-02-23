function [S, stepKPI] = apply_action_and_serve(S, P, L, action)
%APPLY_ACTION_AND_SERVE Serve queues based on action; compute delivered bits
% action.modeU: 0 TN, 1 NTN, 2 SPLIT, 3 DUP (URLLC)
% action.modeE: 0 TN, 1 NTN, 2 SPLIT (eMBB)

stepKPI = struct('URLLC_delivered',0,'URLLC_delay',NaN,'URLLC_success',0,...
                 'eMBB_bitsDL',0,'eMBB_bitsUL',0);

% Available DL/UL rates
RTNdl = L.R_TN_DL;   RNTNdl = L.R_NTN_DL;
RTNul = L.R_TN_UL;   RNTNul = L.R_NTN_UL;

% Convert rates to bits per slot
capTNdl  = RTNdl  * P.dt;
capNTNdl = RNTNdl * P.dt;
capTNul  = RTNul  * P.dt;
capNTNul = RNTNul * P.dt;

if S.HO_TN.active
    capTNdl = 0; capTNul = 0;
end
if S.HO_NTN.active
    capNTNdl = 0; capNTNul = 0;
end

% --- URLLC DL serving ---
if S.qDL_URLLC > 0
    switch action.modeU
        case 0 % TN
            served = min(S.qDL_URLLC, capTNdl);
            S.qDL_URLLC = S.qDL_URLLC - served;
            % if served completes a packet, mark delivery with TN delay
            [S, stepKPI] = deliver_urllc_if_packet_done(S, P, served, L.Dprop_TN, stepKPI);

        case 1 % NTN
            served = min(S.qDL_URLLC, capNTNdl);
            S.qDL_URLLC = S.qDL_URLLC - served;
            [S, stepKPI] = deliver_urllc_if_packet_done(S, P, served, L.Dprop_NTN, stepKPI);

        case 2 % SPLIT 50/50
            servedTN  = min(S.qDL_URLLC, capTNdl/2);
            servedNTN = min(S.qDL_URLLC - servedTN, capNTNdl/2);
            served = servedTN + servedNTN;
            S.qDL_URLLC = S.qDL_URLLC - served;
            % simplistic: assume delivery when full pkt bits served, delay = max(path delays)
            [S, stepKPI] = deliver_urllc_if_packet_done(S, P, served, max(L.Dprop_TN, L.Dprop_NTN), stepKPI);

        case 3 % DUP (send same packet on both; consume both caps)
            % Only meaningful when at least one full pkt can be sent
            servedTN  = min(S.qDL_URLLC, capTNdl);
            servedNTN = min(S.qDL_URLLC, capNTNdl);
            served = max(servedTN, servedNTN); % delivered if either path completes
            S.qDL_URLLC = S.qDL_URLLC - served;
            [S, stepKPI] = deliver_urllc_if_packet_done(S, P, served, min(L.Dprop_TN, L.Dprop_NTN), stepKPI);
        otherwise
            error('Unknown URLLC mode');
    end
end

% --- URLLC UL serving (mirrors DL; just account bits) ---
if S.qUL_URLLC > 0
    switch action.modeU
        case 0
            served = min(S.qUL_URLLC, capTNul); S.qUL_URLLC = S.qUL_URLLC - served;
        case 1
            served = min(S.qUL_URLLC, capNTNul); S.qUL_URLLC = S.qUL_URLLC - served;
        case 2
            servedTN  = min(S.qUL_URLLC, capTNul/2);
            servedNTN = min(S.qUL_URLLC - servedTN, capNTNul/2);
            served = servedTN + servedNTN; S.qUL_URLLC = S.qUL_URLLC - served;
        case 3
            served = min(S.qUL_URLLC, max(capTNul, capNTNul)); S.qUL_URLLC = S.qUL_URLLC - served;
    end
end

% --- eMBB DL serving ---
switch action.modeE
    case 0 % TN
        bits = capTNdl;
    case 1 % NTN
        bits = capNTNdl;
    case 2 % SPLIT
        bits = capTNdl + capNTNdl; % optimistic aggregation for start
    otherwise
        error('Unknown eMBB mode');
end
stepKPI.eMBB_bitsDL = bits;

% --- eMBB UL serving ---
switch action.modeE
    case 0
        bits = capTNul;
    case 1
        bits = capNTNul;
    case 2
        bits = capTNul + capNTNul;
end
stepKPI.eMBB_bitsUL = bits;

% Energy (simple): RX circuitry + UL TX
S.energy_J = S.energy_J + P.energy.rxCircuit_W * P.dt;
% UL TX energy proxy: assume UE transmits if any UL bits served; use max UL power
if (stepKPI.eMBB_bitsUL > 0) || (S.qUL_URLLC > 0)
    % use nominal 23 dBm
    S.energy_J = S.energy_J + P.TN.Pul * P.dt;
end

end

function [S, stepKPI] = deliver_urllc_if_packet_done(S, P, servedBits, propDelay, stepKPI)


if servedBits >= P.URLLC.pktBits
    % deliver earliest undelivered inflight
    idx = find(~[S.URLLC.inflight.delivered], 1, 'first');
    if ~isempty(idx)
        tnow = S.t;
        delay = (tnow - S.URLLC.inflight(idx).genTime) + propDelay;
        S.URLLC.inflight(idx).delivered = true;
        S.URLLC.inflight(idx).delay = delay;

        S.URLLC.rxCount = S.URLLC.rxCount + 1;
        ok = (delay <= P.URLLC.deadline);
        if ok
            S.URLLC.succCount = S.URLLC.succCount + 1;
        else
            S.URLLC.lateCount = S.URLLC.lateCount + 1;
        end
        stepKPI.URLLC_delivered = 1;
        stepKPI.URLLC_delay = delay;
        stepKPI.URLLC_success = ok;
    end
end
end
