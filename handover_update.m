function [S, hoTN, hoNTN] = handover_update(S, P, L)
%HANDOVER_UPDATE Network-controlled HO for TN and NTN using:
% - filtered SINR metric (EWMA) to avoid fast-fading instability
% - hysteresis + TTT
% - execution (interruption) timer

hoTN = false; 
hoNTN = false;

% ---------------------------------------------------------
% 0) Update filtered metrics (RSRP/RSRQ-like filtering)
% ---------------------------------------------------------
% EWMA: filt = alpha*filt + (1-alpha)*instant
S.HO_TN.filtSINR  = P.HO.alpha * S.HO_TN.filtSINR  + (1 - P.HO.alpha) * L.TN.SINR_DL;
S.HO_NTN.filtSINR = P.HO.alpha * S.HO_NTN.filtSINR + (1 - P.HO.alpha) * L.NTN.SINR_DL;

% Use filtered metrics for HO decisions
TNm  = S.HO_TN.filtSINR;
NTNm = S.HO_NTN.filtSINR;

% ---------------------------------------------------------
% 1) TN HO
% ---------------------------------------------------------
if ~S.HO_TN.active
    % Best candidate based on filtered SINR
    [bestMetric, bestBS] = max(TNm);
    servMetric = TNm(S.servBS);

    if bestBS ~= S.servBS && bestMetric > servMetric * P.HO.hyst
        S.HO_TN.ttt_s = S.HO_TN.ttt_s + P.dt;
        S.HO_TN.cand = bestBS;

        if S.HO_TN.ttt_s >= P.HO.TTT_TN_s
            % Start execution
            S.HO_TN.active  = true;
            S.HO_TN.timer_s = P.HO.exec_TN_s;
            S.HO_TN.ttt_s   = 0;
        end
    else
        % Reset TTT if condition not continuously satisfied
        S.HO_TN.ttt_s = 0;
        S.HO_TN.cand  = S.servBS;
    end
else
    % Execution phase
    S.HO_TN.timer_s = S.HO_TN.timer_s - P.dt;
    if S.HO_TN.timer_s <= 0
        % Complete HO
        S.servBS = S.HO_TN.cand;
        S.HO_TN.active  = false;
        S.HO_TN.timer_s = 0;
        hoTN = true;
    end
end

% ---------------------------------------------------------
% 2) NTN HO
% ---------------------------------------------------------
if ~S.HO_NTN.active
    [bestMetric, bestSAT] = max(NTNm);
    servMetric = NTNm(S.servSAT);

    if bestSAT ~= S.servSAT && bestMetric > servMetric * P.HO.hyst
        S.HO_NTN.ttt_s = S.HO_NTN.ttt_s + P.dt;
        S.HO_NTN.cand = bestSAT;

        if S.HO_NTN.ttt_s >= P.HO.TTT_NTN_s
            S.HO_NTN.active  = true;
            S.HO_NTN.timer_s = P.HO.exec_NTN_s;
            S.HO_NTN.ttt_s   = 0;
        end
    else
        S.HO_NTN.ttt_s = 0;
        S.HO_NTN.cand  = S.servSAT;
    end
else
    S.HO_NTN.timer_s = S.HO_NTN.timer_s - P.dt;
    if S.HO_NTN.timer_s <= 0
        S.servSAT = S.HO_NTN.cand;
        S.HO_NTN.active  = false;
        S.HO_NTN.timer_s = 0;
        hoNTN = true;
    end
end

end