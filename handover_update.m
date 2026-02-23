function [S, hoTN, hoNTN] = handover_update(S, P, L)


hoTN = false; hoNTN = false;

% --- TN HO ---
if ~S.HO_TN.active
    % pick best candidate BS based on DL SINR
    [bestSINR, bestBS] = max(L.TN.SINR_DL);
    servSINR = L.TN.SINR_DL(S.servBS);

    if bestBS ~= S.servBS && bestSINR > servSINR * P.HO.hyst
        S.HO_TN.ttt_s = S.HO_TN.ttt_s + P.dt;
        S.HO_TN.cand = bestBS;
        if S.HO_TN.ttt_s >= P.HO.TTT_TN_s
            % start HO execution
            S.HO_TN.active = true;
            S.HO_TN.timer_s = P.HO.exec_TN_s;
            S.HO_TN.ttt_s = 0;
        end
    else
        S.HO_TN.ttt_s = 0;
        S.HO_TN.cand = S.servBS;
    end
else
    % execution phase
    S.HO_TN.timer_s = S.HO_TN.timer_s - P.dt;
    if S.HO_TN.timer_s <= 0
        % complete HO
        S.servBS = S.HO_TN.cand;
        S.HO_TN.active = false;
        S.HO_TN.timer_s = 0;
        hoTN = true;
    end
end

% --- NTN HO ---
if ~S.HO_NTN.active
    [bestSINR, bestSAT] = max(L.NTN.SINR_DL);
    servSINR = L.NTN.SINR_DL(S.servSAT);

    if bestSAT ~= S.servSAT && bestSINR > servSINR * P.HO.hyst
        S.HO_NTN.ttt_s = S.HO_NTN.ttt_s + P.dt;
        S.HO_NTN.cand = bestSAT;
        if S.HO_NTN.ttt_s >= P.HO.TTT_NTN_s
            S.HO_NTN.active = true;
            S.HO_NTN.timer_s = P.HO.exec_NTN_s;
            S.HO_NTN.ttt_s = 0;
        end
    else
        S.HO_NTN.ttt_s = 0;
        S.HO_NTN.cand = S.servSAT;
    end
else
    S.HO_NTN.timer_s = S.HO_NTN.timer_s - P.dt;
    if S.HO_NTN.timer_s <= 0
        S.servSAT = S.HO_NTN.cand;
        S.HO_NTN.active = false;
        S.HO_NTN.timer_s = 0;
        hoNTN = true;
    end
end

end