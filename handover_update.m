function S = handover_update(S, P, L)

% 1) Update filtered rates
if S.filtRateTN == 0
    S.filtRateTN = L.R_TN;
else
    S.filtRateTN = P.HO.alpha * S.filtRateTN + (1 - P.HO.alpha) * L.R_TN;
end

if S.filtRateNTN == 0
    S.filtRateNTN = L.R_NTN;
else
    S.filtRateNTN = P.HO.alpha * S.filtRateNTN + (1 - P.HO.alpha) * L.R_NTN;
end

% 2) If HO interruption is active, count it down
if S.HO.active
    S.HO.interruptTimer = S.HO.interruptTimer - 1;

    if S.HO.interruptTimer <= 0
        S.HO.active = false;
        S.HO.interruptTimer = 0;
    end

    return;
end

% 3) Compare serving and other link
if S.servingLink == 0
    servingMetric = S.filtRateTN;
    otherMetric   = S.filtRateNTN;
    otherLink     = 1;
else
    servingMetric = S.filtRateNTN;
    otherMetric   = S.filtRateTN;
    otherLink     = 0;
end

% Hysteresis rule:
% το άλλο link πρέπει να είναι τουλάχιστον 10% καλύτερο
betterCondition = otherMetric > 1.02 * servingMetric;

% 4) TTT logic
if betterCondition
    if S.HO.candidate == otherLink
        S.HO.timer = S.HO.timer + 1;
    else
        S.HO.candidate = otherLink;
        S.HO.timer = 1;
    end
else
    S.HO.candidate = -1;
    S.HO.timer = 0;
end

% 5) Execute HO if TTT reached
if S.HO.timer >= P.HO.TTT
    S.servingLink = S.HO.candidate;

    S.HO.active = true;
    S.HO.interruptTimer = P.HO.interrupt;

    S.HO.candidate = -1;
    S.HO.timer = 0;
end

end