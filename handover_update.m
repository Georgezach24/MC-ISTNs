function S = handover_update(S, P, L)

for u = 1:P.Nue

    % -----------------------------------
    % 1) Update filtered best TN/NTN rates
    % -----------------------------------
    if S.UE(u).filtRateTN == 0
        S.UE(u).filtRateTN = L.bestTNrate(u);
    else
        S.UE(u).filtRateTN = P.HO.alpha * S.UE(u).filtRateTN + ...
            (1 - P.HO.alpha) * L.bestTNrate(u);
    end

    if S.UE(u).filtRateNTN == 0
        S.UE(u).filtRateNTN = L.bestNTNrate(u);
    else
        S.UE(u).filtRateNTN = P.HO.alpha * S.UE(u).filtRateNTN + ...
            (1 - P.HO.alpha) * L.bestNTNrate(u);
    end

    % -----------------------------------
    % 2) If interruption active, count down
    % -----------------------------------
    if S.UE(u).HO.active
        S.UE(u).HO.interruptTimer = S.UE(u).HO.interruptTimer - 1;

        if S.UE(u).HO.interruptTimer <= 0
            S.UE(u).HO.active = false;
            S.UE(u).HO.interruptTimer = 0;
        end

        continue;
    end

    % -----------------------------------
    % 3) Current serving metric
    % -----------------------------------
    if S.UE(u).servingType == 0
        currentMetric = S.UE(u).filtRateTN;
        otherMetric   = S.UE(u).filtRateNTN;
        candType = 1;
        candID   = L.bestNTNid(u);
    else
        currentMetric = S.UE(u).filtRateNTN;
        otherMetric   = S.UE(u).filtRateTN;
        candType = 0;
        candID   = L.bestTNid(u);
    end

    betterCondition = otherMetric > P.HO.marginRatio * currentMetric;

    % -----------------------------------
    % 4) TTT logic
    % -----------------------------------
    if betterCondition
        if S.UE(u).HO.candidateType == candType && S.UE(u).HO.candidateID == candID
            S.UE(u).HO.timer = S.UE(u).HO.timer + 1;
        else
            S.UE(u).HO.candidateType = candType;
            S.UE(u).HO.candidateID = candID;
            S.UE(u).HO.timer = 1;
        end
    else
        S.UE(u).HO.candidateType = -1;
        S.UE(u).HO.candidateID = -1;
        S.UE(u).HO.timer = 0;
    end

    % -----------------------------------
    % 5) Execute handover
    % -----------------------------------
    if S.UE(u).HO.timer >= P.HO.TTT
        S.UE(u).servingType = S.UE(u).HO.candidateType;
        S.UE(u).servingID   = S.UE(u).HO.candidateID;

        S.UE(u).HO.active = true;
        S.UE(u).HO.interruptTimer = P.HO.interrupt;
        S.UE(u).HO.count = S.UE(u).HO.count + 1;

        S.UE(u).HO.candidateType = -1;
        S.UE(u).HO.candidateID = -1;
        S.UE(u).HO.timer = 0;
    end
end

end