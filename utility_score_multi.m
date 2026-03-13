function [score, info] = utility_score_multi(P, S, L, u, targetType, targetID, loadTN, loadNTN)

% targetType: 0 = TN, 1 = NTN

currentType = S.UE(u).servingType;
currentID   = S.UE(u).servingID;
isStay = (targetType == currentType) && (targetID == currentID);

% --------------------------------------
% 1) Raw rate, delay and energy by domain
% --------------------------------------
if targetType == 0
    rawRate = L.R_TN_all(u, targetID);
    baseDelay_ms = P.domain.TN.baseDelay_ms;
    energyCost = P.domain.TN.energyCost;
    currentLoad = loadTN(targetID);
else
    rawRate = L.R_NTN_all(u, targetID);
    baseDelay_ms = P.domain.NTN.baseDelay_ms;
    energyCost = P.domain.NTN.energyCost;
    currentLoad = loadNTN(targetID);
end

% --------------------------------------
% 2) Expected load if user goes there
% --------------------------------------
if isStay
    expectedLoad = max(currentLoad,1);
else
    expectedLoad = currentLoad + 1;
end

% --------------------------------------
% 3) Effective rate after sharing
% --------------------------------------
effectiveRate = rawRate / expectedLoad;

% --------------------------------------
% 4) Queue pressure and estimated queue delay
% --------------------------------------
if S.UE(u).profile == 1
    % URLLC-oriented: URLLC queue matters more
    queueBits = S.UE(u).Q.urlLC_bits + 0.25*S.UE(u).Q.eMBB_bits;
    W = P.U.URLLC;
else
    % eMBB-oriented: eMBB queue matters more
    queueBits = 0.25*S.UE(u).Q.urlLC_bits + S.UE(u).Q.eMBB_bits;
    W = P.U.eMBB;
end

queueDelay_ms = 1000 * queueBits / max(effectiveRate,1);
totalDelay_ms = baseDelay_ms + queueDelay_ms;

% --------------------------------------
% 5) Normalized metrics
% --------------------------------------
rateNorm   = min(effectiveRate / P.util.rateRef, 2.0);
delayNorm  = min(totalDelay_ms / P.util.delayRef_ms, 2.0);
queueNorm  = min(queueBits / P.util.queueRef_bits, 2.0);
loadNorm   = min(expectedLoad / P.util.loadRef, 2.0);
energyNorm = min(energyCost / P.util.energyRef, 2.0);

% --------------------------------------
% 6) HO penalty
% More expensive if user has just switched recently
% --------------------------------------
if isStay
    hoPenalty = 0;
else
    dwellPenalty = max(0, (P.HO.minDwell - S.UE(u).dwellTime) / P.HO.minDwell);
    hoPenalty = 1 + dwellPenalty;
end

% --------------------------------------
% 7) Final utility score
% --------------------------------------
score = ...
    + W.wRate   * rateNorm ...
    - W.wDelay  * delayNorm ...
    - W.wLoad   * loadNorm ...
    - W.wEnergy * energyNorm ...
    - W.wQueue  * queueNorm ...
    - W.wHO     * hoPenalty;

% optional debug info
info.effectiveRate = effectiveRate;
info.totalDelay_ms = totalDelay_ms;
info.expectedLoad  = expectedLoad;
info.queueBits     = queueBits;
info.hoPenalty     = hoPenalty;

end