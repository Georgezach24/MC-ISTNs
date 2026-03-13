function S = handover_update(S, P, L)

% -----------------------------------
% 0) Snapshot of current node loads
% Exclude users currently under interruption
% -----------------------------------
loadTN  = zeros(P.Nbs,1);
loadNTN = zeros(P.Nsat,1);

for uu = 1:P.Nue
    if S.UE(uu).HO.active
        continue;
    end

    if S.UE(uu).servingType == 0
        loadTN(S.UE(uu).servingID) = loadTN(S.UE(uu).servingID) + 1;
    else
        loadNTN(S.UE(uu).servingID) = loadNTN(S.UE(uu).servingID) + 1;
    end
end

% -----------------------------------
% 1) User-by-user HO update
% -----------------------------------
for u = 1:P.Nue

    % If interruption active, count down
    if S.UE(u).HO.active
        S.UE(u).HO.interruptTimer = S.UE(u).HO.interruptTimer - 1;

        if S.UE(u).HO.interruptTimer <= 0
            S.UE(u).HO.active = false;
            S.UE(u).HO.interruptTimer = 0;
        end

        continue;
    end

    % Increase dwell time while user stays on current serving link
    S.UE(u).dwellTime = S.UE(u).dwellTime + 1;

    currentType = S.UE(u).servingType;
    currentID   = S.UE(u).servingID;

    % -----------------------------------
    % Candidate options
    % - stay on current
    % - best TN
    % - best NTN
    % -----------------------------------
    [scoreStay, ~] = utility_score_multi(P, S, L, u, currentType, currentID, loadTN, loadNTN);
    [scoreTN,   ~] = utility_score_multi(P, S, L, u, 0, L.bestTNid(u), loadTN, loadNTN);
    [scoreNTN,  ~] = utility_score_multi(P, S, L, u, 1, L.bestNTNid(u), loadTN, loadNTN);

    scores = [scoreStay, scoreTN, scoreNTN];
    [bestScore, idx] = max(scores);

    % Decode best target
    if idx == 1
        candType = currentType;
        candID   = currentID;
    elseif idx == 2
        candType = 0;
        candID   = L.bestTNid(u);
    else
        candType = 1;
        candID   = L.bestNTNid(u);
    end

    % Improvement over current
    improvement = bestScore - scoreStay;
    isDifferent = ~(candType == currentType && candID == currentID);

    % -----------------------------------
    % TTT logic based on utility gain
    % -----------------------------------
    if isDifferent && (improvement > P.HO.utilityMargin)

        if S.UE(u).HO.candidateType == candType && S.UE(u).HO.candidateID == candID
            S.UE(u).HO.timer = S.UE(u).HO.timer + 1;
        else
            S.UE(u).HO.candidateType = candType;
            S.UE(u).HO.candidateID   = candID;
            S.UE(u).HO.timer = 1;
        end

    else
        S.UE(u).HO.candidateType = -1;
        S.UE(u).HO.candidateID   = -1;
        S.UE(u).HO.timer = 0;
    end

    % -----------------------------------
    % Execute HO
    % -----------------------------------
    if S.UE(u).HO.timer >= P.HO.TTT
        S.UE(u).servingType = S.UE(u).HO.candidateType;
        S.UE(u).servingID   = S.UE(u).HO.candidateID;

        S.UE(u).HO.active = true;
        S.UE(u).HO.interruptTimer = P.HO.interrupt;
        S.UE(u).HO.count = S.UE(u).HO.count + 1;

        S.UE(u).HO.candidateType = -1;
        S.UE(u).HO.candidateID   = -1;
        S.UE(u).HO.timer = 0;

        % reset dwell time after HO
        S.UE(u).dwellTime = 0;
    end
end

end