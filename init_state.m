function S = init_state(P)

for u = 1:P.Nue
    % User position
    S.UE(u).x = rand * P.areaLen;

    % User speed
    S.UE(u).v = max(0.5, P.UE.vMean + P.UE.vStd*randn());

    % User profile
    if rand < P.profile.urlProb
        S.UE(u).profile = 1;   % URLLC-oriented
    else
        S.UE(u).profile = 0;   % eMBB-oriented
    end

    % Queues
    S.UE(u).Q.urlLC_bits = 0;
    S.UE(u).Q.eMBB_bits = 0;

    % Serving state
    S.UE(u).servingType = 0;   % 0 = TN, 1 = NTN
    S.UE(u).servingID = 1;

    % Handover state
    S.UE(u).HO.candidateType = -1;
    S.UE(u).HO.candidateID = -1;
    S.UE(u).HO.timer = 0;
    S.UE(u).HO.active = false;
    S.UE(u).HO.interruptTimer = 0;
    S.UE(u).HO.count = 0;

    % Time since last HO
    S.UE(u).dwellTime = 0;
end
end



