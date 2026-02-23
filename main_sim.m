clear; clc;

scenarioList = {'A','C'};

for s = 1:numel(scenarioList)
    scenario = scenarioList{s};
    P = default_params(scenario);
    S = init_scenario(P);

    % Satellite x holder (must exist before compute_links)
    S.SATx = zeros(1, numel(P.SAT));

    % Baseline policy state
    lastAction = S.action;

    for k = 1:P.T
        S.t = (k-1)*P.dt;

        % Update geometry + traffic
        S = update_positions(S, P, k);
        S = traffic_step(S, P, k);

        % Compute candidate link metrics (ALL BSs and Sats)
        L = compute_links(S, P);

        % ---------------------------
        % Network-controlled handover
        % ---------------------------
        [S, hoTN, hoNTN] = handover_update(S, P, L);
        hoEvent = hoTN || hoNTN;

        % ---------------------------
        % Baseline policy (replace later with RL)
        % ---------------------------
        action = S.action;

        % URLLC: if either serving link is below SINRmin -> DUP, else best single serving link
        sinrLow = (L.SINR_TN_DL < P.SINRmin) || (L.SINR_NTN_DL < P.SINRmin);
        if sinrLow
            action.modeU = 3; % DUP
        else
            % choose best single (SERVING) link for URLLC
            if L.SINR_TN_DL >= L.SINR_NTN_DL
                action.modeU = 0; % TN
            else
                action.modeU = 1; % NTN
            end
        end

        % eMBB: if both serving links are decent -> SPLIT/aggregate else choose better serving link
        if (L.SINR_TN_DL > P.SINRmin) && (L.SINR_NTN_DL > P.SINRmin)
            action.modeE = 2; % SPLIT/aggregate
        else
            action.modeE = uint8(L.SINR_NTN_DL > L.SINR_TN_DL); % 0 TN, 1 NTN
        end

        actionChanged = (action.modeU ~= lastAction.modeU) || (action.modeE ~= lastAction.modeE);

        % Apply action and serve queues (HO-aware interruption handled inside apply_action_and_serve)
        [S, stepKPI] = apply_action_and_serve(S, P, L, action);

        % Reward (prototype)
        [r, ~] = reward_and_kpis(P, stepKPI, actionChanged, hoEvent); %#ok<NASGU>

        % ---------------------------
        % Logging
        % ---------------------------
        S.log.actionU(k) = uint8(action.modeU);
        S.log.actionE(k) = uint8(action.modeE);

        % serving-link SINR/rates (for easy plotting)
        S.log.SINR_TN_DL(k) = L.SINR_TN_DL;
        S.log.SINR_TN_UL(k) = L.SINR_TN_UL;
        S.log.SINR_NTN_DL(k) = L.SINR_NTN_DL;
        S.log.SINR_NTN_UL(k) = L.SINR_NTN_UL;

        S.log.rate_TN_DL(k) = L.R_TN_DL;
        S.log.rate_TN_UL(k) = L.R_TN_UL;
        S.log.rate_NTN_DL(k) = L.R_NTN_DL;
        S.log.rate_NTN_UL(k) = L.R_NTN_UL;

        % KPIs
        S.log.kpi_URLLC_succ(k) = stepKPI.URLLC_success;
        S.log.kpi_URLLC_delay(k) = stepKPI.URLLC_delay;
        S.log.kpi_eMBB_bitsDL(k) = stepKPI.eMBB_bitsDL;
        S.log.kpi_eMBB_bitsUL(k) = stepKPI.eMBB_bitsUL;

        % HO events
        S.log.hoTN(k) = hoTN;
        S.log.hoNTN(k) = hoNTN;

        % Energy
        S.log.energy_J(k) = S.energy_J;

        % Update state
        S.action = action;
        lastAction = action;
    end

    outFile = sprintf('dataset_%s.mat', scenario);
    save(outFile, 'P', 'S');

    fprintf('Saved %s\n', outFile);
end