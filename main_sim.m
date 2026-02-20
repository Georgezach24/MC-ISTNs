clear; clc;

scenarioList = {'A','C'};

for s = 1:numel(scenarioList)
    scenario = scenarioList{s};
    P = default_params(scenario);
    S = init_scenario(P);

    % Satellite x holder
    S.SATx = zeros(1,2);

    % Simple baseline policy: URLLC uses DUP when either SINR low, else TN.
    % eMBB uses SPLIT when both links decent, else best.
    lastAction = S.action;

    for k = 1:P.T
        S.t = (k-1)*P.dt;

        S = update_positions(S, P, k);
        S = traffic_step(S, P, k);

        L = compute_links(S, P);

        % --- Baseline policy (replace later with RL) ---
        action = S.action;

        sinrLow = (L.SINR_TN_DL < P.SINRmin) || (L.SINR_NTN_DL < P.SINRmin);
        if sinrLow
            action.modeU = 3; % DUP
        else
            % choose best single link for URLLC
            if L.SINR_TN_DL >= L.SINR_NTN_DL
                action.modeU = 0;
            else
                action.modeU = 1;
            end
        end

        % eMBB
        if (L.SINR_TN_DL > P.SINRmin) && (L.SINR_NTN_DL > P.SINRmin)
            action.modeE = 2; % SPLIT/aggregate
        else
            action.modeE = (L.SINR_NTN_DL > L.SINR_TN_DL); % 0 TN, 1 NTN
        end

        actionChanged = (action.modeU ~= lastAction.modeU) || (action.modeE ~= lastAction.modeE);
        hoEvent = false; % will implement next step

        [S, stepKPI] = apply_action_and_serve(S, P, L, action);
        [r, ~] = reward_and_kpis(P, stepKPI, actionChanged, hoEvent);

        % --- Log ---
        S.log.actionU(k) = uint8(action.modeU);
        S.log.actionE(k) = uint8(action.modeE);

        S.log.SINR_TN_DL(k) = L.SINR_TN_DL;
        S.log.SINR_TN_UL(k) = L.SINR_TN_UL;
        S.log.SINR_NTN_DL(k) = L.SINR_NTN_DL;
        S.log.SINR_NTN_UL(k) = L.SINR_NTN_UL;

        S.log.rate_TN_DL(k) = L.R_TN_DL;
        S.log.rate_TN_UL(k) = L.R_TN_UL;
        S.log.rate_NTN_DL(k) = L.R_NTN_DL;
        S.log.rate_NTN_UL(k) = L.R_NTN_UL;

        S.log.kpi_URLLC_succ(k) = stepKPI.URLLC_success;
        S.log.kpi_URLLC_delay(k) = stepKPI.URLLC_delay;
        S.log.kpi_eMBB_bitsDL(k) = stepKPI.eMBB_bitsDL;
        S.log.kpi_eMBB_bitsUL(k) = stepKPI.eMBB_bitsUL;

        S.log.energy_J(k) = S.energy_J;

        % update state
        S.action = action;
        lastAction = action;
    end

    outFile = sprintf('dataset_%s.mat', scenario);
    save(outFile, 'P', 'S');

    fprintf('Saved %s\n', outFile);
end
