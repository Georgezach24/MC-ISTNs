clear; clc;

scenarioList = {'A','C'};

for s = 1:numel(scenarioList)

    scenario = scenarioList{s};
    P = default_params(scenario);
    S = init_scenario(P);

    S.SATx = zeros(1, numel(P.SAT));

    lastAction = S.action;

    for k = 1:P.T
        S.t = (k-1)*P.dt;

        % Update geometry and traffic
        S = update_positions(S, P, k);
        S = traffic_step(S, P, k);

        % Compute all link metrics
        L = compute_links(S, P);

        % Handover update
        [S, hoTN, hoNTN] = handover_update(S, P, L);
        hoEvent = hoTN || hoNTN;

        % -----------------------
        % Baseline Policy
        % -----------------------
        action = S.action;

        sinrLow = (L.SINR_TN_DL < P.SINRmin) || ...
                  (L.SINR_NTN_DL < P.SINRmin);

        if sinrLow
            action.modeU = 3; % DUP
        else
            if L.SINR_TN_DL >= L.SINR_NTN_DL
                action.modeU = 0;
            else
                action.modeU = 1;
            end
        end

        if (L.SINR_TN_DL > P.SINRmin) && ...
           (L.SINR_NTN_DL > P.SINRmin)
            action.modeE = 2;
        else
            action.modeE = uint8(L.SINR_NTN_DL > L.SINR_TN_DL);
        end

        actionChanged = (action.modeU ~= lastAction.modeU) || ...
                        (action.modeE ~= lastAction.modeE);

        [S, stepKPI] = apply_action_and_serve(S, P, L, action);

        [~, ~] = reward_and_kpis(P, stepKPI, actionChanged, hoEvent);

        % Logging
        S.log.actionU(k) = uint8(action.modeU);
        S.log.actionE(k) = uint8(action.modeE);

        S.log.SINR_TN_DL(k) = L.SINR_TN_DL;
        S.log.SINR_NTN_DL(k) = L.SINR_NTN_DL;

        S.log.hoTN(k) = hoTN;
        S.log.hoNTN(k) = hoNTN;

        S.log.servBS(k)  = uint8(S.servBS);
        S.log.servSAT(k) = uint8(S.servSAT);
        S.log.ueX(k)     = S.UE.x;

        S.log.kpi_URLLC_succ(k) = stepKPI.URLLC_success;
        S.log.kpi_URLLC_delay(k) = stepKPI.URLLC_delay;

        lastAction = action;
        S.action = action;
    end

    outFile = sprintf('dataset_%s.mat', scenario);
    save(outFile, 'P', 'S');
    fprintf('Saved %s\n', outFile);

end