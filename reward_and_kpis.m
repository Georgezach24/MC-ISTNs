function [r, kpi] = reward_and_kpis(P, stepKPI, actionChanged, hoEvent)
%REWARD_AND_KPIS Compute reward and KPI summary for RL/dataset

% URLLC: big penalty if fail
if stepKPI.URLLC_delivered == 1
    if stepKPI.URLLC_success == 1
        rU = +1.0;
    else
        rU = -10.0; % late/lost
    end
else
    rU = 0; % no packet delivered this slot
end

% eMBB: normalized throughput reward
rE = 1e-7 * stepKPI.eMBB_bitsDL; % scale to ~0..1

% penalties
penSwitch = -0.1 * double(actionChanged);
penHO = -1.0 * double(hoEvent);  % placeholder for next step
% energy penalty handled via logging or can be included if you want

r = rU + rE + penSwitch + penHO;

kpi = stepKPI;
end
