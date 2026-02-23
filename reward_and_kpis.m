function [r, kpi] = reward_and_kpis(P, stepKPI, actionChanged, hoEvent)
%REWARD_AND_KPIS Simple reward prototype

if stepKPI.URLLC_delivered == 1
    if stepKPI.URLLC_success == 1
        rU = 1.0;
    else
        rU = -10.0;
    end
else
    rU = 0.0;
end

rE = 1e-7 * stepKPI.eMBB_bitsDL;

penSwitch = -0.1 * double(actionChanged);
penHO = -1.0 * double(hoEvent);

r = rU + rE + penSwitch + penHO;
kpi = stepKPI;

end