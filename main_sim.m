clear; clc; close all;

P = default_params();
S = init_state(P);

% Logs
Log.t = (0:P.T-1)*P.dt;
Log.R_TN  = zeros(P.T,1);
Log.R_NTN = zeros(P.T,1);
Log.act   = zeros(P.T,1); % 0=TN, 1=NTN, 2=SPLIT, 3=DUP 

% Example KPIs
Log.bitsDL = zeros(P.T,1);  % delivered bits per slot 
Log.qURLLC = zeros(P.T,1);
Log.qeMBB  = zeros(P.T,1);

for k = 1:P.T
    % 1) Channel / link rates this slot
    L = step_channel(P, S, k);

    % 2) Traffic arrivals
    S = traffic_step(P, S, k);

    % 3) Policy action (baseline: επιλογή καλύτερου rate)
    action = step_policy(P, S, L, k);

    % 4) Serve queues given action + capacities
    [S, stepOut] = step_serve(P, S, L, action);

    % 5) KPI/log
    Log.R_TN(k)  = L.R_TN;
    Log.R_NTN(k) = L.R_NTN;
    Log.act(k)   = action;
    Log.bitsDL(k)= stepOut.bitsDL;
    Log.qURLLC(k)= S.Q.urlLC_bits;
    Log.qeMBB(k) = S.Q.eMBB_bits;
end

plot_results(P, Log);
disp("Simulation Done succesfully!.");