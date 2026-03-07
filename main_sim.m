clear; clc; close all;

P = default_params();
S = init_state(P);

Log.t = (0:P.T-1)*P.dt;
Log.R_TN  = zeros(P.T,1);
Log.R_NTN = zeros(P.T,1);
Log.servingLink = zeros(P.T,1);
Log.bitsDL = zeros(P.T,1);
Log.qURLLC = zeros(P.T,1);
Log.qeMBB  = zeros(P.T,1);
Log.HOactive = zeros(P.T,1);
Log.HOcandidate = zeros(P.T,1);

for k = 1:P.T
    % 1) Channel
    L = step_channel(P, S, k);

    % 2) Handover update
    S = handover_update(S, P, L);

    % 3) Traffic arrivals
    S = traffic_step(P, S, k);

    % 4) Serve using current serving link
    [S, stepOut] = step_serve_handover(P, S, L);

    % 5) Logging
    Log.R_TN(k) = L.R_TN;
    Log.R_NTN(k) = L.R_NTN;
    Log.servingLink(k) = S.servingLink;
    Log.bitsDL(k) = stepOut.bitsDL;
    Log.qURLLC(k) = S.Q.urlLC_bits;
    Log.qeMBB(k) = S.Q.eMBB_bits;
    Log.HOactive(k) = S.HO.active;
    Log.HOcandidate(k) = S.HO.candidate;
end

plot_results(P, Log);
disp("Sim Done.");