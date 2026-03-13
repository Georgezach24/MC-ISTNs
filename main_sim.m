clear; clc; close all;

P = default_params();
S = init_state(P);

% --------------------------
% Logs
% --------------------------
Log.totalBitsDL = zeros(P.T,1);
Log.sumQ_URLLC  = zeros(P.T,1);
Log.sumQ_eMBB   = zeros(P.T,1);
Log.meanBestTN  = zeros(P.T,1);
Log.meanBestNTN = zeros(P.T,1);
Log.numHOactive = zeros(P.T,1);
Log.meanServingType = zeros(P.T,1); % 0 closer to TN, 1 closer to NTN

Log.servingType = zeros(P.T, P.Nue);
Log.servingID   = zeros(P.T, P.Nue);

Log.loadTN = zeros(P.T, P.Nbs);
Log.loadNTN = zeros(P.T, P.Nsat);

Log.profile     = zeros(P.Nue,1);      % 1=URLLC-oriented, 0=eMBB-oriented
Log.bitsDL_user = zeros(P.T, P.Nue);   % served bits per slot per user
Log.qURLLC_user = zeros(P.T, P.Nue);   % URLLC queue per user
Log.qeMBB_user  = zeros(P.T, P.Nue);   % eMBB queue per user
Log.HOcount_user_final = zeros(P.Nue,1);

for u = 1:P.Nue
    Log.profile(u) = S.UE(u).profile;
end

for k = 1:P.T
    % 1) Update channels and user positions
    [S, L] = step_channel(P, S, k);

    % 2) Handover update
    S = handover_update(S, P, L);

    % 3) Traffic arrivals
    S = traffic_step(P, S);

    % 4) Serve traffic
    [S, stepOut] = step_serve(P, S, L);
    for u = 1:P.Nue
        Log.bitsDL_user(k,u) = stepOut.bitsDL_user(u);
        Log.qURLLC_user(k,u) = S.UE(u).Q.urlLC_bits;
        Log.qeMBB_user(k,u)  = S.UE(u).Q.eMBB_bits;
    end

    % 5) Logs
    Log.totalBitsDL(k) = stepOut.totalBitsDL;
    Log.meanBestTN(k)  = mean(L.bestTNrate);
    Log.meanBestNTN(k) = mean(L.bestNTNrate);
    Log.loadTN(k,:) = stepOut.loadTN;
    Log.loadNTN(k,:) = stepOut.loadNTN;

    qU = 0;
    qE = 0;
    numActive = 0;
    meanServ = 0;

    for u = 1:P.Nue
        qU = qU + S.UE(u).Q.urlLC_bits;
        qE = qE + S.UE(u).Q.eMBB_bits;

        if S.UE(u).HO.active
            numActive = numActive + 1;
        end

        meanServ = meanServ + S.UE(u).servingType;

        Log.servingType(k,u) = S.UE(u).servingType;
        Log.servingID(k,u)   = S.UE(u).servingID;
    end

    Log.sumQ_URLLC(k) = qU;
    Log.sumQ_eMBB(k)  = qE;
    Log.numHOactive(k)= numActive;
    Log.meanServingType(k) = meanServ / P.Nue;
end

% Total HO count
totalHOs = 0;
for u = 1:P.Nue
    totalHOs = totalHOs + S.UE(u).HO.count;
end
for u = 1:P.Nue
    Log.HOcount_user_final(u) = S.UE(u).HO.count;
end

summarize_profile_kpis(P, Log);

fprintf('Total Handovers = %d\n', totalHOs);

plot_results(P, Log);