function S = init_scenario(P)
%INIT_SCENARIO Initialize simulation state and preallocate URLLC storage

% UE
S.t = 0;
S.UE.x = P.UE.x0;
S.UE.v = P.UE.v;

% Serving nodes
S.servBS = 1;
S.servSAT = 1;

% Current action
S.action.modeU = 0; % 0 TN, 1 NTN, 2 SPLIT, 3 DUP
S.action.modeE = 0; % 0 TN, 1 NTN, 2 SPLIT

% -----------------------
% Handover state (TN)
% -----------------------
S.HO_TN.active = false;
S.HO_TN.timer_s = 0;
S.HO_TN.ttt_s = 0;
S.HO_TN.cand = S.servBS;
S.HO_TN.filtSINR = ones(1, numel(P.BS));

% -----------------------
% Handover state (NTN)
% -----------------------
S.HO_NTN.active = false;
S.HO_NTN.timer_s = 0;
S.HO_NTN.ttt_s = 0;
S.HO_NTN.cand = S.servSAT;
S.HO_NTN.filtSINR = ones(1, numel(P.SAT));

% -----------------------
% URLLC packet store (preallocated FIFO)
% -----------------------
S.URLLC.genCount  = 0;
S.URLLC.succCount = 0;
S.URLLC.failCount = 0;

S.URLLC.head = 1;  % index of head-of-line active pkt
S.URLLC.tail = 0;  % number of generated pkts

% Estimate max packets based on simulation time and generation rate
% URLLC: 1 pkt per 10ms => 100 pkt/s
genRate = 100; % pkt/s
simTime = P.T * P.dt;

S.URLLC.maxPkts = ceil(simTime * genRate) + 2000; % margin

S.URLLC.queue(S.URLLC.maxPkts,1) = struct( ...
    'genTime', 0, ...
    'delivered', false, ...
    'expired', false, ...
    'delay', NaN);

% Will be set by traffic_step each iteration
S.URLLC.pktBits = 64*8;
S.URLLC.nextGenTime = 0;

% -----------------------
% Logs (preallocate common ones here if you want)
% -----------------------
S.log = struct();

% Optional: preallocate frequently used logs to speed up
S.log.actionU = zeros(P.T,1,'uint8');
S.log.actionE = zeros(P.T,1,'uint8');

S.log.SINR_TN_DL  = zeros(P.T,1);
S.log.SINR_NTN_DL = zeros(P.T,1);

S.log.hoTN  = zeros(P.T,1);
S.log.hoNTN = zeros(P.T,1);

S.log.servBS  = zeros(P.T,1,'uint8');
S.log.servSAT = zeros(P.T,1,'uint8');
S.log.ueX     = zeros(P.T,1);

S.log.kpi_URLLC_succ  = zeros(P.T,1);
S.log.kpi_URLLC_delay = nan(P.T,1);

end