function S = init_scenario(P)
%INIT_SCENARIO Initialize sim state

S.t = 0;

% UE
S.UE.x = P.UE.x0;
S.UE.v = P.UE.v;

% Queues as "bits waiting"
S.qDL_URLLC = 0;
S.qUL_URLLC = 0;
S.qDL_eMBB  = 0;
S.qUL_eMBB  = 0;

% Packet bookkeeping for URLLC delay stats
S.URLLC.nextGenTime = 0;
S.URLLC.inflight = struct('genTime', {}, 'bits', {}, 'delivered', {}, 'delay', {});
S.URLLC.rxCount = 0;
S.URLLC.succCount = 0;
S.URLLC.lateCount = 0;

% Current serving nodes (start with BS1 + SAT1)
S.servBS = 1;
S.servSAT = 1;

% Current action (start conservative)
S.action.modeU = 0; % 0=TN,1=NTN,2=SPLIT,3=DUP
S.action.modeE = 0; % 0=TN,1=NTN,2=SPLIT

% Logs (prealloc basic)
S.log.actionU = zeros(P.T,1,'uint8');
S.log.actionE = zeros(P.T,1,'uint8');

S.log.SINR_TN_DL = zeros(P.T,1);
S.log.SINR_TN_UL = zeros(P.T,1);
S.log.SINR_NTN_DL = zeros(P.T,1);
S.log.SINR_NTN_UL = zeros(P.T,1);

S.log.rate_TN_DL = zeros(P.T,1);
S.log.rate_TN_UL = zeros(P.T,1);
S.log.rate_NTN_DL = zeros(P.T,1);
S.log.rate_NTN_UL = zeros(P.T,1);

S.log.kpi_URLLC_succ = zeros(P.T,1);
S.log.kpi_URLLC_delay = zeros(P.T,1);
S.log.kpi_eMBB_bitsDL = zeros(P.T,1);
S.log.kpi_eMBB_bitsUL = zeros(P.T,1);

S.energy_J = 0;
S.log.energy_J = zeros(P.T,1);

end
