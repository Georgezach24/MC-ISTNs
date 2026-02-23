function P = default_params(scenario)

P.scenario = scenario;
P.dt = 1e-3;            % 1 ms
P.simTime = 30;         % seconds 
P.T = round(P.simTime / P.dt);

% Service QoS
switch upper(scenario)
    case 'A'
        P.URLLC.deadline = 20e-3;     % 20 ms
        P.URLLC.targetRel = 0.9999;
        P.eMBB.softDelay = 200e-3;    % 200 ms
        P.PDCP.Treord = 5e-3;         % reordering timer 
    case 'C'
        P.URLLC.deadline = 10e-3;     % 10 ms
        P.URLLC.targetRel = 0.99999;
        P.eMBB.softDelay = 150e-3;    % 150 ms
        P.PDCP.Treord = 2e-3;
    otherwise
        error('Scenario must be A or C');
end

% --- Handover params (network-controlled) ---
P.HO.hyst_dB = 3;                 % hysteresis
P.HO.hyst = 10^(P.HO.hyst_dB/10); % linear

P.HO.TTT_TN_s  = 40e-3;           % 40 ms
P.HO.TTT_NTN_s = 80e-3;           % 80 ms

P.HO.exec_TN_s  = 10e-3;          % 10 ms interruption
P.HO.exec_NTN_s = 30e-3;          % 30 ms interruption


% Traffic
P.URLLC.pktBytes = 64;
P.URLLC.pktBits  = 8 * P.URLLC.pktBytes;
P.URLLC.period   = 10e-3;            % 1 pkt per 10 ms (100 pps)

P.eMBB.backlogged = true;
P.eMBB.pktBits = 12000;


% Topology (simple 1D line for start)
P.areaLen = 2000;                     % meters

% UE mobility
P.UE.v = 15;                          % m/s
P.UE.x0 = 200;                        % m

% BS positions (2 BS)
P.BS(1).x = 0;
P.BS(2).x = P.areaLen;

% Sat "projection" positions (2 sats) - simple moving points above ground
P.SAT(1).x0 = 300;
P.SAT(2).x0 = 1700;

for i = 1:2
    P.SAT(i).v = 200;        % m/s ground-track proxy
    P.SAT(i).alt = 500e3;    % 500 km
end

% RF / Channel
P.fc = 2e9;                           % 2 GHz (can change)
P.c  = 3e8;
P.N0 = 1e-20;                         % W/Hz (rough), used with W
P.W  = 10e6;                          % 10 MHz

% Tx powers (linear W)
P.TN.Pdl = dbm2watt(46);              % BS DL 46 dBm
P.TN.Pul = dbm2watt(23);              % UE UL 23 dBm

P.NTN.Pdl = dbm2watt(50);             % SAT DL 50 dBm
P.NTN.Pul = dbm2watt(23);             % UE UL 23 dBm

% Pathloss exponents / additional losses
P.TN.alpha = 3.5;
P.NTN.alpha = 2.0;
P.NTN.extraLoss_dB = 10;              % lumped atmospheric + misc

% Fading
P.TN.fading = 'rayleigh';
P.NTN.fading = 'nakagami';
P.NTN.nak_m = 2;

% SINR thresholds for outage
P.SINRmin_dB = -5;
P.SINRmin = 10^(P.SINRmin_dB/10);              % -5 dB threshold (tune later)

% Energy model (very simple)
P.energy.rxCircuit_W = 0.2;           % receive circuitry power
P.energy.switchCost_J = 1e-3;         % per action change (later)
P.energy.hoCost_J = 5e-3;             % per HO event (later)

end

function w = dbm2watt(dbm)
w = 10.^((dbm - 30)/10);
end
