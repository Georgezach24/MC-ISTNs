function P = default_params()

% ---------------------------------
% Time
% ---------------------------------
P.dt = 1e-3;              % 1 ms
P.simTime = 10;           % seconds
P.T = round(P.simTime / P.dt);

% ---------------------------------
% Topology
% ---------------------------------
P.Nue = 12;                % number of users
P.Nbs = 2;                % number of TN base stations
P.Nsat = 2;               % number of NTN nodes
P.areaLen = 2000;         % 1D line [m]

% Fixed node positions (1D for simplicity)
P.BS.pos = [500 1500];
P.SAT.pos = [300 1700];

% ---------------------------------
% Mobility
% ---------------------------------
P.UE.vMean = 5;           % m/s
P.UE.vStd = 1;            % m/s randomization

% ---------------------------------
% Link model
% ---------------------------------
P.link.BW = 20e6;         % Hz
P.link.eta = 0.6;

P.link.snrMin = -10;
P.link.snrMax = 30;

% TN 
P.TN.snr0 = 22;           
P.TN.pathlossCoeff = 0.012;
P.TN.fastSigma = 1.2;

% NTN 
P.NTN.snr0 = 20;
P.NTN.pathlossCoeff = 0.007;
P.NTN.fastSigma = 0.7;

% ---------------------------------
% Traffic
% ---------------------------------
P.traffic.urlLC_mean_bits_per_s = 4e6;
P.traffic.eMBB_mean_bits_per_s  = 15e6;

P.Q.maxBits = 200e6;

% ---------------------------------
% Handover
% ---------------------------------
P.HO.alpha = 0.7;
P.HO.TTT = 12;                % slots
P.HO.interrupt = 8;          % slots
P.HO.marginRatio = 1.02;     

% ---------------------------------
% User profiles
% 1 = URLLC-oriented
% 0 = eMBB-oriented
% ---------------------------------
P.profile.urlProb = 0.4;

% ---------------------------------
% Traffic per profile
% ---------------------------------
P.traffic.URLLC_user.url_bits_per_s = 4e6;
P.traffic.URLLC_user.embb_bits_per_s = 8e6;

P.traffic.eMBB_user.url_bits_per_s = 1e6;
P.traffic.eMBB_user.embb_bits_per_s = 22e6;

% ---------------------------------
% Domain delay / energy proxy
% ---------------------------------
P.domain.TN.baseDelay_ms = 5;
P.domain.NTN.baseDelay_ms = 30;

P.domain.TN.energyCost = 1.0;
P.domain.NTN.energyCost = 1.25;

% ---------------------------------
% Utility normalization
% ---------------------------------
P.util.rateRef = 80e6;         % reference effective rate
P.util.delayRef_ms = 50;       % reference delay
P.util.queueRef_bits = 3e6;    % reference queue
P.util.loadRef = 5;            % reference number of users on node
P.util.energyRef = 1.0;        % reference energy cost

% ---------------------------------
% HO decision parameters
% ---------------------------------
P.HO.utilityMargin = 0.05;     % minimum score gain to consider HO
P.HO.minDwell = 30;            % anti-ping-pong dwell time (slots)

% ---------------------------------
% Utility weights
% ---------------------------------
P.U.URLLC.wRate   = 0.25;
P.U.URLLC.wDelay  = 0.25;
P.U.URLLC.wLoad   = 0.10;
P.U.URLLC.wEnergy = 0.05;
P.U.URLLC.wQueue  = 0.20;
P.U.URLLC.wHO     = 0.15;

P.U.eMBB.wRate   = 0.40;
P.U.eMBB.wDelay  = 0.05;
P.U.eMBB.wLoad   = 0.20;
P.U.eMBB.wEnergy = 0.10;
P.U.eMBB.wQueue  = 0.15;
P.U.eMBB.wHO     = 0.10;
end