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
end