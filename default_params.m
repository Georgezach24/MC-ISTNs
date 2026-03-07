function P = default_params()

% Time
P.dt = 1e-3;          % 1 ms slot
P.simTime = 10;       % seconds
P.T = round(P.simTime / P.dt);

% --- Link models  ---
P.link.useSimpleRateModel = true;

% TN SNR trace parameters (dB)
P.TN.snr_mu = 10;
P.TN.snr_sigma_fast = 1.5;
P.TN.snr_sin_amp = 3;
P.TN.snr_sin_hz  = 0.2;

% NTN SNR trace parameters (dB)
P.NTN.snr_mu = 9.5;
P.NTN.snr_sigma_fast = 1.0;
P.NTN.snr_sin_amp = 3;
P.NTN.snr_sin_hz  = 0.2;

% Clamp
P.link.snr_min = -10;
P.link.snr_max = 30;

% Bandwidth for simple rate model
P.link.BW = 20e6; % 20 MHz

% “Efficiency loss” factor 
P.link.eta = 0.6; % 0..1

% Traffic (bits per slot arrivals)
P.traffic.urlLC_mean_bits_per_s = 2e6; % 2 Mbps offered load 
P.traffic.eMBB_mean_bits_per_s  = 20e6;

% Queue limits 
P.Q.maxBits = 200e6;

% Handover parameters
P.HO.hyst_dB_equiv = 2;   % margin για να θεωρηθεί το άλλο link καλύτερο
P.HO.TTT = 10;            % time-to-trigger σε slots (50 ms αν dt=1ms)
P.HO.interrupt = 10;      % interruption duration σε slots (20 ms)
P.HO.alpha = 0.7;         % filtering coefficient

end