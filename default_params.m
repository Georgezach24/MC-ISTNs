function P = default_params()

% Time
P.dt = 1e-3;          % 1 ms slot
P.simTime = 10;       % seconds
P.T = round(P.simTime / P.dt);

% --- Link models  ---
P.link.useSimpleRateModel = true;

% TN SNR trace parameters (dB)
P.TN.snr_mu = 12;
P.TN.snr_sigma_fast = 2.0;
P.TN.snr_sin_amp = 4;
P.TN.snr_sin_hz  = 0.7;

% NTN SNR trace parameters (dB)
P.NTN.snr_mu = 6;
P.NTN.snr_sigma_fast = 0.7;
P.NTN.snr_sin_amp = 2;
P.NTN.snr_sin_hz  = 0.1;

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

end