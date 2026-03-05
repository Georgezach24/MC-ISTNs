function L = step_channel(P, S, k)

t = (k-1)*P.dt;

% SNR traces (dB)
snrTN = P.TN.snr_mu ...
    + P.TN.snr_sin_amp*sin(2*pi*P.TN.snr_sin_hz*t) ...
    + P.TN.snr_sigma_fast*randn();

snrNTN = P.NTN.snr_mu ...
    + P.NTN.snr_sin_amp*sin(2*pi*P.NTN.snr_sin_hz*t + 1.0) ...
    + P.NTN.snr_sigma_fast*randn();

snrTN  = min(max(snrTN,  P.link.snr_min), P.link.snr_max);
snrNTN = min(max(snrNTN, P.link.snr_min), P.link.snr_max);

L.SNR_TN_dB  = snrTN;
L.SNR_NTN_dB = snrNTN;

% Rate model:
% R = eta * BW * log2(1 + SNR_linear)
snrTN_lin  = 10^(snrTN/10);
snrNTN_lin = 10^(snrNTN/10);

L.R_TN  = P.link.eta * P.link.BW * log2(1 + snrTN_lin);   % bits/s
L.R_NTN = P.link.eta * P.link.BW * log2(1 + snrNTN_lin);  % bits/s

% Convert to per-slot capacities (bits per slot)
L.capTN_bits  = L.R_TN  * P.dt;
L.capNTN_bits = L.R_NTN * P.dt;

end