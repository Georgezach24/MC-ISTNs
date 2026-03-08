function [S, L] = step_channel(P, S, k)

t = (k-1) * P.dt;

L.R_TN_all = zeros(P.Nue, P.Nbs);
L.R_NTN_all = zeros(P.Nue, P.Nsat);

L.bestTNrate = zeros(P.Nue,1);
L.bestTNid   = zeros(P.Nue,1);

L.bestNTNrate = zeros(P.Nue,1);
L.bestNTNid   = zeros(P.Nue,1);

for u = 1:P.Nue
    % -------------------------
    % Update user position
    % -------------------------
    S.UE(u).x = S.UE(u).x + S.UE(u).v * P.dt;

    % keep inside area with simple wrap-around
    if S.UE(u).x > P.areaLen
        S.UE(u).x = S.UE(u).x - P.areaLen;
    end

    x = S.UE(u).x;

    % -------------------------
    % TN rates to all BSs
    % -------------------------
    for b = 1:P.Nbs
        d = abs(x - P.BS.pos(b)) + 1; % avoid zero distance

        snr_dB = P.TN.snr0 ...
            - P.TN.pathlossCoeff * d ...
            + 2*sin(2*pi*0.2*t + 0.4*u + 0.3*b) ...
            + P.TN.fastSigma * randn();

        snr_dB = min(max(snr_dB, P.link.snrMin), P.link.snrMax);
        snr_lin = 10^(snr_dB/10);

        L.R_TN_all(u,b) = P.link.eta * P.link.BW * log2(1 + snr_lin);
    end

    % -------------------------
    % NTN rates to all satellites
    % -------------------------
    for s = 1:P.Nsat
        d = abs(x - P.SAT.pos(s)) + 1;

        snr_dB = P.NTN.snr0 ...
            - P.NTN.pathlossCoeff * d ...
            + 1.5*sin(2*pi*0.05*t + 0.5*u + 0.2*s) ...
            + P.NTN.fastSigma * randn();

        snr_dB = min(max(snr_dB, P.link.snrMin), P.link.snrMax);
        snr_lin = 10^(snr_dB/10);

        L.R_NTN_all(u,s) = P.link.eta * P.link.BW * log2(1 + snr_lin);
    end

    % -------------------------
    % Best TN and NTN candidate
    % -------------------------
    [L.bestTNrate(u), idxTN] = max(L.R_TN_all(u,:));
    [L.bestNTNrate(u), idxNTN] = max(L.R_NTN_all(u,:));

    L.bestTNid(u) = idxTN;
    L.bestNTNid(u) = idxNTN;
end

end