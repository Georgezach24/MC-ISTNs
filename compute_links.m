function L = compute_links(S, P)
%COMPUTE_LINKS Compute DL/UL SINR and rates for serving TN and NTN links

uex = S.UE.x;

% Serving BS and SAT x positions
bsx = P.BS(S.servBS).x;

satIdx = S.servSAT;
satx = S.SATx(satIdx);
alt  = P.SAT(satIdx).alt;   

% Distances 
d_tn  = max(abs(uex - bsx), 1.0);              % meters
d_ntn = sqrt((uex - satx).^2 + alt.^2);        % meters

% Pathloss (linear)
PL_tn  = (d_tn.^P.TN.alpha);
PL_ntn = (d_ntn.^P.NTN.alpha) * db2lin(P.NTN.extraLoss_dB);

% Fading power gains
g_tn  = fading_gain(P.TN.fading, 1, P);
g_ntn = fading_gain(P.NTN.fading, 1, P);

% Received powers (signal)
Pr_tn_DL  = P.TN.Pdl  * g_tn  / PL_tn;
Pr_ntn_DL = P.NTN.Pdl * g_ntn / PL_ntn;

Pr_tn_UL  = P.TN.Pul  * g_tn  / PL_tn;
Pr_ntn_UL = P.NTN.Pul * g_ntn / PL_ntn;

% Noise
N = P.N0 * P.W;

% Interference 
I_tn  = 0;
I_ntn = 0;

SINR_tn_DL  = Pr_tn_DL  / (I_tn  + N);
SINR_ntn_DL = Pr_ntn_DL / (I_ntn + N);
SINR_tn_UL  = Pr_tn_UL  / (I_tn  + N);
SINR_ntn_UL = Pr_ntn_UL / (I_ntn + N);

% Rates (Shannon)
R_tn_DL  = P.W * log2(1 + SINR_tn_DL);
R_ntn_DL = P.W * log2(1 + SINR_ntn_DL);
R_tn_UL  = P.W * log2(1 + SINR_tn_UL);
R_ntn_UL = P.W * log2(1 + SINR_ntn_UL);

% Prop delays
Dprop_tn  = d_tn  / P.c;
Dprop_ntn = d_ntn / P.c;

% Outage flags
out_tn  = (SINR_tn_DL  < P.SINRmin);
out_ntn = (SINR_ntn_DL < P.SINRmin);

% Output struct
L = struct();
L.SINR_TN_DL  = SINR_tn_DL;
L.SINR_NTN_DL = SINR_ntn_DL;
L.SINR_TN_UL  = SINR_tn_UL;
L.SINR_NTN_UL = SINR_ntn_UL;

L.R_TN_DL  = R_tn_DL;
L.R_NTN_DL = R_ntn_DL;
L.R_TN_UL  = R_tn_UL;
L.R_NTN_UL = R_ntn_UL;

L.Dprop_TN  = Dprop_tn;
L.Dprop_NTN = Dprop_ntn;

L.out_TN  = out_tn;
L.out_NTN = out_ntn;

end

function g = fading_gain(type, n, P)
switch lower(type)
    case 'rayleigh'
        h = (randn(n,1) + 1j*randn(n,1))/sqrt(2);
        g = abs(h).^2;

    case 'nakagami'
        m = P.NTN.nak_m;   % assume integer m
        g = 0;
        for i = 1:m
            g = g - log(rand);   % sum of exponentials
        end
        g = g / m;

    otherwise
        g = ones(n,1);
end

g = g(1);
end

function x = db2lin(dB)
x = 10.^(dB/10);
end
