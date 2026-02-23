function L = compute_links(S, P)
%COMPUTE_LINKS Compute DL/UL SINR and rates for ALL candidate TN BSs and NTN sats.
% Returns arrays for candidates and convenience fields for serving links.

uex = S.UE.x;

NB = numel(P.BS);
NS = numel(P.SAT);

N = P.N0 * P.W;   % noise power (W)
I_tn = 0;
I_ntn = 0;

% --- TN candidates ---
SINR_TN_DL = zeros(1,NB);
SINR_TN_UL = zeros(1,NB);
R_TN_DL    = zeros(1,NB);
R_TN_UL    = zeros(1,NB);
Dprop_TN   = zeros(1,NB);

for b = 1:NB
    bsx = P.BS(b).x;

    d_tn = max(abs(uex - bsx), 1.0);
    PL_tn = d_tn.^P.TN.alpha;

    g_tn = fading_gain(P.TN.fading, P, 'TN');

    Pr_tn_DL = P.TN.Pdl * g_tn / PL_tn;
    Pr_tn_UL = P.TN.Pul * g_tn / PL_tn;

    SINR_TN_DL(b) = Pr_tn_DL / (I_tn + N);
    SINR_TN_UL(b) = Pr_tn_UL / (I_tn + N);

    R_TN_DL(b) = P.W * log2(1 + SINR_TN_DL(b));
    R_TN_UL(b) = P.W * log2(1 + SINR_TN_UL(b));

    Dprop_TN(b) = d_tn / P.c;
end

% --- NTN candidates ---
SINR_NTN_DL = zeros(1,NS);
SINR_NTN_UL = zeros(1,NS);
R_NTN_DL    = zeros(1,NS);
R_NTN_UL    = zeros(1,NS);
Dprop_NTN   = zeros(1,NS);

for s = 1:NS
    satx = S.SATx(s);
    alt  = P.SAT(s).alt;

    d_ntn = sqrt((uex - satx).^2 + alt.^2);
    PL_ntn = (d_ntn.^P.NTN.alpha) * db2lin_local(P.NTN.extraLoss_dB);

    g_ntn = fading_gain(P.NTN.fading, P, 'NTN');

    Pr_ntn_DL = P.NTN.Pdl * g_ntn / PL_ntn;
    Pr_ntn_UL = P.NTN.Pul * g_ntn / PL_ntn;

    SINR_NTN_DL(s) = Pr_ntn_DL / (I_ntn + N);
    SINR_NTN_UL(s) = Pr_ntn_UL / (I_ntn + N);

    R_NTN_DL(s) = P.W * log2(1 + SINR_NTN_DL(s));
    R_NTN_UL(s) = P.W * log2(1 + SINR_NTN_UL(s));

    Dprop_NTN(s) = d_ntn / P.c;
end

% Package candidate arrays
L = struct();
L.TN = struct('SINR_DL', SINR_TN_DL, 'SINR_UL', SINR_TN_UL, ...
              'R_DL', R_TN_DL, 'R_UL', R_TN_UL, 'Dprop', Dprop_TN);

L.NTN = struct('SINR_DL', SINR_NTN_DL, 'SINR_UL', SINR_NTN_UL, ...
               'R_DL', R_NTN_DL, 'R_UL', R_NTN_UL, 'Dprop', Dprop_NTN);

% Convenience: serving-link values
L.SINR_TN_DL  = SINR_TN_DL(S.servBS);
L.SINR_TN_UL  = SINR_TN_UL(S.servBS);
L.R_TN_DL     = R_TN_DL(S.servBS);
L.R_TN_UL     = R_TN_UL(S.servBS);
L.Dprop_TN    = Dprop_TN(S.servBS);

L.SINR_NTN_DL = SINR_NTN_DL(S.servSAT);
L.SINR_NTN_UL = SINR_NTN_UL(S.servSAT);
L.R_NTN_DL    = R_NTN_DL(S.servSAT);
L.R_NTN_UL    = R_NTN_UL(S.servSAT);
L.Dprop_NTN   = Dprop_NTN(S.servSAT);

L.out_TN  = (L.SINR_TN_DL < P.SINRmin);
L.out_NTN = (L.SINR_NTN_DL < P.SINRmin);

end

% ---------------- Helpers ----------------

function g = fading_gain(type, P, domain)
switch lower(type)
    case 'rayleigh'
        h = (randn + 1j*randn)/sqrt(2);
        g = abs(h)^2;

    case 'nakagami'
        % Toolbox-free Nakagami-m power gain:
        % |h|^2 ~ Gamma(m, 1/m). For integer m:
        % sum of m exponentials / m
        if strcmpi(domain,'NTN')
            m = P.NTN.nak_m;
        else
            m = 1;
        end
        m = max(1, round(m));
        acc = 0;
        for i = 1:m
            acc = acc - log(rand); % Exp(1)
        end
        g = acc / m;

    otherwise
        g = 1.0;
end
end

function x = db2lin_local(dB)
x = 10.^(dB/10);
end