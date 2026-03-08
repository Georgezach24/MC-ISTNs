function [S, out] = step_serve(P, S, L)

out.bitsDL_user = zeros(P.Nue,1);
out.totalBitsDL = 0;

% -------------------------------------------------
% 1) Count active users per serving node
% -------------------------------------------------
loadTN  = zeros(P.Nbs,1);
loadNTN = zeros(P.Nsat,1);

for u = 1:P.Nue
    if S.UE(u).HO.active
        continue;
    end

    if S.UE(u).servingType == 0
        b = S.UE(u).servingID;
        loadTN(b) = loadTN(b) + 1;
    else
        s = S.UE(u).servingID;
        loadNTN(s) = loadNTN(s) + 1;
    end
end

% Save loads for optional logging
out.loadTN = loadTN;
out.loadNTN = loadNTN;

% -------------------------------------------------
% 2) Serve each user with shared capacity
% -------------------------------------------------
for u = 1:P.Nue

    if S.UE(u).HO.active
        continue;
    end

    if S.UE(u).servingType == 0
        b = S.UE(u).servingID;
        rawCap = L.R_TN_all(u,b) * P.dt;

        if loadTN(b) > 0
            cap = rawCap / loadTN(b);
        else
            cap = rawCap;
        end
    else
        s = S.UE(u).servingID;
        rawCap = L.R_NTN_all(u,s) * P.dt;

        if loadNTN(s) > 0
            cap = rawCap / loadNTN(s);
        else
            cap = rawCap;
        end
    end

    % Serve URLLC first
    serveU = min(S.UE(u).Q.urlLC_bits, cap);
    S.UE(u).Q.urlLC_bits = S.UE(u).Q.urlLC_bits - serveU;
    cap = cap - serveU;

    % Then serve eMBB
    serveE = min(S.UE(u).Q.eMBB_bits, cap);
    S.UE(u).Q.eMBB_bits = S.UE(u).Q.eMBB_bits - serveE;

    bits = serveU + serveE;
    out.bitsDL_user(u) = bits;
    out.totalBitsDL = out.totalBitsDL + bits;
end

end