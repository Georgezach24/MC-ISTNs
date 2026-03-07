function [S, out] = step_serve_handover(P, S, L)

out.bitsDL = 0;

% Αν υπάρχει active HO interruption, δεν μεταδίδεται τίποτα
if S.HO.active
    return;
end

% Capacity based on current serving link
if S.servingLink == 0
    cap = L.capTN_bits;
else
    cap = L.capNTN_bits;
end

% Serve URLLC first
serveU = min(S.Q.urlLC_bits, cap);
S.Q.urlLC_bits = S.Q.urlLC_bits - serveU;
cap = cap - serveU;

% Serve eMBB next
serveE = min(S.Q.eMBB_bits, cap);
S.Q.eMBB_bits = S.Q.eMBB_bits - serveE;

out.bitsDL = serveU + serveE;

end