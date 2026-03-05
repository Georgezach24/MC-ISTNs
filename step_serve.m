function [S, out] = step_serve(P, S, L, action)

out.bitsDL = 0;

% Choose capacity based on action
switch action
    case 0
        cap = L.capTN_bits;
    case 1
        cap = L.capNTN_bits;
    otherwise
        cap = 0;
end

% Serve URLLC first
serveU = min(S.Q.urlLC_bits, cap);
S.Q.urlLC_bits = S.Q.urlLC_bits - serveU;
cap = cap - serveU;

% Then serve eMBB
serveE = min(S.Q.eMBB_bits, cap);
S.Q.eMBB_bits = S.Q.eMBB_bits - serveE;

out.bitsDL = serveU + serveE;

end