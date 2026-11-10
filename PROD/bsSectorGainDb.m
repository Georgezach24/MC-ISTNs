function g = bsSectorGainDb(thetaDeg, phiDeg, ant)
% Κέρδος (dBi) κεραίας τομέα προς (θ, φ): θ ζενιθιακή γωνία (90° = οριζόντια),
% φ αζιμούθιο ως προς την κατεύθυνση του τομέα (TR 38.901 §7.1.3).
%
% Στοιχείο (TR 38.901 Πίν. 7.3-1):
%   A_V(θ) = -min{12((θ-90°)/θ3dB)^2, SLA_V},  θ3dB = 65°, SLA_V = 30 dB
%   A_H(φ) = -min{12(φ/φ3dB)^2, A_m},          φ3dB = 65°, A_m   = 30 dB
%   A(θ,φ) = -min{-(A_V(θ) + A_H(φ)), A_m},    G_E,max = 8 dBi
% Κάθετη συστοιχία M = 10, d_V = 0.5λ, κλίση θ_etilt = 102° (TR 38.901 Πίν. 7.8-1):
%   AF(θ) = |Σ_m w_m exp(j2π(m-1)(d_V/λ)cos θ)|^2,  μέγιστο M (+10 dB) στο θ_etilt.
aV = -min(12*((thetaDeg - 90)/ant.Theta3dBDeg).^2, ant.SlaVDb);
aH = -min(12*(phiDeg/ant.Phi3dBDeg).^2, ant.AmDb);
a  = -min(-(aV + aH), ant.AmDb);
m  = (0:ant.NumVertElements-1);
psi = 2*pi*ant.VertSpacingLambda * (cosd(thetaDeg(:)) - cosd(ant.ElecTiltDeg));
af = abs(sum(exp(1j*psi*m), 2)).^2 / ant.NumVertElements;
g  = ant.ElementGainDbi + a(:) + 10*log10(af);
g  = reshape(g, size(thetaDeg));
end
