function g = bsSectorGainDb(thetaDeg, phiDeg, ant)
% Κέρδος (dBi) κεραίας τομέα σταθμού βάσης προς κατεύθυνση (θ, φ) στο LCS του
% τομέα: θ ζενιθιακή γωνία (90° = οριζόντια), φ αζιμούθιο από την κατεύθυνση
% του τομέα. Χωρίς μηχανική κλίση το LCS συμπίπτει με το GCS εκτός από την
% περιστροφή κατά την κατεύθυνση του τομέα (TR 36.873 §5.1.3, TR 38.901 §7.1.3).
%
% Μοτίβο στοιχείου, TR 36.873 Πίν. 7.1-1 και TR 38.901 Πίν. 7.3-1:
%   A_V(θ) = -min{12((θ-90°)/θ3dB)^2, SLA_V},  θ3dB = 65°, SLA_V = 30 dB
%   A_H(φ) = -min{12(φ/φ3dB)^2, A_m},          φ3dB = 65°, A_m   = 30 dB
%   A(θ,φ) = -min{-(A_V(θ) + A_H(φ)), A_m},    G_E,max = 8 dBi
% Συστοιχία: M = 10 στοιχεία σε κάθετη στήλη, d_V = 0.5λ, όλα σε μία θύρα με
% βάρη w_m = exp(-j2π(m-1)(d_V/λ)cos θ_etilt)/sqrt(M), θ_etilt = 102° (TR 38.901
% Πίν. 7.8-1 για UMa· TR 36.873 Πίν. 7.1-1). Παράγοντας συστοιχίας:
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
