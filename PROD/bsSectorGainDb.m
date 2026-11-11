function g = bsSectorGainDb(thetaDeg, phiDeg, ant)
% Κέρδος (dBi) κεραίας τομέα προς (θ, φ)· θ ζενιθιακή (90° = οριζόντια), φ ως
% προς τον άξονα του τομέα. Στοιχείο: TR 38.901 Πίν. 7.3-1· κάθετη συστοιχία
% M = 10, d_V = 0.5λ, κλίση 102°: Πίν. 7.8-1.
aV = -min(12*((thetaDeg - 90)/ant.Theta3dBDeg).^2, ant.SlaVDb);
aH = -min(12*(phiDeg/ant.Phi3dBDeg).^2, ant.AmDb);
a  = -min(-(aV + aH), ant.AmDb);
m  = (0:ant.NumVertElements-1);
psi = 2*pi*ant.VertSpacingLambda * (cosd(thetaDeg(:)) - cosd(ant.ElecTiltDeg));
af = abs(sum(exp(1j*psi*m), 2)).^2 / ant.NumVertElements;
g  = ant.ElementGainDbi + a(:) + 10*log10(af);
g  = reshape(g, size(thetaDeg));
end
