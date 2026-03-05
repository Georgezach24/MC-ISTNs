function S = traffic_step(P, S, k)

% Mean arrivals per slot
muU = P.traffic.urlLC_mean_bits_per_s * P.dt;
mue = P.traffic.eMBB_mean_bits_per_s  * P.dt;

% Simple Poisson arrivals in bits 
arrU = poissrnd(muU);
arre = poissrnd(mue);

S.Q.urlLC_bits = min(S.Q.urlLC_bits + arrU, P.Q.maxBits);
S.Q.eMBB_bits  = min(S.Q.eMBB_bits  + arre, P.Q.maxBits);

end