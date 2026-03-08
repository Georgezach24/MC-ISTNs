function S = traffic_step(P, S, k)

muU = P.traffic.urlLC_mean_bits_per_s * P.dt;
muE = P.traffic.eMBB_mean_bits_per_s  * P.dt;

for u = 1:P.Nue
    arrU = poissrnd(muU);
    arrE = poissrnd(muE);

    S.UE(u).Q.urlLC_bits = min(S.UE(u).Q.urlLC_bits + arrU, P.Q.maxBits);
    S.UE(u).Q.eMBB_bits  = min(S.UE(u).Q.eMBB_bits + arrE, P.Q.maxBits);
end
end