function S = traffic_step(P, S)

for u = 1:P.Nue
    if S.UE(u).profile == 1
        % URLLC-oriented user
        muU = P.traffic.URLLC_user.url_bits_per_s * P.dt;
        muE = P.traffic.URLLC_user.embb_bits_per_s * P.dt;
    else
        % eMBB-oriented user
        muU = P.traffic.eMBB_user.url_bits_per_s * P.dt;
        muE = P.traffic.eMBB_user.embb_bits_per_s * P.dt;
    end

    % Use poissrnd if available, otherwise Gaussian approx
    arrU = round(muU + sqrt(max(muU,1))*randn());
    arrE = round(muE + sqrt(max(muE,1))*randn());

    arrU = max(arrU,0);
    arrE = max(arrE,0);

    S.UE(u).Q.urlLC_bits = min(S.UE(u).Q.urlLC_bits + arrU, P.Q.maxBits);
    S.UE(u).Q.eMBB_bits  = min(S.UE(u).Q.eMBB_bits + arrE, P.Q.maxBits);
end
end