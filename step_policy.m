function action = step_policy(P, S, L, k)
% 0=TN, 1=NTN, 2=SPLIT, 3=DUP (later)

if L.capTN_bits >= L.capNTN_bits
    action = 0; % TN
else
    action = 1; % NTN
end
end