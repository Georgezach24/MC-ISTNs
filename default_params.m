function P = default_params(scenario)

P.dt = 1e-3;
simTime = 20;
P.T = round(simTime / P.dt);

P.areaLen = 20000;    

P.c = 3e8;
P.W = 20e6;
P.N0 = 1e-20;

% UE
P.UE.v = 20;
P.UE.x0 = P.areaLen/2;

% BS
P.BS(1).x = 0;
P.BS(2).x = P.areaLen;

% Satellites
P.SAT(1).x0 = 0.25*P.areaLen;
P.SAT(2).x0 = 0.75*P.areaLen;

for i=1:2
    P.SAT(i).v = 50;
    P.SAT(i).alt = 500e3;
end

% Channel
P.TN.alpha = 3.5;
P.TN.Pdl = 1;
P.TN.Pul = 1;
P.TN.fading = 'rayleigh';

P.NTN.alpha = 2.2;
P.NTN.Pdl = 1;
P.NTN.Pul = 1;
P.NTN.fading = 'nakagami';
P.NTN.nak_m = 2;
P.NTN.extraLoss_dB = 3;

P.SINRmin = 10^(-5/10);

% Handover
P.HO.hyst_dB = 2;
P.HO.hyst = 10^(P.HO.hyst_dB/10);
P.HO.TTT_TN_s = 0.02;
P.HO.TTT_NTN_s = 0.08;
P.HO.exec_TN_s = 0.005;
P.HO.exec_NTN_s = 0.01;
P.HO.alpha = 0.95;

% URLLC
P.URLLC.pktSize = 64*8;
P.URLLC.lambda = 100;

if scenario == 'A'
    P.URLLC.deadline = 0.02;
else
    P.URLLC.deadline = 0.01;
end

end