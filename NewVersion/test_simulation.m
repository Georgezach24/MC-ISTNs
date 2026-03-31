clc;
clear;

% ------------------ Θέσεις ------------------
% θέση BS
bs_pos = [0, 0;
          0, 30];

% θέση των user
user_pos = [30 60;               
            80 20;
            10 60;
            90 100;
            500 800];              

numUsers = size(user_pos,1); % Αριθμός χρηστών
numBs = size(bs_pos,1);      % Αριθμός BS

% ------------------ Parameters ------------------
simParameters.Carrier = nrCarrierConfig;
simParameters.Carrier.NSizeGrid = 51;           
simParameters.Carrier.SubcarrierSpacing = 30;    
simParameters.Carrier.CyclicPrefix = 'Normal';   

simParameters.CarrierFrequency = 3.5e9;   
simParameters.TxHeight = 25;              
simParameters.TxPower = 40;               
simParameters.RxHeight = 1.5;             
simParameters.RxNoiseFigure = 6;          
simParameters.RxAntTemperature = 290;     

% ------------------ Pathloss model ------------------
simParameters.PathLossModel = '5G-NR';        
simParameters.PathLoss = nrPathLossConfig;
simParameters.PathLoss.Scenario = 'UMa';      
simParameters.PathLoss.EnvironmentHeight = 1; 

% ------------------ Fading channel ------------------
simParameters.DelayProfile = 'TDL-A'; 

if contains(simParameters.DelayProfile,'CDL','IgnoreCase',true)
    channel = nrCDLChannel;
    channel.DelayProfile = simParameters.DelayProfile;
    chInfo = info(channel);
    kFactor = chInfo.KFactorFirstCluster; 
else
    channel = nrTDLChannel;
    channel.DelayProfile = simParameters.DelayProfile;
    chInfo = info(channel);
    kFactor = chInfo.KFactorFirstTap; 
end

simParameters.LOS = kFactor > -Inf;

% ------------------ OFDM info ------------------
waveformInfo = nrOFDMInfo(simParameters.Carrier);

% ------------------ Noise calculation ------------------
kBoltz = physconst('Boltzmann');
NF = 10^(simParameters.RxNoiseFigure/10);
Teq = simParameters.RxAntTemperature + 290*(NF-1); 
N0 = sqrt(kBoltz * waveformInfo.SampleRate * Teq / 2.0);

% ------------------ Bandwidth ------------------
BW_grid = simParameters.Carrier.NSizeGrid * 12 * ...
          simParameters.Carrier.SubcarrierSpacing * 1e3;

% ------------------ Αποθήκευση αποτελεσμάτων ------------------
bestDistanceVec = zeros(numUsers,1);
bestPathLossVec = zeros(numUsers,1);
bestSnrDbVec = zeros(numUsers,1);
capacityMbpsVec = zeros(numUsers,1);
servingBSVec = zeros(numUsers,1);

% Προαιρετικά: αποθήκευση ανά BS
distanceMat = zeros(numUsers,numBs);
pathLossMat = zeros(numUsers,numBs);
snrDbMat = zeros(numUsers,numBs);

% ------------------ Loop για κάθε χρήστη ------------------
for u = 1:numUsers
    
    userBestSNR = -Inf;
    userBestBS = 0;
    userBestDistance = 0;
    userBestPathLoss = 0;
    
    % Loop για κάθε BS
    for b = 1:numBs
        
        % Απόσταση χρήστη - BS b
        d = norm(user_pos(u,:) - bs_pos(b,:));
        distanceMat(u,b) = d;
        
        % Θέσεις πομπού/δέκτη σε 3D
        txPosition = [bs_pos(b,1); bs_pos(b,2); simParameters.TxHeight];
        rxPosition = [user_pos(u,1); user_pos(u,2); simParameters.RxHeight];
        
        % Path loss
        if contains(simParameters.PathLossModel,'5G','IgnoreCase',true)
            pathLoss = nrPathLoss(simParameters.PathLoss, ...
                                  simParameters.CarrierFrequency, ...
                                  simParameters.LOS, ...
                                  txPosition, rxPosition);
        else
            lambda = physconst('LightSpeed') / simParameters.CarrierFrequency;
            pathLoss = fspl(d, lambda);
        end
        
        pathLossMat(u,b) = pathLoss;
        
        % SNR calculation
        fftOccupancy = 12 * simParameters.Carrier.NSizeGrid / waveformInfo.Nfft;
        
        snr_db = (simParameters.TxPower - 30) ...   % dBm -> dBW
                 - pathLoss ...
                 - 10*log10(fftOccupancy) ...
                 - 10*log10(2*N0^2);
        
        snrDbMat(u,b) = snr_db;
        
        fprintf("User %d - BS %d: Distance = %.2f m, PathLoss = %.2f dB, SNR = %.2f dB\n", ...
                u, b, d, pathLoss, snr_db);
        
        % Επιλογή καλύτερου BS με βάση το SNR
        if snr_db > userBestSNR
            userBestSNR = snr_db;
            userBestBS = b;
            userBestDistance = d;
            userBestPathLoss = pathLoss;
        end
    end
    
    % Αποθήκευση καλύτερου BS για τον χρήστη u
    servingBSVec(u) = userBestBS;
    bestDistanceVec(u) = userBestDistance;
    bestPathLossVec(u) = userBestPathLoss;
    bestSnrDbVec(u) = userBestSNR;
    
    % Shannon capacity με βάση το καλύτερο SNR
    snr_lin = 10^(userBestSNR/10);
    capacity = BW_grid * log2(1 + snr_lin);   % bits/s
    capacityMbpsVec(u) = capacity * 1e-6;     % Mbps
    
    fprintf("=> User %d συνδέεται στο BS %d με SNR = %.2f dB\n\n", ...
            u, userBestBS, userBestSNR);
end

% ------------------ Πίνακας αποτελεσμάτων ------------------
userID = (1:numUsers)';

resultsTable = table(userID, servingBSVec, bestDistanceVec, bestPathLossVec, ...
    bestSnrDbVec, capacityMbpsVec, ...
    'VariableNames', {'User','ServingBS','Distance_m','PathLoss_dB','SNR_dB','Capacity_Mbps'});

disp(resultsTable)