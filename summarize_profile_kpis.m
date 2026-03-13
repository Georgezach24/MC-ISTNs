function summarize_profile_kpis(P, Log)

simTime = P.T * P.dt;

urlIdx  = find(Log.profile == 1);   % URLLC-oriented users
embbIdx = find(Log.profile == 0);   % eMBB-oriented users

fprintf('\n================ PROFILE-BASED KPIs ================\n');

% ---------------------------------------------------
% Helper anonymous for safe handling of empty sets
% ---------------------------------------------------
safeMean = @(x) mean(x(:));
safeFrac = @(cond) mean(cond(:));

% ===================================================
% 1) Average queue per profile
% ===================================================
if ~isempty(urlIdx)
    avgQ_url_URLLC = safeMean(Log.qURLLC_user(:,urlIdx));
    avgQ_embb_URLLC = safeMean(Log.qeMBB_user(:,urlIdx));
    avgQ_total_URLLC = safeMean(Log.qURLLC_user(:,urlIdx) + Log.qeMBB_user(:,urlIdx));
else
    avgQ_url_URLLC = NaN; avgQ_embb_URLLC = NaN; avgQ_total_URLLC = NaN;
end

if ~isempty(embbIdx)
    avgQ_url_eMBB = safeMean(Log.qURLLC_user(:,embbIdx));
    avgQ_embb_eMBB = safeMean(Log.qeMBB_user(:,embbIdx));
    avgQ_total_eMBB = safeMean(Log.qURLLC_user(:,embbIdx) + Log.qeMBB_user(:,embbIdx));
else
    avgQ_url_eMBB = NaN; avgQ_embb_eMBB = NaN; avgQ_total_eMBB = NaN;
end

% ===================================================
% 2) Fraction of time on TN / NTN per profile
% servingType: 0=TN, 1=NTN
% ===================================================
if ~isempty(urlIdx)
    fracTN_URLLC  = safeFrac(Log.servingType(:,urlIdx) == 0);
    fracNTN_URLLC = safeFrac(Log.servingType(:,urlIdx) == 1);
else
    fracTN_URLLC = NaN; fracNTN_URLLC = NaN;
end

if ~isempty(embbIdx)
    fracTN_eMBB  = safeFrac(Log.servingType(:,embbIdx) == 0);
    fracNTN_eMBB = safeFrac(Log.servingType(:,embbIdx) == 1);
else
    fracTN_eMBB = NaN; fracNTN_eMBB = NaN;
end

% ===================================================
% 3) HO count per profile
% ===================================================
if ~isempty(urlIdx)
    totalHO_URLLC = sum(Log.HOcount_user_final(urlIdx));
    avgHOperUser_URLLC = mean(Log.HOcount_user_final(urlIdx));
else
    totalHO_URLLC = NaN; avgHOperUser_URLLC = NaN;
end

if ~isempty(embbIdx)
    totalHO_eMBB = sum(Log.HOcount_user_final(embbIdx));
    avgHOperUser_eMBB = mean(Log.HOcount_user_final(embbIdx));
else
    totalHO_eMBB = NaN; avgHOperUser_eMBB = NaN;
end

% ===================================================
% 4) Average served throughput per profile
% ===================================================
if ~isempty(urlIdx)
    totalBits_URLLC = sum(Log.bitsDL_user(:,urlIdx), 'all');
    thr_URLLC_total = totalBits_URLLC / simTime;              % bits/s
    thr_URLLC_perUser = thr_URLLC_total / numel(urlIdx);      % bits/s/user
else
    thr_URLLC_total = NaN; thr_URLLC_perUser = NaN;
end

if ~isempty(embbIdx)
    totalBits_eMBB = sum(Log.bitsDL_user(:,embbIdx), 'all');
    thr_eMBB_total = totalBits_eMBB / simTime;                % bits/s
    thr_eMBB_perUser = thr_eMBB_total / numel(embbIdx);       % bits/s/user
else
    thr_eMBB_total = NaN; thr_eMBB_perUser = NaN;
end

% ===================================================
% Print results
% ===================================================
fprintf('\nUsers by profile:\n');
fprintf('  URLLC-oriented users = %d\n', numel(urlIdx));
fprintf('  eMBB-oriented users  = %d\n', numel(embbIdx));

fprintf('\nAverage queue per profile:\n');
fprintf('  URLLC-profile: avg URLLC queue = %.3f Mbits, avg eMBB queue = %.3f Mbits, avg total queue = %.3f Mbits\n', ...
    avgQ_url_URLLC/1e6, avgQ_embb_URLLC/1e6, avgQ_total_URLLC/1e6);
fprintf('  eMBB-profile : avg URLLC queue = %.3f Mbits, avg eMBB queue = %.3f Mbits, avg total queue = %.3f Mbits\n', ...
    avgQ_url_eMBB/1e6, avgQ_embb_eMBB/1e6, avgQ_total_eMBB/1e6);

fprintf('\nFraction of time on TN / NTN:\n');
fprintf('  URLLC-profile: TN = %.2f%%, NTN = %.2f%%\n', 100*fracTN_URLLC, 100*fracNTN_URLLC);
fprintf('  eMBB-profile : TN = %.2f%%, NTN = %.2f%%\n', 100*fracTN_eMBB, 100*fracNTN_eMBB);

fprintf('\nHandover counts:\n');
fprintf('  URLLC-profile: total HO = %d, avg HO/user = %.2f\n', totalHO_URLLC, avgHOperUser_URLLC);
fprintf('  eMBB-profile : total HO = %d, avg HO/user = %.2f\n', totalHO_eMBB, avgHOperUser_eMBB);

fprintf('\nServed throughput per profile:\n');
fprintf('  URLLC-profile: total = %.3f Mbps, per-user avg = %.3f Mbps/user\n', ...
    thr_URLLC_total/1e6, thr_URLLC_perUser/1e6);
fprintf('  eMBB-profile : total = %.3f Mbps, per-user avg = %.3f Mbps/user\n', ...
    thr_eMBB_total/1e6, thr_eMBB_perUser/1e6);

fprintf('====================================================\n\n');

% ===================================================
% Optional quick comparison plots
% ===================================================
figure;
bar([avgQ_total_URLLC avgQ_total_eMBB] / 1e6);
set(gca, 'XTickLabel', {'Avg Total Queue'});
ylabel('Mbits');
title('Average Total Queue by Profile');
legend({'URLLC-profile','eMBB-profile'}, 'Location', 'best');

figure;
bar([fracTN_URLLC fracTN_eMBB; fracNTN_URLLC fracNTN_eMBB]');
set(gca, 'XTickLabel', {'URLLC-profile','eMBB-profile'});
ylabel('Fraction of time');
title('Time Spent on TN / NTN by Profile');
legend({'TN','NTN'}, 'Location', 'best');

figure;
bar([avgHOperUser_URLLC avgHOperUser_eMBB]);
set(gca, 'XTickLabel', {'Avg HO per User'});
ylabel('Handovers');
title('Average Handover Count per User by Profile');
legend({'URLLC-profile','eMBB-profile'}, 'Location', 'best');

figure;
bar([thr_URLLC_perUser thr_eMBB_perUser] / 1e6);
set(gca, 'XTickLabel', {'Avg Throughput per User'});
ylabel('Mbps/user');
title('Average Served Throughput per User by Profile');
legend({'URLLC-profile','eMBB-profile'}, 'Location', 'best');

end