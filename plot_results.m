function plot_results(P, Log)

t = (0:P.T-1) * P.dt;

figure;
plot(t, Log.meanBestTN/1e6, 'LineWidth', 1.2); hold on;
plot(t, Log.meanBestNTN/1e6, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Mean Best Rate (Mbps)');
title('Mean Best TN and NTN Rates Across Users');
legend('Best TN','Best NTN','Location','best');

figure;
plot(t, Log.sumQ_URLLC/1e6, 'LineWidth', 1.2); hold on;
plot(t, Log.sumQ_eMBB/1e6, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Total Queue (Mbits)');
title('Aggregate Queue Occupancy');
legend('URLLC','eMBB','Location','best');

figure;
avgRate = cumsum(Log.totalBitsDL) ./ ((1:P.T)' * P.dt);
plot(t, avgRate/1e6, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Running Avg Delivered Throughput (Mbps)');
title('Aggregate Running Average Throughput');

figure;
plot(t, Log.numHOactive, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Users in HO');
title('Number of Users Under HO Interruption');

figure;
imagesc(t, 1:P.Nue, Log.servingType');
colorbar;
xlabel('Time (s)');
ylabel('User index');
title('Serving Type per User (0=TN, 1=NTN)');


figure;
plot(t, Log.loadTN(:,1), 'LineWidth', 1.2); hold on;
if size(Log.loadTN,2) > 1
    plot(t, Log.loadTN(:,2), 'LineWidth', 1.2);
end
grid on;
xlabel('Time (s)');
ylabel('Users attached');
title('TN Node Load');
legendStrings = arrayfun(@(x) sprintf('TN BS %d',x), 1:size(Log.loadTN,2), 'UniformOutput', false);
legend(legendStrings, 'Location', 'best');

figure;
plot(t, Log.loadNTN(:,1), 'LineWidth', 1.2); hold on;
if size(Log.loadNTN,2) > 1
    plot(t, Log.loadNTN(:,2), 'LineWidth', 1.2);
end
grid on;
xlabel('Time (s)');
ylabel('Users attached');
title('NTN Node Load');
legendStrings = arrayfun(@(x) sprintf('NTN %d',x), 1:size(Log.loadNTN,2), 'UniformOutput', false);
legend(legendStrings, 'Location', 'best');

end