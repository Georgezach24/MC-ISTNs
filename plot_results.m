function plot_results(P, Log)

t = Log.t;

figure;
plot(t, Log.R_TN/1e6, "-"); hold on;
plot(t, Log.R_NTN/1e6, "-");
grid on;
xlabel("Time (s)");
ylabel("Rate (Mbps)");
title("TN and NTN rates");
legend("TN","NTN","Location","best");

figure;
stairs(t, Log.servingLink, "LineWidth", 1.2);
grid on;
xlabel("Time (s)");
ylabel("Serving Link");
title("Serving link over time");
yticks([0 1]);
yticklabels({'TN','NTN'});

figure;
stairs(t, Log.HOactive, "LineWidth", 1.2);
grid on;
xlabel("Time (s)");
ylabel("HO active");
title("Handover interruption state");

figure;
plot(t, Log.qURLLC/1e6, "-"); hold on;
plot(t, Log.qeMBB/1e6, "-");
grid on;
xlabel("Time (s)");
ylabel("Queue (Mbits)");
title("Queues");
legend("URLLC","eMBB","Location","best");

figure;
plot(t, cumsum(Log.bitsDL)/(t(end)+P.dt)/1e6, "-");
grid on;
xlabel("Time (s)");
ylabel("Avg delivered Mbps");
title("Running average delivered throughput");

end