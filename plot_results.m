function plot_results(P, Log)

t = Log.t;

figure;
plot(t, Log.R_TN/1e6, "-"); hold on;
plot(t, Log.R_NTN/1e6, "-");
grid on; xlabel("Time (s)"); ylabel("Rate (Mbps)");
title("TN vs NTN instantaneous rate");
legend("TN","NTN","Location","best");

figure;
plot(t, Log.qURLLC/1e6, "-"); hold on;
plot(t, Log.qeMBB/1e6, "-");
grid on; xlabel("Time (s)"); ylabel("Queue (Mbits)");
title("Queues");

figure;
plot(t, cumsum(Log.bitsDL)/(t(end)+P.dt)/1e6, "-");
grid on; xlabel("Time (s)"); ylabel("Avg delivered Mbps");
title("Running average delivered throughput");

end