% Script for plotting of laboratory data

clear all
close all

% load('FullRun_Final.mat');
load('CascadeFullRun_1.mat')
% load('Run_NoSync.mat')
% data = data{2};
Data = PreprocessData(data,12*3600); % Black magic that extracts named signals and interpolates any that need it.

t = seconds(Data.time);
t = minutes(t);
Data.p_tau(1) = Data.p_tau(2);
Data.h_tau(1) = Data.h_tau(2);
t_ploss = 8*3600;

%% Figure 1
fullfig
% subplot(2,1,1)
% plot(t,Data.p_tau,'-b')
% hold on
% plot(t,Data.p_ref,'--k')
% xline(t(t_ploss),'r',{'$50\%$ packet loss begins'},'interpreter','latex');
% hold off
% legend('Tank Pressure','Tank Pressure Reference')
% title('Tank Pressure','interpreter','latex')
% xlabel('Time [min]','Interpreter','latex')
% xlim([0 t(end)])
% ylabel('$p_\tau$ [bar]','Interpreter','latex')
% FlipMyFuckingLabel(gca)

% subplot(2,1,2)
plot(t,Data.h_tau,'-b')
hold on
plot(t,Data.h_ref,'--k')
xline(t(t_ploss),'r',{'$50\%$ packet loss begins'},'interpreter','latex');
hold off
legend('Tank Level','Tank Level Reference')
title('Tank Level','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$h_\tau$ [mm]','Interpreter','latex')
FlipMyFuckingLabel(gca)


%% Figure 2

fullfig
subplot(2,1,1)
plot(t,Data.d1,'-b')
hold on
plot(t,Data.d2,'-r')
plot(t,Data.d1_ref,'--k')
hold off
legend('Pump 1','Pump 2','Pump Reference','interpreter','latex','Location','best','FontSize',14)
title('Pump Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

subplot(2,1,2)
plot(t,Data.d1,'-b')
hold on
plot(t,Data.d2,'-r')
plot(t,Data.d1_ref,'--k')
hold off
legend('Pump 1','Pump 2','Pump Reference','interpreter','latex','Location','best','FontSize',14)
title('Pump Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([t(20000) t(20120)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

%% Figure 3
fullfig
subplot(2,1,1)
plot(t,Data.dc_hat,'-r')
hold on
plot(t,Data.dc_kalm,'-b')
hold off
legend('"Disturbance" Measurement','Kalman Filter','interpreter','latex','Location','best','FontSize',14)
title('Disturbance Estimation','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

subplot(2,1,2)
plot(t,Data.dc1,'-b')
hold on
plot(t,Data.dc2,'-r')
hold off
legend('Consumer 1','Consumer 2','interpreter','latex','Location','best','FontSize',14)
title('Consumer Profiles','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

%% Figure 4
fullfig
figure(4)
hold on
plot(t,Data.dc_kalm,'-b')
plot(t,Data.dc_hat,'-r')
hold off
legend('Kalman Filter','"Disturbance" Measurement','interpreter','latex','Location','best','FontSize',14)
title('Leakage Detection','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([t(14250) t(14750)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

%% Export Figures

exportgraphics(figure(1),'LabResultPics/OuterLoopBad.pdf','BackgroundColor','none','ContentType','vector');
exportgraphics(figure(2),'LabResultPics/InnerLoopBad.pdf','BackgroundColor','none','ContentType','vector');
exportgraphics(figure(3),'LabResultPics/DisturbanceEstimationBad.pdf','BackgroundColor','none','ContentType','vector');
exportgraphics(figure(4),'LabResultPics/LeakageDetectionBad.pdf','BackgroundColor','none','ContentType','vector');