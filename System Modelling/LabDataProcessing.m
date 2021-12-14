% Script for plotting of laboratory data

clear; close all;

load('FullRun_Final.mat');
% load('CascadeFullRun_1.mat')
% load('/Users/martin/Documents/Git/Repos/CA7_Project/CA7_Project_Modelling/CascadeFullRun_1.mat')
% load('Run_NoSync.mat')
% data = data{2};
Data = PreprocessData(data,12*3600); % Black magic that extracts named signals and interpolates any that need it.

t = seconds(Data.time);
t = minutes(t);
Data.p_tau(1) = Data.p_tau(2);
Data.h_tau(1) = Data.h_tau(2);
t_leak = 4*3600;
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

% Comment these out if poster version is printed and vice versa
% exportgraphics(figure(1),'LabResultPics/OuterLoopBad.pdf','BackgroundColor','none','ContentType','vector');
% exportgraphics(figure(2),'LabResultPics/InnerLoopBad.pdf','BackgroundColor','none','ContentType','vector');
% exportgraphics(figure(3),'LabResultPics/DisturbanceEstimationBad.pdf','BackgroundColor','none','ContentType','vector');
% exportgraphics(figure(4),'LabResultPics/LeakageDetectionBad.pdf','BackgroundColor','none','ContentType','vector');






%% =======================================================================
% Poster figure version
% =======================================================================
% The below figures are dublicates of the above figures with the figures
% adapted to being included in the poster
close all

% Figure settings 
f_xpos = 100; f_ypos = 500; f_width = 900; f_height = 500;
fig_pos_reg	= [f_xpos f_ypos f_width f_height]; fig_color = 'white';

% Figure 1
f = figure(); f.Position = [f_xpos f_ypos f_width (f_height-150)]; f.Color = fig_color;
figures = f;

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
xline(t(t_leak),'--',{'Leakage occurs'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'$50\%$ packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
hold off
legend('Tank Level','Tank Level Reference','FontSize',14)
title('Tank Level','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylim([350 450])
ylabel('$h_\tau$ [mm]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax = gca;
ax.FontSize = 22; 


%% Figure 2

f = figure(); f.Position = fig_pos_reg; f.Color = fig_color;
figures = [figures f]

subplot(2,1,1)
plot(t,Data.d1,'-b')
hold on
plot(t,Data.d2,'-r')
plot(t,Data.d1_ref,'--k')
xline(t(t_leak),'--',{'Leakage occurs'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);

hold off
legend('Pump 1','Pump 2','Pump Reference','interpreter','latex','Location','northeast','FontSize',15)
title('Pump Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylim([0 0.5])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax = gca;
ax.FontSize = 22; 

subplot(2,1,2)
plot(t,Data.d1,'-b')
hold on
plot(t,Data.d2,'-r')
plot(t,Data.d1_ref,'--k')
hold off
legend('Pump 1','Pump 2','Pump Reference','interpreter','latex','Location','best','FontSize',15)
title('Pump Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([t(20000) t(20120)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax = gca;
ax.FontSize = 22; 

%% Figure 3
f = figure(); f.Position = fig_pos_reg; f.Color = fig_color;
figures = [figures f]

subplot(2,1,1)
plot(t,Data.dc_hat,'-r')
hold on
plot(t,Data.dc_kalm,'-b')
xline(t(t_leak),'--',{'Leakage occurs'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
hold off
legend('"Disturbance" Measurement','Kalman Filter','interpreter','latex','Location','northeast','FontSize',14.5)
title('Disturbance Estimation','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylim([0 0.7])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax1 = gca;
ax1.FontSize = 22; 

subplot(2,1,2)
plot(t,Data.dc1,'-b')
hold on
plot(t,Data.dc2,'-r')
xline(t(t_leak),'--',{'Leakage occurs'}, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
hold off
legend('Consumer 1','Consumer 2','interpreter','latex','Location','best','FontSize',14)
title('Consumer Profiles','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylim([0 0.25])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax2 = gca;
ax2.FontSize = 22;

linkaxes([ax1 ax2], 'x')


%% Figure 4

f = figure(); f.Position = [f_xpos f_ypos f_width (f_height-150)]; f.Color = fig_color;
figures = [figures f]

hold on
plot(t,Data.dc_kalm,'-b')
plot(t,Data.dc_hat,'-r')
hold off
legend('Kalman Filter','"Disturbance" Measurement','interpreter','latex','Location','best','FontSize',15)
title('Leakage Detection','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([t(14250) t(14750)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax = gca;
ax.FontSize = 22; 

%% Figure 5
f = figure(); f.Position = fig_pos_reg; f.Color = fig_color;
figures = [figures f]

subplot(2,1,1)
plot(t,Data.dc_hat,'-r')
hold on
plot(t,Data.dc_kalm,'-b')
xline(t(t_leak),'--',{'Leakage occurs'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
hold off
legend('"Disturbance" Measurement','Kalman Filter','interpreter','latex','Location','northeast','FontSize',14.5)
title('Disturbance Estimation','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([(4*60) (8*60)])
% xlim([0 120])
ylim([0 0.5])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax1 = gca;
ax1.FontSize = 22; 

subplot(2,1,2)
plot(t,Data.dc1,'-b')
hold on
plot(t,Data.dc2,'-r')
xline(t(t_leak),'--',{'Leakage occurs'}, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
xline(t(t_ploss),'--',{'packet loss begins'}, 'LineWidth', 2, 'LabelOrientation', 'horizontal', 'LabelHorizontalAlignment', 'left', 'interpreter','latex','FontSize',17);
hold off
legend('Consumer 1','Consumer 2','interpreter','latex','Location','best','FontSize',14)
title('Consumer Profiles','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
% xlim([0 120])
xlim([(4*60) (6*60)])
ylim([0 0.25])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)
ax2 = gca;
ax2.FontSize = 22;

linkaxes([ax1 ax2], 'x')

%% Export figures

% Just change savepath to whichever fits you!
savepath = '/Users/martin/Documents/Git/Repos/CA7_Writings/CA7_Writings_SciPaper/Poster/Figures'
filename = ["OuterLoop_Poster.pdf"; "InnerLoop_Poster.pdf" ;"DisturbanceEstimation_Poster.pdf"; "LeakageDetection_Poster.pdf"; "DisturbanceEstimationZoom_Poster.pdf"]
for i=1:length(figures)
    f = fullfile(savepath, append(filename(i)))

    exportgraphics(figures(i), f,'Resolution', 10) % Resolution doesn't
	% if export filetype is png
end
