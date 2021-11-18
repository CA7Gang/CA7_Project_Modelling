clear all
close all

%This file requires the "fullfig" functionality from MathWorks community to
%run. Don't say I didn't tell you.


% Load relevant data
Nom = load('NomSysData.mat');
Half = load('HalfNomSysData.mat');
Double = load('DoubleNomSysData.mat');
Out = load('OutDistNomSysData.mat');

%%

t = hours(Out.t);

fullfig
subplot(2,2,1)
plot(t,Nom.x(1,:),'-b')
hold on
plot(t,Half.x(1,:),'color',[0 0.75 0])
plot(t,Double.x(1,:),'-m')
plot(t,Out.x(1,:),'-k')
hold off
legend('Nominal','$200\%$ nominal','$50\%$ nominal','Output disturbance','interpreter','latex','Location','best')
xlabel('Time [hr]','interpreter','latex')
ylabel('$e_{track}$ [bar]','interpreter','latex')
FlipMyFuckingLabel(gca)
ax1 = gca;
hsp1 = get(ax1,'Position');

%

subplot(2,1,2)
plot(t,Nom.y_real,'-b')
hold on
plot(t,Half.y_real,'color',[0 0.75 0])
plot(t,Double.y_real,'-m')
plot(t,Out.y_real,'-k')
plot(t,Out.refval,'--r')
hold off
xlabel('Time [hr]','Interpreter','latex')
ylabel('$p_\tau$ [bar]','Interpreter','latex')
legend('Nominal','$200\%$ nominal','$50\%$ nominal','Output disturbance','Pressure reference','interpreter','latex','Location','best')
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
hsp2 = get(gca, 'Position');                  % Get 'Position' for (2,1,2)
hsp1(3) = hsp2(3);
set(ax1,'Position',hsp1);
FlipMyFuckingLabel(gca)

exportgraphics(figure(1),'LQRTracking.pdf','BackgroundColor','none','ContentType','vector');

%%
fullfig

subplot(3,1,1)
plot(t,Nom.uLQR(1,:),'-b')
hold on
plot(t,Half.uLQR(1,:),'color',[0 0.75 0])
plot(t,Double.uLQR(1,:),'-m')
plot(t,Out.uLQR(1,:),'-k')
hold off
xlabel('Time [hr]','Interpreter','latex')
ylabel('$\Delta u$ [$\frac{m^3}{hr}$]','Interpreter','latex')
legend('Nominal','$200\%$ nominal','$50\%$ nominal','Output disturbance','interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)


subplot(3,1,2)
plot(t(1:end-1),Nom.u(1,:),'b')
hold on
plot(t(1:end-1),Half.u(1,:),'color',[0 0.75 0])
plot(t(1:end-1),Double.u(1,:),'m')
plot(t(1:end-1),Out.u(1,:),'k')
hold off
xlabel('Time [hr]','Interpreter','latex')
ylabel('$u$ [$\frac{m^3}{hr}$]','Interpreter','latex')
legend('Nominal','$200\%$ nominal','$50\%$ nominal','Output disturbance','interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)

%
subplot(3,1,3)
plot(t(1:end-1),Out.dc(1,1:end-1),'-b')
xlabel('Time [hr]','Interpreter','latex')
ylabel('$\delta$ [$\frac{m^3}{hr}$]','Interpreter','latex')
legend('Disturbance','interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)

exportgraphics(figure(2),'LQRControls.pdf','BackgroundColor','none','ContentType','vector');

function FlipMyFuckingLabel(gca)
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    gca.YLabel.PositionMode='auto';
end
