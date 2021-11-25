clear all
close all

load FullSimulation.mat 

for ii = 1:data.numElements
    try
    assignin('base',data.getElementNames{ii},data.getElement(ii)); % Bad practice but I am lazy
    catch 
    foo = data.getElementNames{ii};
    foo = foo(find(~isspace(foo)));
    assignin('base',foo,data.getElement(ii));
    clear foo
    end
end

% Do not emulate the above code for any other purposes unless you know
% exactly what you are doing. Mathworks explicitly recommends AGAINST this
% kind of code design for good reason.
%%

dp = squeeze(d_p.Values.Data);
t = hours(seconds(d_p.Values.Time));
tdp = hours(seconds(dp_ref.Values.Time));

figure(1)
subplot(2,1,1)
plot(t,dp(1,:),'-.b')
hold on
plot(t,dp(2,:),'-.r')
plot(tdp,dp_ref.Values.Data(:,1),'-.k')
plot(tdp,Disturbance.Values.Data,'LineStyle','-.','color',[0 0.75 0])
hold off
xlabel('Time [hr]','Interpreter','latex')
ylabel('Flow [$\frac{m^3}{s}$]','Interpreter','latex')
legend('$d_1$','$d_2$','Flow Reference','Disturbance','Interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)
hold off


subplot(2,1,2)
plot(t,u1.Values.Data,'-.b')
hold on
plot(t,u2.Values.Data,'-.r')
xlabel('Time [hr]','Interpreter','latex')
ylabel('Pump PWM','Interpreter','latex')
legend('$u_1$','$u_2$','Interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)
hold off

figure(2)
plot(t,p_tank.Values.Data,'-.b')
hold on
plot(t,p_ref.Values.Data,'--')
hold off
xlabel('Time [hr]','Interpreter','latex')
ylabel('$p_\tau$ [$bar$]','Interpreter','latex')
legend('Tank Pressure','Pressure Reference','Interpreter','latex','Location','best')
FlipMyFuckingLabel(gca)



exportgraphics(figure(1),'FullSimFlows.pdf','BackgroundColor','none','ContentType','vector');
exportgraphics(figure(2),'FullSimPressures.pdf','BackgroundColor','none','ContentType','vector');

function FlipMyFuckingLabel(gca)
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    gca.YLabel.PositionMode='auto';
end