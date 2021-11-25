clear all
close all

load PumpSimulationData.mat 

for ii = 1:out.logsout.numElements
    try
    assignin('base',out.logsout.getElementNames{ii},out.logsout.getElement(ii)); % Bad practice but I am lazy
    catch 
    foo = out.logsout.getElementNames{ii};
    foo = foo(find(~isspace(foo)));
    assignin('base',foo,out.logsout.getElement(ii));
    clear foo
    end
end

% Do not emulate the above code for any other purposes unless you know
% exactly what you are doing. Mathworks explicitly recommends AGAINST this
% kind of code design for good reason.
%%

subplot(2,1,1)
plot(d1.Values.Time,d1.Values.Data,'b')
hold on
plot(d2.Values.Time,d2.Values.Data,'r')
plot(Ref1.Values.Time,Ref1.Values.Data,'-.b')
plot(Ref2.Values.Time,Ref2.Values.Data,'-.r')
hold off
xlabel('Time [s]','Interpreter','latex')
ylabel('Flow [$\frac{m^3}{s}$]','Interpreter','latex')
legend('$d_1$','$d_2$','Reference $1$','Reference $2$','Interpreter','latex')
FlipMyFuckingLabel(gca)
hold off


subplot(2,1,2)
plot(u1.Values.Time,u1.Values.Data,'b')
hold on
plot(u2.Values.Time,u2.Values.Data,'r')
xlabel('Time [s]','Interpreter','latex')
ylabel('Pump PWM','Interpreter','latex')
legend('$u_1$','$u_2$','Interpreter','latex')
FlipMyFuckingLabel(gca)
hold off


exportgraphics(figure(1),'PumpSimulation.pdf','BackgroundColor','none','ContentType','vector');

function FlipMyFuckingLabel(gca)
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    gca.YLabel.PositionMode='auto';
end