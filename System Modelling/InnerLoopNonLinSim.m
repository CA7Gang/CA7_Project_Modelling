clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1000; % Density of water (kg/L)
v = 1.012*10^-6; % Kinematic viscosity of water
q_mean = 1.5/3600; % Mean assumed flow per second (1.5 cubic meters / hr)

kf = 1.8;

PipeIndex = [2 3 5 6 7 8 9 10 12 13];

PumpIndex = [1 14];

ValveIndex = [4 11];

PipeLen = [5 15 15 15 5 15 5 15 20 5]'; 
PipeDiam = 10^-3.*[25 25 25 25 25 15 25 25 15 15]';

eta = 5*10^-5.*ones(10,1); % Pipe roughness factors


% Pump coefficients
a2 = -[0.0355 0.0355]';
a1 = [0.0004 0.0004]';
a0 = [0.0001 0.0001]';


Kv = ones(1,2);


for ii = 1:length(PipeLen)
Pipes(ii) =  PipeComponent(PipeLen(ii),rho,PipeDiam(ii),g,eta(ii),q_mean,v,kf);
end

for ii = 1:length(a2)
Pumps(ii) = PumpComponent(a0(ii),a1(ii),a2(ii));
end

for ii = 1:length(Kv)
Valves(ii) = ValveComponent('Linear',Kv(ii));
end

for ii = 1:length(PipeIndex)
    Components{PipeIndex(ii)} = @(q,w,OD) Pipes(ii).Lambda(q)+Pipes(ii).mu(q,OD)+Pipes(ii).alpha(q,w);
    Inertias(PipeIndex(ii)) = Pipes(ii).J;
end

for ii = 1:length(PumpIndex)
    Components{PumpIndex(ii)} = @(q,w,OD) Pumps(ii).Lambda(q)+Pumps(ii).mu(q,OD)+Pumps(ii).alpha(q,w);
    Inertias(PumpIndex(ii)) = Pumps(ii).J;
end

for ii = 1:length(ValveIndex)
    Components{ValveIndex(ii)} = @(q,w,OD) Valves(ii).Lambda(q)+Valves(ii).mu(q,OD)+Valves(ii).alpha(q,w);
    Inertias(ValveIndex(ii)) = Valves(ii).J;
end

Inertias = diag(Inertias);



% H = [1 0 1 0 0 0 0;
%     -1 1 0 0 0 0 0;
%     0 0 -1 1 1 0 0;
%     0 -1 0 -1 0 1 0;
%     0 0 0 0 -1 0 -1;
%     0 0 0 0 0 -1 1];

 
H = [1	0	0	0	0	0	0	0	0	0	0	0	0	0;
-1	1	0	0	0	0	0	0	0	0	0	0	0	0;
0	-1	1	0	1	0	0	0	0	0	0	0	0	0;
0	0	-1	1	0	1	0	0	0	0	0	0	0	0;
0	0	0	-1	0	0	0	0	0	0	0	0	0	0;
0	0	0	0	-1	0	1	1	0	0	0	0	0	0;
0	0	0	0	0	-1	-1	0	1	0	0	0	0	0;
0	0	0	0	0	0	0	0	-1	1	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	-1	0	0	0;
0	0	0	0	0	0	0	-1	0	0	1	-1	0	0;
0	0	0	0	0	0	0	0	0	-1	0	1	-1	0;
0	0	0	0	0	0	0	0	0	0	0	0	1	-1;
0	0	0	0	0	0	0	0	0	0	0	0	0	1];


 
% NodeHeights = zeros(size(H,1)-1,1)-h0;
NodeHeights = [0; 0; 0; 0; 0.9; 0; 0; 3; 0.9; 0; 0; 0]-h0; 
p0 = 0;
 
 
fooGraph = GraphModel(H,[3,7],13,[1 13],[5 9],8);

% We need to know where pumps and valves are on the FULL graph before any
% reductions
PumpEdges = [1;14];
ValveEdges = [4;11];

fooSim = HydraulicNetworkSimulation(fooGraph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0,rho,g);

%%

 % Collect the pressure equations in one vector (chords at the top)
                Omega_T = sym(zeros(length(fooGraph.spanT),1));
                Omega_C = sym(zeros(length(fooGraph.chords),1));
                
                for ii = 1:length(fooGraph.spanT)
                    Omega_T(ii) = fooSim.Omega_T(ii);
                end
                for ii = 1:length(fooGraph.chords)
                    Omega_C(ii) = fooSim.Omega_C(ii);
                end
                
                Omegas = [Omega_C;Omega_T];

                syms w1 w2 OD1 OD2 p q3 q7 d1 d5 d9 d8
                w = 66; OD = 0.5; tankflow = 0;
                
                
% Omegas = subs(Omegas,[w1 w2],[w w]);
% Omegas = subs(Omegas,[OD1 OD2],[OD OD]);
% Omegas = subs(Omegas,d6,tankflow);

ResistancePart = fooGraph.Phi*Omegas;
HeightPart = (fooGraph.Psi*(fooSim.NodeHeights))*(fooSim.rho*fooSim.g)/(10^5);
PressurePart = (fooGraph.I*(p-0));

dqdt = fooSim.P*(-ResistancePart+HeightPart+PressurePart);

jacobq = jacobian(dqdt,[q3 q7 d1 d5 d9 d8]);
jacobw = jacobian(dqdt,[w1 w2]);

dqdt = subs(dqdt,[w1 w2],[w w]);
dqdt = subs(dqdt,[OD1 OD2],[OD OD]);
dqdt = subs(dqdt,d8,tankflow);

eqpoint = solve(dqdt == 0);
q0 = struct2array(eqpoint);
pt = q0(4);
q0 = [q0(5:6) q0(1:3) 0]';

w = 100; w = 100;
OD = 0.5; OD = 0.5;

if a1(1) == 0
    A = subs(jacobq,[OD1 OD2 q3 q7 d1 d5 d9 d8],[OD OD q0']);
    B = subs(jacobw,symvar(jacobw),[w w]);
else
    A = subs(jacobq,[OD1 OD2 q3 q7 d1 d5 d9 d8 w1 w2],[OD OD q0' w w]);
    B = subs(jacobw,symvar(jacobw),[q0(3:end)' w w]);
end

eig(A)

% Just to check that f(x0) is actually 0
[dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w w],[OD OD],[q0(1) q0(2)],[q0(3) q0(4) q0(5) 0],q0(6),pt);
double(dqdt)

%% simulering
A = double(A);
B = double(B);
C = [0 0 -1 -1 -1 -1];

ts = 0.025;

LinSys = ss(A,B,C,[]);
LinSys = c2d(LinSys,ts,'Tustin');

Ad = LinSys.A;
Bd = LinSys.B;

s = tf('s');

PI_slow = 1.8033*(s+0.0497)/s;

PI_d = c2d(PI_slow,ts,'tustin');

num = cell2mat(PI_d.Numerator); den = cell2mat(PI_d.Denominator);

clear t

f_p = 1/(3600); % Assumed pump frequency

t = 0:ts:(3600)-ts; % Time vector corresponding to 24 hours
% w = (50*sin(2*pi*f_p*t)); % Pump waveforms
% OD =(0.8+sin(pi*f_p*t+t(end)/2))/2; % Valve waveforms
% OD = OD+0.2;


clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures q_lin pt_lin flowchange linpump
q0 = double(q0);
w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% w1 = w; w2 = w; 
% OD1 = OD; OD2 = OD;
OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
% qc = [0;0];
qc = randn(2,1);
% qc = q0(1:2);
% df = [q0(3:end-1);0];
% df = [0;0;0;0]; 
df = randn(4,1);
pt =0*4*rho*g/10^5;
pt_lin = 0;
% d_t = 0;
d_t = 0;

w0 = 66;

d1ref = 0.2;
d2ref = 0.4;

q_lin(:,1) = randn(6,1);

for ii = 1:length(t)
    [dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w1(ii) w2(ii)],[OD1(ii) OD2(ii)],[qc(1) qc(2)],[df(1) df(2) df(3) df(4)],d_t,pt);
    qc(1) = qc(1)+dqdt(1)*ts;
    qc(2) = qc(2)+dqdt(2)*ts;
    df(1) = df(1)+dqdt(3)*ts; 
    df(2) = df(2)+dqdt(4)*ts; 
    df(3) = df(3)+dqdt(5)*ts; 
    d_t = d_t+dqdt(end)*ts;

    % Must obey basic mass conservation if no leaks are assumed
    masscon = -cumsum([df(1:3);d_t]);
    df(4) = masscon(end);

    e1 = d1ref-df(1);
    e2 = d2ref-df(4);


    w1(ii+1) = w0+filter(num,den,e1);
    w2(ii+1) = w0+filter(num,den,e2);

    pt = pt_new;
    flow(:,ii) = [qc;df;d_t];
    flowL = [qc;df(1:3);d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
    
    q_lin(:,ii+1) = Ad*(flowL-q0)+Bd*([w1(ii);w2(ii)]-w0);
    pt_lin(ii+1) = pt_lin(ii)-0.000096*q_lin(end,ii);
    
end

q_lin = q_lin+q0;

for ii = 1:length(q_lin)
    masscon = cumsum(q_lin(3:end,ii));
    linpump(ii) = -(masscon(end));
end

t = seconds(t);
t = minutes(t);
%% Plots



close all

figure(1)
subplot(2,2,1)
plot(t,flow(1,:),'-b')
hold on
plot(t,flow(2,:),'-k')
plot(t,q_lin(1,1:end-1),'-.r')
plot(t,q_lin(2,1:end-1),'-.g')
hold off
legend('Chord 1','Chord 2','Linearised Chord 1','Linearised Chord 2','Interpreter','latex','Location','best')
title('Chord Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)


subplot(2,2,2)
plot(t,flow(3,:),'-b')
hold on
plot(t,flow(6,:),'-k')
plot(t,q_lin(3,1:end-1),'-.r')
plot(t,linpump(1:end-1),'-.g')
hold off
legend('Pump 1','Pump 2','Linearised Pump 1', 'Linearised Pump 2','Interpreter','latex','Location','best')
title('Pump Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

subplot(2,2,3)
plot(t,flow(4,:),'-b')
hold on
plot(t,flow(5,:),'-k')
plot(t,q_lin(4,1:end-1),'-.r')
plot(t,q_lin(5,1:end-1),'-.g')
hold off
legend('Consumer 1','Consumer 2','Linearised Consumer 1','Linearised Consumer 2','Interpreter','latex','Location','best')
title('Consumer Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 5])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)


subplot(2,2,4)
plot(t,flow(end,:))
hold on
plot(t,q_lin(end,1:end-1),'-.r')
hold off
legend('Tank Flow','Linearised Tank Flow','Interpreter','latex','Location','best')
title('Tank Flows','interpreter','latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$q$ [$\frac{m^3}{hr}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

% figure()
% plot(t,pressures)

figure(2)

subplot(3,1,1)
plot(t,tankpres)
hold on
plot(t,pt_lin(1:end-1),'-.r')
legend('Tank Pressure','Linearised Tank Pressure','Interpreter','latex','Location','best')
title('Tank Pressure','Interpreter','Latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$p_\tau$ [bar]','Interpreter','latex')
FlipMyFuckingLabel(gca)


subplot(3,1,2)
plot(t,OD1)
hold on
plot(t,OD2,'-.r')
hold off
title('Valves','Interpreter','Latex')
legend('Valve 1','Valve 2','Interpreter','latex','Location','best')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$OD$ [$\frac{OD}{OD_{max}}$]','Interpreter','latex')
FlipMyFuckingLabel(gca)

subplot(3,1,3)
plot(t,w1)
hold on
plot(t,w2,'-.r')
hold off
legend('Pump 1','Pump 2','Interpreter','latex','Location','best')
title('Pumps','Interpreter','Latex')
xlabel('Time [min]','Interpreter','latex')
xlim([0 t(end)])
ylabel('$\omega$ [PWM]','Interpreter','latex')
FlipMyFuckingLabel(gca)

%% Export the figures cabron

% exportgraphics(figure(1),'DifferentOPFlows.pdf','BackgroundColor','none','ContentType','vector');
% exportgraphics(figure(2),'DifferentOPPressure.pdf','BackgroundColor','none','ContentType','vector');


%% Export subfigures if you really want to

% axes1 = findall(figure(1),'type','axes');
% axes2 = findall(figure(2),'type','axes');
% 
% for ii = 1:length(axes1)
%     strname = axes1(ii).Title.String;
%     strname(strname == ' ') = [];
%     exportgraphics(axes1(ii),strcat(strname,'.pdf'),'BackgroundColor','none','ContentType','vector');
% end
% 
% for ii = 1:length(axes2)
%     strname = axes2(ii).Title.String;
%     strname(strname == ' ') = [];
%     exportgraphics(axes2(ii),strcat(strname,'.pdf'),'BackgroundColor','none','ContentType','vector');
% end


%% Functions and unused code

function FlipMyFuckingLabel(gca)
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    gca.YLabel.PositionMode='auto';
end
