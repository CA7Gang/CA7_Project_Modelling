clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1; % Density of water (kg/L)
v = 1.012*10^-6; % Kinematic viscosity of water
q_mean = 1.5/3600; % Mean assumed flow per second (1.5 cubic meters / hr)

Reynolds = 25000; 
kf = 1.8;

PipeIndex = [2 4 5 6 7 8 10];

PumpIndex = [1 11];

ValveIndex = [3 9];

PipeLen = [10 20 20 15 10 10 25]'; 
PipeDiam = 10^-3.*[25 25 20 15 25 20 25]';

eta = 5*10^-5.*ones(7,1); % Pipe roughness factors

a2 = -[0.3292 0.3292]';
a1 = [0.038 0.038]';
a0 = [0.0009 0.0009]';

Kv = ones(2,1);


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

H =[1  0  0  0  0  0  0  0  0   0  0;
    -1  1  0  1  0  0  0  0  0   0  0;
     0 -1  1  0  1  0  0  0  0   0  0;
     0  0 -1  0  0  0  0  0  0   0  0;
     0  0  0 -1  0  1  1  0  0   0  0;
     0  0  0  0 -1 -1  0  1  0   0  0;
     0  0  0  0  0  0  0  0 -1   0  0;
     0  0  0  0  0  0 -1  0  1  -1  0;
     0  0  0  0  0  0  0 -1  0   1 -1;
     0  0  0  0  0  0  0  0  0   0  1];
 
% NodeHeights = zeros(size(H,1)-1,1)-h0;
NodeHeights = [0; 0; 0; 0.9; 0; 3; 0.9; 0; 0]-h0; 
p0 = 3;
% p0 = 0;
 
 
fooGraph = GraphModel(H,[2,6],10,[1 10],[4 7],6);

% We need to know where pumps and valves are on the FULL graph before any
% reductions
PumpEdges = [1;11];
ValveEdges = [3;9];

fooSim = HydraulicNetworkSimulation(fooGraph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0,rho,g);

f_p = 2/(24*60); % Assumed pump frequency

t = 0:1:24*60; % Time vector corresponding to 24 minutes
w = 100*(1+square(12*pi*f_p*t,50))/2; % Pump waveforms
OD =(1+cos(pi*f_p*t+t(end)/2))/2; % Valve waveforms
subplot(2,1,1)
plot(t,w)
subplot(2,1,2)
plot(t,OD)
xlim([0 t(end)])
%%

clear flow tankpres df d_t pt w1 w2 OD1 OD2
w1 = w; w2 = w; OD1 = OD; OD2 = OD;
df = [0;0;0;0]; pt = 3; d_t = 0;

for ii = 1:100
    [dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w1(ii) w2(ii)],[OD1(ii) OD2(ii)],[df(1) df(2) df(3) df(4)],d_t,pt,1);
    df(1) = df(1)+dqdt(3); 
    df(2) = df(2)+dqdt(4); 
    df(3) = df(3)+dqdt(5); 
    
    % Must obey mass conservation
%     conflow = cumsum(dqdt);
%     df(4) = -conflow(end);

    df(4) = df(4)+dqdt(6);
    d_t = d_t+dqdt(end);
    pt = pt_new;
    flow(:,ii) = [dqdt(1:2);df;d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
end

%%
close all

figure()
subplot(2,2,1)
plot(t(1):t(length(flow)),flow(1:2,:))
legend('Chord 1','Chord 2')
subplot(2,2,2)
plot(t(1):t(length(flow)),flow(3,:))
hold on
plot(t(1):t(length(flow)),flow(6,:))
hold off
legend('Pump 1','Pump 2')
subplot(2,2,3)
plot(t(1):t(length(flow)),flow(4:5,:))
legend('Consumer 1','Consumer 2')
subplot(2,2,4)
plot(t(1):t(length(flow)),flow(end,:))
legend('Tank')
figure()
plot(t(1):t(length(pressures)),pressures)
figure()
subplot(3,1,1)
plot(t(1):t(length(tankpres)),tankpres)
legend('Tank Pressure')
subplot(3,1,2)
plot(t(1):t(length(tankpres)),OD(1:length(tankpres)))
legend('Valve OD')
subplot(3,1,3)
plot(t(1):t(length(tankpres)),w(1:length(tankpres)))
legend('Pump Speed')





