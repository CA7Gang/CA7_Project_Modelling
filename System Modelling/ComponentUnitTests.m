clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1000; % Density of water (kg/m^3)
v = 1.012*10^-6; % Kinematic viscosity of water
q_mean = 1.5/3600; % Mean assumed flow per second (1.5 cubic meters / hr)
kf = 1.8; % Form factor (see Swamee's WDN book)

PipeIndex = [2 4 5 6 7 8 10];

PumpIndex = [1 11];

ValveIndex = [3 9];

PipeLen = [10 20 20 15 10 10 25]'; 
PipeDiam = 10^-3.*[25 25 20 15 25 20 25]';

eta = 5*10^-5.*ones(7,1); % Pipe roughness factors

% Pump coefficients
a2 = -[0.0355 0.0355]';
a1 = [0.0004 0.0004]'*0;
a0 = [0.0001 0.0001]';

% Valve coefficients
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
Kv = ones(2,1);


% Make all the components

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

% Incidence matrix

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
 
NodeHeights = [0; 0; 0; 0.9; 0; 3; 0.9; 0; 0]-h0; 
% NodeHeights = [0; 0; 0; 0; 0; 3; 0; 0; 0]-h0; 
p0 = 0;

 
% Make the graph model
fooGraph = GraphModel(H,[2,6],10,[1 10],[4 7],6);

% We need to know where pumps and valves are on the FULL graph before any
% reductions
PumpEdges = [1;11];
ValveEdges = [3;9];

% Make the simulation model
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

                syms w1 w2 OD1 OD2 d6 p q2 q6 d1 d4 d7
                w = 66; OD = 0.5; tankflow = 0;
                
                
% Omegas = subs(Omegas,[w1 w2],[w w]);
% Omegas = subs(Omegas,[OD1 OD2],[OD OD]);
% Omegas = subs(Omegas,d6,tankflow);

ResistancePart = fooGraph.Phi*Omegas;
HeightPart = (fooGraph.Psi*(fooSim.NodeHeights))*(fooSim.rho*fooSim.g)/(10^5);
PressurePart = (fooGraph.I*(p-0));

dqdt = fooSim.P*(-ResistancePart+HeightPart+PressurePart);

jacobq = jacobian(dqdt,[q2 q6 d1 d4 d7 d6])
jacobw = jacobian(dqdt,[w1 w2])

dqdt = subs(dqdt,[w1 w2],[w w]);
dqdt = subs(dqdt,[OD1 OD2],[OD OD]);
dqdt = subs(dqdt,d6,tankflow);

eqpoint = solve(dqdt == 0);
q0 = struct2array(eqpoint);
p_t = q0(4);
q0 = [q0(5:6) q0(1:3) 0]';

w = 66; w = 66;
OD = 0.5; OD = 0.5;

A = subs(jacobq,[OD1 OD2 q2 q6 d1 d4 d7 d6],[OD OD q0'])
B = subs(jacobw,symvar(jacobw),[w w])

eig(A)

% Just to check that f(x0) is actually 0
[dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w w],[OD OD],[q0(1) q0(2)],[q0(3) q0(4) q0(5) 0],q0(6),p_t);
double(dqdt)
             
%% Simulation

A = double(A);
B = double(B);

Asaruch = [-0.4146 0 -0.4169 28.6393 0 0;
    0 -0.2409 -0.0219 0 -3.9529 0;
    -0.1943 0.0050 -0.7858 28.6393 8.2352 0;
    0.1206 0 0.4169 -38.0979 0 0;
    0.0806 0.0619 -0.0968 -0.1361 -24.0953 -0.1361;
    -0.0918 -0.1517 -0.1958 8.8483 3.7352 -0.6104;];

Bsaruch = [0.1742 0; 0.0120 0; 0.2363 0; -0.1742 0; -0.0501 -0.0697; -0.0120 -0.1115];

C = [0 0 -1 -1 -1 -1];

ts = 0.025;

LinSys = ss(A,B,C,[]);
LinSys = c2d(LinSys,ts,'Tustin');



f_p = 2/(24*10); % Assumed pump frequency

t = 0:ts:(2*60); % Time vector corresponding to 24 hours
% w = 50*(1+square(12*pi*f_p*t,50))/2+50*(1+square(12*pi*f_p*t,100))/2; % Pump waveforms
% OD =(1+sin(pi*f_p*t+t(end)/2))/2; % Valve waveforms


clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures q_lin pt_lin flowchange LinPump2
q0 = double(q0);
w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% w1 = w; w2 = w; 
% OD1 = OD; OD2 = OD;
OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
qc = [0;0];
% qc = q0(1:2);
% df = [q0(3:end-1);0];
df = [0;0;0;0]; 
pt =0*4*rho*g/10^5;
pt_lin = 0;
% d_t = 0;
d_t = 0;

w0 = 66;

C = [0 0 -1 -1 -1 -1];

q_lin(:,1) = zeros(6,1);

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
    
    flowchange(:,ii) = dqdt;

    pt = pt_new;
    flow(:,ii) = [qc;df;d_t];
    flowL = [qc;df(1:3);d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
    
    q_lin(:,ii+1) = LinSys.A*(flowL-q0)+LinSys.B*([w1(ii);w2(ii)]-w0);
    pt_lin(ii+1) = pt_lin(ii)-0.000096*q_lin(end,ii);
end

% Plots

q_lin = q_lin+q0;

close all

figure()
subplot(2,2,1)
plot(t,flow(1:2,:))
hold on
plot(t,q_lin(1:2,1:end-1),'--')
hold off
legend('Chord 1','Chord 2','LinChord1','LinChord2')
subplot(2,2,2)
plot(t,flow(3,:))
hold on
plot(t,flow(6,:))
plot(t,q_lin(3,1:end-1),'--')
hold off
legend('Pump 1','Pump 2','LinPump1')
subplot(2,2,3)
plot(t,flow(4:5,:))
hold on
plot(t,q_lin(4:5,1:end-1),'--')
hold off
legend('Consumer 1','Consumer 2','LinConsumer1','LinConsumer2')
subplot(2,2,4)
plot(t,flow(end,:))
hold on
plot(t,q_lin(end,1:end-1),'--r')
hold off
legend('Tank','Linearized Tank')
figure()
plot(t,pressures)
figure()
subplot(3,1,1)
plot(t,tankpres)
hold on
plot(t,pt_lin(1:end-1),'--')
legend('Tank Pressure','Linearized Tank Pressure')
subplot(3,1,2)
plot(t,OD1)
hold on
plot(t,OD2)
hold off
legend('Valve 1','Valve 2')
subplot(3,1,3)
plot(t,w1)
hold on
plot(t,w2)
hold off
legend('Pump 1','Pump 2')

