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



ts = 0.025;

f_p = 2/(24*60); % Assumed pump frequency

t = 0:ts:(1*60); % Time vector corresponding to 24 hours
% w = 50*(1+square(12*pi*f_p*t,50))/2+50*(1+square(12*pi*f_p*t,100))/2; % Pump waveforms
% OD =(1+sin(pi*f_p*t+t(end)/2))/2; % Valve waveforms


clear flow tankpres df d_t pt w1 w2 OD1 OD2 pressures
w1 = 66*ones(length(t),1); w2 = 66*ones(length(t),1);
% w1 = w; w2 = w; 
% OD1 = OD; OD2 = OD;
OD1 = 0.5*ones(length(t),1); OD2 = 0.5*ones(length(t),1);
qc = [0;0];
df = [0;0;0;0]; 
pt =0*4*rho*g/10^5;
d_t = 0;

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
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
end
%
close all

figure()
subplot(2,2,1)
plot(t,flow(1:2,:))
legend('Chord 1','Chord 2')
subplot(2,2,2)
plot(t,flow(3,:))
hold on
plot(t,flow(6,:))
hold off
legend('Pump 1','Pump 2')
subplot(2,2,3)
plot(t,flow(4:5,:))
legend('Consumer 1','Consumer 2')
subplot(2,2,4)
plot(t,flow(end,:))
legend('Tank')
figure()
plot(t,pressures)
figure()
subplot(3,1,1)
plot(t,tankpres)
legend('Tank Pressure')
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
figure()
plot(t,flowchange(end,:))


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

                syms w1 w2 OD1 OD2 d6 pt
                w = 66; OD = 0.5; tankflow = 0;
                
                
Omegas = subs(Omegas,[w1 w2],[w w]);
Omegas = subs(Omegas,[OD1 OD2],[OD OD]);
Omegas = subs(Omegas,d6,tankflow);

 ResistancePart = fooGraph.Phi*Omegas;
HeightPart = (fooGraph.Psi*(fooSim.NodeHeights))*(fooSim.rho*fooSim.g)/(10^5);
PressurePart = (fooGraph.I*(pt-0));

dqdt = fooSim.P*(-ResistancePart+HeightPart+PressurePart);

eqpoint = solve(dqdt == 0)
q0 = struct2array(eqpoint);
pt = q0(4);
q0 = [q0(5:6) q0(1:3) 0 0]';
% q0(end-1) = -(q0(3)+q0(4)+q0(5)); % Mass conservation

% Just to check that f(x0) is actually 0
% [dqdt,pbar,pt_new] = fooSim.Model_TimeStep([w w],[OD OD],[q0(1) q0(2)],[q0(3) q0(4) q0(5) 0],q0(6),pt);
% double(dqdt)


syms q w OD p

%Create d_omega wrt qn for chords
for i=1:size(fooSim.r_C,1)
r_q_taylor(i)= diff(fooSim.r_C{i}(q,w,OD),q);
end
%Create d_omega wrt qn for tree
for i=1:size(fooSim.r_T,1)
r_q_taylor(i+2)= diff(fooSim.r_T{i}(q,w,OD),q);
end

%Create d_omega wrt omega for chords
for i=1:size(fooSim.r_C,1)
r_w_taylor(i)= diff(fooSim.r_C{i}(q,w,OD),w);
end
%Create d_omega wrt omega for chords
for i=1:size(fooSim.r_T,1)
r_w_taylor(i+2)= diff(fooSim.r_T{i}(q,w,OD),w);
end

w0 = 66;
OD0 = 0.5;
Q_n = fooSim.Q_n;

q_vector = Q_n*q0;
OD_vector = [0 0 0 OD0 0 0 0 0 OD0 0 0];
w_vector = zeros(11,1);
w_vector(3,1) = w0; w_vector(end,1) = w0;

q0(end-1) = [];
Q_n(:,end-1) = [];

clear r_q r_w A B

% for jj = 1:length(q0)
%     for i = 1:length(r_q_taylor)
%         r_q(i,jj) = subs(r_q_taylor(i),[q,OD,w],[q0(jj)*Q_n(i,jj),OD_vector(i),w_vector(i)]);
%     end
% end

for jj = 1:length(q0)
    for i = 1:length(r_q_taylor)
        r_q(i,jj) = subs(r_q_taylor(i),[q,OD,w],[q0(jj)*Q_n(i,jj),OD_vector(i),w_vector(i)]);
    end
end

w_vector = zeros(11,2);
w_vector(3,1) = w0; w_vector(end,2) = w0;

for jj = 1:length(PumpIndex)
    for i = 1:length(r_w_taylor)
        r_w(i,jj) = subs(r_w_taylor(i),[q,OD,w],[q_vector(i),OD_vector(i),w_vector(i,jj)]);
    end
end


A = double(-fooSim.P*fooGraph.Phi*r_q)
B = double(-fooSim.P*fooGraph.Phi*r_w)
              
%%
Asaruch = [-0.4146 0 -0.4169 28.6393 0 0 0;
           0 -0.2409 -0.0219 0 -3.9529 0 -3.6672;
           -0.1943 0.0050 -0.7858 28.6393 8.2352 0 -1.2224;
           0.1206 0 0.4169 -38.0979 0 0 -5.6719;
           0.0806 0.0619 -0.0968 -0.1361 -24.0953 -0.1361 -2.4448;
           -0.0918 -0.1517 -0.1958 8.8483 3.7352 -0.6104 20.6828]
       
Asaruch(:,end) = []
       
MatrixWeShouldGet = pinv(-fooSim.P*fooGraph.Phi)*Asaruch
           