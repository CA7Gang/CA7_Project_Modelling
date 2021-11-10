clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1000; % Density of water (kg/m^3)
v = 1.012*10^-6; % Kinematic viscosity of water
q_mean = 1.5/3600; % Mean assumed flow per second (1.5 cubic meters / hr)
kf = 1.8; % Form factor (see Swamee's WDN book)

PipeIndex = [2 3 4 5];

PumpIndex = [1];

ValveIndex = [6];

PipeLen = [10 20 20 15]'; 
PipeDiam = 10^-3.*[25 25 20 15]';

eta = 5*10^-5.*ones(length(PipeIndex),1); % Pipe roughness factors

% Pump coefficients
a2 = -[0.0355 0.0355]';
a1 = [0.0004 0.0004]';
a0 = [0.0001 0.0001]';

Kv = ones(length(ValveIndex),1);

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

NodeHeights = [0; 0; 0; 0; 0;]-h0; 
p0 = 0;

H = [1 -1 0 0 0 0; 0 1 -1 0 0 0; 0 1 0 -1 0 0; 0 0 1 0 -1 0; 0 0 0 1 -1 0; 0 0 0 0 1 -1];
H = H';

fooGraph = GraphModel(H,[2],6,[1],[6],4);

% We need to know where pumps and valves are on the FULL graph before any
% reductions
PumpEdges = [1];
ValveEdges = [6];

% Make the simulation model
fooSim = HydraulicNetworkSimulation(fooGraph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0,rho,g);

ts = 0.025;

f_p = 2/(24*60); % Assumed pump frequency

t = 0:ts:(10*60); % Time vector corresponding to 24 hours

clear flow tankpres df d_t pt w1 OD1 pressures
w1 = 50*ones(length(t),1);
OD1 = 0.5*ones(length(t),1);
qc = 0;
df = [0;0]; 
pt = 0;
d_t = 0;

for ii = 1:length(t)
    [dqdt,pbar,pt_new] = fooSim.Model_TimeStep(w1(ii),OD1(ii),qc,[df(1) df(2)],d_t,pt);
    qc = qc+dqdt(1)*ts;
    df(1) = df(1)+dqdt(2)*ts; 
    d_t = d_t+dqdt(end)*ts;

    % Must obey basic mass conservation if no leaks are assumed
    masscon = -cumsum([df(1);d_t]);
    df(2) = masscon(end);
    
    flowchange(:,ii) = dqdt;

    pt = pt_new;
    flow(:,ii) = [qc;df;d_t];
    pressures(:,ii) = double(pbar);
    tankpres(ii) = pt;
end

%%

close all
labels = ["Chord","Pump","Valve","Tank"];
figure()
for ii = 1:size(flow,1)
subplot(2,2,ii)
plot(t,flow(ii,:))
legend(labels(ii));
end

%%

figure()
plot(t,pressures)
figure()
subplot(3,1,1)
plot(t,tankpres)
legend('Tank Pressure')
subplot(3,1,2)
plot(t,OD1)
legend('Valve 1','Valve 2')
subplot(3,1,3)
plot(t,w1)

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

                syms w1 OD1 pt d4
                w = 50; OD = 0.5; tankflow = 0;
                
                
Omegas = subs(Omegas,[w1],[w]);
Omegas = subs(Omegas,[OD1],[OD]);
Omegas = subs(Omegas,[d4],[tankflow]);

ResistancePart = fooGraph.Phi*Omegas;
HeightPart = (fooGraph.Psi*(fooSim.NodeHeights))*(fooSim.rho*fooSim.g)/(10^5);
PressurePart = (fooGraph.I*(pt-0));

dqdt = fooSim.P*(-ResistancePart+HeightPart+PressurePart);

eqpoint = solve(fooSim.P*(-ResistancePart+HeightPart+PressurePart) == 0)
struct2array(eqpoint)
