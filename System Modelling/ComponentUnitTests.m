clear all
close all

g = 9.82; % Gravitational acceleration
h0 = 0; % Height of the reference node
rho = 1; % Density of water (kg/L)
Reynolds = 125000; 
kf = 1.8;

PipeIndex = [2 4 5 6 7 8 10];

PumpIndex = [1 11];

ValveIndex = [3 9];

PipeLen = [10 20 20 15 10 10 25]'; 
PipeDiam = 10e-3.*[25 25 20 15 25 20 25]';
PipeHeights = zeros(7,1)-h0;
eta = 0.05.*ones(7,1); % Pipe roughness factors

a2 = -[0.0367 0.0367]';
a0 = [7.335e-5 7.335e-5]';

Kv = ones(2,1);


for ii = 1:length(PipeLen)
Pipes(ii) =  PipeComponent(PipeLen(ii),rho,PipeDiam(ii),g,eta(ii),Reynolds,kf,PipeHeights(ii));
end

for ii = 1:length(a2)
Pumps(ii) = PumpComponent(a0(ii),0,a2(ii));
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
 

fooGraph = GraphModel(H,[2,6],10,[1 10],[4 7],6);

PumpEdges = [1;11];
ValveEdges = [3;9];



fooSim = HydraulicNetworkSimulation(fooGraph,Components,PumpEdges,ValveEdges,Inertias);