clear all
close all
Length = 10; rho = 1; Area = 0.5; Diameter = 0.025; g = 9.82; eta = 0.05; Reynolds = 125000; kf = 1.8; h = -3;


fooPipe = PipeComponent(Length,rho,Area,Diameter,g,eta,Reynolds,kf,h); % Every object property except mus and dps should be nonzero

fooValve = ValveComponent('foo',1);

fooPump = PumpComponent(2,4);

H = [1 0 1 0 0 0 0;
    -1 1 0 0 0 0 0;
    0 0 -1 1 1 0 0;
    0 -1 0 -1 0 1 0;
    0 0 0 0 -1 0 -1;
    0 0 0 0 0 -1 1];

fooGraph = GraphModel(H,[1,4],4,[1 6],[2 5]);

for ii = 1:numel(fooGraph.edges)
    Components.Pipes(ii) = fooPipe;
end

fooSim = HydraulicNetworkSimulation(fooGraph,Components);