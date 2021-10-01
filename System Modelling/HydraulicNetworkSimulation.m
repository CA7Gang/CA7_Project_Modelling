classdef HydraulicNetworkSimulation
    %HydraulicNetworkSimulation This class creates a hydraulic network
    %simulation setup based on a provided graph model of the network
    %topology and a number of component characteristics
    
    properties
        Graph % Graph model of the hydraulic network
        Components % List of system components (hydraulic resistances)
        C % Incidence matrix for components - assigns components to correct edge
        q % Symbolic variable for flows in edges
        lambda % Variable containing expression for pressure drop across edges
    end
    
    methods
        function obj = HydraulicNetworkSimulation(Graph,Components)
            %HydraulicNetworkSimulation Constructs the simulation model
            obj.Graph = Graph;
            obj.Components = Components;
            
            for ii = 1:length(obj.Components.Pipes)
                    obj.lambda(ii) = Components.Pipes(ii).Lambda;
            end
            
            syms q [numel(obj.Graph.edges) 1]
            obj.q = q;
            
            obj.lambda = obj.lambda*abs(obj.q)*obj.q;
            
            
           
            
            
            
            
            
            
            
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

