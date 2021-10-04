classdef HydraulicNetworkSimulation
    %HydraulicNetworkSimulation This class creates a hydraulic network
    %simulation setup based on a provided graph model of the network
    %topology and a number of component characteristics
    
    properties
        Graph % Graph model of the hydraulic network
        Components % List of system components (hydraulic resistances)
        C % Incidence matrix for components - assigns components to correct edge
        
        q % Symbolic variable for flows in edges
        q_C % Symbolic variable for chordal flows
        q_T % Symbolic variable for tree flows
        
        p % Symbolic variable for nodal pressures
        pdot % Symbolic variable for nodal pressure derivative
        
        d % Symbolic variable for nodal demands
        d_p % Symbolic variable for nodal production
        d_c % Symbolic variable for nodal consumption
        d_t % Symbolic variable for tank demand (bilateral)
        
        r_all % Variable containing expression for hydraulic resistance across all edges
        r_T % Variable containing expression for hydraulic resistance across spanning tree edges
        r_C % Variable containing expression for hydraulic resistance across chordal edges
        
        lambda_C % Symbolic variable containing chordal pressure drops
        lambda_T % Symbolic variable containing spanning tree pressure drops
    end
    
    methods
        function obj = HydraulicNetworkSimulation(Graph,Components)
            %HydraulicNetworkSimulation Constructs the simulation model
            obj.Graph = Graph;
            obj.Components = Components;
            
            for ii = 1:length(obj.Components.Pipes)
                    obj.r_all(ii,1) = Components.Pipes(ii).Lambda;
            end
            
            syms q [numel(obj.Graph.edges) 1]
            syms d [numel(obj.Graph.vertices) 1]
            syms pdot [numel(obj.Graph.vertices) 1]
            syms p [numel(obj.Graph.vertices) 1]
            syms lambda_C
            syms lambda_T
            
            for ii = 1:length(obj.Graph.chords)
                obj.r_C(ii,1) = obj.r_all(obj.Graph.chords(ii));
            end
            
             for ii = 1:length(obj.Graph.spanT)
                obj.r_T(ii,1) = obj.r_all(obj.Graph.spanT(ii));
            end
            
            obj.q = q;
            obj.d = d;
            obj.pdot = pdot;
            obj.p = p;           
            
            [obj.q_C, obj.q_T, obj.d_p, obj.d_c, obj.d_t] = ParseGraphInfo(obj);
            [obj.lambda_C, obj.lambda_T] = ComputePressureDrops(obj);
            q_C = KCL(obj);
            
            
            
                
            
            
            end
            
            
            function [q_C,q_T,d_p,d_c,d_t] = ParseGraphInfo(obj)
                for ii = 1:length(obj.Graph.chords)
                    q_C(ii,1) = obj.q(obj.Graph.chords(ii));
                end
                
                for ii = 1:length(obj.Graph.producers)
                    d_p(ii,1) = obj.d(obj.Graph.producers(ii));
                end
                
                for ii = 1:length(obj.Graph.consumers)
                    d_c(ii,1) = obj.d(obj.Graph.consumers(ii));
                end
                
                for ii = 1:length(obj.Graph.vref)
                    d_t(ii,1) = obj.d(obj.Graph.vref(ii));
                end
                
                H_bar_C = obj.Graph.H_bar_C; H_bar_T = obj.Graph.H_bar_T;
                M_bar_c = obj.Graph.M_bar_c; M_bar_p = obj.Graph.M_bar_p;
                M_bar_t = obj.Graph.M_bar_t;
                
                % Make expression for q_T term by term
                part_chords = -inv(H_bar_T)*H_bar_C*q_C;
                part_producers = inv(H_bar_T)*M_bar_p*d_p;
                part_consumers = inv(H_bar_T)*M_bar_c*d_c;
                part_tank = inv(H_bar_T)*M_bar_t*d_t;
                
                % Combine the terms
                q_T = part_chords + part_producers + part_consumers + part_tank;      
            end
        
            function [lambdaC, lambdaT] = ComputePressureDrops(obj)
                lambdaC = obj.r_C.*abs(obj.q_C).*obj.q_C;
                lambdaT = obj.r_T.*abs(obj.q_T).*obj.q_T;
            end
            
            function q_C = KCL(obj)
                
                eqn = obj.lambda_C-obj.Graph.H_bar_C'*inv(obj.Graph.H_bar_T)'*obj.lambda_T == 0;
                q_C = solve(eqn,obj.q_C);
            end
    end
end

