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
        q_n % Symbolic variable for final composite state vector
        
        p % Symbolic variable for nodal pressures
        pdot % Symbolic variable for nodal pressure derivative
        
        d % Symbolic variable for nodal demands
        d_f % Symbolic variable for non-zero nodal flows
        d_t % Symbolic variable for tank demands (bilateral)
        
        r_all % Variable containing expression for hydraulic resistance across all edges
        r_T % Variable containing expression for hydraulic resistance across spanning tree edges
        r_C % Variable containing expression for hydraulic resistance across chordal edges
        
        w % Symbolic variable representing pump speeds
        OD % Symbolic variable representing valve opening degrees
        
        PumpEdges % Variable containing indices of edges with pumps
        ValveEdges % Variable containing indices of edges with valves
        Inertias % Variable containing pipe inertias
        
        PF % Symbolic variable containing pressure functions
        
    end
    
    methods
        function obj = HydraulicNetworkSimulation(Graph,Components,PumpEdges,ValveEdges,Inertias)
            %HydraulicNetworkSimulation Constructs the simulation model
            obj.Graph = Graph;
            obj.Components = Components;
            obj.ValveEdges = ValveEdges;
            obj.PumpEdges = PumpEdges;
            obj.Inertias = Inertias;
            
            
            syms q [numel(obj.Graph.edges) 1]
            syms d [numel(obj.Graph.vertices) 1]
            syms pdot [numel(obj.Graph.vertices) 1]
            syms p [numel(obj.Graph.vertices) 1]
            syms w [numel(obj.Graph.producers) 1]
            syms OD [numel(obj.Graph.consumers) 1]
            
            
            for ii = 1:length(obj.Graph.chords)
                obj.r_C{ii,1} = obj.Components{obj.Graph.chords(ii)};
            end
            
             for ii = 1:length(obj.Graph.spanT)
                obj.r_T{ii,1} = obj.Components{obj.Graph.spanT(ii)};
            end
            
            obj.q = q;
            obj.d = d;
            obj.pdot = pdot;
            obj.p = p;
            obj.w = sym(zeros(numel(obj.Graph.edges),1));
            
            for ii = 1:length(obj.PumpEdges)
                obj.w(obj.PumpEdges(ii)) = w(ii);
            end
           
            obj.OD = sym(zeros(numel(obj.Graph.edges),1));
            
            for ii = 1:length(obj.ValveEdges)
                obj.OD(obj.ValveEdges(ii)) = OD(ii);
            end
            
            [obj.q_C, obj.q_T, obj.d_f, obj.d_t,obj.q_n] = ParseGraphInfo(obj);
            obj.PF = ComputePressureDrops(obj);
            obj.q_C = KVL(obj);       
              
            end
            
            
            function [q_C,q_T,d_f,d_t,q_n] = ParseGraphInfo(obj)
                for ii = 1:length(obj.Graph.chords)
                    q_C(ii,1) = obj.q(obj.Graph.chords(ii));
                end
                
                pc = sort([obj.Graph.consumers obj.Graph.producers]);
                
                for ii = 1:length(pc)
                    d_f(ii,1) = obj.d(pc(ii));
                end

                for ii = 1:length(obj.Graph.tanks)
                    d_t(ii,1) = obj.d(obj.Graph.tanks(ii));
                end
                
                H_bar_T = obj.Graph.H_bar_T; H_bar_C = obj.Graph.H_bar_C;
                G_bar = obj.Graph.G_bar; F_bar = obj.Graph.F_bar;
                
                q_T = -inv(H_bar_T)*H_bar_C*q_C + inv(H_bar_T)*F_bar*d_f + H_bar_T*G_bar*d_t;
                
                q_n = [q_C;d_f;d_t];          
            end
        
            function PF = ComputePressureDrops(obj)
                for ii = 1:length(obj.r_T)
                    PF(ii) = obj.r_T{ii}(obj.q_T(ii),obj.w(obj.Graph.spanT(ii)),obj.OD(obj.Graph.spanT(ii)));
                end

            end
            
%             function q_C = KVL(obj)
%                
%                 eqn = obj.lambda_C-obj.Graph.H_bar_C'*inv(obj.Graph.H_bar_T)'*obj.lambda_T == 0; % KVL equation for the network
%                 eqn = subs(eqn,[obj.d(1) obj.d(2) obj.d(3) obj.d(4) obj.d(5) obj.d(6)],[2 -1 0 0 -2 1]); % Just example values for now. Will need som modification.
%                 q_C = vpasolve(eqn,obj.q_C); % Solve with vpasolve since not polynomial, doesn't have closed-form solution.
%             end
    end
end

