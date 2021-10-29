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
        qdot % Symbolic variable for derivate of flows in edges
        
        qC_eq % Equation describing KCL
        
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
        
        Omega_T % Symbolic variable containing pressure functions of spanning tree
        Omega_C % Symbolic variable containing pressure functions of chords
        
        NodeHeights % Geodesic heights of the non-reference nodes
        p0 % Initial tank pressures
        
    end
    
    methods
        function obj = HydraulicNetworkSimulation(Graph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0)
            %HydraulicNetworkSimulation Constructs the simulation model
            obj.Graph = Graph;
            obj.Components = Components;
            obj.ValveEdges = ValveEdges;
            obj.PumpEdges = PumpEdges;
            obj.Inertias = Inertias;
            obj.NodeHeights = NodeHeights;
            obj.p0 = p0;
            
            
            syms q [numel(obj.Graph.edges) 1]
            syms d [numel(obj.Graph.vertices) 1]
            syms pdot [numel(obj.Graph.vertices) 1]
            syms qdot [numel(obj.Graph.edges) 1]
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
            obj.qdot = qdot;
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
            [obj.Omega_T, obj.Omega_C] = ComputePressureDrops(obj);    
            obj.qC_eq = KCL(obj);
              
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
        
            function [Omega_T,Omega_C] = ComputePressureDrops(obj)
                for ii = 1:length(obj.r_T)
                    Omega_T(ii) = obj.r_T{ii}(obj.q_T(ii),obj.w(obj.Graph.spanT(ii)),obj.OD(obj.Graph.spanT(ii)));
                end  
                for ii = 1:length(obj.r_C)
                    Omega_C(ii) = obj.r_C{ii}(obj.q_C(ii),obj.w(obj.Graph.chords(ii)),obj.OD(obj.Graph.chords(ii)));
                end
            end
            
            function eqn = KCL(obj)
                HbarT = obj.Graph.H_bar_T; HbarC = obj.Graph.H_bar_C;
                eqn = obj.Omega_C'-HbarC'*inv(HbarT)'*obj.Omega_T' == 0;
            end
        
            function [dqdt,pbar,dptdt] = Model_TimeStep(obj,w,OD,df,d_t,pt)
                for ii = 1:length(obj.PumpEdges)
                    PumpSubs(ii) = obj.w(obj.PumpEdges(ii));
                end
                for ii = 1:length(obj.ValveEdges)
                ValveSubs(ii) = obj.OD(obj.ValveEdges(ii));
                end
                
                Omega_T = subs(obj.Omega_T,[PumpSubs ValveSubs],[w OD]);
                
                for ii = 1:length(obj.d_f)
                    dfvars(ii) = obj.d_f(ii);
                end
                for ii = 1:length(obj.d_t)
                    dtvars(ii) = obj.d_t(ii);
                end
                for ii = 1:length(obj.q_C)
                    qCvars(ii) = obj.q_C(ii);
                end
                
                KCLeq = subs(obj.qC_eq,[dfvars dtvars],[df d_t]);
                q_C = vpasolve(KCLeq,[obj.q_C(1);obj.q_C(2)]);
                
                names = fieldnames(q_C);
                for ii = 1:length(names)
                    qc(ii) = getfield(q_C,names{ii});
                end
                
                Omega_T = subs(Omega_T,[qCvars dfvars dtvars],[qc df d_t])
                Omega_C = subs(obj.Omega_C, [qCvars],[qc])
                
                Omegas = zeros(max(obj.Graph.edges),1);
                
                for ii = 1:length(obj.Graph.spanT)
                    Omegas(obj.Graph.spanT(ii)) = Omega_T(ii);
                end
                for ii = 1:length(obj.Graph.chords)
                    Omegas(obj.Graph.chords(ii)) = Omega_C(ii)
                end

                
                pbar = inv(obj.Graph.H_bar_T)'*Omega_T' - obj.NodeHeights;
                
%                 for ii = 1:length(obj.Graph.spanT)
%                     J(ii) = obj.Inertias(obj.Graph.spanT(ii));
% %                 end
%                 J = diag(J);
                
                P = pinv(obj.Graph.Phi*obj.Inertias*obj.Graph.Phi');
                
                qdot = -P*obj.Graph.Phi*Omegas+P*(obj.Graph.Psi*(obj.NodeHeights))+P*(obj.Graph.I*(pt-0))
                
                
            
             end
            
            

    end
end

