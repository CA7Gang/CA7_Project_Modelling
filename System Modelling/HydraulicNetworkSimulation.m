classdef HydraulicNetworkSimulation
    %HydraulicNetworkSimulation This class creates a hydraulic network
    %simulation setup based on a provided graph model of the network
    %topology and a number of component characteristics
    
    properties
        Graph % Graph model of the hydraulic network
        Components % List of system components (hydraulic resistances)
        C % Incidence matrix for components - assigns components to correct edge
        P % Inverse of Phi J Phi^T matrix that we don't want to have to invert every time
        
        g % Gravitational acceleration constant
        rho % Density of the working fluid
        
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
        function obj = HydraulicNetworkSimulation(Graph,Components,PumpEdges,ValveEdges,Inertias,NodeHeights,p0,rho,g)
            %HydraulicNetworkSimulation Constructs the simulation model
            obj.Graph = Graph;
            obj.Components = Components;
            obj.ValveEdges = ValveEdges;
            obj.PumpEdges = PumpEdges;
            obj.Inertias = Inertias;
            obj.NodeHeights = NodeHeights;
            obj.p0 = p0;
            obj.P = pinv(obj.Graph.Phi*obj.Inertias*obj.Graph.Phi');
            
            obj.rho = rho;
            obj.g = g;
            
            % Define the symbolic variables we'll need later
            syms q [numel(obj.Graph.edges) 1]
            syms d [numel(obj.Graph.vertices) 1]
            syms pdot [numel(obj.Graph.vertices) 1]
            syms qdot [numel(obj.Graph.edges) 1]
            syms p [numel(obj.Graph.vertices) 1]
            syms w [numel(obj.Graph.producers) 1]
            syms OD [numel(obj.Graph.consumers) 1]
            
            % Get the chord resistance functions
            for ii = 1:length(obj.Graph.chords)
                obj.r_C{ii,1} = obj.Components{obj.Graph.chords(ii)};
            end
            
            % Get the spanning tree resistance functions
            for ii = 1:length(obj.Graph.spanT)
                obj.r_T{ii,1} = obj.Components{obj.Graph.spanT(ii)};
            end
            
            % Map the symbolic variables onto the simulation object
            obj.q = q;
            obj.d = d;
            obj.pdot = pdot;
            obj.qdot = qdot;
            obj.p = p;
            obj.w = sym(zeros(numel(obj.Graph.edges),1));
            obj.OD = sym(zeros(numel(obj.Graph.edges),1));
            
            % Get the edges that have pumps on them and map the symbolic
            % pump speeds in
            for ii = 1:length(obj.PumpEdges)
                obj.w(obj.PumpEdges(ii)) = w(ii);
            end
           
            
            % Get the edges that have valves on them and map the symbolic
            % opening degrees in
            
            for ii = 1:length(obj.ValveEdges)
                obj.OD(obj.ValveEdges(ii)) = OD(ii);
            end
            
            % Collect the 
            [obj.q_C, obj.q_T, obj.d_f, obj.d_t,obj.q_n] = ParseGraphInfo(obj);
            [obj.Omega_T, obj.Omega_C] = ComputePressureDrops(obj);    
            obj.qC_eq = KCL(obj);
              
            end
            
            
            function [q_C,q_T,d_f,d_t,q_n] = ParseGraphInfo(obj)
                
                % Identify the chord flows
                for ii = 1:length(obj.Graph.chords)
                    q_C(ii,1) = obj.q(obj.Graph.chords(ii));
                end
                
                % Make a sorted (ascending) vector of open nodes
                pc = sort([obj.Graph.consumers obj.Graph.producers]);
                
                % Identify the flows at the open nodes in ascending order
                for ii = 1:length(pc)
                    d_f(ii,1) = obj.d(pc(ii));
                end
                
                % Identify the flows at the tank-connected nodes
                for ii = 1:length(obj.Graph.tanks)
                    d_t(ii,1) = obj.d(obj.Graph.tanks(ii));
                end
                
                % Map object properties into local variables so we can
                % actually read stuff
                H_bar_T = obj.Graph.H_bar_T; H_bar_C = obj.Graph.H_bar_C;
                G_bar = obj.Graph.G_bar; F_bar = obj.Graph.F_bar;
                
                % Use the identify for spanning tree flows based on
                % q_C,d_f,d_t
                q_T = -inv(H_bar_T)*H_bar_C*q_C + inv(H_bar_T)*F_bar*d_f + H_bar_T*G_bar*d_t;
                
                % Define the vector of all free flows
                q_n = [q_C;d_f;d_t];          
            end
        
            function [Omega_T,Omega_C] = ComputePressureDrops(obj)
                % Find the pressure drop functions across the spanning tree
                % as a function of the tree flows
                for ii = 1:length(obj.r_T)
                    Omega_T(ii) = obj.r_T{ii}(obj.q_T(ii),obj.w(obj.Graph.spanT(ii)),obj.OD(obj.Graph.spanT(ii)));
                end  
                % Find the pressure drop functions across the chords as a
                % function of the chordal flows
                for ii = 1:length(obj.r_C)
                    Omega_C(ii) = obj.r_C{ii}(obj.q_C(ii),obj.w(obj.Graph.chords(ii)),obj.OD(obj.Graph.chords(ii)));
                end
            end
            
            function eqn = KCL(obj)
                % Make the KCL equation for the graph so we can use it
                % elsewhere
                HbarT = obj.Graph.H_bar_T; HbarC = obj.Graph.H_bar_C;
                eqn = obj.Omega_C'-HbarC'*inv(HbarT)'*obj.Omega_T' == 0;
            end
        
            function [dqdt,pbar,pt_new] = Model_TimeStep(obj,w,OD,df,d_t,pt_old,ts)
                
                % Figure out which edges are pumps and valves so we can
                % substitute in OD and w
                for ii = 1:length(obj.PumpEdges)
                    PumpSubs(ii) = obj.w(obj.PumpEdges(ii));
                end
                
                for ii = 1:length(obj.ValveEdges)
                    ValveSubs(ii) = obj.OD(obj.ValveEdges(ii));
                end
                
                % Substitute OD and w into the pressure loss expression for
                % the spanning tree
                Omega_T = subs(obj.Omega_T,[PumpSubs ValveSubs],[w OD]);
                
                % Find the correct elements for the nodal demands
                for ii = 1:length(obj.d_f)
                    dfvars(ii) = obj.d_f(ii);
                end
                for ii = 1:length(obj.d_t)
                    dtvars(ii) = obj.d_t(ii);
                end
                for ii = 1:length(obj.q_C)
                    qCvars(ii) = obj.q_C(ii);
                end
                
                % Substitute in the nodal demands and find the chord flows
                % numerically
                KCLeq = subs(obj.qC_eq,[dfvars dtvars],[df d_t]);
                q_C = vpasolve(KCLeq,[obj.q_C(1);obj.q_C(2)]);
                
                % Extract the values of q_C from the solution
                names = fieldnames(q_C);
                for ii = 1:length(names)
                    qc(ii) = getfield(q_C,names{ii});
                end
                
                % Substitute the actual flows into the pressure equations
                Omega_T = subs(Omega_T,[qCvars dfvars dtvars],[qc df d_t]);
                Omega_C = subs(obj.Omega_C, [qCvars],[qc]);
                
                
                % Collect the pressure equations in one vector
                Omegas = zeros(max(obj.Graph.edges),1);
                
                for ii = 1:length(obj.Graph.spanT)
                    Omegas(obj.Graph.spanT(ii)) = Omega_T(ii);
                end
                for ii = 1:length(obj.Graph.chords)
                    Omegas(obj.Graph.chords(ii)) = Omega_C(ii);
                end

                % Find the non-reference node pressures
                pbar = inv(obj.Graph.H_bar_T)'*Omega_T' - obj.NodeHeights*(obj.rho*obj.g)/(10^5);
                
                % Get the change in flows
                ResistancePart = obj.Graph.Phi*Omegas;
                HeightPart = (obj.Graph.Psi*(obj.NodeHeights))*(obj.rho*obj.g)/(10^5);
                PressurePart = (obj.Graph.I*(pt_old-0))*1/10^5;
%                 dqdt = -P*obj.Graph.Phi*Omegas+P*(obj.Graph.Psi*(obj.NodeHeights))+P*(obj.Graph.I*(pt_old-0));
                dqdt = obj.P*(-ResistancePart+HeightPart+PressurePart); % Note that P is the inverse of Phi J Phi^T
                

                % Find out where the tanks are
                tankstartindex = numel(obj.Graph.chords)+numel(obj.Graph.consumers)+numel(obj.Graph.producers);
                
                % Get the new tank pressures
                pt_new = pt_old - 0.000096*dqdt(tankstartindex+1:end)*ts; % Need to introduce code here to get the actual tank constant from the graph model
                
               
            
             end
            
            

    end
end

