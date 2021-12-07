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
        Q_n % Coefficient matrix for the free flows

        
        qC_eq % Equation describing KCL
        
        p % Symbolic variable for nodal pressures

        
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

        PumpIndex % Indices of pumps (ordered) for use in simulation
        ValveIndex % Indices of valves (ordered) for use in simulation
        
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
            
            % Collect all the resistance functions in one vector
            obj.r_all = cell(numel(obj.Graph.spanT)+numel(obj.Graph.chords),1);
            
             for ii = 1:length(obj.r_T)
                    obj.r_all{(obj.Graph.spanT(ii))} = obj.r_T{ii};
             end  
             
             for ii = 1:length(obj.r_C)
                    obj.r_all{obj.Graph.chords(ii)} = obj.r_C{ii};
             end

            
            % Map the symbolic variables onto the simulation object
            obj.q = q;
            obj.d = d;
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
            
            % Collect the graph flows and compute the pressure drops
            [obj.q_C, obj.q_T, obj.d_f, obj.d_t,obj.q_n,obj.Q_n] = ParseGraphInfo(obj);
            [obj.Omega_T, obj.Omega_C] = ComputePressureDrops(obj);  
            
            % Solve the static equation to find the steady-state chord flow equations
            obj.qC_eq = KCL(obj);

            for ii = 1:length(obj.PumpEdges)
                obj.PumpIndex(ii) = find(obj.Graph.spanT == obj.PumpEdges(ii));
            end
            
            for ii = 1:length(obj.ValveEdges)
                obj.ValveIndex(ii) = find(obj.Graph.spanT == obj.ValveEdges(ii));
            end


             % Create a function file to speed up simulation

              syms q [numel(obj.Graph.edges) 1];
              syms w [numel(obj.Graph.edges) 1];
              syms OD [numel(obj.Graph.edges) 1];

             ResFunIndices = [obj.Graph.chords obj.Graph.spanT];

             for ii = 1:length(obj.Graph.edges)
                    Omegas(ii) = obj.r_all{ResFunIndices(ii)}(q(ii),w(ii),OD(ii));
             end

             symOmegas = sym(Omegas);

             matlabFunction(symOmegas,'File','fooFile','Vars',{q w OD});

            end

            
            
            function [q_C,q_T,d_f,d_t,q_n,Q_n] = ParseGraphInfo(obj)
                
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
                
                % Use the identity for spanning tree flows based on
                % q_C,d_f,d_t
                q_T = -H_bar_T\H_bar_C*q_C + H_bar_T\F_bar*d_f + H_bar_T\G_bar*d_t;
                
                % Define the vector of all free flows
                q_n = [q_C;d_f;d_t];        
                
                Q_t = nan(length(obj.Graph.spanT),length(q_n));
                
                for jj = 1:length(q_n)
                    for ii = 1:length(q_T)
                        Q_t(ii,jj) = diff(q_T(ii),q_n(jj));
                    end
                end
                
                Q_c = nan(length(obj.Graph.chords));
                
                 for jj = 1:length(q_n)
                    for ii = 1:length(q_C)
                        Q_c(ii,jj) = diff(q_C(ii),q_n(jj));
                    end
                 end
                
                Q_n = [Q_c;Q_t];
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
        
            function [dqdt,pbar,pt_new] = Model_TimeStep(obj,w,OD,qc,df,d_t,pt_old)
                
                % Make vectors describing the mapping of flows, pump speeds
                % and opening degrees into the system equations
                
                q_vec = obj.Q_n*[qc';df';d_t];
                w_vec = zeros(length(obj.Graph.spanT),1);
                o_vec = zeros(length(obj.Graph.spanT),1);
                
                
                spanT = obj.Graph.spanT;
                chords = obj.Graph.chords;
                
%                 for ii = 1:length(obj.PumpEdges)
%                     PumpIndex(ii) = find(spanT == obj.PumpEdges(ii));
%                 end
%                 
%                 for ii = 1:length(obj.ValveEdges)
%                     ValveIndex(ii) = find(spanT == obj.ValveEdges(ii));
%                 end
                
                for ii = 1:length(obj.PumpEdges)
                    w_vec(obj.PumpIndex(ii),ii) = 1;
                end
                
                for ii = 1:length(obj.ValveEdges)
                    o_vec(obj.ValveIndex(ii),ii) = 1;
                end
                
                w_vec = [zeros(size(chords,2)); w_vec];
                
                o_vec = [zeros(size(chords,2)); o_vec];
                
                w_vec = w_vec*w';
                o_vec = o_vec*OD';
                
 
                % Use the edge pressure equations to find the pressures in
                % the full graph and spanning tree.
                
                omegaIndices = [chords spanT];
                
                
%                 for ii = 1:length(obj.Graph.edges)
%                     Omegas(ii) = obj.r_all{omegaIndices(ii)}(q_vec(ii),w_vec(ii),o_vec(ii));
%                 end
                Omegas = fooFile(q_vec,w_vec,o_vec);
                
                Omega_T = Omegas(size(chords,2)+1:end);
               
                
    
                
                % Find the non-reference node pressures
                pbar = obj.Graph.Hinv_bar_T'*Omega_T' - obj.NodeHeights*(obj.rho*obj.g)/(10^5);
                
                % Get the change in flows
                ResistancePart = -obj.Graph.Phi*Omegas';
                HeightPart = (obj.Graph.Psi*(obj.NodeHeights))*(obj.rho*obj.g)/(10^5);
                PressurePart = (obj.Graph.I*(pt_old-0));
                
                dqdt = obj.P*(ResistancePart+HeightPart+PressurePart); % Note that P is the inverse of Phi J Phi^T
                

                % Find out where the tanks are
%                 tankstartindex = numel(obj.Graph.chords)+numel(obj.Graph.consumers)+numel(obj.Graph.producers)+1;
                
                % Get the new tank pressures
                pt_new = pt_old - 0.000096*d_t; % Need to introduce code here to get the actual tank constant from the graph model
%                 pt_new = pt_old - 0.0000096*d_t;
             end
            
            

    end
end

