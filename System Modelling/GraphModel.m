classdef GraphModel
    %GraphModel Graph Theory-based model of a hydraulic network
    
    properties
        chords % Chords of graph theory model
        vref % Model reference vertex
        spanT % Model spanning tree
        edges % Edges of the model
        vertices % Vertices of the model
        consumers % Open vertices that consume flow
        producers % Open vertices that produce flow
        tanks % Vertices connected to tanks
       
        
        H % Object incidence matrix
        H_T % Spanning tree incidence matrix 
        H_C % Chord incidence matrix
        H_bar % Reduced object incidence matrix (without reference vertex)
        H_bar_T % Reduced spanning tree incidence matrix
        H_bar_C % Reduced chord incidence matrix
        Hinv_bar_T % Inverted spanning tree matrix (for efficiency)
        B % Loop matrix
        F_bar % Non-tank demand matrix
        G_bar % Tank demand matrix
        Phi % Flow block matrix in final equation
        Psi % Pressure block matrix in final equation
        I % Identity matrix for pressures in final equation
        
    end
    
    methods
        function obj = GraphModel(H,chords,vref,producers,consumers,tanks)
            %GraphModel Construct a graph theory model of a hydraulic
            %network
            if nargin > 0
                
                % Assign all the basic graph information
                obj.H = H;
                obj.chords = chords;
                obj.vref = vref;
                obj.consumers = consumers;
                obj.producers = producers;
                obj.tanks = tanks;
                [obj.edges,obj.vertices] = GetGraphProperties(obj);
                obj.spanT = GetSpanningTree(obj);
                
                % Form partitions based on graph info
                [obj.H_T,obj.H_C] = GetHPartitions(obj,obj.H);
                
                % Form reduced matrix and its partitions
                obj.H_bar = obj.H; 
                obj.H_bar(vref,:) = [];
                [obj.H_bar_T,obj.H_bar_C] = GetHPartitions(obj,obj.H_bar);

                obj.Hinv_bar_T = inv(obj.H_bar_T);
                
                % Make the loop matrix
                numC = numel(obj.chords); % Total number of chords 
                obj.B = [eye(numC,numC), -obj.H_bar_C'*inv(obj.H_bar_T')];
%                 obj.B = [eye(numC,numC), -obj.H_bar_C'*pinv(obj.H_bar_T')]; % Probably more robust implementation
%                 
                % Make the demand matrices
                [obj.F_bar,obj.G_bar] = GetDemandMatrices(obj);
                
                % Make the final block matrices
                [obj.Phi, obj.Psi, obj.I] = MakeBlockMatrices(obj);
            end     
                       
        end
        
        function [edges,vertices] = GetGraphProperties(obj)
            %GetGraphProperties Gets the vertices and edges of the graph
            %from the H matrix.
            edges = 1:1:size(obj.H,2);
            vertices = 1:1:size(obj.H,1);
        end
        
        function [H_T,H_C] = GetHPartitions(obj,H)
            H_T = [];  H_C = [];
            % Form tree partition via spanning tree indices
            for ii = 1:length(obj.spanT)
                H_T = [H_T H(:,obj.spanT(ii))];
            end
            % Form chord partition via chord indices
            for ii = 1:length(obj.chords)
                H_C = [H_C H(:,obj.chords(ii))];
            end
        end
        
        function spanT = GetSpanningTree(obj)
            spanT = [];
            notChords = ismember(obj.edges,obj.chords); % Make a vector that contains all the elements of the edges vector that are NOT chords
            for ii = 1:length(notChords)
                if notChords(ii) == 0 % Check if the current edge is a chord
                    spanT = [spanT obj.edges(ii)]; % If not a chord, append to the spanning tree vector
                else
                end
            end
        end
        
        function [F, G] = GetDemandMatrices(obj)
            F = zeros(numel(obj.vertices),numel(obj.producers)+numel(obj.consumers));
            G = zeros(numel(obj.vertices),numel(obj.tanks));
            
            pc = sort([obj.producers obj.consumers]); % Make a descending-sorted list of the producer and consumer nodes
            
            for ii = 1:length(pc)
                F(pc(ii),ii) = 1; % Assign ones to the locations in F that have non-zero nodal demand (producers and consumers)
            end
            
            for ii = 1:length(obj.tanks)
                G(obj.tanks(ii)) = 1; % Assign ones to the locations in G that are tanks
            end
            
            % Kill the row corresponding to the reference node in both F
            % and G
            G(obj.vref,:) = [];
            F(obj.vref,:) = [];
        end
          
        function [Phi,Psi,I] = MakeBlockMatrices(obj)
            numC = numel(obj.chords); % Get the number of chords
            numFlows = numel(obj.consumers)+numel(obj.producers); % Get the correct total number of non-zero nodal demands
            numTanks = numel(obj.tanks); % Get the number of tanks 
            numNodes = numel(obj.vertices); % Get the total number of nodes
            
            Phi = [eye(numC,numC),-obj.H_bar_C'*inv(obj.H_bar_T)'; % Collect the flow equations in one big matrix
            zeros(numFlows,numC),obj.F_bar'*inv(obj.H_bar_T)';
            zeros(numTanks,numC),obj.G_bar'*inv(obj.H_bar_T)'];

            Psi = [zeros(numC,numNodes-1);obj.F_bar';obj.G_bar']; % Collect the height equations in one big matrix
            
            I = [zeros(numC,numTanks); zeros(numFlows,numTanks); eye(numTanks)]; % Collect the pressure equations in one big matrix
            
            % Make a sorted vector of consumers and producers
            pc = sort([obj.consumers obj.producers]);
            pc = [pc obj.tanks];
                        
            % Identify the flow that needs to die (reference flow) assuming
            % chord flows are always listed at the top
            somebodyKillMePlease = find(pc == obj.vref)+numC;
            
            % Ceterum censeo reference fluxus delendam esse
            Psi(somebodyKillMePlease,:) = [];
            Phi(somebodyKillMePlease,:) = [];
            I(somebodyKillMePlease,:) = [];
        end
    end
end

