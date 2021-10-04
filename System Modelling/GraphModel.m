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
       
        
        H % Object incidence matrix
        H_T % Spanning tree incidence matrix 
        H_C % Chord incidence matrix
        H_bar % Reduced object incidence matrix (without reference vertex)
        H_bar_T % Reduced spanning tree incidence matrix
        H_bar_C % Reduced chord incidence matrix
        B % Loop matrix
        M_bar_c % Consumer matrix
        M_bar_p % Producer matrix
        M_bar_t % Tank matrix
        
    end
    
    methods
        function obj = GraphModel(H,chords,vref,producers,consumers)
            %GraphModel Construct a graph theory model of a hydraulic
            %network
            if nargin > 0
                
                % Assign all the basic graph information
                obj.H = H;
                obj.chords = chords;
                obj.vref = vref;
                obj.consumers = consumers;
                obj.producers = producers;
                [obj.edges,obj.vertices] = GetGraphProperties(obj);
                obj.spanT = GetSpanningTree(obj);
                
                % Form partitions based on graph info
                [obj.H_T,obj.H_C] = GetHPartitions(obj,obj.H);
                
                % Form reduced matrix and its partitions
                obj.H_bar = obj.H; 
                obj.H_bar(vref,:) = [];
                [obj.H_bar_T,obj.H_bar_C] = GetHPartitions(obj,obj.H_bar);
                
                % Make the loop matrix
                numC = numel(obj.chords); % Total number of chords 
%                 obj.B = [eye(numC,numC), -obj.H_bar_C'*inv(obj.H_bar_T')];
                obj.B = [eye(numC,numC), -obj.H_bar_C'*pinv(obj.H_bar_T')]; % Probably more robust implementation
                
                % Make the demand matrices
                [obj.M_bar_c obj.M_bar_p obj.M_bar_t] = GetDemandMatrices(obj);
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
            notChords = ismember(obj.edges,obj.chords);
            for ii = 1:length(notChords)
                if notChords(ii) == 0
                    spanT = [spanT obj.edges(ii)];
                else
                end
            end
        end
        
        function [Mc, Mp, Mt] = GetDemandMatrices(obj)
            Mc = zeros(numel(obj.vertices),numel(obj.consumers)); % Has vertices-ref rows and consumers columns
            Mp = zeros(numel(obj.vertices),numel(obj.producers)); 
            Mt = zeros(numel(obj.vertices),numel(obj.vref));
            
            for ii = 1:numel(obj.consumers)
                Mc(obj.consumers(ii),ii) = 1;   
            end
            for ii = 1:numel(obj.producers)
                Mp(obj.producers(ii),ii) = 1; 
            end
            Mc(obj.vref,:) = [];
            Mp(obj.vref,:) = [];
            Mt(obj.vref,:) = [];
            
        end
        
    end
end

