classdef GraphApprox
    %GRAPHAPPROX Graph Approximation of SG
    %   ???

    properties (SetAccess = private)
        level     % Level of fineness (this is m in \Gamma_m)
        vertices  % Row vector of vertices in the graph
        edges     % Row vector of edges in the graph
    end

    methods

        function self = GraphApprox(level)
            %GRAPHAPPROX Construct an instance of this class
            %   Input should be an int; no error checking (yet?)

            self.level = level;

            % Construct V_0
            Z = zeros(1, level+1);  % Off-by-one ugh.
            self.vertices = [Vertex(Z) Vertex(Z+1) Vertex(Z+2)];
            for i = 2:level+1
                % Construct V_{i-1} from V_{i-2}
                new_vertices = Vertex.empty(0, 3^(i-1));
                count = 1;  % To index new_vertices
                for vertex = self.vertices
                    address = vertex.get_address();
                    % Assign addresses in a cyclic fashion
                    new_address = [address(1:i-1) ...
                                   mod(address(i)+1, 3) ...
                                   address(i+1:end)];
                    new_vertices(count) = Vertex(new_address);
                    count = count + 1;
                end
                self.vertices = [self.vertices new_vertices];  % Save updates
            end

            % Construct edges (HOW?)
            self.edges = 0;

        end

        function level = get_level(self)
            %GET_LEVEL Returns level
            level = self.level;
        end

        function vertices = get_vertices(self)
            %GET_VERTICES Returns vertices
            vertices = self.vertices;
        end

        function edges = get_edges(self)
            %GET_EDGES Returns edges
            edges = self.edges;
        end

    end

end