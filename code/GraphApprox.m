classdef GraphApprox
    %GRAPHAPPROX Graph Approximation of SG

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
            
            % Construct vertices
            % Construct V_0
            Z = zeros(1, level+1);  % Off-by-one ugh.
            self.vertices = [Vertex(Z) Vertex(Z+1) Vertex(Z+2)];
            for i = 2:level+1
                % Construct V_{i-1} from V_{i-2}
                new_vertices = Vertex.empty(0, 3^(i-1));
                count = 1;  % To index new_vertices

                % Instead of iterating over the vertices in V_{i-2},
                % we iterate over the i-2 cells instead
                % Cell addressing is unique, so no problems here
                cell_addresses = CartesianProd([0, 1, 2], i-2);
                for cell_index = 1:3^(i-2)
                    if i == 2
                        cell = [];  % Annoying edge case
                    else
                        cell = cell_addresses(cell_index, :);
                    end
                    for head = 0:2
                        address = [head cell];
                        % Pad the end of address with the vertex bit
                        address(end+1:level+1) = head;
                        % Assign addresses in a cyclic fashion
                        new_address = [address(1:i-1) ...
                                       mod(address(i)+1, 3) ...
                                       address(i+1:end)];
                        new_vertices(count) = Vertex(new_address);
                        count = count + 1;
                    end
                end

                self.vertices = [self.vertices new_vertices];  % Save updates
            end

            % Construct edges
            % Two vertices share an edge iff they are in the same m-cell
            % Can postpone this, because we can just count neighbors then
            % divide by 2
            self.edges = 0;

        end

        function draw(self, side)
            %DRAW Draws the graph in xy-plane
            %   side controls the side length of the bounding triangle
            %   Assumes that a figure window is already open and hold is on
            % Draw vertices
            for vertex = self.vertices
                vertex.draw(side)
            end
        end

    end

end
