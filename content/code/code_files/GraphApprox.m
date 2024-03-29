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
            self.edges = zeros([0,2]);
            
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
            for v = 1:(3^(level+1)+3)/2
                [nnbrs, ~] = size(self.vertices(v).neighbors);
                for n = 1:nnbrs
                    nbr = GraphApprox.lookup(self.vertices, level, self.vertices(v).neighbors(n, 1:end));
                    self.edges = [self.edges; v, nbr];
                end
            end
          

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
    
    methods(Static)
        
        function idx = lookup(g,level,address)
            for i = 1: (3^(level+1) + 3)/2
                if (all(g(i).address == address))
                    idx = i;
                end
            end
        end
        
    end
    
end
