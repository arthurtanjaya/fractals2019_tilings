classdef Vertex
    %VERTEX For vertices in \Gamma_m
    
    properties (SetAccess = private)
        level      % Level of fineness (this is m in \Gamma_m)
        address    % A word of length m (row vector)
        neighbors  % A vector of addresses of its neighbors
        xycoords   % Coordinates in the Cartesian plane
    end
    
    methods(Static)
        
        function point = Primary(point)
            %GET_PRIMARY Returns primary address
            %   Passthrough if primary already
            %   Twin to get_secondary()
            %   Given $F_wq_i$ with $\left|w\right| = m$, the other address
            %   can be found using the following procedure.
            %   Let $n$ be the least integer such that $w_k = i$ for all
            %   $n < k \le m$. Then $j = w_n$ and $w' = (w_1, \dots,
            %   w_{n-1}, i, j, \dots, j)$ (with $\left|w'\right| = m$).
            q = point(1);
            for i = length(point):-1:2
                if point(i) < q  % Different and currently secondary
                    point(1) = point(i);
                    point(i) = q;
                    point(i+1:end) = point(1);
                    return
                elseif point(i) > q  % Different but already primary
                    return
                end
            end
        end
        
        function point = Secondary(point)
            %GET_SECONDARY Returns secondary address
            %   Passthrough if secondary already
            %   Twin to get_primary()
            %   Given $F_wq_i$ with $\left|w\right| = m$, the other address
            %   can be found using the following procedure.
            %   Let $n$ be the least integer such that $w_k = i$ for all
            %   $n < k \le m$. Then $j = w_n$ and $w' = (w_1, \dots,
            %   w_{n-1}, i, j, \dots, j)$ (with $\left|w'\right| = m$).
            
            q = point(1);
            for i = length(point):-1:2
                if point(i) > q  % Different and currently primary
                    point(1) = point(i);
                    point(i) = q;
                    point(i+1:end) = point(1);
                    return
                elseif point(i) < q  % Different but already secondary
                    return
                end
            end
        end
        
    end
    
    methods
        
        function self = Vertex(address)
            %VERTEX Construct an instance of this class
            %   Input should be a valid address; no error checking (yet?)
            self.level = length(address)-1;  % Off-by-one ugh.
            self.address = address;
            self.address = self.get_primary();
            self = self.set_neighbors();
            self = self.set_xycoords();
        end
        
        function prim = get_primary(self)
            
            point = self.address;
            prim = self.Primary(point);
            
        end
        
        function sec = get_secondary(self)
            
            point = self.address;
            sec = self.Secondary(point);
            
        end
        
        function self = set_neighbors(self)
            %SET_NEIGHBORS Updates addresses of neighbors of self
            
            if all(self.address == self.address(1))
                % Case 1: self is a boundary point
                % Then it has only 2 neighbors
                if self.address(1) == 0
                    self.neighbors = [[self.address(1:end-1) 1];
                        [self.address(1:end-1) 2]];
                elseif self.address(1) == 1
                    self.neighbors = [[self.address(1:end-1) 0];
                        [self.address(1:end-1) 2]];
                else  % self.address(1) == 2
                    self.neighbors = [[self.address(1:end-1) 0];
                        [self.address(1:end-1) 1]];
                end
            else
                % Case 2: self is not a boundary point
                % Then it has 4 neighbors
                % 2 of them come from one cell
                primary = self.get_primary();
                if primary(end) == 0
                    self.neighbors = [[primary(1:end-1) 1];
                        [primary(1:end-1) 2]];
                elseif primary(end) == 1
                    self.neighbors = [[primary(1:end-1) 0];
                        [primary(1:end-1) 2]];
                else  % primary(end) == 2
                    self.neighbors = [[primary(1:end-1) 0];
                        [primary(1:end-1) 1]];
                end
                % The other 2 come from another cell
                secondary = self.get_secondary();
                if secondary(end) == 0
                    self.neighbors = [self.neighbors;
                        [secondary(1:end-1) 1];
                        [secondary(1:end-1) 2]];
                elseif secondary(end) == 1
                    self.neighbors = [self.neighbors;
                        [secondary(1:end-1) 0];
                        [secondary(1:end-1) 2]];
                else  % secondary(end) == 2
                    self.neighbors = [self.neighbors;
                        [secondary(1:end-1) 0];
                        [secondary(1:end-1) 1]];
                end
            end
            
            old = self.neighbors;
            new = [];
            [height,~] = size(old);
            for i = 1:height
                new = [new; Vertex.Primary(old(i,1:end))];
            end
            self.neighbors = new;
        end
        
        function self = set_xycoords(self)
            %SET_NEIGHBORS Updates addresses of neighbors of self
            % Start at a boundary point
            q = [[0.5, 3^0.5/2];
                [0, 0];
                [1, 0]];
            self.xycoords = q(self.address(1)+1, :);  % Off-by-one again
            % Apply IFS based on address of vertex
            self.xycoords = ApplyIFS_2D(self.xycoords, self.address(2:end));
        end
        
        function draw(self, side)
            %DRAW Draws the vertex in xy-plane
            %   side controls the side length of the bounding triangle
            %   Assumes that a figure window is already open and hold is on
            coords = self.xycoords * side;
            plot(coords(1), coords(2), '.r')
        end
        
    end
    
end
