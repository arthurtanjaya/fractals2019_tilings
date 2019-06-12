classdef Vertex
    %VERTEX For vertices in \Gamma_m

    properties (SetAccess = private)
        level      % Level of fineness (this is m in \Gamma_m)
        address    % A word of length m (row vector)
        neighbors  % A vector of addresses of its neighbors
    end

    methods

        function self = Vertex(address)
            %VERTEX Construct an instance of this class
            %   Input should be a valid address; no error checking (yet?)
            self.level = length(address)-1;  % Off-by-one ugh.
            self.address = address;
            self.address = self.get_primary();
            self = self.set_neighbors();
        end

        function level = get_level(self)
            %GET_LEVEL Returns level
            level = self.level;
        end

        function address = get_address(self)
            %GET_ADDRESS Returns address
            address = self.address;
        end

        function neighbors = get_neighbors(self)
            %GET_NEIGHBORS Returns neighbors
            neighbors = self.neighbors;
        end

        function point = get_primary(self)
            %GET_PRIMARY Returns primary address
            %   Passthrough if primary already
            %   Twin to get_secondary()
            %   Adapted from https://github.com/seraphinalee/fractals

            point = self.address;  % TODO Is this a shallow or deep copy???
            % Looks at the first values...
            q = point(1);
            % For each listed map...
            for i = 2:self.level+1
                % Looks for the first point that departs from q
                if point(i) ~= q
                    % and sees what direction it goes in,
                    % correcting if necessary
                    if isequal(point(i-1:i), [1 0])
                        point(i-1:i) = [0 1];
                        point(1:i-1) = point(i-1);
                    elseif isequal(point(i-1:i), [2 1])
                        point(i-1:i) = [1 2];
                        point(1:i-1) = point(i-1);
                    elseif isequal(point(i-1:i), [0 2])
                        point(i-1:i) = [2 0];
                        point(1:i-1) = point(i-1);
                    end
                    return
                end
            end

        end

        function point = get_secondary(self)
            %GET_SECONDARY Returns secondary address
            %   Passthrough if secondary already
            %   Twin to get_primary()
            %   Adapted from https://github.com/seraphinalee/fractals

            point = self.address;  % TODO Is this a shallow or deep copy???
            % Looks at the first values...
            q = point(1);
            % For each listed map...
            for i = 2:self.level+1
                % Looks for the first point that departs from q
                if point(i) ~= q
                    % and sees what direction it goes in,
                    % correcting if necessary
                    if isequal(point(i-1:i), [0 1])
                        point(i-1:i) = [1 0];
                        point(1:i-1) = point(i-1);
                    elseif isequal(point(i-1:i), [1 2])
                        point(i-1:i) = [2 1];
                        point(1:i-1) = point(i-1);
                    elseif isequal(point(i-1:i), [2 0])
                        point(i-1:i) = [0 2];
                        point(1:i-1) = point(i-1);
                    end
                    return
                end
            end

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

        end

        function draw(self, side)
        %DRAW Draws the vertex in xy-plane
        %   side controls the side length of the bounding triangle
        %   Assumes that a figure window is already open and hold is on

        % Start at a boundary point
        q = [[0.5, 3^0.5/2];
             [0, 0];
             [1, 0]];
        coord = q(self.address(1)+1, :) * side;  % Off-by-one again

        % Apply IFS based on address of vertex
        coord = ApplyIFS_2D(coord, self.address(2:end));

        % Plot!
        plot(coord(1), coord(2), '.')

        end

    end

end
