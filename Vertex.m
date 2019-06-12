classdef Vertex
    %VERTEX For vertices in \Gamma_m
    %   ???

    properties
        level      % Level of fineness (this is m in \Gamma_m)
        address    % A word of length m (row vector)
        neighbors  % A vector of addresses of its neighbors
    end

    methods

        function self = Vertex(address)
            %VERTEX Construct an instance of this class
            %   Input should be a valid address; no error checking (yet?)

            self.level = length(address);
            self.address = address;
            self.neighbors = self.get_neighbors();

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
            for i = 2:self.level
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
            for i = 2:self.level
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

        function neighbors = get_neighbors(self)
            %GET_NEIGHBORS Returns addresses of neighbors of self
            %   ???

            if all(self.address == self.address(1))
                % Case 1: self is a boundary point
                % Then it has only 2 neighbors
                if self.address(1) == 0
                    neighbors = [[self.address(1:end-1) 1];
                                 [self.address(1:end-1) 2]];
                elseif self.address(1) == 1
                    neighbors = [[self.address(1:end-1) 0];
                                 [self.address(1:end-1) 2]];
                else  % self.address(1) == 2
                    neighbors = [[self.address(1:end-1) 0];
                                 [self.address(1:end-1) 1]];
                end
            else
                % Case 2: self is not a boundary point
                % Then it has 4 neighbors
                % 2 of them come from one cell
                primary = self.get_primary();
                if primary(end) == 0
                    neighbors = [[primary(1:end-1) 1];
                                 [primary(1:end-1) 2]];
                elseif primary(end) == 1
                    neighbors = [[primary(1:end-1) 0];
                                 [primary(1:end-1) 2]];
                else  % primary(end) == 2
                    neighbors = [[primary(1:end-1) 0];
                                 [primary(1:end-1) 1]];
                end
                % The other 2 come from another cell
                secondary = self.get_secondary();
                if secondary(end) == 0
                    neighbors = [neighbors;
                                 [secondary(1:end-1) 1];
                                 [secondary(1:end-1) 2]];
                elseif secondary(end) == 1
                    neighbors = [neighbors;
                                 [secondary(1:end-1) 0];
                                 [secondary(1:end-1) 2]];
                else  % secondary(end) == 2
                    neighbors = [neighbors;
                                 [secondary(1:end-1) 0];
                                 [secondary(1:end-1) 1]];
                end
            end

        end

    end

end
