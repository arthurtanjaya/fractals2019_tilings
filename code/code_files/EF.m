classdef EF
    
    methods (Static)
        
        function is_eigen = EigenfunctionCheck(level, eigenvalue, eigenfunction)
            % EIGENFUNCTIONCHECK returns true if eigenfunction is an
            % eigenfunction of the laplacian with the given eigenvalue.
            % Eigenfunction must be a row vector of length equal to the
            % number of vertices in Gamma_level.
            gamma = GraphApprox(level);
            npts = length(gamma.vertices);
            lapl = GLaplacian(gamma, eigenfunction);
            laplpar = lapl(4:npts);
            is_eigen = norm(laplpar+eigenvalue*eigenfunction(4:npts)) < 0.001;
        end
        
        function ok = EigenfunctionsCheck(efs, evs)
            %EIGENFUNCTIONSCHECK returns true if each row of efs matrix is
            %an eigenfunction of the laplacian with the eigenvalue equal to
            %the corresponding entry of evs (a column vector).
            [n, plevel] = size(efs);
            level = int8( log(2*plevel-3)/log(3) - 1 );
            arr = zeros(1, n);
            for i = 1:n
                arr(i) = EF.EigenfunctionCheck(level, evs(i, 1), efs(i, :));
            end
            ok = all(arr);
        end
        
        function eigenfunctions = Eigenfunctions(level, eigenvalue)
            % EIGENFUNCTIONS(l, e) is the list of Dirichlet eigenfunctions with
            % eigenvalue e on SG approximation of level l
            npts = length(GraphApprox(level).vertices);
            syms x [1 npts];
            laplvar = GLaplacian(GraphApprox(level), x);
            laplparvar = laplvar(4:npts);
            eigenfunctions = solve([laplparvar == -sym(eigenvalue)*x(4:npts), x(1:3) == [0 0 0], norm(x)^2 ~= 0], x, 'ReturnCondition', true);
        end
        
        function [func, list] = Shrink(old, level, vertex)
            %SHRINK reproduces old - a function on Gamma_n - on the 1-cell
            %adjacent to the given vertex of Gamma_{n+1}
            %(where level must be equal to n+1)
            graphNew = GraphApprox(level).vertices;
            graphOld = GraphApprox(level-1).vertices;
            func = zeros([1,(3^(level+1)+3)/2]);
            list = zeros([1,(3^(level+1)+3)/2]);
            for i = 1 : (3^level+3)/2
                oldAddress = graphOld(i).address;
                newAddress = Vertex.Primary([oldAddress(1), vertex, oldAddress(2:end)]);
                j = GraphApprox.lookup(graphNew, level, newAddress);
                func(j) = old(i);
                list(i) = j;
            end
        end
        
        function new = ShrinkAll(old, level, vertex)
            %SHRINKALL applies Shrink to several functions (represented as
            % columns of the matrix old). All functions are shrunk to the
            % same 1-cell given by vertex.
            length = (3^(level+1)+3)/2;
            [~, num] = size(old{1});
            for i = 1:num
                [func, list] = EF.Shrink(old{1}(1:end,i)', level, vertex);
                new{1}(1:length,i) = func(1:length)';
            end
            new{2} = old{2};
            new{3} = list(old{3});
        end
        
        function [func, list] = Reflect(old, level, vertices)
            graph = GraphApprox(level).vertices;
            func = zeros([1,(3^(level+1)+3)/2]);
            list  = zeros([1,(3^(level+1)+3)/2]);
            for i = 1 : (3^(level+1)+3)/2
                oldAddress = graph(i).address;
                newAddress = oldAddress;
                newAddress(oldAddress == vertices(1)) = vertices(2);
                newAddress(oldAddress == vertices(2)) = vertices(1);
                newAddress = Vertex.Primary(newAddress);
                j = GraphApprox.lookup(graph, level, newAddress);
                func(j) = old(i);
                list(i) = j;
            end
        end
        
        function [func,list] = ExtendLocal(old, level, vertices)
            [func1, list] = EF.Shrink(old, level, vertices(1));
            func2 = EF.Shrink(EF.Reflect(old, level-1,vertices), level, vertices(2));
            func = func1-func2;
        end
        
        function new = ExtendAll(old, level, vertices)
            length = (3^(level+1)+3)/2;
            [~, num] = size(old{1});
            for i = 1:num
                [func, list] = EF.ExtendLocal(old{1}(1:end,i)', level, vertices);
                new{1}(1:length,i) = func(1:length)';
            end
            new{2} = old{2};
            new{3} = list(old{3});
        end
    end
end

