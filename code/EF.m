classdef EF

    methods (Static)

        function is_eigen = EigenfunctionCheck(level, eigenvalue, eigenfunction)
            % EIGENFUNCTION is a Dirichlet eigenfunction of EIGENVALUE
            % Returns a boolean
            gamma = GraphApprox(level);
            npts = length(gamma.vertices);
            lapl = GLaplacian(gamma, eigenfunction);
            laplpar = lapl(4:npts);
            is_eigen = norm(laplpar+eigenvalue*eigenfunction(4:npts)) < 0.001;  % Degree of accuracy
        end

        function ok = EigenfunctionsCheck(efs, evs)
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

        function EFandEV = EigenfunctionsAll(level)
            % EFandEV finds all eigenvalues and eigenfunctions of the
            % /level-approximation of SG
            npts = length(GraphApprox(level).vertices);
            syms x [1 npts];
            syms ev;
            laplvar = GLaplacian(GraphApprox(level), x);
            laplparvar = laplvar(4:npts);
            EFandEV = solve([laplparvar == -ev*x(4:npts), x(1:3) == [0 0 0], norm(x)^2 ~= 0], [x, ev], 'ReturnCondition', true);
        end
        
        function new = Symmetric(old, level)
            if level == 2
                new = old - old([2,1,3,4,6,5,10,12,11,7,9,8,13,15,14],1);
            else
                new = [];
            end
        end
        
        function new = SymmetricAll(old, level, varargin)
            new{1} = old{1};
            if nargin == 3
                list = varargin{1};
            else
                [~, num] = size(old{1});
                list = 1:num;
            end
            for i = list
                new{1}(1:end,i) = EF.Symmetric(old{1}(1:end,i), level);
            end
            new{2} = old{2};
            new{3} = old{3};
        end
    end
end

