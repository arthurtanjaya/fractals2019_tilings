classdef EF
    
    methods(Static)
        
        function eigenfunctionCheck = EigenfunctionCheck(level,eigenvalue,eigenfunction)
            %EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
            npts = length(GraphApprox(level).vertices);
            lapl = GLaplacian(GraphApprox(level),eigenfunction);
            laplpar = lapl(4:npts);
            eigenfunctionCheck = norm(laplpar+eigenvalue*eigenfunction(4:npts)) < 0.001;
        end
        
        function eigenfunctions = Eigenfunctions(level,eigenvalue)
            %EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
            npts = length(GraphApprox(level).vertices);
            syms x [1 npts];
            laplvar = GLaplacian(GraphApprox(level),x);
            laplparvar = laplvar(4:npts);
            eigenfunctions = solve([laplparvar == -sym(eigenvalue)*x(4:npts), x(1:3) == [0,0,0], norm(x)^2~= 0], x, 'ReturnCondition', true);
        end 
        
        function EFandEV = EigenfunctionsAll(level)
            %EFandEV finds all eigenvalues and eigenfunctions of the
            %/level-approximation of SG
            npts = length(GraphApprox(level).vertices);
            syms x [1 npts];
            syms ev;
            laplvar = GLaplacian(GraphApprox(level),x);
            laplparvar = laplvar(4:npts);
            EFandEV = solve([laplparvar == -ev*x(4:npts), x(1:3) == [0,0,0], norm(x)^2~=0], [x, ev], 'ReturnCondition', true);
        end
            
    end
end