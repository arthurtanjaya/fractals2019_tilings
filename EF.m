classdef EF

methods(Static)     
function eigenfunction = Eigenfunction(level,eigenvalue)
%EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
npts = length(GraphApprox(level).vertices);
syms x [1 npts];
laplvar = GLaplacian(GraphApprox(level),x);
laplparvar = laplvar(4:npts);
eigenfunction = solve([laplparvar == -eigenvalue*x(4:npts), x(1:3) == [0,0,0]], x, 'ReturnCondition', true);
end

function eigenfunctionCheck = EigenfunctionCheck(level,eigenvalue,eigenfunction)
%EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
npts = length(GraphApprox(level).vertices);
lapl = GLaplacian(GraphApprox(level),eigenfunction);
laplpar = lapl(4:npts);
eigenfunctionCheck = isequal(laplpar,-eigenvalue*eigenfunction(4:npts));
end
end
end