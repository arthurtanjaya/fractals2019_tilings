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

function orthCheck = OrthogonalityCheck(level,functions)
orthCheck = [];
    V = GraphApprox(level);
[nfns,~] = size(functions);
for i = 1 : nfns
    for j = 1 : (i-1)
         orthCheck = [orthCheck,GVertexIP(V, functions(i,:), functions(j,:))];
    end
end
end


function orthB = FindOrth(level,efBasis)
dim = length(efBasis.parameters);
syms coefs [dim dim];
qBasis = [];
x = [efBasis.x1,efBasis.x2,efBasis.x3,efBasis.x4,efBasis.x5,efBasis.x6];
for i = 1 : dim
    qBasis = [qBasis; subs(x, efBasis.parameters, coefs(i,:))];
end
qBasis
orthB = solve([EF.OrthogonalityCheck(level,qBasis) == zeros, rank(coefs)== dim], coefs);
end
end
end