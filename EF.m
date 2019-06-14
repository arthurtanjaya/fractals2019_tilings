classdef EF

methods(Static)     
function eigenfunction = Eigenfunction(level,eigenvalue)
%EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
npts = length(GraphApprox(level).vertices);
syms x [1 npts];
laplvar = GLaplacian(GraphApprox(level),x);
laplparvar = laplvar(4:npts);
eigenfunction = solve([laplparvar == -eigenvalue*x(4:npts), x(1:3) == [0,0,0], norm(x)~= 0], x, 'ReturnCondition', true);
end

function EFandEV = EigenfunctionsAll(level)
%EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
npts = length(GraphApprox(level).vertices);
syms x [1 npts];
syms ev;
laplvar = GLaplacian(GraphApprox(level),x);
laplparvar = laplvar(4:npts);
EFandEV = solve([laplparvar == -ev*x(4:npts), x(1:3) == [0,0,0], norm(x)^2~=0], [x, ev], 'ReturnCondition', true);
end

function eigenfunctionCheck = EigenfunctionCheck(level,eigenvalue,eigenfunction)
%EIGENFUNCTION is a Dirichlet eigenfunction of /eigenvalue
npts = length(GraphApprox(level).vertices);
lapl = GLaplacian(GraphApprox(level),eigenfunction);
laplpar = lapl(4:npts);
eigenfunctionCheck = isequal(laplpar,-eigenvalue*eigenfunction(4:npts));
end

function orthArray = OrthogonalityCheck(functions)
orthArray = [];
[nfns,~] = size(functions);
for i = 1 : nfns
    for j = 1 : (i-1)
         orthArray = [orthArray, dot(functions(i,:),functions(j,:))];
    end
    orthArray = [orthArray, (norm(functions(i,:)))^2-1];
end
end

function orthB = FindOrth(efBasis)
par = length(efBasis.parameters);
dim = length(efBasis.ev);
syms coefs [dim par];
x = [];
qBasis = [];
i = 1;
while isfield(efBasis,strcat('x',int2str(i)))
    x = [x, getfield(efBasis,strcat('x',int2str(i)))];
    i = i + 1;
end
for j = 1 : dim
    qBasis = [qBasis; subs(x(j,:), efBasis.parameters, coefs(j,:))];
end
orthB = solve(EF.OrthogonalityCheck(qBasis) == zeros([1, (dim * (dim + 1)/2)]));
end

function orthB = FindOrthTwo(efBasis)
par = length(efBasis.parameters);
dim = length(efBasis.ev);

    function orth = F(coefs)
    x = [];
    qBasis = [];
    i = 1;
    while isfield(efBasis,strcat('x',int2str(i)))
        x = [x, getfield(efBasis,strcat('x',int2str(i)))];
        i = i + 1;
    end
    for j = 1 : dim
        qBasis = [qBasis; subs(x(j,:), efBasis.parameters, coefs(j,:))];
    end
    orth = double(EF.OrthogonalityCheck(qBasis));
    end
    

fun = @F;
orthB = fsolve(fun, ones([dim, par]));
end
end
end