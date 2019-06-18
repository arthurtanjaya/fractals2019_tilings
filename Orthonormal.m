classdef Orthonormal
    
    methods(Static)
        function orthArray = OrthogonalityCheck(functions)
            [nfns,~] = size(functions);
            for i = 1 : nfns
                for j = 1 : (i-1)
                    orthArray((i^2-i)/2+j) = dot(functions(i,:),functions(j,:));
                end
                orthArray((i^2+i)/2) = norm(functions(i,:))^2-1;
                endS
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
            transpose(Orthonormal.OrthogonalityCheck(qBasis) == zeros([1, (dim * (dim + 1)/2)]))
            orthB = solve(simplify(Orthonormal.OrthogonalityCheck(qBasis) == zeros([1, (dim * (dim + 1)/2)])), coefs);
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
        
        function orthB = FindOrth3(data)
            dim = length(data.eigenvalues);
            syms P [dim dim];
            indep = ones([dim, dim]);
            acc = [];
            for i = 1:dim
                if ismember(data.eigenvalues(i), acc)
                    for j = 1 : i
                        if data.eigenvalues(i) ~= data.eigenvalues(j)
                            indep(i,j) = 0;
                            indep(j,i) = 0;
                        end
                    end
                else
                    for j = 1 : dim
                        indep(j,i) = 0;
                    end
                    indep(i,i) = 1;
                    acc = [acc, data.eigenvalues(i)];
                end
            end
            indep(3,2) = 0;
            new = (data.x)*(P.*indep);
            eqa = Orthonormal.OrthogonalityCheck(new') == zeros([1, (dim ^2 + dim)/2]);
            reshape(eqa, [(dim^2+dim)/2, 1])
            sol = solve(eqa);
            resP = zeros([dim dim]);
            for i = 1 : dim
                for j = 1 : dim
                    if isfield(sol, strcat('P',int2str(i),'_',int2str(j)))
                        resP(i,j) = sol.(strcat('P',int2str(i),'_',int2str(j)));
                    else
                        resP(i, j) = 0;
                    end
                end
            end
            orthB = (data.x)*(resP.*indep);
        end
        
        function [orth, eigenvalues] = FindSpectrum(basis, tile, limit)
            digits(5)
            dim = length(basis.eigenvalues);
            sdim = length(tile);
            syms P [dim dim];
            assume(P, 'real')
            indep = ones([dim, dim]);
            for i = 1:dim
                for j = 1 : i
                    if basis.eigenvalues(i) ~= basis.eigenvalues(j)
                        indep(i,j) = 0;
                        indep(j,i) = 0;
                    end
                end
            end
            new = (basis.x)*(P.*indep);
            for n = 1:limit
                chosen = randperm(dim, sdim);
                subnew = new(tile, chosen);
                eqa = simplify(vpa(subnew'*subnew) == eye(sdim));
                sol = solve(reshape(eqa, [1, sdim^2]));
                resP = zeros([dim dim]);
                for i = 1 : dim
                    for j = 1 : dim
                        if isfield(sol, strcat('P',int2str(i),'_',int2str(j))) && not(isempty(sol.(strcat('P',int2str(i),'_',int2str(j)))))
                            resP(i,j) = subs(sol.(strcat('P',int2str(i),'_',int2str(j)))(1,1), P, zeros([dim dim]));
                        else
                            resP(i, j) = 0;
                        end
                    end
                end
                if resP == zeros([dim dim])
                    continue
                else
                    orthB = (basis.x)*(resP.*indep);
                    orth = orthB(1:end,chosen);
                    eigenvalues = basis.eigenvalues(chosen);
                    if logical(norm(double(orth(tile, 1:end)'*orth(tile, 1:end)) - eye(sdim)) > 0.01)
                        double(orth(tile, 1:end)'*orth(tile, 1:end))
                        "not o.n."
                    elseif not(EF.EigenfunctionsCheck(orth', eigenvalues'))
                        "not basis"
                    else
                        "OK"
                    end
                    return
                end
            end
            orth = [];
            eigenvalues = [];
        end
    end
end

