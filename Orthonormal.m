classdef Orthonormal
    
    methods(Static)
        function orthArray = OrthogonalityCheck(functions)
            [nfns,~] = size(functions);
            for i = 1 : nfns
                for j = 1 : (i-1)
                    orthArray((i^2-i)/2+j) = dot(functions(i,:),functions(j,:));
                end
                orthArray((i^2+i)/2) = norm(functions(i,:))^2-1;
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
            sol = solve(eqa)
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
            resP.*indep
            orthB = (data.x)*(resP.*indep);
        end
    end
end

