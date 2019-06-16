classdef Orthonormal
    
    methods(Static)
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

