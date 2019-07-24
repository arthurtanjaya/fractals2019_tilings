classdef SpectralDecimation
    
    methods(Static)
        function evlist = Eigenvalues(level)
            %EIGENVLAUES(level) gives a list of eigenvalues
            if level == 1
                evlist = [2, 5];
            else
                evlist = [5, 6];
                prev = SpectralDecimation.Eigenvalues(level - 1);
                for i = 1 : length(prev)
                    ev = prev(i);
                    evlist = [evlist, (5+(25-4*ev)^(1/2))/2, (5-(25-4*ev)^(1/2))/2];
                end
            end
        end
        
        function new = Initial(level)
            %INITIAL(level) gives the list of initial eigenfunctions and
            %the corresponding eigenvalues (which are 2, 5, and 6).
            %Eigenfunctions are represented as column vectors, eigenvalues
            %are stored in a row vector.
            if level == 1
                twos = EF.Eigenfunctions(level, 2);
                dim1 = length(twos.parameters);
                for i = 1 : (3^(level+1) + 3)/ 2
                    for j = 1 : dim1
                        new.x(i,j) = subs(twos.(strcat('x',int2str(i))), twos.parameters, [zeros([1,j-1]), 1, zeros([1, dim1-j])]);
                    end
                end
                fives = EF.Eigenfunctions(level, 5);
                dim2 = length(fives.parameters);
                for i = 1 : (3^(level+1) + 3)/ 2
                    for j = 1 : dim2
                        new.x(i,dim1 + j) = subs(fives.(strcat('x',int2str(i))), fives.parameters, [zeros([1,j-1]), 1, zeros([1, dim2-j])]);
                    end
                end
                new.eigenvalues = [2*ones([1,dim1]), 5*ones([1,dim2])];
                
            else
                fives = EF.Eigenfunctions(level, 5);
                dim1 = length(fives.parameters);
                for i = 1 : (3^(level+1) + 3)/ 2
                    for j = 1 : dim1
                        new.x(i,j) = subs(fives.(strcat('x',int2str(i))), fives.parameters, [zeros([1,j-1]), 1, zeros([1, dim1-j])]);
                    end
                end
                sixes = EF.Eigenfunctions(level, 6);
                dim2 = length(sixes.parameters);
                for i = 1 : (3^(level+1) + 3)/ 2
                    for j = 1 : dim2
                        new.x(i, dim1 + j) = subs(sixes.(strcat('x',int2str(i))), sixes.parameters, [zeros([1,j-1]), 1, zeros([1, dim2-j])]);
                    end
                end
                new.eigenvalues = [5*ones([1,dim1]), 6*ones([1,dim2,1])];
            end
        end
        
        function expand = Continued(level, eigenvalue, old)
            %CONTINUED(l, e, old) gives the list of continued
            %eigenfunctions with eigenvalue e 
            %that result from eigenfunction old on SG
            %approximation of level l-1.
            expand = zeros([1, (3^(level+1)+3)/2]);
            SG = GraphApprox(level).vertices;
            SGprev = GraphApprox(level - 1).vertices;
            for j = 1 : (3^(level+1) + 3)/2
                if SG(j).address(1) == SG(j).address(end)
                    oldaddress = SG(j).address(1:end-1);
                    expand(1,j) = old(GraphApprox.lookup(SGprev,level-1,oldaddress));
                else
                    nbrs = SG(j).neighbors;
                    fx1x2 = [];
                    y1y2 = [];
                    for l = 1 : 4
                        nbr = nbrs(l,:);
                        if nbr(1) == nbr(end)
                            fx1x2 = [fx1x2, old(GraphApprox.lookup(SGprev,level-1,nbr(1:end-1)))];
                        else
                            y1y2 = [y1y2, GraphApprox.lookup(SG,level,nbr)];
                        end
                    end
                    x3b = intersect(SG(y1y2(1)).neighbors, SG(y1y2(2)).neighbors, 'rows');
                    for m = 1 : 2
                        t = x3b(m,:);
                        if t(1) == t(end)
                            fx3 = old(GraphApprox.lookup(SGprev,level-1,t(1:end-1)));
                        end
                    end
                    fy0 = ((4 - eigenvalue) * (fx1x2(1) + fx1x2(2)) + 2 * fx3)/((2 - eigenvalue)*(5 - eigenvalue));
                    expand(1,j) = fy0;
                end
            end
        end
        
        function result = Expand(level, old, varargin)
            result{1} = [];
            result{2} = [];
            result{3} = old{3};
            if nargin == 2 
            for i = 1 : length(old{2})
                ev = old{2}(i);
                evplus = (5 + (25 - 4 * ev)^(1/2)) / 2;
                evminus = (5 - (25 - 4 * ev)^(1/2)) / 2;
                newplus = SpectralDecimation.Continued(level, evplus, old{1}(1:end,i)');
                newminus = SpectralDecimation.Continued(level, evminus, old{1}(1:end,i)');
                if evplus ~= 2 && evplus ~= 5
                    result{2} = [result{2}, evplus];
                    result{1} = [result{1}, newplus'];
                end
                if evminus ~= 2 && evminus ~= 5
                    result{2} = [result{2}, evminus];
                    result{1} = [result{1}, newminus'];
                end
            end
            else
                seq = varargin{1};
                for i = 1 : length(seq)
                    ev = old{2}(i);
                    if seq(i) >= 0
                        evnew = (5 + (25 - 4 * ev)^(1/2)) / 2;
                    else
                        evnew = (5 - (25 - 4 * ev)^(1/2)) / 2;
                    end
                    new = SpectralDecimation.Continued(level, evnew, old{1}(1:end,i)');
                    if evnew ~= 2 && evnew ~= 5
                        result{2} = [result{2}, evnew];
                        result{1} = [result{1}, new'];
                    end
                end
            end
        end
        
        function result = Eigenfunctions(level)
            result = SpectralDecimation.Initial(level);
            if level > 1
                old = SpectralDecimation.Eigenfunctions(level - 1);
                for i = 1 : length(old.eigenvalues)
                    ev = old.eigenvalues(i);
                    evplus = (5 + (25 - 4 * ev)^(1/2)) / 2;
                    evminus = (5 - (25 - 4 * ev)^(1/2)) / 2;
                    newplus = SpectralDecimation.Continued(level, evplus, old.x(1:end,i)');
                    newminus = SpectralDecimation.Continued(level, evminus, old.x(1:end,i)');
                    if evplus ~= 2 && evplus ~= 5
                        result.eigenvalues = [result.eigenvalues, evplus];
                        result.x = [result.x, newplus'];
                    end
                    if evminus ~= 2 && evminus ~= 5
                        result.eigenvalues = [result.eigenvalues, evminus];
                        result.x = [result.x, newminus'];
                    end
                    
                end
            end
        end
    end
end