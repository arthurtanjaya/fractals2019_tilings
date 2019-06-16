classdef SpectralDecimation
    
    methods(Static)
        function evlist = Eigenvalues(level)
            %EIGENVLAUES(level) gives a list of eigenvalues on level /level
            %bla
            if level == 1
                evlist = [2, 5];
            else
                evlist = [5, 6];
                prev = EF.Eigenvalues(level - 1);
                for i = 1 : length(prev)
                    ev = prev(i);
                    evlist = [evlist, (5 + (25 - 4 * ev)^(1/2)) / 2, (5 - (25 - 4 * ev)^(1/2)) / 2];
                end
            end
        end
        
        function new = Initial(level)
            if level == 1
                for i = 1 : (3^(level+1) + 3)/ 2
                    new.(strcat('x',int2str(i))) = [EF.Eigenfunctions(level, 2).(strcat('x',int2str(i))); EF.Eigenfunctions(level, 5).(strcat('x',int2str(i)))];
                end
                
            else
                for i = 1 : (3^(level+1) + 3)/ 2
                    new.(strcat('x',int2str(i))) = [EF.Eigenfunctions(level, 5).(strcat('x',int2str(i))); EF.Eigenfunctions(level, 6).(strcat('x',int2str(i)))];
                end
            end
        end
        
        function expand = Continued(level, eigenvalue, old)
            expand = zeros([1, (3^(level+1)+3)/2]);
            SG = GraphApprox(level).vertices;
            SGprev = GraphApprox(level - 1).vertices;
            function idx = lookup(on,address)
                if on == 0
                    for i = 1: (3^level + 3)/2
                        if (all(SGprev(i).address == address))
                            idx = i;
                        end
                    end
                else
                    for i = 1: (3^(level+1) + 3)/2
                        if (all(SG(i).address == address))
                            idx = i;
                        end
                    end
                end
            end
            for j = 1 : (3^(level+1) + 3)/2
                if SG(j).address(1) == SG(j).address(end)
                    oldaddress = SG(j).address(1:end-1);
                    expand(1,j) = old(lookup(0,oldaddress));
                else
                    nbrs = SG(j).neighbors;
                    fx1x2 = [];
                    y1y2 = [];
                    for l = 1 : 4
                        nbr = nbrs(l,:);
                        if nbr(1) == nbr(end)
                            fx1x2 = [fx1x2, old(lookup(0,nbr(1:end-1)))];
                        else
                            y1y2 = [y1y2, lookup(1,nbr)];
                        end
                    end
                    x3b = intersect(SG(y1y2(1)).neighbors, SG(y1y2(2)).neighbors, 'rows');
                    for m = 1 : 2
                        t = x3b(m,:);
                        if t(1) == t(end)
                            fx3 = old(lookup(0,t(1:end-1)));
                        end
                    end
                    fy0 = ((4 - eigenvalue) * (fx1x2(1) + fx1x2(2)) + 2 * fx3)/((2 - eigenvalue)*(5 - eigenvalue));
                    expand(1,j) = fy0;
                end
            end
        end
    end
end