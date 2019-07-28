classdef Orthonormal
    
    methods(Static)
        function check = Test(spectrum, details)
            %TEST(s,d) returns 1 if the set of functions s (formatted as a result
            %of functions below) is an orthogonal spectrum with the given
            %eigenvalues
            mat = double(spectrum{1}(spectrum{3}, 1:end)'*spectrum{1}(spectrum{3}, 1:end));
            %orth = norm(mat - eye(length(spectrum{3}))) < 0.1;
            orth = norm(mat - diag(diag(mat))) < 0.1 && rank(mat) == length(mat);
            eigen = EF.EigenfunctionsCheck(spectrum{1}', spectrum{2}');
            if details
                disp(["orthogonal:",orth])
                disp(mat)
                disp(["eigenfunctions:", eigen]) 
            end
            check = orth && eigen;
        end
        
        %the spectrum searching functions below have the same
        %specifications as FindSpectrum
        %choose set of eigenvalues randomly, solve
        function result = FindSpectrum(basis, tile, limit, numerical, varargin)
            %FINDSPECTRUM(b, t, l, n) returns an orthogonal spectrum of the
            %subset (i.e., set of points) t of the  SG approximations given by b.
            %l specifies the number of attempts to solve for the spectrum.
            %If n is set to true, equations are solved numerically.
            %The output is a cell array where the first entry conatins the
            %matrix of eigenfunctions, the second enrty is the row vector
            %of eigenvalues, the third entry is t.
            w = waitbar(0,mat2str(tile));
            if numerical
                digits(4)
            end
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
                waitbar(n/limit,w);
                chosen = randperm(dim, sdim);
                subnew = new(tile, chosen);
                conditions = [];
                if nargin == 5
                    syms u [dim+3 1];
                    for k = 1:sdim
                        conditions = [conditions,subs(varargin{1}, u, new(1:end,chosen(k)))];
                    end
                end
                if numerical
                    eqa = vpa(subnew'*subnew == eye(sdim));
                    sol = vpasolve([reshape(eqa, [1, sdim^2]), conditions]);
                else
                    eqa = vpa(simplify(subnew'*subnew == eye(sdim)));
                    sol = solve([reshape(eqa, [1, sdim^2]),conditions]);
                end
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
                    result = {orth, eigenvalues, tile};
                    digits(32)
                    close(w);
                    return
                end
            end
            digits(32)
            result = {};
            close(w);
        end
        
        %same but with parallel
        function result = FindSpectrumPar(basis, tile, limit,numerical)
            if numerical
                digits(4)
            else
            digits(32)
            end
            results = [];
            dim = length(basis.eigenvalues);
            sdim = length(tile);
            indep = ones([dim, dim]);
            for i = 1:dim
                for j = 1 : i
                    if basis.eigenvalues(i) ~= basis.eigenvalues(j)
                        indep(i,j) = 0;
                        indep(j,i) = 0;
                    end
                end
            end
            syms P [dim dim];
            assume(P, 'real')
            new = (basis.x)*(P.*indep);
            tilenew = new(tile, 1:end);
            eqa = vpa(simplify(tilenew'*tilenew == eye(dim)));
            parfor n = 1:limit
                chosen = randperm(dim, sdim);
                subeqa = eqa(chosen, chosen);
                if numerical
                    sol = vpasolve(reshape(subeqa, [1, sdim^2]));
                else
                    sol = solve(reshape(subeqa, [1, sdim^2]));
                end
                resP = zeros([dim dim]);
                for i = 1 : dim
                    for j = 1 : dim
                        if isfield(sol, strcat('P',int2str(i),'_',int2str(j))) && not(isempty(sol.(strcat('P',int2str(i),'_',int2str(j)))))
                            vary = symvar(sol.(strcat('P',int2str(i),'_',int2str(j)))(1,1))
                            l = length(vary)
                            resP(i,j) = subs(sol.(strcat('P',int2str(i),'_',int2str(j)))(1,1), vary, zeros([1 l]));%sol.(strcat('P',int2str(i),'_',int2str(j)))(1,1)
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
                    result = {orth, eigenvalues, tile};
                    results = [results,{result}];
                end
            end
            digits(32)
            if not(isempty(results))
                result = results(1);
                result = result{1};
            else
                result = {};
            end
        end
        
        %try adding vetors one by one choosing random eigenvalue each time
        function result = Spectrum(basis, tile, limit)
            [npts, ~] = size(basis.x);
            dim = length(tile);
            result{3} = tile;
            for i = 1:limit
                evunique = unique(basis.eigenvalues);
                result{1} = zeros([npts, dim]);
                result{2} = zeros([1 dim]);
                curr = 0;
                while not(isempty(evunique))
                    idx = randi(length(evunique));
                    ev = evunique(idx);
                    subbasis = basis.x(1:end, basis.eigenvalues == ev);
                    tilebasis = basis.x(tile, basis.eigenvalues == ev);
                    [~, num] = size(tilebasis);
                    syms v [num 1]
                    new = tilebasis*v;
                    if curr > 0
                        sol = solve([simplify(result{1}(tile,1:curr)'* new == zeros([curr,1])), simplify(norm(new)^2 == 1)], v);
                    else
                        sol = vpasolve(simplify(norm(new)^2 == 1), v);
                    end
                    if isempty(sol)
                        evunique(idx) = [];
                    else
                        resV = zeros([num 1]);
                        if num == 1
                            resV = sol;
                        elseif isempty(sol.v1)
                            evunique(idx) = [];
                            continue
                        else
                            for j = 1 : num
                                if isfield(sol, strcat('v',int2str(j)))
                                    resV(j,1) = sol.(strcat('v',int2str(j)))(1,1);
                                end
                            end
                        end
                        result{1}(1:end, curr+1) = subbasis*resV;
                        result{2}(curr+1) = ev;
                        curr = curr + 1;
                    end
                    if result{2}(dim) ~= 0
                        return
                    end
                end
            end
        end
        
        %same but using null instead of solve
        function result = Spectrum2(basis, tile, limit, varargin)
            [npts, ~] = size(basis.x);
            dim = length(tile);
            result{3} = tile;
            if nargin > 3 
                evs = varargin{1};
            end
            for i = 1:limit
                if nargin > 3
                    evs = evs(randperm(length(evs)));
                end
                evunique = unique(basis.eigenvalues);
                result{1} = zeros([npts, dim]);
                result{2} = zeros([1 dim]);
                curr = 0;
                counter = 0;
                while not(isempty(evunique) || (nargin > 3 && counter == length(evs)))
                    if nargin > 3
                        ev = evs(curr+1);
                        evactual = basis.eigenvalues(abs(basis.eigenvalues - ev)<0.01);
                        ev = evactual(1);
                    else
                        idx = randi(length(evunique));
                        ev = evunique(idx);
                    end
                    subbasis = basis.x(1:end, basis.eigenvalues == ev);
                    tilebasis = basis.x(tile, basis.eigenvalues == ev);
                    [~, num] = size(tilebasis);
                    syms v [num 1]
                    new = tilebasis*v;
                    if curr > 0
                        [A,~] = equationsToMatrix(result{1}(tile,1:curr)'* new == zeros([curr,1]), v);
                        sol = null(A);
                        if isempty(sol)
                            evunique(evunique == ev) = [];
                        else
                            [~, nsp] = size(sol);
                            order = randperm(nsp);
                            for k = 1:nsp
                                vec = subbasis*sol(1:end,order(k));
                                if vec(tile,1) == zeros([dim 1])
                                    if k == nsp
                                        evunique(evunique == ev) = [];
                                    end
                                    continue
                                else
                                    result{1}(1:end, curr+1) = vec;
                                    result{2}(curr+1) = ev;
                                    curr = curr + 1;
                                    break
                                end
                            end
                        end
                    else
                        %{
                        sol = vpasolve(simplify(norm(new)^2 == 1), v);
                        resV = zeros([num 1]);
                        if num == 1
                            resV = sol;
                        elseif isempty(sol.v1)
                            evunique(evunique == ev) = [];
                            continue
                        else
                            for j = 1 : num
                                if isfield(sol, strcat('v',int2str(j)))
                                    resV(j,1) = sol.(strcat('v',int2str(j)))(1,1);
                                end
                            end
                        end
                        %}
                        resV = zeros([num 1]);
                        resV(1,1) = 1;
                        result{1}(1:end, curr+1) = subbasis*resV;
                        result{2}(curr+1) = ev;
                        curr = curr + 1;
                       
                    end
                    if result{2}(dim) ~= 0
                        return
                    end
                    counter = counter + 1;
                end
            end
        end
    end
end

