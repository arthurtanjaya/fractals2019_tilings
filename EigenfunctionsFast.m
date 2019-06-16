EigenfunctionsFast(level)
%DEPRICATED
%EFandEV finds all  eigenfunctions of the
%/level-approximation of SG
evs = EF.Eigenvalues(level);
npts = length(GraphApprox(level).vertices);
for i = 1 : npts
    new.(strcat('x',int2str(i))) = sym.empty(length(evs),0);
end
new.ev = evs;
syms x [1 npts];
for i = 1 : length(evs)
    ev = evs(i);
    tempEFs = EF.Eigenfunctions(level, ev);
    for l = 1 : npts
        res = tempEFs.(strcat('x',int2str(l)));
        update = new.(strcat('x',int2str(l)));
        update(i,1) = res(1,1);
        new.(strcat('x',int2str(l))) = update;
    end
end

