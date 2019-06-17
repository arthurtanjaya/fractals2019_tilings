function lapl = GLaplacian(gamma, u)
%GENERGYIP Computes the energy of two function
%gamma is a GRAPHAPPROX
% u,v are vectors where the i'th entry corresponds to the value of the
% function on the i'th vertex
pts = gamma.vertices;
npts = length(pts);
lapl = [];
for i = 1:npts
    laplP = 0;
    nbrs = pts(i).neighbors;
    for j = 1:npts
        if ismember(pts(j).address,nbrs,'rows')
            laplP = laplP + (u(j)-u(i));
        end
    end
    lapl = [lapl,laplP];
end
end
