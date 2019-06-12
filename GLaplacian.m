function lapl = GLaplacian(gamma, u)
%GENERGYIP Computes the energy of two function
%gamma is a GRAPHAPPROX
% u,v are vectors where the i'th entry corresponds to the value of the
% function on the i'th vertex
pts = gamma.vertices;
npts = length(pts);
lapl = zeros([1,npts]);
for i = 1:npts
    laplP = 0;
    for j = 1:npts
        if ismember(get_address(pts(j)),pts(i).get_neighbors)
            laplP = laplP + (u(j)-u(i));
        end
        lapl(i) = laplP;
    end
end
end
