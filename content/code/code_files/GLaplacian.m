function lapl = GLaplacian(gamma, u)
%GLAPLACIAN Computes the energy of two function
% gamma is a GraphApprox
% u is a vector where the i'th entry corresponds to the value of the
% function on the i'th vertex

pts = gamma.vertices;
npts = length(pts);
lapl = sym(zeros(1, npts));
for i = 1:npts
    nbrs = pts(i).neighbors;
    for j = 1:npts
        if ismember(pts(j).address, nbrs, 'rows')
            lapl(i) = lapl(i) + u(j) - u(i);
        end
    end
end

end
