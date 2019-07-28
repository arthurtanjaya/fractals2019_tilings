function energy = GEnergyIP(gamma, u, v)
%GENERGYIP Computes the energy of two functions
% gamma is a GraphApprox
% u, v are vectors where the i'th entry corresponds to the value of the
% function on the i'th vertex

pts = gamma.vertices;
npts = length(pts);
energy = 0;
for i = 1:npts
    nbrs = pts(i).neighbors;
    for j = 1:i-1
        if ismember(pts(j).address, nbrs, 'rows') || ismember(pts(j).get_secondary(), nbrs, 'rows')
            energy = energy + (u(i)-u(j))*(v(i)-v(j));
        end
    end
end

end
