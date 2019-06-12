function energy = GEnergyIP(gamma, u, v)
%GENERGYIP Computes the energy of two function
%gamma is a GRAPHAPPROX
% u,v are vectors where the i'th entry corresponds to the value of the
% function on the i'th vertex
pts = gamma.vertices;
energy = 0;
npts = length(pts);
for i = 1:npts
    for j = 1:(i-1)
        if ismember(get_address(pts(j)),pts(i).get_neighbors)
            energy = energy + (u(i)-u(j))*(v(i)-v(j));
        end
    end
end

end
