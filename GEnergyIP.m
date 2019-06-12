function energy = GEnergyIP(u, v)
%GENERGYIP Computes the energy inner product of two functions on \Gamma_m
%   The answer is not renormalized
%   No error checking yet; we assume u, v are properly formatted
%   This is currently WRONG

energy = 0;
[nrows, ncols] = size(u);
for i = 1:nrows
    for j = 1:ncols
        energy = energy + u(i, j) * v(i, j);
    end
end

end
