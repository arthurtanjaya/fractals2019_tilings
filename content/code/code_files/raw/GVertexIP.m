function inner_product = GVertexIP(u, v)
%GVERTEXIP Computes the vertex inner product of two functions on \Gamma_m
%   The answer is not renormalized
%   No error checking yet; we assume u, v are properly formatted

inner_product = 0;
[nrows, ncols] = size(u);
for i = 1:nrows
    for j = 1:ncols
        inner_product = inner_product + u(i, j) * v(i, j);
    end
end

end
