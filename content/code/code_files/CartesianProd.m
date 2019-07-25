function product = CartesianProd(vector, n)
%CARTESIANPROD Returns the n-fold Cartesian product of vector
%   If n == 0 then this function returns {[]}
% Taken from https://www.mathworks.com/matlabcentral/answers/282777-n-fold-cartesian-product

if n == 0
    product = {[]};
else
    product = repmat({vector}, 1, n);
    [product{1:n}] = ndgrid(product{:});
    product = reshape(cat(n, product{:}), [], n);
end

end
