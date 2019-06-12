function product = CartesianProd(vector, n)
%CARTESIANPROD Returns the n-fold Cartesian product of vector
% Taken from https://www.mathworks.com/matlabcentral/answers/282777-n-fold-cartesian-product

product = repmat({vector}, 1, n);
[product{1:n}] = ndgrid(product{:});
product = reshape(cat(n, product{:}), [], n);

end
