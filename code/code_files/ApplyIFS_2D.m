function coord = ApplyIFS_2D(coord, word)
%APPLYIFS_2D applies F_{word} to coord, in 2D Cartesian coordinates
%   Note that word should not include the first bit of a vertex address

% Define parameters
r = 1/2;              % Contraction ratio
q = [[0.5, 3^0.5/2];
     [0, 0];
     [1, 0]];         % Fixed points of respective contractions
F = @(x, i) r*(x-q(i+1, :)) + q(i+1, :);  % Contractions

for w = word(end:-1:1)
    coord = F(coord, w);
end
