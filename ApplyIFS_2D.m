function coord = ApplyIFS_2D(coord, word)
%APPLYIFS_2D Applies F_{word} to coord, in 2D Cartesian coordinates
%   Note that word should not include the first bit of a vertex address

% Define parameters
N = 3;                % Number of contractions
r = 1/2;              % Contraction ratio
q = [[0.5, 3^0.5/2];
     [0, 0];
     [1, 0]];         % Fixed points of respective contractions
F = cell(1, 3);       % Contractions
for i = 1:N
    F{i} = @(x) r*(x-q(i)) + q(i);
end

for i = word
    coord = F{i+1}(coord);  % This annoying 0/1-indexing issue...
end
