function dists = setdists(tree)
%SETDISTS Summary of this function goes here
%   Detailed explanation goes here
    dists = zeros(height(tree), 1);
    for i = 2:height(tree)
        dists(i, 1) = calcdists(tree, i);
    end
end

