function temp = helper(temp, inds, i)

    k = repelem(i, size(temp, 1));
    temp(:, end + 1) = {i};
    temp
end
