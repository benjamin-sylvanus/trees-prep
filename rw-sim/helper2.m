function [inds] = helper2 (bsize, ci)

    try
        inds = sub2ind(bsize, ci(:, 1), ci(:, 2), ci(:, 3));
    catch
    end

end
