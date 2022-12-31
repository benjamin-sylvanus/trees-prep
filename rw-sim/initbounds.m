function [b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, voxelScale)

    swc = tree;
    swc(1, "Parent") = {1};
    X = swc.X + 2; Y = swc.Y + 2; Z = swc.Z + 2;
    dim = [X Y Z];

    % Threshold Radii
    [~, ~, bin] = histcounts(swc.Radii, 5);
    avg = mean(swc.Radii(bin == 1));
    swc.Radii(swc.Radii < avg) = avg;
    r = swc.Radii;
    ranges = [min(dim - r); max(dim + r)]; % (min - ri) (max + ri)

    % Define Voxel Size
    VSIZE = voxelScale * min(dists(2:end));

    % Translate X,Y,Z
    tran = ((dim - ranges(1, :)) ./ VSIZE) + 3;

    % Update Bounds
    boundSize = ceil((ranges(2, :) - ranges(1, :)) ./ VSIZE) + 6;

    % Radii Scale
    R = swc.Radii ./ (VSIZE);

    % Set Tree XYZ to Translated
    % min(swc.tree.Z-swc.tree.Radii);
    swc.X = tran(:, 1); swc.Y = tran(:, 2); swc.Z = tran(:, 3);
    swc.Radii = R;
    % min(swc.tree.Z-swc.tree.Radii)

    [b, pairs] = calcBounds(swc{:, :});
    addpath(genpath(['/Users/benjaminsylvanus/Documents/GitHub/' ...
                     'SparseMatrixGenerator/ndSparse_G4_2021_03_16']));
    b = b.a;

    % Input [r1;r2] [x1 y1 z1; x2 y2 z2]
    % Output [min(ri-ci) max(ri-ci)]
    function [b, pairs] = calcBounds(swc)
        b = struct();
        pairs = zeros(size(swc, 1), 2);

        for i = 1:size(swc, 1)
            % extract parent id
            pair = swc(swc(i, 6), :);

            pairs(i, 1) = swc(i, 1);
            pairs(i, 2) = pair(1, 1);

            if (pairs(i, 2) == 0)
                pairs(i, 2) = pairs(i, 1);
            end

            b.a(i, :) = pairBounds(swc(i, :), pair);
        end

        function b = pairBounds(p1, p2)
            % ri: pair radii
            % pxi: pair coords

            ri = [p1(1, 5); p2(1, 5)];
            pxi = [p1(1, 2:4) - 1; p2(1, 2:4) + 1];
            b = {[round(min(pxi - ri)); ceil(max(pxi + ri))]};
        end

    end

end
