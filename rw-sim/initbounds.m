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
    tran = ((dim - ranges(1, :)) ./ VSIZE) + 1;

    % Update Bounds
    boundSize = ceil((ranges(2, :) - ranges(1, :)) ./ VSIZE) + 2;

    % Radii Scale
    R = swc.Radii ./ (VSIZE);

    % Set Tree XYZ to Translated
    % min(swc.tree.Z-swc.tree.Radii);
    swc.X = tran(:, 1); swc.Y = tran(:, 2); swc.Z = tran(:, 3);
    swc.Radii = R;
    % min(swc.tree.Z-swc.tree.Radii)

    [b, pairs] = calcBounds(swc, VSIZE);
    addpath(genpath(['/Users/benjaminsylvanus/Documents/GitHub/' ...
                     'SparseMatrixGenerator/ndSparse_G4_2021_03_16']));
    b = b.a;

    % Input [r1;r2] [x1 y1 z1; x2 y2 z2]
    % Output [min(ri-ci) max(ri-ci)]
    function [b, pairs] = calcBounds(swc, vs)
        b = struct();
        t = swc;
        mr = vs;
        pairs = [];

        for i = 1:numel(t.X)
            pair = swc(swc{i, "Parent"}, :);
            pairs(i, 1) = t{i, "NodeId"};
            pairs(i, 2) = pair{1, "NodeId"};

            if (pairs(i, 2) == 0)
                pairs(i, 2) = pairs(i, 1);
            end

            b.a(i, :) = pairBounds(t(i, :), pair, vs);
        end

        function b = pairBounds(p1, p2, m)
            % ri: pair radii
            % pxi: pair coords

            ri = [p1.Radii; p2.Radii];
            pxi = [p1{1, 2:4}; p2{1, 2:4}];
            b = {[round(min(pxi - ri)); ceil(max(pxi + ri))]};
            ba = [min(pi - (ri + m)), max(pi + (ri + m))];
        end

    end

end