function [LUT, B, pairs, boundSize, swc, memoized_distance,A,vsize,ranges] = prepSim(voxelscale)
    % addpath(genpath("../../treestoolbox"));
    addpath(genpath("./"));

    % addpath("/Users/bensylvanus/Library/Application Support/MathWorks/MATLAB Add-Ons/" + ...
    % "Collections/random unit vector generator");

    tree = read_swc('democell10ele.swc');

    tree{:,2:5} = tree{:,2:5}.*1e-3;
    
    % clc;
    tic;
    % clearvars;
    % load("tree.mat", "tree");
    dists = zeros(height(tree), 1);

    for i = 2:height(tree)
        dists(i, 1) = calcdists(tree, i);
    end

    toc;

    tic;
    [b, swc, boundSize, pairs, vsize,ranges] = initbounds(tree, dists, voxelscale);

    toc;
    tic;
    [A, indicies, t2, LUT] = generateLUT(boundSize, b);
    toc;

    A = A(~cellfun('isempty', A));

    % get pairs from A
    B = cell(size(A, 1), 2);

    % extract node id [child, parent] -> linear index in B
    % B(i) = {[children],[parents]}
    for i = 1:length(B)
        pairlist = A{i};
        ps = pairs(pairlist, :);

        % pairlist => [child, parent]
        children = ps(:, 1);
        parents = ps(:, 2);

        B{i, 1} = children;
        B{i, 2} = parents;
    end

    % extract node id [child, parent] -> linear index in B
    % B(i) = {[children],[parents]}
    memoized_distance = zeros(length(pairs), 3);

    for i = 1:length(pairs)
        % extract child and parent ids
        baseid = pairs(i, 1); targetid = pairs(i, 2);

        p1 = swc{baseid, 2:5};
        p2 = swc{targetid, 2:5};

        x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
        x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

        memoized_distance(i, 3) = (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2;
        memoized_distance(i, 2) = targetid;
        memoized_distance(i, 1) = baseid;
    end

end
