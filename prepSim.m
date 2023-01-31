function ...
        [ ...
         LUT, ...
         B, ...
         pairs, ...
         boundSize, ...
         swc, ...
         memoized_distance, ...
         A, ...
         vsize, ...
         ranges, ...
         b ...
     ] ...
        = prepSim(voxelscale, filename, scale)
    addpath(genpath("./"));
    tree = read_swc(filename);

    % neuroglancer scale [8nm 8nm 33nm]
    % size is [515892;356400;5293]

    % nm * 1e-3 = um
    % [x/8 y/8 z/33]*1e-3
    % * any scale nm -> scalenm = [sx sy sz]
    % * scalenm -> passed params
    %     tree{:,2:4} = tree{:,2:4}; % units [xyz]: nm;
    %     tree{:,2:5} = tree{:,2:5}.*1e-3; % units [xyz&r]: um;

    % scale should convert swc coords -> um;
    tree{:, 2:5} = tree{:, 2:5} .* scale; % units [xyz&r]: um;

    dists = setdists(tree);

    tic;
    [b, swc, boundSize, pairs, vsize, ranges] = initbounds(tree, dists, voxelscale);
    toc;

    tic;
    [A, LUT] = generateLUT(boundSize, b);
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

%     mainLoop(swc{:,:},LUT,b,pairs);
end
