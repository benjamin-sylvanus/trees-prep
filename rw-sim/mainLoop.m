function [A] = mainLoop(tree, A, pairBounds, pairs)
    swc = tree{:, :};

    tic;

    for i = 1:size(pairs, 1)

        Bound = pairBounds{i, 1};

        Sx = Bound(1, 1); Sy = Bound(1, 2); Sz = Bound(1, 3);

        Nx = Bound(2, 1); Ny = Bound(2, 2); Nz = Bound(2, 3);

        [y0, x0, z0] = meshgrid(Sy:Ny, Sx:Nx, Sz:Nz);

        p1 = swc(pairs(i, 1), 2:5);

        pair = swc(pairs(i, 2), :);

        p2 = pair(1, 2:5);

        x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);

        x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

        pos = swc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, Nx, Ny, Nz);

        pos_fill = A(Sx:Nx, Sy:Ny, Sz:Nz);

        if sum(pos, "all") > 0
            is = isosurface(y0, x0, z0, pos, 0);
            sum(pos, "all");
            p = patch('Faces', is.faces, 'Vertices', is.vertices);
            p.FaceColor = "green";
            p.FaceAlpha = 0.1;
            p.EdgeColor = "none";
        end

        A(Sx:Nx, Sy:Ny, Sz:Nz) = pos_fill | pos;

    end

    A(Sx:Nx, Sy:Ny, Sz:Nz) = pos_fill;
end
