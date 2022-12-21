function [A, poses] = mainLoop(tree, A, pairBounds, pairs)
    A = A;
    swc = tree;

    ecount = 0; prevElapse = 0; tic;

    for i = 1:size(pairs, 1)
        Bound = pairBounds{i, 1};

        Sx = Bound(1, 1); Sy = Bound(1, 2); Sz = Bound(1, 3);

        Nx = Bound(2, 1) + 1; Ny = Bound(2, 2) + 1; Nz = Bound(2, 3) + 1;

        [y0, x0, z0] = meshgrid(Sy:Ny, Sx:Nx, Sz:Nz);
        % pos_fill = ndSparse.build(range(Bound)+2);

        p1 = swc{pairs(i, 1), 2:5};

        pair = swc(pairs(i, 2), :);

        p2 = pair{1, 2:5};

        x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);

        x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

        pos = swc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, Nx, Ny, Nz);

        pos_fill = A(Sx:Nx, Sy:Ny, Sz:Nz);

        Q = find(pos == 1);
        [ix, iy, iz] = ind2sub(size(pos), Q);
        mxyz = mean([iy, ix, iz]);

        poses(i, :) = mxyz + [Sy Sx Sz];

        hold on;

        try
            is = isosurface(pos, 0);
            sum(pos, "all");
            is.vertices = is.vertices + [Sy Sx Sz];
            p = patch('Faces', is.faces, 'Vertices', is.vertices);
            p.FaceColor = 'green';
            p.FaceAlpha = 0.3;
            p.EdgeColor = 'none';
            view(3);
            axis tight;

%             [Y, X, Z]= sphere(16);
%             x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
%             X2 = X * r1; Y2 = Y * r1; Z2 = Z * r1;
%             h = surf(Y2 + y1,X2 + x1,Z2 + z1);
%             h.FaceAlpha=0.05;
%             h.EdgeColor="none";
%             h.FaceColor = 'blue';
%         
%             x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
%             X2 = X * r2; Y2 = Y * r2; Z2 = Z * r2;
%             h = surf(Y2 + y2,X2 + x2,Z2 + z2);
%             h.FaceAlpha = 0.05;
%             h.EdgeColor="none";
%             h.FaceColor = 'red';
        catch ME
        end

        A(Sx:Nx, Sy:Ny, Sz:Nz) = pos_fill | pos;

    end

    A(Sx:Nx, Sy:Ny, Sz:Nz) = pos_fill;
end
