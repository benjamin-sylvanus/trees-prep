function vox = float2vox(pos)
    fpos = floor(pos);
    cpos = ceil(pos);

    ci = enumerate_vox(cpos, fpos);

    vox = ci;

    function ci = enumerate_vox(cpos, fpos)
        xi = fpos(1); yi = fpos(2); zi = fpos(3);
        xe = cpos(1); ye = cpos(2); ze = cpos(3);
        rx = xi:xe;
        ry = yi:ye;
        rz = zi:ze;

        ci = distributeelem(rx, ry, rz);

        function ci = distributeelem(rx, ry, rz)
            csize = prod([size(rx, 2), size(ry, 2), size(rz, 2)]);
            cx = repmat(rx', length(ry), 1);
            cy = repelem(ry, length(rx));
            cx_cy = [cx, cy'];
            cxy = repmat(cx_cy, length(rz), 1);
            cz = repelem(rz, size(cx_cy, 1));
            scxy = [size(cxy)];
            scz = [size(cz)];
            ci = [cxy, cz'];
        end

        % if a steps previous voxels trace to inside and these voxels
        % aren't repr in curr then an error will occur.  enumerate all
        % min-max voxel combinations

    end

end
