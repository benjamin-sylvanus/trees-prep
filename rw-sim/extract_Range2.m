function [cx, cy, cz, ci] = extract_Range2(sub)
    xi = sub(1, 1);
    yi = sub(1, 2);
    zi = sub(1, 3);
    xe = sub(2, 1);
    ye = sub(2, 2);
    ze = sub(2, 3);
    rx = xi:xe;
    ry = yi:ye;
    rz = zi:ze;

    r = {rx, ry, rz};

    [cx, cy, cz, ci] = distributeelem(rx, ry, rz);

    function [cx, cy, cz, ci] = distributeelem(rx, ry, rz)
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

end
