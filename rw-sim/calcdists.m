function d = calcdists(tree, i)

    cxyz = tree{i, ["X", "Y", "Z"]};
    parentId = tree{i, "Parent"};
    pxyz = tree{parentId, ["X", "Y", "Z"]};
    coord_1 = cxyz;
    coord_2 = pxyz;
    difference = abs((coord_2 - coord_1) .^ 2);
    d_squared = sum(difference);
    d = d_squared;
    d = sqrt(d_squared);
end
