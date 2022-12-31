function tree = read_swc(filename)
    fid = fopen(filename);
    A = textscan (fid, '%s', 'delimiter', '\n');
    A = A{1};
    fclose (fid);
    swc = [];

    for counter = 1:length (A)

        if ~isempty (A{counter}) % allow empty lines in between
            % allow comments: lines starting with #:
            if ~strcmp (A{counter} (1), '#')
                swc0 = sscanf (A{counter}, '%f')';
                swc = [swc; swc0]; %#ok<AGROW>
            end

        end

    end

    NodeID = swc(:, 1); Coords = swc(:, 3:5);
    Radii = swc(:, 6); Parents = swc(:, 7);

    tree = table(NodeID, Coords(:, 1), Coords(:, 2), Coords(:, 3), ...
        Radii, Parents, 'VariableNames', ...
        {'NodeId', 'X', 'Y', 'Z', 'Radii', 'Parent'});

    save("tree.mat", "tree");
end
