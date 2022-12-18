function swc = read_t(filename)
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
                swc = [swc; swc0];
            end

        end

    end

end
