function [A, indicies, t2, LUT] = generateLUT(B, b)
    A = cell(0, 0);
    tic;
    LUT = zeros(B);
    [cx, cy, cz, ci] = cellfun(@extract_Range2, b, 'UniformOutput', false);

    bsx = repmat(B(1), size(b)); bsy = repmat(B(2), size(b)); bsz = repmat(B(3), size(b));

    bsize = [bsx(:), bsy(:), bsz(:)];
    bsize = num2cell(bsize, 2);
    inds = ci;
    inds = cellfun(@helper2, bsize, ci, 'UniformOutput', false);

    indicies = 0;
    t2 = 0;
    la = 0;

    for i = 1:length(inds)
        idx = [inds{i, 1}];
        sub = LUT(idx);
        idy = find(sub < 1);
        j = [1:length(idy)]';

        nextid = la + j;

        lj = length(j);

        if length(j) > 0

            while true
                la = nextid(length(j));

                % Add size to array
                if la > size(A, 1)
                    A = [A; cell(100000, 1)];
                end

                if la < size(A, 1)
                    break;
                end

            end

        end

        sub(idy(j)) = deal(nextid);
        idn = find(sub > 0);

        for j = 1:length(idx)
            A(sub(j), 1) = {[A{sub(j), 1}; i]};
        end

        LUT(idx) = sub;
    end

end
