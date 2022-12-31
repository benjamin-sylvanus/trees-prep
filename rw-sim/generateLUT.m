function [A, indicies, t2, LUT] = generateLUT(B, b)
    A = cell(0, 0);
    tic;
    LUT = zeros(B);
    [ci] = cellfun(@extract_Range, b, 'UniformOutput', false);

    bsx = repmat(B(1), size(b));
    bsy = repmat(B(2), size(b));
    bsz = repmat(B(3), size(b));

    bsize = num2cell([bsx(:), bsy(:), bsz(:)], 2);
    inds = cellfun(@cf_sub2ind, bsize, ci, 'UniformOutput', false);

    indicies = 0;
    t2 = 0;
    la = 0;

    for i = 1:length(inds)
        idx = [inds{i, 1}];
        sub = LUT(idx);
        idy = find(sub < 1);
        j = (1:length(idy))';

        nextid = la + j;

        if ~isempty(j)

            while true
                la = nextid(length(j));

                % Add size to array
                if la > size(A, 1)
                    A = [A; cell(500000, 1)]; %#ok<AGROW>
                end

                if la < size(A, 1)
                    break;
                end

            end

        end

        sub(idy(j)) = deal(nextid);

        for j = 1:length(idx)
            A(sub(j), 1) = {[A{sub(j), 1}; i]};
        end

        LUT(idx) = sub;
    end

    function inds = cf_sub2ind(bsize, ci)
        inds = sub2ind(bsize, ci(:, 1), ci(:, 2), ci(:, 3));
    end

end
