classdef random_walker_sim
    %RANDOM_WALKER_SIM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        lookup_table;
        index_array;
        pairs;
        swc;
        step_size;
        randomwalker;
        perm_prob;
        boundSize;
        rwpath = [];

    end

    methods (Access = public)

        function obj = random_walker_sim(LUT, index_array, pairs, ...
                bounds, swc, step_size, perm_prob)
            %RANDOM_WALKER_SIM Construct an instance of this class
            %   Detailed explanation goes here
            obj.lookup_table = LUT;
            obj.index_array = index_array;
            obj.pairs = pairs;
            obj.swc = swc;
            obj.step_size = step_size;
            obj.perm_prob = perm_prob;
            obj.boundSize = bounds;

            obj.randomwalker = randomwalker(1, obj.step_size, bounds, ...
                obj.swc, obj.lookup_table, obj.index_array, obj.pairs);
        end

        function obj = eventloop(obj, iter)
            i = 1;
            obj.rwpath = zeros(iter, 3);

            while i <= iter
                % set random next position
                obj.randomwalker = obj.randomwalker.setnext();

                % get positions
                current = obj.randomwalker.curr;
                next = obj.randomwalker.next;

                % check positions

                inside = obj.checkpos(current, next, obj.swc, obj.lookup_table, obj.index_array, obj.pairs);

                if inside
                    %                     disp(obj.randomwalker.curr);
                    obj.randomwalker.curr = obj.randomwalker.next;
                    %                     disp(obj.randomwalker.curr);
                else
                    %                     disp("not inside")
                end

                obj.rwpath(i, :) = obj.randomwalker.curr(:)';
                i = i + 1;
            end

        end

        function inside = checkpos(obj, curr, next, swc, LUT, A, pairs)
            cvoxes = float2vox(curr);
            nvoxes = float2vox(next);

            if all(all(cvoxes <= obj.boundSize')) && all(all(nvoxes < obj.boundSize')) && all(all(cvoxes > 0)) && all(all(nvoxes > 0))

                cvoxes = float2vox(curr);

                %             cvoxes = unique(cvoxes);

                cvoxes = sub2ind(obj.boundSize, cvoxes(1, :), cvoxes(2, :), cvoxes(3, :));
                % checkvoxes correct
                cindicies = LUT(cvoxes);

                nvoxes = float2vox(next);
                %             nvoxes = unique(nvoxes,"stable");
                try
                    nvoxes = sub2ind(obj.boundSize, nvoxes(1, :), nvoxes(2, :), nvoxes(3, :));

                catch ME
                    disp(ME);
                end

                % checkvoxes correct
                nindicies = LUT(nvoxes);

                indicies = [nindicies; cindicies];

                indicies = obj.combineinds(cindicies, nindicies);

                % check lookup-table using subscripted index of [x y z]
                %             indicies = LUT(pos(1), pos(2), pos(3));

                x0 = curr(1); y0 = curr(2); z0 = curr(3);
                nx0 = next(1); ny0 = next(2); nz0 = next(3);

                currinside = false;
                nextinside = false;

                if (indicies ~= 0)

                    % get pairs from A
                    pairlist = A{indicies};
                    ps = pairs(pairlist, :);

                    % pairlist => [child, parent]
                    children = ps(:, 1);
                    parents = ps(:, 2);

                    % for each pair: check if point is inside
                    for i = 1:length(children)

                        % get base and target ids
                        baseid = children(i); targetid = parents(i);
                        p1 = swc{baseid, 2:5};
                        p2 = swc{targetid, 2:5};
                        x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                        x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
                        currinside = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2) | currinside;
                        nextinside = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2) | nextinside;
                    end

                    if currinside && nextinside
                        %                     disp("INSIDE: ");
                        inside = 1;
                    elseif currinside && ~nextinside
                        %                     disp("current in, next out");
                        inside = 0;
                    else
                        disp("OUTSIDE: ")
                        inside = 0;
                    end

                else
                    inside = 1;
                end

            else
                disp("OUTSIDE RANGE");
                inside = 0;
            end

        end

        function inds = combineinds(~, cind, nind)

            c = isnumeric(cind);
            n = isnumeric(nind);

            if c && n
                inds = [cind; nind];
            elseif c && ~n
                inds = cind;
            elseif ~c && n
                inds = nind;
            else
                inds = 0;
            end

        end

    end

end
