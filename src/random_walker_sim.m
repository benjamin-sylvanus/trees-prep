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
        rwpath;
        steplen;
        logical_alloc;
        voxels;
        currstate;

    end

    methods (Access = public)

        function obj = random_walker_sim(LUT, index_array, pairs, ...
                bounds, swc, step_size, perm_prob, iter, init_in)
            %RANDOM_WALKER_SIM Construct an instance of this class
            %   Detailed explanation goes here

            % initialize vars
            obj.lookup_table = LUT;
            obj.index_array = index_array;
            obj.pairs = pairs;
            obj.swc = swc;
            obj.step_size = step_size;
            obj.perm_prob = perm_prob;
            obj.boundSize = bounds;
            obj.logical_alloc = zeros(1, 1, "logical");
            obj.voxels = cell(iter, 2);
            obj.currstate = init_in;
            obj.randomwalker = randomwalker(init_in, obj.step_size, bounds, ...
                obj.swc, obj.lookup_table, obj.index_array, obj.pairs, iter);
        end

        function obj = eventloop(obj, iter)
            i = 1;
            obj.rwpath = zeros(iter, 3);
            % obj.steplen = zeros(iter, 1);
            obj.randomwalker.step = obj.step_size;

            while i <= iter

                % init step_size
                % obj.steplen(i) = obj.step_size;
                % obj.randomwalker.step = obj.step_size;
                % set random next position

                obj.randomwalker = obj.randomwalker.setnext(i);

                % get positions
                current = obj.randomwalker.curr;
                next = obj.randomwalker.next;

                % check positions
                [obj, inside] = obj.checkpos(current, next, obj.swc, ...
                obj.lookup_table, obj.index_array, obj.pairs, i);

                % inside -> update vars
                if inside
                    obj.randomwalker.curr = next;
                    obj.currstate = true;
                else
                    % random walker was outside or crossed boundary
                    if (rand < obj.perm_prob)
                        obj.randomwalker.curr = next;
                        obj.currstate = false;
                    end

                    % TODO implement this outcome:
                    % ^ LIMIT STEP WHEN:
                    % * crosses out
                    % ^ ENABLE STEP WHEN:
                    % * crosses in
                    % * remains out
                end

                obj.rwpath(i, :) = obj.randomwalker.curr(:)';
                i = i + 1;
            end

        end

        function [obj, inside] = checkpos(obj, curr, next, swc, LUT, A, pairs, i)

            indicies = obj.preprocesses(curr, next);

            if (indicies ~= 0)

                [obj, currinside, nextinside] = check_connections(obj, 0, indicies, A, swc, curr, next);

                inside = obj.insideLogic(currinside, nextinside);
            else
                inside = 0;
            end

        end

        function inds = combineinds(~, cind, nind)

            cidx = cind > 0;
            nidx = nind > 0;

            tc = cind(cidx);
            tn = nind(nidx);

            c = ~isempty(tc);
            n = ~isempty(tn);

            if c && n
                inds = [tc; tn];
            elseif c && ~n
                inds = tc;
            elseif ~c && n
                inds = tn;
            else
                inds = 0;
            end

        end

        function inds = combineind(~, nind)

            nidx = nind > 0;

            tn = nind(nidx);

            n = ~isempty(tn);

            if n
                inds = tn;
            else
                inds = 0;
            end

        end

        function [obj, currinside, nextinside] = check_connections(obj, flag, indicies, A, swc, curr, next)
            %             x0 = curr(1); y0 = curr(2); z0 = curr(3);
            nx0 = next(1); ny0 = next(2); nz0 = next(3);

            if islogical(obj.currstate)
                currinside = obj.currstate;
            else
                currinside = false;
            end

            nextinside = false;

            children = A{indicies, 1};
            parents = A{indicies, 2};

            %             % get pairs from A
            %             pairlist = A{indicies};
            %             ps = obj.pairs(pairlist, :);
            %
            %             % pairlist => [child, parent]
            %             children = ps(:, 1);
            %             parents = ps(:, 2);

            % for each pair: check if point is inside
            for i = 1:length(children)
                % get base and target ids
                % TODO extract values to store in A rather than at runtime
                baseid = children(i); targetid = parents(i);
                p1 = swc(baseid, 2:5);
                p2 = swc(targetid, 2:5);

                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

                tcn = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2, obj.logical_alloc);

                % inside ith connection or already inside
                % currinside = tcn(1) | currinside;
                nextinside = tcn | nextinside;

                % if both positions are inside:
                % break loop to stop extra iterations
                if currinside && nextinside
                    break;
                end

            end

        end

        function inside = checkone_connection(obj, indicies, A, swc, pos)
            x0 = pos(1); y0 = pos(2); z0 = pos(3);
            inside = false;

            children = A{indicies, 1};
            parents = A{indicies, 2};

            %             % get pairs from A
            %             pairlist = A{indicies};
            %             ps = obj.pairs(pairlist, :);
            %
            %             % pairlist => [child, parent]
            %             children = ps(:, 1);
            %             parents = ps(:, 2);

            % for each pair: check if point is inside
            for i = 1:length(children)
                % get base and target ids
                baseid = children(i); targetid = parents(i);
                p1 = swc(baseid, 2:5);
                p2 = swc(targetid, 2:5);
                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
                inside = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2, obj.logical_alloc) | inside;
            end

        end

        function logit(obj, connections, A)
            % get pairs from A
            pairlist = A{connections};
            ps = obj.pairs(pairlist, :);

            % pairlist => [child, parent]
            children = ps(:, 1); parents = ps(:, 2);
            allnodes = unique([children; parents]);

            fid = fopen(fullfile('YourLogFile.txt'), 'a');

            if fid == -1
                error('Cannot open log file.');
            end

            for i = 1:length(allnodes)
                fprintf(fid, "%d\n", allnodes(i));
            end

            fclose(fid);
        end

        function linearArray = rowOp(obj, row)

        end

        function indicies = preprocesses(obj, curr, next)
            LUT = obj.lookup_table;
            % convert float to voxel
            nvoxes = float2vox(next);
            nvoxes = nvoxes';
            nv = all(nvoxes > 0, 2) & all(nvoxes < obj.boundSize, 2);
            NVOXES = nvoxes(nv, :)';
            %NVOXES are always valid indicies
            NVOXES = sub2ind(obj.boundSize, NVOXES(1, :), ...
            NVOXES(2, :), NVOXES(3, :));
            nindicies = LUT(NVOXES);
            indicies = obj.combineind(nindicies');
        end

        function inside = insideLogic(obj, currinside, nextinside)
            % both positions within a connection
            if currinside && nextinside
                inside = 1;
                % disp("Was just inside")

                % current pos in, next pos out
            elseif currinside && ~nextinside
                % disp("next_outside")
                inside = 0;

                % current pos out, next pos in
            elseif ~currinside && nextinside

                disp("RW Outside NextInside: ");
                inside = 0;

                % both positions outside
            elseif ~currinside && ~nextinside
                disp("RW Outside: ");
                %                     scatter3(curr(2), curr(1), curr(3));

                inside = 0;
            else
                disp("OUTSIDE: ")
                inside = 0;
            end

        end

    end

end
