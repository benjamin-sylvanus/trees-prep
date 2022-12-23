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
                bounds, swc, step_size, perm_prob, iter)
            %RANDOM_WALKER_SIM Construct an instance of this class
            %   Detailed explanation goes here
            obj.lookup_table = LUT;
            obj.index_array = index_array;
            obj.pairs = pairs;
            obj.swc = swc;
            obj.step_size = step_size;
            obj.perm_prob = perm_prob;
            obj.boundSize = bounds;
            obj.logical_alloc = zeros(2, 1, "logical");
            obj.voxels = cell(iter, 2);

            obj.randomwalker = randomwalker(1, obj.step_size, bounds, ...
                obj.swc, obj.lookup_table, obj.index_array, obj.pairs, iter);
        end

        function obj = eventloop(obj, iter)
            i = 1;
            obj.rwpath = zeros(iter, 3);
            obj.steplen = zeros(iter, 1);

            while i <= iter

                obj.steplen(i) = obj.step_size;

                obj.randomwalker.step = obj.step_size;
                % set random next position
                obj.randomwalker = obj.randomwalker.setnext(i);

                % get positions
                current = obj.randomwalker.curr;
                next = obj.randomwalker.next;

                % check positions

                [obj, inside] = obj.checkpos(current, next, obj.swc, ...
                obj.lookup_table, obj.index_array, obj.pairs, i);

                if inside
                    obj.randomwalker.curr = next;
                    obj.currstate = true;
                else
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
            cvoxes = float2vox(curr);
            nvoxes = float2vox(next);

            zero_cv = cvoxes > 0;
            max_cv = cvoxes < obj.boundSize;
            zero_nv = nvoxes > 0;
            max_nv = nvoxes < obj.boundSize;

            cv = (zero_cv & max_cv);
            nv = (zero_nv & max_nv);

            % combine columns
            combcv = cv(:, 1) & cv(:, 2) & cv(:, 3);
            combnv = nv(:, 1) & nv(:, 2) & nv(:, 3);

            CVOXES = cvoxes(combcv, :)';
            NVOXES = nvoxes(combnv, :)';

            % CVOXES and NVOXES are always valid indicies

            %             if all(cvoxes <= obj.boundSize) && ...
            %                     all(nvoxes <= obj.boundSize) && ...
            %                     all(cvoxes > 0) && ...
            %                     all(nvoxes > 0)

            CVOXES = sub2ind(obj.boundSize, CVOXES(1, :), ...
            CVOXES(2, :), CVOXES(3, :));
            cindicies = LUT(CVOXES);

            NVOXES = sub2ind(obj.boundSize, NVOXES(1, :), ...
                NVOXES(2, :), NVOXES(3, :));
            nindicies = LUT(NVOXES);

            cindicies = cindicies(cindicies > 0);

            indicies = obj.combineinds(cindicies', nindicies');

            if i > 2
                indicies = obj.combineinds(indicies, obj.voxels{i - 1, 1}');
                indicies = obj.combineinds(indicies, obj.voxels{i - 1, 2}');
            end

            obj.voxels{i, 1} = CVOXES;
            obj.voxels{i, 2} = NVOXES;

            if (indicies ~= 0)

                [obj, currinside, nextinside] = check_connections(obj, 0, indicies, A, swc, curr, next);

                % both positions within a connection
                if currinside && nextinside
                    inside = 1;
                    %                     disp("Was just inside")

                    % current pos in, next pos out
                elseif currinside && ~nextinside
                    %                     disp("next_outside")
                    inside = 0;

                    if i > 1
                        obj.voxels{i, 1} = obj.voxels{i - 1, 1};
                        obj.voxels{i, 2} = obj.voxels{i - 1, 2};
                    end

                    % current pos out, next pos in
                elseif ~currinside && nextinside

                    %                     disp("RW Outside NextInside: ");
                    inside = 0;

                    % both positions outside
                elseif ~currinside && ~nextinside
                    disp("RW Outside: ");
                    scatter3(curr(2), curr(1), curr(3));

                    inside = 0;
                else
                    disp("OUTSIDE: ")
                    inside = 0;
                end

            else
                inside = 0;
            end

            %             else
            %                 % if the step exits the bounds: check rw inside cell
            %                 if all(all(cvoxes <= obj.boundSize')) && all(all(cvoxes > 0))
            %                     cvoxes = float2vox(curr);
            %                     cvoxes = sub2ind(obj.boundSize, cvoxes(1, :), ...
            %                         cvoxes(2, :), cvoxes(3, :));
            %                     % checkvoxes correct
            %                     cindicies = LUT(cvoxes);
            %                     indicies = obj.combineinds(cindicies, cindicies);
            %
            %                     if (indicies ~= 0)
            %                         currinside = checkone_connection(obj, indicies, A, swc, curr);
            %
            %                         if currinside
            %                             % if there was a conn but rw is inside cell
            %                             % prevent step.
            %                             inside = 0;
            %                         else
            %                             % if there was a connection and rw is outside
            %                             % then step is allowed
            %                             inside = 1;
            %                         end
            %
            %                     else
            %                         % if the were no connections, rw is outside and
            %                         % should be allowed to step.
            %                         inside = 1;
            %                     end
            %
            %                 else
            %
            %                     % TODO : currently ignoring the chance that some of the voxels are within.
            %                     % If a single voxel crosses bounds of LUT we throw out step.
            %                     % ~ SOLUTION remove filter voxes that are outside bounds.
            %
            %                     inside = 0;
            %                     disp("NEEDS FIXING: randomwalker is inside range but step goes out");
            %                 end
            %
            %             end

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

        function [obj, currinside, nextinside] = check_connections(obj, flag, indicies, A, swc, curr, next)
            x0 = curr(1); y0 = curr(2); z0 = curr(3);
            nx0 = next(1); ny0 = next(2); nz0 = next(3);

            if islogical(obj.currstate)
                currinside = obj.currstate;
            else
                currinside = false;
            end

            nextinside = false;

            % get pairs from A
            pairlist = A{indicies};
            ps = obj.pairs(pairlist, :);

            % pairlist => [child, parent]
            children = ps(:, 1);
            parents = ps(:, 2);

            % for each pair: check if point is inside
            % TODO: store previous inside or outside bool: this will reduce the number of connections checked
            for i = 1:length(children)
                % get base and target ids
                baseid = children(i); targetid = parents(i);
                p1 = swc(baseid, 2:5);
                p2 = swc(targetid, 2:5);

                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

                tcn = pointbasedswc2v([x0 nx0], [y0 ny0], [z0 nz0], x1, x2, y1, y2, z1, z2, r1, r2, obj.logical_alloc);

                % inside ith connection or already inside
                currinside = tcn(1) | currinside;
                nextinside = tcn(2) | nextinside;

                % if both positions are inside:
                % break loop to stop extra iterations
                if currinside && nextinside
                    % adaptive stepsize
                    cdists = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
                    cdists = min([r1, r2, cdists]);
                    obj.step_size = cdists;

                    break;
                end

            end

        end

        function inside = checkone_connection(obj, indicies, A, swc, pos)
            x0 = pos(1); y0 = pos(2); z0 = pos(3);
            inside = false;

            % get pairs from A
            pairlist = A{indicies};
            ps = obj.pairs(pairlist, :);

            % pairlist => [child, parent]
            children = ps(:, 1);
            parents = ps(:, 2);

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

    end

end
