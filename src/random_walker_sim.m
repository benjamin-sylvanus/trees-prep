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
        steplen;

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
            obj.steplen = zeros(iter, 1);
            

            while i <= iter
                obj.steplen(i) = obj.step_size;
                

                % set random next position
                obj.randomwalker.step = obj.step_size;
                obj.randomwalker = obj.randomwalker.setnext();

                % get positions
                current = obj.randomwalker.curr;
                next = obj.randomwalker.next;

                % check positions

                [obj,inside] = obj.checkpos(current, next, obj.swc, ...
                obj.lookup_table, obj.index_array, obj.pairs);

                if inside
                    %                     disp(obj.randomwalker.curr);
                    obj.randomwalker.curr = obj.randomwalker.next;
                    %                     disp(obj.randomwalker.curr);
                else
                    % disp("Curr Inside")
                end

                obj.rwpath(i, :) = obj.randomwalker.curr(:)';
                i = i + 1;
            end

        end

        function [obj,inside] = checkpos(obj, curr, next, swc, LUT, A, pairs)
            cvoxes = float2vox(curr);
            nvoxes = float2vox(next);

            if all(all(cvoxes <= obj.boundSize')) && ...
                    all(all(nvoxes < obj.boundSize')) && ...
                    all(all(cvoxes > 0)) && ...
                    all(all(nvoxes > 0))

                cvoxes = float2vox(curr);

                %             cvoxes = unique(cvoxes);

                cvoxes = sub2ind(obj.boundSize, cvoxes(1, :), ...
                cvoxes(2, :), cvoxes(3, :));

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

                indicies = obj.combineinds(cindicies, nindicies);

                % check lookup-table using subscripted index of [x y z]
                %             indicies = LUT(pos(1), pos(2), pos(3));

                if (indicies ~= 0)

                    [obj,currinside, nextinside] = check_connections(obj, 0, indicies, A, swc, curr, next);

                    if currinside && nextinside
                        %                     disp("INSIDE: ");
                        inside = 1;
                    elseif currinside && ~nextinside
                        
                        %                     disp("current in, next out");
                        inside = 0;
                    elseif ~currinside && ~nextinside
                        obj.logit(indicies,A);
                        scatter3(curr(2),curr(1),curr(3));

%                         disp("RW Outside: ");
                        inside = 0;
                    elseif ~currinside && nextinside
                        obj.logit(indicies,A);
                        disp("RW Outside NextInside: ");
                        inside = 0;
                    else
                        disp("OUTSIDE: ")
                        inside = 0;
                    end

                else
                    inside = 0;
                end

            else
                % if the step exits the bounds: check rw inside cell
                if all(all(cvoxes <= obj.boundSize')) && all(all(cvoxes > 0))
                    cvoxes = float2vox(curr);
                    cvoxes = sub2ind(obj.boundSize, cvoxes(1, :), ...
                        cvoxes(2, :), cvoxes(3, :));
                    % checkvoxes correct
                    cindicies = LUT(cvoxes);
                    indicies = obj.combineinds(cindicies, cindicies);

                    if (indicies ~= 0)
                        currinside = checkone_connection(obj, indicies, A, swc, curr);

                        if currinside
                            % if there was a conn but rw is inside cell
                            % prevent step.
                            inside = 0;
                        else
                            % if there was a connection and rw is outside
                            % then step is allowed
                            inside = 1;
                        end

                    else
                        % if the were no connections, rw is outside and
                        % should be allowed to step.
                        inside = 1;
                    end

                else
                    inside = 0;
                    disp("NEEDS FIXING: randomwalker is inside range but step goes out");
                end

                % disp("OUTSIDE RANGE");

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

        function [obj, currinside, nextinside] = check_connections(obj, flag, indicies, A, swc, curr, next)
            x0 = curr(1); y0 = curr(2); z0 = curr(3);
            nx0 = next(1); ny0 = next(2); nz0 = next(3);

            currinside = false;
            nextinside = false;

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
                tcn = pointbasedswc2v([x0 nx0], [y0 ny0], [z0 nz0], x1, x2, y1, y2, z1, z2, r1, r2);
%                 if tcn(1) && tcn(2)
%                    cdists = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
%                    cdists = min([r1,r2,cdists]);
% %                    fid = fopen(fullfile('distslog.txt'), 'a');
% %                     if fid == -1
% %                       error('Cannot open log file.');
% %                     end
% %                     fprintf(fid,"%f\n",cdists);
% %                     fclose(fid);
%                     obj.step_size = cdists;
%                 end
                currinside = tcn(1) | currinside;
                nextinside = tcn(2) | nextinside;
                if currinside && nextinside
                   cdists = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
                   cdists = min([r1,r2,cdists]);
                   obj.step_size = cdists/5;
%                     disp("early exit: ");
                    break;
                end

%                 nextinside = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2) | nextinside;
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
                inside = pointbasedswc2v(x0, y0, z0, x1, x2, y1, y2, z1, z2, r1, r2) | inside;
            end

        end

        function logit(obj, connections,A)
            % get pairs from A
            pairlist = A{connections};
            ps = obj.pairs(pairlist, :);

            % pairlist => [child, parent]
            children = ps(:, 1);
            parents = ps(:, 2);
            allnodes = unique([children;parents]);
            fid = fopen(fullfile('YourLogFile.txt'), 'a');
            if fid == -1
              error('Cannot open log file.');
            end
%             fprintf(fid,"%s\n",datetime(now,0));
            for i = 1:length(allnodes)
                fprintf(fid,"%d\n",allnodes(i));
            end
%             fprintf(fid,"-------------------\n");
            fclose(fid);

        end

    end

end
