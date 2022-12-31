classdef random_walker_sim
    %RANDOM_WALKER_SIM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        lookup_table;
        index_array;
        pairs;
        swc;
        step_size;
        perm_prob;
        boundSize;
        rwpath;
        steplen;
        logical_alloc;
        voxels;
        currstate;
        particles;
        particle_num;

    end

    methods (Access = public)

        function obj = random_walker_sim(LUT, index_array, pairs, ...
                bounds, swc, step_size, perm_prob, iter, init_in, particle_num)
            %RANDOM_WALKER_SIM Constructor
            %   Initialize object for Monte Carlo Simulations

            % init vars
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
            obj.particle_num = particle_num;
            obj.particles = obj.init_particles(iter);
        end

        function obj = eventloop(obj, iter)
            i = 1;
            obj.currstate = zeros(obj.particle_num, 1, "logical");
            obj.currstate(:) = true;
            obj.rwpath = zeros(iter, obj.particle_num, 3);
            obj.randomwalker.step = obj.step_size;
            % fid = fopen("log.txt","a+");
            fid = fopen("X.txt", "a+");

            while i <= iter

                if mod(i, 1000) == 0
                    toc;
                    tic;
                end

                lim = obj;

                for j = 1:lim.particle_num

                    obj.particles{j} = obj.particles{j}.setnext(i);

                    % get positions
                    current = obj.particles{j}.curr;
                    next = obj.particles{j}.next;

                    % check positions
                    [obj, inside] = obj.checkpos(current, next, obj.swc, obj.index_array, j);

                    % inside -> update vars
                    if inside
                        obj.particles{j}.curr = next;
                        obj.currstate(j) = true;
                    else
                        % random walker was outside or crossed boundary
                        if (rand < obj.perm_prob)
                            obj.particles(j).curr = next;
                            obj.currstate(j) = false;
                        end

                        % TODO implement this outcome:
                        % ^ LIMIT STEP WHEN:
                        % * crosses out
                        % ^ ENABLE STEP WHEN:
                        % * crosses in
                        % * remains out
                    end

                    obj.rwpath(i, j, :) = obj.particles{j}.curr(:)';
                end

                i = i + 1;
            end

            fclose(fid);

        end

        function [obj, inside] = checkpos(obj, curr, next, swc, A, j)
            % extract {child,parent} Ids of segments near location of randomwalker
            indicies = obj.preprocesses(curr, next);

            % if there are elements of the cell near the random walker
            if (indicies ~= 0)

                % check if the random walker is inside these connections
                [obj, currinside, nextinside] = check_connections(obj, 0, indicies, A, swc, curr, next, j);

                % determine whether boundary was crossed
                % using states of current position and next
                inside = insideLogic(currinside, nextinside);
            else
                % if no elements are near the randomwalker then no boundaries were crossed
                inside = 0;
            end

        end

        function inds = combineind(~, nind)

            % extract non-zero elements
            nidx = nind > 0; tn = nind(nidx);

            % logicals for non empty index arrays
            n = ~isempty(tn);

            if n
                % if there are non-zero elements return them
                inds = tn;
            else
                % return 0
                inds = 0;
            end

        end

        function [obj, currinside, nextinside] = check_connections(obj, ~, indicies, A, swc, ~, next, j)
            % initialize current and next states
            currinside = obj.currstate(j); nextinside = false;

            % extract child and parent ids
            children = A{indicies, 1}; parents = A{indicies, 2};

            % set xyz for inside-outside calculation
            nx0 = next(1); ny0 = next(2); nz0 = next(3);

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
                nextinside = tcn | nextinside;

                % if both positions are inside: break
                if currinside && nextinside
                    break;
                end

            end

        end

        function indicies = preprocesses(obj, ~, next)
            % convert float to voxel
            nvoxes = float2vox(next)';

            % logical array for elements within range 0 < n < boundSize
            nv = all(nvoxes > 0, 2) & all(nvoxes < obj.boundSize, 2);

            % extract elements within range
            NVOXES = nvoxes(nv, :)';

            % convert subscript index to linear index
            % ^ (faster than using sub index to get values from lookup table)
            NVOXES = sub2ind(obj.boundSize, NVOXES(1, :), NVOXES(2, :), NVOXES(3, :));

            % extract [child,parent] node ids
            nindicies = obj.lookup_table(NVOXES);

            % return non-zero indicies or return 0
            indicies = obj.combineind(nindicies');
        end

        function particles = init_particles(obj, iter)
            % create n-randomwalkers
            particles = cell(obj.particle_num, 1);
            tic;

            for i = 1:obj.particle_num
                % initialize randomwalker
                particles{i} = randomwalker(true, obj.step_size, obj.boundSize, ...
                obj.swc, obj.lookup_table, obj.index_array, obj.pairs, iter);
            end

            fprintf("Particle Init: %f seconds\n", toc);
        end

        function [particle, currstate] = cellgap(obj, particle, currstate, i)
            % get positions
            current = particle.curr;
            next = particle.setnext(i);

            % check positions
            [~, inside] = obj.checkpos2(current, next, obj.swc, obj.index_array, currstate);

            if inside
                particle.curr = next;
                currstate = true;
            else
                % random walker was outside or crossed boundary
                if (rand < obj.perm_prob)
                    particle.curr = next;
                    currstate = false;
                end

                % TODO implement this outcome:
                % ^ LIMIT STEP WHEN:
                % * crosses out
                % ^ ENABLE STEP WHEN:
                % * crosses in
                % * remains out
            end

        end

        function obj = eventloop2(obj, iter)
            i = 1;
            obj.currstate = zeros(obj.particle_num, 1, "logical");
            obj.currstate(:) = true;
            obj.rwpath = zeros(iter, obj.particle_num, 3);
            % fid = fopen("log.txt","a+");
            %             fid = fopen("X.txt", "a+");

            while i <= iter

                if mod(i, 1000) == 0
                    toc;
                    % fprintf(fid, "I: %d\t%f\n", i, toc);
                    tic;
                end

                eles = num2cell(repelem(i, obj.particle_num, 1));
                [res, states] = cellfun(@obj.cellgap, obj.particles, num2cell(obj.currstate), eles, "UniformOutput", false);
                obj.particles = res;
                obj.currstate = cell2mat(states);

                for j = 1:obj.particle_num
                    obj.rwpath(i, j, :) = obj.particles{j}.curr(:)';
                    %                     fprintf(fid, "%d ", obj.particles{j}.curr(1));
                end

                %                 fprintf(fid,"\n");

                i = i + 1;
            end

            %             fclose(fid);

        end

        function [obj, inside] = checkpos2(obj, curr, next, swc, A, currstate)
            % extract {child,parent} Ids of segments near location of randomwalker
            indicies = obj.preprocesses(curr, next);

            % if there are elements of the cell near the random walker
            if (indicies ~= 0)

                % check if the random walker is inside these connections
                [obj, currinside, nextinside] = check_connections2(obj, 0, indicies, A, swc, curr, next, currstate);

                % determine whether boundary was crossed
                % using states of current position and next
                inside = obj.insideLogic(currinside, nextinside);
            else
                % if no elements are near the randomwalker then no boundaries were crossed
                inside = 0;
            end

        end

        function [obj, currinside, nextinside] = check_connections2(obj, ~, indicies, A, swc, ~, next, currstate)
            % initialize current and next states
            currinside = currstate; nextinside = false;

            % extract child and parent ids
            children = A{indicies, 1}; parents = A{indicies, 2};

            % set xyz for inside-outside calculation
            nx0 = next(1); ny0 = next(2); nz0 = next(3);

            swc1 = swc(children, 2:5);
            swc2 = swc(parents, 2:5);

            % for each pair: check if point is inside
            for i = 1:length(children)
                % get base and target ids
                % TODO extract values to store in A rather than at runtime
                %                 baseid = children(i); targetid = parents(i);
                %                 p1 = swc(baseid, 2:5);
                %                 p2 = swc(targetid, 2:5);
                p1 = swc1(i, :);
                p2 = swc2(i, :);

                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

                tcn = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2, obj.logical_alloc);

                % inside ith connection or already inside
                nextinside = tcn | nextinside;

                % if both positions are inside: break
                if currinside && nextinside
                    break;
                end

            end

        end

    end

    methods (Static)

        function inside = insideLogic(currinside, nextinside)
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

                inside = 0;
            else
                disp("OUTSIDE: ")
                inside = 0;
            end

        end

    end

end
