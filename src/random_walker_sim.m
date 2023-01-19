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
        currstate;
        particles;
        particle_num;
        memoized_distance;
        chunkSize;
        path;
    end

    methods (Access = public)

        function obj = random_walker_sim(LUT, index_array, pairs, ...
                bounds, swc, step_size, perm_prob, iter, init_in, ...
                particle_num, memoized_distance)
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
            obj.logical_alloc = false;
            obj.currstate = init_in;
            obj.particle_num = particle_num;
            obj.chunkSize = 10;
            obj.particles = obj.init_particles(iter, obj.chunkSize);
            obj.memoized_distance = memoized_distance;
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

        function particles = init_particles(obj, iter, chunk)
            % create n-randomwalkers
            particles = cell(obj.particle_num, 1);
            tic;

            for i = 1:obj.particle_num
                % initialize randomwalker
                particles{i} = randomwalker(true, obj.step_size, obj.boundSize, ...
                obj.swc, obj.lookup_table, obj.index_array, obj.pairs, iter, chunk);
            end

            fprintf("Particle Init: %f seconds\n", toc);
        end

        function [particle, currstate] = cellgap(obj, particle, currstate, i, chunk)

            % flag would be arr(5);
            % state is      arr(4);

            state = particle(4);
            flag = particle(5);

            if (flag)
                particle(5) = false;
            else
                % get positions
                current = particle(1:3);

                next = particle.setnext(chunk);

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
                    else
                        particle.flag = true;
                    end

                    % TODO implement this outcome:
                    % ^ LIMIT STEP WHEN:
                    % * crosses out
                    % ^ ENABLE STEP WHEN:
                    % * crosses in
                    % * remains out

                end

            end

        end

        function obj = eventloop2(obj, iter)
            i = 1;
            chunk_iter = 1;
            obj.currstate = zeros(obj.particle_num, 1, "logical");
            obj.currstate(:) = true;
            [path, ~, matObj, initPosObj] = obj.setsimulationpath(iter, obj.particle_num);
            obj.path = path;
            obj.rwpath = zeros(iter, obj.particle_num, 3);

            while i <= iter

                % create arrays of repeated i and and chunk_iter
                % for args of cellfun
                eles = num2cell(repelem(i, obj.particle_num, 1));
                chunks = num2cell(repelem(chunk_iter, obj.particle_num, 1));

                % ~ for all particles determine state:
                % ~ res(i) represents state of particle(i) where state is inside 1 | outside 0
                [res, states] = cellfun(@obj.cellgap, obj.particles, num2cell(obj.currstate), ...
                eles, chunks, "UniformOutput", false);

                % Update Particles with result
                obj.particles = res;

                % Update Particle States with result
                obj.currstate = cell2mat(states);

                % Record position of each particle
                for j = 1:obj.particle_num
                    % * rwpath stores values from 1:chunkSize
                    obj.rwpath(chunk_iter, j, :) = obj.particles{j}.curr(:)';
                end

                % When a chunk is complete
                if chunk_iter == obj.chunkSize
                    toc;

                    fprintf('Writing %d of %d\n\n', i, iter);

                    % Write data in memory to to matfile
                    matObj.data((i - (obj.chunkSize - 1):i), :, :) = obj.rwpath;
                    fprintf("Initializing Random Positions: \n");
                    tic;
                    % New random vectors for each particle
                    for k = 1:obj.particle_num
                        % initialize randomwalker
                        obj.particles{k} = obj.particles{k}.init_rands();
                    end

                    toc;
                    fprintf("-------------------------------\n");
                    % Reset chunk_iter
                    chunk_iter = 0;
                    tic;
                end

                chunk_iter = chunk_iter + 1; i = i + 1;
            end

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
                p1 = swc1(i, :);
                p2 = swc2(i, :);

                x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
                x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

                dist = obj.memoized_distance(children(i));

                if r1 > r2
                    tcn = pointbasedswc2v(nx0, ny0, nz0, x1, x2, y1, y2, z1, z2, r1, r2, obj.logical_alloc, dist);

                else
                    tcn = pointbasedswc2v(nx0, ny0, nz0, x2, x1, y2, y1, z2, z1, r2, r1, obj.logical_alloc, dist);
                end

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

        % Todo Actually implement a file system this is just a temp method for Chunk testing
        function [path, fileName, matObj, initPosObj] = setsimulationpath(iter, particle_num)
            % Get the folders containing sim
            folders = dir("./simulations/sim*");
            % Set the next folder name
            next = numel(folders);
            str = sprintf("simulations/sim%d", next);

            % Create Folder
            mkdir(str);
            path = str;

            % Filename of Simulation results
            fileName = str + '/DATA.mat';
            initialPos = str + '/initialPos.mat';

            % Init matfile with memory alloc of (step_num,particle_num,[x,y,z])
            % Enables write to portion of matfile without loading entire file into memory
            matObj = matfile(fileName);
            matObj.Properties.Writable = true;
            matObj.data(iter, particle_num, 3) = 0;
            initPosObj = matfile(initialPos);
            initPosObj.Properties.Writable = true;
            initPosObj.data(particle_num, 3) = 0;

        end

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
