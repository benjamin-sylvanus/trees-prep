function obj = eventloop(obj, iter)

    i = 1;
    chunk_iter = 1;
    obj.currstate = zeros(obj.particle_num, 1, "logical");
    obj.currstate(:) = true;
    [path, ~, matObj, initPosObj] = obj.setsimulationpath(iter, obj.particle_num);
    obj.path = path;
    obj.rwpath = zeros(iter, obj.particle_num, 3);
    PARTS = [obj.particles{:}]';
    PARTS = [[PARTS(:).curr]', [obj.currstate(:)], [PARTS(:).flag]'];
    particles = PARTS;
    parvars = struct();
    parvars.swc = obj.swc;
    parvars.index_array = obj.index_array;
    %     parvars.memoized_distance = obj.memoized_distance;
    parvars.limit = obj.particle_num;
    parvars.boundSize = obj.boundSize;
    parvars.lookup_table = obj.lookup_table;
    parvars.step = obj.step_size;

    PARVARS = parallel.pool.Constant(parvars);
    data = matObj.data;
    data(:) = 0;

    initPosObj.data(:, :) = particles(:, 1:3);

    parfor i = 1:obj.particle_num
        tic;

        rwpath = zeros(iter, 5);

        for j = 1:iter
            particles(i, :) = cellgap2(particles(i, :), ...
                PARVARS.Value.swc, ...
                PARVARS.Value.index_array, ...
                PARVARS.Value.boundSize, ...
                PARVARS.Value.lookup_table, ...
                PARVARS.Value.step); %#ok<PFOUS>
            rwpath(j, :) = particles(i, :);
        end

        t = toc;
        fprintf("Completed particle %d in %f seconds\n", i, t);
        % Write data in memory to to matfile
        data(:, i, :) = reshape(rwpath(:, 1:3), [size(data(:, i, :))]);
    end

    matObj.data(:, :, :) = data(:, :, :);

    %     while i <= iter
    %
    %         % create arrays of repeated i and and chunk_iter
    %         % for args of cellfun
    %         %         eles = num2cell(repelem(i, obj.particle_num, 1));
    %         %         chunks = num2cell(repelem(chunk_iter, obj.particle_num, 1));
    %
    %         % ~ for all particles determine state:
    %         % ~ res(i) represents state of particle(i) where state is inside 1 | outside 0
    %         %         [res, states] = cellfun(@cellgap2, particles, ...
    %         %         eles, chunks, "UniformOutput", false);
    %
    %         t = tic;
    %         parfor j = 1:PARVARS.Value.limit
    %             particles(j, :) = cellgap2(particles(j, :), ...
    %                 PARVARS.Value.swc, ...
    %                 PARVARS.Value.index_array, ...
    %                 PARVARS.Value.boundSize, ...
    %                 PARVARS.Value.lookup_table);
    %         end
    %         fprintf("parfor: %f\t seconds\n",toc(t));
    %
    %         % Update Particles with result
    %         %         obj.particles = res;
    %
    %         % Update Particle States with result
    %         %         obj.currstate = cell2mat(states);
    %
    %         % Record position of each particle
    %         for j = 1:obj.particle_num
    %             % * rwpath stores values from 1:chunkSize
    %             %             obj.rwpath(chunk_iter, j, :) = obj.particles{j}.curr(:)';
    %             obj.rwpath(chunk_iter, j, :) = particles(j, 1:3);
    %         end
    %
    %         % When a chunk is complete
    %         if chunk_iter == obj.chunkSize
    %             toc;
    %
    %             fprintf('Writing %d of %d\n\n', i, iter);
    %
    %             % Write data in memory to to matfile
    %             matObj.data((i - (obj.chunkSize - 1):i), :, :) = obj.rwpath;
    %             fprintf("Initializing Random Positions: \n");
    %             tic;
    %             % New random vectors for each particle
    %             for k = 1:obj.particle_num
    %                 % initialize randomwalker
    %                 obj.particles{k} = obj.particles{k}.init_rands();
    %             end
    %
    %             toc;
    %             fprintf("-------------------------------\n");
    %             % Reset chunk_iter
    %             chunk_iter = 0;
    %             tic;
    %         end
    %
    %         chunk_iter = chunk_iter + 1; i = i + 1;
    %     end

end
