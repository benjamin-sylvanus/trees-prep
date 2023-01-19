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

end
