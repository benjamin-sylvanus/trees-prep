function obj = eventloop(obj, iter)
    [path, ~, matObj, initPosObj] = obj.setsimulationpath(iter, ...
        obj.particle_num);
    obj.path = path;
    i = 1;
    chunk_iter = 1000;
    obj.currstate = zeros(obj.particle_num, 1, "logical");
    obj.currstate(:) = true;
    %     obj.rwpath = zeros(iter, obj.particle_num, 3);
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
    parvars.perm_prob = obj.perm_prob;

    PARVARS = parallel.pool.Constant(parvars);
    %     data = matObj.data;
    %     data(:) = 0;

    initPosObj.data(:, :) = particles(:, 1:3);

    for batch = 1:chunk_iter:obj.particle_num
        chunk = batch:batch + chunk_iter - 1;
        data = matObj.data(:, chunk, :);
        data(:) = 0;
        parties = particles(chunk, :);
        t0 = tic;

        parfor i = 1:numel(chunk)
            rwpath = zeros(iter, 5);
            swc = PARVARS.Value.swc;
            index_array = PARVARS.Value.index_array;
            boundSize = PARVARS.Value.boundSize;
            lookup_table = PARVARS.Value.lookup_table;
            step = PARVARS.Value.step;
            perm_prob = PARVARS.Value.perm_prob;

            rands = random_unit_vector(3, iter);

            for j = 1:iter
                crand = rands(:, j);
                parties(i, :) = cellgap2(parties(i, :), ...
                    swc, ...
                    index_array, ...
                    boundSize, ...
                    lookup_table, ...
                    step, ...
                    perm_prob, ...
                    crand);
                rwpath(j, :) = parties(i, :);
            end

            % Write data in memory to to matfile

            data(:, i, :) = reshape(rwpath(:, 1:3), [size(data(:, i, :))]);
        end

        matObj.data(:, chunk, :) = data(:, :, :);
        t = toc(t0);
        fprintf("Completed batch (%d:%d)x%d steps in %f seconds\n", ...
            batch, batch + chunk_iter - 1, iter, t);
    end

end
