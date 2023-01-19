function setresults(path)

    % data of size (step_num, particle_num, 3);
    load(path + "/DATA.mat");

    % TODO get correct time step this is a temp value.
    tstep = 2e-4;

    % time (ms) of each step;
    t = tstep * (1:size(data, 1))';

    % load the initial positions of the particles;
    initpos = load(path + "/initialpos.mat").data;
    ds = size(data);

    % reshape matrix of init positions;
    initpos = reshape(initpos, 1, ds(2), ds(3));

    % repeat copies of init x,y,z for each step: for matrix ops
    pos_repeated = repmat(initpos, ds(1), 1, 1);

    % ({xi,yi,zi}-{x0,y0,z0})^2
    % calculate displacement;
    dis = (data - pos_repeated) .^ 2;

    % isolate x y z vectors;
    dispX = dis(:, :, 1); dispY = dis(:, :, 2); dispZ = dis(:, :, 3);

    % calculate mean x y z: column-wise
    mdx = mean(dispX, 2); mdy = mean(dispY, 2); mdz = mean(dispZ, 2);

    % Diffusivity defined as sqrt(mean displacement)/2*t
    diffuse = @(x, y) (sqrt(x)) ./ (2 * y);
    diffx = diffuse(mdx, t); diffy = diffuse(mdy, t); diffz = diffuse(mdz, t);
    resultsPath = path + "/results";

    mkdir(resultsPath);
    resultsObj = matfile(resultsPath + "/results.mat");
    resultsObj.Properties.Writable = true;
    resultsObj.t = t;
    resultsObj.diffx = diffx;
    resultsObj.diffy = diffy;
    resultsObj.diffz = diffz;

    clearvars -except resultsPath path;

    load(resultsPath + "/results.mat");
    close all;

    f = figure("Name", "Diffusivity X Y Z");

    hold on;

    plot(diffx, t);
    plot(diffy, t);
    plot(diffz, t);

    savefig(f, resultsPath + "/" + (f.Name));

    clearvars -except path resultsPath; close all;
    load(path + "/DATA.mat");
    dataunique = cell(size(data, 2), 1);

    for i = 1:size(data, 2)
        q = data(:, i, :);
        q = reshape(q, size(data, [1, 3]));
        size(q);
        dataunique{i} = unique(q, "rows", "stable");
    end

    f = figure("Name", "Particle Paths");
    hold on;
    axis equal;
    step = 5;

    for i = 1:length(dataunique)
        rwpath = dataunique{i};
        h = plot3([rwpath(1:step:end, 2)], ...
            [rwpath(1:step:end, 1)], ...
            [rwpath(1:step:end, 3)]);
    end

    savefig(f, resultsPath + "/" + (f.Name));
