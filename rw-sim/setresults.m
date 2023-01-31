function setresults(path, vsize, tstep)

    % data of size (step_num, particle_num, 3);
    load(path + "/DATA.mat"); % unitless

    %{
    % TODO get correct time step this is a temp value.
    % fig 4 LUT intrinsic diff

%     fiermens_lee_cookbook;
    % step size: 3d: sqrt(6*intrind0* tstep);
    % all in micrometer and millisecond.
    % t = stepno * time_step;
    % ((min length^2)/6)/d0 = tstep;
    % 10% min radius. (filtered nodes)

    % d = 0.1*min(filtered_radius);
    % or
    % d = (min length^2) d/6/d0

    % tstep = d/6/D0;
    % D0 is 2um^2 / ms

    % nm->um = nm*1e3;
    % swc .* 1e-3;

    % t = step_no * tstep;
    % % step size: 3d: sqrt(6 * d0 * tstep);

    % Nstep = optional
    % Npar = optional

    % Nstep = 1e3 (minimal);
    % Npar = 1e4 (minimal);

    % Use above equations to define tstep after
    % diff  -> 1/t;
    % slope = r^2/5;

    %}

    % time (ms) of each step;
    t = (tstep) .* (1:size(data, 1))'; % ms

    % load the initial positions of the particles;
    %     initpos = load(path + "/initialpos.mat").data;
    initpos = load(path + "/initialpos.mat").data; % unitless
    ds = size(data);

    % reshape matrix of init positions;
    initpos = reshape(initpos, 1, ds(2), ds(3));

    % repeat copies of init x,y,z for each step: for matrix ops
    pos_repeated = repmat(initpos, ds(1), 1, 1);

    % ({xi,yi,zi}-{x0,y0,z0})^2
    % calculate displacement;

    dis_1_2 = (data - pos_repeated) .* vsize; % (unitless - unitless) * um
    dis = dis_1_2 .^ 2; % um^2

    % isolate x y z vectors;
    dispX = (dis(:, :, 1)); % um^2
    dispY = (dis(:, :, 2)); % um^2
    dispZ = (dis(:, :, 3)); % um^2

    % calculate mean x y z: column-wise
    mdx = mean(dispX, 2); mdy = mean(dispY, 2); mdz = mean(dispZ, 2);

    % Diffusivity defined as mean displacement/2*t
    diffuse = @(x, y) (x) ./ (2 * y); % (um^2/ms);
    diffx = diffuse(mdx, t); % (um^2/ms);
    diffy = diffuse(mdy, t); % (um^2/ms);
    diffz = diffuse(mdz, t); % (um^2/ms);

    dis_1_2k = (data - pos_repeated) .* vsize; % (unitless - unitless) * um

    dis_Kn = dis_1_2k .^ 4; % um^4
    dispX_Kn = (dis_Kn(:, :, 1)); % um^4
    dispY_Kn = (dis_Kn(:, :, 2)); % um^4
    dispZ_Kn = (dis_Kn(:, :, 3)); % um^4

    mdx_kn = mean(dispX_Kn, 2);
    mdy_kn = mean(dispY_Kn, 2);
    mdz_kn = mean(dispZ_Kn, 2);

    dis_Kd = dis_1_2k .^ 2; % um^2
    dispX_Kd = (dis_Kd(:, :, 1)); % um^2
    dispY_Kd = (dis_Kd(:, :, 2)); % um^2
    dispZ_Kd = (dis_Kd(:, :, 3)); % um^2

    mdx_kd = mean(dispX_Kd, 2);
    mdy_kd = mean(dispY_Kd, 2);
    mdz_kd = mean(dispZ_Kd, 2);

    % Kurt defined as (mean (displacement)^4)/(mean displacement^2)^2;
    kurt = @(x, y) ((x) ./ (y .^ 2))-3; % units
    kurtx = kurt(mdx_kn, mdx_kd); % units
    kurty = kurt(mdy_kn, mdy_kd); % units
    kurtz = kurt(mdz_kn, mdz_kd); % units


    resultsPath = path + "/results";

    mkdir(resultsPath);
    resultsObj = matfile(resultsPath + "/results.mat");
    resultsObj.Properties.Writable = true;
    resultsObj.t = t;
    resultsObj.diffx = diffx;
    resultsObj.diffy = diffy;
    resultsObj.diffz = diffz;
    resultsObj.kurtx = kurtx;
    resultsObj.kurty = kurty;
    resultsObj.kurtz = kurtz;

    clearvars -except resultsPath path;

    load(resultsPath + "/results.mat");

    fprintf("|(Diffusivity t0 / D0)-1|: %f %f %f\n", ...
    abs(([diffx(1) diffy(1) diffz(1)]./2) - 1));

    fprintf("Kurtosis t0: %f %f %f\n",[kurtx(1);kurty(1); kurtz(1)]);

    close all;

    f = figure("Name", "Diffusivity X Y Z");

    hold on;


    % (um^2/ms)/(1/ms) -> (um^2/ms) * (ms) -> slope: (um^2)
    % max(y) ~ D0;
    % slope r^5/2
    scatter(1 ./ t, diffx); % (1/ms,um^2/ms)
    scatter(1 ./ t, diffy); % (1/ms,um^2/ms)
    scatter(1 ./ t, diffz); % (1/ms,um^2/ms)

    ax = gca;
    xlabel("$\frac{1}{t}(\frac{1}{ms})$", 'fontsize', 14, "Interpreter", "latex");
    ylabel("$\frac{\textrm{mean displacement}^2}{2t}$ $(\frac{\mu m^2}{ms})$", 'fontsize', 14, "Interpreter", "latex");

    legend(["$\frac{(x-x_0)^2}{2t}$", "$\frac{(y-y_0)^2}{2t}$", "$\frac{(z-z_0)^2}{2t}$"], "Interpreter", "latex", FontSize = 20);

    savefig(f, resultsPath + "/" + (f.Name));

    f = figure("Name", "Kurtosis X Y Z");

    hold on;
    % (um^2/ms)/(1/ms) -> (um^2/ms) * (ms) -> slope: (um^2)
    % max(y) ~ D0;
    % slope r^5/2
    scatter(1 ./ t, kurtx); % (1/ms,um^2/ms)
    scatter(1 ./ t, kurty); % (1/ms,um^2/ms)
    scatter(1 ./ t, kurtz); % (1/ms,um^2/ms)
    savefig(f, resultsPath + "/" + (f.Name));

    %     clearvars -except path resultsPath; close all;
    %     load(path + "/DATA.mat");
    %     dataunique = cell(size(data, 2), 1);
    %
    %     for i = 1:size(data, 2)
    %         q = data(:, i, :);
    %         q = reshape(q, size(data, [1, 3]));
    %         size(q);
    %         dataunique{i} = unique(q, "rows", "stable");
    %     end
    %
    %     f = figure("Name", "Particle Paths");
    %     hold on;
    %     axis equal;
    %     step = 100;
    %
    %     for i = 1:length(dataunique)
    %         rwpath = dataunique{i};
    %         h = plot3([rwpath(1:step:end, 2)], ...
    %             [rwpath(1:step:end, 1)], ...
    %             [rwpath(1:step:end, 3)]);
    %     end
    %
    %     savefig(f, resultsPath + "/" + (f.Name));
