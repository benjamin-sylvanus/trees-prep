function setresults(path,vsize,tstep)

    % data of size (step_num, particle_num, 3);
    load(path + "/DATA.mat");

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

    % time (ms) of each step;
    t = tstep * (1:size(data, 1))';

    % load the initial positions of the particles;
    initpos = load(path + "/initialpos.mat").data;
    ds = size(data);

    % reshape matrix of init positions;
    initpos = reshape(initpos.*vsize, 1, ds(2), ds(3));

    % repeat copies of init x,y,z for each step: for matrix ops
    pos_repeated = repmat(initpos.*vsize, ds(1), 1, 1);

    % ({xi,yi,zi}-{x0,y0,z0})^2
    % calculate displacement;
    dis = (data - pos_repeated) .^ 2;

    % isolate x y z vectors;
%     dispX = (dis(:, :, 1).*vsize)./8;
%     dispY = (dis(:, :, 2).*vsize)./8;
%     dispZ = (dis(:, :, 3).*vsize)./33;
    dispX = (dis(:, :, 1).*vsize);
    dispY = (dis(:, :, 2).*vsize);
    dispZ = (dis(:, :, 3).*vsize);

    % calculate mean x y z: column-wise
    mdx = mean(dispX, 2); mdy = mean(dispY, 2); mdz = mean(dispZ, 2);

    % Diffusivity defined as mean displacement/2*t
    diffuse = @(x, y) (x) ./ (2 * y);
    % micrometer^2/ ms
    % d^2/t
    diffx = diffuse(mdx, t); diffy = diffuse(mdy, t); diffz = diffuse(mdz, t);
    resultsPath = path + "/results";

%     figure()
% 
%     plot(sqrt(t),diffx);

    % xlim([0,1]);

    % use micron
    % scale into micrometer and ms
    % 2um^2 / ms

    
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

%     plot(diffx, 1/t);
%     plot(diffy, 1/t);
%     plot(diffz, 1/t);
    scatter(1./t,diffx);
    scatter(1./t,diffy);
    scatter(1./t,diffz);

    ax = gca;
    xlabel("$\frac{1}{t}(\frac{1}{ms})$",'fontsize',14,"Interpreter","latex");
    ylabel("$\frac{\textrm{mean displacement}^2}{2t}$ $(\frac{\mu m^2}{ms})$",'fontsize',14,"Interpreter","latex");

    legend(["$\frac{(x-x_0)^2}{2t}$","$\frac{(y-y_0)^2}{2t}$","$\frac{(z-z_0)^2}{2t}$"],"Interpreter","latex",FontSize=20);


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
    step = 1;

    for i = 1:length(dataunique)
        rwpath = dataunique{i};
        h = plot3([rwpath(1:step:end, 2)], ...
            [rwpath(1:step:end, 1)], ...
            [rwpath(1:step:end, 3)]);
    end

    savefig(f, resultsPath + "/" + (f.Name));
