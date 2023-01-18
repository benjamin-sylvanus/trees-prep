%%
clc; clearvars; close all;
[LUT, B, pairs, boundSize, swc, memoized_distance] = prepSim();

swcmat = swc{:, :};
clear sim;
clc;
% at least 1e3 steps
% at least 1e4 particles
step_num = 1e4;
particle_num = 1e3;
% dependent on geometry and diffusion time
step_size = 1;
perm_prob = 0;
init_in = true;
tic;

sim = random_walker_sim(LUT, B, pairs, boundSize, swcmat, step_size, ...
    perm_prob, step_num, init_in, particle_num, memoized_distance(:, 3));

sim = eventloop(sim, step_num);

toc;

% fid=fopen('/Users/benjaminsylvanus/CLionProjects/testclion/data/LUT.bin','w');
% fwrite(fid,lut,"uint16","b");
%%

%%
clearvars;
load("simulations/sim5/DATA.mat");
tstep = 2e-4;
t = tstep * (1:size(data, 1))';

initpos = load("simulations/sim5/initialpos.mat").data;

pos_repeated = repelem(initpos, size(data, 1), 1);
pos_repeated = reshape(pos_repeated, size(data));

displacement = (data - pos_repeated) .^ 2;

displaceX = displacement(:, :, 1);
displaceY = displacement(:, :, 2);
displaceZ = displacement(:, :, 3);

meanDispX = mean(displaceX, 2);
meanDispY = mean(displaceY, 2);
meanDispZ = mean(displaceZ, 2);

diffx = (sqrt(meanDispX)) ./ (2 * t);
diffy = (sqrt(meanDispY)) ./ (2 * t);
diffz = (sqrt(meanDispZ)) ./ (2 * t);

close all;

figure("Name", "Diffusivity X");
plot(diffx, t);
figure("Name", "Diffusivity Y");
plot(diffy, t);
figure("Name", "Diffusivity Z");
plot(diffz, t);

figure("Name", "Diffusivity X Y Z");
hold on;
plot(diffx, t);

plot(diffy, t);

plot(diffz, t);

% (r,c)
% r: particle; c: steps;

% sqrt(mean(displacement))/2t
%%
clearvars;
load("simulations/sim5/DATA.mat");
dataunique = cell(size(data, 2), 1);

for i = 1:size(data, 2)
    q = data(:, i, :);
    q = reshape(q, size(data, [1, 3]));
    size(q);
    dataunique{i} = unique(q, "rows", "stable");
end

hold on;
axis equal;
step = 1;

for i = 1:length(dataunique)
    rwpath = dataunique{i};
    h = plot3([rwpath(1:step:end, 2)], ...
        [rwpath(1:step:end, 1)], ...
        [rwpath(1:step:end, 3)]);
end
