%%
clc; clearvars; close all;
[LUT, B, pairs, boundSize, swc, memoized_distance] = prepSim();

swcmat = swc{:, :};
clear sim;
clc;
% at least 1e3 steps
% at least 1e4 particles 
step_num = 1e3;
particle_num = 1e4;
% dependent on geometry and diffusion time
step_size = 0.2;
perm_prob = 0;
init_in = true;
tic;

sim = random_walker_sim(LUT, B, pairs, boundSize, swcmat, step_size, ...
    perm_prob, step_num, init_in, particle_num, memoized_distance(:, 3));

sim = eventloop(sim,step_num);

toc;

% fid=fopen('/Users/benjaminsylvanus/CLionProjects/testclion/data/LUT.bin','w');
% fwrite(fid,lut,"uint16","b");
%%

%%
clearvars;
load("simulations/sim12/DATA.mat");
tstep = 2e-4;
t = tstep * (1:size(data,1));
initpos = load("simulations/sim12/initialpos.mat").data;

pos_repeated = repelem(initpos,size(data,1),1);
pos_repeated  =reshape(pos_repeated, size(data));

displacement = data - pos_repeated;

displaceX = displacement(:,:,1);
displaceY = displacement(:,:,2);
displaceZ = displacement(:,:,3);

x_example = displaceX(1,:);
y_example = displaceY(1,:);
z_example = displaceZ(1,:);



sqrt(mean(x_example))
xres = zeros(size(x_example));
yres = zeros(size(y_example));
zres = zeros(size(z_example));

for i = 1:size(x_example,2)
    xres(i) = sqrt(mean(abs(x_example(1:i))));
    yres(i) = sqrt(mean(abs(y_example(1:i))));
    zres(i) = sqrt(mean(abs(z_example(1:i))));
end

diffx = xres./(2*t);
diffy = yres./(2*t);
diffz = zres./(2*t);

clf;
hold on;
plot(diffx,t);
plot(diffy,t);
plot(diffz,t);
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
