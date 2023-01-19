%%
clc; clearvars; close all;
[LUT, B, pairs, boundSize, swc, memoized_distance] = prepSim();

swcmat = swc{:, :};
clear sim;
clc;
% at least 1e3 steps
% at least 1e4 particles
step_num = 1e2;
particle_num = 1e5;
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

path = sim.path;
% path = "./simulations/sim5";
clearvars -except path;
setresults(path);
