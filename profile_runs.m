clc; clearvars; close all;
tic;

addpath(genpath("../../treestoolbox"));
addpath(genpath("./"));
% addpath("/Users/bensylvanus/Library/Application Support/MathWorks/MATLAB Add-Ons/" + ...
% "Collections/random unit vector generator");

tree = read_swc('exampleTree.swc');

% clc;
tic;
% clearvars;
% load("tree.mat", "tree");
dists = zeros(height(tree), 1);

for i = 2:height(tree)
    dists(i, 1) = calcdists(tree, i);
end

toc;

tic;
[b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, 1);

toc;
tic;
[A, indicies, t2, LUT] = generateLUT(boundSize, b);
toc;

A = A(~cellfun('isempty', A));
sizes = cellfun('size', A, 1);

% get pairs from A
B = cell(size(A, 1), 2);

% extract node id [child, parent] -> linear index in B
% B(i) = {[children],[parents]}
for i = 1:length(B)
    pairlist = A{i};
    ps = pairs(pairlist, :);

    % pairlist => [child, parent]
    children = ps(:, 1);
    parents = ps(:, 2);

    B{i, 1} = children;
    B{i, 2} = parents;
end

%%
% TODO only reasonable improvement would be vectorized rw to avoid for loop
% TODO and to convert sub2ind calculations to custom function
close all;
figure();
hold on;
[~] = mainLoop(swc, zeros(boundSize), b, pairs);
axis equal;
% view(3);
hold on;

clc;
%{
% TODO implement time step determine how to walk concurrently
X.time_step = 2e-4; % Time of each step (ms)
X.step_num = 5e5; % # step
X.particle_num = 1e3; % # particle
%}
%%
clc;
step_num = 1e3;
particle_num = 1e3;
tic;
sim = random_walker_sim(LUT, B, pairs, boundSize, swc{:, :}, 1, 0, step_num, true, particle_num);

sim = sim.eventloop2(step_num);
% f = @() sim.eventloop2(step_num);
% t = timeit(f);
toc;

rwpath = sim.rwpath;
hold on;

for i = 1:sim.particle_num
    h = plot3([rwpath(:, i, 2)], [rwpath(:, i, 1)], [rwpath(:, i, 3)]);
end
