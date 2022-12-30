clc; clearvars; close all;
tic;

addpath(genpath("../../treestoolbox"));
addpath(genpath("./"));
addpath("/Users/bensylvanus/Library/Application Support/MathWorks/MATLAB Add-Ons/" + ...
"Collections/random unit vector generator");

tree = read_t('exampleTree.swc');

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
[b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, 0.8);

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
[~, ~] = mainLoop(swc, zeros(boundSize), b, pairs);
axis equal;
% view(3);
hold on;

tic;
clc;
iter = 1000;
plts = zeros(30, iter, 3);
tic;

%{
% TODO implement time step determine how to walk concurrently
X.time_step = 2e-4; % Time of each step (ms)
X.step_num = 5e5; % # step
X.particle_num = 1e3; % # particle
%}

for i = 1:200
    tic;
    sim = random_walker_sim(LUT, B, pairs, boundSize, swc{:, :}, 0.5, 0, iter, true);
    sim = sim.eventloop(iter);
    rwpath = sim.rwpath;
    %     rwpath = unique(rwpath, "rows", "stable");

    plts(i, :, :) = rwpath;
    %     h = plot3([rwpath(1:1000:end, 2)], [rwpath(1:1000:end, 1)], [rwpath(1:1000:end, 3)]);
    toc;
    tim = toc;
    fprintf("%12.0f\t%f\n", size(unique(rwpath, "rows"), 1), tim);
end

% axis equal;
toc;

for i = 1:200
    rwpath = plts(i, :, :);
    rwpath = reshape(rwpath, size(rwpath, [2 3]));
    rwpath = unique(rwpath, "rows", "stable");
    h = plot3([rwpath(1:10:end, 2)], [rwpath(1:10:end, 1)], [rwpath(1:10:end, 3)]);
end
