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
sizesB = cellfun("size",B,1);

% extract node id [child, parent] -> linear index in B
% B(i) = {[children],[parents]}
memoized_distance = zeros(length(pairs), 3);

for i = 1:length(pairs)
    % extract child and parent ids
    baseid = pairs(i, 1); targetid = pairs(i, 2);

    p1 = swc{baseid, 2:5};
    p2 = swc{targetid, 2:5};

    x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
    x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);

    memoized_distance(i, 3) = (x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2;
    memoized_distance(i, 2) = targetid;
    memoized_distance(i, 1) = baseid;
end

%%
% TODO only reasonable improvement would be vectorized rw to avoid for loop
% TODO and to convert sub2ind calculations to custom function
swcmat = swc{:, :};
%%
close all;
figure();
hold on;
[~] = mainLoop(swcmat, zeros(boundSize), b, pairs);
% axis equal;
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

swcmat = swc{:, :};

clc;
% at least 1e3 steps
% at least 1e4 particles
step_num = 1e3;
particle_num = 1e4;
% dependent on geometry and diffusion time
step_size = 1;
perm_prob = 0;
init_in = true;
tic;

sim = random_walker_sim(LUT, B, pairs, boundSize, swcmat, step_size, ...
    perm_prob, step_num, init_in, particle_num, memoized_distance(:, 3));

sim = sim.eventloop2(step_num);

toc;

% rwpath = sim.rwpath;
% % hold on;
% %
% % for i = 1:sim.particle_num
% %     h = plot3([rwpath(:, i, 2)], [rwpath(:, i, 1)], [rwpath(:, i, 3)]);
% % end
% %
% % f = @() sim.eventloop2(step_num);
% % t = timeit(f);
%%
clearvars;
load("simulations/sim1/matFileOfPositions.mat");

dataunique = cell(size(data, 2), 1);

for i = 1:size(data, 2)
    q = data(:, i, :);
    q = reshape(q, size(data, [1, 3]));
    size(q);
    dataunique{i} = unique(q, "rows", "stable");
end

% data = unique(data,"stable");

hold on;
axis equal;
step = 5;

for i = 1:length(dataunique)
    rwpath = dataunique{i};
    h = plot3([rwpath(1:step:end, 2)], ...
        [rwpath(1:step:end, 1)], ...
        [rwpath(1:step:end, 3)]);
end
