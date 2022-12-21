clc; clearvars; close all;
tic;
addpath(genpath("./"));
swc = read_t('exampleTree.swc');
NodeID = swc(:, 1); Coords = swc(:, 3:5);
Radii = swc(:, 6); Parents = swc(:, 7);

tree = table(NodeID, Coords(:, 1), Coords(:, 2), Coords(:, 3), ...
    Radii, Parents, 'VariableNames', ...
    {'NodeId', 'X', 'Y', 'Z', 'Radii', 'Parent'});


save("tree.mat", "tree");

clearvars; clc;
tic;
load("tree.mat", "tree");
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

%{
change swc to matrix


%}
%%
close all;
figure();
hold on;
[~, ~] = mainLoop(swc, zeros(boundSize), b, pairs);
axis equal
%%
clf
tic;
clc;
for i = 1:20
    tic;
    sim = random_walker_sim(LUT, A, pairs, boundSize, swc{:,:}, 1, 0);
    sim = sim.eventloop(10000);
    rwpath = sim.rwpath;
    hold on;
    h = plot3(rwpath(:, 2) + 1, rwpath(:, 1) + 1, rwpath(:, 3) + 1);
    if mod(i,5) == 0
        toc;
    end
    fprintf("%12.0f\n",size(unique(rwpath,"rows"),1));
    pause(0.1);
end
toc;

axis equal;