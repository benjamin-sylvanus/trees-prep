%%
zz_test_tmp;

% addpath(genpath('/Users/benjaminsylvanus/Desktop/random_unit_vector'));
%%
clc; clearvars; close all;
tic;
addpath(genpath("./"));
swc = read_t('exampleTree.swc');
NodeID = swc(:, 1); Coords = swc(:, 3:5);
Radii = swc(:, 6); Parents = swc(:, 7);

tree = table(NodeID, Coords(:, 1), Coords(:, 2), Coords(:, 3), ...
    Radii, Parents, 'VariableNames', ...
    {'NodeId', 'X', 'Y', 'Z', 'Radii', 'Parent'});

g = graph(NodeID(2:end), Parents(2:end), 1, string(1:length(Parents)));

%%
sg = shortestpathtree(g, [243, 252, 257, 260, 267], 267);
swc(unique(str2double(sg.Edges.EndNodes)));
save("tree.mat", "tree");
unique(str2double(sg.Edges.EndNodes));
sg1 = subgraph(g, unique(str2double(sg.Edges.EndNodes)));
newnames = string(1:height(sg1.Nodes))';
sg1.Nodes.Name = newnames;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
tic;
load("tree.mat", "tree");
dists = zeros(height(tree), 1);

for i = 2:height(tree)
    dists(i, 1) = calcdists(tree, i);
end

toc;

tic;
[b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, 1);

toc;
% addpath(genpath(['/Users/benjaminsylvanus/Documents/GitHub/' ...
%         'SparseMatrixGenerator/ndSparse_G4_2021_03_16']));
tic;
[A, indicies, t2, LUT] = generateLUT(boundSize, b);
toc;

A = A(~cellfun('isempty', A));
sizes = cellfun('size', A, 1);

%%

%%%%%% IF THE RW CROSSES OUTSIDE BOUNDS OF
% LUT THE CURRPOINT STILL NEEDS TO BE CHECKED

close all;
figure();
hold on;
[~, ~] = mainLoop(swc, zeros(boundSize), b, pairs);
axis equal
%%
tic;

for i = 1:2
    tic;
    sim = random_walker_sim(LUT, A, pairs, boundSize, swc, 1, 0);
    sim = sim.eventloop(10000);
    rwpath = sim.rwpath;
    hold on;
    h = plot3(rwpath(:, 2) + 1, rwpath(:, 1) + 1, rwpath(:, 3) + 1);
    h.Color = 'k';
    toc;
end

toc;

axis equal;

%% To PLOT
% close all;
figure();
hold on;
[C, poses] = mainLoop(swc, zeros(boundSize), b, pairs);
%%
step = VSIZE / 1000;
% % Direction is homogenous sampling of cos(theta) phi
% % https://www.ni.com/docs/en-US/bundle/labview/page/gmath/3d_cartesian_coordinate_rotation_euler.html
iter = 1000;
ps = 1;

for k = 1:10

    if ps > size(b, 1)
        ps = 1;
    end

    p = zeros(iter, 3);
    prng = [{14:18}; {8:11}; {14:17}];
    pos = poses(ps, :)';

    for i = 1:iter
        [pos, updated] = updateposition(step, pos);
        p(i, :) = pos;
        pos = updated;
    end

    h = plot3(p(:, 1) + 1, p(:, 2) + 1, p(:, 3) + 1);
    h.Color = [rand; rand; rand];
    ps = ps + 1;
end
