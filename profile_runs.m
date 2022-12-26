clc; clearvars; close all;
tic;
addpath(genpath("../../treestoolbox"));
addpath(genpath("./"));
addpath(" / Users / bensylvanus / Library / Application Support / MathWorks / MATLAB Add - Ons / " + ...
"Collections/random unit vector generator");

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
[b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, 0.45);

toc;
tic;
[A, indicies, t2, LUT] = generateLUT(boundSize, b);
toc;

A = A(~cellfun('isempty', A));
sizes = cellfun('size', A, 1);

% get pairs from A
B = cell(size(A,1),2);
for i = 1:length(B)
    pairlist = A{i};
    ps = pairs(pairlist, :);
    
    % pairlist => [child, parent]
    children = ps(:, 1);
    parents = ps(:, 2);
    
    B{i,1}=children;
    B{i,2}=parents;
end

%{
change swc to matrix

%}
%%
% TODO only reasonable improvement would be vectorized rw to avoid for loop
% TODO and to convert sub2ind calculations to custom function
close all;
figure();
hold on;
[~, ~] = mainLoop(swc, zeros(boundSize), b, pairs);
axis equal;
view(3)
%%
hold on;
tic;
clc;
iter = 10000;

for i = 1:10
    tic;
    sim = random_walker_sim(LUT, B, pairs, boundSize, swc{:, :}, 1, 0, iter,true);
    sim = sim.eventloop(iter);
    rwpath = sim.rwpath;
    rwpath = unique(rwpath, "rows", "stable");

%     h = plot3(rwpath(:, 2), rwpath(:, 1), rwpath(:, 3));
    h = plot3([rwpath(1:100:end, 2)], [rwpath(1:100:end, 1)], [rwpath(1:100:end, 3)]);

    tim = toc;
    fprintf("%12.0f\t%f\n", size(unique(rwpath, "rows"), 1), tim);
    %     pause(0.1);
end

toc;




