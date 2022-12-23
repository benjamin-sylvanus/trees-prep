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
[b, swc, boundSize, pairs, VSIZE] = initbounds(tree, dists, 3);

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
axis equal;
view(3)

hold on;
tic;
clc;
iter = 5000;
for i = 1:5
    tic;
    sim = random_walker_sim(LUT, A, pairs, boundSize, swc{:,:}, 0.1, 0,iter);
    sim = sim.eventloop(iter);
    rwpath = sim.rwpath;
    rwpath = unique(rwpath,"rows","stable");
    h = plot3(rwpath(:, 2) , rwpath(:, 1) , rwpath(:, 3) );
    tim = toc;
    fprintf("%12.0f\t%f\n",size(unique(rwpath,"rows"),1),tim);
%     pause(0.1);
end
toc;

axis equal;
