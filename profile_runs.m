%% Section 1
clc; clearvars; close all;
scale = 1;

[LUT, B, pairs, boundSize, swc, memoized_distance,A,vsize,ranges,b] = ...
    prepSim(0.3,"democell.swc",1); 

clc;
% dependent on geometry and diffusion time
D0 = 2; % 2um^2/ms; 
d = (min(swc.Radii).*0.05)*vsize; % swc is in voxel scale 
tstep = ((d*d)/6)/D0;
step_size = d./vsize; 
perm_prob = 0.0000;
init_in = true;

% TODO Enable diff units of vsize and add conversion in this units table
% TODO Apply units to results function 
C = {"D0","d","tstep","step_size", "vsize";
    D0,d,tstep,step_size,vsize;
    str2symunit('mcm^2/ms'), str2symunit('um'), str2symunit('ms'),"1",str2symunit('mcm')};
S = cell2sym(C); V = S(1:2,1:5);
V(2,1) = S(2,1)* S(3,1); V(2,2) = S(2,2)* S(3,2);
V(2,3) = S(2,3)* S(3,3); V(2,4) = S(2,4)* S(3,4);
V(2,5) = S(2,5)* S(3,5); clear S C; V


% function to reverse scaling & translation done in prepsim.
reverse_scale = @(swc,vs,original_range) original_range(1,:)+(vs.*(swc-2));
reversed =reverse_scale(swc{:,2:4},vsize,ranges);

% function to reverse just scaling (voxel->init)
prior_voxelized = @(swc,vs) vs*(swc{:,2:5});
pv = prior_voxelized(swc,vsize);

swcmat = swc{:, :};
clear sim;

% at least 1e3 steps
% at least 1e4 particles
step_num = 1e4; 
particle_num = 1e4;

% points = zeros(length(A),1);
% uniqs{1,1} = A{1};
% for i = 1:length(A)
%     ele = [A{i}];
%     b = false;
%     for j = 1:length(uniqs)
%          if numel(uniqs{j}) == numel(ele)
%             if (all(uniqs{j}==ele))
%                 b=true;
%                 points(i) = j;
%             end
%         end
%     end
%     if ~b
%         uniqs{j+1,1} = ele;
%         points(i) = j+1;
%     end
% end

%% Section 2
% initialize sim
tic; 
sim = random_walker_sim(LUT, B, pairs, boundSize, swcmat, step_size, ...
    perm_prob, step_num, init_in, particle_num, memoized_distance(:, 3));
toc; 

tic;
    s = swc.Radii(1)-step_size;
    x0 = swc.X(1); y0 = swc.X(1); z0 = swc.Z(1);
    ps = sim.particles;
    
for i = 1:particle_num
    n = 1; m = 1;
    theta = 2*pi*rand(n,m);
    v = rand(n,m);
    phi = acos((2.*v)-1);
    r = rand(n,m).^(1/3);
    x = x0 +  s*r .* sin(phi) .* cos(theta);
    y = y0 +  s*r .* sin(phi) .* sin(theta);
    z = z0 +  s*r .* cos(phi);
    ps{i}.curr = [x,y,z]';
end

sim.particles = ps;
toc;
% save("randominits","sim");
% clear sim;
% load("randominits");

%% Section 3
tic;
sim = eventloop(sim, step_num);
toc;


%% Section 4

% fid=fopen('/Users/benjaminsylvanus/CLionProjects/testclion/data/LUT.bin','w');
% fwrite(fid,lut,"uint16","b");
path = sim.path;
% path = "./simulations/sim19";
% clearvars -except path vsize tstep
% ;
setresults(path,vsize,tstep);

% refline((swc.Radii(1).*vsize)^5/2);
clear data;
load(sim.path + "/initialPos.mat");
figure();
scatter3(data(1:1000:end,1),data(1:1000:end,2),data(1:1000:end,3));
axis equal;
hold on;
mainLoop(swc{:,:},LUT,b,pairs);

%%
path = sim.path;
load(path + "/DATA.mat");
dataunique = cell(size(data, 2), 1);


for i = 1:size(data, 2)
    q = data(:, i, :);
    q = reshape(q, size(data, [1, 3]));
    size(q);
    dataunique{i} = unique(q, "rows", "stable");
end

f = figure("Name", "Particle Paths");
hold on;
axis equal;
step = 100;
step_part = 100;
for i = 1:step_part:length(dataunique)
    rwpath = dataunique{i};
    h = plot3([rwpath(1:step:end, 2)], ...
        [rwpath(1:step:end, 1)], ...
        [rwpath(1:step:end, 3)]);
end

% 
% for i = 1:size(data, 2)
%     q = data(:, i, :);
%     q = reshape(q, size(data, [1, 3]));
%     size(q);
%     dataunique{i} = unique(q, "rows", "stable");
% end
% 
% f = figure("Name", "Particle Paths");
% hold on;
% axis equal;
% step_iter = 2;
% step_part = 1000;
% for i = 1:step_part:length(dataunique)
%     rwpath = dataunique{3:5};
%     list = rwpath(1:step_iter:end,:);
% 
%     for j = 1:length(list)-1
%         pause(0.001)
%         curr = list(j,:);
%         next = list(j+1,:);
%         h = plot3([curr(1),next(1)], ...
%         [curr(2),next(2)], ...
%         [curr(3),next(3)]);
%     end
% %         h = plot3([rwpath(1:step_iter:end, 2)], ...
% %         [rwpath(1:step_iter:end, 1)], ...
% %         [rwpath(1:step_iter:end, 3)]);
% %     for j = 1:step_iter-1
% %             pause(0.01)
% %     h = plot3([rwpath(j*step_iter:j*step_iter+1, 2)], ...
% %         [rwpath(j*step_iter:j*step_iter+1, 1)], ...
% %         [rwpath(j*step_iter:j*step_iter+1, 3)]);
% %     end
% end
%%
% theres a big problem with large simulation sizes:
% parfor disables writing directly to matobj
% large sims exceed memory (16G) 
% ~ Potential Solution
% Could chunk with for loop 
