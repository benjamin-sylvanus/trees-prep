
clc; clearvars; close all;

[LUT, B, pairs, boundSize, swc, memoized_distance,A,vsize,ranges] = prepSim(0.1);

reverse_scale = @(swc,vs,original_range) original_range(1,:)+(vs.*(swc-2));

reversed =reverse_scale(swc{:,2:4},vsize,ranges);
% just voxel scale 

prior_voxelized = @(swc,vs) vs*(swc{:,2:5});
pv = prior_voxelized(swc,vsize);

swcmat = swc{:, :};
clear sim;
clc;
% at least 1e3 steps
% at least 1e4 particles
step_num = 1e3;
particle_num = 1e5;
% dependent on geometry and diffusion time
D0 = 2; % 2um^2/ms;
d = min(swc.Radii)*0.5;
tstep = d/6/D0;
step_size = d;
perm_prob = 0.0000;
init_in = true;

points = zeros(length(A),1);
uniqs{1,1} = A{1};
for i = 1:length(A)
    ele = [A{i}];
    b = false;
    for j = 1:length(uniqs)
         if numel(uniqs{j}) == numel(ele)
            if (all(uniqs{j}==ele))
                b=true;
                points(i) = j;
            end
        end
    end
    if ~b
        uniqs{j+1,1} = ele;
        points(i) = j+1;
    end
end


%%

tic;
sim = random_walker_sim(LUT, B, pairs, boundSize, swcmat, step_size, ...
    perm_prob, step_num, init_in, particle_num, memoized_distance(:, 3));
toc;
%%
tic;
sim = eventloop(sim, step_num);
toc;
%%
% fid=fopen('/Users/benjaminsylvanus/CLionProjects/testclion/data/LUT.bin','w');
% fwrite(fid,lut,"uint16","b");
path = sim.path;
% path = "./simulations/sim11";
clearvars -except path vsize tstep;
setresults(path,vsize,tstep);

%%
% theres a big problem with large simulation sizes:
% parfor disables writing directly to matobj
% large sims exceed memory (16G) 
% ~ Potential Solution
% Could chunk with for loop 
