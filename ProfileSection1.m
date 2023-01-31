%% PROFILE SECTION 1

clc; clearvars; close all;
[LUT, B, pairs, boundSize, swc, memoized_distance, A, vsize, ranges] = ...
    prepSim(0.08, "democell.swc", 1);

% dependent on geometry and diffusion time
D0 = 2; % 2um^2/ms;
d = (min(swc.Radii) .* 0.1) * vsize; % swc is in voxel scale
tstep = ((d * d) / 6) / D0;
step_size = d ./ vsize;
perm_prob = 0.0000;
init_in = true;

% TODO Enable diff units of vsize and add conversion in this units table
% TODO Apply units to results function
C = {"D0", "d", "tstep", "step_size", "vsize";
     D0, d, tstep, step_size, vsize;
     str2symunit('mcm^2/ms'), str2symunit('um'), str2symunit('ms'), "1", str2symunit('mcm')};
S = cell2sym(C); V = S(1:2, 1:5);
V(2, 1) = S(2, 1) * S(3, 1); V(2, 2) = S(2, 2) * S(3, 2);
V(2, 3) = S(2, 3) * S(3, 3); V(2, 4) = S(2, 4) * S(3, 4);
V(2, 5) = S(2, 5) * S(3, 5); clear S C; V

% function to reverse scaling & translation done in prepsim.
reverse_scale = @(swc, vs, original_range) original_range(1, :) + (vs .* (swc - 2));
reversed = reverse_scale(swc{:, 2:4}, vsize, ranges);

% function to reverse just scaling (voxel->init)
prior_voxelized = @(swc, vs) vs * (swc{:, 2:5});
pv = prior_voxelized(swc, vsize);

swcmat = swc{:, :};
clear sim;

% at least 1e3 steps
% at least 1e4 particles
step_num = 1e4;
particle_num = 1e3;

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
