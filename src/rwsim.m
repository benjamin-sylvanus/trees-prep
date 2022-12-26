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

% embeddings -> adda openai
% least semantic correlation to user query.
% use cos similarity

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

%%

% function [A,poses] = mainLoop(tree,A,pairBounds,pairs)
%         A=A;
%         swc = tree;
%
%         ecount = 0; prevElapse = 0; tic;
%
%         for i = 1:size(pairs,1)
%                 Bound = pairBounds{i,1};
%
%                 Sx = Bound(1,1); Sy = Bound(1,2); Sz = Bound(1,3);
%
%                 Nx = Bound(2,1) + 1; Ny = Bound(2,2) + 1; Nz = Bound(2,3) + 1;
%
%                 [y0,x0,z0]=meshgrid(Sy:Ny,Sx:Nx,Sz:Nz);
% %                 pos_fill = ndSparse.build(range(Bound)+2);
%
%                 p1 = swc{pairs(i,1),2:5};
%
%                 pair = swc(pairs(i,2),:);
%
%                 p2 = pair{1,2:5};
%
%                 x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
%
%                 x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
%
%                 pos = swc2v(x0,y0,z0,x1,x2,y1,y2,z1,z2,r1,r2,Nx,Ny,Nz);
%
%                 pos_fill = A(Sx:Nx,Sy:Ny,Sz:Nz);
%
%                 Q = find(pos==1);
%                 [ix iy iz]=ind2sub(size(pos),Q);
%                 mxyz = mean([iy,ix,iz]);
%
%                 poses(i,:)=mxyz + [Sy Sx Sz];
%
%
%
%                 hold on;
%                 try
%                 is = isosurface(pos,0);
%                 sum(pos,"all");
%                 is.vertices = is.vertices + [Sy Sx Sz];
%                 p = patch('Faces',is.faces,'Vertices',is.vertices);
%                 p.FaceColor='green';
%                 p.FaceAlpha=0.3;
%                 p.EdgeColor='none';
%                 view(3);
%                 axis tight;
%                 catch ME
%                 end
%
%                A(Sx:Nx,Sy:Ny,Sz:Nz) = pos_fill | pos;
%
%             end
%
%
%             A(Sx:Nx,Sy:Ny,Sz:Nz) = pos_fill;
% end

% Replaced Calulation for r comparison

% function pos = swc2v(x0,y0,z0,x1,x2,y1,y2,z1,z2,r1,r2,Nx,Ny,Nz)
%         t = ( (x0-x1)*(x2-x1) + (y0-y1)*(y2-y1) + (z0-z1)*(z2-z1) ) ./...
%             ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
%         x = x1 + (x2-x1)*t;
%         y = y1 + (y2-y1)*t;
%         z = z1 + (z2-z1)*t;
%
%         list1 = (x-x1).*(x-x2) + (y-y1).*(y-y2) + (z-z1).*(z-z2) <0;
%         list2 = ~list1;
%
%         dist2 = (x0-x).^2 + (y0-y).^2 + (z0-z).^2;
%
%
%         %     r = r1 + sqrt((x-x1).^2 + (y-y1).^2 + (z-z1).^2) /...
%         %         sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2) * (r2-r1);
%
%         %     r = ( c + r2 ) / (sqrt ( 1 - ( |r1-r2 | / l ) )
%
%         %     c = ( |r1 - r2| * l ) / L
%
%             if r2 > r1
%                 pos = swc2v(x0,y0,z0,x2,x1,y2,y1,z2,z1,r2,r1,Nx,Ny,Nz);
%             else
%
%             rd = abs( r1 - r2 );
%
%             % distance from orthogonal vector to p2
%             l = sqrt( (x-x2).^2 + (y-y2).^2 + (z-z2).^2 );
%
%             % distance from p1 -> p2
%             L = sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
%
%             c = (rd * l) ./ sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
%             r = (c + r2) ./ sqrt(1 - ( (rd / L) .^2) );
%
%
%             pos1 = dist2<=(r.^2); % smaller in one line and less than and equal
%             pos2 = ( ( (x0-x1).^2 + (y0-y1).^2 + (z0-z1).^2 ) <= (r1^2) ) | ...
%                    ( ( (x0-x2).^2 + (y0-y2).^2 + (z0-z2).^2 ) <= (r2^2) );
%             pos = ndSparse.build(size(x0)); % use false
%
%
%             pos(list1) = pos1(list1);
%             pos(list2) = pos2(list2);
%             end
%     end

% function [pos,pos2] = updateposition(step,pos)
%         vect = random_unit_vector(3,1);
%         delta = vect * step;
%         pos2 = pos + delta;
% end

%% Functions

%  function swc = read_t(filename)
%         fid = fopen(filename);
%         A = textscan (fid, '%s', 'delimiter', '\n');
%         A = A{1};
%         fclose (fid);
%         swc = [];
%         for counter  = 1 : length (A)
%             if ~isempty (A{counter})  % allow empty lines in between
%                 % allow comments: lines starting with #:
%                 if ~strcmp (A{counter} (1), '#')
%                     swc0   = sscanf (A{counter}, '%f')';
%                     swc    = [swc; swc0];
%                 end
%             end
%         end
% end

% function d = calcdists(tree,i)
%
%         cxyz = tree{i,["X","Y","Z"]};
%         parentId = tree{i,"Parent"};
%         pxyz = tree{parentId,["X","Y","Z"]};
%         coord_1 = cxyz;
%         coord_2 = pxyz;
%         difference = abs((coord_2-coord_1).^2);
%         d_squared = sum(difference);
%         d = d_squared;
%         d = sqrt(d_squared);
% end

% function [b,swc,boundSize,pairs,VSIZE] = initbounds(tree,dists,voxelScale)
%
%         swc = tree;
%         swc(1,"Parent") = {1};
%         X = swc.X+2; Y = swc.Y+2; Z = swc.Z+2;
%         dim = [X Y Z];
%
%         % Threshold Radii
%         [~,~,bin] = histcounts(swc.Radii, 5);
%         avg = mean(swc.Radii(bin==1));
%         swc.Radii(swc.Radii < avg)=avg;
%         r = swc.Radii;
%         ranges = [min(dim - r); max(dim+r)]; % (min - ri) (max + ri)
%
%         % Define Voxel Size
%         VSIZE =  voxelScale*min(dists(2:end));
%
%         % Translate X,Y,Z
%         tran = ((dim-ranges(1,:))./VSIZE)+1;
%
%         % Update Bounds
%         boundSize=ceil((ranges(2,:) - ranges(1,:))./VSIZE)+2;
%
%         % Radii Scale
%         R = swc.Radii./(VSIZE);
%
%         % Set Tree XYZ to Translated
%         % min(swc.tree.Z-swc.tree.Radii);
%         swc.X = tran(:,1); swc.Y = tran(:,2); swc.Z = tran(:,3);
%         swc.Radii = R;
%         % min(swc.tree.Z-swc.tree.Radii)
%
%         [b,pairs] = calcBounds(swc,VSIZE);
%         addpath(genpath(['/Users/benjaminsylvanus/Documents/GitHub/' ...
%                 'SparseMatrixGenerator/ndSparse_G4_2021_03_16']));
%         b = b.a;
% end

% Input [r1;r2] [x1 y1 z1; x2 y2 z2]
% Output [min(ri-ci) max(ri-ci)]

% function [b,pairs] = calcBounds(swc,vs)
%         b = struct();
%         t = swc;
%         mr = vs;
%         pairs = [];
%         for i = 1:numel(t.X)
%                 pair = swc(swc{i,"Parent"},:);
%                 pairs(i,1) = t{i,"NodeId"};
%                 pairs(i,2) = pair{1,"NodeId"};
%                 if (pairs(i,2)==0)
%                         pairs(i,2) = pairs(i,1);
%                 end
%                 b.a(i,:) = pairBounds(t(i,:),pair,vs);
%         end
%
%         function b = pairBounds(p1, p2, m)
%                     % ri: pair radii
%                     % pxi: pair coords
%
%                     ri = [p1.Radii; p2.Radii];
%                     pxi = [p1{1,2:4}; p2{1,2:4}];
%                     b = {[round(min(pxi-ri)); ceil(max(pxi+ri))]};
%                     ba = [min(pi-(ri+m)), max(pi+(ri+m))];
%         end
% end

% function [A,indicies,t2,LUT] = generateLUT(B,b)
%         A = cell(0,0);
%         tic;
%         LUT =zeros(B);
%         [cx,cy,cz,ci]  = cellfun(@extract_Range2, b,'UniformOutput',false);
%
%         bsx = repmat(B(1),size(b));
%         bsy = repmat(B(2),size(b));
%         bsz = repmat(B(3),size(b));
%
%         bsize = [bsx(:),bsy(:),bsz(:)];
%         bsize = num2cell(bsize,2);
%         inds = ci;
%         inds = cellfun(@helper2,bsize,ci,'UniformOutput',false);
%
%         indicies = 0;
%         t2 = 0;
%         la=0;
%
%         for i = 1:length(inds)
%                 idx = [inds{i,1}];
%                 sub = LUT(idx);
%                 idy = find(sub < 1);
%                 j = [1:length(idy)]';
%
%                 nextid = la+j;
%
%                 lj = length(j);
%
%                 if length(j) > 0
%                         while true
%                                 la = nextid(length(j));
%                                 if la > size(A,1)
%                                         A = [A;cell(100000,1)];
%                                 end
%                                 if la < size(A,1)
%                                         break;
%                                 end
%                         end
%                 end
% %                 A =  [A;cell(lj,1)];
%                 sub(idy(j)) = deal(nextid);
%
% %                 for j = 1:length(idy)
% %                         nextid = size(A,1)+1;
% %                         A = [A;{[]}];
% %                         sub(idy(j))=nextid;
% %                 end
%                 idn = find(sub>0);
%                 for j = 1:length(idx)
%                         A(sub(j),1) = {[A{sub(j),1};i]};
%                 end
%                 LUT(idx)=sub;
%         end
% end

% function [cx,cy,cz,ci] =extract_Range2(sub)
%         xi = sub(1,1);
%         yi = sub(1,2);
%         zi = sub(1,3);
%         xe = sub(2,1);
%         ye = sub(2,2);
%         ze = sub(2,3);
%         rx = xi:xe;
%         ry = yi:ye;
%         rz = zi:ze;
%
%         r = {rx,ry,rz};
%
%        [cx,cy,cz,ci] = distributeelem(rx,ry,rz);
%         function [cx,cy,cz,ci] = distributeelem(rx,ry,rz)
%                 csize = prod([size(rx,2),size(ry,2),size(rz,2)]);
%                 cx = repmat(rx',length(ry),1);
%                 cy = repelem(ry,length(rx));
%                 cx_cy = [cx, cy'];
%                 cxy = repmat(cx_cy,length(rz),1);
%                 cz = repelem(rz,size(cx_cy,1));
%                 scxy = [size(cxy)];
%                 scz = [size(cz)];
%                 ci = [cxy,cz'];
%         end
% end

% function temp = helper(temp,inds,i)
%
%         k = repelem(i,size(temp,1));
%         temp(:,end+1)={i};
%         temp
% end

% function [inds] = helper2 (bsize,ci)
% try
%         inds = sub2ind(bsize,ci(:,1),ci(:,2),ci(:,3));
% catch
% end
% end

% %%
% rw = randomwalker(1,1,boundSize,swc, LUT, A, pairs);
%
% % to verify inside visually
% inds = pairs(A{LUT(rw.curr(1), rw.curr(2),rw.curr(3))},:);
%
% % pairlist => [child, parent]
% children = inds(:, 1);
% parents = inds(:, 2);
% [X, Y, Z]= sphere(16);
% clf;
% hold on;
%
% % for each pair: check if point is inside
% for i = 1:length(children)
%
%     % get base and target ids
%     baseid = children(i); targetid = parents(i);
%
%     % p1 = swc{pairs(i, 1), 2:5};
%     p1 = swc{baseid, 2:5}; p2 = swc{targetid, 2:5};
%
%     x1 = p1(1); y1 = p1(2); z1 = p1(3); r1 = p1(4);
%     X2 = X * r1; Y2 = Y * r1; Z2 = Z * r1;
%     h = surf(X2 + x1,Y2 + y1,Z2 + z1);
%     h.FaceAlpha=0.05;
%     h.EdgeColor="none";
%     h.FaceColor = 'red';
%
%     x2 = p2(1); y2 = p2(2); z2 = p2(3); r2 = p2(4);
%     X2 = X * r2; Y2 = Y * r2; Z2 = Z * r2;
%     h = surf(X2 + x2,Y2 + y2,Z2 + z2);
%     h.FaceAlpha = 0.05;
%     h.EdgeColor="none";
%     h.FaceColor = 'red';
% end
% X2 = X * 0.1; Y2 = Y * 0.1; Z2 = Z * 0.1;
% h = surf(X2 + rw.curr(1),Y2 + rw.curr(2),Z2 + rw.curr(3));
% h.FaceColor = "blue";
