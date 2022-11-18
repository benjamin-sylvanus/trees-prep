addpath("/Users/benjaminsylvanus/Documents/CoopMRI/treestoolbox-master/")

d =  '/Users/benjaminsylvanus/Downloads/SWC';

start_trees; clc;

cd( '/Users/benjaminsylvanus/Downloads/SWC');

dswc = dir("*.swc");
names = {dswc.name}'
load_tree(names{1});
t = {};
for i = 1:length(dswc)
    t{i,1} = load_tree(names{i  },' ');
end
cd ('/Users/benjaminsylvanus/Desktop/layer-analysis/')

%%
close all;
figure();
hold on;
for i = 1:100

    hold on;
    plot_tree(t{i,1});

end


%%
lens = len_tree(t{1,1});
tic;
rt = resample_tree(t{1,1},100,'-s -d');        
toc; 
% plot_tree(t{1,1});
%%
figure();
hist(lens,5);
figure;
hist(len_tree(rt),5);

%%
close all;
for i = 1:5
    tic; 
    figure();
    rt = resample_tree(t{i,1},100,'-s -d');     
    toc; 
    figure();
    title(names{i});
    figure();
    s1 = stats_tree(t{i,1});
    figure();
    s2 = stats_tree(rt);
end

