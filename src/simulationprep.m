clearvars;
addpath(genpath('/Users/bensylvanus/Desktop/code/work/'));
d = '/Users/bensylvanus/Desktop/code/work/SWC';

start_trees; clc;

cd('/Users/bensylvanus/Desktop/code/work/SWC');

dswc = dir("*.swc");
names = {dswc.name}';
load_tree(names{1});
t = {};

for i = 1:length(dswc)
    t{i, 1} = load_tree(names{i}, ' '); %#ok<SAGROW> 
end

cd ('/Users/bensylvanus/Desktop/code/work/trees-prep/')

%%

clearvars; start_trees;

load("t.mat");

%%
close all;
rt = resample_tree(t{33}, 2000, '-s');
% rt.D(rt.D>4000) = 4000;
figure
hist(rt.D); 
%%
swc_tree(rt, 'exampleTree.swc');
%%
lens = len_tree(t{1, 1});
tic;
rt = resample_tree(t{1, 1}, 100, '-s -d');
toc;
% plot_tree(t{1,1});
%%
figure();
hist(lens, 5);
figure;
hist(len_tree(rt), 5);
