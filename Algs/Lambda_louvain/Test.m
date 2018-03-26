%% Test Lambda Louvain algorithm
load('../Code/graphfiles/A300lfr.mat')
addpath('../Lambda_louvain/')
n = size(A,1);
lam = .3;
w = ones(n,1);
c = lambda_louvain(A,lam,w);
max(c)