addpath('Algs')
addpath('Algs/Lambda_louvain/')

load Networks/karate.mat

% Get adjacency matrix of a graph and set lambda
A = Problem.A;
lam = .6;
n = size(A,1);

%% Exactly solve Lambda-CC using Gurobi Software
% Requires installation of Gurobi. 
% Will only work for small graphs (constraint matrix with O(n^3) rows is
% explicitly generated)

[C,cOpt,lcc_objective] = lam_cc(A,lam);

fprintf('Graph has %d nodes, lam = %f. \nOptimal LamCC clustering has %d clusters and an objective score of %f \n',size(A,1),lam,max(cOpt),lcc_objective)

%% Run GrowCluster and GrowClique

cGC = GrowCluster(A,lam);
GCobj = lamCCobj(A,lam,cGC);
fprintf('GrowCluster forms %d clusters, objective is %f \n',max(cGC),GCobj) 


% k = number of random cliques GrowClique grows at each step
k = 2;
cCl = GrowClique(A,k);
Clobj = cluster_deletion_objective(A,cCl);
fprintf('GrowClique forms %d cliques. Cluster Deletion score is %d \n',max(cCl),Clobj)

%% Run Lambda-Louvain (Standard LambdaCC, not degree-weighted)

% itlim: the number of times lambda_louvain goes through all the nodes and
% performs locally optimal moves before collapsing into supernodes.
itlim = 10;
cLL = lambda_louvain(A,lam,ones(n,1),10);
LLobj = lamCCobj(A,lam,cLL);
fprintf('Lambda-Louvain forms %d clusters, objective is %f \n',max(cLL),LLobj) 

%% Run LP-based algorithm, ThreeLP

% Get LP-relaxation
[D, LPbound] = lam_cc_LPrelax(A,lam);  

% Apply rounding procedure t times and output best result
t = 1;
[cLP3,LP3obj] = ThreeLP_round(A,D,t);
fprintf('ThreeLP forms %d clusters, objective is %f \n',max(cLP3),LP3obj) 

%% Run TwoCD

[Dist,Bound] = CD_LPrelax(A);
cCD2 = CD_lp_round(Dist);
CD2obj = cluster_deletion_objective(A,cCD2);

fprintf('TwoCD forms %d cliques, Cluster Deletion score is %d \n',max(cCD2),CD2obj)
