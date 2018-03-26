%% Main Genes experiment 
% Load processed microarray data (see Get_Microarray_Data.m)

load MatFiles/GetMicroData.mat
n = size(G1,1);

Acor = corrcoef(G1'); % get correlation coefficients

    
%% Threshold to form a graph
threshold = .9;
A = Acor > threshold; 
A = A-diag(diag(A));
m = nnz(A)/2;

connected = find(sum(A) > 0);
noconnections = numel((find(sum(A) == 0)));

fprintf('With a threshold of %f, there are %d edges\n',threshold, m);
fprintf('There are %d nodes with zero edges \n',noconnections)
fprintf('The final graph has %d nodes \n',n-noconnections)
fprintf('Average degree is %f among nodes of nonzero degree \n',m/(n-noconnections))

G = A(connected,connected);
SGD = gene_sgdid(connected,:);
ORF = gene_orf(connected,:);


%% Cluster the graph with our algorithm CD4

addpath('../Algs')
addpath('MatFiles')
load T131.mat
[Dist,Bound] = CD_LPrelax(G,T131);

c = CD_lp_round(Dist);
lamCCobj(G,0,c)

save('Genes131opt','c')
%% We can cluster with other algorithms as well

cGC = GrowClique(G,10);
lamCCobj(G,0,cGC)