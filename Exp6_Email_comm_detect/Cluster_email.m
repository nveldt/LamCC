load Email

A = spones(A);
A = A-diag(diag(A));
currdir = pwd;

cd('../../LFR_lcc/Code')
n = size(A,1);
k = 42;
addpath('../Community_BGLL_Matlab/')
addpath('Algs/')
addpath('../Lambda_louvain/')
addpath('../find_densest_subgraph-master/matlab_wrapper_loadlibrary/')
addpath('../metis-5.0.2/metismex/')
addpath('/Users/nateveldt/GitHubRepos/LFR_lcc/quick/quick-source')


[I,J] = find(A);
EdgeList = [I,J];

%%
scores = [];
for i = 1:1
n = size(A,1);

% Louvain
COMTY = cluster_jl_cpp(triu(full(A)),1,1,0,1);
cLou = COMTY.COM{end};


% Infomap
cInfo = Run_InfoMap(A,1);


%for k = 2:100
% Metisk
k = 27;
cMet = metismex('PartGraphRecursive',sparse(A),k);

% Graclus
k = 13;
cGrac = Run_Graclus(A,EdgeList,k);


scores = [scores; ARI(cMet+1,truth) ARI(cGrac+1,truth) ARI(cLou,truth) ARI(cInfo,truth)]

end

%% Recursive maximum quasi-clique

% den = .55;
% minsize = 50;
% cRMQ = recursive_maximum_quasiclique(A,den,minsize);

%% Grow-Cluster and Lambda-Louvain
addpath('/Users/nateveldt/GitHubRepos/LFR_lcc/Lambda_louvain')
% Lams = .1:.005:.14;
% numlam = numel(Lams);
% numlams = numlam;
% GCclus = zeros(n,numlam);
% LLclus = zeros(n,numlam);
% GCs = zeros(numlam,1);
% LLs = zeros(numlam,1);
% limit = 2;
% for ii = 1:numlam
%     cGrow = GrowCluster(A,Lams(ii),1);
%     cLL = lambda_louvain(A,Lams(ii),ones(n,1),2);
%     GCclus(:,ii) = cGrow;
%     LLclus(:,ii) = cLL;
% end

cd(currdir)
%% Compute the lambda-louvain clustering a different way
score = [];
for i = 1:20

addpath('~/GitHubRepos/GenLouvain-2.1/')
Lams = .00005:.00005:.0003;
Lams = 10e-5;
numlam = numel(Lams);
LLclus2 = zeros(n,numlam);
LL2s = zeros(numlam,1);
e = ones(1,n);
d = full(sum(A,1));
limit = 1000;
for ii = 1:numlam
    %cLL = lambda_louvain(A,Lams(ii),ones(n,1),5);
    lam = Lams(ii);
    B = @(i) A(:,i) - lam*d'*d(i);
    B = full(A-lam*d'*d);
    cLL = iterated_genlouvain(B);
    LLclus2(:,ii) = cLL;
    ari = ARI(cLL,truth);
    sns = Get_sns(A,cLL);
    score = [score; ari];
    fprintf('ARI = %.10f, Lam = %.10f \n',ari,lam)
    fprintf('LL: %.10f  %.10f  %.10f \n',min(sns), mean(sns), max(sns))
end

snst = Get_sns(A,truth);
fprintf('Tr: %f  %f  %f \n',min(snst), mean(snst), max(snst))
end
median(score)

%% NMI and ARI computations

Clusterings = [cMet'+1 cGrac+1 cLou' cInfo LLclus2];

numlam = numel(Lams);
numalgs = size(Clusterings,2);
ariMat = zeros(numalgs,1);
nmiMat = zeros(numalgs,1);


%%
for i = 1:numalgs
    ariMat(i) = ARI(Clusterings(:,i),truth);
    nmiMat(i) = NMI(Clusterings(:,i),truth);
end

plot(Lams,ariMat(1:end))
hold on
plot(Lams,nmiMat(1:end))