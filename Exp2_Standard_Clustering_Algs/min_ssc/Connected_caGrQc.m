addpath('../Community_BGLL_Matlab/')
addpath('Algs/')
addpath('../Lambda_louvain/')
addpath('../find_densest_subgraph-master/matlab_wrapper_loadlibrary/')
addpath('../metis-5.0.2/metismex/')
addpath('caGrQc_exps')
%%
Lams = [logspace(-5,-1,10), .15:.1:.95];
numlam = numel(Lams);
numlams = numlam;
numtimes = 1;
Louobj = zeros(numlams,numtimes);
Metisobj = Louobj;
Graclusobj = Louobj;
Denobj = Louobj;
Infoobj = Louobj;

i = 1;
addpath('~/data/MatlabGraphs/')

%% caGcQc
graph = 'ca-GrQc';
load(graph);
A = Problem.A;
A = A-diag(diag(A));
A = StandardizeFully(A);
n = size(A,1);
[I,J] = find(A);
EdgeList = [I J];
Kvals = [2:10 15:3:30 50 100 150 200 300 400 500 750 1000 1250 1500 1750 2000 2250 2500];


% Louvain
COMTY = cluster_jl_cpp(triu(full(A)),1,1,0,1);
cLou = COMTY.COM{end};
Louobj(:,i) = LamRange_objs(A,cLou,Lams);

% Infomap
cInfo = Run_InfoMap(A,1);
Infoobj(:,i) = LamRange_objs(A,cInfo,Lams);

% Recursive densest subgraph
cDen = meng_recursive_dense_subgraph(A);
Denobj(:,i) = LamRange_objs(A,cDen,Lams);

% Metis
cMet = zeros(n,numel(Kvals));
for ii = 1:numel(Kvals)
    kk = Kvals(ii);
    c = metismex('PartGraphRecursive',sparse(A),kk);
    cMet(:,ii) = c+1;
end

AllObjsMet = zeros(numel(Kvals),numlam);
for ii = 1:numel(Kvals);
    AllObjsMet(ii,:) = LamRange_objs(A,cMet(:,ii),Lams);
end
[Metisobj(:,i), metisLamclus] = min(AllObjsMet);

% Graclus
cGrac = Run_Graclus(A,EdgeList,Kvals);
AllObjsGrac = zeros(numel(Kvals),numlam);
for ii = 1:numel(Kvals);
    AllObjsGrac(ii,:) = LamRange_objs(A,cGrac(:,ii),Lams);
end
[Graclusobj(:,i), gracLamclus] = min(AllObjsGrac);

save(strcat('clusterings/Nonlcc_connected_',graph))

