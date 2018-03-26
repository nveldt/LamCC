%% Solve the GC objective with a real network

addpath('/homes/lveldt/GitHubRepos')
addpath('~/data/MatlabGraphs/')
addpath('caGrQc_exps')

clear

fprintf('Now for ca-GrQc\n')
graph = 'ca-GrQc';

load AcaGrQc.mat

A = StandardizeFully(A);


n = size(A,1);
Lams = [logspace(-5,-1,10), .15:.1:.95];
numlam = numel(Lams);

% Grow-Cluster and Lambda-Louvain
GCs = zeros(numlam,1);
LLs = zeros(numlam,1);
GCclus = zeros(n,numlam);
LLclus = zeros(n,numlam);
for ii = 1:numlam
    cGrow = GrowCluster(A,Lams(ii),1);
    fprintf('All done with GC \n')
    cLL = lambda_louvain(A,Lams(ii),ones(n,1),3);
    GCs(ii) = lamCCobj(A,Lams(ii),cGrow);
    LLs(ii) = lamCCobj(A,Lams(ii),cLL);
    GCclus(:,ii) = cGrow;
    LLclus(:,ii) = cLL;
    fprintf('Done with Lambda = %f \n',Lams(ii));
end

% Trivial Clusterings
Giant = LamRange_objs(A,ones(n,1),Lams);
Singletons = LamRange_objs(A,(1:n)',Lams);

save(strcat('clusterings/LccAlgs_connected_',graph))