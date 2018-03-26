%% Running the new lambda louvain algorithm
addpath('../Lambda_louvain/')
graph = 'scom10k15';
load(strcat('graphfiles/',graph))

Lams = [logspace(-6,-1,10), .15:.1:.95];

n = size(A,1);
numLam = numel(Lams);
w = ones(n,1);
LLs = zeros(numLam,1);

FullLLlist = zeros(n,numLam);
for i = 1:numLam
    tic
    FullLLlist(:,i) = lambda_louvain(A,lam,w);
    LLs(i) = lamCCobj(A,Lams(i),FullLLlist(:,i));
    fprintf('Lambda Louvain with lam = %f took %f seconds \n',Lams(i),toc)
    save(strcat('clusterings/LL_PseudoPlot_',graph),'Lams','LLs','FullLLlist','n')
end
