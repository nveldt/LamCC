function [C,BestObj] = GrowCluster(A,lam,times)
% GROWCLUSTER(A,LAM): A variation on the MergeCluster algorithm
% where one cluster is formed at a time and then removed from the graph.
%
% Procedure outline:
% 1. If any nodes still left in graph, select a pivot
% 2. Run GrowClusterStep to greedily grow a cluster around pivot
% 3. Remove the cluster from the current graph and go back to 1.
%
% Run this procedure for "times" instantiations, and output the one with
% the best correlation clustering objective score

n = size(A,1);
BestObj = n^3;

if nargin < 3
    times = 1;
end

for t = 1:times
    
    Vinds = (1:n)';
    C = [];
    Canswer = zeros(n,1);
    ClusNum = 1;
    while numel(Vinds >0)

        [clusLocal,~] = GrowClusterStep(A(Vinds,Vinds),lam); % grab a cluster

        % the indices in clusLocal are the indices of the subgraph
        % A(Vinds,Vinds). We must get the original indices of these nodes

        clus = Vinds(clusLocal);

        Canswer(clus) = ClusNum;
        ClusNum = ClusNum + 1;
        Vinds = setdiff(Vinds,clus);    % updates Vinds
    end

    if times > 1
    Obj = lamCCobj(A,lam,Canswer);
    
            if Obj < BestObj
                BestObj = Obj;
                BestC = Canswer;
            end
    else
        BestC = Canswer;
    end
end

C = Canswer;

end