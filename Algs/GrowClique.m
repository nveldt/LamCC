function c = GrowClique(A,times)
% GROWCLIQUE: this is equivalent to running the Grow-Cluster algorithm when
% the chosen value of lambda is so high that we are guaranteed to form
% clusters that are so dense they are just cliques.
%
% This code is basically just a wrapper for GrowCliqueStep
% Coded by Nate Veldt on 5-4-2017
%

n = size(A,1);

if nargin < 2
    times = 1;
end
c = zeros(n,1);
    
    Vinds = (1:n)';
    ClusNum = 1;
    while numel(Vinds >0)
        SubA = A(Vinds,Vinds);
        clusLocal = GrowCliqueStep2(SubA,times); % grab a cluster

        % the indices in clusLocal are the indices of the subgraph
        % A(Vinds,Vinds). We must get the original indices of these nodes
        clus = Vinds(clusLocal);
        
        c(clus) = ClusNum;
        ClusNum = ClusNum + 1;
        Vinds = setdiff(Vinds,clus);    % updates Vinds
    end

end