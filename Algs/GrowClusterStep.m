function [clus,indicator] = GrowClusterStep(A,lam,seed)
% CLUS = GROWCLUSTERSTEP(A): This performs one step (finds one cluster) in
% the GrowCluster algorithm:
%
%   1. randomly select a node in the graph with adjacency matrix A
%   2. randomly add one of its edges greedily (which improves CC objective)
%   3. For every node v outside current cluster, compute
%                       del_v = cut(C,v) - lam*|C|
%   4. If there exists del_v > 0, add to cluster C the node v with
%       maximum del_v score. If not, end the algorithm.
%   5. Return to step 3 and repeat
%
%
% Inputs:
%       A = adjacency matrix for a graph
%       lam = parameter which controls how much we add
%       seed = optimal seed parameter, telling you where to begin merging
%               in the graph
%
% Coded by Nate Veldt on 2-03-17

n = size(A,1);

if nargin < 3
    seed = randi(n);    % select a random seed if one is not given
end
    
clus = seed;            % current cluster made of starter/seed/pivot node   

% We want to efficiently keep track of del_v scores

sizeC = numel(seed);    % keep |C| updated
cutC = A(:,clus);       % keep cut(C,v) updated
if size(cutC,2) > 1
    cutC = sum(cutC,2);
end
    
assert(size(cutC,2) == 1);   % this should be a vector


del = cutC - lam*sizeC;
m = max(del);


while m > 0
    canAdd = find(del == m);    % take all possible nodes we could add
    j = randi(numel(canAdd));   % select one at random
    toAdd = canAdd(j);
    
    clus = [clus;toAdd];        % add it to the current cluster
    clus = unique(clus);
    sizeC = sizeC + 1;          % our cluster grew by one
    
    newNeighbors = A(:,toAdd);  % find all neighbors of the new node
                                % these are the nodes whose connection to
                                % the cluster C increased by one;
                                
    cutC = cutC + newNeighbors;
    cutC(clus) = 0;             % make sure we don't add back in nodes that
                                % are already in the cluster
    
    del = cutC - lam*sizeC;
    m = max(del);
end
indicator = zeros(n,1);
indicator(clus) = 1;

end