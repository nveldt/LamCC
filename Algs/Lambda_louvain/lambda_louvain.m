function c = lambda_louvain(A,lam,w,itlim)
% c = lambda_louvain(A,lam): A heuristic algorithm for approximating the
% lambdaCC objective, in the style of the Louvain algorithm. The mechanics
% of the algorithm are the same as the traditional Louvain algorithm, but
% instead of moving nodes based on changes in the modularity objective, we
% move around nodes based on changes in the lambdaCC objective.
%
% The algorithm is recursive. A can be weighted.
%
% This can be generalized later, but currently this algorithm is meant to
% approximate the optimal standard LambdaCC objective (not degree-weighted)
%
% Inputs:
%   * A: adjacency matrix for an undirected network
%   * lam: a value between 0 and 1
%   * w: node weights (basically the number of original nodes inside the
%        supernode)
%   * itlim: limit the number of times you go through and move nodes before
%       collapsing them into a smaller graph and moving on.

n = size(A,1);
assert(lam<1)
assert(lam>0)

if nargin < 4
    itlim = n;
end

if nargin < 3
    w = ones(n,1);
end

% Get the updated clustering based on one step 
cnew = ll_step(A,lam,w,itlim);

if max(cnew) == n
    c = cnew;       % do nothing, there was no change in the last step
else
    fprintf('New level \n')
    % compress the graph based on the new clustering
    [A,w] = compress(A,cnew,w);
    
    % Get the full solution for the compressed graph
    cComp = lambda_louvain(A,lam,w,itlim);
    
    % And convert it back to size
    c = cComp(cnew);
end
    
