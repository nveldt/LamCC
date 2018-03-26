function [Obj,noncliques] = cluster_deletion_objective(A,c)
% Given a sparse graph A and a clustering vector c, this computes the total
% number of edges cut by the clustering in c
% This also checks to make sure the clusters are cliques.
%
% In general, A will be sparse.
%
% Coded by Nate Veldt on June 14, 2017


n = size(A,1);
A = A-diag(diag(A));

assert(issparse(A) == 1)
assert(size(c,2) == 1)

Obj = 0;
noncliques = [];
% Get the number of clusters
numClus = max(c);

for i = 1:numClus
    
    % Get the cluster
    S = find(c == i);
    AS = A(:,S);
    % Check that it's a clique
    b = numel(S);
    B = AS(S,:);
    nnzB = nnz(B);
    if nnz(B) ~= b^2-b
        fprintf('Cluster %d is not a clique \n',i);
        nonclique = [nonclique; i];
    end
    % Find the number of edges that are cut
    vol = full(sum(nonzeros(AS)));
    edges = nnzB;
    cut = vol-edges;
    
    Obj = Obj + cut;
end
Obj = Obj/2;

