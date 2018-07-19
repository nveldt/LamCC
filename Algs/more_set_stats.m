function [snc,ssc,cut,vol,edges,cond] = more_set_stats(A,S,volA)

if nargin < 3
    volA = sum(nonzeros(A));
end

if numel(S) == size(A,1)
    % then we have an indicator vector
    S = find(S);
    AS = A(:,S);
else
    % then we have a subset
    assert(min(S) >= 1)
    assert(max(S) <= size(A,1))
    AS = A(:,S);
end

vol = full(sum(nonzeros(AS)));
edges = full(sum(nonzeros(AS(S,:))));
cut = vol-edges;

cond = cut/min(vol,volA-vol);

snc = cut/(vol*(volA-vol));

n = size(A,1);
s = numel(S);
ssc = cut/(s*(n-s));

end