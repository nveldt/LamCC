function c = renumber(c)
% Take in a cluster with labels that are somewhat all over the place, and
% renumber them so that max(c) = number of unique clusters

map = unique(c);
numclus = numel(map);
n = numel(c);

% Map old cluster names to new cluster names
% Choose convention that the first element is alwasy in cluster 1
for i = 1:n
    newClus = find(map == c(i));
    c(i) = newClus;    
end


end