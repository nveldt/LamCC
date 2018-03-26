%% Load a graph and check the scaled normalized cut scores for its ground truth clusters

load YouTube_7truths
truth = GroundTruthVecs(:,7);

truth = fix_clus_vec(truth);

%%
volA = sum(nonzeros(A));
next = 1;
SnCuts = zeros(numel(unique(truth)),1);
for i = 1:max(truth)
    S = find(truth == i);
    
    snc = more_set_stats(A,S,volA);
    SnCuts(i) = snc;
    
end

min(SnCuts)
max(SnCuts)
%% For Email it seems like the best is 1e-4
    

n = size(A,1);
sns = zeros(n,1);
for i = 1:n
    sns(i) = SnCuts(truth(i));
end

min(sns)
median(sns)
max(sns)

S = find(truth == 3748);

AS = A(:,S);
vol = full(sum(nonzeros(AS)));
edges = full(sum(nonzeros(AS(S,:))));
cut = vol-edges;
cut