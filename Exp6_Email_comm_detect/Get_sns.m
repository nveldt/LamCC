function sns = Get_sns(A,c)
% Get a vector of the scaled normalized cut value of the cluster that node
% i is in

volA = sum(nonzeros(A));
SnCuts = zeros(max(c),1);

for i = 1:max(c)
    S = find(c == i);
    
    snc = more_set_stats(A,S,volA);
    SnCuts(i) = snc;
end

n = size(A,1);
sns = zeros(n,1);
for i = 1:n
    sns(i) = SnCuts(c(i));
end

end