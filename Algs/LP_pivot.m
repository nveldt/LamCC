function [cBest,obj] = LP_pivot(A,Dist,k)
% LPPIVOT: an implementation of the 2.5-approximation for probability
% constraints correlation clustering by Ailon et al.
% Run it k times. 
% This takes as input the distances from the LP relaxation, already assumed
% to be solved.
%
% Since we are using it for lambdaCC, we do take lambda as input and output
% the clustering out of the k roundings that gives the best objective score

obj = nnz(A);
for i = 1:k
    c = OneRun(Dist);
    
    o = CCminDisagreeObj(A,c);
    
    if o < obj
        obj = o;
        cBest = c;
    end
end


end

function c = OneRun(Dist)
% Runs the LP pivot algorithm just once

n = size(Dist,1);
Vinds = (1:n)';
C = [];

while numel(Vinds >0)
    i = randi(numel(Vinds));    % select random index from Vinds
    pivot = Vinds(i);           % index for pivot
    
    % Now go through all unclustered items, and add all the objects with
    % probability 1-dij, which you get from the distance matrix.
    Cluster = [];
    for j = 1:numel(Vinds)
        candidate = Vinds(j);
        dist = Dist(pivot,candidate);
        prob = 1-dist;
        
        if rand(1) < prob
            Cluster = [Cluster; candidate];           
        end
    end
    
    NewCinds = [pivot; Cluster];  % new cluster
    ci = zeros(n,1);
    ci(NewCinds) = 1;           % Forms a cluster indicator vector
    C = [C ci];
    Vinds = setdiff(Vinds,NewCinds);    % updates Vinds
end
m = size(C,2);
v = 1:m;
c = C*v';
C;
end