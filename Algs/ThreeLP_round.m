function [cBest,obj] = ThreeLP_round(A,Dist,k)
% Given the distance matrix Dist from solving the LP relaxtion for
% LambdaCC, this gives the randomized rounding scheme for ThreeLP, which
% selects a node uniformly at random and clusters it with all unclusterd
% neighbors that have distance less than three to the pivot.

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
    Dist(pivot,pivot) = 0;
    
    % Take everything that is within 1/3 and isn't already in another
    % cluster
    Close = find(Dist(pivot,:) < 1/3);
    Cluster = intersect(Vinds,Close);
    
    NewCinds = [pivot; Cluster];  % new cluster
    ci = zeros(n,1);
    ci(NewCinds) = 1;           % Forms a cluster indicator vector
    C = [C ci];
    Vinds = setdiff(Vinds,NewCinds);    % updates Vinds
end
m = size(C,2);
v = 1:m;
c = C*v';

end