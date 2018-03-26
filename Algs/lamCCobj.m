function Obj = lamCCobj(A,lam,c)
% Given a sparse graph A and a parameter lambda, and a clustering vector c,
% this quickly computes the lambda-CC objective corresponding to c. 
%
% Coded by Nate Veldt on May 5, 2017. Updated 9-05-17
%
n = size(A,1);
if issparse(A) == 0
    A = sparse(A);
end

Obj = lam*(n*(n-1)/2 - nnz(A)/2);

numClus = max(c);

for i = 1:numClus
    ClusterI = find(c == i);
    ObjI = clusterScore(A,ClusterI,lam,n);
    
    Obj = Obj + .5*ObjI;
end


end

function y = clusterScore(A,S,lam,n)
% This returns cut(S) - lam*|S|*|V\S|
AS = A(:,S);
vol = full(sum(nonzeros(AS)));
edges = full(sum(nonzeros(AS(S,:))));
cut = vol-edges;
sizeS = numel(S);
y = cut - lam*sizeS*(n-sizeS);

end