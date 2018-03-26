function Obj = CCminDisagreeObj(A,c)
% Obj = CCminDisagreeObj(A,c)
% Returns the weight of disagreements correlation clustering 
% objective based for the clustering c, 
% given the coefficients from the matrix A, a signed network.
%
% c can either be a cluster vector or an indicator matrix
%
% This code is slow since it iteratively visits all pairs of nodes. For
% sparse matrices we can do better than computing the objective score in
% O(n^2) time.

[n,t] = size(c);
if t > 1 && max(max(c)) ==1
    % then we were given an indicator matrix
    cOld = c;
    c = zeros(n,1);
    for i = 1:t
        c = c + i*cOld(:,i);
    end
end

Obj = 0;

A = full(A); 
for j=1:n-1
    for i=(j+1):n
        if A(i,j) > 0 && c(i) ~= c(j)
            Obj = Obj + A(i,j);
      
        end
        
        if A(i,j) < 0 && c(i) == c(j)
            Obj = Obj - A(i,j);
        end            
    end
end