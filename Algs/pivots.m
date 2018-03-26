function [cBest,obj] = pivots(A,k)
% Run the pivot algorithm k times on the matrix A and take the best result
obj = nnz(A);
for i = 1:k
    [~,c] = PivotCC(A);
    
    o = CCminDisagreeObj(A,c);
    
    if o < obj
        obj = o;
        cBest = c;
    end
end