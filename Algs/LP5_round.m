function [cBest,obj] = LP5_round(A,Dist,k)
% This is the rounding scheme for the LP5 algorithm outlined in the paper.
%
% This repeatedly calls the gamma_delta_rounding scheme k times and returns
% the best objective value.

gamma = 2/5;
delta = 1/5;
obj = nnz(A);
for i = 1:k
    [~,c] = gamma_delta_rounding(Dist,gamma,delta);
    
    o = CCminDisagreeObj(A,c);
    
    if o < obj
        obj = o;
        cBest = c;
    end
end