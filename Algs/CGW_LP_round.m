function [cBest,obj] = CGW_LP_round(A,Dist,k)
% Standarnd LP rounding scheme to accomplish 4-approximation by Charikar,
% Guruswami, and Wirth.
%
% This repeatedly calls the gamma_delta_rounding scheme k times and returns
% the best objective value.
%

gamma = 1/2;
delta = 1/4;

obj = nnz(A);
for i = 1:k
    [~,c] = gamma_delta_rounding(Dist,gamma,delta);
    
    o = CCminDisagreeObj(A,c);
    
    if o < obj
        obj = o;
        cBest = c;
    end
end