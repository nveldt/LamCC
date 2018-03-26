function Pr = ProbabilitySharringAttributes(A,info)
% Pr = ProbabilitySharringAttributes(A,info): given a social network (from
% the FaceBook 100 dataset for example) and a list of attributes given in
% info, return a vector Pr which gives the probability that two random
% people in the network A share attribute k for each attribute in 'info'.
%
% A = n x n unweighted adjacency matrix
% info = n x k matrix, where k is the number of attributes

assert(size(A,1) == size(info,1))

[n,k] = size(info);
Pr = zeros(k,1);
TotalPairs = nchoosek(n,2);     % total number of pairs of people

for j = 1:k
    attribute = info(:,j);
    
    NumSame = 0;
    for a = 1:n-1
        for b = a+1:n
            if attribute(a) == attribute(b) && attribute(a)>0
                NumSame = NumSame + 1;
            end     
        end
    end
    Pr(j) = NumSame/TotalPairs;
           
end
   

end