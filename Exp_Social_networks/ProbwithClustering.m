function [PrCgivenK, PrKgivenC,PrK] = ProbwithClustering(A,info,C)
% C is a clustering, A is a social network unweighted adjacency matrix, and
% info is a set of attributes. Compute two sets of probabilities:
% PrCgivenK: probability that two people are in the same clustering given
% that they share an attribute
% PrKgivenC: probability that two people share an attribute given they are
% in the same clustering.

PrK = ProbabilitySharingAttributes(A,info); 
% PrK(i) give the probability that two random people share attribute i

PrC = ProbabilitySharingAttributes(A,C);
% gives the probability that two randomly chosen people are in the same
% cluster

assert(size(A,1) == size(info,1))

[n,k] = size(info);
PrCnK = zeros(k,1);             % probability of sharing both a specific 
                                % attribute and the same cluster
                                
TotalPairs = nchoosek(n,2);     % total number of pairs of people

for j = 1:k
    attribute = info(:,j);
    
    Share = 0;
    for a = 1:n-1
        for b = (a+1):n
            if attribute(a) == attribute(b) && C(a) == C(b) && attribute(a) > 0
                Share = Share + 1;
            end     
        end
    end
    PrCnK(j) = Share/TotalPairs;
end

PrCgivenK = PrCnK./PrK;
PrKgivenC = PrCnK/PrC;
   
end