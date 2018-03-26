function [B,wB] = compress(A,c,w)
% Given a network A and a clustering c, compress/boil down/agglomerate
% nodes in the same cluster into a new weighted network B

numclus = max(c);

% w is the weights of former nodes
wB = zeros(numclus,1);
B = sparse(numclus,numclus);
for i = 1:numclus-1
    % the new ith supernode has weight equal 
    %  to the sum of weights of the previous nodes
    wB(i) = sum(w(c == i));
    for j = i+1:numclus
        B(i,j) = sum(sum(A(c == i,c == j)));
    end
end
wB(end) = sum(w(c == numclus));
B = B+B';


end