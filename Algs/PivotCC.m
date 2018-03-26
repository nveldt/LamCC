function [C,cPiv] = PivotCC(A)
% Pivot Algorithm for correlation clustering, the algorithm developed by
% Ailon et. al.
% Also known as KwikCluster
%
% Reference:
% Nir Ailon, Moses Charikar, and Alantha Newman. 
% Aggregating inconsistent information: ranking and clustering. 
% Journal of the ACM (JACM), 55(5):23, 2008.

n = size(A,1);
c = ones(n,1);

Vinds = (1:n)';
C = [];

while numel(Vinds >0)
    i = randi(numel(Vinds));    % select random index from Vinds
    pivot = Vinds(i);           % index for pivot
    scores = A(pivot,Vinds);    % gives a list of similarities with pivot
    similar = find(scores >0);  % gives a list of those that are similar
    Similar = Vinds(similar);   % gets original indices
    
    NewCinds = [pivot; Similar];  % new cluster
    ci = zeros(n,1);
    ci(NewCinds) = 1;           % Forms a cluster indicator vector
    C = [C ci];
    Vinds = setdiff(Vinds,NewCinds);    % updates Vinds
end

cPiv = zeros(n,1);
for t = 1:size(C,2)
    cPiv = cPiv + t*C(:,t);
end

end
