function BestClus = GrowCliqueStep(A,times)
% GROWCLIQUESTEP: given a network A, this algorithm searches for a maximal
% random clique in the network for the specificed number of "times", and
% then returns the largest of them
%
% Coded by Nate Veldt on 5-4-2017

BestClus = 0;
BestSize = 0;
n = size(A,1);
for i = 1:times
    clus = randi(n);       % current cluster made of starter/seed/pivot node   

    % We want to efficiently keep track of del_v scores

    sizeC = 1;    % keep |C| updated
    cutC = A(:,clus);       % keep cut(C,v) updated    
    assert(size(cutC,2) == 1);   % this should be a vector


    Neigh = cutC;                 % We can merge anything that neighbors 
                                  % EVERYTHING in cluster C, and isn't already
                                  % in C.

    while nnz(Neigh) > 0

        canAdd = find(Neigh);    % take all possible nodes we could add
        j = randi(numel(canAdd));     % select one at random
        toAdd = canAdd(j);

        clus = [clus;toAdd];        % add it to the current cluster
        sizeC = sizeC + 1;          % our cluster grew by one

        newNeighbors = A(:,toAdd);  % find all neighbors of the new node

        Neigh = Neigh.*newNeighbors; % now only consider neighbors that are also neighbors
                                            % of the new node "toAdd"
        Neigh(toAdd) = 0;

    end

    if sizeC > BestSize
        BestSize = sizeC;
        BestClus = clus;
    end

end

end