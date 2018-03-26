function c = WorkCliqueClustering(B)
% WORKCLIQUECLUSTERING: form a clustering of collaborators by grouping them
% in the work with the most collaborators, then remove this work and its
% collaborators from matrix B do this recursively until all collaborators
% are clustered in some clique.
%
% B is a bipartite matrix storing colalborator-work connections
%
% Instead of collaborators we can talk about actors, and movies instead of
% works.
%
% This is inefficient code, intended simply to find the solution without
% regard for runtime.

[m,n] = size(B);

% We keep an updated list of the indices we are considering 
ActorsLeft = (1:m)';
MoviesLeft = (1:n)';
clusterNumber = 1;
c = zeros(m,1);
C = [];
while numel(ActorsLeft > 0)
    
    % We have a bipartite matrix, and the indices that corresponds to each

    % Find the "fullest" movie remaining
    [~,wheremax] = max(sum(B,1));

    % themax is the largest number of actors in any remaining movie
    % wheremax is the index in the "remaining B" of "themax". We must change to
    % get the index from the original dataset:

    MovInd = MoviesLeft(wheremax); % index of movie in original dataset

    newActInds = find(B(:,wheremax)); % the indices of actors in this movie that haven't already been clustered
                                    % but these indices correspond to
                                    % indices in the ActorsLeft list

    ActInds = ActorsLeft(newActInds);   % these are the original indices

    % Update B
    [mb,~] = size(B);
    B = B(setdiff(1:mb,newActInds),[1:wheremax-1, wheremax+1:end]);


    c(ActInds) = clusterNumber;
    clusterNumber = clusterNumber + 1;
    if mod(clusterNumber,500) == 0
        fprintf('ActorsLeft = %d, numclus = %d \n',numel(ActorsLeft),clusterNumber);
    end
    % Update which movies and actors are left
    ActorsLeft = setdiff(ActorsLeft,ActInds);
    MoviesLeft = setdiff(MoviesLeft,MovInd);

end