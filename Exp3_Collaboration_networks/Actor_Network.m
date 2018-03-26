%% Load the original data

clear
load NotreDame_actors.mat
AcMo= Problem.A;

%% Remove movies with more than 50 actors
B = AcMo;

NumCollaborators = sum(B,1);
NoBigWorks = find(NumCollaborators < 51);
NoEmptyWorks = find(NumCollaborators > 0);

MainWorks = intersect(NoBigWorks,NoEmptyWorks);
B = B(:,MainWorks);

%% Now remove any actor not in a movie
Collaborators = find(sum(B,2) > 0);

B = B(Collaborators,:);

%% And check that there are no empty rows or columns

NumCollaborators = sum(B,1);

assert(numel(find(NumCollaborators == 0)) == 0)

NumMovies = sum(B,2);
assert(numel(find(NumMovies == 0)) == 0)

%% Turn this into an actor-actor collaboration network
A = spones(B*B');
A = A-diag(diag(A));

text = 'Actors dataset, removed movies with more than 50 actors and then removed any actor not in a movie.';
save('ActorsMovies50','A','B','text')

%% Find the Movie-Clique Solution
c = WorkCliqueClustering(B);
cluster_deletion_objective(A,c)

%% Find the GrowClique solution
c = GrowClique(A,500);
cluster_deletion_objective(A,c)

%% As long as you have run the Julia code first, see the objective
% obtained by running the NP-hard 2-approximation

Clustering_txt2mat(char('ActorMovie50_2approx.txt'))
load ActorMovie50_2approx.txt.mat
Greedy2Score = cluster_deletion_objective(A,c)