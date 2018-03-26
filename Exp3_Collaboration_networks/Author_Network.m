%% Load the original data

clear
load dblp_2007_08.mat
AuPa= A';               % Author-Paper bipartite graph

%% Remove papers with over 50 authors, or less than 2 authors
B = AuPa;

%%
NumCollaborators = sum(B,1);
NoBigWorks = find(NumCollaborators < 51);
NoEmptyWorks = find(NumCollaborators > 1);

MainWorks = intersect(NoBigWorks,NoEmptyWorks);
B = B(:,MainWorks);
size(B)
%% Now remove any author that has no papers
Collaborators = find(sum(B,2) > 0);

B = B(Collaborators,:);
size(B)

%% And check that there are no empty rows or columns

NumCollaborators = sum(B,1);

assert(numel(find(NumCollaborators == 0)) == 0)

NumWorks = sum(B,2);
assert(numel(find(NumWorks == 0)) == 0)

%% Turn this into an actor-actor collaboration network
A = spones(B*B');
A = A-diag(diag(A));

text = 'DBLP authors dataset, kept papers with between 2 and 50 authors, and then removed any author with no papers.';
%save('AuthorPaper50','A','B','text')

%% Find the Work-Clique Solution
addpath('../Algs/')
c = WorkCliqueClustering(B);
cluster_deletion_objective(A,c)

%% Find the GrowClique solution
c = GrowClique(A,500);
cluster_deletion_objective(A,c)

%% As long as you have run the Julia code first, see the objective
% obtained by running the NP-hard 2-approximation

Clustering_txt2mat(char('AuthorPaper50_2approx.txt'))
load AuthorPaper50_2approx.txt.mat
Greedy2Score = cluster_deletion_objective(A,c)