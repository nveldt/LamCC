%% Read in EdgeList

file = fopen('email-Eu-core.txt','r');
formatSpec = '%d %d';
sizeA = [2 Inf];

Edges = fscanf(file,formatSpec,sizeA)' +1 ;
fclose(file)
%% Form the graph
n = 1005;
A = sparse(Edges(:,1),Edges(:,2),ones(numel(Edges(:,1)),1),n,n);

%%
A = A+A';
d = sum(A,2);

%% Read in the communities
file = fopen('email-Eu-core-department-labels.txt');
Labels = fscanf(file,formatSpec,sizeA)' + 1;

truth = Labels(:,2);

save('Email','A','truth')