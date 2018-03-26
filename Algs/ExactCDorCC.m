function [C,c,obj] = ExactCDorCC(A,solveCC,T)
% [C,c,obj] = ExactCDorCC(A,T):
% Can be used to solve the exact Correlation Clustering objective
% or the exact cluster deletion objective on an unweighted, undirected
% network (instead of positive and negative edges we just have edges and
% non-edges.
%
% The correlation clustering objective on an unweighted, undirected,
% unsigned network is:
%
%       min \sum_{i < j} (1-A_{ij}) + 2*A_{ij}*x_{ij} - x_{ij}
%
% (minimizing disagreements, not maximizing agreements)
% which in vector notation is
%
%       min e'*(e-vecA) + + 2*vecA'*x - e'*x
%
% Where vecA is the linearzation of the upper trianglular part
% of the adjacency matrix (of length N = n(n+1)/2), e = ones(N,1)
% and x is the length N binary vector indicating "distances" between two
% nodes. x_{ij} = 0 if i and j are clustered together.
%
% A is the adjacency matrix for 
% Set OnlyCliques = 0 if you want the cluster deletion objective

n = size(A,1);
N = n*(n-1)/2; %number of variables x_ij
if nargin < 3
    T = triangle_inequality(n);
end

assert(all(nonzeros(A) == 1))

% Vectorize the matrix, but only for entris i < j
rA = zeros(N,1);
cA = zeros(N,1);
vecA = zeros(N,1);
next = 0;
for j = 2:n
    for i = 1:j-1
        next = next+1;
        vecA(next) = A(i,j);
        rA(next) = i;
        cA(next) = j;
    end
end

% This is the part of the objective without the first constant out front
e = ones(N,1);
RealObj = 2*vecA - e;

clear model
model.obj = full(RealObj);

% The constant offset = e'(e-vecA) = N - sum(vecA);
model.objcon = N - sum(vecA);
b = zeros(size(T,1),1);
sense = '<';

% Constraints matrix:
% if you want clique clusters, add an extra constraint
if solveCC == 0
    Con = [T;-speye(N)];
    rhs = [b;vecA-1];
else
    Con = T;
    rhs = b;
end

model.A = Con;
model.rhs = rhs;
model.sense = sense;
model.vtype = repmat('B',1,N);
model.modelsense = 'min';

clear params;

% Have strict tolerance constraints
params.outputflag = 1;
%params.OptimalityTol = 1e-09;

if solveCC
    %if we're doing correlation clustering, we gain nothing by presolving
    params.presolve = 0; 
end

params.MIPGap = 0;
result = gurobi(model, params);

ExactClustering = result.x;

obj = result.objval;

G = zeros(n,n);

for i = 1:N
     connect = ExactClustering(i);
     G(rA(i),cA(i)) = 1-connect;
end
G = G + G';

% Performing pivot on G gives connected components
[C,c] = PivotCC(G);