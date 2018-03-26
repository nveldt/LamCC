function [C,c,obj] = lam_cc(A,lam,T)
% LAM_CC Solve the lambda parameter correlation clustering objective
%
% This solves correlation clustering with a weight-matrix of
%   A - lambda*ones
% which gives a positive benefit of 1-lambda for each edge
% and a benefit of lambda for breaking each non-edge.
%
% [C,c,obj] = lam_cc(A,lambda)
% 
% Returns the cluster vector c,
% and the objective value obj. 
n = size(A,1);
if nargin < 3
    T = triangle_inequality(n);
end


assert(lam > 0)
assert(lam < 1) 
assert(all(nonzeros(A) == 1))

A = A - diag(diag(A));
tA = triu(A-lam,1);
n = size(A,1);

[rA,cA,vecA] = find(tA);

% Now that we have vectorized A, we can deal with this fixed ordering of
% the elements of A, with the mapping defined above.
wp = (1-lam)*(vecA>0);
wn = lam*(vecA<0);

RealObj = wp-wn; % The objective we minimize here is \sum_{i,j} A_{ij} d_{ij}
                 % If you actually want to count number of agreements, add
                 % sum(wp) to the output of Gurobi

% Run Gurobi using the whole matrix A

N = n*(n-1)/2; %number of variables x_ij

clear model
model.obj = RealObj; 

b = zeros(size(T,1),1);
sense = '<';

model.A = T;
model.rhs = b;
model.sense = sense;
model.vtype = repmat('B',1,N);
model.modelsense = 'min';

clear params;
params.outputflag = 0;
params.OptimalityTol = 1e-09;
%params.FeasibilityTol = 1e-09;
params.MIPGap = 1e-08;
result = gurobi(model, params);

ExactClustering = result.x;

% Add a constant to the objective to actually get weight of disagreements
obj = result.objval + lam*(N-nnz(A)/2);

G = zeros(n,n);

for i = 1:N
     connect = ExactClustering(i);
     G(rA(i),cA(i)) = 1-connect;
end
G = G + G';

% Performing pivot on G gives connected components
[C,c] = PivotCC(G);



