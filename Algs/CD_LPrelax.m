function [Dist,Bound] = CD_LPrelax(A,T)
% CD_LPRELAX: return the distances obtained from solving the LP
% relaxation of the cluster deletion objective
%
% This returns:
% * Dist: "distances" between nodes, useful for forming
%   clusterings via LP rounding solutions
% * Bound: the objective value found by the LP relaxation, useful for
%   giving a lower bound on the optimal objective

n = size(A,1);
if nargin < 2
    T = triangle_inequality(n);
end

assert(all(nonzeros(A) == 1))

% Vectorize the matrix, but only for entris i < j
N = nchoosek(n,2);
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
% Force distance to be 1 for non-edges
Con = [T;-speye(N)];
rhs = [b;vecA-1];

model.A = Con;
model.rhs = rhs;
model.sense = sense;
model.vtype = repmat('B',1,N);
model.modelsense = 'min';

clear params;

% Have strict tolerance constraints
params.outputflag = 1;
params.OptimalityTol = 1e-9;    % minimum value for tolerance
params.FeasibilityTol = 1e-9;   % also the minimum value for tolerance
result = gurobi(model, params);

xij = result.x;

Bound = result.objval;

Dist = sparse(rA,cA,xij,n,n);   % list of distances
Dist = Dist + Dist';