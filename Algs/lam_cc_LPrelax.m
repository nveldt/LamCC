function [Dist,Bound] = lam_cc_LPrelax(A,lam,T)
% LAM_CC_LPRELAX: return the distances obtained from solving the LP
% relaxation of the lambda-CC objective
%
% This returns:
% * Dist: "distances" between nodes, useful for forming
%   clusterings via LP rounding solutions
% * Bound: the objective value found by the LP relaxation, useful for
%   giving a lower bound on the optimal objective

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

% Add a constant to the objective to actually get weight of disagreements
model.objcon = lam*(N-nnz(A)/2);

% In order to relax, we introduce constraints xij \in (0,1
b = zeros(size(T,1),1);
T = [T;speye(N)];
b = [b;ones(N,1)];
sense = '<';

model.A = T;

model.rhs = b;
model.sense = sense;
model.vtype = repmat('C',1,N);
model.modelsense = 'min';

clear params;
params.outputflag = 0;
params.OptimalityTol = 1e-9;    % minimum value for tolerance
params.FeasibilityTol = 1e-9;   % also the minimum value for tolerance
result = gurobi(model, params);

xij = result.x;

Bound = result.objval;

Dist = sparse(rA,cA,xij,n,n);   % list of distances
Dist = Dist + Dist';




