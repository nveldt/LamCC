function L = laplacian(A,n)
% LAPLACIAN Compute the (combinatorial) Laplacian matrix for the matrix A.
%
% L = laplacian(A) returns the combinatorial Laplacian matrix for a 
% matrix A.  This function works for both a dense A, a sparse A, and an
% operator A (passed as a function handle).  The return type is the same as
% the input type.  
%
% Formally, L = D - A where D = diag(A*e) is the diagonal
% matrix of degrees.  This formulation handles weighted graphs in the 
% natural way, i.e. the "degree" is the row-sum.
%
% If A is a function handle, then you MUST pass the additional argument "n"
% to denote the size.
%
% Example:
%  % Put in Laplacian eigenmaps example here...

% History
% :2011-03-10: Initial coding

if isa(A,'function_handle')
    e = ones(n,1);
    fhand = true;
    d = A(e);
else
    fhand = false;
    d = sum(A,2);
end

if fhand
    d = full(d); 
    L = @(x) d.*x - A(x);
else
    L = diag(d) - A;
end
