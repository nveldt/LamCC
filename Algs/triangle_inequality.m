function C = triangle_inequality(n)
%  TRIANGLE_INEQUALITY(n): Returns the constraints matrix that encodes the constraints:
% d_ij <= d_ik + d_jk
% for every disjoint pair of nodes i and j in an undirected graph 
% with n nodes
%
% Given an adjacency matrix A, we choose the linear ordering on the d_ij 
% defined by the 'find' function,
% i.e. we count of the d_ij by first counting through columns and then
% counting through rows.

G = ones(n) - eye(n);
G = triu(G);

%p now is the number of dij's, not the number of edges
[r,c,~] = find(G);
p = nnz(G);
pp = 1:p;

% The opposite direction mapping
MatToLin = sparse(r,c,pp',n,n);
MatToLin = MatToLin + MatToLin'; % make it symmetric to give us 
                                 % both indexings 

% 3 non-zeros per row, this many rows, 
% (n choose 2) pairs ij, n-2 pairs ijk with k distinct
nzs = 3*(n*(n-1)/2*(n-2)); 

ei = zeros(nzs,1);
ej = zeros(nzs,1);
ev = zeros(nzs,1);

row = 0;
count = 0;
for i=(1:n)
    for j=(i+1):n
        for k=1:n
            if k==i || k == j
                continue
            end
            ij = MatToLin(i,j);
            ik = MatToLin(i,k);
            jk = MatToLin(j,k);
            
            row = row + 1;
            
            count = count + 1;
            ei(count) = row;
            ej(count) = ij;
            ev(count) = 1;
            
            count = count + 1;
            ei(count) = row;
            ej(count) = jk;
            ev(count) = -1;
            
            count = count + 1;
            ei(count) = row;
            ej(count) = ik;
            ev(count) = -1;
        end
    end
end

C = sparse(ei,ej,ev, n*(n-1)/2*(n-2), n*(n-1)/2);

end