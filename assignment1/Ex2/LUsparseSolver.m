function [x,lambda] = LUsparseSolver(n,u,d0)
% This problem uses LU factorization for solving the quadratic program in
% problem 2
[K,d] = KKTmatrix(n,u,d0);
K = sparse(K);
[L,U,p] = lu(K,'vector');
z = U\(L\d(p));

x = z(1:n+1);
lambda = z(n+2:end);
end
