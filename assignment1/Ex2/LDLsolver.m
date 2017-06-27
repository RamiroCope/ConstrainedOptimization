function [x,lambda] = LDLsolver(n,u,d0)
% This function uses LDL factorization for solving the quadratic program in
% problem 2

[K,d] = KKTmatrix(n,u,d0);

z = zeros(size(K,1),1);
[L,D,p] = ldl(K,'lower','vector');
z(p) = L'\(D\(L\d(p)));
x = z(1:n+1);
lambda = z(n+2:end);
end