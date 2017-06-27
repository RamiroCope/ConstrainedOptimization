function [x,lambda] = RangeSpaceSolver(n,u,d0)
% This problem uses the range-space method for solving the quadratic program in
% problem 2
[H,g,A,b] = matrixForm(n,u,d0);
L = chol(H,'lower');
v = L'\(L\g);
La = chol(A'*(L'\(L\A)),'lower');
lambda = La'\(La\(b+A'*v));
x = L'\(L\(A*lambda-g));
end
