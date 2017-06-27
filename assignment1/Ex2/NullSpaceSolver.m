function [x,lambda] = NullSpaceSolver(n,u,d0)
% This problem uses the null-space method for solving the quadratic program in
% problem 2
[H,g,A,b] = matrixForm(n,u,d0);

[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1);
Q2 = Q(:,m1+1:end);
R = Rbar(1:m1,1:m1);
xy = R'\b;
xz = (Q2'*H*Q2)\(-Q2'*(H*Q1*xy+g));
x = Q1*xy + Q2*xz;
lambda = R\Q1'*(H*x+g);
end