function [H,g,A,b] = matrixForm(n,u,d0)
% This function constructs H,g,A and b to express problem 2 as a quadratic
% program of the form 1/2x'Hx + g'x s.t. A'x = b

H = eye(n+1);
g = -u*ones(n+1,1);
At = [diag(-ones(n,1)) + diag(ones(n-1,1),-1) zeros(n,1)];
At(n,n+1) = -1;
At(1,n) = 1;
A = At';
b = zeros(n,1);
b(1) = -d0;

end