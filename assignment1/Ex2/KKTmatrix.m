function [K,d] = KKTmatrix(n,u,d0)
% This function constructs the KKT matrix for solving the quadratic 
% program in problem 2

[H,g,A,b] = matrixForm(n,u,d0);
K = [H -A; A' zeros(length(b))];
d = [-g; b];

end