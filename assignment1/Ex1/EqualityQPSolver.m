function [x,lambda] = EqualityQPSolver(H,g,A,b)

if isempty(A) && isempty(b)
    x = H\-g;
    lambda = [];
else
m = size(H,1);
z = [H -A; A' zeros(length(b))]\[-g; b];
x = z(1:m);
lambda = z(m+1:end);
end

end