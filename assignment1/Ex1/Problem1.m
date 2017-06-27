clear all
%% Solve Quadratic Problem
H = [6 2 1; 2 5 2; 1 2 4];
g = [-8; -3; -3];

A = [1 0 1; 0 1 1]';
b = [3; 0];

[x,lambda] = EqualityQPSolver(H,g,A,b);

%% Random convex quadratic programs
for j = 1:100
p = 1;
while p ~= 0
H = rand(3);
[~,p] = chol(H);
end
xrand = rand(3,1);
A = rand(3,2);

b = A'*xrand;
g = -H*xrand;


[xt,lambdat] = EqualityQPSolver(H,g,A,b);
if norm(xt-xrand) > 10e-4
    disp('The solution is not correct')
end
end
%% Sensivity
H = [6 2 1; 2 5 2; 1 2 4];
g = [-8; -3; -3];

A = [1 0 1; 0 1 1]';
b = [3; 0];

cx = A;
Wxx = H;
Wxp = [eye(3); zeros(2,3)];
cp = [zeros(3,2); -eye(2)];
[Spx,Splambda] = EqualityQPSensitivity(Wxp,cp,Wxx,cx);
epsilon = 10e-4;

% Sensitivity testing
for i = 1:3
    gt = g;
    bt = b;

    gt(i) = g(i) + epsilon;

    [xt,lambdat] = EqualityQPSolver(H,gt,A,bt);
    Spxt(i,:) = (xt - x)./epsilon;
    Splambdat(i,:) = (lambdat - lambda)./epsilon;
end


for i = 1:2
    gt = g;
    bt = b;

    bt(i) = b(i) + epsilon;

    [xt,lambdat] = EqualityQPSolver(H,gt,A,bt);
    Spxt(i+3,:) = (xt - x)./epsilon;
    Splambdat(i+3,:) = (lambdat - lambda)./epsilon;
end

if norm(Spxt - Spx) > 10e-4 || norm(Splambdat - Splambda) > 10e-4
    disp('The solution is not correct')
end

%% Dual program
t = cputime;
[x,lamdba] = EqualityQPdual(H,g,A,b);
