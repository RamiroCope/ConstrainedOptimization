%
% Test LP solver
% Original from John, 
% modified for Assignment 2 
%
m = 100;
n = 1000;
seed = 12345;
rng(seed, 'twister'); % updated so new Matlab versions don't complain
% Generate random LP problems
fprintf('\nGenerating random LP problem (N&W 14.15 p419-420)\n');
A = randn(m,n);
x = zeros(n,1);
x(1:m) = abs(rand(m,1));
lambda = zeros(n,1);
lambda(m+1:n) = abs(rand(n-m,1));
mu = rand(m,1);
g = A'*mu + lambda;
b = A*x;
% Choose the starting point (x0, lambda0, s0) 
% with the components set to large positive values.
scalef = 100; % "large values" -> scale factor
x_init = scalef*ones(n,1);
lambda_init = scalef*ones(size(lambda));
mu_init = scalef*ones(size(mu));
% Run
fprintf('\nRunning primal dual interior point LP solver...\n');
[xlp,converged,mulp,lambdalp,iter,tol] = LPSolver(g,A,b,x_init,lambda_init,mu_init);
fprintf('Iterations: %d\n', iter);
% Check the computed vs "real" is the algorithm converged
if converged
    fprintf('Converged!\nErrors (solver tolerance %e)\n', tol);
    fprintf('x:      %e\n', max(abs(xlp-x)));
    fprintf('mu:     %e\n', max(abs(mulp-mu)));
    fprintf('lambda: %e\n', max(abs(lambdalp-lambda)));
end
