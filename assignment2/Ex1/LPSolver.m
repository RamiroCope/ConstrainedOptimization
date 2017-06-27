function [x,Converged,mu,lambda,iter,tol] = LPSolver(g,A,b,x,lambda_init,mu_init)
% LPIPPD   Primal-Dual Interior-Point LP Solver
%
%          min  g'*x
%           x
%          s.t. A x  = b      (Lagrange multiplier: mu)
%                 x >= 0      (Lagrange multiplier: lambda)
%
% Syntax: [x,Converged,mu,lambda,iter,tol] = LPsolver(g,A,b,x,lambda_init,mu_init)
%

% Created: 04.12.2007
% Author : John Bagterp J�rgensen
%          IMM, Technical University of Denmark
% -------------------------------------------------------------------------
% Original modified for Assignment 2 of Constrained Optimization (9.5.2017)
% Uses simplified parts from Assignment 1 PDPCIPA QP

% Initialize
n = length(A);
lambda = lambda_init;
mu = mu_init;
tol = 1e-9;
tolL = tol;
tolA = tol;
tols = tol;
eta = 0.99;
maxiter = 100;
iter = 0;
% Compute residuals
[rL, rA, rC, s] = computeResiduals(g,A,b,x,mu,lambda,n);
% Convergence check
Converged = checkConvergence(rL,rA,s,tolL,tolA,tols);

while ~Converged && (iter < maxiter)
    iter = iter+1;
    % Hessian (save factorization for speed)
    xdivlambda = x./lambda;
    H = A*diag(xdivlambda)*A';
    L = chol(H,'lower');
    % Affine step -> solve
    [dx,dmu,dlambda] = stepSolve(x,lambda,xdivlambda,A,L,rA,rL,rC);
    % Affine step -> step length
    [alpha, beta] = stepLength(x,dx,lambda,dlambda);
    % Center parameter
    xAff = x + alpha*dx;
    lambdaAff = lambda + beta*dlambda;
    sAff = sum(xAff.*lambdaAff)/n;      % duality gap for affine step
    sigma = (sAff/s)^3;                 % centering parameter 
    % Center+Corrector step to follow tractory to the primal dual solution set
    rC = rC + dx.*dlambda - sigma*s;
    [dx,dmu,dlambda] = stepSolve(x,lambda,xdivlambda,A,L,rA,rL,rC);
    % Step length
    [alpha, beta] = stepLength(x,dx,lambda,dlambda);
    % Take step - update from (45) 
    x = x + eta*alpha*dx;
    mu = mu + eta*beta*dmu;
    lambda = lambda + eta*beta*dlambda;
    % Compute residuals
    [rL, rA, rC, s] = computeResiduals(g,A,b,x,mu,lambda,n);
    % Check if converged
    Converged = checkConvergence(rL,rA,s,tolL,tolA,tols);
end
% Signal the user that the algorithm didn't converge
if ~Converged 
    throw(MException('dtu:copt', 'Did not converge!'));
end

end

function [converged] = checkConvergence(rL,rA,s,tolL,tolA,tols)
% Stopping criteria from (48)
converged = (norm(rL) <= tolL) && (norm(rA) <= tolA) && (abs(s) <= tols);
end

function [rL,rA,rC,s] = computeResiduals(g,A,b,x,mu,lambda,n)
% Following eqns (46)
rL = g - A'*mu - lambda;
rA = A*x - b;
rC = x.*lambda;
s = sum(rC)/n;
end

function [alpha,beta] = stepLength(x,dx,lambda,dlambda)
% Find alpha_aff and beta_aff using (35) and (36)
idx_a = find(dx < 0);
alpha = min(1, min(-x(idx_a)./dx(idx_a)));
idx_b = find(dlambda < 0);
beta = min(1, min(-lambda(idx_b)./dlambda(idx_b)));
end

function [dx,dmu,dlambda] = stepSolve(x,lambda,xdivlambda,A,L,rA,rL,rC)
% Solve for dx, dmu, dlambda
rhs = -rA + A*((x.*rL + rC)./lambda);
dmu = L'\(L\rhs);
dx = xdivlambda.*(A'*dmu) - (x.*rL + rC)./lambda;
dlambda = -(rC+lambda.*dx)./x;
end
    