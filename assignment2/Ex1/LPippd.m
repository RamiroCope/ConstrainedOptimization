function [x,Converged,mu,lambda,iter] = LPippd(g,A,b,x,lambda_init,mu_init)
% LPIPPD   Primal-Dual Interior-Point LP Solver
%
%          min  g'*x
%           x
%          s.t. A x  = b      (Lagrange multiplier: mu)
%                 x >= 0      (Lagrange multiplier: lambda)
%
% Syntax: [x,info,mu,lambda,iter] = LPippd(g,A,b,x)
%
%         info = true   : Converged
%              = false  : Not Converged

% Created: 04.12.2007
% Author : John Bagterp J�rgensen
%          IMM, Technical University of Denmark
% -------------------------------------------------------------------------
% Modified for Assignment 2 of Constrained Optimization

%
[m,n]=size(A);

maxit = 100;
tolL = 1.0e-9;
tolA = 1.0e-9;
tols = 1.0e-9;

eta = 0.99;

lambda = lambda_init; % ones(n,1);
mu = mu_init; % zeros(m,1);

% Compute residuals
[rL, rA, rC, s] = computeResiduals(g,A,b,x,mu,lambda,n);

% Converged
Converged = checkConvergence(rL,rA,s,tolL,tolA,tols);
%    
iter = 0;
while ~Converged && (iter<maxit)
    iter = iter+1;
    
    % ====================================================================
    % Form and Factorize Hessian Matrix
    % ====================================================================
    xdivlambda = x./lambda;
    H = A*diag(xdivlambda)*A';
    L = chol(H,'lower');
    
    % ====================================================================
    % Affine Step
    % ====================================================================
    % Solve
    tmp = (x.*rL + rC)./lambda;
    rhs = -rA + A*tmp;
    
    dmu = L'\(L\rhs);
    dx = xdivlambda.*(A'*dmu) - tmp;
    dlambda = -(rC+lambda.*dx)./x;
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(dlambda < 0.0);
    beta = min([1.0; -lambda(idx,1)./dlambda(idx,1)]);
    
    % ====================================================================
    % Center Parameter
    % ====================================================================
    xAff = x + alpha*dx;
    lambdaAff = lambda + beta*dlambda;
    sAff = sum(xAff.*lambdaAff)/n;
    
    sigma = (sAff/s)^3;
    tau = sigma*s;    

    % ====================================================================
    % Center-Corrector Step
    % ====================================================================
    rC = rC + dx.*dlambda - tau;
    
    tmp = (x.*rL + rC)./lambda;
    rhs = -rA + A*tmp;
    
    dmu = L'\(L\rhs);
    dx = xdivlambda.*(A'*dmu) - tmp;
    dlambda = -(rC+lambda.*dx)./x;
    
    % Step length
    idx = find(dx < 0.0);
    alpha = min([1.0; -x(idx,1)./dx(idx,1)]);
    
    idx = find(dlambda < 0.0);
    beta = min([1.0; -lambda(idx,1)./dlambda(idx,1)]);

    % ====================================================================
    % Take step 
    % ====================================================================
    x = x + (eta*alpha)*dx;
    mu = mu + (eta*beta)*dmu;
    lambda = lambda + (eta*beta)*dlambda;
    
    % ====================================================================
    % Residuals and Convergence
    % ====================================================================
    % Compute residuals
    rL = g - A'*mu - lambda;    % Lagrangian gradient
    rA = A*x - b;               % Equality Constraint
    rC = x.*lambda;             % Complementarity
    s = sum(rC)/n;              % Duality gap

    % Converged
    Converged = checkConvergence(rL,rA,s,tolL,tolA,tols);
end

% Return solution
if ~Converged
    x=[];
    mu=[];
    lambda=[];
end

end

function [ converged ] = checkConvergence(rL,rA,s,tolL,tolA,tols)
converged = (norm(rL,inf) <= tolL) && (norm(rA,inf) <= tolA) && (abs(s) <= tols);
end

function [rL,rA,rC,s] = computeResiduals(g,A,b,x,mu,lambda,n)
rL = g - A'*mu - lambda;    % Lagrangian gradient
rA = A*x - b;               % Equality Constraint
rC = x.*lambda;             % Complementarity
s = sum(rC)/n;              % Duality gap
end
    