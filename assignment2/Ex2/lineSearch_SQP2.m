function [x,lambda,info] = local_SQP(objfun,const,x,lambda,tol,maxit)
% The following function solves the nonlinear optimization problem in
% exercise 18.3 in N&W.

% initial evaluation of the function and the constraints
[f,g,H] = objfun(x);
[c,A,d2c1,d2c2,d2c3] = const(x);

dxL = g-A'*lambda; % gradient of the Lagrangian  
F = [dxL; c]; %% KKT conditions
k = 1;
n = length(x);
x(:,k) = x;
% While the KKT conditions are not satisfied and the number of iterations
% is less than the threshold
while ((k < maxit) && (norm(F,'inf') > tol))
    % compute hessian of the Lagrangian
    d2L = H - lambda(1)*d2c1 - lambda(2)*d2c2 - lambda(3)*d2c3;
    % solve the system to get the Newton's descent direction
    z = [d2L -A'; A zeros(3)]\[-dxL; -c];
    p = z(1:n);
    plambda = z(n+1:end);
    % Backtraking line search method
    sigma = 1; % Hessian positive definite
    rho = 0.5; % rho in (0,1)
    mu = (g'*p + sigma/2*p'*d2L*p)/((1-rho)*norm(c,1));
    alpha = 1; % full step
    tau = 0.8; % contraction factor
    
    % check sufficient decrease condition
    fn = objfun(x(:,k)+alpha*p);
    cn = const(x(:,k)+alpha*p);
    eta = 0.4;
    while fn+mu*norm(cn,1) > f+mu*norm(c,1)+eta*alpha*(g'*p-mu*norm(c,1))
        alpha = tau*alpha; % update alpha 
        fn = objfun(x(:,k)+alpha*p);
        cn = const(x(:,k)+alpha*p);
    end
    f = fn;
    alphaval(k) = alpha; % Store alfa value in every iteration
    % obtain next point in the iteration sequence
    x(:,k+1) = x(:,k) + alpha*p;
    % update lambda value
    lambda = lambda + alpha*plambda;
    
    k = k+1;
    % evaluate function and constraints on the new point
    [~,g,H] = objfun(x(:,k));
    [c,A,d2c1,d2c2,d2c3] = const(x(:,k));
    dxL = g-A'*lambda; % gradient of the Lagrangian
    % check KKT conditions
    F = [dxL;
        c];
end
info.iter = k; % return number of iterations
info.seq = x; % return iteration sequence
info.alpha = alphaval; % return step lenght sequence
x = x(:,k); % return last point as solution
end

