function [x,lambda,info] = lineSearch_SQP(objfun,const,x,lambda,tol,maxit)

% The following function solves the nonlinear programming problem in
% exercise 18.3 in N&W using the damped BFGS approximation to the Hessian
% of the Lagrangian.

% initial evaluation of the function and the constraints
[f,g] = objfun(x);
[c,A] = const(x);

dxL = g-A'*lambda; % gradient of the Lagrangian  
F = [dxL; c];

k = 1;
n = length(x);
x(:,k) = x;

B = eye(n); % Hessian approximation initial value

while ((k < maxit) && (norm(F,'inf') > tol))
    
    % solve the system to get the search direction
    z = [B -A'; A zeros(3)]\[-dxL; -c];
    p = z(1:n);
    plambda = z(n+1:end);
    
    % Backtraking line search method
    sigma = 1; % Hessian positive definite
    rho = 0.1; % rho in (0,1)
    mu = (g'*p + sigma/2*p'*B*p)/((1-rho)*norm(c,1));
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
    [f,g] = objfun(x(:,k));
    [c,A] = const(x(:,k));
    dxLn = g-A'*lambda; % new gradient of the Lagrangian
    F = [dxLn; c];
    
    % damped BFGS updating
    s = x(:,k) - x(:,k-1);
    y = dxLn - dxL;
    dxL = dxLn;
    
    if s'*y >= 0.2*s'*B*s 
        theta = 1;
    else
        theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
    end
    r = theta*y+(1-theta)*B*s; % ensure possitive definiteness
    
    B = B - B*(s*s')*B/(s'*B*s) + r*r'/(s'*r); % B updating

    % ensure possitive definiteness
end
% store number of iterations and iteration sequence
info.iter = k; 
info.seq = x;
info.alpha = alphaval; % return step lenght sequence
x = x(:,k); % return last point as solution
end