function [x,lambda,info] = BFGS_SQP(objfun,const,x,lambda,tol,maxit)

% The following function solves the nonlinear programming problem in
% exercise 18.3 in N&W using the damped BFGS approximation to the Hessian
% of the Lagrangian.

% initial evaluation of the function and the constraints
[~,g] = objfun(x);
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
    % obtain next point in the iteration sequence
    p = z(1:n);
    x(:,k+1) = x(:,k) + p;
    % update lambda value
    lambda = lambda + z(n+1:end);
    k = k+1;
    % evaluate function and constraints on the new point
    [~,g] = objfun(x(:,k));
    [c,A] = const(x(:,k));
    dxLn = g-A'*lambda; % new gradient of the Lagrangian
    % check KKT conditions
    F = [dxLn; c];
    % damped BFGS updating
    s = x(:,k) - x(:,k-1);
    y = dxLn - dxL;
    dxL = dxLn;
     % ensure possitive definiteness
    if s'*y >= 0.2*s'*B*s 
        theta = 1;
    else
        theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
    end
    r = theta*y+(1-theta)*B*s;
    % update Bs
    B = B - (B*(s*s')*B)/(s'*B*s) + r*r'/(s'*r); 
end
% store number of iterations and iteration sequence
info.iter = k; 
info.seq = x;
x = x(:,k); % return last point as solution
end