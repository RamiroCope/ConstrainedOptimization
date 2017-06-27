function [x,lambda,info] = ineqSQP(x,lambda,tol,maxit)

% Initializing info storage variables
subQPiter = zeros(1,maxit); 
subQPfval = zeros(1,maxit);
fcalls=1;
k=1;
n=length(x);
x(:,k) = x; % storing solution iterates

% Initializing problem
B = eye(2);
[~,g,~] = objfun(x); [cineq,dcineq,~,~] = constraint(x);
F = [g - dcineq'*lambda; cineq];
 
QPoptions = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'none');

while ((k < maxit) && (norm(F(1:2),'inf') > tol))   %&& (norm(cineq,'inf') < tol)) 
  
    [f,g,~] = objfun(x(:,k)); fcalls=fcalls+1;
    [cineq,dcineq,~,~] = constraint(x(:,k));
    
    % Defining variables
    gL = g - dcineq'*lambda;   %gradient of Lagrangian
    
    % solve system to get search direction
    [incremental, fval,~,output, lambda_qp] = quadprog(B,g,-dcineq, cineq, [], [], [], [], [], QPoptions);
    subQPiter(k) = output.iterations;
    subQPfval(k) = fval;
    plambda = lambda_qp.ineqlin - lambda;
     
    % Backtraking line search method
    sigma = 1; % Hessian positive definite
    rho = 0.2; % rho in (0,1)
    mu = (g'*incremental + sigma/2*incremental'*B*incremental)/((1-rho)*norm(cineq,1));
    alpha = 1; % full step
    tau = 0.8; % contraction factor
    
    % check sufficient decrease condition
    fn = objfun(x(:,k)+alpha*incremental); fcalls=fcalls+1;
    cn = constraint(x(:,k)+alpha*incremental);
    eta = 0.2;
    while fn+mu*norm(cn,1) > f+mu*norm(cineq,1)+eta*alpha*(g'*incremental-mu*norm(cineq,1))
        alpha = tau*alpha; % update alpha 
        fn = objfun(x(:,k)+alpha*incremental); fcalls=fcalls+1;
        cn = constraint(x(:,k)+alpha*incremental);
    end
    f = fn;
    
    % update lambda value
    lambda = lambda + alpha*plambda;
    gL = g - dcineq'*lambda;   %gradient of Lagrangian
    
    % evaluate with new x(k+1)
    x(:,k+1) = x(:,k) + alpha*incremental;
    [f,g,~] = objfun(x(:,k+1)); fcalls=fcalls+1;
    [cineq,dcineq,~,~] = constraint(x(:,k+1));
    
    % compute gradient of lagrangian with new x(k+1)
    gL_new = g - dcineq'*lambda;
    
    % update Hessian matrix by modified BFGS procedure
    p = x(:,k+1) - x(:,k);
    q = gL_new - gL;
    B = update_hessian(B, p, q);
    
    F = [gL_new; cineq];
    k = k+1; 
    
end
info.iter = k;
info.seq = x;
info.fcalls = fcalls;
info.subQPiter = subQPiter;
info.subQPfval = subQPfval;
end



