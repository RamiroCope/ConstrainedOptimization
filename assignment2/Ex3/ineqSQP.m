function [x,lambda,info] = ineqSQP(x,lambda,tol,maxit)

% Initializing info storage variables
subQPiter = zeros(1,maxit); 
subQPfval = zeros(1,maxit);
k=1;
fcalls=1;
n=length(x);
x(:,k) = x; % storing solution iterates

% Initial evaluation of the problem and constraints
[~,g,~] = objfun(x);                 %gradient of original function
[cineq,dcineq,~,~] = constraint(x);  %constraints and its gradient

B = eye(2);                          %Hessian of L to be aproximated
F = [g - dcineq'*lambda; cineq];     %1st Order KKT conditions

QPoptions = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'none');

% Checking max iterations and KKT conditions
while ((k < maxit) && (norm(F(1:2),'inf') > tol) )    
    
    [~,g,~] = objfun(x(:,k)); fcalls=fcalls+1;
    [cineq,dcineq,~,~] = constraint(x(:,k)); %c(x) and A(x) in report notation
    
    % Defining variables
    gL = g - dcineq'*lambda;    %gradient of Lagrangian
    
    % Solving QP subproblem: to get Newton's descent direction
    [z, fval,~,output, lambda] = quadprog(B,g,-1*(dcineq), cineq, [], [], [], [], [], QPoptions);
    % Storing SubQP info
    subQPiter(k) = output.iterations;
    subQPfval(k) = fval;
    
    % Obtain next point in iteration sequence
    incremental = z(1:n); % p in report's notation
    x(:,k+1) = x(:,k) + incremental;
    
    % Update lambda
    lambda = lambda.ineqlin; %only when using quadprog
    
    % Evaluate objective function and constraints with new point, x(k+1)
    [~,g,~] = objfun(x(:,k+1)); fcalls=fcalls+1;
    [cineq, dcineq] = constraint(x(:,k+1));
    
    % Compute gradient of lagrangian with new x(k+1)
    gL_new = g - dcineq'*lambda;
    
    % Update Hessian matrix by damped BFGS procedure
    p = x(:,k+1) - x(:,k);
    q = gL_new - gL;
    B = update_hessian(B, p, q);
    
    % Update KKT conditions
    F = [gL_new; cineq];
    k = k+1; 
    
    
end
% storing convergence stats
info.iter   = k;
info.seq    = x;
info.fcalls = fcalls;
info.subQPiter = subQPiter;
info.subQPfval = subQPfval;
end


