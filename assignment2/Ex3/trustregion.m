function [ x,lambda,info ] = trustregion( x,lambda,tol,maxit )
%TRUSTREGION Summary of this function goes here
%   Detailed explanation goes here

% Initializing info storage variables
subQPiter = zeros(1,maxit); 
subQPfval = zeros(1,maxit);
fcalls=1;
k=1;
n=length(x);
x(:,k) = x; % storing solution iterates hg

% Padding Hessian approx. with zeros to
% accommodate trust-region formulation 
zero = [0 0; 0 0];
B = [eye(2) zero;
    zero   zero];
mu = 40; dk = 1;
[~,g,~] = objfun(x); [cineq,dcineq,~,~] = constraint(x);
g = [g; mu; mu];
F = [g(1:2) - dcineq'*lambda; cineq];

QPoptions = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'none');
      
while ((k < maxit) && (norm(F(1:2),'inf') > tol) )   %&& (norm(cineq,'inf') < tol)) 
   
    [~,g,~] = objfun(x(:,k)); fcalls=fcalls+1;
    g = [g; mu; mu];
    [cineq,dcineq,~,~] = constraint(x(:,k));
    cineq=[cineq;0;0];
    dcineq=[dcineq zero;
            zero   zero];
    % Defining variables
    gL = g(1:2) - dcineq(1:2,1:2)'*lambda;    %gradient of Lagrangian
    
    % Solving QP subproblem
    LB = [-dk -dk 0 0]';
    UB = [dk dk inf inf]';
    [z, fval,~,output, lambda] = quadprog(B,g,-1*(dcineq), cineq, [], [], LB, UB, [], QPoptions);
    subQPiter(k) = output.iterations;
    subQPfval(k) = fval;
    
    incremental = z(1:n);
    x(:,k+1) = x(:,k) + incremental;
    t = z(n+1:end);
    lambda = lambda.ineqlin(1:2); %only when using quadprog
    
    % evaluate with new x(k+1)
    [~,g,~] = objfun(x(:,k+1)); fcalls=fcalls+1;
    [cineq, dcineq] = constraint(x(:,k+1));
    
    % compute gradient of lagrangian with new x(k+1)
    gL_new = g - dcineq'*lambda;
    
    % update Hessian matrix by modified BFGS procedure
    p = x(:,k+1) - x(:,k);
    q = gL_new - gL;
    B(1:2,1:2) = update_hessian(B(1:2,1:2), p, q);
    
    F = [gL_new; cineq];
    k = k+1; 
    
end
info.iter = k;
info.seq = x;
info.fcalls = fcalls;
info.subQPiter = subQPiter;
info.subQPfval = subQPfval;
end









