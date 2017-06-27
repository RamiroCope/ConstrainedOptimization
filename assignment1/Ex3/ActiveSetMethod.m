function [x,lambda,info] = ActiveSetMethod(H,g,A,b,xk)
% This function solves strictly convex quadratic programs with inequality
% constraints. 
% The function returns:
% x: solution to the optimization problem.
% lambda: lagrange multipliers at the solution
% info.iter: number of iterations
% info.xs: iteration sequence
% info.ws: working set in each iteration

m = size(A,1);
lambda = zeros(m,1);

% Initial working set
w = (A*xk == 0);

stop = false;
k = 1;
while ~ stop
    ws(:,k) = w;
    xs(:,k) = xk;
    Ak = A(w,:);
    bk = b(w);
    gk = H*xk+g;
    if isempty(Ak) && isempty(bk)
        pk = EqualityQPSolver(H,gk,[],[]);
    else
        pk = EqualityQPSolver(H,gk,Ak',zeros(size(Ak,1),1));
    end
    if norm(pk) < 10e-4
        lambdak = Ak'\(H*xk +g);
        if sum(lambdak < 0) == 0
            x = xk;
            lambda(w) = lambdak;
            stop = true;
        else
           [~,idx] = min(lambdak);
           wk = ones(length(lambdak),1);
           wk(idx) = 0;
           w(w == 1) = wk;
        end
    else
        Ai = A(~w,:);
        bi = b(~w,:);
        i = Ai*pk < 0;
        alpha = ones(length(bi),1);
        alpha(i) = (bi(i)-Ai(i,:)*xk)./(Ai(i,:)*pk);
        [alphak,idx] = min(alpha);
        if alphak < 1
            xk = xk + alphak*pk;
            wk = zeros(length(alpha),1);
            wk(idx) = 1;
            w(w == 0) = wk;
        else
            xk = xk + pk;
        end
    end
    k = k+1;
end
info.iter = k;
info.xs = xs;
info.ws = ws;
end