clear 
close all
%% Equality Constrained SQP
% Initial values
x = [-1.8;1.7;1.9;-0.8;-0.8];
lambda = [1; 1; 1];
tol = 1e-3;
maxit = 200;
% solve nonlinear optimization problem
[x1,lambda1,info1] = local_SQP(@objfun,@const,x,lambda,tol,maxit);
% solve using a damped BFGS approximation to the Hessian of the Lagrangian
[x2,lambda2,info2] = BFGS_SQP(@objfun,@const,x,lambda,tol,maxit);
% solve using a damped BFGS approcimation and a line search method
[x3,lambda3,info3] = lineSearch_SQP(@objfun,@const,x,lambda,tol,maxit);
% just line search
[x4,lambda4,info4] = lineSearch_SQP2(@objfun,@const,x,lambda,tol,maxit);

%% Convergence semilogys
err1 = (info1.seq - x1);
err1_norm = diag(sqrt(err1'*err1));
semilogy(1:info1.iter,err1_norm,'-o')
hold on
err2 = (info2.seq - x2);
err2_norm = diag(sqrt(err2'*err2));
semilogy(1:info2.iter,err2_norm,'-o')
err3 = (info3.seq - x3);
err3_norm = diag(sqrt(err3'*err3));
semilogy(1:info3.iter,err3_norm,'-o')
err4 = (info4.seq - x4);
err4_norm = diag(sqrt(err4'*err4));
semilogy(1:info4.iter,err4_norm,'-o')
xlabel('iteration')
ylabel('2 norm of error')
legend('Local SQP','BFGS','BFGS + linesearch','linesearch')
axis([1 6 1e-5 1])
print('convergence2','-dpng')

%% Iteration sequence
%% Generate INFO 2 iteration sequence table
clear input;
clc;
% Add data, change as needed
input.data = info2.seq'; % make sure to have this correctly formatted
input.tablePositioning = 'h';
input.tableColLabels = {'$x_1$','$x_2$','$x_3$','$x_4$','$x_5$'};
input.tableColumnAlignment = 'l';
input.tableCaption = 'Iteration sequence for BFGS approximation';
input.tableLabel = 'ex2:bfgs'; % prepends "table" -> table:MyTableLabel
% Generate tex output
texOut = latexTable(input);

%% Generate INFO 3 iteration sequence table
clear input;
clc;
% Add data, change as needed
input.data = [info4.seq' [info4.alpha 0]']; % make sure to have this correctly formatted
input.tablePositioning = 'h';
input.tableColLabels = {'$x_1$','$x_2$','$x_3$','$x_4$','$x_5$','$\alpha$'};
input.tableColumnAlignment = 'l';
input.tableCaption = 'Iteration sequence for local SQP and line search';
input.tableLabel = 'ex2:localSQPls'; % prepends "table" -> table:MyTableLabel
% Generate tex output
texOut = latexTable(input);