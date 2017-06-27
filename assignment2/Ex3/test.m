clear all; close all; 


x0=[1,3 ; -5,5; 8,4]'; %3 starting points  
n=size(x0,2);    %number of starting points
lambda=[1 ; 1];  %initial lambda
tol=1e-5;
maxit=100; 
rep=30;         %number of repetitions for TicToc
T1=zeros(rep,n); %TicToc stats for ineqSQP
T2=zeros(rep,n); %TicToc stats for linesearch_ineqSQP
T3=zeros(rep,n); %TicToc stats for trustregion_ineqSQP
e1 = cell(1,n); e2 = cell(1,n); e3 = cell(1,n); %store errors for each algo
solution = [3;2];

% Structs to store info for each algorithm
ISQP{n}=[]; LS_ISQP{n}=[];  TR_ISQP{n}=[];
X1{n} = []; L1{n} = []; I1{n} = [];
X2{n} = []; L2{n} = []; I2{n} = [];
X3{n} = []; L3{n} = []; I3{n} = [];

%Contour plot
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours 
plot(x0(1,:),x0(2,:),'k.','markersize',70) %change desired starting point

% Iterating algorithms over 3 starting points
for j=1:n
    
    for i=1:1%rep
        
        % Algo 1: local IQP
        T=tic;
        [x,lambda,info] = ineqSQP(x0(:,j),lambda,tol,maxit);
        T1(i,j)=toc(T);
        X1{j} = x;  L1{j} = lambda;  I1{j} = info;
        ISQP{j} = struct('x',X1(j),'lambda',L1(j),'info',I1(j));

        % Algo 2: linesearch IQP
        T=tic;
        [x,lambda,info] = linesearch_ineqSQP(x0(:,j),lambda,tol,20);
        T2(i,j)=toc(T);
        X2{j} = x;  L2{j} = lambda;  I2{j} = info;
        LS_ISQP{j} = struct('x',X2(j),'lambda',L2(j),'info',I2(j));

        % Algo 3: trust-region IQP
        T=tic;
        [x,lambda,info] = trustregion(x0(:,j),lambda,tol,maxit);
        T3(i,j)=toc(T);
        X3{j} = x;  L3{j} = lambda;  I3{j} = info;
        TR_ISQP{j} = struct('x',X3(j),'lambda',L3(j),'info',I3(j));
    
    end%for i (repetitions)
    
        
    e1{j} = sqrt(sum( (([X1{j}(1,:)-3; X1{j}(2,:)-2]).^2) )); %algo1 error/iter
    e2{j} = sqrt(sum( (([X2{j}(1,:)-3; X2{j}(2,:)-2]).^2) )); %algo2 error/iter
    e3{j} = sqrt(sum( (([X3{j}(1,:)-3; X3{j}(2,:)-2]).^2) )); %algo3 error/iter
    
end%for j (starting points)

    h1 = plot(X1{1}(1,:),X1{1}(2,:),'k-o','linewidth',2);
    h2 = plot(X2{1}(1,:),X2{1}(2,:),'r-o','linewidth',2);
    h3 = plot(X3{1}(1,:),X3{1}(2,:),'b-o','linewidth',2);
    
    legend([h1,h2,h3],'local IQP','backtracking line search IQP', 'trust-region IQP','location','southeast')
    
    %feeding iterations to error_plot
    iter1 = [I1{1}.iter;I1{2}.iter;I1{3}.iter];
    iter2 = [I2{1}.iter;I2{2}.iter;I2{3}.iter];
    iter3 = [I3{1}.iter;I3{2}.iter;I3{3}.iter];
    error_plot(e1,e2,e3,iter1,iter2,iter3,n,maxit);
    
    % Timing stats for each algo
    T1=mean(T1); T2=mean(T2); T3=mean(T3);
    
    
    % Clean Workspace
%clear i I1 I2 I3 j L1 L2 L3 lambda T x X1 X2 X3 info 

%% Convergence Plots

% Local IQP Plot
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours
plot(x0(1,:),x0(2,:),'k.','markersize',70) 
h1 = plot(X1{1}(1,:),X1{1}(2,:),'m-o','linewidth',2);
h2 = plot(X1{2}(1,:),X1{2}(2,:),'r-o','linewidth',2);
h3 = plot(X1{3}(1,:),X1{3}(2,:),'b-o','linewidth',2);
legend([h1,h2,h3],'Starting Point 1','Starting Point 2', 'Starting Point 3','location','northeast')
% show constraints
x1lim = 14; x1 = -x1lim:0.01:x1lim;
x2lim = 10; x2 = -x1lim:0.01:x2lim;
x1f = -10:0.01:14;
%plot(x1, x1.^2 + 4.*x1 + 4, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
%plot(x1, 0.4*x1, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% set axis as needed
axis([-5 14 -5 8])
hold off

%% Line search IQP Plot
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours
plot(x0(1,:),x0(2,:),'k.','markersize',70) 
h1 = plot(X2{1}(1,:),X2{1}(2,:),'m-o','linewidth',2);
h2 = plot(X2{2}(1,:),X2{2}(2,:),'r-o','linewidth',2);
h3 = plot(X2{3}(1,:),X2{3}(2,:),'b-o','linewidth',2);
legend([h1,h2,h3],'Starting Point 1','Starting Point 2', 'Starting Point 3','location','northeast')
% show constraints
x1lim = 10; x1 = -x1lim:0.01:x1lim;
x2lim = 10; x2 = -x1lim:0.01:x2lim;
x1f = -10:0.01:10;
%plot(x1, x1.^2 + 4.*x1 + 4, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
%plot(x1, 0.4*x1, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% set axis as needed
axis([-5 8 -5 5])
hold off

%% Trust Region IQP Plot
fig = figure('Position', [100, 0, 1000,1000]);
clf;
hold on
himmelblauContours
plot(x0(1,:),x0(2,:),'k.','markersize',70) 
h1 = plot(X3{1}(1,:),X3{1}(2,:),'m-o','linewidth',2);
h2 = plot(X3{2}(1,:),X3{2}(2,:),'r-o','linewidth',2);
h3 = plot(X3{3}(1,:),X3{3}(2,:),'b-o','linewidth',2);
legend([h1,h2,h3],'Starting Point 1','Starting Point 2', 'Starting Point 3','location','northeast')
% show constraints
x1lim = 10; x1 = -x1lim:0.01:x1lim;
x2lim = 10; x2 = -x1lim:0.01:x2lim;
x1f = -10:0.01:10;
%plot(x1, x1.^2 + 4.*x1 + 4, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
%plot(x1, 0.4*x1, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% set axis as needed
axis([-5 8 -5 5])
hold off




%% LATEX TABLES (for iteration seq of all algorithms)

% LOCAL IQP
%Starting Point 1
matrix = I1{1}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
                  '12', '13', '14', '15', '16', '17', '18'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'localseq1.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
   
   %%
%Starting Point 2
matrix = I1{2}.seq(:,1:12)';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
                  '12'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'localseq2.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');

   %%
%Starting Point 3
matrix = I1{3}.seq(1:17)';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
                  '12', '13', '14', '15', '16', '17'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'localseq13.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');


%% LINESEARCH
%Starting Point 1
matrix = I2{1}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'lsseq1.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');

   %%
%Starting Point 2
matrix = I2{2}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'lsseq2.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
%%
%Starting Point 3
matrix = I2{3}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
                  '12', '13', '14', '15', '16'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'lsseq3.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');

%% TRUST REGION
%Starting Point 1
matrix = I3{1}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'trseq1.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
%%
%Starting Point 2
matrix = I3{2}.seq(:,1:9)';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'trseq2.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
%%
%Starting Point 3
matrix = I3{3}.seq';
   rowLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', ...
                  '12', '13', '14'};
   columnLabels = { '$x_1$', '$x_2$'};
   matrix2latex(matrix, 'trseq3.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'small');
   
%% LATEX Tables with TicToc, Fcalls, etc.
