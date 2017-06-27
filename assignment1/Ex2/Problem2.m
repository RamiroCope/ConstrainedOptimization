close all
clear
% This script will be used for testing the computing performance of 6
% different QP solvers

u = 0.2;
d0 = 1;

%% DENSE SOLVERS

%% LU solver

n = linspace(10,1000,10);
e1 = zeros(10,1); 

t = cputime;
for j = 1:1000
[x,lambda] = LUsolver(n(1),u,d0);
end
e1(1) = (cputime-t)/1000;

for i = 2:4
t = cputime;
for j = 1:100
[x,lambda] = LUsolver(n(i),u,d0);
end
e1(i) = (cputime -t)/100;
end

for i = 5:10
t = cputime;
[x,lambda] = LUsolver(n(i),u,d0);
e1(i) = (cputime -t);
end
%% LDL solver
e2 = zeros(10,1);

t = cputime;
for j = 1:1000
[x,lambda] = LDLsolver(n(1),u,d0);
end
e2(1) = (cputime-t)/1000;

for i = 2:4
t = cputime;
for j = 1:100
[x,lambda1] = LDLsolver(n(i),u,d0);
end
e2(i) = (cputime-t)/100;
end

for i = 5:10
t = cputime;
[x,lambda1] = LDLsolver(n(i),u,d0);
e2(i) = (cputime -t);
end

%% Null-Space solver
e3 = zeros(10,1);

t = cputime;
for j = 1:1000
[x,lambda] = NullSpaceSolver(n(1),u,d0);
end
e3(1) = (cputime-t)/1000;

for i = 2:4
t = cputime;
for j = 1:100
[x,lambda] = NullSpaceSolver(n(i),u,d0);
end
e3(i) = (cputime-t)/100;
end

for i = 5:10
t = cputime;
[x,lambda] = NullSpaceSolver(n(i),u,d0);
e3(i) = (cputime-t);
end
%% Range-Space solver
e4 = zeros(10,1);

t = cputime;
for j = 1:1000
[x,lambda] = RangeSpaceSolver(n(1),u,d0);
end
e4(1) = (cputime-t)/1000;

for i = 2:4
t = cputime;
for j = 1:100
[x,lambda] = RangeSpaceSolver(n(i),u,d0);
end
e4(i) = (cputime-t)/100;
end

for i = 4:10
t = cputime;
[x,lambda] = RangeSpaceSolver(n(i),u,d0);
e4(i) = (cputime-t);
end
%% CPUtime plot
subplot(2,1,1)
loglog(n,e1,n,e2,n,e3,n,e4,'linewidth',1.5)
hold on
plot(n,n.^2./10^6,'k--')
legend('LU','LDL','Null-Space','Range-Space','location','northwest')
subplot(2,1,2)
plot(n,e1,n,e2,n,e3,n,e4,'linewidth',1.5)
legend('LU','LDL','Null-Space','Range-Space','location','northwest')
%% SPARSE SOLVERS

%% Sparsity pattern
K = KKTmatrix(100,0.2,1);
figure
spy(K);
legend('non-zero elements','location','southeast')

%% LU sparse solver

e5 = zeros(10,1); 

t = cputime;
for j = 1:1000
[x,lambda] = LUsparseSolver(n(1),u,d0);
end
e5(1) = (cputime-t)/1000;

for i = 2:4
t = cputime;
for j = 1:100
[x,lambda] = LUsparseSolver(n(i),u,d0);
end
e5(i) = (cputime -t)/100;
end

for i = 5:10
t = cputime;
[x,lambda] = LDLsparseSolver(n(i),u,d0);
e5(i) = (cputime -t);
end

%% LDL sparse solver

e6 = zeros(10,1); 

t = cputime;
for j = 1:1000
[x,lambda] = LDLsparseSolver(n(1),u,d0);
end
e6(1) = (cputime-t)/1000;

for i = 2:5
t = cputime;
for j = 1:100
[x,lambda] = LDLsparseSolver(n(i),u,d0);
end
e6(i) = (cputime -t)/100;
end

for i = 6:10
t = cputime;
[x,lambda] = LDLsparseSolver(n(i),u,d0);
e6(i) = (cputime -t);
end

%% CPUtime plot
figure
subplot(2,1,1)
loglog(n,e3,n,e4,n,e5,n,e6,'linewidth',1.5)
legend('Null-Space','Range-Space','LU sparse','LDL sparse','location','northwest')
hold on
plot(n,n.^2./10^6,'k--')
subplot(2,1,2)
plot(n,e1,n,e2,n,e3,n,e4,n,e5,n,e6,'linewidth',1.5)
legend('LU','LDL','Null-Space','Range-Space','LU sparse','LDL sparse','location','northwest')