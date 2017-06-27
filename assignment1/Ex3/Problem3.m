clear;
close all;

H = [2 0; 0 2];
g = [-2 ; -5];
A = [1 -1 -1 1 0;-2 -2 2 0 1]';
b = [-2; -6; -2; 0; 0];
%% 3.1 : Contour plot with constraints
x=-2:0.05:5;
y=-2:0.05:5;
[X Y] = meshgrid(x,y);

q=(X-1).^2+(Y-2.5).^2;
contourf(X,Y,q,linspace(0,30,50))
axis([-1 4.5 -1 3])
colorbar;

yc1=x/2+1;
yc2=-x/2+3;
yc3=x/2-1;
xc4=zeros(size(x));
yc5=zeros(size(y));
hold on
%     plot(x,yc1);
%     plot(x,yc2);
%     plot(x,yc3);
%     plot(xc4,y);
%     plot(x,yc5);
    fill([x x(end) x(1)],[yc1 y(end) y(end)],[0.7 0.7 0.7],'facealpha',0.3, 'facecolor', 'black')
    fill([x x(end) x(1)],[yc2 y(end) y(end)],[0.7 0.7 0.7],'facealpha',0.3, 'facecolor', 'black')
    fill([x x(end) x(1)],[yc3 y(1) y(1)],[0.7 0.7 0.7],'facealpha',0.3, 'facecolor', 'black')
    fill([xc4 x(1) x(1)],[y y(end) y(1)],[0.7 0.7 0.7],'facealpha',0.3, 'facecolor', 'black')
    fill([x x(end) x(1)],[yc5 y(1) y(1)],[0.7 0.7 0.7],'facealpha',0.3, 'facecolor', 'black')

X = quadprog(H,g,-A,-b);
hold on
plot(X(1),X(2),'ro', 'MarkerSize', 6)

%% 3.8 : Find a feasible point using linprog

[m,n] = size(A);
xe = [-1; -1];
ganma = sign(A*xe - b);
Af = [A diag(ganma); zeros(m,n) eye(m)];
bf = [b; zeros(m,1)]; 
f = [zeros(n,1); ones(m,1)];
z = linprog(f,-Af,-bf);
xk = z(1:n);
[x,lambda,info] = ActiveSetMethod(H,g,A,b,xk);

plot(info.xs(1,:),info.xs(2,:),'r--','LineWidth',1.5)

%% 3.10 %%  (16.47) Big M Penalty Method: L-inf norm


Hm = [2 0 0;...
      0 2 0;...
      0 0 1e-3];
gm = [g ; big_M];
xm = [-1 ; -1 ; v];
Aineq_m =   [A ones(size(A,1),1); 0 0 1;]; 
bineq_m = [b ; 0];

% Implementation:
[x,lambda,info] = ActiveSetMethod(Hm,gm,Aineq_m,bineq_m,xm);

% Summary
plot(info.xs(1,:),info.xs(2,:),'r--','LineWidth',1.5)



%% 3.10 %%  (16.47) Big M Penalty Method: L-1 norm

% NOTE: The objective function should be minimized wrt x AND eta
%       - eta should equal 0 if global solution is found
%       - set eta large enough such that all constraints are satisfied
%       - if during solving eta converges to positive number,
%         increase big_M until eta converges to 0
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%  Setting up Hyperparameters
eta = 200 ;
big_M = 1; 

%  Setting up Variables
Hm = [2 0 0 0 0 0 0;...
      0 2 0 0 0 0 0;...
      0 0 1 0 0 0 0;...
      0 0 0 1 0 0 0;...
      0 0 0 0 1 0 0;...
      0 0 0 0 0 1 0;...
      0 0 0 0 0 0 1];
gm = [g ; big_M; big_M; big_M; big_M; big_M;];
xm = [-1 ; -1 ; eta; eta; eta; eta; eta;];
one = ones(5,1);
unos = ones(6,1);
Aineq_m = [A one; 0 0 1;];
Aineq_m = [Aineq_m unos unos unos unos];
bineq_m = [b ; 0; 0; 0; 0; 0;];

% Implementation:
[x,lambda,info] = ActiveSetMethod(Hm,gm,Aineq_m,bineq_m,xm);

% Summary
x
info.iter 
info.ws
info.xs

plot(info.xs(1,:),info.xs(2,:),'r--','LineWidth',1.5)
% Plot with iterations
% RM: If anyone gets a chance, that would be great! : )
%  but maybe they're not needed
%  Right Here! 
%
%
%% (16.48) Missing 

%  Setting up Hyperparameters
v = 200 ;
big_M = 0.005; 

%  Setting up Variables
e =  [1 1 1 1 1];  % <-- FIGURE OUT HOW e GOES BELOW IN G AND H!
                   %     Hm, gm and xm might look different than below
Hm = [2 0 0 0 0 0 0;...
      0 2 0 0 0 0 0;...
      0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0];
gm = [g ; big_M; big_M; big_M; big_M; big_M;];
xm = [0 ; 0 ; v; v; v; v; v;];
one = ones(5,1);
unos = ones(6,1);
Aineq_m =   [A one; 0 0 1;]; 
Aineq_m = [Aineq_m unos unos unos unos];
bineq_m = [b ; 0; 0; 0; 0; 0;];

% Implementation:
[x,lambda,info] = ActiveSetMethod(Hm,gm,Aineq_m,bineq_m,xm);

% Summary
x
info.iter
info.ws
info.xs

% Plot with iterations
% RM: If anyone gets a chance, that would be great! : )
%     but maybe they're not needed
%  Right Here!               



