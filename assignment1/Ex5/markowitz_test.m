G = [   2.3     .93     .62     .74     -.23    0;
        0.93    1.4     0.22    0.56    0.26    0;
        .62     .22     1.8     .78     -0.27   0;
        .74     .56     .78     3.4     -0.56   0;
        -0.23   0.26    -0.27   -0.56   2.6     0;
        0       0       0       0       0       0;
];

r_f = 2;
returns = [15.1, 12.5, 14.7, 9.02, 17.68, r_f]; 
R = 15;

% Quadprog form
H = G; % covariance matrix
g = zeros(6,1);
A = [ones(1,6); returns]';
b = [1; R]; % target R
C = eye(6);
d = zeros(6,1);

% just testing
x = 0*ones(6,1);
y = ones(2,1);
z = ones(6,1);
s = ones(6,1);

tic
[x_ours, xsteps] = pcipa( H, g, A, b, C, d, x, y, z, s );
toc
% x_qp = 0.1655 0.1365 0.3115 0.0266 0.3352 0.0247


%% Setting up variables according to quadprog() documentation
% This constraint allows us to specify our target return
A1 = returns;
b1 = R;

% This constraints makes portfolio weights (x) equal to 1
A2 = [1,1,1,1,1,1];
b2 = 1;

% Combining the above 2 equality constraints into 1; to make quadprog happy 
Aeq = [A1; A2];
beq = [b1; b2];

% Inequality constraint to disallow short-selling
Aineq = -eye(6);
bineq = zeros(6,1); 

% Lower and Upper Bounds on x
LB = zeros(6,1); 
UB = [1; 1; 1; 1; 1; 1];

%optimal port. weights found by quadprog
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIterations', 200, 'Display', 'none');
tic
[x_qp, fval, exitflag, output, lambda] = quadprog( H, [], Aineq, bineq, Aeq, beq, LB, UB, ones(6,1) ,options ); 
toc

disp((x_qp - x_ours)')
disp(xsteps)

%% Plotting Markowitz
G = [   2.3     .93     .62     .74     -.23    0;
        0.93    1.4     0.22    0.56    0.26    0;
        .62     .22     1.8     .78     -0.27   0;
        .74     .56     .78     3.4     -0.56   0;
        -0.23   0.26    -0.27   -0.56   2.6     0;
        0       0       0       0       0       0;
];
asset_var2 = diag(G);
asset_std2 = sqrt(asset_var2);
r_f = 2;
R = 15;
returns = [15.1, 12.5, 14.7, 9.02, 17.68, r_f]; 

% Quadprog form
H = G; % covariance matrix
g = zeros(6,1);
A = [ones(1,6); returns]';
b = [1; R]; % target R
C = eye(6);
d = zeros(6,1);
% We iterate on the seq of target returns
xset = zeros(6,100);
b1 = linspace(2,25,100); %seq of target returns
for i=1:100
    [x1,x1steps] = pcipa( H, g, A, [1; b1(i)], C, d, x, y, z, s );
    xset(: ,i) = x1;
end
figure;
subplot(1,2,1)
hold on;
xsum = sum(xset,1);
plot(b1, xset); 
title('Optimal Portfolios');
legend('asset 1', 'asset 2', 'asset 3', 'asset 4', 'asset 5', 'rfs');
xlabel('Target Return (%)');
ylabel('Asset Weights'); hold off; 
axis([2 17.5 0 1]);

% Obtaining Variance of each optimal porfolio in xset to plot eff frontier

%Calculating variance of all optimal porfolios (found in xset)
opt_port_var = zeros(100);
opt_port_std = zeros(100);
opt_port_return= zeros(100); %should be the same as b1
for i=1:length(xset)
    opt_port_var(i) = xset(:,i)'*G*xset(:,i);
    opt_port_std(i) = sqrt(opt_port_var(i));
    
    %Calculating returns to verify with b1 that they're correct
    opt_port_return(i) = returns*xset(:,i);
end

subplot(1,2,2)
plot(opt_port_return, opt_port_std)
title('Efficient Frontier');hold on;
scatter(returns, asset_std2)
plot(15, port_std, 'k +' );
xlabel('Return (%)');
ylabel('Risk (std)');  hold off;   
xlim([0 17.7])
