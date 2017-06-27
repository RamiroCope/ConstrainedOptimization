clear all;
clc; 

% RUN SECTION BY SECTION : variable names might repeat
%% Given Variables: Asset returns, target portfolio return, assets covariance
returns = [15.1, 12.5, 14.7, 9.02, 17.68]; 
R = 10; 
covariance = [2.3 .93 .62 .74 -.23;
              0.93 1.4 0.22 0.56 0.26;
              .62  .22 1.8  .78  -0.27;
              .74 .56 .78 3.4 -0.56;
              -0.23 0.26 -0.27 -0.56 2.6];  
asset_var = diag(covariance);    % variance of individual asssets
asset_std = sqrt(asset_var);     % std. dev of individual asssets

%% Setting up variables according to quadprog() documentation
H = covariance;
f = []; %f is not used, hence [] the empty brace

% This constraint allows us to specify our target return
A1 = returns;
b1 = R;

% This constraints makes portfolio weights (x) equal to 1
A2 = [1,1,1,1,1];
b2 = 1;

% Combining the above 2 equality constraints into 1; to make quadprog happy 
Aeq = [A1; A2];
beq = [b1; b2];

% Inequality constraint to disallow short-selling
Aineq = -eye(5);
bineq = zeros(5,1); 

% Lower and Upper Bounds on x
LB = zeros(5,1); 
UB = [1; 1; 1; 1; 1];
%% Solving
%optimal port. weights found by quadprog
x = quadprog( H, f, Aineq, bineq, Aeq, beq, LB, UB ); 

%Calculating portfolio variance (i.e., inherent risk of portfolio)
port_var = x'*covariance*x; % --> 2.0923
port_std = sqrt(port_var);  % --> 1.4465 

%Calculating porfolios expected return (should be 10%)
port_expected_return = returns*x; 

%% Efficient Frontier (without Risk Free Security)
% We minimize the variance for a sequence of target returns
b1 = linspace(2,25,1000); %seq of target returns

% We iterate on the seq of target returns
xset = zeros(5,1000);
for i=1:1000
    beq = [b1(i); b2];
    xset(: ,i) = quadprog( H, f, Aineq, bineq, Aeq, beq, LB, UB);
end
%% Plot
plot(b1, xset, 'LineWidth',5); hold on;
title('Optimal Portfolios');
legend('asset 1', 'asset 2', 'asset 3', 'asset 4', 'asset 5');
xlabel('Target Return (%)');
ylabel('Asset Weights');
axis([9 17.5 0 1]); hold off; 

%% Obtaining Variance of each optimal porfolio in xset to plot eff frontier

%Calculating variance of all optimal porfolios (found in xset)
opt_port_var = zeros(1000,1);
opt_port_std = zeros(1000,1);
opt_port_return= zeros(1000,1); %should be the same as b1
for i=1:length(xset)
    opt_port_var(i) = xset(:,i)'*covariance*xset(:,i);
    opt_port_std(i) = sqrt(opt_port_var(i));
    
    %Calculating returns to verify with b1 that they're correct
    opt_port_return(i) = returns*xset(:,i);
end

%% Plot
plot(opt_port_return(305:682), opt_port_std(305:682),  'LineWidth',5)
title('Efficient Frontier');hold on;
xlabel('Return (%)');
ylabel('Risk (std)'); 
scatter(returns, asset_std, 'filled');
axis([0 20 0 2]); hold off;

%% Adding Risk Free Security

returns = [15.1, 12.5, 14.7, 9.02, 17.68, 2]; 
R = 15; 
covariance = [2.3 .93 .62 .74 -.23, 0;
              0.93 1.4 0.22 0.56 0.26, 0;
              .62  .22 1.8  .78  -0.27, 0;
              .74 .56 .78 3.4 -0.56, 0;
              -0.23 0.26 -0.27 -0.56 2.6, 0;          
               0    0   0   0    0  0;]; 
asset_var2 = diag(covariance);
asset_std2 = sqrt(asset_var2);

           
%% Setting up variables according to quadprog() documentation
H = covariance;
f = []; %f is not used, hence [] the empty brace

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
%% Solving
%optimal port. weights found by quadprog
x = quadprog( H, f, Aineq, bineq, Aeq, beq, LB, UB ); 

%Calculating portfolio variance (i.e., inherent risk of portfolio)
port_var = x'*covariance*x; 
port_std = sqrt(port_var);  

%Calculating porfolios expected return (should be 10%)
port_expected_return = returns*x; 

%% Efficient Frontier (with Risk Free Security)
% We minimize the variance for a sequence of target returns
b1 = linspace(2,25,1000); %seq of target returns

% We iterate on the seq of target returns
xset = zeros(6,1000);
for i=1:1000
    beq = [b1(i); b2];
    xset(: ,i) = quadprog( H, f, Aineq, bineq, Aeq, beq, LB, UB);
end

%% Plot
plot(b1, xset, 'LineWidth',5); hold on;
title('Optimal Portfolios');
legend('asset 1', 'asset 2', 'asset 3', 'asset 4', 'asset 5', 'rfs');
xlabel('Target Return (%)');
ylabel('Asset Weights');
axis([2 17.5 0 1]);
vline(15, 'k+: '); hold off; 

%% Obtaining Variance of each optimal porfolio in xset to plot eff frontier

%Calculating variance of all optimal porfolios (found in xset)
opt_port_var = zeros(1000,1);
opt_port_std = zeros(1000,1);
opt_port_return= zeros(1000,1); %should be the same as b1
for i=1:length(xset)
    opt_port_var(i) = xset(:,i)'*covariance*xset(:,i);
    opt_port_std(i) = sqrt(opt_port_var(i));
    
    %Calculating returns to verify with b1 that they're correct
    opt_port_return(i) = returns*xset(:,i);
end

%% Plot
plot(opt_port_return(1:685), opt_port_std(1:685), 'LineWidth',5)
title('Efficient Frontier with Risk Free Asset');hold on;
scatter(returns, asset_std2, 'filled')
xlabel('Return (%)');
ylabel('Risk (std)');
plot(15, port_std, 'k o' , 'MarkerSize',20);
hold off;           
