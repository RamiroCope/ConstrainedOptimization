close all

%% Rewrite the problem in the quadratic form
H = [ 2 0;
      0 2 ];
g = [-2; -5];
% inequality constraints
C = [-1 2; 1 2; 1 -2];
d = [2; 6; 2];

%% Calculate the minimum using the implemented algorith
x = [0;0];
y = 0.5*ones(0,1);
z = 0.5*ones(3,1);
s = 0.5*ones(3,1);
A = zeros(2,0);
b = zeros(0,1);

%% Check with Matlab's `quadprog`
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'MaxIterations', 200);
tic
xqp = quadprog( H, g, C, d, [], [], zeros(2,1), [], [], options ); 
toc

%% Our implementation of interior point algorithm
%  (arguments need to be a bit modified compared to Matlab's quadprog 
%   to match the definition's on the slides)
tic
[xopt, xsteps] = pcipa(H,g,A,b,-1*C',-1*d,x,y,z,s);
toc
%xsteps = xsteps(:,all(~isnan(xsteps))); % remove NaN

%% Plot the function and the constraints
[X1,X2] = meshgrid(0:0.01:6, 0:0.01:3);
Q = X1.^2 - 2.*X1 + X2.^2 - 5.*X2 + 6.25;

figure
hold on
contourf(X1,X2,Q,50), colorbar
%contour(X1,X2,Q,20,'ShowText','on','LineWidth',2), colorbar

scatter(xopt(1), xopt(2))
plot(xsteps(1,:), xsteps(2,:),'-','LineWidth',2,'Color','red')

patch('Faces',[1 2 3], 'Vertices', [2 0; 6 2; 6 0], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [0 1; 6 4; 0 4], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [0 3; 6 0; 6 3], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')

xlim([0 6]); ylim([0 3]);
xlabel('x_1'); ylabel('x_2');

%% Check the results
format long
abs(xopt - xqp)'
format short