%% Plot the function and the constraints
[X1,X2] = meshgrid(0:0.01:6, 0:0.01:3);
Q = X1.^2 - 2.*X1 + X2.^2 - 5.*X2 + 6.25;

figure
hold on
contourf(X1,X2,Q,50), colorbar

scatter(1.4, 1.7)

xsteps = [2 0; 1 0; 1 1.5; 1.4 1.7];
plot(xsteps(:,1), xsteps(:,2),'x--','LineWidth',2,'Color','red')

patch('Faces',[1 2 3], 'Vertices', [2 0; 6 2; 6 0], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [0 1; 6 4; 0 4], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [0 3; 6 0; 6 3], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')

xlim([0 6]); ylim([0 3]);
xlabel('x_1'); ylabel('x_2');
