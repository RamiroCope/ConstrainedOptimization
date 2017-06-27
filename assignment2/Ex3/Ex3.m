%% Plot the problem given in Ex3-
x2lim = 5; x2 = -x1lim:0.01:x2lim;
[X1,X2] = meshgrid(x1, x2);

% https://en.wikipedia.org/wiki/Himmelblau%27s_function
F = (X1.^2 + X2 - 11).^2 + (X1 + X2.^2 - 7).^2;

figure(1)
fprintf('Plotting...\n');
clf
axis([-x1lim x1lim -x1lim x2lim])
hold on
% Objective function
fprintf('Adding contour lines...\n');
[~,h] = contour(X1,X2,F, 'LevelList', ([
    0:2:10 10:10:100 100:20:200
]), 'LineWidth', 1.4); colorbar; % , 'ShowText', 'on'
% Constraints
fprintf('Adding constraints...\n');
plot(x1, x1.^2 + 4.*x1 + 4, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
plot(x1, 0.4*x1, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-')
x1f = -x1lim:0.01:x1lim;
fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none')
patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.4, 'EdgeColor', 'none')
hold off
%title('Himmelblau''s function: $f(\mathbf{x}) = (x_1^2 + x_2 - 11)^2 + (x_1 + x_2^2 - 7)^2$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$x_1$', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('$x_2$', 'FontSize', 14, 'Interpreter', 'latex')
set(gca, 'FontSize', 14)
set(gca, 'XTick', -x1lim:x1lim); set(gca, 'YTick', -x2lim:x2lim); 
grid on
fprintf('Done! Wait couple of seconds to see the rendered plot.\n');