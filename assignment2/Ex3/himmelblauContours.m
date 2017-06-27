    function himmelblauContours

    Himmelblau = @(x1,x2)(x1.^2+x2-11).^2 + (x1+x2.^2-7).^2;

    % Compute the function values
    x1 = linspace(-10,10, 400); % linspace(-5,5, 100);
    x2 = linspace(-10,10, 400); % linspace(-5,5, 100);
    [X1,X2]=meshgrid(x1,x2);
    F = Himmelblau(X1,X2);
    

    % Make contour plot
    v = [0:2:10 10:10:100 100:20:200]; hold on
    contour(X1,X2,F,v,'linewidth',2);
    
  

    xlabel('x_1','Fontsize',14)
    ylabel('x_2','Fontsize',14)
    set(gca,'fontsize',14);
    set(gcf, 'Color', 'w');
    colorbar
    axis image


    %Xmin = [3, 2; -3.77931, -3.28319; -2.80512, 3.13131; 3.58443, -1.84813];
    %Xsaddle = [-3.07303, -0.081353; -0.127961, -1.95371; 0.0866775, 2.88425; 3.38515, 0.0738519];
    %Xmax = [-0.270845, -0.923039];

    % Plot stationary points
    %plot(Xmin(:,1),Xmin(:,2),'b.','markersize',40)
    %plot(Xsaddle(:,1),Xsaddle(:,2),'k.','markersize',40)
    %plot(Xmax(:,1), Xmax(:,2),'r.','markersize',40)
    %set(gcf, 'Color', 'w');
    %export_fig '../img/himmelblauContours.png'