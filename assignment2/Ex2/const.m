% The following function evaluates the constraints in exercise 18.3
% in N&W and computes the gradient and the Hessian at a specific point.

function [c,dc,d2c1,d2c2,d2c3] = const(x)

c = [x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-10;
     x(2)*x(3)-5*x(4)*x(5);
     x(1)^3+x(2)^3+1];

dc = [2*x(1) 0 3*x(1)^2;
      2*x(2) x(3) 3*x(2)^2;
      2*x(3) x(2) 0;
      2*x(4) -5*x(5) 0;
      2*x(5) -5*x(4) 0]';
  
d2c1 = diag(2*ones(5));
d2c2 = zeros(5);
d2c2(2,3) = 1; d2c2(3,2) = 1; d2c2(4,5) = -5; d2c2(5,4) = -5;
d2c3 = zeros(5);
d2c3(1,1) = 6*x(1);
d2c3(2,2) = 6*x(2);

end