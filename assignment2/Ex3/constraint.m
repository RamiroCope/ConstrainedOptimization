function [ cineq,dcineq,d2cineq1,d2cineq2 ] = constraint( x )
%CONSTRAINT evalutates the equality and inequality constraints for a given x.


cineq = [x(1)^2 + 4*x(1) - x(2) + 4;
          -4*x(1) + 10*x(2)];
  
dcineq = [ 2*x(1)+4   -4;
           -1         10]';

d2cineq1 = [2  0;
            0  0];

d2cineq2 = [0  0;
            0  0];
end

