 
function [f,g,H] = objfun(x)

%Evaluates the objective function, its first order and second order
%derivatives for a supplied x vector. 

x1=x(1);
x2=x(2);
f = (x1.^2+x2-11).^2 + (x1+x2.^2-7).^2;
g = [4*(x1.^2+x2-11).*x1+2*(x1+x2.^2-7); 2*(x1.^2+x2-11)+4*(x1+x2.^2-7).*x2];
H = [12*x1^2+4*x2-42,       4*(x1+x2);
      4*(x1+x2),      4*x1+12*x2.^2-26]; 


%{
With Symbolic Tool
    % dummy x for testing, and extracting supplied x
    x = [1;1];
    x1 = x(1,1); x2 = x(2,1);
    
    % creating symbolic math: obj function
    syms x_1 x_2 
    %func = exp(x_1*x_2*x_3*x_4*x_5) - 0.5*(x_1^3 + x_2^3 +1)^2
    func = (x_1^2 + x_2 - 11)^2 + (x_
    x_1=x2; x_2=x2; 
    
    % returning solution to obj function 
    f = eval(func);
    
    % creating symbolic math: objective function
    syms x_1 x_2 x_3 x_4 x_5
    df1 = diff(func, x_1);
    df2 = diff(func, x_2);
    df3 = diff(func, x_3);
    df4 = diff(func, x_4);
    df5 = diff(func, x_5);
    
    % returning 1st order solution
    df = zeros(5,1);
    x_1=x2; x_2=x2; x_3=x3; x_4=x4; x_5=x5;
    df(1,1) = eval(df1);
    df(2,1) = eval(df2);
    df(3,1) = eval(df3);
    df(4,1) = eval(df4);
    df(5,1) = eval(df5);
      
      %}