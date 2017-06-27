function [x,lambda] = EqualityQPdual(H,g,A,b)
AtH = A'/H;
lambda = (AtH*A)\(AtH*g+b);
x = H\(A*lambda-g);

end