function [ B_new ] = update_hessian(B, p, q )
%UNTITLED4 Updates Hessian matrix by the modified BFGS procedure
%   Detailed explanation goes here

if p'*q >= 0.2*p'*B*p
    theta = 1;
    
elseif p'*q < 0.2*p'*B*p
    theta = 0.8*p'*B*p / (p'*B*p - p'*q);    
end

r = theta*q + (1-theta)*(B*p);

nom = B*(p*p')*B; 
denom = (p'*B*p);
B_new = B - (nom / denom) + r*r'/(p'*r); 

end


