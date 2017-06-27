% Implementation from http://etd.dtu.dk/thesis/220437/ep08_19.pdf
% used to check ours algo vs quadprog vs this
function [x_stop ,y_stop ,z_stop ,s_stop ,k] = test(x,y,z,s,G,g,C,d,A,b) 

eta = 0.95;
[mA,nA] = size(A); 
[mC,nC] = size(C);
e = ones(nC,1);
rL =G*x+g-A*y-C*z; 
rA = -A'*x + b;
rC = -C'*x + s + d;
rsz = s.*z;
mu = sum(z.*s)/nC;

k=0;
maxk = 200;
epsL=1e-10; epsA=1e-10; epsC=1e-10; epsmu=1e-10;


while (k<=maxk && norm(rL)>=epsL && norm(rA)>=epsA && norm(rC)>=epsC && abs(mu)>=epsmu)

    lhs = [G,-A,-C;
        -A' , spalloc(nA,nA,0) , spalloc(nA,nC,0);
        -C' , spalloc(nC, nA , 0 ) , sparse(-diag(s./z)) 
        ] ;
    [L,D,P] = ldl(lhs);
    rhs = [-rL; -rA; -rC+rsz ./ z ] ;
    dxyz_a = P*(L'\(D\(L\(P'*rhs))));
    dx_a = dxyz_a(1:length(x));
    dy_a = dxyz_a(length(x)+1:length(x)+length(y));
    dz_a = dxyz_a(length(x)+length(y)+1:length(x)+length(y)+length(z)); 
    ds_a = -((rsz+s.*dz_a)./z);
    %Compute alpha aff 
    alpha_a = 1;
    idx_z = find(dz_a<0); 
    if (isempty(idx_z)==0)
        alpha_a = min(alpha_a ,min(-z(idx_z)./dz_a(idx_z))); 
    end
    idx_s = find(ds_a<0); 
    if (isempty(idx_s)==0)
        alpha_a = min(alpha_a,min(-s(idx_s)./ds_a(idx_s))); 
    end
    disp(sprintf('[test] Alpha affine is %f', alpha_a));
    %Compute the affine duality gap
    mu_a = ((z+alpha_a*dz_a)'*(s+alpha_a*ds_a))/nC;
    %Compute the centering parameter 
    sigma = (mu_a/mu)^3;
    
    %Solve system
    rsz = rsz + ds_a.*dz_a - sigma*mu*e;
    rhs = [-rL;-rA;-rC+rsz ./ z ] ;
    dxyz = P*(L'\(D\(L\(P'*rhs))));
    dx = dxyz(1:length(x));
    dy = dxyz(length(x)+1:length(x)+length(y));
    dz = dxyz(length(x)+length(y)+1:length(x)+length(y)+length(z)); 
    ds = -((rsz+s.*dz)./z);
    %Compute alpha
    alpha = 1;
    idx_z = find(dz<0);
    if ( isempty ( idx_z )==0)
        alpha = min(alpha,min(-z(idx_z)./dz(idx_z))); 
    end
    idx_s = find(ds<0);
    if (isempty(idx_s)==0)
        alpha = min(alpha,min(-s(idx_s)./ds(idx_s))); 
    end
    disp(sprintf('[test] Alpha is %f', alpha));
    %Update x,z,s
    x = x + eta*alpha*dx;
    y = y + eta*alpha*dy; 
    z = z + eta*alpha*dz; 
    s =s + eta*alpha*ds;
    k = k+1;
    %Update rhs
    rL =G*x+g-A*y-C*z; 
    rA = -A'*x + b;
    rC = -C'*x + s + d;
    rsz = s.*z;
    mu = sum(z.*s)/nC;
end
%Output
x_stop = x; y_stop = y; z_stop=z; s_stop=s;


end

