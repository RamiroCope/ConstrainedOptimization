function [ x_r, xiter ] = pcipa( H, g, A, b, C, d, x, y, z, s )
% Predictor-Corrector Interior-Point Algorithm

mc = length(s);
S = diag(s); 
Z = diag(z); 
e = ones(mc, 1); 
% a = 0.999;
% ahat = 1;


%Initial step
rL = H*x+g-A*y-C*z;
rA = b-A'*x;
rC = s+d-C'*x;
rsz = s.*z;
Hhat = H+C*(inv(S)*Z)*C';
HA = [Hhat -A; -A' zeros(size(A,2))];
rhatL = rL - C*((inv(S)*Z))*(rC-inv(Z)*rsz);
[L,D,p] = ldl(HA,'vector');
bl = [-rhatL; -rA];
xx(p) = L'\(D\(L\bl(p)));
dxaff = xx(1:end-size(A,2));
dyaff = xx(end-size(A,2)+1:end);
dzaff = -((inv(S)*Z))*C'*dxaff'+((inv(S)*Z))*(rC-inv(Z)*rsz); 
dsaff = -inv(Z)*rsz-(inv(Z)*S)*dzaff;
z = max(1, abs(z+dzaff));
s = max(1, abs(s+dsaff));


% Compute residuals
rL = H*x+g-A*y-C*z;
rA = b-A'*x;
rC = s+d-C'*x;
rsz = s.*z;
mu = (z'*s)/mc;


iter = 0; % current iteration number
maxiter = 200;
AdvancedConvCheck = true;
xiter = x;
% epsilon = 1e-8;
% mu0 = 1e-8;
eta = 0.995; % or 0.95 not sure if mistake or on purpose (slides say 0.995)
epsL=1e-8; epsA=1e-8; epsC=1e-8; epsmu=1e-8;
if size(A,2) == 0
    epsL = 1e-2;
end
Stop = check_convergence();

while ~Stop && (iter <= maxiter)
    iter = iter + 1;
    % Precompute
    S=diag(s); Z=diag(z);
    invStimesZ = diag(z./s);
    
    % Compute Hhat
    Hhat = H+C*(invStimesZ)*C';
    % Compute LDL-factorization of [Hhat -A; -A' 0]
    HA = [Hhat -A; -A' zeros(size(A,2))];
    
    % ## Affine Direction
    rhatL = rL - C*((z./s).*(rC-rsz./z));% rL - C*(invStimesZ)*(rC-inv(Z)*rsz);
    % Solve HA*[dxaff; dyaff] = [-rhatL; -rA] using LDL (Ax=b)
    [L,D,p] = ldl(HA,'vector');
    bl = [-rhatL; -rA];
    xx(p) = L'\(D\(L\bl(p)));
    dxaff = xx(1:end-size(A,2));
    dyaff = xx(end-size(A,2)+1:end);
    % Compute dzaff
    dzaff = -(invStimesZ)*C'*dxaff' + (z./s).*(rC-rsz./z);%-(invStimesZ)*C'*dxaff' + (invStimesZ)*(rC-inv(Z)*rsz);%%(C'*(-(z./s).*dxaff'))
    % Compute dsaff
    dsaff = -((rsz+s.*dzaff)./z);%-inv(Z)*rsz - inv(Z)*S*dzaff;%%
    
    % ### Compute the larget aaff 
    %     s.t. z+aaff*dzaff >= 0 and s+aaff*dsaff >= 0
    aaff = 1.0;
    dzs = find(dzaff < 0.0); 
    dss = find(dsaff < 0.0);
    if (~isempty(dzs))
        aaff = min(1, min(-z(dzs)./dzaff(dzs))); 
    end
    if (~isempty(dss))
        aaff = min(aaff, min(-s(dss)./dsaff(dss))); 
    end
    %disp(sprintf('Alpha affine is %f', aaff));
    
    % ## Duality gap and centering parameter
    % affine duality gap
    muaff = ((z+aaff*dzaff)'*(s+aaff*dsaff))/mc;
    % centering param
    sigma = (muaff/mu)^3;
    
    % ## Affine-centering-correction direction
    rhatsz = rsz + dsaff.*dzaff - sigma*mu*e; 
    rhatL = rL - C*((z./s).*(rC-rhatsz./z));%rhatL = rL - C*(invStimesZ)*(rC-inv(Z)*rhatsz);
    % Solve the system again
    bl = [-rhatL; -rA];
    xx(p) = L'\(D\(L\bl(p)));
    dx = xx(1:end-size(A,2));
    dy = xx(end-size(A,2)+1:end); 
    % Compute dzaff
    dz = -(invStimesZ)*C'*dxaff' + (z./s).*(rC-rhatsz./z);%-(invStimesZ)*C'*dxaff' + (invStimesZ)*(rC-inv(Z)*rhatsz);%%(C'*(-(z./s).*dx'))
    % Compute dsaff
    ds = -((rhatsz+s.*dz)./z);%-inv(Z)*rsz - inv(Z)*S*dzaff;%%
    
    % ### Compute the largest alpha
    %     s.t. z+alpha*dz >= 0 and s+alpha*ds >= 0    
    alpha = 1;
    dzs = find(dz < 0.0); 
    dss = find(ds < 0.0);
    if (~isempty(dzs))
        alpha = min(1, min(-z(dzs)./dz(dzs))); 
    end 
    if (~isempty(dss))
        alpha = min(alpha, min(-s(dss)./ds(dss))); 
    end
    %disp(sprintf('Alpha is %f', alpha));
    
    % ## Update iteration
    ahat = eta*(alpha); 
    x = x+ahat*dx';
    xiter = [xiter x];
    y = y+ahat*dy';
    z = z+ahat*dz;
    s = s+ahat*ds;
    S = diag(s); 
    Z = diag(z);
    % residuals
    rL = H*x+g-A*y-C*z;
    rA = b-A'*x;
    rC = s+d-C'*x;
    rsz = s.*z;
    mu = (z'*s)/mc;
    
    % ## Check for convergence conditions
    Stop = check_convergence();
    
end % end while

x_r = x; y_r = y; s_r = s; z_r = z;

%%% Helper functions %%%
    function [converges] = check_convergence()
        AdvancedConvCheck = true;
        
        %for markowitz when close to singularity for the last iterations
        if (isnan(max(H*x+g-A*y-C*z))) || (max(rL) > 1000000)
            converges = true;
            return
        end
        
        if AdvancedConvCheck
            arL = (max(H*x+g-A*y-C*z) <= 1e-2*max(1,norm([H g A C], 'inf')));
            if size(rA,1) == 0
                arA = true;
            else
                arA = (max(b-A'*x) <= 1e-8*max(1,norm([A' b], 'inf')));
            end
            arC = (max(d+s-C'*x) <= 1e-8*max(1,norm([eye(size(C,2)) d C'], 'inf')));
            amu = (mu <= 1e-8*10^(-2)*1e-4);
            %disp(sprintf('\nChecking covergence %d (rL=%f, rA=%f, rC=%f, mu=%f)\n', iter, max(H*x+g-A*y-C*z), max(b-A'*x), max(d+s-C'*x), (mu)));

            converges = arL && arA && arC && amu;
        else
            %disp(sprintf('\nChecking covergence %d (rL=%f, rA=%f, rC=%f, mu=%f)\n', iter, norm(rL), norm(rA), norm(rC), abs(mu)));
            if size(rA,1) == 0
                convrA = true;
            else
                convrA = norm(rA)<=epsA;
            end
            
            converges = (norm(rL)<=epsL) && (convrA) && (norm(rC)<=epsC) && (abs(mu)<=epsmu);
        end
    end
    
end

