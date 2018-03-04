% Newton-Raphson solver for strain due to stress-loading, using J2
% isotropic hardening model
function [SigN1,EpsN1,KN1,res_hist] = NewtonRaphson(SigStar,SigN,EpsN,Ce,H,G,KN,tol)
%% Step 1: Initialize k=0
SigK = SigN;                     % Sig_{n+1}^0 = Sig_{n}
EpsK = EpsN;                     % Eps_{n+1}^0 = Eps_{n}

% Initial residual
r0 = SigStar - SigK;
normr0 = voigtnorm(r0);

% Solve for initial increment in strain
[~,~,~,~,CTO] = J2(SigK,Ce,H,G,KN);
dEps = CTO\r0;
EpsK = EpsK + dEps;

% iteration
k = 1;   

res_hist = 999*ones(100,1);
res_hist(1) = 1;
res = 1000;
while res > tol
    % Step 2: Compute Sig + CTO with J2 model
    [Sig_out, K_out,~,~,CTO] = RadialReturn(SigN,Ce,H,G,dEps,KN);
   
    SigK = Sig_out;

    % Step 3: Compute residual
    rk = SigStar - SigK;
    
    % Step 4: Compute dEps
    dEK = CTO\rk;
    
    
    % Step 5: Update Eps
    EpsK = EpsK + dEK;
    dEps = dEps + dEK;

    % Check convergence
    res = voigtnorm(rk)/normr0;
    res_hist(k+1) = res;
   
    if (k > 100)
        error('Error: Too many iterations')
       
    end
    k = k + 1;
end
% Output final forms
EpsN1 = EpsK;
SigN1 = SigK;
KN1 = K_out;

end