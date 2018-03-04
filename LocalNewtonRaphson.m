% Local Newton-Raphson solver for Drucker-Prager material model
function [xout,Jout] = LocalNewtonRaphson(Evetr,Esetr,xN,Lambda,K,G,mu,tol)
% xk(1) = Eps_{v}^{e}
% xk(2) = Eps_{s}^{e}
% xk(3) = dLam

%% Step 1: Initialize k=0
xk = xN;
Lambdak = Lambda;

% Initial residual
rk = 999*ones(3,1);
[alphak,betak] = Hardening(mu,Lambdak);

rk(1) = xk(1) - Evetr + xk(3)*betak;
rk(2) = xk(2) - Esetr + xk(3);
rk(3) = 3*G*xk(2) + alphak*K*xk(1);
normr0 = norm(rk);

% Solve for initial increment in strain
Jk = Jacobian(mu,K,G,alphak,betak,Lambdak,xk);
dxk = -inv(Jk)*rk;
xk = xk + dxk;

% iteration
k = 1;   
res = 1000;
while res > tol
    
    % Update Lambda and compute residual
    Lambdak = Lambda + xk(3);
    [alphak,betak] = Hardening(mu,Lambdak);
    rk(1) = xk(1) - Evetr + xk(3)*betak;
    rk(2) = xk(2) - Esetr + xk(3);
    rk(3) = 3*G*xk(2) + alphak*K*xk(1);
    
    % Check convergence
    res = norm(rk)/normr0;
   
    % Compute Jacobian
    Jk = Jacobian(mu,K,G,alphak,betak,Lambdak,xk);
    
    % Compute increment for x
    dxk = -inv(Jk)*rk;
    
    % Update x
    xk = xk + dxk;
    
    if (k > 10)
        error('Error: Too many iterations')
    end
    k = k + 1;
end
% Output final forms
xout = xk;
Jout = Jk;
end