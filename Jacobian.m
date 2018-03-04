% Calculate local jacobian for Drucker-Prager model
% mu(1) = alpha0
% mu(2) = a
% mu(3) = k
% mu(4) = beta0
function [J] = Jacobian(mu,K,G,alpha,beta,Lambda,xk)
J = zeros(3,3);

% dalpha / dlambda
if Lambda > 0
    dAlph = mu(2)*mu(3)*(mu(3)-Lambda)/sqrt(mu(3)*Lambda)/(mu(3)+Lambda)^2;
else
    dAlph = 0;
end

J(1,1) = 1;
J(1,2) = 0;
J(1,3) = beta + xk(3)*dAlph;
J(2,1) = 0;
J(2,2) = 1;
J(2,3) = 1;
J(3,1) = K*alpha;
J(3,2) = 3*G;
J(3,3) = dAlph*K*xk(1);


end