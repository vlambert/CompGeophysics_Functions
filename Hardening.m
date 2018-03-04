% Calculate hardening parameters for Drucker-Prager model
% mu(1) = alpha0
% mu(2) = a
% mu(3) = k
% mu(4) = beta0
function [alpha,beta] = Hardening(mu,Lambda)
alpha = mu(1) + 2*mu(2)*sqrt(mu(3)*Lambda)/(mu(3)+Lambda);
beta = alpha - mu(4);
end