% Radial return function for 3D plasticity with isotropic hardening
function [Sig_out,nhat,dLambda,CTO] = DruckerPrager(SigStar,SigN,EpsN,dEps,Lambda,Ce,Ceinv,K,G,mu,tol)

Epsl  = EpsN + dEps;         % trial elastic strain
Evetr = sum(Epsl(1:3));               % Ev^e
eetr  = deviatoric(Epsl);             % e^e
Esetr = sqrt(2/3)*voigtnorm(eetr);    % Es^e
nhat  = eetr/voigtnorm(eetr);         

% Check Yield function
p = K*Evetr;
q = 3*G*Esetr;
[alpha,~] = Hardening(mu,Lambda);

Ftr = q + alpha*p;
Ftr

if (Ftr <= 0 ) % elastic
   dSig    = Ce(:,1)*dEps(1)
   dEps(2) = Ceinv(2,:)*dSig;
   dEps(3) = Ceinv(3,:)*dSig;
   dEps
   Sigout = SigN + dSig;
   EpsEout = EpsN+dEps; 
   Lambda_out = Lambda;
   dEpsPout = 0*voigtI2;
else        
    
end