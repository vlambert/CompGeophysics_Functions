function [Sigout,EpsE_out,Lambda_out,dEps_out,res_hist] = GlobalNewtonRaphson(SigStar,SigN,EpsEN,dEps,Lambda,Ce,Ceinv,K,G,mu,tol)

Eetr  = EpsEN + dEps;         % trial elastic strain
Evetr = sum(Eetr(1:3));               % Ev^e
eetr  = deviatoric(Eetr);             % e^e
Esetr = sqrt(2/3)*voigtnorm(eetr);    % Es^e
nhat  = eetr/voigtnorm(eetr);         

% Check Yield function
p = K*Evetr;
q = 3*G*Esetr;
[alpha,~] = Hardening(mu,Lambda);

Ftr = q + alpha*p

if (Ftr <= 0 ) % elastic
   dSig    = Ce(:,1)*dEps(1)
   dEps(2) = Ceinv(2,:)*dSig;
   dEps(3) = Ceinv(3,:)*dSig;
   EpsE_out = EpsEN+dEps;
   dEps_out = dEps;
   Sigout = SigN + dSig
   Lambda_out = Lambda;
   res_hist = 999;
else        
    % For local NR
    xN = zeros(3,1);
    xN(1) = sum(EpsEN(1:3));                        % E_v^e n
    xN(2) = sqrt(2/3)*voigtnorm(deviatoric(EpsEN)); % E_s^e n
    
    Sig0 = p*voigtI2 + sqrt(2/3)*q*nhat;
    R0 = SigStar - Sig0(2:3);
    normR0 = norm(R0);
    
    Rk = 999*ones(2,1);
    l = 1; % global iteration number
    res_hist = 999*ones(10,1);
    res_hist(1) = 1;
    res = 1000;
    while res > tol
        Eetr
        % Obtain Jacobian and strains with Local NR
        [xout,Jout] = LocalNewtonRaphson(Evetr,Esetr,xN,Lambda,K,G,mu,tol);
        
        % Output strains
        Evel  = xout(1);
        Esel  = xout(2);
        dLaml = xout(3);

        % Calculate Stress
        pl = K*Evel;
        ql = 3*G*Esel;
        Sigl = pl*voigtI2 + sqrt(2/3)*ql*nhat;

        % Calculate CTO
        Ji = inv(Jout);
        dEvedE = Ji(1,1)*voigtI2 + sqrt(2/3)*Ji(1,2)*nhat;
        dEsedE = Ji(2,1)*voigtI2 + sqrt(2/3)*Ji(2,2)*nhat;
        dndE   = sqrt(2/3)*(1/Esetr) *(voigtI4 - 1/3*(voigtI2*voigtI2')-nhat*nhat');

        CTO = K*voigtI2*dEvedE' + sqrt(2/3)*(3*G*nhat*dEsedE' + ql*dndE);
        
        % Compute residual
        Rk = SigStar - Sigl(2:3);

        % Step 4: Compute dEps
        dEl = CTO(2:3,2:3)\Rk;

        % Step 5: Update dEps 
        dEps(2:3) = dEps(2:3) + dEl;  
        
        % Update lambda and Eetr
        Lambda = Lambda+dLaml;
        [~,beta]=Hardening(mu,Lambda);
        Eetr = Eetr - dLaml*(1/3*beta*voigtI2 + sqrt(3/2)*nhat);
        
        % Compute elastic strain components
        Evetr = sum(Eetr(1:3));               % Ev^e
        eetr  = deviatoric(Eetr);             % e^e
        Esetr = sqrt(2/3)*voigtnorm(eetr);    % Es^e
        nhat  = eetr/voigtnorm(eetr);       

        % Check convergence
        res = norm(Rk)/normR0
        res_hist(l) = res;

        if (l > 10)
            error('Error: Too many iterations')

        end
        l = l + 1;
    end

    % Output final forms
    Lambda_out = Lambda + dLaml;
    Sigout = Sigl
    dEps_out  = dEps;
    EpsE_out = Eetr;
end
