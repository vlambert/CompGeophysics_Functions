% Global Newton-Raphson solver for loading in Drucker-Prager model
function [Sigout,EpsE_out,Lambda_out,dEps_out,res_hist] = GlobalNewtonRaphson(SigStar,EpsEN,dEps,Lambda,Ce,K,G,mu,tol)
% Initial elastic strain for local NR search
xN = zeros(3,1);
xN(1) = sum(EpsEN(1:3));                        % E_v^e n
xN(2) = sqrt(2/3)*voigtnorm(deviatoric(EpsEN)); % E_s^e n

% Assign initial trial elastic strain and compute components
Eetr  = EpsEN + dEps;                 % trial elastic strain
Evetr = sum(Eetr(1:3));               % Ev^e
eetr  = deviatoric(Eetr);             % e^e
Esetr = sqrt(2/3)*voigtnorm(eetr);    % Es^e
nhat  = eetr/voigtnorm(eetr);         

% Initial stress residual
p = K*Evetr;
q = 3*G*Esetr;
Sig0 = p*voigtI2 + sqrt(2/3)*q*nhat;
R0 = SigStar - Sig0(2:3);
normR0 = norm(R0);

% Begin iteration for strain increment
l = 1; % global iteration number
res_hist = 999*ones(10,1);
res = 1000;
res_hist(1) = 1;
while res > tol
    % Check yield function
    p = K*Evetr;
    q = 3*G*Esetr;
    [alpha,~]=Hardening(mu,Lambda);
    Ftr = q + alpha*p;

    % If elastic then CTO and stress are simple
    if (Ftr <= 0 ) % elastic
        CTO = Ce;
        dLaml = 0;
        Sigl  = p*voigtI2 + sqrt(2/3)*q*nhat;
    % If plastic then compute plastic strain component with local NR search
    else        
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
    end
    % Compute stress residual
    %Rl = SigStar - Sigl(2:3); % axisymmetric 
    Rl = SigStar - Sigl(3);   % plane strain

    % Step 4: Compute dEps
    dEl = CTO(2:3,2:3)\Rl; % axisymmetric solve for dEps2 and dEps3
    %dEl = CTO(3,3)\Rl; % plane strain solve for dEps3

    %  Update total strain increment
    dEps(2:3) = dEps(2:3) + dEl; % axisymmetric
    %dEps(3) = dEps(3) + dEl; % plane strain
    
    % Update Eetr
    Eetr = EpsEN + dEps;

    % Compute elastic strain components
    Evetr = sum(Eetr(1:3));               % Ev^e
    eetr  = deviatoric(Eetr);             % e^e
    Esetr = sqrt(2/3)*voigtnorm(eetr);    % Es^e
    nhat  = eetr/voigtnorm(eetr);       

    % Check convergence
    res = norm(Rl)/normR0;
    res_hist(l+1) = res;

    if (l > 10)
        error('Error: Too many iterations')
    end
    l = l + 1;
end
    % Output final forms
    Lambda_out = Lambda + dLaml;
    [~,beta] = Hardening(mu,Lambda_out);
    Sigout = Sigl;
    dEps_out  = dEps;
    EpsE_out = Eetr - dLaml*(1/3*beta*voigtI2 + sqrt(3/2)*nhat);
end