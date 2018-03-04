% Radial return function for 3D plasticity with isotropic hardening
function [Sig_out,K_out,nhat,dLambda,CTO] = RadialReturn(Sig_in,Ce,H,G,dEps,K)
% Calculate trial stress    
Sig_trial = Sig_in + Ce*dEps;

% Determine deviatoric component
Dev_trial = deviatoric(Sig_trial);

% Check Yield function
Ftr = voigtnorm(Dev_trial) - K;

if (Ftr <= 0 ) % elastic
   dLambda = 0;
   nhat = 0*Dev_trial;
   Sig_out = Sig_trial;
   K_out   = K;            % no change in yielding condition
   CTO = Ce;
else           % plastic
   dLambda = Ftr/(2*G+H);
   nhat = Dev_trial/voigtnorm(Dev_trial);  
   Sig_out = Sig_trial - 2*G*dLambda*nhat;
   K_out    = K + H*dLambda;
   CTO = Ce - (4*G*G/(2*G+H))*(nhat*nhat');% - ...
        %(4*G*G*dLambda / voigtnorm(Dev_trial)) *...
        %(voigtI4 - 1/3 * (voigtI2*voigtI2') - (nhat*nhat'));
end

% Check that stress is acceptable
fcheck = voigtnorm(deviatoric(Sig_out)) - K_out;
if ((fcheck - 1.e-12) > 0)
   error('Error. \nF is positive'); 
end

end