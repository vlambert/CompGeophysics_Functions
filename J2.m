% Checks elastic-plastic state with J2 model
function [Sig_out,K_out,nhat,dLambda,CTO] = J2(Sig_in,Ce,H,G,K)
% Determine deviatoric component
Dev_trial = deviatoric(Sig_in);

% Check Yield function
Ftr = voigtnorm(Dev_trial) - K;

if (Ftr <= 0 ) % elastic
   dLambda = 0;
   nhat = 0*Dev_trial;
   Sig_out = Sig_in;
   K_out   = K;            % no change in yielding condition
   CTO = Ce;
else           % plastic
   dLambda = Ftr/(2*G+H);
   nhat = Dev_trial/voigtnorm(Dev_trial);  
   Sig_out = Sig_in - 2*G*dLambda*nhat;
   K_out    = K + H*dLambda;
   CTO = Ce - (4*G*G/(2*G+H))*(nhat*nhat') - ...
        (4*G*G*dLambda / voigtnorm(Dev_trial)) *...
        (voigtI4 - 1/3 * (voigtI2*voigtI2') - (nhat*nhat'));
end

% Check that stress is acceptable
fcheck = voigtnorm(deviatoric(Sig_out)) - K_out;
if ((fcheck - 1.e-12) > 0)
   error('Error. \nF is positive'); 
end

end