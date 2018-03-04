clear all; close all;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%    Numerical Integraion for Drucker-Prager infinitesimal       %
%                      plasticity model                          %
%                Ae 223, Valere Lambert, 2017                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Linear elastic parameters
E  = 25000;          % Young's Modulus (kPa)
nu = 0.3;            % Poisson's ratio
G  = E/(2*(1+nu));   % Shear modulus (kPa)
K  = E/(3*(1-2*nu)); % Bulk modulus (kPa)

% Plastic parameters
mu = zeros(4,1); %PIV
mu(1) = 0.7;     % alpha0
mu(2) = 0.25;    % a
mu(3) = 0.1;     % k
mu(4) = 0.7;     % beta0

% alpha = alpha0 + 2*a*sqrt(k*Lambda)/(k+Lambda)
% beta  = alpha - beta0

% Stiffness Matrix 
Ce = E/(1+nu)/(1-2*nu) * ...
     [  1-nu  nu   nu     0           0          0;
         nu  1-nu  nu     0           0          0;
         nu   nu  1-nu    0           0          0;
         0    0    0     (1-2*nu)/2   0          0;
         0    0    0      0          (1-2*nu)/2  0;
         0    0    0      0           0         (1-2*nu)/2];

Ceinv = inv(Ce);
   
% Initial stress conditions
p0 = -50; % mean normal stress (kPa)
q0 = 0;   % deviatoric stress (kPa)

%% % % % % % % % % % % % % % % % % % % %
%             Model Loading            %
% % % % % % % % % % % % % % % % % % % %%

% Applied strain increment
dEps    = zeros(6,1);
dEps(1) = -0.002;      % -0.2% vertical strail increment

% Initial strain (assumed elastic)
Eps0 = p0/(3*K);
n = ceil((-0.03 + Eps0)/dEps(1)); % number of time station

Eps   = zeros(6,n);               % Strain tensor (history)
Eps(1:3,1) = Eps0*ones(3,1);
EpsE  = Eps;                      % Elastic strain tensor

% Stress tensor and initial condition
Sigma = zeros(6,n);       % Stress tensor (history)
Sigma(:,1) = p0*voigtI2;  

% Plastic Multiplier
Lambda = zeros(1,n);

% Stress constrains
SigStar=Sigma(2:3,1);    % axisymmetric - hold Sig2 and Sig3 fixed
%SigStar = Sigma(3,1);   % plane strain - hold sig3 fixed

% Residuals for time dEps1 = -1%, -2%, -3%
res1 = 999*ones(100,1);
res2 = 999*ones(100,1);
res3 = 999*ones(100,1);

tolerance = 1e-10;
for step = 2:n
    disp(step);
    
    % Call Global Newton Raphson Solver for Strain and Stress increments
    [Sig_out,EpsE_out,Lambda_out,dEps_out,res] = GlobalNewtonRaphson(SigStar,EpsE(:,step-1),dEps,...
                        Lambda(step-1),Ce,K,G,mu,tolerance);
    Sigma(:,step) = Sig_out;
    Eps(:,step)   = Eps(:,step-1)+dEps_out; % Total Strain
    EpsE(:,step)  = EpsE_out;               % Elastic Strain
    Lambda(step)  = Lambda_out;

    % Get residuals dEps1 = -1%, -2%, -3%
    if(abs(Eps(1,step) + 0.01) < 1e-3 )
       res1 = res; 
    elseif(abs(Eps(1,step) + 0.02) < 1e-3 )
       res2 = res;
    elseif(abs(Eps(1,step) + 0.03) < 1e-3)
       res3 = res;
    end
    
end
%%

% Compute strain and stress components 
Ev = sum(Eps(1:3,:),1);  % Volumetric strain magnitude
e = zeros(size(Eps));    % Deviatoric strain  
Es = zeros(size(Ev));    % Deviatoric strain magnitude
q = zeros(size(Ev));     % Deviatoric stress
for i = 1:n
    e(:,i) = deviatoric(Eps(:,i));
    Es(i)  = sqrt(2/3)*voigtnorm(e(:,i));
    q(i)   = sqrt(3/2)*voigtnorm(deviatoric(Sigma(:,i)));
end
p = 1/3 * sum(Sigma(1:3,:),1); % mean normal stress history

% Deviatoric Strain magnitude versus deviatoric strain
figure(1);clf;
plot(Es,q,'k:','LineWidth',2);
xlabel('\epsilon_{s}')
ylabel('q')
set(gca,'FontSize',14)
grid on

% Volumetric strain magnitude versus deviatoric strain magnitude
figure(2);clf;
plot(Es,Ev,'k:','LineWidth',2);
xlabel('\epsilon_{s}')
ylabel('\epsilon_{v}')
set(gca,'FontSize',14)
grid on

% Convergence figure
figure(3);clf;
k1 = find(res1 ~=999);
k2 = find(res2 ~=999);
k3 = find(res3 ~=999);
semilogy((1:length(k1))-1,res1(k1),'k:','LineWidth',2); hold on;
semilogy((1:length(k2))-1,res2(k2),'r-.','LineWidth',2);
semilogy((1:length(k3))-1,res3(k3),'b:','LineWidth',2);
xlabel('Iteration Step');
ylabel('| r | / | r_0 |');
set(gca,'FontSize',14)
legend('\Delta \epsilon_{1} = -1%','\Delta \epsilon_{1} = -2%','\Delta \epsilon_{1} = -3%');
legend boxoff
grid on
axis square
