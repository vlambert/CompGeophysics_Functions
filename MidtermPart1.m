clear all; close all;
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%    Solver for 3D continuum problems incorporating isotropic    %
%    linear elasticity and elasto-plastic von Mises yielding     %
%                                                                %
%                Ae 223, Valere Lambert, 2017                    %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% Linear elastic parameters
E  = 300;          % Young's Modulus (MPa)
nu = 0.45;         % Poisson's ratio
G  = E/(2*(1+nu)); % Shear modulus (MPa)

% Stiffness Matrix 
Ce = E/(1+nu)/(1-2*nu) * ...
     [  1-nu  nu   nu     0           0          0;
         nu  1-nu  nu     0           0          0;
         nu   nu  1-nu    0           0          0;
         0    0    0     (1-2*nu)/2   0          0;
         0    0    0      0          (1-2*nu)/2  0;
         0    0    0      0           0         (1-2*nu)/2];

Ceinv = inv(Ce);

% Plasticity parameters
ET = 0.1*E;                % Tangent modulus (Isotropic hardening)
ET = 0;                     % Tangent modulus (Perfect plasticity)
Kappa = 0.35;               % size of initial yield surface (MPa)
H     = 2/3 * E*ET/(E-ET) ; % Hardening parameter

%% % % % % % % % % % % % % % % % % % % %
%      Part (a) Uniaxial Loading       %
% % % % % % % % % % % % % % % % % % % %%

% Imposed Strain Tensor (Voigt notation)
Eps = zeros(6,1);
Eps(1) = 0.00005; % E11
Eps(2) = 0; % E22
Eps(3) = 0; % E33
Eps(4) = 0; % E12
Eps(5) = 0; % E13
Eps(6) = 0; % E23

Eps(4:6) = 2*Eps(4:6); % symmetric tensor

% Number of time steps
n = 1000;

Sigma = zeros(6,n);      % Stress tensor (history)
S     = zeros(6,n);      % Deviatoric stress (history)

Epshist   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;
dEpshist = zeros(6,n);
dSighist = zeros(6,n);
Sigouthist = zeros(6,n);
dEps = Eps;
for step = 2:n
    % Strain increment
    Epshist(:,step) = Epshist(:,step-1)+dEps;
    dEpshist(:,step) = dEps;
    % Calculate Stress
    [Sig_out,K_out,nhat,dLambda] = RadialReturn(Sigma(:,step-1),Ce,H,G,dEps,K(step-1));
    K(step)       = K_out;
    Sigouthist(:,step) = Sig_out;
    % Axial stress
    dSig = Ce(:,1)*dEps(1) - dLambda*nhat;  % Calculate stress increment from applied dE11
    dEps(2) = Ceinv(2,:)*dSig + dLambda*nhat(2); % Calculate resulting dE22 = dE33
    dEps(3) = Ceinv(3,:)*dSig + dLambda*nhat(3);
    dSighist(:,step) = dSig;
    Sigma(:,step) = Sigma(:,step-1) + dSig;
    S(:,step)     = deviatoric(Sigma(:,step));
    
end
%%
normS = sum(S.*S,1).^(0.5);

% Axial components of stress versus E11
figure(1);clf;
subplot(2,2,1);
plot(Epshist(1,:),Sigma(1,:),'k');
xlabel('\epsilon_{11}');
ylabel('\sigma_{11}');

axis square

subplot(2,2,2);
plot(Epshist(1,:),Sigma(2,:),'k');
xlabel('\epsilon_{11}');
ylabel('\sigma_{22}');

axis square

subplot(2,2,3);
plot(Epshist(1,:),Sigma(3,:),'k');
xlabel('\epsilon_{11}');
ylabel('\sigma_{33}');

axis square

subplot(2,2,4);
plot(Epshist(1,:),normS,'k');
xlabel('\epsilon_{11}');
ylabel('|S|');

axis square

% Visual check of axial strains
figure(2);clf;
plot(1:n,Epshist(1,:),'k-','LineWidth',2); hold on;
plot(1:n,Epshist(2,:),'r-','LineWidth',2);
plot(1:n,Epshist(3,:),'b-','LineWidth',2);
xlabel('Load step');
ylabel('Cumulative Strain');
legend('Location','northwest','\epsilon_{11}','\epsilon_{22}','\epsilon_{22}');
legend boxoff
set(gca,'FontSize',14)

Sigma(:,end)
return
%% % % % % % % % % % % % % % % % % % % %
%      Part (b) Simple Shear           %
% % % % % % % % % % % % % % % % % % % %%


% Imposed Strain Tensor (Voigt notation)
Eps = zeros(6,1);
Eps(1) = 0; % E11
Eps(2) = 0; % E22
Eps(3) = 0; % E33
Eps(4) = 0.01; % E12
Eps(5) = 0; % E13
Eps(6) = 0; % E23

Eps(4:6) = 2*Eps(4:6); % symmetric tensor

% Trials for different number of time steps
% Number of time steps
n = 10;

Sigma10 = zeros(6,n);      % Stress tensor (history)

Epshist10   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

dEps = Eps/n;
for step = 2:n
    % Strain increment
    Epshist10(:,step) = Epshist10(:,step-1)+dEps;
    
    % Calculate Stress
    [Sig_out,K_out,nhat,dLambda] = RadialReturn(Sigma10(:,step-1),Ce,H,G,dEps,K(step-1));
    Sigma10(:,step) = Sig_out;
    K(step)       = K_out;
   
end

% Number of time steps
n = 100;

Sigma100 = zeros(6,n);      % Stress tensor (history)

Epshist100   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

dEps = Eps/n;
for step = 2:n
    % Strain increment
    Epshist100(:,step) = Epshist100(:,step-1)+dEps;
    
    % Calculate Stress
    [Sig_out,K_out,nhat,dLambda] = RadialReturn(Sigma100(:,step-1),Ce,H,G,dEps,K(step-1));
    Sigma100(:,step) = Sig_out;
    K(step)       = K_out;
   
end

% Number of time steps
n = 1000;

Sigma1000 = zeros(6,n);      % Stress tensor (history)

Epshist1000   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

dEps = Eps/n;
for step = 2:n
    % Strain increment
    Epshist1000(:,step) = Epshist1000(:,step-1)+dEps;
    
    % Calculate Stress
    [Sig_out,K_out,nhat,dLambda] = RadialReturn(Sigma1000(:,step-1),Ce,H,G,dEps,K(step-1));
    Sigma1000(:,step) = Sig_out;
    K(step)       = K_out;
   
end
%%
figure(3);clf;
plot(Epshist10(4,:)/2,Sigma10(4,:),'k--','LineWidth',2);hold on
plot(Epshist100(4,:)/2,Sigma100(4,:),'r--','LineWidth',2);
plot(Epshist1000(4,:)/2,Sigma1000(4,:),'b--','LineWidth',2);
xlabel('\epsilon_{12}');
ylabel('\sigma_{12}');
title('Simple Shear Loading');
xlim([min(Epshist1000(4,:)/2) max(Epshist1000(4,:)/2)]);
ylim([ 0 max(Sigma1000(4,:))+0.5]);
legend('Location','southeast','n = 10','n = 100','n = 100');
legend boxoff
set(gca,'FontSize',14)

Sigma1000(:,end)

%% % % % % % % % % % % % % % % % % % % %
%      Part (c) 3D Loading             %
% % % % % % % % % % % % % % % % % % % %%

% Imposed Strain Tensor (Voigt notation)
Eps = zeros(6,1);
Eps(1) = 0.05;  % E11
Eps(2) = -0.03; % E22
Eps(3) = -0.03; % E33
Eps(4) = 0.01;  % E12
Eps(5) = 0.025;  % E13
Eps(6) = 0.00;  % E23

Eps(4:6) = 2*Eps(4:6); % symmetric tensor

% Number of time steps
n = 10000;

Sigma = zeros(6,n);      % Stress tensor (history)
S     = zeros(6,n);      % Deviatoric stress (history)

Epshist   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

dEps = Eps/n;
for step = 2:n
    % Strain increment
    Epshist(:,step) = Epshist(:,step-1)+dEps;
    
    % Calculate Stress
    [Sig_out,K_out,nhat,dLambda] = RadialReturn(Sigma(:,step-1),Ce,H,G,dEps,K(step-1));
    Sigma(:,step) = Sig_out;
    K(step)       = K_out;
    S(:,step)     = deviatoric(Sig_out);
   
end

normEps = (sum(Epshist.*Epshist,1)).^(0.5);

figure(4);clf;
subplot(2,3,1);
plot(normEps,Sigma(1,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{11} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
xlim([ 0 0.02]);
axis square

subplot(2,3,2);
plot(normEps,Sigma(2,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{22} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
axis square

subplot(2,3,3);
plot(normEps,Sigma(3,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{33} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
axis square

subplot(2,3,4);
plot(normEps,Sigma(4,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{12} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
axis square

subplot(2,3,5);
plot(normEps,Sigma(5,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{13} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
axis square

subplot(2,3,6);
plot(normEps,Sigma(6,:),'k');
xlabel('|\epsilon|')
ylabel('\sigma_{23} (MPa)')
set(gca,'FontSize',14)
xlim([0 max(normEps)]);
axis square

figure(5);clf;
plot(1:n,Epshist(1,:),'k-','LineWidth',2); hold on;
plot(1:n,Epshist(2,:),'r-','LineWidth',2);
plot(1:n,Epshist(3,:),'b-','LineWidth',2);
plot(1:n,Epshist(4,:)/2,'m-','LineWidth',2); 
plot(1:n,Epshist(5,:)/2,'y-','LineWidth',2);
plot(1:n,Epshist(6,:)/2,'g-','LineWidth',2);
xlabel('Load step');
ylabel('Cumulative Strain');
legend('Location','northwest','\epsilon_{11}','\epsilon_{22}','\epsilon_{33}',...
                              '\epsilon_{12}','\epsilon_{13}','\epsilon_{23}');
legend boxoff
set(gca,'FontSize',14)


Sigma(:,end)
