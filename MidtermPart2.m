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
%ET = 0;                     % Tangent modulus (Perfect plasticity)
Kappa = 0.35;               % size of initial yield surface (MPa)
H     = 2/3 * E*ET/(E-ET) ; % Hardening parameter


%% % % % % % % % % % % % % % % % % % % %
%     Part (c) Uniaxial Loading        %
% % % % % % % % % % % % % % % % % % % %%

% Time
dt = 0.05; % time step (arbitrary units)
t = (0:dt:1)';
n = length(t);

Eps   = zeros(6,n);      % Strain tensor (history)
Sigma = zeros(6,n);      % Stress tensor (history)
S     = zeros(6,n);      % Deviatoric stress (history)

K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

% Residuals for time 0.2 and 0.9
res2 = zeros(100,1);
res9 = zeros(100,1);

SigAxial = [1; 0; 0; 0; 0; 0]; % Sig11 uniaxial loading 

tolerance = 1e-8;
for step = 2:n
    disp(step);
    % Applied load at time step
    SigStar = 0.5*t(step) * SigAxial; % MPa
    
    % Solve for Stress and Strain
    [Sig_out,Eps_out,K_out,res] = NewtonRaphson(SigStar,Sigma(:,step-1),Eps(:,step-1),...
                        Ce,H,G,K(step-1),tolerance);
    Sigma(:,step) = Sig_out;
    Eps(:,step)   = Eps_out;
    K(step) = K_out;
    
    % Get residuals for t = 0.2 and t = 0.9
    if(t(step) == 0.2)
       res2 = res; 
    elseif(t(step) == 0.9)
       res9 = res;
    end
    
end
%%
% Axial stress versus strain
figure(1);clf;
plot(Eps(1,:),Sigma(1,:),'k','LineWidth',2);
xlabel('\epsilon_{11}')
ylabel('\sigma_{11} (MPa)')
set(gca,'FontSize',14)
xlim([ 0 max(Eps(1,:))]);
axis square


figure(2);clf;
k2 = find(res2 ~=999);
k9 = find(res9 ~=999);
m2 = (res2(max(k2)) - res2(max(k2)-1));
m9 = (res9(max(k9)) - res9(max(k9)-1));
semilogy((1:length(k2))-1,res2(k2),'k-','LineWidth',2); hold on;
semilogy((1:length(k9))-1,res9(k9),'r-','LineWidth',2);
xlabel('Iteration k');
ylabel('| r | / | r_0 |');
set(gca,'FontSize',14)
xlim([0 max(length(k9),length(k2))-1])
ylim([min(res9(k9)) 1])
legend(sprintf('t = 0.2, m = %0.3g',m2),sprintf('t = 0.9, m = %0.3g',m9));
legend boxoff
axis square



%% % % % % % % % % % % % % % % % % % % %
%  Part (d) Asymmetric Compression     %
% % % % % % % % % % % % % % % % % % % %%

% Time
dt = 0.05; % time step (arbitrary units)
t = (0:dt:1)';
n = length(t);

Eps   = zeros(6,n);      % Strain tensor (history)
Sigma = zeros(6,n);      % Stress tensor (history)
P     = zeros(1,n);      % Pressure (history)
S     = zeros(1,n);      % Norm Deviatoric stress (history)

Epshist   = zeros(6,n);   % Strain (history)
K = zeros(1,n);           % Yield stress  (history)
K(1) = Kappa;

% Residuals time 0.2 and 0.9
res2d = zeros(100,1);
res9d = zeros(100,1);

tolerance = 1e-8;
for step = 2:n
    disp(step);
    % Applied load at time step
    if(t(step) <= 0.4)
        SigStar = [-t(step); -t(step); -t(step); 0; 0; 0]; % MPa
    else
        SigStar = [-t(step); -0.4; -0.4; 0; 0; 0]; % MPa
    end
    
    % Solve for Stress and Strain
    [Sig_out,Eps_out,K_out,res] = NewtonRaphson(SigStar,Sigma(:,step-1),Eps(:,step-1),...
                        Ce,H,G,K(step-1),tolerance);
    Sigma(:,step) = Sig_out;
    Eps(:,step)   = Eps_out;
    K(step) = K_out;
    Dev = deviatoric(Sig_out);
    S(step) = voigtnorm(Dev);
    P(step) = 1/3 * sum(Sig_out(1:3));
    
    % Get residuals for t = 0.2 and t = 0.9
    if(t(step) == 0.2)
       res2d = res; 
    elseif(t(step) == 0.9)
       res9d = res;
    end
    
end
%%
% Axial stress versus strain
figure(3);clf;
subplot(1,2,1);
plot(P,sqrt(3/2)*S,'k','LineWidth',2);
ylabel('Q (MPa)')
xlabel('P (MPa)')
set(gca,'FontSize',14)
xlim([min(P) max(P)])
axis square
subplot(1,2,2)
plot(Eps(1,:),sqrt(3/2)*S,'k','LineWidth',2);
xlabel('\epsilon_{11}')
ylabel('Q (MPa)')
xlim([min(Eps(1,:)) max(Eps(1,:))])
set(gca,'FontSize',14)
axis square

%%
eps = 1e-18;
figure(4);clf;
k2d = find(res2d ~=999);
k9d = find(res9d ~=999);
m2d = (res2d(max(k2d)) - res2d(max(k2d)-1));
m9d = (res9d(max(k9d)) - res9d(max(k9d)-1));
semilogy((1:length(k2d))-1,max(res2d(k2d),eps),'k-','LineWidth',2); hold on;
semilogy((1:length(k9d))-1,res9d(k9d),'r-','LineWidth',2);
xlabel('Iteration k');
ylabel('| r | / | r_0 |');
set(gca,'FontSize',14)
xlim([0 max(length(k9d),length(k2d))-1])
ylim([min(res9d(k9d)) 1])
legend(sprintf('t = 0.2, m = %0.3g',m2d),sprintf('t = 0.9, m = %0.3g',m9d));
legend boxoff
axis square

