
%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
% Evaluates the slip history on a fault in antiplane   %
% strain under the rate- and state-dependent friction  %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%


%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%            P H Y S I C A L   M O D E L               %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% rigidity (MPa)
G=30e3;

% stress-interaction function
s12h=@(x2,x3,y2,y3,W) G*( ...
    -(x3-y3)./((x2-y2).^2+(x3-y3).^2)+(x3+y3)./((x2-y2).^2+(x3+y3).^2) ...
    +(x3-y3-W)./((x2-y2).^2+(x3-y3-W).^2)-(x3+y3+W)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

s13h=@(x2,x3,y2,y3,W) G*( ...
    (x2-y2)./((x2-y2).^2+(x3-y3).^2)-(x2-y2)./((x2-y2).^2+(x3+y3).^2) ...
    -(x2-y2)./((x2-y2).^2+(x3-y3-W).^2)+(x2-y2)./((x2-y2).^2+(x3+y3+W).^2) ...
    )/2/pi;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                        M E S H                       %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

Locking_Depth=35e3; %m
M=150;
dz=Locking_Depth/M;
% top of slip patch
y3=(0:M-1)'*dz;
% width (down dip) of slip patch
W=ones(M,1)*dz;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%             S T R E S S   K E R N E L S              %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%
ss.K=zeros(M,M);
for k=1:M
    % we evaluate the stress at the center 
    % of the slip patches
    ss.K(:,k)=s12h(0,y3+dz/2,0,y3(k),W(k));
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         F R I C T I O N   P A R A M E T E R S        %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

% default friction properties (velocity-weakening friction)
% static friction coefficient
ss.mu0=0.6*ones(size(y3));
% frictional parameters
ss.a=1e-2*ones(size(y3));
ss.b=0.006*ones(size(y3));
% normal stress (MPa)
ss.sigma=100.0*ones(size(y3));
% characteristic weakening distance (m)
ss.L=0.01*ones(size(y3));
% plate velocity (m/s)
ss.Vpl=1e-9*ones(size(y3));
% reference slip rate (m/s)
ss.Vo=1e-6*ones(size(y3));
% shear speed (m/s)
ss.Vs=3e3*ones(size(y3)); 

% Top 2km and bottom 5km are velocity strengthening
top=floor(2e3/(Locking_Depth/M));
bottom=ceil((Locking_Depth-5e3)/(Locking_Depth/M));
ss.b(top:bottom)=0.014*ones(length(top:bottom),1);

hstar=G.*ss.L./(ss.b-ss.a)./ss.sigma;
Ti = 5*(ss.b-ss.a).*ss.sigma.*0.5.*(y3(bottom)-y3(top))./(G*ss.Vpl);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%         N U M E R I C A L   S O L U T I O N          %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

ss.dgf=3;
Y0=zeros(M*ss.dgf,1);
Y0(1:ss.dgf:end)=zeros(M,1);
Y0(2:ss.dgf:end)=(ss.mu0+(ss.a-ss.b).*log(ss.Vpl./ss.Vo)).*ss.sigma+G*ss.Vpl./(2*ss.Vs); 
Y0(3:ss.dgf:end)=log(ss.Vo./ss.Vpl); 

% initialize the function handle with 
% set constitutive parameters
yp=@(t,y) odefunAntiplane(t,y,ss);

% solve the system
options=odeset('Refine',1,'RelTol',3e-14,'InitialStep',1e-7);
[t,Y]=ode45(yp,[0 3e10],Y0,options);  % time - 0 to 3e10 seconds 

% compute the instantaneous derivative
Yp=zeros(size(Y));
for i=1:length(t)
    Yp(i,:)=yp(t(i),Y(i,:)');
end

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                      %
%                    F I G U R E S                     %
%                                                      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %%

figure(1);clf;set(gcf,'name','Time evolution')

subplot(2,1,1);cla;
pcolor(t/3.15e7,y3/1e3,log10(Yp(:,1:ss.dgf:end)')), shading flat
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
title(h,'Velocity (m/s)')
xlabel('Time (yr)')
ylabel('Depth (km)');

subplot(2,1,2);cla;
plot(t/3.15e7,log10(Yp(:,(M/2-1)*ss.dgf+1)'))
xlabel('Time (yr)')
ylabel('Velocity (m/s) log10')
title('Time series at the fault center')


figure(2);clf;set(gcf,'name','Evolution with time steps')

subplot(2,1,1);cla;
pcolor((1:length(t)),y3/1e3,log10(Yp(:,1:ss.dgf:end)')), shading flat
set(gca,'YDir','reverse');
h=colorbar('Location','NorthOutside');
title(h,'Velocity (m/s)')
xlabel('Time Step')
ylabel('Depth (km)');

subplot(2,1,2);cla;
plot((1:length(t)),log10(Yp(:,(M/2-1)*ss.dgf+1)'))
xlabel('Step')
ylabel('Velocity (m/s) log10')
title('Evolution at the fault center')



