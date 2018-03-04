clear all; close all;
%% % % % % % % % % % % % % % % % % % % % %
%      Simple 2D Finite Difference       %
%       Thermal Diffusion Model          %
%                                        %
%        Valere Lambert, 2017            %
% % % % % % % % % % % % % % % % % % % % %%

%% % % % % % % % % % % % % % % % % % % % %
%             Create 2D Grid             %
% % % % % % % % % % % % % % % % % % % % %%
% Here we consider a 2D uniform grid with 
% even spacing in both direction (dx/dy),
% however this is not a requirement as we 
% calculate the grid spacing for each grid
% therefore x and y can be made non-uniform

dx = 0.5; % m
dy = dx;      
x = (0:dx:10)';
y = (0:dy:10)';
% Number of grids in each direction
nx = length(x);
ny = length(y);

% Grid distances
dXL=zeros(nx*ny,1);
dXR=zeros(nx*ny,1);
dYL=zeros(nx*ny,1);
dYR=zeros(nx*ny,1);

for i=1:length(dXL)
    indexX = mod(i,nx);
    indexY = floor((i-1)/(nx))+1;
    if indexX==1                       %left 
        dXL(i)=1/( x(indexX+1)-x(indexX) )^2; 
        dXR(i)=1/( x(indexX+1)-x(indexX) )^2; 
    elseif indexX==0                   %right
        dXL(i)=1/( x(end)-x(end-1) )^2;
        dXR(i)=1/( x(end)-x(end-1) )^2;
    else                               %middle
        dXL(i)=1/( x(indexX)-x(indexX-1) )^2;  
        dXR(i)=1/( x(indexX+1)-x(indexX) )^2;
    end
    
    if indexY==1                       %top
        dYL(i)= 1/( y(indexY+1)-y(indexY)  )^2;
        dYR(i)= 1/( y(indexY+1)-y(indexY) )^2;
    elseif indexY==ny          %bottom
        dYL(i)= 1/( y(indexY)-y(indexY-1) )^2;  
        dYR(i)= 1/( y(indexY)-y(indexY-1) )^2; 
    else                               %middle
        dYL(i)= 1/( y(indexY)-y(indexY-1) )^2;
        dYR(i)= 1/( y(indexY+1)-y(indexY) )^2;
    end
end
dXL=reshape(dXL,nx,ny);
dXR=reshape(dXR,nx,ny);
dYL=reshape(dYL,nx,ny);
dYR=reshape(dYR,nx,ny);

%% % % % % % % % % % % % % % % % % % % % %
%        Set Model Conditions            %
% % % % % % % % % % % % % % % % % % % % %%

% Thermal diffusivity
Kappa = 2.3e-5; % m^2/s (Iron)

% Heating source
HeatZoneX = find(x > 5 & x< 8);
HeatZoneY = 1; % surface
Q = zeros(nx,ny);
Q(HeatZoneX,HeatZoneY) = Q(HeatZoneX,HeatZoneY) + 1e-4; % Degrees Celsius / s
Q = Q(:);

% Initial conditions (Degrees Celsius)
T0 = zeros(nx,ny); 
% domain = find(x >= 2 & x < 8);
% T0(domain,domain) = ...
%     T0(domain,domain) + 2;

% Boundary conditions
top    = zeros(nx,1);
bottom = zeros(nx,1);
left   = zeros(ny,1);
right  = zeros(ny,1);

% Timing
dt = 0.1;       % Time step (s)
tmax = 60*60*20; % Duration (s) - 2 hours of simulation of time
t = (0:dt:tmax)'; 

% Create vector for temperature history
Te = zeros(length(T0(:)),length(t));
Tph = zeros(length(T0(:)),length(t));
% Te = T0 at t=0
Te(:,1) = T0(:);

%% % % % % % % % % % % % % % % % % % % % %
%      Solve for Thermal Evolution       %
% % % % % % % % % % % % % % % % % % % % %%
tic
% dT/dt = Kappa * Grad^2 T + dQ/dt
for tstep = 2:length(t)
    % Temperature from previous time step
    te = Te(:,tstep-1);

    % % dT / dt
    te=reshape(te,[nx,ny]);
    Tv=[top,te,bottom];
    % % Vertical Diffusion
    Cz=dYR.*(Tv(:,3:end)-Tv(:,2:end-1)) + dYL.*(Tv(:,1:end-2)-Tv(:,2:end-1));
    % % Horizontal Diffusion
    Th=[left';te;right'];
    Cx=dXR.*(Th(3:end,:)-Th(2:end-1,:)) + dXL.*(Th(1:end-2,:)-Th(2:end-1,:));
    Tp=reshape(Cx+Cz,[nx*ny,1]);

    % Update temperature with diffusion and heating terms
    Tph(:,tstep) = Tp;
    Te(:,tstep) = te(:) + Kappa*Tp*dt + Q*dt;
end
toc
%% % % % % % % % % % % % % % % % % % % % %
%    Plotting for Thermal Profiles       %
% % % % % % % % % % % % % % % % % % % % %%
Tef = reshape(Te(:,end),nx,ny);

figure(1);clf;
subplot(1,2,1);
pcolor(x,y,T0');
set(gca,'Ydir','reverse')
title('Initial conditions')
xlabel('x (m)');
ylabel('y (m)');
c=colorbar;
ylabel(c,'Temperature (Celsius)')
axis square

subplot(1,2,2)
pcolor(x,y,Tef');
set(gca,'Ydir','reverse')
title('Thermal Profile after 2 hours')
xlabel('x (m)');
ylabel('y (m)');
c=colorbar;
ylabel(c,'Temperature (Celsius)')
axis square




