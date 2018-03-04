clear all; close all;
%% % % % % % % % % % % % % % % % % % % % %
%      Simple 2D Finite Difference       %
%       Thermal Diffusion Model          %
%                                        %
%        Valere Lambert, 2017            %
% % % % % % % % % % % % % % % % % % % % %%

%% % % % % % % % % % % % % % % % % % % % %
%      Solve for Thermal Evolution       %
% % % % % % % % % % % % % % % % % % % % %%

% Create 2D grid (Assume that grid spacing is even in both directions)
dx = 0.5; % m
dy = dx;      
x = (0:dx:10)';
y = (0:dy:10)';
% Number of grids in each direction
nx = length(x);
ny = length(y);

% Thermal diffusivity
Kappa = 2.3e-5; % m^2/s (Iron)

% Initial conditions
T0 = zeros(nx,ny);
domain = find(x >= 2 & x < 8);
T0(domain,domain) = ...
    T0(domain,domain) + 2;

% Timing
dt = 0.1;       % Time step (s)
tmax = 60*60*2; % Duration (s) - 2 hours of simulation of time
t = (0:dt:tmax)'; 

% Create vector for temperature history
Te = zeros(length(T0(:)),length(t));
Tph = zeros(length(T0(:)),length(t));
% Te = T0 at t=0
Te(:,1) = T0(:);

%% % % % % % % % % % % % % % % % % % % % %
%      Solve for Thermal Evolution       %
% % % % % % % % % % % % % % % % % % % % %%
% dT/dt = Kappa * Grad^2 T
for tstep = 2:length(t)
    % Initialize temperature rate of change dT/dt
    Tp = zeros(size(T0(:)));
    % Temperature from previous time step
    te = Te(:,tstep-1);
    % Boundary condition (0 <= Co <= 0.5)
    % Co = 0    No heat loss from system
    % Co = 0.5  Absorbing boundary
    Co=0.;
    Bound = Co/((Co-1));
    % Loop over cells to solve for thermal gradient
    for i=1:length(te)
        indexX = mod(i,nx);
        indexY = floor((i-1)/nx)+1;
        if indexX==1        % Left boundary
            Cx = ( te(i+1) - te(i)   ) / ( dx )^2+...
                 (Bound * te(i) ) / ( dx )^2;   
        elseif indexX==0    % Right boundary
            Cx = -( te(i)  - te(i-1) )/ ( dx )^2+...
                  (Bound * te(i)) / ( dx)^2;
        else                % Central grids
            Cx = -( te(i)  - te(i-1) )/ ( dx )^2+...
                  ( te(i+1)- te(i)   )/ ( dx )^2;
        end

        if indexY==1        % Top boundary 
            Cy = ( te(i+nx)- te(i)      ) / ( dy )^2+...
                 (Bound * te(i) ) / ( dy )^2;   
        elseif indexY==ny   % Bottom boundary
            Cy = -( te(i)   - te(i-nx)   ) / ( dy )^2+...
                  (Bound * te(i)) / ( dy)^2;
        else                % Central grids
            Cy = ( te(i+nx)- te(i)      ) / ( dy   )^2-...
                 ( te(i)   - te(i-nx)   ) / ( dy )^2;
        end
        Tp(i)=Cx+Cy; 
    end
    % Update temperature
    Tph(:,tstep) = Tp;
    Te(:,tstep) = te + Kappa*Tp*dt;
end
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


