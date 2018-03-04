function [yp] = odefunAntiplane(~,y,ss)
% function ODEFUNANTIPLANE describes the evolution of the ordinary
% differential equation y' = f(t,y), where the state 
% vector y is 
%
%        /        s          \
%    y = |       tau         |
%        \ log(theta Vo / L) /
%
% Instead of integrating numerically the againg law
%
%    d theta / dt = 1 - V theta / L
%
% as is, we operate the following change of variable
%
%    phi = ln (theta Vo / L)
%
% Considering the following relationships
%
%    d phi = d theta / theta
%
%    d theta = d phi theta = d phi exp(phi) L / Vo
%
%    d theta / dt = d phi exp(phi) / dt L / Vo
%
%    d phi exp(phi) / dt L / Vo = 1 - V theta / L = 1 - V exp(phi) / Vo
%
% we obtain the evolution law for the new variable
% 
%    d phi / dt = ( Vo exp(-phi) - V ) / L
%
% rigidity (MPa)
G=30e3;

% shear stress
tau=y(2:ss.dgf:end);

% state variable
th=y(3:ss.dgf:end);

% norm of velocity
% V=2*ss.Vo.*exp(...
%     (tau-ss.mu0.*ss.sigma-ss.b.*ss.sigma.*th)./(ss.a.*ss.sigma));

% velocity with radiation damping for approximation of inertial terms

V=(2.*ss.Vs.*ss.a.*ss.sigma./G).*...
    Lambert_W(G*ss.Vo./(2*ss.Vs.*ss.a.*ss.sigma).*...
    exp((tau-ss.mu0.*ss.sigma-ss.sigma.*ss.b.*th)./(ss.sigma.*ss.a)));

% initialize yp
yp=zeros(size(y));

% velocity
yp(1:ss.dgf:end)=V;

% shear stress rate of change
yp(2:ss.dgf:end)=ss.K*(V-ss.Vpl);

% rate of state
yp(3:ss.dgf:end)=(ss.Vo.*exp(-th)-V)./ss.L;

end
