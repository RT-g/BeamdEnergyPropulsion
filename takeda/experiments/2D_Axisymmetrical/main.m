%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               basic MUSCL solver for Euler system equations
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                         U_t + F(U)_x + G(U)_r + H(U) = 0,
%
% MUSCL based numerical schemes extend the idea of using a linear
% piecewise approximation to each cell by using slope limited left and
% right extrapolated states. This results in the following high
% resolution, TVD discretisation scheme.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Refs:
%   [1] Toro, E. F., "Riemann Solvers and Numerical Methods for Fluid
%   Dynamics" Springer-Verlag, Second Edition, 1999.
%   [2] Balsara, Dinshaw S. "A two-dimensional HLLC Riemann solver for
%   conservation laws: Application to Euler and magnetohydrodynamic flows."
%   Journal of Computational Physics 231.22 (2012): 7476-7503.
%   [3] Einfeldt, Bernd. "On Godunov-type methods for gas dynamics." SIAM
%   Journal on Numerical Analysis 25.2 (1988): 294-318.
%   [4] Kurganov, Alexander, and Eitan Tadmor. "Solution of two-dimensional
%   Riemann problems for gas dynamics without Riemann problem solvers."
%   Numerical Methods for Partial Differential Equations 18.5 (2002): 584-608.
%   [5] Vides, Jeaniffer, Boniface Nkonga, and Edouard Audit. "A simple
%   two-dimensional extension of the HLL Riemann solver for gas dynamics."
%   (2014).
%
% coded by Ryota Takeda, 2021.07.01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all; clc;
global Param % Values specific to gas

%% Parameters
Settings;

% Discretize spatial domain
Lx=2.5*1e-3; dx=Lx/Param.hp.nx; xc=dx/2:dx:Lx;
Lr=2.5*1e-3; dr=Lr/Param.hp.nr; rc=dr/2:dr:Lr;
[x,r] = meshgrid(xc,rc);

% Set Initial Condition
[Q0,a0] = setIC2d(x,r,Param.hp.IC);

% Run scheme
switch Param.hp.fluxMth
case ChakaravarthyOsher
    fluxMth = @ChakaravarthyOsher; % 2nd-order HLLE2d (working on it)
    O = 2;
otherwise, error('flux assamble not available');
end

% Set q-array & adjust grid for ghost cells
% 開放端の場合は精度に応じて必要となるが、壁面の場合は要らない。
[nx,nr,in,jn] = setGhostCells(Param.hp.nx,Param.hp.nr,O);
q0=zeros(nx,nr,4); q0(in,jn,1:4)=Q0;

% Boundary Conditions in ghost cells
q0 = setBC2d(q0,nx,nr,O);

% Discretize time domain
dt0 = Param.hp.CFL*min(dx,dr)/a0;

% Initialize parpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj); parpool('local',2); end



% Configure figure
if plotFig
    figure(1);
    subplot(2,2,1); [~,h1]=contourf(x,r,r0); axis('square'); xlabel('x'); ylabel('r'); title('\rho');
    subplot(2,2,2); [~,h2]=contourf(x,r,u0); axis('square'); xlabel('x'); ylabel('r'); title('u_x');
    subplot(2,2,3); [~,h3]=contourf(x,r,v0); axis('square'); xlabel('x'); ylabel('r'); title('u_r');
    subplot(2,2,4); [~,h4]=contourf(x,r,p0); axis('square'); xlabel('x'); ylabel('r'); title('p');
end

%% Solver Loop

% Load GC
q=q0; t=0; it=0; dt=dt0; a=a0; Param.LV.x_laser0=0; dq=zeros(nx,nr,4);

tic
while t < tEnd
    % caluculate S, Uionz, x at the center
    Power_laser = getPowerLaser(t); %変数をtから定数にすると定常解が出せる
    Param.LV.S_laser0 = Param.LV.R_peak * Power_laser/4/Param.LV.W_G/Param.LV.W_T * 1e-3; % GW/m^2 波頭のレーザー強度
    u_ionz0 = getSWVelocity(0);
    Param.LV.x_laser0 = Param.LV.x_laser0 + u_ionz0 * dt; %m Ionized wave front これはξ座標とも一致

    % caluculate ΔQ
    dq = fluxMth(q,dq,dt,dx,dr,nx,nr);
    q = q + dq;

    q = setBC2d(q,nx,nr,O);

	% Compute flow properties
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho;
    p=(Param.GC.gamma-1)*rho.*(E-0.5*(u.^2+v.^2)); c=sqrt(Param.GC.gamma*p./rho);

    % Update dt and time
    vn=sqrt(u.^2+v.^2); lambda1=vn+c; lambda2=vn-c;
    a = max(abs([lambda1(:);lambda2(:)]));
    dt=CFL*min(dx/a,dr/a); if t+dt>tEnd; dt=tEnd-t; end
	t=t+dt; it=it+1;

    % Plot figure
    if plotFig && rem(it,2) == 0
        set(h1,'ZData',rho(in,jn));
        set(h2,'ZData',u(in,jn));
        set(h3,'ZData',v(in,jn));
        set(h4,'ZData',p(in,jn));
        drawnow
    end
end
cputime = toc;

% Remove ghost cells
q=q(in,jn,1:4); nx=nx-2; nr=nr-2;

% compute flow properties
rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho; p=(Param.GC.gamma-1)*rho.*(E-0.5*(u.^2+v.^2));

%% Calculation of flow parameters
c = sqrt(Param.GC.gamma*p./rho);   % Speed of sound
Mx = u./c; Mr = v./c; U = sqrt(u.^2+v.^2); M = U./c;
p_ref = 101325;         % Reference air pressure (N/m^2)
rho_ref= 1.225;           % Reference air density (kg/m^3)
s = 1/(Param.GC.gamma-1)*(log(p/p_ref)+Param.GC.gamma*log(rho_ref./rho));
                        % Entropy w.r.t reference condition
ss = log(p./rho.^Param.GC.gamma);  % Dimensionless Entropy
r_x = rho.*u;             % Mass Flow rate per unit area
r_r = rho.*v;             % Mass Flow rate per unit area
e = p./((Param.GC.gamma-1)*rho);   % internal Energy

%% Final plot
figure(2); offset=0.05; n=22; % contour lines
s1=subplot(2,3,1); contour(x,r,rho,n); axis('square'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); contour(x,r,U,n); axis('square'); xlabel('x(m)'); ylabel('Velocity Magnitud (m/s)');
s3=subplot(2,3,3); contour(x,r,p,n); axis('square'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); contour(x,r,ss,n);axis('square'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); contour(x,r,M,n); axis('square'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); contour(x,r,e,n); axis('square'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'MUSCL with genuinely 2D HLL fluxes');
