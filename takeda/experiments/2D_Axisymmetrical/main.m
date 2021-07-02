% Title : 2D Axisymmetric CFD, general coordinate
% Time step: 1st
clear; close all; clc;
global Param % Values specific to gas

%% Parameters
% Hyperparameter (Computation Values)
CFL     = 0.50;     % CFL number. 1以下にする 0.5~0.99くらい;
tEnd    = 0.05;     % Final time;
nx      = 50;      % Number of cells/Elements in x;
nr      = 50;      % Number of cells/Elements in r;
IC      = 'constant';       % 19 IC cases are available;
gas     = 'air';    % selected gas, air/argon/helium;
method  = 4;	    % 1:Dim by Dim, 2:HLLE2d 1st-order, 3:HLLE2d 2nd-order;
plotFig = true;     % true:visualize evolution;
Param.hp.coodinate = 'general'; %coordinates, general/cartesian
Param.hp.fluxMth ='ChakravarthyOsher';  % HLLE1d, HLLE2d;
Param.hp.limiter = 'mm';
Param.hp.theta_1 = 1;
Param.hp.theta_2 = 0.5;
Param.hp.kappa = 1/3;
Param.hp.bcomp = 2; %圧縮パラメータ, (1,4]
Param.hp.eta = 0.1; %圧縮パラメータ, (1,4]

% Laser Values
setLaserConstants

% Discretize spatial domain
Lx=2.5*1e-3; dx=Lx/nx; xc=dx/2:dx:Lx;
Lr=2.5*1e-3; dr=Lr/nr; rc=dr/2:dr:Lr;
[x,r] = meshgrid(xc,rc);

% Set Initial Condition
[r0,u0,v0,p0] = setEulerInitialCondition2d(x,r,IC,gas);
E0 = p0./((Param.GC.gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
c0 = sqrt(Param.GC.gamma*p0./r0);                    % Speed of sound
Q0 = cat(3, r0, r0.*u0, r0.*v0, r0.*E0);    % initial state

% Set q-array & adjust grid for ghost cells
nx=nx+2; nr=nr+2; q0=zeros(nx,nr,4); q0(2:nx-1,2:nr-1,1:4)=Q0;

% Boundary Conditions in ghost cells
q0(:,1,:)=q0(:,2,:); q0(:,nr,:)=q0(:,nr-1,:);   % Natural BCs
q0(1,:,:)=q0(2,:,:); q0(nx,:,:)=q0(nx-1,:,:);   % Natural BCs

% Discretize time domain
vn = sqrt(u0.^2+v0.^2); lambda1=vn+c0; lambda2=vn-c0;
a0 = max(abs([lambda1(:);lambda2(:)]));
dt0=CFL*min(dx./a0,dr./a0);

% Initialize parpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj); parpool('local',2); end

% Run scheme
switch method
    case 1, MUSCL_EulerRes2d = @MUSCL_EulerRes2d_v0; % Do HLLE1d Dim by Dim
    case 2, MUSCL_EulerRes2d = @MUSCL_EulerRes2d_v1; % 1st-order HLLE2d (working on it)
    case 3, MUSCL_EulerRes2d = @MUSCL_EulerRes2d_v2; % 2nd-order HLLE2d (working on it)
    case 4, MUSCL_EulerRes2d = @MUSCL_HLLE; % 2nd-order HLLE2d (working on it)
    otherwise, error('flux assamble not available');
end

% Configure figure
in=2:nx-1; jn=2:nr-1; % internal indexes
if plotFig
    figure(1);
    subplot(2,2,1); [~,h1]=contourf(x,r,r0); axis('square'); xlabel('x'); ylabel('r'); title('\rho');
    subplot(2,2,2); [~,h2]=contourf(x,r,u0); axis('square'); xlabel('x'); ylabel('r'); title('u_x');
    subplot(2,2,3); [~,h3]=contourf(x,r,v0); axis('square'); xlabel('x'); ylabel('r'); title('u_r');
    subplot(2,2,4); [~,h4]=contourf(x,r,p0); axis('square'); xlabel('x'); ylabel('r'); title('p');
end

%% Solver Loop

% Load GC
q=q0; t=0; it=0; dt=dt0; a=a0; Param.LV.x_laser0=0; dq_ex=zeros(nx,nr,4);

tic
while t < tEnd

    % RK2 1st step
    Power_laser = getPowerLaser(t);
    Param.LV.S_laser0 = Param.LV.R_peak * Power_laser/4/Param.LV.W_G/Param.LV.W_T * 1e-3; % GW/m^2 波頭のレーザー強度
    u_ionz0 = getSWVelocity(0);
    Param.LV.x_laser0 = Param.LV.x_laser0 + u_ionz0 * dt; %m Ionized wave front これはξ座標とも一致

    dq = MUSCL_EulerRes2d(q,dq_ex,dt,dx,dr,nx,nr);
    dq_ex = dq;
    q = q + dq;

    q(:,1,:)=q(:,2,:); q(:,nr,:)=q(:,nr-1,:);   % Natural BCs
    q(1,:,:)=q(2,:,:); q(nx,:,:)=q(nx-1,:,:);   % Natural BCs


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
