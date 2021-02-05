clear
close all

% Computation Values
tic

nx = 250; % Number of position segments
ny = 250; 
nt = 200; % Number of time segments
ng = 10; % Graph
dx = 0.01e-2; % Segment length
dy = 0.01e-2;
CFL = 0.2; % Courant–Friedrichs–Lewy condition



% Physical Values

f = 170e9;
mu_0 = 1.2566370614e-6; 
epsilon_0 = 8.854187817e-12; 
c = 1./sqrt(mu_0*epsilon_0);
lambda_0 = c/f;
u_ionz0 = 0.004187e3*1254^1.16*1.2; %Helium
P_0 = 1.013e5*1;
R_air = 2077.1;% 287.1 / 2077.1
T_0 = 293;
rho_0 = P_0/R_air/T_0;
u_0 = 0;
v_0 = 0;
gamma = 1.66;
sigma = 0.000864;
I = 0;
hr = 5; % Heating Region 

% Initial Values

P = zeros(nx,ny);
u = zeros(nx,ny);
v = zeros(nx,ny);
rho = zeros(nx,ny);
H = zeros(nx,ny);

P(:,:) = P_0;
u(:,:) = u_0;
v(:,:) = v_0;
rho(:,:) = P_0/R_air/T_0;


Ttemp_x = zeros(1,nx);
Ttemp_y = zeros(1,ny);


% Vector of states
U1 = ones(nx,ny).*rho; % Mass
U2 = ones(nx,ny).*rho.*u; % Momentum x
U3 = ones(nx,ny).*rho.*v; %Momentum y
U4 = ones(nx,ny).*P/(gamma-1)+1/2.*rho.*(u.^2+v.^2); % Energy


% New calculated vector of states
U1_cal = zeros(nx,ny);
U2_cal = zeros(nx,ny);
U3_cal = zeros(nx,ny);
U4_cal = zeros(nx,ny);

% Vector of flux
F1 = zeros(nx,ny);
F2 = zeros(nx,ny);
F3 = zeros(nx,ny);
F4 = zeros(nx,ny);

G1 = zeros(nx,ny);
G2 = zeros(nx,ny);
G3 = zeros(nx,ny);
G4 = zeros(nx,ny);

% Vector of flux between cells
F1_half = zeros(nx,ny);
F2_half = zeros(nx,ny);
F3_half = zeros(nx,ny);
F4_half = zeros(nx,ny);

G1_half = zeros(nx,ny);
G2_half = zeros(nx,ny);
G3_half = zeros(nx,ny);
G4_half = zeros(nx,ny);

% Calculation

t = 0;
t_laser = 0;
x_laser = 0;
umax_x = u_ionz0;
umax_y = u_ionz0;




for n1 = 1:nt
    dt = CFL*(dx/umax_x + dy/umax_y);
    t = t+dt;
    t_laser = t_laser+dt;
    
    Power_laser = 1254*exp(-t_laser*1e6/0.619);%524.2
    u_ionz = 0.004187e3*Power_laser^1.16; % Helium
    x_laser = x_laser+u_ionz*dt;
    y_laser = dy*ny/2;
    
    for n2 = 2:nx-1
        x = n2*dx;
        
        for n3 = 2:ny-1
            y = n3*dy;
            % Input Power Definition
            if x-hr*dx < x_laser && x > x_laser
                ita = 0.1;
                w = Power_laser*1e9/hr/dx*exp(-1/2*((y - y_laser)/sigma)^2);
            else
                w = 0;
            end
            

        
            % Euler Calculation
            U1_cal(n2,n3) = U1(n2,n3) - dt/dx*(F1_half(n2,n3)-F1_half(n2-1,n3)) - dt/dy*(G1_half(n2,n3)-G1_half(n2,n3-1));
            U2_cal(n2,n3) = U2(n2,n3) - dt/dx*(F2_half(n2,n3)-F2_half(n2-1,n3)) - dt/dy*(G2_half(n2,n3)-G2_half(n2,n3-1));
            U3_cal(n2,n3) = U3(n2,n3) - dt/dx*(F3_half(n2,n3)-F3_half(n2-1,n3)) - dt/dy*(G3_half(n2,n3)-G3_half(n2,n3-1));
            U4_cal(n2,n3) = U4(n2,n3) - dt/dx*(F4_half(n2,n3)-F4_half(n2-1,n3)) - dt/dy*(G4_half(n2,n3)-G4_half(n2,n3-1)) + w*dt;
        
        end
    end
    
    % Boundary Conditions

    
    % Left
    for i = 2:nx-1
        Ttemp_x(i) = (gamma-1)*(U4_cal(i,2)-(1/2)*(U2_cal(i,2).^2+U3_cal(i,2).^2)./U1_cal(i,2))./U1_cal(i,2)/R_air; % Temperature to set the condition
        U1_cal(i,1) = P_0/R_air./Ttemp_x(i);
        U2_cal(i,1) = U2_cal(i,2)./U1_cal(i,2).*U1_cal(i,1);
        U3_cal(i,1) = U3_cal(i,2)./U1_cal(i,2).*U1_cal(i,1);
        U4_cal(i,1) = U4_cal(i,2)./U1_cal(i,2).*U1_cal(i,1);
    end
    
    
    
    % Right
    for i = 2:nx-1
        Ttemp_x(i) = (gamma-1)*(U4_cal(i,ny-1)-(1/2)*(U2_cal(i,ny-1).^2+U3_cal(i,ny-1).^2)./U1_cal(i,ny-1))./U1_cal(i,ny-1)/R_air; % Temperature to set the condition
        U1_cal(i,ny) = P_0/R_air./Ttemp_x(i);
        U2_cal(i,ny) = U2_cal(i,ny-1)./U1_cal(i,ny-1).*U1_cal(i,ny);
        U3_cal(i,ny) = U3_cal(i,ny-1)./U1_cal(i,ny-1).*U1_cal(i,ny);
        U4_cal(i,ny) = U4_cal(i,ny-1)./U1_cal(i,ny-1).*U1_cal(i,ny);
    end
    
    % Top
    
    for i = 2:ny-1
        Ttemp_y(i) = (gamma-1)*(U4_cal(nx-1,i)-(1/2)*(U2_cal(nx-1,i).^2+U3_cal(nx-1,i).^2)./U1_cal(nx-1,i))./U1_cal(nx-1,i)/R_air; % Temperature to set the condition
        U1_cal(nx,i) = P_0/R_air./Ttemp_y(i);
        U2_cal(nx,i) = U2_cal(nx-1,i)./U1_cal(nx-1,i).*U1_cal(nx,i);
        U3_cal(nx,i) = U3_cal(nx-1,i)./U1_cal(nx-1,i).*U1_cal(nx,i);
        U4_cal(nx,i) = U4_cal(nx-1,i)./U1_cal(nx-1,i).*U1_cal(nx,i);
    end
    
    
    % Plate
    U1_cal(1,:) = U1_cal(2,:);
    U2_cal(1,:) = 0;
    U3_cal(1,:) = 0;
    U4_cal(1,:) = U4_cal(2,:);
    
    % Corners
    U1_cal(nx,1) = (U1_cal(nx,2) + U1_cal(nx-1,1))/2;
    U2_cal(nx,1) = (U2_cal(nx,2) + U2_cal(nx-1,1))/2;
    U3_cal(nx,1) = (U3_cal(nx,2) + U3_cal(nx-1,1))/2;
    U4_cal(nx,1) = (U4_cal(nx,2) + U4_cal(nx-1,1))/2;
    
    U1_cal(nx,ny) = (U1_cal(nx-1,ny) + U1_cal(nx,ny-1))/2;
    U2_cal(nx,ny) = (U2_cal(nx-1,ny) + U2_cal(nx,ny-1))/2;
    U3_cal(nx,ny) = (U3_cal(nx-1,ny) + U3_cal(nx,ny-1))/2;
    U4_cal(nx,ny) = (U4_cal(nx-1,ny) + U4_cal(nx,ny-1))/2;
    
    % Update the matrix
    U1 = U1_cal;
    U2 = U2_cal;
    U3 = U3_cal;
    U4 = U4_cal;

    rho = U1;
    u = U2./rho;
    v = U3./rho;
    P = (gamma-1)*(U4-(1/2)*rho.*(u.^2+v.^2));
    H = (U4+P)./rho;
    
    
    F1 = U2;
    F2 = rho.*u.^2+P;
    F3 = rho.*u.*v;
    F4 = (U4+P).*u;
    
    G1 = U3;
    G2 = rho.*u.*v;
    G3 = rho.*v.^2+P;
    G4 = (U4+P).*v;
    
    umax_x = u_ionz;
    umax_y = u_ionz;
    
    for n2 = 1:nx-1
        
        for n3 = 1:ny-1
        
            % Calculation of average values between cells
            
            % Direction x
            
            % Expressing rho, u, v and H in the basis of w (Roe's paper)
            rho_av_x = sqrt(rho(n2,n3)*rho(n2+1,n3));
            u_av_x = (sqrt(rho(n2,n3))*u(n2,n3)+sqrt(rho(n2+1,n3))*u(n2+1,n3))/(sqrt(rho(n2,n3))+sqrt(rho(n2+1,n3)));
            v_av_x = (sqrt(rho(n2,n3))*v(n2,n3)+sqrt(rho(n2+1,n3))*v(n2+1,n3))/(sqrt(rho(n2,n3))+sqrt(rho(n2+1,n3)));
            H_av_x = (sqrt(rho(n2,n3))*H(n2,n3)+sqrt(rho(n2+1,n3))*H(n2+1,n3))/(sqrt(rho(n2,n3))+sqrt(rho(n2+1,n3)));
            
            a_av_x = sqrt((gamma-1)*(H_av_x-1/2*(u_av_x^2 + v_av_x^2))); %Square root of the derivative of P by rho, used in the eigenvalues of the Jacobian matrix
            
            % Eigenvalues of the Jacobian matrix
            lambda_x(1) = u_av_x-a_av_x;
            lambda_x(2) = u_av_x;
            lambda_x(3) = u_av_x+a_av_x;
            lambda_x(4) = u_av_x;
            
            if (umax_x<max(abs(lambda_x)))
                umax_x = max(abs(lambda_x));
            end
            
            %R Matrix constructed from eigenvectors in w-space
            
            R_x = zeros(4,4);
            R_x(1,:) = 1;
            R_x(1,4) = 0;
            R_x(2,:) = lambda_x;
            R_x(2,4) = 0;
            R_x(3,:) = v_av_x;
            R_x(3,4) = 1;
            R_x(4,1) = H_av_x-u_av_x*a_av_x;
            R_x(4,2) = 1/2*(u_av_x^2 + v_av_x^2);
            R_x(4,3) = H_av_x+u_av_x*a_av_x;
            R_x(4,4) = v_av_x;
            
            
            % Solving the linear problem
            
            A_x = zeros(4,4);
            for i =1:4
                A_x(i,i) = abs(lambda_x(i));
            end
            
            RAR_inv_x = R_x*A_x/R_x;
            
            S_x = zeros(4,1);
            S_x(1) = RAR_inv_x(1,1)*(U1(n2+1,n3)-U1(n2,n3)) + RAR_inv_x(1,2)*(U2(n2+1,n3)-U2(n2,n3)) + RAR_inv_x(1,3)*(U3(n2+1,n3)-U3(n2,n3)) + RAR_inv_x(1,4)*(U4(n2+1,n3)-U4(n2,n3));
            S_x(2) = RAR_inv_x(2,1)*(U1(n2+1,n3)-U1(n2,n3)) + RAR_inv_x(2,2)*(U2(n2+1,n3)-U2(n2,n3)) + RAR_inv_x(2,3)*(U3(n2+1,n3)-U3(n2,n3)) + RAR_inv_x(2,4)*(U4(n2+1,n3)-U4(n2,n3));
            S_x(3) = RAR_inv_x(3,1)*(U1(n2+1,n3)-U1(n2,n3)) + RAR_inv_x(3,2)*(U2(n2+1,n3)-U2(n2,n3)) + RAR_inv_x(3,3)*(U3(n2+1,n3)-U3(n2,n3)) + RAR_inv_x(3,4)*(U4(n2+1,n3)-U4(n2,n3));
            S_x(4) = RAR_inv_x(4,1)*(U1(n2+1,n3)-U1(n2,n3)) + RAR_inv_x(4,2)*(U2(n2+1,n3)-U2(n2,n3)) + RAR_inv_x(4,3)*(U3(n2+1,n3)-U3(n2,n3)) + RAR_inv_x(4,4)*(U4(n2+1,n3)-U4(n2,n3));
            
            F1_half(n2,n3) = 1/2*(F1(n2+1,n3)+F1(n2,n3)-S_x(1));
            F2_half(n2,n3) = 1/2*(F2(n2+1,n3)+F2(n2,n3)-S_x(2));
            F3_half(n2,n3) = 1/2*(F3(n2+1,n3)+F3(n2,n3)-S_x(3));
            F4_half(n2,n3) = 1/2*(F4(n2+1,n3)+F4(n2,n3)-S_x(4));
            
            
            % Direction y
            
            % Expressing rho, u, v and H in the basis of w (Roe's paper)
            rho_av_y = sqrt(rho(n2,n3)*rho(n2,n3+1));
            u_av_y = (sqrt(rho(n2,n3))*u(n2,n3)+sqrt(rho(n2,n3+1))*u(n2,n3+1))/(sqrt(rho(n2,n3))+sqrt(rho(n2,n3+1)));
            v_av_y = (sqrt(rho(n2,n3))*v(n2,n3)+sqrt(rho(n2,n3+1))*v(n2,n3+1))/(sqrt(rho(n2,n3))+sqrt(rho(n2,n3+1)));
            H_av_y = (sqrt(rho(n2,n3))*H(n2,n3)+sqrt(rho(n2,n3+1))*H(n2,n3+1))/(sqrt(rho(n2,n3))+sqrt(rho(n2,n3+1)));
            
            a_av_y = sqrt((gamma-1)*(H_av_y-1/2*(u_av_y^2 + v_av_y^2))); %Square root of the derivative of P by rho, used in the eigenvalues of the Jacobian matrix
            
            % Eigenvalues of the Jacobian matrix
            lambda_y(1) = v_av_y-a_av_y;
            lambda_y(2) = v_av_y;
            lambda_y(3) = v_av_y+a_av_y;
            lambda_y(4) = v_av_y;
            
            if (umax_y<max(abs(lambda_y)))
                umax_y = max(abs(lambda_y));
            end
            
            %R Matrix constructed from eigenvectors in w-space
            
            R_y = zeros(4,4);
            R_y(1,:) = 1;
            R_y(1,4) = 0;
            R_y(2,:) = u_av_y;
            R_y(2,4) = -1;
            R_y(3,:) = lambda_y;
            R_y(3,4) = 0;
            R_y(4,1) = H_av_y-v_av_y*a_av_y;
            R_y(4,2) = 1/2*(u_av_y^2 + v_av_y^2);
            R_y(4,3) = H_av_y+v_av_y*a_av_y;
            R_y(4,4) = -u_av_y;

            % Solving the linear problem
            
            A_y = zeros(4,4);
            for i =1:4
                A_y(i,i) = abs(lambda_y(i));
            end
            
            RAR_inv_y = R_y*A_y/R_y;
            
            S_y = zeros(4,1);
            S_y(1) = RAR_inv_y(1,1)*(U1(n2,n3+1)-U1(n2,n3)) + RAR_inv_y(1,2)*(U2(n2,n3+1)-U2(n2,n3)) + RAR_inv_y(1,3)*(U3(n2,n3+1)-U3(n2,n3)) + RAR_inv_y(1,4)*(U4(n2,n3+1)-U4(n2,n3));
            S_y(2) = RAR_inv_y(2,1)*(U1(n2,n3+1)-U1(n2,n3)) + RAR_inv_y(2,2)*(U2(n2,n3+1)-U2(n2,n3)) + RAR_inv_y(2,3)*(U3(n2,n3+1)-U3(n2,n3)) + RAR_inv_y(2,4)*(U4(n2,n3+1)-U4(n2,n3));
            S_y(3) = RAR_inv_y(3,1)*(U1(n2,n3+1)-U1(n2,n3)) + RAR_inv_y(3,2)*(U2(n2,n3+1)-U2(n2,n3)) + RAR_inv_y(3,3)*(U3(n2,n3+1)-U3(n2,n3)) + RAR_inv_y(3,4)*(U4(n2,n3+1)-U4(n2,n3));
            S_y(4) = RAR_inv_y(4,1)*(U1(n2,n3+1)-U1(n2,n3)) + RAR_inv_y(4,2)*(U2(n2,n3+1)-U2(n2,n3)) + RAR_inv_y(4,3)*(U3(n2,n3+1)-U3(n2,n3)) + RAR_inv_y(4,4)*(U4(n2,n3+1)-U4(n2,n3));
            
            G1_half(n2,n3) = 1/2*(G1(n2,n3+1)+G1(n2,n3)-S_y(1));
            G2_half(n2,n3) = 1/2*(G2(n2,n3+1)+G2(n2,n3)-S_y(2));
            G3_half(n2,n3) = 1/2*(G3(n2,n3+1)+G3(n2,n3)-S_y(3));
            G4_half(n2,n3) = 1/2*(G4(n2,n3+1)+G4(n2,n3)-S_y(4));
        end
    end
    
    if (mod(n1,ng) == 0)
        I = I+1;
        Per = n1/nt*100;
        disp(Per)
        figure(I);
        view(135,45)
        sur = surface(1:ny,1:nx,P);
        set(sur,'LineStyle','none')
        title('Pressure Colormap');
        ylabel('Position x');
        xlabel('Position y');
        axis ([0 inf 0 inf 0 50000000]);
        M(I) = getframe(); %#ok<SAGROW>
        close
    end
    figure(100)
    tiledlayout(4,1)
    
    nexttile
    plot(1:nx,rho(:,ny/2))
    
    nexttile
    plot(1:nx,P(:,ny/2))
    
    nexttile
    plot(1:nx,F4_half(:,ny/2))
    
    nexttile
    plot(1:nx,G4_half(:,ny/2))
%     
%     nexttile
%     plot(1:nt,power_array)
% 
%     nexttile
%     plot(1:nx,U4(:,ny/2))
        
end

% figure
% 
% figure(2)
% 
% contourf(1:nx,1:ny,P)
% title('Pressure Colormap')
% xlabel('Position x')
% ylabel('Position y')
% 
% figure(1)
% 
% sur = surface(1:ny,1:nx,P);
% set(sur,'LineStyle','none')
% title('Pressure Colormap')
% ylabel('Position x')
% xlabel('Position y')

% figure(2)
% 
% surface(1:ny,1:nx,u)
% title('X Speed Colormap')
% ylabel('Position x')
% xlabel('Position y')
% 
% figure(3)
% 
% surface(1:ny,1:nx,v)
% title('Y Speed Colormap')
% ylabel('Position x')
% xlabel('Position y')
% 
% figure(4)
% 
% surface(1:ny,1:nx,rho)
% title('Density Colormap')
% ylabel('Position x')
% xlabel('Position y')

toc

movie(M,15)
