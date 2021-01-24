% Title : 1D CFD (Two-temperatures model, T_vib and T_gas)
% Author: Kuniyoshi Tabata
% Created date：2006115

% clear all

%% Time step and spatial grid
% Spatial grid
    N = 600;
    dx = 0.5e-3;                              % 空間刻み：5 mm

% Time step
    N_t = 20000;

% CFL condition
    CFL = 0.01;


%% Physicsal constatnts

    m_e = 9.10938356e-31;                      % Electron mass
    q = 1.60217662e-19;                        % Elementary charge
    e = 8.85418782e-12;                        % Permittivity
    u = 4.0*3.1415e-7;                         % Permeability
    c = (1.0/(e*u))^0.5;                       % Speed of light
    kB = 1.380649e-23;                         % Boltzmann's constant
    Avogadro = 6.02e+23;                       % Avogadro's constant
    h = 6.62607015e-34;                        % Planck's constant

    m_O2 = 32/Avogadro*1e-3;                   % N2 mass
    m_N2 = 28/Avogadro*1e-3;                   % O2 mass
    ratio_N2 = 0.8;                            % Mass ratio of N2
    ratio_O2 = 0.2;                            % Mass ratio of O2
    m_air = m_N2*ratio_N2 + m_O2*ratio_O2;     % Air mass

    R_air = kB/m_air;                          % Gas constant
    
    e_vib_N2 = 0.288;        % 窒素分子の振動特性エネルギー [eV](2016_Electron-Vibrational energy exchange), 3380 K(実在気体流れの力学)
    e_vib_O2 = 0.196;        % 酸素分子の振動特性エネルギー [eV](2016_Electron-vibration relaxation in oxygen plasmas), 2230 K(実在気体流れの力学)

    
%% 電離，電子励起に関する定数 / Contants of ionziation and excitation
% 窒素分子に関する定数 / Constatnts of nitrogen molecules
    epsilon_ion_N2 = 15.58;    % 電離エネルギー / Ionization energy [eV]
    epsilon_ex_N2(1) = 6.17;   % A(v=0-4), 電子励起準位のエネルギー [eV]
    epsilon_ex_N2(2) = 7.0;    % A(v=5-9)
    epsilon_ex_N2(3) = 7.35;   % B
    epsilon_ex_N2(4) = 7.36;   % W
    epsilon_ex_N2(5) = 7.80;   % A(v>10)
    epsilon_ex_N2(6) = 8.16;   % B'
    epsilon_ex_N2(7) = 8.40;   % a'
    epsilon_ex_N2(8) = 8.55;   % a
    epsilon_ex_N2(9) = 8.89;   % w
    epsilon_ex_N2(10) = 11.03; % C
    epsilon_ex_N2(11) = 11.88; % E
    epsilon_ex_N2(12) = 12.25; % a"
    epsilon_ex_N2(13) = 13;    % sum of singlets


% 酸素分子に関する定数
    epsilon_ion_O2 = 12.07;   % 電離エネルギー / Ionization energy [eV]
    epsilon_ex_O2(1) = 0;     % 基底準位
    epsilon_ex_O2(1) = 13.118;% M
    epsilon_ex_O2(1) = 12.464;% H
    epsilon_ex_O2(1) = 9.918; % F
    epsilon_ex_O2(1) = 9.592; % E
    epsilon_ex_O2(1) = 6.119; % B
    epsilon_ex_O2(1) = 4.490; % c
    epsilon_ex_O2(1) = 4.340; % A
    epsilon_ex_O2(1) = 4.256; % C
    epsilon_ex_O2(1) = 1.627; % b
    epsilon_ex_O2(1)= 0.977;  % a
    

%% Numdrical conditions

% Ambient conditions
    P_ambient = 101300;                        % Ambient pressure [Pa]
    T_ambient = 293;                           % Ambient temperature [K] 
    n_ambient_air = P_ambient/(kB*T_ambient);  % Ambient neutral particles dendity [1/m^3]
    rho_ambient = m_air*n_ambient_air;         % Mass density of Air
    u_ambient = 0;                             % Flow velocity [m/s]
    gamma_gas = 7/5;                           % Specific heat ratio

% Propagation velocity and peak intensity
    S_peak = 1.9;                              % Peak power density [GW/m^2]
    U_ion = 419*S_peak-14.9;                 % Propagation velocity of ionization front [m/s]
    
% Length of ionization front propagation
    length = dx * N;
    
% Thruster configuration
    eta = 0.486225;                            % Heat absorption efficiency
    phi = 1;%0.4;                              % Plasma space occupancy
    omega_0 = 20.4e-3;                         % Beam radius [mm]
    d = 60e-3;                                 % Thruster's diameter [mm]
    
% Milliemter-wave and Plasma
    f = 170e+9;                                % Frequency
    omega = 2*pi*f;                            % Angular frequency
    lambda = c/f;                              % Wavelength
    
% Absorption layer
    u_m = 5.3e+9*P_ambient*760/101300;         % Collisional frequency between electrons and neutral particles
    n_cutoff = m_e*e/q^2*(omega^2+u_m^2);      % Cutoff density
    n_e = n_cutoff;
    r_abs = n_e/n_cutoff;
    q_abs = u_m/omega;
    k_0 = omega/c;
    l_abs = 1/(k_0/2^0.5*(r_abs-1+((1-r_abs^2)+(q_abs*r_abs)^2)^0.5)^0.5);   % Absorption depth (170 GHzの場合0.2mm程度)


%% Initial conditions
%--------------------------- Pre-allocation -----------------------------%
	n_e_d = zeros(1,N);
    
    rho_av = zeros(1,N);
    n_air = zeros(1,N);
    u_av = zeros(1,N);
    H_av = zeros(1,N);
    gamma_av = zeros(1,N);
    c_av = zeros(1,N);
    
    rho = zeros(1,N);
    u = zeros(1,N);
    u2 = zeros(1,N);                    % 電離波面固定系に対する2次側の流速
    P = zeros(1,N);
    T_gas = zeros(1,N);
    T_vib = zeros(1,N);
    H = zeros(1,N);
    gamma = zeros(1,N);                 % Heat capacity ratio => IMPORTANT

    
    Q = zeros(4,N);
    
    E = zeros(4,N);
    E_half = zeros(4,N);
    
    Q_VT = zeros(1,N);
    Q_ela = zeros(1,N);
    
    L = zeros(3,3,N);
    R = zeros(3,3,N);
    R_inv = zeros(3,3,N);
    
    w = zeros(1,N);
    
    d_shock = zeros(1,2);
    d_ion = zeros(1,2);
    
    % Output for graphix
    Energy_gas = zeros(1,N);
    Energy_pressure = zeros(1,N);
    Energy_vib = zeros(1,N);
    Energy_flow = zeros(1,N);
    
    
%------------------------ Initial conditions ----------------------------%
    
    rho(:) = rho_ambient;                                        % Mass density
    n_air(:) = n_ambient_air;                                    % Particle density
    u(:) = u_ambient;                                            % Flow velocity
    P(:) = P_ambient;                                            % Pressure of heavy particles
    T_gas(:) = T_ambient;                                        % Translational and rotational temperature
    T_vib(:) = T_ambient;                                        % Vibrational temperature
    e_gas = 5/2*R_air.*T_gas;                                    % Specific translational and rotational energy
    e_vib = R_air*q*e_vib_N2/kB./(exp(q*e_vib_N2/kB./T_vib)-1);  % Specific vibrational energy
    H = e_gas + 0.5*u.^2 + P./rho;                               % Specific enthalpy
    gamma = (e_gas+e_vib+R_air*T_gas)./(e_gas+e_vib);            % Specific heat ratio
   
    Q(1,:) = rho;
    Q(2,:) = rho.*u;
    Q(3,:) = rho.*e_gas + 0.5*rho.*u.^2;
    Q(4,:) = rho.*e_vib;
    
    E(1,:) = Q(2,:);
    E(2,:) = rho.*u.^2 + P;
    E(3,:) = (Q(3,:)+P).*u;
    E(4,:) = Q(4,:).*u;
    
    i = 1:N;                                  % x axis set
    X = i*dx;
    
    
%% ----------------------- MAIN Program ------------------------- %%
 % Time count
    t = 0;
    u_max = U_ion;

for step = 0:N_t
    
   %% Termination conditon
   
    if U_ion *t < dx*N 
   
        
   %% Time step update
    
    dt = CFL*dx/u_max;
    t = t+dt;
    
    if mod(t,1e-5) < dt
        t*1e+3  % 10 usごとに時間を表示[ms]
    end
    
   %% The placement of shock wave => MSCのときこの方法はうまくいかないので改善が必要
   %if mod(t,1e-6) < dt
   %    d_shock = [d_shock; t*1e+6 max(find(P == max(P))*dx)];
   %    d_ion = [d_ion; t*1e+6 max(find(T_vib == max(T_vib))*dx)];
   %end
   
   
   %% Calculation of values between cells => Roe average
        
        rho_av(1:N-1) = (rho(1:N-1).*rho(2:N)).^0.5;
        u_av(1:N-1) = (rho(1:N-1).^0.5.*u(1:N-1) + rho(2:N).^0.5.*u(2:N))./(rho(1:N-1).^0.5+rho(2:N).^0.5);
        H_av(1:N-1) = (rho(1:N-1).^0.5.*H(1:N-1) + rho(2:N).^0.5.*H(2:N))./(rho(1:N-1).^0.5+rho(2:N).^0.5);
        gamma_av(1:N-1) = (gamma(1:N-1).*gamma(2:N)).^0.5;
        c_av(1:N-1) = ((gamma_gas-1).*(H_av(1:N-1)-0.5*u_av(1:N-1).^2)).^0.5;
       
        rho_av(N) = rho_av(N-1);
        u_av(N) = u_av(N-1);
        H_av(N) = H_av(N-1);
        gamma_av(N) = gamma_av(N-1);
        c_av(N) = c_av(N-1);
    
    
   %% Calculation of matrix, R
        
        R(1,:,:) = 1;
        R(2,1,:) = u_av-c_av;
        R(2,2,:) = u_av;
        R(2,3,:) = u_av+c_av;
        R(3,1,:) = H_av - u_av.*c_av;
        R(3,2,:) = 0.5*u_av.^2;
        R(3,3,:) = H_av + u_av.*c_av;
    
    
   %% Calculation of inverse matrix, R^(-1)
        
        B = (gamma_gas-1)./c_av.^2;
        A = B.*u_av.^2/2;
    
        R_inv(1,1,:) = 0.5 * ( A + u_av./c_av );
        R_inv(1,2,:) = 0.5 * ( -1*B.*u_av - 1./c_av);
        R_inv(1,3,:) = 0.5 * B;
    
        R_inv(2,1,:) = 1-A;
        R_inv(2,2,:) = B.*u_av;
        R_inv(2,3,:) = -1*B;
    
        R_inv(3,1,:) = 0.5 * ( A - u_av./c_av );
        R_inv(3,2,:) = 0.5 * ( -1*B.*u_av + 1./c_av);
        R_inv(3,3,:) = 0.5 * B;
    
    
   %% Calculation of diagonal matrix
        
        L(1,1,:) = u_av - c_av;
        L(2,2,:) = u_av;
        L(3,3,:) = u_av + c_av;
    
    if u_max < max(abs(L),[],'all')
        u_max = max(abs(L),[],'all');
    end
    
    
   %% Calculation of numerical flux between cells
    
    for i = 1:N-1
        E_half(1:3,i) = 0.5 * ( E(1:3,i)+E(1:3,i+1) - R(1:3,1:3,i)*abs(L(1:3,1:3,i))*R_inv(1:3,1:3,i)*(Q(1:3,i+1)-Q(1:3,i)));
        E_half(4,i) = 0.5 * ( E(4,i+1)+E(4,i) - abs(u(i)).*(Q(4,i+1)-Q(4,i)));
    end
        E_half(1:3,N) = E_half(1:3,N-1);
        E_half(4,N) = E_half(4,N-1);
    
    
   %% Energy source term => 吸収長を波長程度で与える
   
   % Heat input by ionization front propagation
        % Source term
        for i = 1:N
            if t*U_ion < length
                if (i-1)*dx <= t*U_ion && i*dx >= t*U_ion
                    w(i) = eta * 2 *(omega_0/d)^2 * S_peak*1e+9 /dx;
                else
                    w(i) = 0;
                end
            else
                w(:) = 0;
            end            
        end
    
    % Vibrational-translational relaxation energy
        % k_2 の計算
            k_VT = 7.8e-18.*T_gas.*exp(-218./(T_gas.^(1/3))+690./T_gas)./(1-exp(-1*q*e_vib_N2./(kB.*T_gas)));
        
        % 振動エネルギーの計算
            E_v_1 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_vib))-1); % [J]
            E_v_3 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_gas))-1); % [J]
            
        % エネルギーソース項の計算
            Q_VT = n_air.*n_air.*k_VT.*(1-exp(-1*q*e_vib_N2./(kB.*T_gas))).*(E_v_1-E_v_3); % [J/m^3]
            
    % 運動量輸送エネルギー / Elastic collision energy transfer
    % とりあえずサハの式から電子数密度を算出する．
            n_e_d = (n_air*4.82e+21.*T_vib.^1.5.*exp(-11600*epsilon_ion_N2./T_vib)).^0.5;
            Q_ela = 3*m_e/m_air*kB.*5.3e+9.*P*760/101300.*(T_vib-T_gas).*1e+21*P_ambient/101300;            % 全領域で1e+21にしているので，緩和が速い可能性がある．
    

   %% Solve Euler equation
    % Time evolution
       % Souce term
        Q_electron = w/phi - Q_VT - Q_ela;
        Q_gas = Q_VT + Q_ela;
        
       % Time evolution
        Q(1,2:N-1) = Q(1,2:N-1) - dt/dx * (E_half(1,2:N-1) - E_half(1,1:N-2));
        Q(2,2:N-1) = Q(2,2:N-1) - dt/dx * (E_half(2,2:N-1) - E_half(2,1:N-2));
        Q(3,2:N-1) = Q(3,2:N-1) - dt/dx * (E_half(3,2:N-1) - E_half(3,1:N-2)) + Q_gas(2:N-1) * dt;          % 並進・回転エネルギー保存則
        Q(4,2:N-1) = Q(4,2:N-1) - dt/dx * (E_half(4,2:N-1) - E_half(4,1:N-2)) + Q_electron(2:N-1) * dt;     % 振動エネルギー保存則, 空間占有率を考慮して入熱
        
    % Boundary conditions (Left): Closed boundary
        Q(1,1) = Q(1,2);
        Q(2,1) = 0;
        Q(3,1) = Q(3,2);
        Q(4,1) = Q(4,2);
        
    % Boundary conditions (Right): Open boundary
        T_temp = (gamma_gas-1)/R_air/Q(1,N-1)*(Q(3,N-1)-0.5*Q(2,N-1)^2/Q(1,N-1));
        Q(1,N) = P_ambient/R_air/T_temp;
        Q(2,N) = Q(2,N-1)/Q(1,N-1)*Q(1,N);
        Q(3,N) = Q(3,N-1)/Q(1,N-1)*Q(1,N);
        Q(4,N) = Q(4,N-1)/Q(1,N-1)*Q(1,N);
       
        
    %% Update matrix
        
        rho = Q(1,:);
        n_air = rho/m_air;
        u = Q(2,:)./rho;
        u2 = U_ion-u;
        M2 =  u2./c_av;
        e_gas = Q(3,:)./rho - 0.5.*u.^2;
        e_vib = Q(4,:)./rho;
        T_gas = e_gas/(5/2*R_air);
        T_vib = e_vib_N2*q/kB./log(R_air*e_vib_N2*q/kB./e_vib+1);
        P = rho.*R_air.*T_gas;
        H = (Q(3,:)+P)./rho;
        gamma = (Q(3,:)+Q(4,:)+P)./(Q(3,:)+Q(4,:));
        
        
        E(1,:) = Q(2,:);
        E(2,:) = rho.*u.^2 + P;
        E(3,:) = (Q(3,:)+P).*u;
        E(4,:) = Q(4,:).*u;
        
    end
    
end

   % Output for Graphix
    Energy_pressure = P - P_ambient;
    Energy_flow = 0.5*rho.*u.^2;
    Energy_gas = Q(3,:) - 0.5*rho_ambient.*u.^2 - rho.*5/2*R_air.*T_ambient;
    Energy_vib = Q(4,:) - rho_ambient*R_air*q*e_vib_N2/kB/(exp(q*e_vib_N2/kB/T_ambient)-1);
    
    Energy_1 = Energy_pressure;
    Energy_2 = Energy_1 + Energy_flow;
    Energy_3 = Energy_2 + Energy_gas;
    Energy_4 = Energy_3 + Energy_vib;