% Title : 2D CFD (One-temperature model)
% Time step: 1st explicit, Spatial step: MUSCL method, minmod limiter
% Author: Kuniyoshi Tabata
% Created date：200825

% limiter_1_3はMUSCL法ではなく，1次精度の計算
% limiter_1_2はMUSCL法を用いた方法
% limiter_1_1はヤコビアンが間違っている<=宇宙研・谷口君が指摘してくれた．

% E_halfの計算=> u からu_avに変更
% 入熱量の計算 => totalで半分しか入っていない可能性

% 200929 衝撃波固定系で斜め衝撃波関係式が成立するかどうか

% すでにある関数と同じ変数を定義するとエラーを吐くことがあるので注意．maxなど

% 201023 衝撃波固定系で計算を開始．まだ，流れを上流から流入させるだけの計算になっている．収束なども確認できるコードとなっている．
% 201203 VT緩和項の改善．dissociationに関連する加熱項，fast gas heating
% に関連する項を加えた．さらに，計算領域を小さくした．プラズマ構造の挿入．
% 201207 計算領域を100 mmに直した

% 構造ありなしで変更する点．①structure == の選択　②ビーム直径 27 or 28mm ③Plasam space occupancy 

 clear all


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
    m_O = 16/Avogadro*1e-3;                    % N2 mass
    m_N = 14/Avogadro*1e-3;                    % O2 mass
    ratio_N2 = 0.8;                            % Ratio of N2
    ratio_O2 = 0.2;                            % Ratio of O2
    %ratio_N = 1e-4;                            % Ratio of N
    %ratio_O = 1e-4;                            % Ratio of O
    m_air = m_N2*ratio_N2 + m_O2*ratio_O2;     % Air mass

    R_air = kB/m_air;                          % Gas constant
    gamma = 7/5;                               % Heat capacity ratio
    
    e_vib_N2 = 0.288;                          % 窒素分子の振動特性エネルギー [eV](2016_Electron-Vibrational energy exchange), 3380 K(実在気体流れの力学)
    e_vib_O2 = 0.196;                          % 酸素分子の振動特性エネルギー [eV](2016_Electron-vibration relaxation in oxygen plasmas), 2230 K(実在気体流れの力学)
    e_D = 9.76;                                % 窒素分子の解離エネルギー [eV]

    
    
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

    
    
%% Millimeter-wave and plasma structure
% Millimeter-wave
    f = 170e+9;                                     % Beam frequency
    omega = 2*pi*f;                                 % Beam angular frequency
    lambda = c/f;                                   % Beam wavelength
    w_beam = 40.8e-3;                               % Beam diameter of flat-top beam

% Comb-shaped structure
    w_plasma = 27e-3;                               % Plasma width of comb-shaped structure (拡散構造：ここをw_beamとする．)
    pitch = 0;                                      % Pitch between plasmoids
    phi = 1.0;                                      % Plasma space occupancy
    %N_beam = 31;%31;
    %N_plasma = 31;%31;
    %N_pitch = 0;
    
% Thruster configuration
    d_tube = 100e-3;                                % Tube diameter [mm]
    l_tube = 100e-3;                                % Tube length [mm]
    
% Absorption layer
    P_ambient = 101300;
    u_m = 5.3e+9*P_ambient*760/101300;              % Collisional frequency between electrons and neutral particles
    n_cutoff = m_e*e/q^2*(omega^2+u_m^2);           % Cutoff density
    n_e = n_cutoff;
    r_abs = n_e/n_cutoff;
    q_abs = u_m/omega;
    k_0 = omega/c;
    %l_abs = 1/(k_0/2^0.5*(r_abs-1+((1-r_abs^2)+(q_abs*r_abs)^2)^0.5)^0.5);   % Absorption depth (170 GHzの場合0.2mm程度)

% Absorption rate of millimeter-wave energy
    eta = 0.5;%1.0;                                 % Heat absorption efficiency
    l_abs = 1.0e-3;%lambda/2;                       % Absorption length
    
% Propagation velocity and peak intensity
    S_peak = 2.2;                                   % Peak power density [GW/m^2]
    S_local = eta*S_peak/phi;                       % Local beam intensity corrected by energy concentration [GW/m^2]
    U_ioniz = 419*S_peak-14.9;                      % Propagation velocity of ionization front [m/s]



%% Time step and spatial grid
% Spatial grid
    dx = 1.0e-3;                                    % Spatial grid in x direction
    dy = 1.0e-3;                                    % Spatial grid in y direction
    Nx = round(l_tube/dx);                          % Grid number in x direction
    Ny = round(d_tube/dy);                          % Grid number in y direction
 
% Time step
    dt = 0;                                         % Time step
    Nt = 50000;                                     % Maximum time step

% CFL condition
    CFL = 0.2;                                      % Courant number
    
% 
    N_plasma = w_plasma/dy;                         % Plasma region in y direction
    N_plasma_x = Nx*0.9;                            % Plasma region in x direction
    N_pitch = round(pitch/dy);                      % Pitch between plasmoids
    
    N_abs = round(l_abs/dx);              % Region of milliemter-wave absorption
    if N_abs == 0
        N_abs = 1;
    end

% MUSCL method order
    mu_muscl = 1/3;
    omega_limiter = 1.0;                       % Compression parameter for limiter function


%% Numdrical conditions
% Ambient conditions
    P_ambient = 101300;                             % Ambient pressure [Pa]
    T_ambient = 293;                                % Ambient temperature [K] 
    n_ambient_air = P_ambient/(kB*T_ambient);       % Ambient neutral particles dendity [1/m^3]
    rho_ambient = m_air*n_ambient_air;              % Mass density of Air
    u_ambient = -U_ioniz;                           % Flow velocity (x direction) [m/s]
    v_ambient = 0;                                  % Flow velocity (y direction) [m/s]
    


%% Initial conditions
%--------------------------- Pre-allocation -----------------------------%

    n_e_d = zeros(1,1,Ny,Nx);
    n_air = zeros(1,1,Ny,Nx);
    
	rho_av = zeros(1,1,Ny,Nx);
    u_av = zeros(1,1,Ny,Nx);
    v_av = zeros(1,1,Ny,Nx);
    H_av = zeros(1,1,Ny,Nx);
    c_av = zeros(1,1,Ny,Nx);
    
    rho = zeros(1,1,Ny,Nx);
    u = zeros(1,1,Ny,Nx);
    v = zeros(1,1,Ny,Nx);
    P = zeros(1,1,Ny,Nx);
    T_gas = zeros(1,1,Ny,Nx);
    T_vib = zeros(1,1,Ny,Nx);
    e_gas = zeros(1,1,Ny,Nx);
    e_vib = zeros(1,1,Ny,Nx);
    H = zeros(1,1,Ny,Nx);
    
    T_temp = zeros(1,1,Ny);
    
    b2 = zeros(1,1,Ny,Nx);
    b1 = zeros(1,1,Ny,Nx);
    
    Q = zeros(5,1,Ny,Nx);
    
    Q_L = zeros(5,1,Ny,Nx);
    Q_R = zeros(5,1,Ny,Nx);
    Delta_minus = zeros(5,1,Ny,Nx);
    Delta_plus  = zeros(5,1,Ny,Nx);
    
    sign_1 = zeros(5,1,Ny,Nx);
    sign_2 = zeros(5,1,Ny,Nx);
    sign_3 = zeros(5,1,Ny,Nx);
    sign_4 = zeros(5,1,Ny,Nx);
    
    SIGN_1 = zeros(5,1,Ny,Nx);
    SIGN_2 = zeros(5,1,Ny,Nx);
    
    minmod_1 = zeros(5,1,Ny,Nx);
    minmod_2  = zeros(5,1,Ny,Nx);
    
    E = zeros(5,1,Ny,Nx);
    E_half = zeros(5,1,Ny,Nx);
    E_L = zeros(5,1,Ny,Nx);
    E_R = zeros(5,1,Ny,Nx);
    
    F = zeros(5,1,Ny,Nx);
    F_half = zeros(5,1,Ny,Nx);
    F_L = zeros(5,1,Ny,Nx);
    F_R = zeros(5,1,Ny,Nx);
    
    L = zeros(4,4,Ny,Nx);
    R = zeros(4,4,Ny,Nx);
    R_inv = zeros(4,4,Ny,Nx);
    
    W = zeros(1,1,Ny,Nx);
    E_Td = zeros(1,1,Ny,Nx);
    eta_vib = zeros(1,1,Ny,Nx);
    eta_gas = zeros(1,1,Ny,Nx);
    
    Q_VT = zeros(1,1,Ny,Nx);
    Q_ela = zeros(1,1,Ny,Nx);
    Q_electronic = zeros(1,1,Ny,Nx);
    
    X = zeros(1,Nx);
    Y = zeros(1,Ny);

    
%--------------------------- For graphix -----------------------------%

    rho_plot = zeros(Ny,Nx);
    T_gas_plot = zeros(Ny,Nx);
    T_vib_plot = zeros(Ny,Nx);
    P_plot = zeros(Ny,Nx);
    u_plot = zeros(Ny,Nx);
    v_plot = zeros(Ny,Nx);
    W_plot = zeros(Ny,Nx);
    a_plot = zeros(Ny,Nx);          % 音速
    M_plot = zeros(Ny,Nx);
    
    Nx_min = 1;
    N_max = zeros(1,Ny);
    x_p = zeros(1,Ny);
    x_v = zeros(1,Ny);
    
    u_rel = zeros(Ny,Nx);
    v_rel = zeros(Ny,Nx);
    M1 = zeros(1,Ny);
    M2_plot = zeros(Ny,Nx);
    M2_plot_theory = zeros(Ny,Nx);
    Error = zeros(Ny,Nx);
    
    P_closed = zeros(1,Nx);
    T_gas_closed = zeros(1,Nx);
    T_vib_closed = zeros(1,Nx);
    rho_closed = zeros(1,Nx);
    u_closed = zeros(1,Nx);
    v_closed = zeros(1,Nx);
    
    P_wall = zeros(1,2);
    Distance = zeros(1,Ny+1);
    Norm = zeros(1,7);
    
    Q_VT_plot = zeros(Ny,Nx);
    Q_D_plot = zeros(Ny,Nx);
    
    filename_excel_1 = 'Pressure_wall';
    filename_excel_2 = 'Distance';
    filename_excel_3 = 'Norm';
    
%------------------------ Initial conditions ----------------------------%
    
    rho(:) = rho_ambient;                                 % Mass density
    n_air(:) = n_ambient_air;
    n_N2 = ratio_N2*n_air;
    n_O2 = ratio_O2*n_air;
    u(:) = u_ambient;                                     % Flow velocity (x direction)
    v(:) = v_ambient;                                     % Flow velocity (y direction)
    P(:) = P_ambient;                                     % Pressure of heavy particles
    T_gas(:) = T_ambient;                                 % Translational and rotational temperature
    T_vib(:) = T_ambient;                                 % Vibrational temperature
    e_gas(:) = 5/2*R_air*T_ambient;
    e_vib(:) = 1./rho.* ( n_N2*q*e_vib_N2./(exp(q*e_vib_N2/kB./T_vib)-1) + n_O2*q*e_vib_O2./(exp(q*e_vib_O2/kB./T_vib)-1)) ;  % Specific vibrational energy
    H(:) = gamma/(gamma-1)*R_air*T_gas+0.5*(u.^2+v.^2);   % Enthalpy
   
    Q(1,1,:,:) = rho;
    Q(2,1,:,:) = rho.*u;
    Q(3,1,:,:) = rho.*v;
    Q(4,1,:,:) = P/(gamma-1)+0.5*rho.*(u.^2+v.^2);
    Q(5,1,:,:) = rho.*e_vib;
    
    E(1,1,:,:) = Q(2,1,:,:);
    E(2,1,:,:) = rho.*u.^2 + P;
    E(3,1,:,:) = rho.*u.*v;
    E(4,1,:,:) = (Q(4,1,:,:)+P).*u;
    E(5,1,:,:) = Q(5,1,:,:).*u;
    
    F(1,1,:,:) = Q(3,1,:,:);
    F(2,1,:,:) = rho.*u.*v;
    F(3,1,:,:) = rho.*v.^2 + P;
    F(4,1,:,:) = (Q(4,1,:,:)+P).*v;
    F(5,1,:,:) = Q(5,1,:,:).*v;
    
    i = 1:Nx;                                  % x axis set
    j = 1:Ny;                                  % y axis set
    X = i*dx;
    Y = j*dy;
   
% Relaxation time for reactions O2 = [O2-O2 O2-N2 O2-O O2-N], N2 = [N2-N2 N2-O2 N2-O]
    mu_O2 = Avogadro * 1e+3 * [m_O2*m_O2/(m_O2+m_O2) m_O2*m_N2/(m_O2+m_N2) m_O2*m_O/(m_O2+m_O) m_O2*m_N/(m_O2+m_N)];
    mu_N2 = Avogadro * 1e+3 * [m_N2*m_N2/(m_N2+m_N2) m_N2*m_O2/(m_N2+m_O2) m_N2*m_O/(m_N2+m_O)];
    
    A_O2 = 1.16e-3 * mu_O2.^(1/2) * (e_vib_O2*q/kB)^(4/3);
    A_N2 = 1.16e-3 * mu_N2.^(1/2) * (e_vib_N2*q/kB)^(4/3);
    
    phi_O2 = [m_O2/(m_O2+m_N2+m_O+m_N) m_N2/(m_O2+m_N2+m_O+m_N) m_O/(m_O2+m_N2+m_O+m_N) m_N/(m_O2+m_N2+m_O+m_N)];
    phi_N2 = [m_N2/(m_N2+m_O2+m_O) m_O2/(m_N2+m_O2+m_O) m_O/(m_N2+m_O2+m_O)];
    
    t_VT_O2 = zeros(1,1,Ny,Nx);   % Relaxation time
    t_VT_N2 = zeros(1,1,Ny,Nx);
    
    Q_VT_O2 = zeros(1,1,Ny,Nx);   % Energy transfer
    Q_VT_N2 = zeros(1,1,Ny,Nx);
    
    Q_vib = zeros(1,1,Ny,Nx);
    Q_gas = zeros(1,1,Ny,Nx);
    
    
%% ----------------------- MAIN Program ------------------------- %%
 % Time count
    k = 0;
    t_initial = 0;
    t = t_initial;
    u_max = U_ioniz;
    dt = 0.5*lambda/u_max;  % 着火点が壁から1λ離れているとする．
    
    ii = 1;
    
    step = 1
    
%% ----------------------- MAIN Program ------------------------- %%
for step = 1:Nt
  
    step
    t

    
   %% Fractional Time Step
   
   for k = 1:3
       
    if k == 2
        dt = CFL / max(abs(u_av)/dx+abs(v_av)/dy+c_av*(1/dx^2+1/dy^2)^0.5,[],'all');
    end
   
    if k == 1 | k == 3       % x方向の時間発展 
   
   %% Calculation of values between cells
        
        rho_av(1,1,:,1:Nx-1) = ( rho(1,1,:,1:Nx-1) .* rho(1,1,:,2:Nx) ).^0.5;
        u_av(1,1,:,1:Nx-1) = (rho(1,1,:,1:Nx-1).^0.5.*u(1,1,:,1:Nx-1) + rho(1,1,:,2:Nx).^0.5.*u(1,1,:,2:Nx))./(rho(1,1,:,1:Nx-1).^0.5+rho(1,1,:,2:Nx).^0.5);
        v_av(1,1,:,1:Nx-1) = (rho(1,1,:,1:Nx-1).^0.5.*v(1,1,:,1:Nx-1) + rho(1,1,:,2:Nx).^0.5.*v(1,1,:,2:Nx))./(rho(1,1,:,1:Nx-1).^0.5+rho(1,1,:,2:Nx).^0.5);
        H_av(1,1,:,1:Nx-1) = (rho(1,1,:,1:Nx-1).^0.5.*H(1,1,:,1:Nx-1) + rho(1,1,:,2:Nx).^0.5.*H(1,1,:,2:Nx))./(rho(1,1,:,1:Nx-1).^0.5+rho(1,1,:,2:Nx).^0.5);
        c_av(1,1,:,1:Nx-1) = ((gamma-1)*(H_av(1,1,:,1:Nx-1)-0.5*(u_av(1,1,:,1:Nx-1).^2+v_av(1,1,:,1:Nx-1).^2))).^0.5;
    
        rho_av(1,1,:,Nx) = rho_av(1,1,:,Nx-1);
        u_av(1,1,:,Nx) = u_av(1,1,:,Nx-1);
        v_av(1,1,:,Nx) = v_av(1,1,:,Nx-1);
        H_av(1,1,:,Nx) = H_av(1,1,:,Nx-1);
        c_av(1,1,:,Nx) = c_av(1,1,:,Nx-1);
    
    
   %% Calculation of matrix, R
    
        R(1,1,:,:) = 1; 
        R(2,1,:,:) = u_av-c_av;
        R(3,1,:,:) = v_av;
        R(4,1,:,:) = H_av-u_av.*c_av;
        
        R(1,2,:,:) = 1;
        R(2,2,:,:) = u_av;
        R(3,2,:,:) = v_av;
        R(4,2,:,:) = (u_av.^2+v_av.^2)/2;
        
        R(1,3,:,:) = 1;
        R(2,3,:,:) = u_av+c_av;
        R(3,3,:,:) = v_av;
        R(4,3,:,:) = H_av+u_av.*c_av;
        
        R(1,4,:,:) = 0;
        R(2,4,:,:) = 0;
        R(3,4,:,:) = 1;
        R(4,4,:,:) = v_av;

        
    %% Calculation of inverse matrix, R^(-1)
        
        b2 = (gamma-1)./c_av.^2;
        b1 = b2.*(u_av.^2+v_av.^2)/2;
    
        R_inv(1,1,:,:) = (b1+u_av./c_av)/2;
        R_inv(2,1,:,:) = 1-b1;
        R_inv(3,1,:,:) = (b1-u_av./c_av)/2;
        R_inv(4,1,:,:) = -v_av;
        
        R_inv(1,2,:,:) = (-b2.*u_av-1./c_av)/2;
        R_inv(2,2,:,:) = b2.*u_av;
        R_inv(3,2,:,:) = (-b2.*u_av+1./c_av)/2;
        R_inv(4,2,:,:) = 0;
        
        R_inv(1,3,:,:) = -b2.*v_av/2;
        R_inv(2,3,:,:) = b2.*v_av;
        R_inv(3,3,:,:) = -b2.*v_av/2;
        R_inv(4,3,:,:) = 1;
        
        R_inv(1,4,:,:) = b2/2;
        R_inv(2,4,:,:) = -b2;
        R_inv(3,4,:,:) = b2/2;
        R_inv(4,4,:,:) = 0;
    
    
   %% Calculation of diagonal matrix
        
        L(1,1,:,:) = u_av - c_av;
        L(2,2,:,:) = u_av;
        L(3,3,:,:) = u_av + c_av;
        L(4,4,:,:) = u_av;

   %% Calculation of numerical flux between cells
    
    for j = 1:Ny
        for i = 1:Nx-1
            E_half(1:4,1,j,i) = 0.5 * ( E(1:4,1,j,i)+E(1:4,1,j,i+1) - R(1:4,1:4,j,i)*abs(L(1:4,1:4,j,i))*R_inv(1:4,1:4,j,i)*(Q(1:4,1,j,i+1)-Q(1:4,1,j,i)));
            E_half(5,1,j,i) = 0.5 * ( E(5,1,j,i+1)+E(5,1,j,i) - abs(u_av(1,1,j,i)).*(Q(5,1,j,i+1)-Q(5,1,j,i)));
        end
        E_half(:,:,j,Nx) = E_half(:,:,j,Nx-1);
    end
    
    
   %% Heat input by ionization front propagation % 変更必要
    
    % Source term
        W(:,:,:,:) = 0;
        
        structure = 0;  % 0: no structure, gaussian, 1: no structure, flat-top, 2: structure, flat-top
        time = lambda/U_ioniz;
        
        % No plasma structure
        if structure == 0
            for i = N_plasma_x-N_abs+1 : N_plasma_x
                for j = 1:N_plasma
                    W(1,1,Ny/2-(N_plasma-1)/2-1+j,i) = S_local*1e+9 /(N_abs*dx) * exp(-2*(((Ny/2-(N_plasma-1)/2)-1+j-Ny/2)*dx)^2/(w_beam/2)^2); % Gaussian beam
                end
            end
            %for i = N_plasma_x-N_abs+1 : N_plasma_x
            %    for j = 1:N_plasma
            %        W(1,1,Ny/2-(N_plasma-1)/2-1+j,i) = S_local*1e+9 / (N_abs*dx);  % Flat-top beam
            %    end
            %end
        
        elseif structure == 1
            if t <= time
                for i = N_plasma_x-N_abs+1 : N_plasma_x
                    W(1,1,38:61,i) = S_local*1e+9 / (N_abs*dx)*t/time;  % Flat-top beam
                end
            else
                for i = N_plasma_x-N_abs+1 : N_plasma_x
                    W(1,1,38:61,i) = S_local*1e+9 / (N_abs*dx);  % Flat-top beam
                end
            end
        elseif structure == 2
            if t <= time
                for i = N_plasma_x-N_abs+1 : N_plasma_x
                    for j = Ny/2-N_plasma/2:Ny/2+N_plasma/2
                        if mod(j,4)==0
                            W(1,1,j,i) = S_local*1e+9 / (N_abs*dx)*t/time  * exp(-2*((j-Ny/2)*dx)^2/(w_beam/2)^2);  % Flat-top beam
                            W(1,1,j+1,i) = S_local*1e+9 / (N_abs*dx)*t/time  * exp(-2*((j-Ny/2)*dx)^2/(w_beam/2)^2);  % Flat-top beam
                        end
                    end
                end
             else
                for i = N_plasma_x-N_abs+1 : N_plasma_x
                    for j = Ny/2-N_plasma/2:Ny/2+N_plasma/2
                        if mod(j,4)==0
                            W(1,1,j,i) = S_local*1e+9 / (N_abs*dx)  * exp(-2*((j-Ny/2)*dx)^2/(w_beam/2)^2);  % Flat-top beam
                            W(1,1,j+1,i) = S_local*1e+9 / (N_abs*dx)  * exp(-2*((j-Ny/2)*dx)^2/(w_beam/2)^2);  % Flat-top beam
                        end
                    end
                end
             end
        end
        
        % Townsend number calc.
        %     E_Td(:,:,:,:) = 0;
        %for i = N_plasma_x-N_abs+1 : N_plasma_x
        %    for j = 1:N_plasma
        %        E_Td(1,1,Ny/2-(N_plasma-1)/2-1+j,i) = ( 377* S_local*1e+9 * exp(-2*(((Ny/2-(N_plasma-1)/2)-1+j -Ny/2)*dx)^2/(w_beam/2)^2))^0.5/n_ambient_Air*1e+21;
        %        E_Td(1,1,Ny/2-(N_plasma-1)/2-1+j,i) = E_Td(1,1,Ny/2-(N_plasma-1)/2-1+j,N_plasma_x) / (1 + (omega/u_m)^2)^0.5;
        %    end
        %end
            E_Td = (377 * (W *phi/eta*(N_abs*dx))).^0.5 /n_ambient_air * 1e+21;
            E_Td = E_Td / (1 + (omega/u_m)^2)^0.5;
        
        % Eta_vib calc. => Hancoの修正版，自分で計算したもの，fast gas heatingを含む
        eta_vib(:,:,:,:) = 0;
        eta_gas(:,:,:,:) = 0;
        for j = 1:Ny
            for i = 1:Nx
                if E_Td(1,1,j,i) < 100
                    eta_vib(1,1,j,i) = 0.0000012216 * E_Td(1,1,j,i)^4 - 0.0002533235 * E_Td(1,1,j,i)^3 + 0.0096548449 * E_Td(1,1,j,i)^2 + 0.0107151002 * E_Td(1,1,j,i) + 94.9440890255 ;
                else
                    eta_vib(1,1,j,i) = -0.0000130955 * E_Td(1,1,j,i)^3 + 0.0090790427 * E_Td(1,1,j,i)^2 - 2.2469804542 * E_Td(1,1,j,i) + 208.4585917809 ; 
                end
            
                if E_Td(1,1,j,i) < 100
                    eta_gas(1,1,j,i) = -0.0000003598 * E_Td(1,1,j,i)^4 + 0.0000738890 * E_Td(1,1,j,i)^3 - 0.0026803169 * E_Td(1,1,j,i)^2 - 0.0106442568 * E_Td(1,1,j,i) + 1.3570686444 ;
                else
                    eta_gas(1,1,j,i) = 0.0000040816 * E_Td(1,1,j,i)^3 - 0.0028546398 * E_Td(1,1,j,i)^2 + 0.6946713755 * E_Td(1,1,j,i) - 33.6177117334 ; 
                end
            end
        end
        
        %{
        eta_vib(:,:,:,:) = 0;
        for j = 1:Ny
            if E_Td(1,1,j,70) < 100
                eta_vib(1,1,j,70) = 0.000000028435*E_Td(1,1,j,70)^4 - 0.0000047365*E_Td(1,1,j,70)^3 + 0.00011513*E_Td(1,1,j,70)^2 + 0.0027434*E_Td(1,1,j,70)+0.88884;
            else
                eta_vib(1,1,j,70) = 3874.3*E_Td(1,1,j,70)^(-1.9855);
            end
        end
        %}
        
        %{
    % Vibrational-translational relaxation energy
        % k_2 の計算
            k_VT = 7.8e-18.*T_gas.*exp(-218./(T_gas.^(1/3))+690./T_gas)./(1-exp(-1*q*e_vib_N2./(kB.*T_gas)));
            k_VT_2 = 2.3e-19 * exp(-1280./T_gas) + 2.7e-17 * exp(-10840./T_gas);
            t_VT = ((1-exp(-1*q*e_vib_N2./(kB.*T_gas))) .* (k_VT.*n_air + k_VT_2.*1e-4.*n_air)).^(-1);   % N原子が1e-4の割合で存在することを仮定
        
        % 振動エネルギーの計算
            E_v_1 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_vib))-1); % [J]
            E_v_3 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_gas))-1); % [J]
            
        % エネルギーソース項の計算
            Q_VT = n_air.*(E_v_1-E_v_3)./t_VT; % [J/m^3]
        %}
        %{
    % Vibrational-translational relaxation energy (PATTERN 2, 後半は酸素原子による反応，Hancoと同じ)
        % 緩和時間の計算
            t_VT = 1.0/(7e-10 * exp(-141*T_gas.^(-1/3)) .* n_air * 1e-6 + 1e-4 .* n_air * 5e-18 .* exp(-128/T_gas.^(0.5)));
            
        % 振動エネルギーの計算
            E_v_1 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_vib))-1); % [J]
            E_v_3 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_gas))-1); % [J]
            
        % エネルギーソース項の計算
            Q_VT = n_air.*(E_v_1-E_v_3)./t_VT; % [J/m^3]
          %}
        %{
    % Vibrational-translational relaxation energy
        % k_2 の計算
            k_VT = 7.8e-18.*T_gas.*exp(-218./(T_gas.^(1/3))+690./T_gas)./(1-exp(-1*q*e_vib_N2./(kB.*T_gas)));
        
        % 振動エネルギーの計算
            E_v_1 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_vib))-1); % [J]
            E_v_3 = q*e_vib_N2./(exp(q*e_vib_N2./(kB*T_gas))-1); % [J]
            
        % エネルギーソース項の計算
            Q_VT = n_air.*n_air.*k_VT.*(1-exp(-1*q*e_vib_N2./(kB.*T_gas))).*(E_v_1-E_v_3); % [J/m^3]
        %}
        
    % V-T 緩和 (PATTERN 3.)
                % Energy calculation
                    E_vib_N2    = n_N2 * q*e_vib_N2./(exp(q*e_vib_N2/(kB*T_vib))-1);
                    E_vib_N2_eq = n_N2 * q*e_vib_N2./(exp(q*e_vib_N2/(kB*T_gas))-1);
                    E_vib_O2    = n_O2 * q*e_vib_O2./(exp(q*e_vib_O2/(kB*T_vib))-1);
                    E_vib_O2_eq = n_O2 * q*e_vib_O2./(exp(q*e_vib_O2/(kB*T_gas))-1);
                
                % Relaxation time
                    t_VT_O2 = sum(phi_O2.*(1/(P_ambient/101300) .* exp(A_O2.*(T_gas.^(-1/3)-0.015*mu_O2.^(1/4))-18.42)).^(-1));
                    t_VT_O2 = (t_VT_O2).^(-1);
                    t_VT_N2 = sum(phi_N2.*(1/(P_ambient/101300) .* exp(A_N2.*(T_gas.^(-1/3)-0.015*mu_N2.^(1/4))-18.42)).^(-1));
                    t_VT_N2 = (t_VT_N2).^(-1);
                    
                % Source term
                    Q_VT_O2 = (E_vib_O2 - E_vib_O2_eq) ./ t_VT_O2;
                    Q_VT_N2 = (E_vib_N2 - E_vib_N2_eq) ./ t_VT_N2;
                    
                    
            % Vibrational energy loss due to strong non-equilibrium // 
                % Reaction rate
                    k_51 = 4.98e-9.*T_gas.^(-1.5).*exp(-113260./T_gas);
                    k_52 = k_51.*ratio_N2/(ratio_N)^2;
                    
                    alpha_D = (m_N/(m_N2+m_N))^2;
                    n_D = -1.5;
                    L_D = 2*(1-alpha_D)/(pi^2*alpha_D^(3/4)) * (T_gas/113200).^(1.5-n_D) .* ( 1 + 7*(1-alpha_D)*(1+(alpha_D)^0.5)*T_gas/(2*113200));
                    T_a = alpha_D * T_vib + (1-alpha_D) * T_gas;
                    Z_D = (1-exp(-3354/T_vib))./(1-exp(-3354/T_gas)) .* (1-L_D) .* exp(-113200*(1/T_vib-1/T_gas))...
                            + L_D .* exp(-113200*(1/T_a-1/T_gas));

                % Source term
                    Q_D = q*e_D*(Z_D.*k_51.*n_N2.^2);
                    for j = 1:Ny             % 陽解法で解く際は必要．発散してしまう
                        for i= 1:Nx
                            if Q_D(1,1,j,i) > 1e+11
                                Q_D(1,1,j,i) = 1e+11;
                            end
                        end
                    end

    
    % 運動量輸送エネルギー / Elastic collision energy transfer
    % とりあえずサハの式から電子数密度を算出する．
            %n_e_d = (n_air*4.82e+21.*T_vib.^1.5.*exp(-11600*epsilon_ion_N2./T_vib)).^0.5;
            Q_ela = 3*m_e/m_air*kB.*5.3e+9.*P*760/101300.*(T_vib-T_gas).*1e+21*P_ambient/101300;            % 全領域で1e+21にしているので，緩和が速い可能性がある．
            
            
    % SUMMARY
        Q_vib = eta_vib/100 .* W - Q_VT_N2 - Q_VT_O2 - Q_ela;% - 2*Q_D;
        Q_gas = eta_gas/100 .* W + Q_VT_N2 + Q_VT_O2 + Q_ela;% + Q_D;
            
        
    %% Solve Euler equation
    % Time evolution
        %if k == 1
            Q(1:3,:,:,2:Nx-1) = Q(1:3,:,:,2:Nx-1) - 0.5*dt/dx * (E_half(1:3,:,:,2:Nx-1) - E_half(1:3,:,:,1:Nx-2));
            Q(4,:,:,2:Nx-1) = Q(4,:,:,2:Nx-1) - 0.5*dt/dx * (E_half(4,:,:,2:Nx-1) - E_half(4,:,:,1:Nx-2)) + Q_gas(1,1,:,2:Nx-1) * 0.5*dt;
            Q(5,:,:,2:Nx-1) = Q(5,:,:,2:Nx-1) - 0.5*dt/dx * (E_half(5,:,:,2:Nx-1) - E_half(5,:,:,1:Nx-2)) + Q_vib(1,1,:,2:Nx-1) * 0.5*dt;
        %elseif k == 3
        %    Q(:,:,:,2:Nx-1) = Q(:,:,:,2:Nx-1) - 0.5*dt/dx * (E_half(:,:,:,2:Nx-1) - E_half(:,:,:,1:Nx-2));
        %end
    
   
    % Boundary conditions (Left): Closed boundary
        T_temp_1(1,1,:) = (gamma-1)/R_air/Q(1,1,:,2).*(Q(4,1,:,2)-0.5*(Q(2,1,:,2).^2+Q(3,1,:,2).^2)./Q(1,1,:,2));
        Q(1,1,:,1) = Q(1,1,:,2);
        Q(2,1,:,1) = Q(2,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(3,1,:,1) = Q(3,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(4,1,:,1) = Q(4,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(5,1,:,1) = Q(5,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
    
    % Boundary conditions (Right): Open boundary
        T_temp_1(1,1,:) = (gamma-1)/R_air/Q(1,1,:,Nx-1).*(Q(4,1,:,Nx-1)-0.5*(Q(2,1,:,Nx-1).^2+Q(3,1,:,Nx-1).^2)./Q(1,1,:,Nx-1));
        Q(1,1,:,Nx) = P_ambient/R_air./T_temp_1;
        Q(2,1,:,Nx) = Q(2,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(3,1,:,Nx) = Q(3,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(4,1,:,Nx) = Q(4,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(5,1,:,Nx) = Q(5,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        
    % Boundary conditions (Upper): Closed boundary
        T_temp_2(1,1,:) = (gamma-1)/R_air/Q(1,1,2,:).*(Q(4,1,2,:)-0.5*(Q(2,1,2,:).^2+Q(3,1,2,:).^2)./Q(1,1,2,:));
        Q(1,1,1,:) = Q(1,1,2,:);
        Q(2,1,1,:) = Q(2,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(3,1,1,:) = Q(3,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(4,1,1,:) = Q(4,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(5,1,1,:) = Q(5,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
    
    % Boundary conditions (Lower): Open boundary
        T_temp_2(1,1,:) = (gamma-1)/R_air/Q(1,1,Ny-1,:).*(Q(4,1,Ny-1,:)-0.5*(Q(2,1,Ny-1,:).^2+Q(3,1,Ny-1,:).^2)./Q(1,1,Ny-1,:));
        Q(1,1,Ny,:) = Q(1,1,Ny-1,:);
        Q(2,1,Ny,:) = Q(2,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(3,1,Ny,:) = Q(3,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(4,1,Ny,:) = Q(4,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(5,1,Ny,:) = Q(5,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        
        

        
    %% Update matrix
        
        rho = Q(1,1,:,:);
        n_air = rho/m_air;
        n_N2 = ratio_N2*n_air;
        n_O2 = ratio_O2*n_air;
        u = Q(2,1,:,:)./rho;
        v = Q(3,1,:,:)./rho;
        P = (gamma-1)*(Q(4,1,:,:)-0.5.*rho.*(u.^2+v.^2));
        H = (Q(4,1,:,:)+P)./rho;
        T_gas = P./(R_air.*rho);
        e_gas = Q(4,1,:,:)./rho - 0.5.*(u.^2+v.^2);
        e_vib = Q(5,1,:,:)./rho;
        
        Temp_E_vib = Q(5,1,:,:);
        
        % Newton法による振動温度の算出
          error = 1e-9;
          kmax = 50;
          for j = 1:Ny
              for i= 1:Nx
                  x = 300;
                  for k = 1:kmax
                      E_v(1,1,j,i) =  n_N2(1,1,j,i) * q*e_vib_N2/(exp(q*e_vib_N2/(kB*x))-1) + n_O2(1,1,j,i) * q*e_vib_O2/(exp(q*e_vib_O2/(kB*x))-1) - Temp_E_vib(1,1,j,i);
                      E_v_d(1,1,j,i) = n_N2(1,1,j,i) * (q*e_vib_N2/kB)/x^(2) * q*e_vib_N2 * exp((q*e_vib_N2/kB)/x) / (exp(q*e_vib_N2/kB/x)-1)^(2) ...
                              + n_O2(1,1,j,i) * (q*e_vib_O2/kB)/x^(2) * q*e_vib_O2 * exp((q*e_vib_O2/kB)/x) / (exp(q*e_vib_O2/kB/x)-1)^(2);    % 微分した関数
                      x = x - E_v(1,1,j,i)/E_v_d(1,1,j,i);
                    
                      if abs(E_v(1,1,j,i)) < error
                          break;
                      end
                  end
                  T_vib(1,1,j,i) = x;
              end
          end
        
        
        E(1,1,:,:) = Q(2,1,:,:);
        E(2,1,:,:) = rho.*u.^2 + P;
        E(3,1,:,:) = rho.*u.*v;
        E(4,1,:,:) = (Q(4,1,:,:)+P).*u;
        E(5,1,:,:) = Q(5,1,:,:).*u;
        
        F(1,1,:,:) = Q(3,1,:,:);
        F(2,1,:,:) = rho.*u.*v;
        F(3,1,:,:) = rho.*v.^2 + P;
        F(4,1,:,:) = (Q(4,1,:,:)+P).*v;
        F(5,1,:,:) = Q(5,1,:,:).*v;
   
        
   
   elseif k == 2        % y方向の時間発展
       
   %% Calculation of values between cells
        
        rho_av(1,1,1:Ny-1,:) = ( rho(1,1,1:Ny-1,:) .* rho(1,1,2:Ny,:) ).^0.5;
        u_av(1,1,1:Ny-1,:) = (rho(1,1,1:Ny-1,:).^0.5.*u(1,1,1:Ny-1,:) + rho(1,1,2:Ny,:).^0.5.*u(1,1,2:Ny,:))./(rho(1,1,1:Ny-1,:).^0.5+rho(1,1,2:Ny,:).^0.5);
        v_av(1,1,1:Ny-1,:) = (rho(1,1,1:Ny-1,:).^0.5.*v(1,1,1:Ny-1,:) + rho(1,1,2:Ny,:).^0.5.*v(1,1,2:Ny,:))./(rho(1,1,1:Ny-1,:).^0.5+rho(1,1,2:Ny,:).^0.5);
        H_av(1,1,1:Ny-1,:) = (rho(1,1,1:Ny-1,:).^0.5.*H(1,1,1:Ny-1,:) + rho(1,1,2:Ny,:).^0.5.*H(1,1,2:Ny,:))./(rho(1,1,1:Ny-1,:).^0.5+rho(1,1,2:Ny,:).^0.5);
        c_av(1,1,1:Ny-1,:) = ((gamma-1)*(H_av(1,1,1:Ny-1,:)-0.5*(u_av(1,1,1:Ny-1,:).^2+v_av(1,1,1:Ny-1,:).^2))).^0.5;
    
        rho_av(1,1,Ny,:) = rho_av(1,1,Ny-1,:);
        u_av(1,1,Ny,:) = u_av(1,1,Ny-1,:);
        v_av(1,1,Ny,:) = v_av(1,1,Ny-1,:);
        H_av(1,1,Ny,:) = H_av(1,1,Ny-1,:);
        c_av(1,1,Ny,:) = c_av(1,1,Ny-1,:);
    
    
   %% Calculation of matrix, R
    
        R(1,1,:,:) = 1; 
        R(2,1,:,:) = u_av;
        R(3,1,:,:) = v_av-c_av;
        R(4,1,:,:) = H_av-v_av.*c_av;
        
        R(1,2,:,:) = 1;
        R(2,2,:,:) = u_av;
        R(3,2,:,:) = v_av;
        R(4,2,:,:) = (u_av.^2+v_av.^2)/2;
        
        R(1,3,:,:) = 1;
        R(2,3,:,:) = u_av;
        R(3,3,:,:) = v_av+c_av;
        R(4,3,:,:) = H_av+v_av.*c_av;
        
        R(1,4,:,:) = 0;
        R(2,4,:,:) = 1;
        R(3,4,:,:) = 0;
        R(4,4,:,:) = u_av;

        
    %% Calculation of inverse matrix, R^(-1)
        
        b2 = (gamma-1)./c_av.^2;
        b1 = b2.*(u_av.^2+v_av.^2)/2;
    
        R_inv(1,1,:,:) = (b1+v_av./c_av)/2;
        R_inv(2,1,:,:) = 1-b1;
        R_inv(3,1,:,:) = (b1-v_av./c_av)/2;
        R_inv(4,1,:,:) = -u_av;
        
        R_inv(1,2,:,:) = -b2.*u_av/2;
        R_inv(2,2,:,:) = b2.*u_av;
        R_inv(3,2,:,:) = -b2.*u_av/2;
        R_inv(4,2,:,:) = 1;
        
        R_inv(1,3,:,:) = (-b2.*v_av-1./c_av)/2;
        R_inv(2,3,:,:) = b2.*v_av;
        R_inv(3,3,:,:) = (-b2.*v_av+1./c_av)/2;
        R_inv(4,3,:,:) = 0;
        
        R_inv(1,4,:,:) = b2/2;
        R_inv(2,4,:,:) = -b2;
        R_inv(3,4,:,:) = b2/2;
        R_inv(4,4,:,:) = 0;
    
    
   %% Calculation of diagonal matrix
        
        L(1,1,:,:) = v_av - c_av;
        L(2,2,:,:) = v_av;
        L(3,3,:,:) = v_av + c_av;
        L(4,4,:,:) = v_av;

    
   %% Calculation of numerical flux between cells
  
    for i = 1:Nx
        for j = 1:Ny-1
            F_half(1:4,:,j,i) = 0.5 * ( F(1:4,:,j+1,i)+F(1:4,:,j,i) - R(1:4,1:4,j,i)*abs(L(1:4,1:4,j,i))*R_inv(1:4,1:4,j,i)*(Q(1:4,:,j+1,i)-Q(1:4,:,j,i)));
            F_half(5,:,j,i) = 0.5 * ( F(5,:,j+1,i)+F(5,:,j,i) - abs(v_av(1,1,j,i)).*(Q(5,:,j+1,i)-Q(5,:,j,i)));
        end
        F_half(:,:,Ny,i) = F_half(:,:,Ny-1,i);
    end
  
    
   %% Solve Euler equation
    % Time evolution       
        Q(1:3,:,2:Ny-1,:) = Q(1:3,:,2:Ny-1,:) - dt/dy * (F_half(1:3,:,2:Ny-1,:) - F_half(1:3,:,1:Ny-2,:));
        Q(4,:,2:Ny-1,:) = Q(4,:,2:Ny-1,:) - dt/dy * (F_half(4,:,2:Ny-1,:) - F_half(4,:,1:Ny-2,:));%+ Q_gas(1,1,2:Ny-1,:) * 0.5*dt;
        Q(5,:,2:Ny-1,:) = Q(5,:,2:Ny-1,:) - dt/dy * (F_half(5,:,2:Ny-1,:) - F_half(5,:,1:Ny-2,:));%+ Q_electron(1,1,2:Ny-1,:) * 0.5*dt;
        
    % Boundary conditions (Left): Closed boundary
        T_temp_1(1,1,:) = (gamma-1)/R_air/Q(1,1,:,2).*(Q(4,1,:,2)-0.5*(Q(2,1,:,2).^2+Q(3,1,:,2).^2)./Q(1,1,:,2));
        Q(1,1,:,1) = Q(1,1,:,2);
        Q(2,1,:,1) = Q(2,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(3,1,:,1) = Q(3,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(4,1,:,1) = Q(4,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
        Q(5,1,:,1) = Q(5,1,:,2)./Q(1,1,:,2).*Q(1,1,:,2);
    
    % Boundary conditions (Right): Open boundary
        T_temp_1(1,1,:) = (gamma-1)/R_air/Q(1,1,:,Nx-1).*(Q(4,1,:,Nx-1)-0.5*(Q(2,1,:,Nx-1).^2+Q(3,1,:,Nx-1).^2)./Q(1,1,:,Nx-1));
        Q(1,1,:,Nx) = P_ambient/R_air./T_temp_1;
        Q(2,1,:,Nx) = Q(2,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(3,1,:,Nx) = Q(3,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(4,1,:,Nx) = Q(4,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        Q(5,1,:,Nx) = Q(5,1,:,Nx-1)./Q(1,1,:,Nx-1).*Q(1,1,:,Nx);
        
    % Boundary conditions (Upper): Closed boundary
        T_temp_2(1,1,:) = (gamma-1)/R_air/Q(1,1,2,:).*(Q(4,1,2,:)-0.5*(Q(2,1,2,:).^2+Q(3,1,2,:).^2)./Q(1,1,2,:));
        Q(1,1,1,:) = Q(1,1,2,:);
        Q(2,1,1,:) = Q(2,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(3,1,1,:) = Q(3,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(4,1,1,:) = Q(4,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
        Q(5,1,1,:) = Q(5,1,2,:)./Q(1,1,2,:).*Q(1,1,2,:);
    
    % Boundary conditions (Lower): Open boundary
        T_temp_2(1,1,:) = (gamma-1)/R_air/Q(1,1,Ny-1,:).*(Q(4,1,Ny-1,:)-0.5*(Q(2,1,Ny-1,:).^2+Q(3,1,Ny-1,:).^2)./Q(1,1,Ny-1,:));
        Q(1,1,Ny,:) = Q(1,1,Ny-1,:);
        Q(2,1,Ny,:) = Q(2,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(3,1,Ny,:) = Q(3,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(4,1,Ny,:) = Q(4,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        Q(5,1,Ny,:) = Q(5,1,Ny-1,:)./Q(1,1,Ny-1,:).*Q(1,1,Ny,:);
        
        
    %% Update matrix
        
        rho = Q(1,1,:,:);
        n_air = rho/m_air;
        n_N2 = ratio_N2*n_air;
        n_O2 = ratio_O2*n_air;
        u = Q(2,1,:,:)./rho;
        v = Q(3,1,:,:)./rho;
        P = (gamma-1)*(Q(4,1,:,:)-0.5.*rho.*(u.^2+v.^2));
        H = (Q(4,1,:,:)+P)./rho;
        T_gas = P./(R_air.*rho);
        e_gas = Q(4,1,:,:)./rho - 0.5.*(u.^2+v.^2);
        e_vib = Q(5,1,:,:)./rho;
        
        E_vib = Q(5,1,:,:);
        
        % Newton法による振動温度の算出
          error = 1e-9;
          kmax = 50;
          for j = 1:Ny
              for i= 1:Nx
                  x = 300;
                  for k = 1:kmax
                      E_v(1,1,j,i) =  n_N2(1,1,j,i) * q*e_vib_N2/(exp(q*e_vib_N2/(kB*x))-1) + n_O2(1,1,j,i) * q*e_vib_O2/(exp(q*e_vib_O2/(kB*x))-1) - E_vib(1,1,j,i);
                      E_v_d(1,1,j,i) = n_N2(1,1,j,i) * (q*e_vib_N2/kB)/x^(2) * q*e_vib_N2 * exp((q*e_vib_N2/kB)/x) / (exp(q*e_vib_N2/kB/x)-1)^(2) ...
                              + n_O2(1,1,j,i) * (q*e_vib_O2/kB)/x^(2) * q*e_vib_O2 * exp((q*e_vib_O2/kB)/x) / (exp(q*e_vib_O2/kB/x)-1)^(2);    % 微分した関数
                      x = x - E_v(1,1,j,i)/E_v_d(1,1,j,i);
                    
                      if abs(E_v(1,1,j,i)) < error
                          break;
                      end
                  end
                  T_vib(1,1,j,i) = x;
              end
          end


   
   end
        
   end  %% the end of Fractuional Time Step
   
   
   %% Update: Time step / 時間刻み幅の更新
        t = t + dt;
        
        
   %% Output for Graphix
    %% L2ノルムを取り出す
    if mod(step,20) == 0
    L2_1 = 0;
    L2_2 = 0;
    L2_3 = 0;
    L2_4 = 0;
    L2_5 = 0;
    L2_6 = 0;
    for i = 1:Nx
        for i = 1:Ny
            L2_1 = L2_1+(rho_plot(j,i)-rho(1,1,j,i))^2;
            L2_2 = L2_2+(T_gas_plot(j,i)-T_gas(1,1,j,i))^2;
            L2_3 = L2_3+(T_vib_plot(j,i)-T_vib(1,1,j,i))^2;
            L2_4 = L2_4+(P_plot(j,i)-P(1,1,j,i))^2;
            L2_5 = L2_5+(u_plot(j,i)-u(1,1,j,i))^2;
            L2_6 = L2_6+(v_plot(j,i)-v(1,1,j,i))^2;
        end
    end
    L2_1 = L2_1^0.5;
    L2_2 = L2_2^0.5;
    L2_3 = L2_3^0.5;
    L2_4 = L2_4^0.5;
    L2_5 = L2_5^0.5;
    L2_6 = L2_6^0.5;
    
         Norm = [Norm; step L2_1 L2_2 L2_3 L2_4 L2_5 L2_6];
    end
    
    for j = 1:Ny
        for i = 1:Nx
            rho_plot(j,i) = rho(1,1,j,i);
            T_gas_plot(j,i) = T_gas(1,1,j,i);
            T_vib_plot(j,i) = T_vib(1,1,j,i);
            P_plot(j,i) = P(1,1,j,i);
            u_plot(j,i) = u(1,1,j,i);
            v_plot(j,i) = v(1,1,j,i);
            W_plot(j,i) = W(1,1,j,i);
            gamma_plot(j,i) = (rho_plot(j,i)*e_gas(1,1,j,i)+rho_plot(j,i)*e_vib(1,1,j,i)+P(1,1,j,i))/(rho_plot(j,i)*e_gas(1,1,j,i)+rho_plot(j,i)*e_vib(1,1,j,i));
            a_plot(j,i) = (gamma_plot(j,i)*P_plot(j,i)/rho_plot(j,i))^0.5;
            M_plot(j,i) = (u_plot(j,i)^2+v_plot(j,i)^2)^0.5/a_plot(j,i);
            Q_VT_plot(j,i) = Q_VT_O2(1,1,j,i)+Q_VT_N2(1,1,j,i);
            Q_D_plot(j,i) = Q_D(1,1,j,i);
        end
    end
    
    
          
    % 衝撃波部分の取り出し（反射衝撃波のそり落とし）
    for i = Nx:-1:1
        if abs(P_plot(Ny,i)-P_ambient)/P_ambient > 2
            Nx_min = i;
            break
        end
    end
    
    % それぞれのyで圧力が最大になるxを抽出する(このやり方だと衝撃波の先端をとらえることができない)
    for j = 1:Ny
        for i = Nx_min:Nx
            if P_plot(j,i) > N_max
                N_max = P_plot(j,i);
                x_p(j) = i*dx;
            end
        end
        N_max = 0;
    end 
    %x_p = smoothdata(x_p);   % dataの平滑化
    
    % それぞれのyで振動温度が最大になるxを抽出する
    N_max = 0;
    for j = 1:Ny
        for i = Nx_min:Nx
            if T_vib_plot(j,i) > N_max
                N_max = T_vib_plot(j,i);
                x_v(j) = i*dx;
            end
        end
        N_max = 0;
    end 
    %x_v = smoothdata(x_v);   % dataの平滑化
    
    % 圧力の最大値と振動温度の最大値の位置が一致するとき，前方の衝撃波形成領域のx位置を抜き出す．
    if x_v(Ny/2) == x_p(Ny/2)
        for j = 1:Ny
            for i = Nx:-1:2
                if P_plot(j,i-1) > P_plot(j,i)*1.5
                    x_p(j) = (i-1)*dx;
                    break;
                end
            end
        end
    end
    
    d_vp = x_v-x_p;
    if mod(step,20) == 0
         Distance = [Distance; step d_vp];
    end
    
    
     P_closed(:) = P_plot(50,:);
     T_gas_closed(:) = T_gas_plot(50,:);
     T_vib_closed(:) = T_vib_plot(50,:);
     rho_closed(:) = rho_plot(50,:);
     u_closed(:) = u_plot(50,:);
     v_closed(:) = v_plot(50,:);
     
     if mod(t,1e-6) < dt
         P_wall = [P_wall; (t-t_initial)*1e+6 mean(P(1,1,:,1))];
     end
     
     %xlswrite(filename_excel_1,P_wall);
     
    
    
     
    %% Output for images
    
    if mod(t, 1e-5) < dt  %mod(step,10) == 0
        fig = figure('visible', 'off');
        
        figure_pressure_exp3(X,Y,P_plot);
        filename = [sprintf('Sim_data/Images/Pressure/Pressure_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_density_exp3(X,Y,rho_plot);
        filename = [sprintf('Sim_data/Images/Density/Density_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_temperature_exp3(X,Y,T_gas_plot);
        filename = [sprintf('Sim_data/Images/Temperature/Temperature_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_vib_temperature_exp3(X,Y,T_vib_plot);
        filename = [sprintf('Sim_data/Images/Vib_Temperature/Vib_Temperature_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_input_power_exp3(X,Y,W_plot);
        filename = [sprintf('Sim_data/Images/Power/Power_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        
        figure_HSHS_exp3(X,Y,P_plot,T_vib_plot);
        filename = [sprintf('Sim_data/Images/HSHS/HSHS_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        
        figure_u_exp3(X,Y,u_plot);
        filename = [sprintf('Sim_data/Images/Velocity_u/Velocity_u_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_v_exp3(X,Y,v_plot);
        filename = [sprintf('Sim_data/Images/Velocity_v/Velocity_v_%d.jpg', ii)];        
        SaveImage(fig,500,400,filename);
        
        figure_Mach_exp3(X,Y,M_plot);
        filename = [sprintf('Sim_data/Images/Mach/Mach_%d.jpg', ii)];        
        SaveImage(fig,600,400,filename);
 
        %{
        Pressure_plasma(X,P_closed);
        filename = [sprintf('Sim_data/Images/1D/Pressure/Pressure_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
        
        Temperature_plasma(X,T_gas_closed);
        filename = [sprintf('Sim_data/Images/1D/Temperature/Temperature_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
        
        Vib_Temperature_plasma(X,T_vib_closed);
        filename = [sprintf('Sim_data/Images/1D/Vib_Temperature/Vib_Temperature_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
        
        Density_plasma(X,rho_closed);
        filename = [sprintf('Sim_data/Images/1D/Density/Density_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
        
        U_plasma(X,u_closed);
        filename = [sprintf('Sim_data/Images/1D/Velocity_u/Velocity_u_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
        
        V_plasma(X,rho_closed);
        filename = [sprintf('Sim_data/Images/1D/Velocity_v/Velocity_v_%d.jpg', ii)];        
        SaveImage(fig,600,450,filename);
   %}
        
        ii = ii+1;
    end
    
    %%　圧力のピーク位置と電離波面の位置関係を横軸時間でプロットする．
    
     
end

    xlswrite(filename_excel_2,Distance);
    xlswrite(filename_excel_3,Norm);