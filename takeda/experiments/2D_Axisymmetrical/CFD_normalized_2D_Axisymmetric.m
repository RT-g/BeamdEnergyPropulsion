% 2次元軸対称流れCFD
% 一般化座標
% auther = Ryota Takeda

clear, close all
% clearvars

% Movie
base_dir = fullfile(pwd, string(datetime(now, 'ConvertFrom', 'datenum', 'Format', 'yyMMddHH')));
mkdir(base_dir, 'Pressure');
mkdir(base_dir, 'Temperature');
mkdir(base_dir, 'U');
mkdir(base_dir, 'V');

tic
vp=VideoWriter(fullfile(base_dir,'movie_Pressure.avi'));
vt=VideoWriter(fullfile(base_dir,'movie_Temperature.avi'));
vu=VideoWriter(fullfile(base_dir,'movie_Velocity_u.avi'));
vv=VideoWriter(fullfile(base_dir,'movie_Velocity_v.avi'));
open(vp);
open(vt);
open(vu);
open(vv);

% Computation Values
dx = 0.01e-3; % m Segment length, xi(ξ) direction
dh = 0.01e-3; % m Segment length, eta(η) direction

nx = 250; % Number of position segments, xは進展方向
nh = 180; %2.5e-3/dh; %全長2.5 mmになるように設定 
nt = 2000; % Number of time segments
ng = 40;%40; % Graphの分割数を決める値
CFL = 0.7; % クーラン数 1以下にする 0.5~0.99くらいで，小さくしすぎると進みが遅い
select_gas = 'air';
theta_1 = 1;
theta_2 = 0.5;
phi = 1/3;
beta_comp = 2; %圧縮パラメータ, (1,4]

% Physical Values

% Laser Values
A_G = 0.224517656;
B_G = 0.77548;
sigma_G1 = 0.84473; %mm
sigma_G2 = 1.75302; %mm
A_T = 0.764525962;
B_T = 0.235474467;
sigma_T1 = 1.557528296;
sigma_T2 = 4.050355336;
lambda = 10.6; %単位はum
M2_G = 15;
M2_T = 21;
W_G0 = 1.7;%単位はmm
W_T0 = 2.0;
R_peak = 2.13;
a_s = (B_G*1e6/sigma_G2^2 + B_T*1e6/sigma_T2^2);

% Values specific to gas
filename = '../../data/value.csv';
opts = detectImportOptions(filename);
gas_table = readtable(filename, 'ReadRowNames', true);
m = gas_table(select_gas,'m').Variables;
gamma = gas_table(select_gas, 'gamma').Variables; % 本当はガンマ関数があるのでこの変数名はよくない
eta_trans = gas_table(select_gas, 'eta_trans').Variables;
slope = gas_table(select_gas, 'a').Variables;
intercept = gas_table(select_gas, 'b').Variables;
slope_low = gas_table(select_gas, 'a_low').Variables;
intercept_low = gas_table(select_gas, 'b_low').Variables;
a = gas_table(select_gas, 'speed_of_sound').Variables;
R_0 = 8314.0;
R = R_0/m;

% Common Values for gases
u_ionz0 = slope*524.4^intercept*1.2*1e3; %m/s
P_0 = 1.013e5*1;
T_0 = 293;
rho_0 = P_0/R/T_0;
u_0 = 0;
v_0 = 0;

l = 0.2e-3; %m 加熱長さ レーザーは0.2 mm
hr = 20; %Heating region

% Initial Values

P = zeros(nx,nh);
Pplateau = zeros(nt,nh);
u = zeros(nx,nh);
v = zeros(nx,nh);
T = zeros(nx,nh);
rho = zeros(nx,nh);
H = zeros(nx,nh);
t_list = zeros(nt,1);
x_list = (1:nx) * dx * 1e3; % mm
r_list = (1:nh) * dh * 1e3; % mm
uionz_list = zeros(nt,1);
xionz = zeros(nh,1);
t_txt = zeros(nt,1);
S_laser = zeros(1,nh);

P(:,:) = P_0;
Pplateau(:,:) = P_0;
u(:,:) = u_0;
v(:,:) = v_0;
rho(:,:) = P_0/R/T_0;
T(:,:) = T_0;

% Vector of states
Q1 = ones(nx,nh).*rho; % Mass
Q2 = ones(nx,nh).*rho.*u; % Momentum x
Q3 = ones(nx,nh).*rho.*v; %Momentum r
Q4 = ones(nx,nh).*P/(gamma-1)+1/2.*rho.*(u.^2+v.^2); % Energy

% New calculated vector of states
Q1_cal = ones(nx,nh).*rho; % Mass
Q2_cal = zeros(nx,nh);
Q3_cal = zeros(nx,nh);
Q4_cal = ones(nx,nh).*P/(gamma-1);

% diff vector of states
dQ_1 = zeros(nx,nh); 
dQ_2 = zeros(nx,nh);
dQ_3 = zeros(nx,nh);
dQ_4 = zeros(nx,nh);

rhs1 = zeros(nx,nh);
rhs2 = zeros(nx,nh);
rhs3 = zeros(nx,nh);
rhs4 = zeros(nx,nh);

dQ1_1 = zeros(nx,nh); 
dQ1_2 = zeros(nx,nh);
dQ1_3 = zeros(nx,nh);
dQ1_4 = zeros(nx,nh);

dQ2_1 = zeros(nx,nh); 
dQ2_2 = zeros(nx,nh);
dQ2_3 = zeros(nx,nh);
dQ2_4 = zeros(nx,nh);

dQ3_1 = zeros(nx,nh); 
dQ3_2 = zeros(nx,nh);
dQ3_3 = zeros(nx,nh);
dQ3_4 = zeros(nx,nh);

dQ4_1 = zeros(nx,nh); 
dQ4_2 = zeros(nx,nh);
dQ4_3 = zeros(nx,nh);
dQ4_4 = zeros(nx,nh);

dQ5_1 = zeros(nx,nh); 
dQ5_2 = zeros(nx,nh);
dQ5_3 = zeros(nx,nh);
dQ5_4 = zeros(nx,nh);

% Vector of flux
E1 = zeros(nx,nh); % 0
E2 = P; % P
E3 = zeros(nx,nh); % 0
E4 = zeros(nx,nh); % 0

F1 = zeros(nx,nh); %0
F2 = zeros(nx,nh); %0
F3 = P; %P
F4 = zeros(nx,nh); %0

H1 = zeros(nx,nh); %0
H2 = zeros(nx,nh); %0
H3 = zeros(nx,nh); %0
H4 = zeros(nx,nh); %0

W = zeros(nx,nh);

% Vector of flux between cells
E1_half = zeros(nx,nh);
E2_half = ones(nx,nh).*rho.*u.^2+P; % P
E3_half = zeros(nx,nh);
E4_half = zeros(nx,nh);

F1_half = zeros(nx,nh);
F2_half = zeros(nx,nh);
F3_half = ones(nx,nh).*rho.*v.^2+P; %P
F4_half = zeros(nx,nh);

% Vector of coordinate transformation
x_x = zeros(nx,nh); % ξ_x
x_r = zeros(nx,nh); % ξ_r
h_x = zeros(nx,nh); % η_x
h_r = zeros(nx,nh); % η_r
r = zeros(1,nh); % r, catersian coodinate
% Vector of Jacobian
J = ones(nx,nh);

% ADI
T_x = zeros(4,4);
T_h = zeros(4,4);

% Calculation
%eta = eta_trans
eta = 0.05; % eta_trans ではどうやら値が大きすぎる
t = 0;
dt = dx/u_ionz0; %s すぐに淘汰される値なので意味はない。初期値
x_laser0 = 0;
umax = u_ionz0;
vmax = u_ionz0;
I = 0;

for n1 = 1:nt
    % 時間刻み幅の設定。ここ次第で計算がいくらでも長くなる
    % レーザー強度の関数の形状に応じてクーラン数を調整
    if (t*1e6<0.085)
        CFL=0.7;
    elseif (t*1e6<0.125)
        CFL=0.8;
    else
        CFL=0.9; 
    end
    
    dt = CFL * min(dx/umax,  dh/vmax); %s, 計算位置>波面位置となるようにdtを決定
    
    t = t+dt; %s
    t_list(n1,1) = t * 1e6; %us

    % レーザー強度の時間減衰(10J) 単位はMW
    Power_laser = 1; 
    % Power_laser = 8.15*exp(-0.866*t*1e6); 
    % 正確に再現すると温度などが高く出すぎる
%     if (t*1e6<0.085)
%         Power_laser = 287.9222823*t*1e6+0.0005175756469;
%     elseif (t*1e6<0.125)
%         Power_laser =-476.3277981*t*1e6+66.1043297;
%     else
%         Power_laser = 8.15*exp(-0.866*t*1e6); 
%     end  
    
    % 進展方向のレーザー強度の変化を考慮しない場合のビーム半径
    W_G = W_G0*1e-3; %m
    W_T = W_T0*1e-3; %m

    % 波頭の値
    S_laser0 = R_peak * Power_laser / 4 / W_G / W_T * 1e-3; % GW/m^2 波頭のレーザー強度
    u_ionz_line1 = slope * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept * a; %m/s 波頭の速度
    u_ionz_line3 = slope_low * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept_low * a;
    if (u_ionz_line1 > u_ionz_line3)
        u_ionz0 = u_ionz_line1; % Line1, 松井さん博論
    else
        u_ionz0 = u_ionz_line3; % Line3, 松井さん博論
    end
    x_laser0 = x_laser0 + u_ionz0 * dt; %m Ionized wave front

    % η方向の計算
    for n3 = 1:nh        
        h = n3*dh; %m, general coordinate

        % 速度の調整
        % ガウシアン分布/トップハット分布を仮定 二次元極座標のためx,yをそれぞれr/sqrt(2)としている。
        % line1
        b_s_line1 = intercept/(1-intercept);
        r_line1 = -(2/(3*a_s^2*b_s_line1^2*h + sqrt(9*a_s^4*b_s_line1^4*h^2 + 4*a_s^3*b_s_line1^3)))^(1/3) + ((3*a_s^2*b_s_line1^2*h + sqrt(9*a_s^4*b_s_line1^4*h^2 + 4*a_s^3*b_s_line1^3))/2)^(1/3)/a_s/b_s_line1; %m, Cartesian coordinates
        G_r_line1 = A_G*exp(-2*(r_line1*1e3/sqrt(2)/sigma_G1)^4)+B_G*exp(-2*(r_line1*1e3/sqrt(2)/sigma_G2)^2);
        T_r_line1 = A_T*exp(-2*(r_line1*1e3/sqrt(2)/sigma_T1)^4)+B_T*exp(-2*(r_line1*1e3/sqrt(2)/sigma_T2)^2);
        cos_line1 = (G_r_line1*T_r_line1)^b_s_line1;
        u_ionz_line1 = slope * (S_laser0 * cos_line1 * 1e9 / rho_0 / a ^ 3) ^ intercept * a; %m/s 波頭の速度
        % line3
        b_s_line3 = intercept_low/(1-intercept_low);
        r_line3 = -(2/(3*a_s^2*b_s_line3^2*h + sqrt(9*a_s^4*b_s_line3^4*h^2 + 4*a_s^3*b_s_line3^3)))^(1/3) + ((3*a_s^2*b_s_line3^2*h + sqrt(9*a_s^4*b_s_line3^4*h^2 + 4*a_s^3*b_s_line3^3))/2)^(1/3)/a_s/b_s_line3; %m, Cartesian coordinates
        G_r_line3 = A_G*exp(-2*(r_line3*1e3/sqrt(2)/sigma_G1)^4)+B_G*exp(-2*(r_line3*1e3/sqrt(2)/sigma_G2)^2);
        T_r_line3 = A_T*exp(-2*(r_line3*1e3/sqrt(2)/sigma_T1)^4)+B_T*exp(-2*(r_line3*1e3/sqrt(2)/sigma_T2)^2);
        cos_line3 = (G_r_line3*T_r_line3)^b_s_line3;
        u_ionz_line3 = slope_low * (S_laser0 * cos_line3 * 1e9 / rho_0 / a ^ 3) ^ intercept_low * a;
        % 速度が大きい方を取得
        if  (u_ionz_line1 > u_ionz_line3)
            u_ionz = u_ionz_line1; % Line1, 松井さん博論
            b = intercept;
            b_s = b_s_line1;
            r(n3) = r_line1;
            S_laser(n3) = S_laser0 * cos_line1;
        else
            u_ionz = u_ionz_line3; % Line3, 松井さん博論
            b = intercept_low;
            b_s = b_s_line3;
            r(n3) = r_line3;
            S_laser(n3) = S_laser0 * cos_line3;
        end
        % Jacobian
        J(:,n3) = (1-a_s*r(n3)^2)^(-b_s);
        x_x(:,n3) = 1;
        x_r(:,n3) = 2*a_s*b_s*r(n3) * (1-a_s*r(n3)^2)^(-2*b_s-1) * ((1-a_s*r(n3)^2)^(-2*b_s)-1)^(-0.5);
        h_x(:,n3) = 0;
        h_r(:,n3) = (1-a_s*r(n3)^2)^(-b_s);

        % ξ方向の計算
        for n2 = 1:nx
            x = n2*dx; %m, general coordinate
            % 加熱領域をx_laser0よりも前にしてしまうと計算が壊れるので実際に考えられる値よりも少し進めたほうがいい?
            if x > x_laser0-l && x < x_laser0 % ξが現実の長さを表しているわけではないので実際は少し補正が必要だが、ここでは無視
                % W(n2,n3) = eta * S_laser(n3) / l * 1e9; %W/m3
                W(n2,n3) = eta * S_laser0 / l * 1e9; %W/m3
            else
                W(n2,n3) = 0;
            end
        end
    end
    % Beam Warming Calculation ビームウォーミング法
    Q1_cal = Q1 + dQ_1;
    Q2_cal = Q2 + dQ_2;
    Q3_cal = Q3 + dQ_3;
    Q4_cal = Q4 + dQ_4;

    % Boundary Conditions    
    % Left(回転軸)
    for i = 2:nx-1
        Q1_cal(i,1) = Q1_cal(i,2); 
        Q2_cal(i,1) = Q2_cal(i,2); 
        Q3_cal(i,1) = 0; %中心ではur=0のはず
        Q4_cal(i,1) = Q4_cal(i,2);
    end
    
    % Right(流出)
    for i = 2:nx-1
        Q1_cal(i,nh) = Q1_cal(i,nh-1);
        Q2_cal(i,nh) = Q2_cal(i,nh-1);
        Q3_cal(i,nh) = Q3_cal(i,nh-1);
        Q4_cal(i,nh) = Q4_cal(i,nh-1);
    end
    
    % Top(光軸方向)流出
    for i = 2:nh-1
        Q1_cal(nx,i) = Q1_cal(nx-1,i);
        Q2_cal(nx,i) = Q2_cal(nx-1,i);
        Q3_cal(nx,i) = Q3_cal(nx-1,i);
        Q4_cal(nx,i) = Q4_cal(nx-1,i);
    end
    
    % Plate(アルミ板)
    Q1_cal(1,:) = Q1_cal(2,:);
    Q2_cal(1,:) = 0;
    Q3_cal(1,:) = 0;
    Q4_cal(1,:) = Q4_cal(2,:);
    
    % Corners
    Q1_cal(nx,1) = (Q1_cal(nx,2) + Q1_cal(nx-1,1))/2;
    Q2_cal(nx,1) = (Q2_cal(nx,2) + Q2_cal(nx-1,1))/2;
    Q3_cal(nx,1) = (Q3_cal(nx,2) + Q3_cal(nx-1,1))/2;
    Q4_cal(nx,1) = (Q4_cal(nx,2) + Q4_cal(nx-1,1))/2;
    
    Q1_cal(nx,nh) = (Q1_cal(nx-1,nh) + Q1_cal(nx,nh-1))/2;
    Q2_cal(nx,nh) = (Q2_cal(nx-1,nh) + Q2_cal(nx,nh-1))/2;
    Q3_cal(nx,nh) = (Q3_cal(nx-1,nh) + Q3_cal(nx,nh-1))/2;
    Q4_cal(nx,nh) = (Q4_cal(nx-1,nh) + Q4_cal(nx,nh-1))/2;
    

    % 各行列の更新
    
    Q1 = Q1_cal;
    Q2 = Q2_cal;
    Q3 = Q3_cal;
    Q4 = Q4_cal;

    rho = J.*Q1;
    u = J.*Q2./rho; %ξ方向速度
    v = J.*Q3./rho; %η方向速度
    e = J.*Q4;
    P = (gamma-1)*(J.*Q4-(1/2)*rho.*(u.^2+v.^2));
    Pplateau(n1,:) = P(1,:);
    % disp(Pplateau)
    H = (e+P)./rho;
    T = P./rho./R; % Temperature
    U = x_x.*u + x_r.*v;
    V = h_x.*u + h_r.*v;

    E1 = rho.*U./J;
    E2 = (rho.*u.*U + x_x.*P)./J;
    E3 = (rho.*U.*v + x_r.*P)./J;
    E4 = (e+P).*U./J;
    
    F1 = rho.*V./J;
    F2 = (rho.*V.*u + h_x.*P)./J;
    F3 = (rho.*v.*V + h_r.*P)./J;
    F4 = (e+P).*V./J;
    
    H1 = rho.*v./J;
    H2 = rho.*u.*v./J;
    H3 = rho.*v.^2./J;
    H4 = (e+P).*v./J;

    umax = max(u, [], 'all');
    vmax = max(v, [], 'all');

    for n3 = 2:nh-2
        for n2 = 2:nx-2
            % ヤコビアンの逆
            x_xi = h_r(n2,n3)/J(n2,n3);
            x_h = -x_r(n2,n3)/J(n2,n3);
            r_xi = -h_x(n2,n3)/J(n2,n3);
            r_h = x_x(n2,n3)/J(n2,n3);
            g11 = x_xi^2 + r_xi^2;
            g22 = x_h^2 + r_h^2;
            g12 = x_xi*x_h + r_xi*r_h;
            G0 = sqrt(g11/g22);
            
            % ξ方向の計算
            % セル間の平均値の計算
            % Expressing rho, u, v and H in the basis of w (Roe's paper)
            D_x = sqrt(rho(n2+1,n3)/rho(n2,n3));
            u_av_x = (u(n2,n3)+D_x*u(n2+1,n3))/(1+D_x); 
            v_av_x = (v(n2,n3)+D_x*v(n2+1,n3))/(1+D_x); 
            H_av_x = (H(n2,n3)+D_x*H(n2+1,n3))/(1+D_x); 
            a_av_x = sqrt((gamma-1)*(H_av_x-1/2*(u_av_x.^2 + v_av_x.^2))); %Square root of the derivative of P br rho, used in the eigenvalues of the Jacobian matrix

            % ヤコビアン行列の固有値計算 
            U_av_x = u_av_x*x_x(n2,n3) + v_av_x*x_r(n2,n3);
            lambda_x(1) = U_av_x;
            lambda_x(2) = U_av_x;
            lambda_x(3) = U_av_x + a_av_x*sqrt(x_x(n2,n3)^2+x_r(n2,n3)^2);
            lambda_x(4) = U_av_x - a_av_x*sqrt(x_x(n2,n3)^2+x_r(n2,n3)^2);

            Lambda_x_p = diag((lambda_x + abs(lambda_x))/2);
            Lambda_x_m = diag((lambda_x - abs(lambda_x))/2);

            N = zeros(4,4);
            N(1,1) = 1;
            N(2,1) = -u_av_x;
            N(2,2) = 1;
            N(3,1) = -v_av_x;
            N(3,3) = 1;
            N(4,1) = (gamma-1)*(u_av_x.^2 + v_av_x.^2)/2;
            N(4,2) = -(gamma-1)*u_av_x;
            N(4,3) = -(gamma-1)*v_av_x;
            N(4,4) = gamma-1;

            C_x = zeros(4,4);
            C_x(1,1) = 1;
            C_x(1,4) = -1/a_av_x^2;
            C_x(2,2) = r_h;
            C_x(2,3) = -x_h;
            C_x(2,4) = sqrt(g22)/a_av_x;
            C_x(3,2) = x_h;
            C_x(3,3) = r_h;
            C_x(4,2) = -r_h;
            C_x(4,3) = x_h;
            C_x(4,4) = sqrt(g22)/a_av_x;

            T_x = C_x * N;

            dQ_x1 = [Q1(n2+2,n3);Q2(n2+2,n3);Q3(n2+2,n3);Q4(n2+2,n3)] - [Q1(n2+1,n3);Q2(n2+1,n3);Q3(n2+1,n3);Q4(n2+1,n3)];
            dQ_x2 = [Q1(n2+1,n3);Q2(n2+1,n3);Q3(n2+1,n3);Q4(n2+1,n3)] - [Q1(n2,n3);Q2(n2,n3);Q3(n2,n3);Q4(n2,n3)];
            dQ_x3 = [Q1(n2,n3);Q2(n2,n3);Q3(n2,n3);Q4(n2,n3)] - [Q1(n2-1,n3);Q2(n2-1,n3);Q3(n2-1,n3);Q4(n2-1,n3)];

            dW_x1 = T_x\dQ_x1;
            dW_x2 = T_x\dQ_x2;
            dW_x3 = T_x\dQ_x3;

            dE1 = T_x*[minmod(Lambda_x_m(1)*dW_x1(1), beta_comp*Lambda_x_m(1)*dW_x2(1));...
                       minmod(Lambda_x_m(2)*dW_x1(2), beta_comp*Lambda_x_m(2)*dW_x2(2));...
                       minmod(Lambda_x_m(3)*dW_x1(3), beta_comp*Lambda_x_m(3)*dW_x2(3));...
                       minmod(Lambda_x_m(4)*dW_x1(4), beta_comp*Lambda_x_m(4)*dW_x2(4))];
            dE2 = T_x*[minmod(Lambda_x_m(1)*dW_x2(1), beta_comp*Lambda_x_m(1)*dW_x1(1));...
                       minmod(Lambda_x_m(2)*dW_x2(2), beta_comp*Lambda_x_m(2)*dW_x1(2));...
                       minmod(Lambda_x_m(3)*dW_x2(3), beta_comp*Lambda_x_m(3)*dW_x1(3));...
                       minmod(Lambda_x_m(4)*dW_x2(4), beta_comp*Lambda_x_m(4)*dW_x1(4))];
            dE3 = T_x*[minmod(Lambda_x_p(1)*dW_x2(1), beta_comp*Lambda_x_p(1)*dW_x3(1));...
                       minmod(Lambda_x_p(2)*dW_x2(2), beta_comp*Lambda_x_p(2)*dW_x3(2));...
                       minmod(Lambda_x_p(3)*dW_x2(3), beta_comp*Lambda_x_p(3)*dW_x3(3));...
                       minmod(Lambda_x_p(4)*dW_x2(4), beta_comp*Lambda_x_p(4)*dW_x3(4))];
            dE4 = T_x*[minmod(Lambda_x_p(1)*dW_x3(1), beta_comp*Lambda_x_p(1)*dW_x2(1));...
                       minmod(Lambda_x_p(2)*dW_x3(2), beta_comp*Lambda_x_p(2)*dW_x2(2));...
                       minmod(Lambda_x_p(3)*dW_x3(3), beta_comp*Lambda_x_p(3)*dW_x2(3));...
                       minmod(Lambda_x_p(4)*dW_x3(4), beta_comp*Lambda_x_p(4)*dW_x2(4))];

            E1_half(n2,n3) = 1/2*(E1(n2+1,n3) + E1(n2,n3)) - 1/2*(dE3(1)-dE2(1)) - (1-phi)/4*dE1(1) - (1+phi)/4*dE2(1) + (1+phi)/4*dE3(1) + (1-phi)/4*dE4(1);
            E2_half(n2,n3) = 1/2*(E2(n2+1,n3) + E2(n2,n3)) - 1/2*(dE3(2)-dE2(2)) - (1-phi)/4*dE1(2) - (1+phi)/4*dE2(2) + (1+phi)/4*dE3(2) + (1-phi)/4*dE4(2);
            E3_half(n2,n3) = 1/2*(E3(n2+1,n3) + E3(n2,n3)) - 1/2*(dE3(3)-dE2(3)) - (1-phi)/4*dE1(3) - (1+phi)/4*dE2(3) + (1+phi)/4*dE3(3) + (1-phi)/4*dE4(3);
            E4_half(n2,n3) = 1/2*(E4(n2+1,n3) + E4(n2,n3)) - 1/2*(dE3(4)-dE2(4)) - (1-phi)/4*dE1(4) - (1+phi)/4*dE2(4) + (1+phi)/4*dE3(4) + (1-phi)/4*dE4(4);
            
            % η方向の計算
            % セル間の平均値の計算
            % Expressing rho, u, v and H in the basis of w (Roe's paper)
            D_h = sqrt(rho(n2+1,n3)/rho(n2,n3));
            u_av_h = (u(n2,n3)+D_h*u(n2+1,n3))/(1+D_h); 
            v_av_h = (v(n2,n3)+D_h*v(n2+1,n3))/(1+D_h); 
            H_av_h = (H(n2,n3)+D_h*H(n2+1,n3))/(1+D_h); 
            a_av_h = sqrt((gamma-1)*(H_av_h-1/2*(u_av_h.^2 + v_av_h.^2))); %Square root of the derivative of P br rho, used in the eigenvalues of the Jacobian matrix

            % ヤコビアン行列の固有値計算 
            V_av_h = u_av_h*h_x(n2,n3) + v_av_h*h_r(n2,n3);
            lambda_h(1) = V_av_h;
            lambda_h(2) = V_av_h;
            lambda_h(3) = V_av_h + a_av_h*sqrt(h_x(n2,n3)^2+h_r(n2,n3)^2);
            lambda_h(4) = V_av_h - a_av_h*sqrt(h_x(n2,n3)^2+h_r(n2,n3)^2);

            Lambda_h_p = diag((lambda_h + abs(lambda_h))/2);
            Lambda_h_m = diag((lambda_h - abs(lambda_h))/2);

            N = zeros(4,4);
            N(1,1) = 1;
            N(2,1) = -u_av_h;
            N(2,2) = 1;
            N(3,1) = -v_av_h;
            N(3,3) = 1;
            N(4,1) = (gamma-1)*(u_av_h.^2 + v_av_h.^2)/2;
            N(4,2) = -(gamma-1)*u_av_h;
            N(4,3) = -(gamma-1)*v_av_h;
            N(4,4) = gamma-1;

            C_h = zeros(4,4);
            C_h(1,1) = 1;
            C_h(1,4) = -1/a_av_h^2;
            C_h(2,2) = x_xi;
            C_h(2,3) = r_xi;
            C_h(3,2) = -r_xi;
            C_h(3,3) = x_xi;
            C_h(3,4) = sqrt(g11)/a_av_h;
            C_h(4,2) = r_xi;
            C_h(4,3) = -x_xi;
            C_h(4,4) = sqrt(g11)/a_av_h;

            T_h = C_h * N;

            dQ_h1 = [Q1(n2,n3+2);Q2(n2,n3+2);Q3(n2,n3+2);Q4(n2,n3+2)] - [Q1(n2,n3+1);Q2(n2,n3+1);Q3(n2,n3+1);Q4(n2,n3+1)];
            dQ_h2 = [Q1(n2,n3+1);Q2(n2,n3+1);Q3(n2,n3+1);Q4(n2,n3+1)] - [Q1(n2,n3);Q2(n2,n3);Q3(n2,n3);Q4(n2,n3)];
            dQ_h3 = [Q1(n2,n3);Q2(n2,n3);Q3(n2,n3);Q4(n2,n3)] - [Q1(n2,n3-1);Q2(n2,n3-1);Q3(n2,n3-1);Q4(n2,n3-1)];

            dW_h1 = T_h\dQ_h1;
            dW_h2 = T_h\dQ_h2;
            dW_h3 = T_h\dQ_h3;

            dF1 = T_h*[minmod(Lambda_h_m(1)*dW_h1(1), beta_comp*Lambda_h_m(1)*dW_h2(1));...
                       minmod(Lambda_h_m(2)*dW_h1(2), beta_comp*Lambda_h_m(2)*dW_h2(2));...
                       minmod(Lambda_h_m(3)*dW_h1(3), beta_comp*Lambda_h_m(3)*dW_h2(3));...
                       minmod(Lambda_h_m(4)*dW_h1(4), beta_comp*Lambda_h_m(4)*dW_h2(4))];
            dF2 = T_h*[minmod(Lambda_h_m(1)*dW_h2(1), beta_comp*Lambda_h_m(1)*dW_h1(1));...
                       minmod(Lambda_h_m(2)*dW_h2(2), beta_comp*Lambda_h_m(2)*dW_h1(2));...
                       minmod(Lambda_h_m(3)*dW_h2(3), beta_comp*Lambda_h_m(3)*dW_h1(3));...
                       minmod(Lambda_h_m(4)*dW_h2(4), beta_comp*Lambda_h_m(4)*dW_h1(4))];
            dF3 = T_h*[minmod(Lambda_h_p(1)*dW_h2(1), beta_comp*Lambda_h_p(1)*dW_h3(1));...
                       minmod(Lambda_h_p(2)*dW_h2(2), beta_comp*Lambda_h_p(2)*dW_h3(2));...
                       minmod(Lambda_h_p(3)*dW_h2(3), beta_comp*Lambda_h_p(3)*dW_h3(3));...
                       minmod(Lambda_h_p(4)*dW_h2(4), beta_comp*Lambda_h_p(4)*dW_h3(4))];
            dF4 = T_h*[minmod(Lambda_h_p(1)*dW_h3(1), beta_comp*Lambda_h_p(1)*dW_h2(1));...
                       minmod(Lambda_h_p(2)*dW_h3(2), beta_comp*Lambda_h_p(2)*dW_h2(2));...
                       minmod(Lambda_h_p(3)*dW_h3(3), beta_comp*Lambda_h_p(3)*dW_h2(3));...
                       minmod(Lambda_h_p(4)*dW_h3(4), beta_comp*Lambda_h_p(4)*dW_h2(4))];

            F1_half(n2,n3) = 1/2*(F1(n2,n3+1) + F1(n2,n3)) - 1/2*(dF3(1)-dF2(1)) - (1-phi)/4*dF1(1) - (1+phi)/4*dF2(1) + (1+phi)/4*dF3(1) + (1-phi)/4*dF4(1);
            F2_half(n2,n3) = 1/2*(F2(n2,n3+1) + F2(n2,n3)) - 1/2*(dF3(2)-dF2(2)) - (1-phi)/4*dF1(2) - (1+phi)/4*dF2(2) + (1+phi)/4*dF3(2) + (1-phi)/4*dF4(2);
            F3_half(n2,n3) = 1/2*(F3(n2,n3+1) + F3(n2,n3)) - 1/2*(dF3(3)-dF2(3)) - (1-phi)/4*dF1(3) - (1+phi)/4*dF2(3) + (1+phi)/4*dF3(3) + (1-phi)/4*dF4(3);
            F4_half(n2,n3) = 1/2*(F4(n2,n3+1) + F4(n2,n3)) - 1/2*(dF3(4)-dF2(4)) - (1-phi)/4*dF1(4) - (1+phi)/4*dF2(4) + (1+phi)/4*dF3(4) + (1-phi)/4*dF4(4);

            T_hT_x = [[1 0 0 0];...
                      [0 J(n2,n3)/2/g22 g12/g22 -J(n2,n3)/g22];...
                      [0 G0/2-g12/2/g22 J(n2,n3)/g22 G0/2+g12/2/g22];...
                      [0 G0/2+g12/2/g22 -J(n2,n3)/g22 G0/2-g12/2/g22]];

            % ソースタームのヤコビ行列
            omega = gamma*e(n2,n3)/rho(n2,n3);
            phi2 = (gamma-1)/2*(u(n2,n3)^2+v(n2,n3)^2);
            C = [[0 0 1 0];...
                 [-u(n2,n3)*v(n2,n3) v(n2,n3) u(n2,n3) 0];...
                 [-v(n2,n3)^2 0 2*v(n2,n3) 0];...
                 [v(n2,n3)*(2*phi2-omega) -(gamma-1)*u(n2,n3)*v(n2,n3) omega-phi2-(gamma-1)*v(n2,n3)^2 gamma*v(n2,n3)]];

            % dQの計算            
            rhs1(n2,n3) = -dt/(1+theta_2)*((E1_half(n2,n3)-E1_half(n2-1,n3))/dx + (F1_half(n2,n3)-F1_half(n2,n3-1))/dh + H1(n2,n3)/r(n3)) + theta_2/(1+theta_2)*dQ_1(n2,n3);
            rhs2(n2,n3) = -dt/(1+theta_2)*((E2_half(n2,n3)-E2_half(n2-1,n3))/dx + (F2_half(n2,n3)-F2_half(n2,n3-1))/dh + H2(n2,n3)/r(n3)) + theta_2/(1+theta_2)*dQ_2(n2,n3);
            rhs3(n2,n3) = -dt/(1+theta_2)*((E3_half(n2,n3)-E3_half(n2-1,n3))/dx + (F3_half(n2,n3)-F3_half(n2,n3-1))/dh + H3(n2,n3)/r(n3)) + theta_2/(1+theta_2)*dQ_3(n2,n3);
            rhs4(n2,n3) = -dt/(1+theta_2)*((E4_half(n2,n3)-E4_half(n2-1,n3))/dx + (F4_half(n2,n3)-F4_half(n2,n3-1))/dh + H4(n2,n3)/r(n3)) + theta_2/(1+theta_2)*dQ_4(n2,n3) + theta_1/(1+theta_2)*dt*W(n2,n3);%/r(n3)

            dQ1 = T_x\[rhs1(n2,n3);rhs2(n2,n3);rhs3(n2,n3);rhs4(n2,n3)];
            dQ1_1(n2,n3) = dQ1(1);
            dQ1_2(n2,n3) = dQ1(2);
            dQ1_3(n2,n3) = dQ1(3);
            dQ1_4(n2,n3) = dQ1(4);

            dQ1_m = [dQ1_1(n2-1,n3);dQ1_2(n2-1,n3);dQ1_3(n2-1,n3);dQ1_4(n2-1,n3)];
            dQ1_p = [dQ1_1(n2+1,n3);dQ1_2(n2+1,n3);dQ1_3(n2+1,n3);dQ1_4(n2+1,n3)];
            dQ2 = (eye(4)+theta_1/(1+theta_2)*dt/2/dx*diag(abs(lambda_x)))\(dQ1+theta_1/(1+theta_2)*dt/2/dx*(Lambda_x_p*dQ1_m-Lambda_x_m*dQ1_p));
            dQ2_1(n2,n3) = dQ2(1);
            dQ2_2(n2,n3) = dQ2(2);
            dQ2_3(n2,n3) = dQ2(3);
            dQ2_4(n2,n3) = dQ2(4);

            dQ3 = T_hT_x\dQ2;
            dQ3_1(n2,n3) = dQ3(1);
            dQ3_2(n2,n3) = dQ3(2);
            dQ3_3(n2,n3) = dQ3(3);
            dQ3_4(n2,n3) = dQ3(4);

            dQ3_m = [dQ3_1(n2,n3-1);dQ3_2(n2,n3-1);dQ3_3(n2,n3-1);dQ3_4(n2,n3-1)];
            dQ3_p = [dQ3_1(n2,n3+1);dQ3_2(n2,n3+1);dQ3_3(n2,n3+1);dQ3_4(n2,n3+1)];
            dQ4 = (eye(4)+theta_1/(1+theta_2)*dt/2/dh*diag(abs(lambda_h)))\(dQ3+theta_1/(1+theta_2)*dt/2/dh*(Lambda_h_p*dQ3_m-Lambda_h_m*dQ3_p));
            dQ4_1(n2,n3) = dQ4(1);
            dQ4_2(n2,n3) = dQ4(2);
            dQ4_3(n2,n3) = dQ4(3);
            dQ4_4(n2,n3) = dQ4(4);

            dQ5 = T_h\dQ4;
            dQ5_1(n2,n3) = dQ5(1);
            dQ5_2(n2,n3) = dQ5(2);
            dQ5_3(n2,n3) = dQ5(3);
            dQ5_4(n2,n3) = dQ5(4);

            dQ = (eye(4)+theta_1/(1+theta_2)*dt/r(n3)*C)\dQ5;
            dQ_1(n2,n3) = dQ(1);
            dQ_2(n2,n3) = dQ(2);
            dQ_3(n2,n3) = dQ(3);
            dQ_4(n2,n3) = dQ(4);
        end
    end
    
    if (mod(n1,ng) == 0)
        I = I+1;
        Per = n1/nt*100;
        t_num = int64(nt/n1);
        t_txt(t_num,1)=t*1e6;
        disp([num2str(Per),'%']) % 進捗
        disp(['time: ',num2str(t*1e6),' us']) % 時刻
        f0 = figure;
        f1 = figure;
        f2 = figure;
        f3 = figure;
        f4 = figure;

        % make Support Line Matrix
        figure(f0);
        plot(r_list, xionz*1e3);
        title('ionized wave propagation');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        xlim([0 2.5]);
        ylim([0 10]);

        % make p_t movie
        figure(f1);
        view(135,45)
        sur = surface(r_list, x_list, P / 10^5,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('Pressure Colormap /atm');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([1 50]);
        frame = getframe(gcf);
        writeVideo(vp,frame);
        
        % % make T_t movie
        % figure(f2);
        % view(135,45)
        % sur = surface(r_list, x_list, T,'FaceAlpha',0.5);
        % set(sur,'LineStyle','none')
        % title('Temperature Colormap /K');
        % ylabel('Position z /mm');
        % xlabel('Position r /mm');
        % zlim([1 5*1e4]);
        % frame = getframe(gcf);
        % writeVideo(vt,frame);
        
        % make u_t movie
        figure(f3);
        view(135,45)
        sur = surface((1:nh) *dh * 10^3, x_list, u/1e3,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('u Colormap / km/s');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([0 4]);
        frame = getframe(gcf);
        writeVideo(vu,frame);
        
        % make v_t movie
        figure(f4);
        view(135,45)
        sur = surface(r_list, x_list, v/1e3,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('v Colormap / km/s');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([0 3]);
        frame = getframe(gcf);
        writeVideo(vv,frame);
        
        % make csv file
        pressure_filename = fullfile(base_dir, 'Pressure', append('Pressure_',num2str(I),'.csv'));
        temperature_filename = fullfile(base_dir, 'Temperature', append('Temperature_',num2str(I),'.csv'));
        u_filename = fullfile(base_dir, 'U', append('U_',num2str(I),'.csv'));
        v_filename = fullfile(base_dir, 'V', append('V_',num2str(I),'.csv'));
        writematrix(P/1e5,pressure_filename)
        writematrix(T,temperature_filename)
        writematrix(u/1e3,u_filename)
        writematrix(v/1e3,v_filename)
        %M(I) = getframe(); %#ok<SAGROW>
        
        %close
    end
        
end
close(vp);
close(vt);
close(vu);
close(vv);

writematrix(t_list,fullfile(base_dir,'time_scale.csv'))
writematrix(x_list,fullfile(base_dir,'x_scale.csv'))
writematrix(r_list,fullfile(base_dir,'r_scale.csv'))
writematrix(Pplateau/1e5,fullfile(base_dir,'plateau_pressure_data.csv'))

toc


