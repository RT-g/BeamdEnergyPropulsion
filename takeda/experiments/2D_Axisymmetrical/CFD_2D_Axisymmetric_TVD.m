% 2次元軸対称流れCFD
% auther = Ryota Takeda

clear
close all
clearvars

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
dx = 0.01e-3; % m Segment length 
dr = 0.005e-3; % m Segment length % dr の分解能が大事っぽい
dt = 0.1; %s すぐに淘汰される値なので意味はない。初期値
nx = 1000; % Number of position segments, xは進展方向
nr = 1000; %2.5e-3/dr; %全長2.5 mmになるように設定 
nt = 5000; % Number of time segments
ng = 40;%40; % Graphの分割数を決める値
CFL = 0.5; % クーラン数 1以下にする 0.5~0.99くらいで，小さくしすぎると進みが遅い
delta = 0.1; %エントロピー補正量？大体10%
select_gas = 'air';

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

P = zeros(nx,nr);
Pplateau = zeros(nt,nr);
u = zeros(nx,nr);
v = zeros(nx,nr);
T = zeros(nx,nr);
rho = zeros(nx,nr);
H = zeros(nx,nr);
t_list = zeros(nt,1);
x_list = (1:nx) * dx * 1e3; % mm
r_list = (1:nr) * dr * 1e3; % mm
uionz_list = zeros(nt,1);
xionz = zeros(nr,1);
t_txt = zeros(nt/ng,1);

P(:,:) = P_0;
Pplateau(:,:) = P_0;
u(:,:) = u_0;
v(:,:) = v_0;
rho(:,:) = P_0/R/T_0;

T(:,:) = T_0;
T_x = zeros(1,nx);
T_r = zeros(1,nr);
T_x(1,:) = T_0;
T_r(1,:) = T_0;

% Vector of states
U1 = ones(nx,nr).*rho; % Mass
U2 = ones(nx,nr).*rho.*u; % Momentum x
U3 = ones(nx,nr).*rho.*v; %Momentum r
U4 = ones(nx,nr).*P/(gamma-1)+1/2.*rho.*(u.^2+v.^2); % Energy

% New calculated vector of states
U1_cal = ones(nx,nr).*rho; % Mass
U2_cal = zeros(nx,nr);
U3_cal = zeros(nx,nr);
U4_cal = ones(nx,nr).*P/(gamma-1);

% Vector of flux
F1 = ones(nx,nr).*rho.*u; % 0
F2 = ones(nx,nr).*rho.*u.^2+P; % P
F3 = ones(nx,nr).*rho.*u.*v; % 0
F4 = (U4+P).*u; % 0

G1 = ones(nx,nr).*rho.*v; %0
G2 = ones(nx,nr).*rho.*v.*u; %0
G3 = ones(nx,nr).*rho.*v.^2+P; %P
G4 = (U4+P).*v; %0

% Vector of flux between cells
F1_half = zeros(nx,nr);
F2_half = ones(nx,nr).*rho.*u.^2+P; % P
F3_half = zeros(nx,nr);
F4_half = zeros(nx,nr);

G1_half = zeros(nx,nr);
G2_half = zeros(nx,nr);
G3_half = ones(nx,nr).*rho.*v.^2+P; %P
G4_half = zeros(nx,nr);

% Calculation
%eta = eta_trans
eta = 0.05; % eta_trans ではどうやら値が大きすぎる
t = 0;
x_laser0 = 0;
umax = u_ionz0;
vmax = u_ionz0;
I = 0;

for n1 = 1:nt
    % 時間刻み幅の設定。ここ次第で計算がいくらでも長くなる    
    dt = CFL * min(dx/umax,  dr/vmax); %s, 計算位置>波面位置となるようにdtを決定
    dtdx = dt/dx;
    dtdr = dt/dr;

    t = t+dt; %s
    t_list(n1,1) = t * 1e6; %us

    % レーザー強度の時間減衰(10J) 単位はMW
    Power_laser = 8.15*exp(-0.866*t*1e6); 
    % 正確に再現すると温度などが高く出すぎる
    % if (t*1e6<0.085)
    %     Power_laser = 287.9222823*t*1e6+0.0005175756469;
    % elseif (t*1e6<0.125)
    %     Power_laser =-476.3277981*t*1e6+66.1043297;
    % else
    %     Power_laser = 8.15*exp(-0.866*t*1e6); 
    % end  
    
    % 進展方向のレーザー強度の変化を考慮しない場合のビーム半径
    W_G = W_G0*1e-3; %m
    W_T = W_T0*1e-3; %m

    % 波頭の値
    S_laser0 = R_peak * Power_laser / 4 / W_G / W_T * 1e-3; % GW/m^2 波頭のレーザー強度
    u_ionz_line1 = slope * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept * a; %m/s 波頭の速度
    u_ionz_line3 = slope_low * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept_low * a;
    if (u_ionz_line1 > u_ionz_line3)
        u_ionz = u_ionz_line1; % Line1, 松井さん博論
        b = intercept;
    else
        u_ionz = u_ionz_line3; % Line3, 松井さん博論
        b = intercept_low;
    end
    x_laser0 = x_laser0 + u_ionz * dt; %m 電離波面の波頭進展位置=Laserによって加熱される最初の位置
    x_laser = x_laser0; %m 横方向の加熱位置を決めるための処理。初期値はその時刻における波頭位置とする

    % 半径方向の計算
    for n3 = 2:nr-1
        r = n3*dr; %m

        % ガウシアン分布/トップハット分布を仮定 二次元極座標のためx,yをそれぞれr/sqrt(2)としている。
        G_r = A_G*exp(-2*(r*1e3/sqrt(2)/sigma_G1)^4)+B_G*exp(-2*(r*1e3/sqrt(2)/sigma_G2)^2);
        T_r = A_T*exp(-2*(r*1e3/sqrt(2)/sigma_T1)^4)+B_T*exp(-2*(r*1e3/sqrt(2)/sigma_T2)^2);
        S_laser = R_peak * Power_laser/4/W_G/W_T*G_r*T_r*1e-3; %局所レーザー強度, GW/m2
        x_laser = x_laser - sqrt((S_laser0/S_laser) ^(2*b/(1-b))-1)*dr;  %m 松井さんD論より、横方向の波面位置を積分により導出。最初の方は外側は負の値            
        if x_laser > 0
            xionz(n3,1) = x_laser; % 電離波面の進展位置。加熱箇所を可視化する。
            break_r = r;
        end

        % 進展方向の計算
        for n2 = 2:nx-1
            x = n2*dx; %m
            % % その位置におけるビーム半径の計算。大した距離ではないのであまり計算する必要はない。
            % W_G = sqrt((W_G0*1e-3)^2+M2_G^2*(lambda*1e-6)^2*x^2/pi^2/(W_G0*1e-3)^2); %m
            % W_T = sqrt((W_T0*1e-3)^2+M2_T^2*(lambda*1e-6)^2*x^2/pi^2/(W_T0*1e-3)^2); %m
            
            % % 波頭の値
            % S_laser0 = R_peak *Power_laser/4/W_G/W_T*1e-3; % GW/m^2 波頭のレーザー強度
            % u_ionz = a*(S_laser0)^b*10^3; %m/s 波頭の速度
            % x_laser0 = x_laser0+u_ionz*dt; %電離波面の波頭進展位置=Laserによって加熱される最初の位置
            % x_laser = x_laser0; %横方向の加熱位置を決めるための処理。初期値はその時刻における波頭位置とする
        
            % Input Power Definition
            % 加熱領域をx_laserよりも前にしてしまうと計算が壊れるので実際に考えられる値よりも少し進めたほうがいい?
            if x > x_laser-l && x < x_laser
                w = eta * S_laser / l * 1e9; %W/m3
            else
                w = 0;
            end

            % Euler Calculation 有限体積法、一次後退差分(風上差分)
            rn = abs(r);
            U1_cal(n2,n3) = U1(n2,n3) - dtdx*(F1_half(n2,n3)-F1_half(n2-1,n3)) - dtdr*(G1_half(n2,n3)-G1_half(n2,n3-1))-dt/rn*G1(n2,n3);
            U2_cal(n2,n3) = U2(n2,n3) - dtdx*(F2_half(n2,n3)-F2_half(n2-1,n3)) - dtdr*(G2_half(n2,n3)-G2_half(n2,n3-1))-dt/rn*G2(n2,n3) ;
            U3_cal(n2,n3) = U3(n2,n3) - dtdx*(F3_half(n2,n3)-F3_half(n2-1,n3)) - dtdr*(G3_half(n2,n3)-G3_half(n2,n3-1))+dt/rn*(P(n2,n3)-G3(n2,n3));
            U4_cal(n2,n3) = U4(n2,n3) - dtdx*(F4_half(n2,n3)-F4_half(n2-1,n3)) - dtdr*(G4_half(n2,n3)-G4_half(n2,n3-1))-dt/rn*G4(n2,n3) + w*dt;
            
        end
    end    
    

    % Boundary Conditions
    
    % Left(回転軸)
    for i = 2:nx-1
        %T_x(i) = (gamma-1)*(U4_cal(i,2)-(1/2)*(U2_cal(i,2).^2+U3_cal(i,2).^2)./U1_cal(i,2))./U1_cal(i,2)/R; % Temperature to set the condition
        U1_cal(i,1) = U1_cal(i,2); 
        U2_cal(i,1) = U2_cal(i,2); 
        U3_cal(i,1) = 0; %中心ではur=0のはず
        U4_cal(i,1) = U4_cal(i,2);
    end
    
    % Right(流出)
    for i = 2:nx-1
        % T_x(i) = (gamma-1)*(U4_cal(i,nr-1)-(1/2)*(U2_cal(i,nr-1).^2+U3_cal(i,nr-1).^2)./U1_cal(i,nr-1))./U1_cal(i,nr-1)/R; % Temperature to set the condition
        U1_cal(i,nr) = U1_cal(i,nr-1);%P_0/R./T_x(i);
        U2_cal(i,nr) = U2_cal(i,nr-1);%U2_cal(i,nr-1)./U1_cal(i,nr-1).*U1_cal(i,nr);
        U3_cal(i,nr) = U3_cal(i,nr-1);%U3_cal(i,nr-1)./U1_cal(i,nr-1).*U1_cal(i,nr);
        U4_cal(i,nr) = U4_cal(i,nr-1);%U4_cal(i,nr-1)./U1_cal(i,nr-1).*U1_cal(i,nr);
    end
    
    % Top(光軸方向)流出
    for i = 2:nr-1
        % T_r(i) = (gamma-1)*(U4_cal(nx-1,i)-(1/2)*(U2_cal(nx-1,i).^2+U3_cal(nx-1,i).^2)./U1_cal(nx-1,i))./U1_cal(nx-1,i)/R; % Temperature to set the condition
        U1_cal(nx,i) = U1_cal(nx-1,i); %P_0/R./T_r(i);
        U2_cal(nx,i) = U2_cal(nx-1,i);%./U1_cal(nx-1,i).*U1_cal(nx,i);
        U3_cal(nx,i) = U3_cal(nx-1,i);%./U1_cal(nx-1,i).*U1_cal(nx,i);
        U4_cal(nx,i) = U4_cal(nx-1,i);%./U1_cal(nx-1,i).*U1_cal(nx,i);
    end
    
    % Plate(アルミ板)
    U1_cal(1,:) = U1_cal(2,:);
    U2_cal(1,:) = 0;
    U3_cal(1,:) = 0;
    U4_cal(1,:) = U4_cal(2,:);
    
    % Corners
    U1_cal(nx,1) = (U1_cal(nx,2) + U1_cal(nx-1,1))/2;
    U2_cal(nx,1) = (U2_cal(nx,2) + U2_cal(nx-1,1))/2;
    U3_cal(nx,1) = (U3_cal(nx,2) + U3_cal(nx-1,1))/2;
    U4_cal(nx,1) = (U4_cal(nx,2) + U4_cal(nx-1,1))/2;
    
    U1_cal(nx,nr) = (U1_cal(nx-1,nr) + U1_cal(nx,nr-1))/2;
    U2_cal(nx,nr) = (U2_cal(nx-1,nr) + U2_cal(nx,nr-1))/2;
    U3_cal(nx,nr) = (U3_cal(nx-1,nr) + U3_cal(nx,nr-1))/2;
    U4_cal(nx,nr) = (U4_cal(nx-1,nr) + U4_cal(nx,nr-1))/2;
    

    % 各行列の更新
    U1 = U1_cal;
    U2 = U2_cal;
    U3 = U3_cal;
    U4 = U4_cal;
    
    rho = U1;
    u = U2./rho; %軸方向速度
    v = U3./rho; %半径方向速度
    P = (gamma-1)*(U4-(1/2)*rho.*(u.^2+v.^2));
    Pplateau(n1,:) = P(1,:);
    % disp(Pplateau)
    H = (U4+P)./rho;
    T = P./rho./R; % Temperature

    F1 = U2;
    F2 = rho.*u.^2+P;
    F3 = rho.*u.*v;
    F4 = (U4+P).*u;
    
    G1 = U3;
    G2 = rho.*u.*v;
    G3 = rho.*v.^2+P;
    G4 = (U4+P).*v;
        
    umax = max(u, [], 'all');
    vmax = max(v, [], 'all');
    

    % セル間の平均値の計算 ゴドノフ法
    for n2 = 2:nx-2
        for n3 = 2:nr-2
            if n3*dr > break_r + 20*dr ||  n2*dx > x_laser0 + 20 * dx
                F1_half(n2,n3) = 1/2*(F1(n2+1,n3)+F1(n2,n3));
                F2_half(n2,n3) = 1/2*(F2(n2+1,n3)+F2(n2,n3));
                F3_half(n2,n3) = 1/2*(F3(n2+1,n3)+F3(n2,n3));
                F4_half(n2,n3) = 1/2*(F4(n2+1,n3)+F4(n2,n3));
                G1_half(n2,n3) = 1/2*(G1(n2,n3+1)+G1(n2,n3));
                G2_half(n2,n3) = 1/2*(G2(n2,n3+1)+G2(n2,n3));
                G3_half(n2,n3) = 1/2*(G3(n2,n3+1)+G3(n2,n3));
                G4_half(n2,n3) = 1/2*(G4(n2,n3+1)+G4(n2,n3));
                continue
            end
            % Direction x
            % Solving the Riemann problem
            % 2次精度風上TVD - Harten-Yee TVD Scheme

            rsx = riemannsolver(gamma,rho(n2:n2+1,n3),u(n2:n2+1,n3), v(n2:n2+1,n3),H(n2:n2+1,n3),P(n2:n2+1,n3),U1(n2:n2+1,n3),U2(n2:n2+1,n3),U3(n2:n2+1,n3),U4(n2:n2+1,n3),'F');
            rslx = riemannsolver(gamma,rho(n2-1:n2,n3),u(n2-1:n2,n3),v(n2-1:n2,n3),H(n2-1:n2,n3),P(n2-1:n2,n3),U1(n2-1:n2,n3),U2(n2-1:n2,n3),U3(n2-1:n2,n3),U4(n2-1:n2,n3),'F');
            rsrx = riemannsolver(gamma,rho(n2+1:n2+2,n3),u(n2+1:n2+2,n3),v(n2+1:n2+2,n3),H(n2+1:n2+2,n3),P(n2+1:n2+2,n3),U1(n2+1:n2+2,n3),U2(n2+1:n2+2,n3),U3(n2+1:n2+2,n3),U4(n2+1:n2+2,n3),'F');
            alpha_x= cell2mat(rsx(1));
            R_x= cell2mat(rsx(2));
            lambda_x = cell2mat(rsx(3));
            alpha_lx = cell2mat(rslx(1));
            alpha_rx = cell2mat(rsrx(1));
            ph_x = zeros(4,1);

            for ni = 1:4
                % 線形場と非線形場で分けたほうがいい説もある。非線形はminmod, 線形はsuperbee
                if ismember(ni, 1:2)
                    g_lx = superbee(alpha_lx(ni), alpha_x(ni));
                    g_rx = superbee(alpha_x(ni), alpha_rx(ni));
                else
                    g_lx = minmod(alpha_lx(ni), alpha_x(ni));
                    g_rx = minmod(alpha_x(ni), alpha_rx(ni));
                end
                if alpha_x(ni) == 0
                    beta_x = 0;
                else
                    beta_x = tvd_sigma(lambda_x(ni),dtdx,delta)*(g_rx-g_lx)/alpha_x(ni);
                end
                ph_x(ni) = tvd_sigma(lambda_x(ni),dtdx,delta)*(g_lx+g_rx) - enthalpy_fix(lambda_x(ni)+beta_x, delta)*alpha_x(ni);
            end

            if (umax<max(abs(lambda_x)))
                umax = max(abs(lambda_x));
            end
            tvd_x = R_x*ph_x;
            % tvd_x = -R_x*diag(abs(lambda_x))*alpha_x; % ここをTVDスキームに変えるだけでTVDにできる

            F1_half(n2,n3) = 1/2*(F1(n2+1,n3)+F1(n2,n3)+tvd_x(1));
            F2_half(n2,n3) = 1/2*(F2(n2+1,n3)+F2(n2,n3)+tvd_x(2));
            F3_half(n2,n3) = 1/2*(F3(n2+1,n3)+F3(n2,n3)+tvd_x(3));
            F4_half(n2,n3) = 1/2*(F4(n2+1,n3)+F4(n2,n3)+tvd_x(4));
            

            % Direction r
            % % Solving the linear problem
            % % 2次精度風上TVD
            rsr = riemannsolver(gamma,rho(n2,n3:n3+1),u(n2,n3:n3+1), v(n2,n3:n3+1),H(n2,n3:n3+1),P(n2,n3:n3+1),U1(n2,n3:n3+1),U2(n2,n3:n3+1),U3(n2,n3:n3+1),U4(n2,n3:n3+1),'G');
            rslr = riemannsolver(gamma,rho(n2,n3-1:n3),u(n2,n3-1:n3),v(n2,n3-1:n3),H(n2,n3-1:n3),P(n2,n3-1:n3),U1(n2,n3-1:n3),U2(n2,n3-1:n3),U3(n2,n3-1:n3),U4(n2,n3-1:n3),'G');
            rsrr = riemannsolver(gamma,rho(n2,n3+1:n3+2),u(n2,n3+1:n3+2),v(n2,n3+1:n3+2),H(n2,n3+1:n3+2),P(n2,n3+1:n3+2),U1(n2,n3+1:n3+2),U2(n2,n3+1:n3+2),U3(n2,n3+1:n3+2),U4(n2,n3+1:n3+2),'G');
            alpha_r= cell2mat(rsr(1));
            R_r= cell2mat(rsr(2));
            lambda_r = cell2mat(rsr(3));
            alpha_lr = cell2mat(rslr(1));
            alpha_rr = cell2mat(rsrr(1));
            ph_r = zeros(4,1);

            for ni = 1:4
                if ismember(ni, 1:2)
                    g_lr = superbee(alpha_lr(ni), alpha_r(ni));
                    g_rr = superbee(alpha_r(ni), alpha_rr(ni));
                else
                    g_lr = minmod(alpha_lr(ni), alpha_r(ni));
                    g_rr = minmod(alpha_r(ni), alpha_rr(ni));
                end
                if alpha_r(ni) == 0
                    beta_r = 0;
                else
                    beta_r = tvd_sigma(lambda_r(ni),dtdr,delta)*(g_rr-g_lr)/alpha_r(ni);
                end
                ph_r(ni) = tvd_sigma(lambda_r(ni),dtdr,delta)*(g_lr+g_rr) - enthalpy_fix(lambda_r(ni)+beta_r, delta)*alpha_r(ni);
            end

            if (vmax<max(abs(lambda_r)))
                vmax = max(abs(lambda_r));
            end
            tvd_r = R_r*ph_r;
            % tvd_r = -R_r*diag(abs(lambda_r))*alpha_r; % ここをTVDスキームに変えるだけでTVDにできる

            G1_half(n2,n3) = 1/2*(G1(n2,n3+1)+G1(n2,n3)+tvd_r(1));
            G2_half(n2,n3) = 1/2*(G2(n2,n3+1)+G2(n2,n3)+tvd_r(2));
            G3_half(n2,n3) = 1/2*(G3(n2,n3+1)+G3(n2,n3)+tvd_r(3));
            G4_half(n2,n3) = 1/2*(G4(n2,n3+1)+G4(n2,n3)+tvd_r(4));
        end
    end

    % Boundary Conditions

    F1_half(1,:) = F1_half(2,:);
    F2_half(1,:) = F2_half(2,:);
    F3_half(1,:) = F3_half(2,:);
    F4_half(1,:) = F4_half(2,:);

    G1_half(1,:) = G1_half(2,:);
    G2_half(1,:) = G2_half(2,:);
    G3_half(1,:) = G3_half(2,:);
    G4_half(1,:) = G4_half(2,:);

    F1_half(:,1) = F1_half(:,2);
    F2_half(:,1) = F2_half(:,2);
    F3_half(:,1) = F3_half(:,2);
    F4_half(:,1) = F4_half(:,2);

    G1_half(:,1) = G1_half(:,2);
    G2_half(:,1) = G2_half(:,2);
    G3_half(:,1) = G3_half(:,2);
    G4_half(:,1) = G4_half(:,2);

    
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
        ylim([0 2.5]);

        % make p_t movie
        figure(f1);
        view(135,45)
        sur = surface(r_list, x_list, P / 10^5,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('Pressure Colormap /atm');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([1 100]);
        frame = getframe(gcf);
        writeVideo(vp,frame);
        
        % make T_t movie
        figure(f2);
        view(135,45)
        sur = surface(r_list, x_list, T,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('Temperature Colormap /K');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([1 1*1e5]);
        frame = getframe(gcf);
        writeVideo(vt,frame);
        
        % make u_t movie
        figure(f3);
        view(135,45)
        sur = surface((1:nr) *dr * 10^3, x_list, u/1e3,'FaceAlpha',0.5);
        set(sur,'LineStyle','none')
        title('u Colormap / km/s');
        ylabel('Position z /mm');
        xlabel('Position r /mm');
        zlim([0 5]);
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

writematrix(t_txt,fullfile(base_dir,'time_scale.csv'))
writematrix(x_list,fullfile(base_dir,'x_scale.csv'))
writematrix(r_list,fullfile(base_dir,'r_scale.csv'))
writematrix(Pplateau/1e5,fullfile(base_dir,'plateau_pressure_data.csv'))

disp(['time: ',num2str(t*1e6),' us']) % 時刻

toc


