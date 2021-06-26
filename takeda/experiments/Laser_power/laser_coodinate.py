from math import exp, sqrt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, gridspec

dx = 0.01e-3  # m Segment length
dh = 0.005e-3  # m Segment length  # dr の分解能が大事っぽい
dt = 0.1  # s すぐに淘汰される値なので意味はない。初期値
nz = 1000  # Number of position segments, xは進展方向
nr = 1000  # 2.5e-3/dr  #全長2.5 mmになるように設定
nt = 1000  # Number of time segments
ng = 40  # 40  # Graphの分割数を決める値
CFL = 0.7  # クーラン数 1以下にする 0.5~0.99くらいで，小さくしすぎると進みが遅い
eta = 0.1  # 加熱効率

A_G = 0.224517656
B_G = 0.77548
sigma_G1 = 0.84473  # mm
sigma_G2 = 1.75302  # mm
A_T = 0.764525962
B_T = 0.235474467
sigma_T1 = 1.557528296
sigma_T2 = 4.050355336
lamda = 10.6  # 単位はum
M2_G = 15
M2_T = 21
W_G0 = 1.7  # 単位はmm
W_T0 = 2.0
R_peak = 2.13
hr = 20  # Heating region
R_0 = 8314.0

# 気体種固有の値
select_gas = 'air'
gas_dict = pd.read_csv('../../data/value.csv', header=0, index_col=0).to_dict(orient='index')[select_gas]
m = gas_dict['m']
gamma_gas = gas_dict['gamma']
eta_trans = gas_dict['eta_trans']
a = gas_dict['a']
b = gas_dict['b']
a_low = gas_dict['a_low']
b_low = gas_dict['b_low']
ss = gas_dict['speed_of_sound']
rho = gas_dict['rho']
R = R_0 / m
t_terminate = 2.3

u_ionz0 = a * 524.4 ** b * 1.2 * 1e3  # m/s
P_0 = 1.013e5
T_0 = 293
rho_0 = P_0 / R / T_0
u_0 = 0
v_0 = 0
gamma = gamma_gas
I = 0
l = 0.2e-3  # m 加熱長さ レーザーは0.2 mm
eta = 0.1
t = 0
x_laser0 = 0
umax_x = u_ionz0
umax_r = u_ionz0

t_list = []
uionz_list = []
uionz_low_list = []
p_list = []
s_list = []
x_list = []

fig, (axL, axR) = plt.subplots(ncols=2, figsize=(10, 4))

# 時間刻み幅の設定。ここ次第で計算がいくらでも長くなる
dt = 1e-8  # s
intercept = 0.46
intercept_low = 0.3
a_s = (1 / sigma_G2**2 + 1 / sigma_T2**2)
slope = 0.16
slope_low = 0.68

for n1 in range(nt):
    t = t + dt  # s
    t_list.append(t * 1e6)
    # レーザー強度の時間減衰(10J) 単位はMW
    if t * 1e6 < 0.085:
        Power_laser = 287.9222823 * t * 1e6 + 0.0005175756469
    elif t * 1e6 < 0.125:
        Power_laser = -476.3277981 * t * 1e6 + 66.1043297
    else:
        Power_laser = 8.15 * exp(-0.866 * t * 1e6)
    p_list.append(Power_laser)
    # 進展方向のレーザー強度の変化を考慮しない場合のビーム半径
    W_G = W_G0 * 1e-3  # m
    W_T = W_T0 * 1e-3  # m
    # 波頭の値
    S_laser0 = R_peak * Power_laser / 4 / W_G / W_T * 1e-3  # GW/m ** 2 波頭のレーザー強度
    u_ionz_line1 = a * (S_laser0 * 1e9 / rho / (ss**3)) ** b * ss  # m/s 波頭の速度
    u_ionz_line3 = a_low * (S_laser0 * 1e9 / rho / (ss**3)) ** b_low * ss  # m/s 波頭の速度
    if u_ionz_line1 > u_ionz_line3:
        u_ionz = u_ionz_line1
        beta = b
        t_500 = t
    else:
        u_ionz = u_ionz_line3
        beta = b_low
    # u_ionz = a * S_laser0 ** b * 1e3  # m/s 無次元化していない波頭の速度
    x_laser0 = x_laser0 + u_ionz * dt  # m 電離波面の波頭進展位置=Laserによって加熱される最初の位置
    x_laser = x_laser0  # m 横方向の加熱位置を決めるための処理。初期値はその時刻における波頭位置とする
    uionz_list.append(u_ionz * 1e-3)
    # uionz_low_list.append(u_ionz_low * 1e-3)
    s_list.append(S_laser0)
    # x_list.append(x_laser0 * 1e3)
    r_list = []
    x_lateral_list = []
    s_lateral_list = []
    # 進展方向の計算
    for n3 in range(nz):
        h = n3 * dh  # m
        # 半径方向の計算
        # line1
        b_s_line1 = intercept / (1 - intercept)
        r_line1 = -(2 / (3 * a_s ** 2 * b_s_line1 ** 2 * h + sqrt(9 * a_s ** 4 * b_s_line1 ** 4 * h**2 + 4 * a_s ** 3 * b_s_line1 ** 3)))**(1 / 3) + (2 / (3 * a_s ** 2 * b_s_line1 ** 2 * h + sqrt(9 * a_s ** 4 * b_s_line1 ** 4 * h**2 + 4 * a_s ** 3 * b_s_line1 ** 3)))**(-1 / 3)/a_s/b_s_line1
        G_r_line1 = A_G*exp(-2*(r_line1*1e3/sqrt(2)/sigma_G1)**4)+B_G*exp(-2*(r_line1*1e3/sqrt(2)/sigma_G2)**2)
        T_r_line1 = A_T*exp(-2*(r_line1*1e3/sqrt(2)/sigma_T1)**4)+B_T*exp(-2*(r_line1*1e3/sqrt(2)/sigma_T2)**2)
        cos_line1 = (G_r_line1*T_r_line1)**b_s_line1
        u_ionz_line1 = slope * (S_laser0 * cos_line1 * 1e9 / rho_0 / a ** 3) ** intercept * a
        # line3
        b_s_line3 = intercept_low/(1-intercept_low)
        r_line3 = -(2/(3*a_s**2*b_s_line3**2*h + sqrt(9*a_s**4*b_s_line3**4*h**2 + 4*a_s**3*b_s_line3**3)))**(1/3) + (2/(3*a_s**2*b_s_line3**2*h + sqrt(9*a_s**4*b_s_line3**4*h**2 + 4*a_s**3*b_s_line3**3)))**(-1/3)/a_s/b_s_line3
        G_r_line3 = A_G*exp(-2*(r_line3*1e3/sqrt(2)/sigma_G1)**4)+B_G*exp(-2*(r_line3*1e3/sqrt(2)/sigma_G2)**2)
        T_r_line3 = A_T*exp(-2*(r_line3*1e3/sqrt(2)/sigma_T1)**4)+B_T*exp(-2*(r_line3*1e3/sqrt(2)/sigma_T2)**2)
        cos_line3 = (G_r_line3*T_r_line3)**b_s_line3
        u_ionz_line3 = slope_low * (S_laser0 * cos_line3 * 1e9 / rho_0 / a ** 3) ** intercept_low * a
        if (u_ionz_line1 > u_ionz_line3):
            b = intercept
            b_s = b_s_line1
            r = r_line1
        else:
            b = intercept_low
            b_s = b_s_line3
            r = r_line3
        if n3 % 25 ==0:
            r_list = []
            x_list = []
            for n2 in range(nr):
                xi = n2 * dx  # m
                r_list.append(r)
                # print((1-a_s*r**2)**(-2*b_s))
                x_list.append(xi-sqrt((1-a_s*r**2)**(-2*b_s)-1))

                # r_minus = list(map(lambda x: x * -1, r_list))
                # r_list = list(reversed(r_minus)) + r_list
                # x_lateral_list = list(reversed(x_lateral_list)) + x_lateral_list
                # s_lateral_list = list(reversed(s_lateral_list)) + s_lateral_list
            plt.plot(r_list, x_list, label=f't = {round(t*1e6,2)}')
                # axR.plot(r_list, s_lateral_list, label=f't = {round(t*1e6,2)}')
plt.show()
axL.set_xlabel('r')
axL.set_ylabel('x')
axL.set_xlim(-2, 2)
axL.set_ylim(-4, 0)
axR.set_xlabel('r')
axR.set_ylabel('S')
axR.set_xlim(-2, 2)
axL.set_aspect('equal')
# axR.set_ylim(0, 500)
axL.legend(loc='upper right')  # 凡例
axR.legend(loc='upper right')  # 凡例
plt.savefig('波面形状とレーザープロファイル履歴.png')
fig.show()

# figure()でグラフを表示する領域をつくり，figというオブジェクトにする．
fig = plt.figure(figsize=(10, 8))

# add_subplot()でグラフを描画する領域を追加する．引数は行，列，場所
spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[1, 4])
ax1 = fig.add_subplot(4, 5, 1)
ax2 = fig.add_subplot(4, 5, 6)
ax3 = fig.add_subplot(4, 5, 11)
ax4 = fig.add_subplot(4, 5, 16)
ax5 = fig.add_subplot(spec[1])

c1, c2, c3, c4 = "blue", "green", "red", "black"  # 各プロットの色
l1, l2, l3, l4, l5, l6 = "power", "intensity", "velocity", "x", 'line1', 'line3'  # 各ラベル

ax1.plot(t_list, p_list, color=c1, label=l1)
ax2.plot(t_list, s_list, color=c2, label=l2)
ax3.plot(t_list, uionz_list, color=c3, label=l3)
ax4.plot(t_list, x_list, color=c4, label=l4)
ax5.plot(s_list, uionz_list, color=c1, label=l5)
# ax5.plot(s_list, uionz_low_list, color=c2, label=l6)
ax1.legend(loc='upper right')  # 凡例
ax2.legend(loc='upper right')  # 凡例
ax3.legend(loc='upper right')  # 凡例
ax4.legend(loc='upper right')  # 凡例
ax5.legend(loc='upper right')  # 凡例

ax1.vlines([2.3], 0, 1000, 'yellow', linestyles='dashed')
ax2.vlines([2.3], 0, 1000, 'yellow', linestyles='dashed')
ax3.vlines([2.3], 0, 1000, 'yellow', linestyles='dashed')
ax4.vlines([2.3], 0, 1000, 'yellow', linestyles='dashed')
ax1.vlines([t_500*1e6], 0, 1000, 'red', linestyles='dashed')
ax2.vlines([t_500*1e6], 0, 1000, 'red', linestyles='dashed')
ax3.vlines([t_500*1e6], 0, 1000, 'red', linestyles='dashed')
ax4.vlines([t_500*1e6], 0, 1000, 'red', linestyles='dashed')
ax5.set_xscale('log')
ax5.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
ax5.set_yscale('log')
ax5.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
ax5.vlines([500], 0.1, 10, 'red', linestyles='dashed')
intencity_terminate = 8.15 * exp(-0.866 * 2.3) *  R_peak / 4 / W_G / W_T * 1e-3
ax5.vlines([intencity_terminate], 0.1, 10, 'yellow', linestyles='dashed')

ax1.set_xlim(0, 5)
ax1.set_ylim(0, 25)
ax2.set_xlim(0, 5)
ax2.set_ylim(0, 1000)
ax3.set_xlim(0, 5)
ax3.set_ylim(0, 8)
ax4.set_xlim(0, 5)
ax4.set_ylim(0, 16)
ax5.set_xlim(20, 2000)
ax5.set_ylim(0.3, 10)
# fig.tight_layout()  # レイアウトの設定
plt.savefig('laser.png')
plt.show()
