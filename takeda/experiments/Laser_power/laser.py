from math import exp
import matplotlib.pyplot as plt

dx = 0.01e-3  # m Segment length
dr = 0.005e-3  # m Segment length  # dr の分解能が大事っぽい
dt = 0.1  # s すぐに淘汰される値なので意味はない。初期値
nx = 250  # Number of position segments, xは進展方向
nr = 500  # 2.5e-3/dr  #全長2.5 mmになるように設定
nt = 5000  # Number of time segments
ng = 40  # 40  # Graphの分割数を決める値
CFL = 0.7  # クーラン数 1以下にする 0.5~0.99くらいで，小さくしすぎると進みが遅い
eta = 0.1  # 加熱効率
select_gas = 'air'

helium_a = 0.004519362
helium_b = 1.10603162
argon_a = 1.18
argon_b = 0.21
air_a = 0.354446826
air_b = 0.379696856
m_air = 28.966
m_argon = 40.0
m_helium = 4.0
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

R_0 = 8314.0

if select_gas == 'helium':
    m = m_helium
    gamma_gas = 1.67  # 単原子分子
    eta_trans = 0.54
    a = helium_a
    b = helium_b
    hr = 20  # Heating region
elif select_gas == 'argon':
    m = m_argon
    gamma_gas = 1.67  # 単原子分子
    eta_trans = 0.42
    a = argon_a
    b = argon_b
    hr = 20  # Heating region
elif select_gas == 'air':
    m = m_air
    gamma_gas = 1.4  # 2原子分子
    eta_trans = 0.48
    a = air_a
    b = air_b
    hr = 20  # Heating region

R = R_0 / m

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
p_list = []
s_list = []
x_list = []

for n1 in range(nt):
    # 時間刻み幅の設定。ここ次第で計算がいくらでも長くなる

    dt = 1e-8  # s
    t = t + dt  # s
    t_list.append(t * 1e6)
    # レーザー強度の時間減衰(10J) 単位はMW
    # Power_laser = 8.15 * exp(-0.866 * t * 1e6)
    # 正確に再現すると温度などが高く出すぎる
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
    u_ionz = a * (S_laser0) ** b * 1e3  # m/s 波頭の速度
    x_laser0 = x_laser0 + u_ionz * dt  # m 電離波面の波頭進展位置=Laserによって加熱される最初の位置
    x_laser = x_laser0  # m 横方向の加熱位置を決めるための処理。初期値はその時刻における波頭位置とする
    uionz_list.append(u_ionz * 1e-3)
    s_list.append(S_laser0)
    x_list.append(x_laser0 * 1e3)


#     # 進展方向の計算
#     for n2 in range(nz):
#         x = n2 * dx  # m
# #          # その位置におけるビーム半径の計算
# #          W_G = sqrt((W_G0*1e-3) ** 2+M2_G ** 2*(lamda*1e-6) ** 2*x ** 2/pi ** 2/(W_G0*1e-3) ** 2) # m
# #          W_T = sqrt((W_T0*1e-3) ** 2+M2_T ** 2*(lamda*1e-6) ** 2*x ** 2/pi ** 2/(W_T0*1e-3) ** 2) # m
# #          # 波頭の値
# #          S_laser0 = R_peak *Power_laser/4/W_G/W_T*1e-3 #  GW/m ** 2 波頭のレーザー強度
# #          u_ionz = a*(S_laser0) ** b*10 ** 3 # m/s 波頭の速度
# #          x_laser0 = x_laser0+u_ionz*dt # 電離波面の波頭進展位置=Laserによって加熱される最初の位置
# #          x_laser = x_laser0 # 横方向の加熱位置を決めるための処理。初期値はその時刻における波頭位置とする

#         # 半径方向の計算
#         for n3 in range(nr):
#             r = n3 * dr  # m
#             # ガウシアン分布/トップハット分布を仮定 二次元極座標のためx,yをそれぞれr/sqrt(2)としている。
#             G_r = A_G * exp(-2 * (r * 1e3 / sqrt(2) / sigma_G1) ** 4) + B_G * exp(-2 * (r * 1e3 / sqrt(2) / sigma_G2) ** 2)
#             T_r = A_T * exp(-2 * (r * 1e3 / sqrt(2) / sigma_T1) ** 4) + B_T * exp(-2 * (r * 1e3 / sqrt(2) / sigma_T2) ** 2)
#             S_laser = R_peak * Power_laser / 4 / W_G / W_T * G_r * T_r * 1e-3  # 局所レーザー強度, GW/m2
#             x_laser = x_laser - sqrt((S_laser0 / S_laser) ** (2 * b / (1 - b)) - 1) * dr  # m 松井さんD論より、横方向の波面位置を積分により導出。最初の方は外側は負の値
#             # disp(x_laser*1e3)

#             #  Input Power Definition
#             #  加熱領域をx_laserよりも前にしてしまうと計算が壊れるので実際に考えられる値よりも少し進めたほうがいい
#             #  x_laserが0以下の時は加熱領域が0とする
#             if x_laser > 0:
#                 if x > x_laser - l and x < x_laser:
#                     w = eta * S_laser / l * 1e9  # W/m3
#                 else:
#                     w = 0
#                 end
#             else
#                 w = 0

#figure()でグラフを表示する領域をつくり，figというオブジェクトにする．
fig = plt.figure()

#add_subplot()でグラフを描画する領域を追加する．引数は行，列，場所
ax1 = fig.add_subplot(4, 1, 1)
ax2 = fig.add_subplot(4, 1, 2)
ax3 = fig.add_subplot(4, 1, 3)
ax4 = fig.add_subplot(4, 1, 4)

c1,c2,c3,c4 = "blue","green","red","black"      # 各プロットの色
l1,l2,l3,l4 = "power","intensity","velocity","x"   # 各ラベル

ax1.plot(t_list, p_list, color=c1, label=l1)
ax2.plot(t_list, s_list, color=c2, label=l2)
ax3.plot(t_list, uionz_list, color=c3, label=l3)
ax4.plot(t_list, x_list, color=c4, label=l4)
ax1.legend(loc = 'upper right') #凡例
ax2.legend(loc = 'upper right') #凡例
ax3.legend(loc = 'upper right') #凡例
ax4.legend(loc = 'upper right') #凡例
ax1.set_xlim(0, 5)
ax1.set_ylim(0, 25)
ax2.set_xlim(0, 5)
ax2.set_ylim(0, 1000)
ax3.set_xlim(0, 5)
ax3.set_ylim(0, 8)
ax4.set_xlim(0, 5)
ax4.set_ylim(0, 16)
fig.tight_layout()  # レイアウトの設定
plt.show()
