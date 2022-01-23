import os
import sys
import warnings
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
from numba import jit
from PIL import Image, ImageDraw

from limitter import minmod

warnings.filterwarnings('ignore')


dt = 2e-10
nmax = 25000
theta_1 = 1  # 0~1, 0で陽解法、0.5でCrank Nicorson, 1で陰解法
theta_2 = 0.5  # 時間についての差分精度。0で1次, 0.5で2次

gamma = 1.4

m = 28.966
R_0 = 8314.0
R = R_0 / m
T_0 = 293

p_0 = 1.013e5
T_0 = 293
rho_0 = p_0 / R / T_0
u_0 = 0
v_0 = 0

PI = 1.013e5
RHOI = PI / R / T_0
UI = 0.0
VI = 0.0

PE = 1.013e4
RHOE = PE / R / T_0
UE = 0.0
VE = 0.0

jmax = 1000
kmax = 350
dx = 0.01e-3
dy = 0.01e-3
xmin, xmid, xmax = 0.0, 1.0, dx * (jmax - 1)
ymin, ymid, ymax = 0.0, 0.5, dy * (kmax - 1)
x = np.linspace(xmin, xmax, jmax)
y = np.linspace(ymin, ymax, kmax)
X, Y = np.meshgrid(x, y, indexing='ij')

# dx = (xmax - xmin) / (jmax - 1)
# dy = (ymax - ymin) / (kmax - 1)

dtdx = dt / dx
dtdy = dt / dy

slope = 0.16
intercept = 0.46
slope_low = 0.68
intercept_low = 0.30

A_G = 0.224517656
B_G = 0.77548
sigma_G1 = 0.84473  # mm
sigma_G2 = 1.75302  # mm
A_T = 0.764525962
B_T = 0.235474467
sigma_T1 = 1.557528296
sigma_T2 = 4.050355336
lmd = 10.6  # 単位はum
M2_G = 15
M2_T = 21
W_G0 = 1.7  # mm
W_T0 = 2.0
R_peak = 2.13
ss = 347  # 音速
l = 0.2e-3  # m 加熱長さ レーザーは0.2 mm
eta = 0.48


# ### Roeスキームによる計算
def init(shocktube=False):
    """
    Q の初期条件を決定
    :return Q:
    """
    Q = np.zeros([jmax, kmax, 4])
    if shocktube:
        Q[x <= xmid, :, 0] = RHOI
        Q[x <= xmid, :, 1] = RHOI * UI
        Q[x <= xmid, :, 2] = RHOI * VI
        Q[x <= xmid, :, 3] = (PI / (gamma - 1.0) + 0.5 * RHOI * (UI**2 + VI**2))

        Q[x > xmid, :, 0] = RHOE
        Q[x > xmid, :, 1] = RHOE * UE
        Q[x > xmid, :, 2] = RHOE * VE
        Q[x > xmid, :, 3] = (PE / (gamma - 1.0) + 0.5 * RHOE * (UE**2 + VE**2))
    else:
        Q[:, :, 0] = rho_0
        Q[:, :, 1] = 0
        Q[:, :, 2] = 0
        Q[:, :, 3] = (p_0 / (gamma - 1.0))
    return Q


def laser_intensity(nnum=nmax, deltat=dt):
    """
    :return S_laser: nmax * kmax
    :return x_laser: nmax * kmax
    :return u_ionz: nmax * kmax
    """
    response = {}
    x_laser = np.zeros([nnum, kmax])
    t = np.arange(0, nnum) * deltat
    # Power_laser = 8.15 * np.exp(-0.866 * t * 1e6)  # 松井
    Power_laser = 6.092231655 * np.exp(-0.7663806156 * t * 1e6)  # 菅村

    W_G = W_G0 * 1e-3
    W_T = W_T0 * 1e-3
    S_laser0 = R_peak * Power_laser / 4 / W_G / W_T * 1e-3  # 波頭のレーザー強度
    G_y = A_G * np.exp(-2 * (y * 1e3 / np.sqrt(2) / sigma_G1)**4) + B_G * np.exp(-2 * (y * 1e3 / np.sqrt(2) / sigma_G2)**2)
    T_y = A_T * np.exp(-2 * (y * 1e3 / np.sqrt(2) / sigma_T1)**4) + B_T * np.exp(-2 * (y * 1e3 / np.sqrt(2) / sigma_T2)**2)
    S_ratio = G_y * T_y  # S_laser / S_laser0
    S_laser = S_laser0[:, np.newaxis] * S_ratio[np.newaxis, :]  # 局所レーザー強度, GW/m2
    u_ionz = np.where((S_laser * 1e9 / rho_0 / ss**3)**(intercept - intercept_low) * slope / slope_low >= 1,
                      slope * (S_laser * 1e9 / rho_0 / ss**3) ** intercept * ss,
                      slope_low * (S_laser * 1e9 / rho_0 / ss**3) ** intercept_low * ss)
    for n in range(1, nnum):
        x_laser[n, 0] = x_laser[n - 1, 0] + u_ionz[n, 0] * deltat
        for k in range(1, kmax):
            # ルートの中身が正になっている間は進展位置を計算する。マイナスになるところには電離波面は存在しない
            if S_ratio[k] ** (-2 * intercept / (1 - intercept)) >= 1:
                # 切片はLine1のもので統一する仮定
                x_laser[n, k] = x_laser[n, k - 1] - np.sqrt(S_ratio[k] ** (-2 * intercept / (1 - intercept)) - 1) * dy  # m 電離波面の波頭進展位置=Laserによって加熱される最初の位置
            else:
                break
    response['S_laser'] = S_laser
    response['x_laser'] = x_laser
    response['u_ionz'] = u_ionz
    return response


def laser_propagation():
    """
    電離波面の進展とレーザー強度の減衰をgif化する
    """
    response = laser_intensity(nnum=50, deltat=1e-7)
    x_laser = response['x_laser']
    s_laser = response['S_laser']
    os.makedirs('propagation', exist_ok=True)

    fig = plt.figure(figsize=(8, 8))
    font_dict = dict(style="italic", size=16)
    for n in range(50):
        ax = fig.add_subplot(1, 2, 1)
        ax.plot(y * 1e3, x_laser[n] * 1e3)
        ax.set_ylim([0, 15])
        ax.set_xlim([0, 7.5])
        ax.text(4.0, 0.8, f"t={round(n*1e-1,1)} us", fontdict=font_dict)

        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(y * 1e3, s_laser[n])
        ax2.set_ylim([0, 1000])
        ax2.set_xlim([0, 5])
        ax2.text(4.0, 0.8, f"t={round(n*1e-1,1)} us", fontdict=font_dict)
        fig.savefig(f'propagation/{n}.png')
        plt.clf()
        plt.close()

    images = []
    for n in range(0, 50):
        im = Image.open(f'propagation/{n}.png')
        images.append(im)

    images[0].save('Propagation.gif', save_all=True, append_images=images[1:50], duration=0.1, loop=0)


def calc_CFL(Q):
    rho, rhou, rhov, e = Q[:, :, 0], Q[:, :, 1], Q[:, :, 2], Q[:, :, 3]

    u = rhou / rho
    v = rhov / rho
    p = (gamma - 1.0) * (e - 0.5 * rho * u ** 2)

    c = np.sqrt(gamma * p / rho)
    sp_u = c + np.abs(u)
    sp_v = c + np.abs(v)
    return round(np.max(sp_u) * dtdx, 3), round(np.max(sp_v) * dtdy, 3)


def Qtoq(Q):
    """
    保存量ベクトルQから基本量ベクトルqを形成する
    :param Q: 保存量ベクトルrho, rhou, e
    :return q: 基本量ベクトルrho, u, p
    """
    q = np.zeros([jmax, kmax, 4])
    rho, rhou, rhov, e = Q[:, :, 0], Q[:, :, 1], Q[:, :, 2], Q[:, :, 3]
    u = rhou / rho
    v = rhov / rho

    q[:, :, 0] = rho
    q[:, :, 1] = u
    q[:, :, 2] = v
    q[:, :, 3] = (gamma - 1.0) * (e - 0.5 * rho * (u**2 + v**2))  # p

    return q


def qtoQ(q):
    """
    基本量ベクトルqから保存量ベクトルQを形成する
    :param q: 基本量ベクトル
    :return Q: 保存量ベクトル
    """
    Q = np.zeros([jmax, kmax, 4])
    rho, u, v, p = q[:, :, 0], q[:, :, 1], q[:, :, 2], q[:, :, 3]
    e = p / (gamma - 1.0) + 0.5 * rho * (u**2 + v**2)

    Q[:, :, 0] = rho
    Q[:, :, 1] = rho * u
    Q[:, :, 2] = rho * v
    Q[:, :, 3] = e

    return Q


def MUSCL(Q, order: int = 2, kappa: float = 1 / 3, limit=minmod):
    """
    空間一次精度で求められた数値流速を代入することにより高次精度の数値流束を得る。
    保存変数ではなく基本変数で内挿しあとから戻した方が安定になりやすい。
    :param Q: 保存変数
    :param order: 精度のオーダー
    :param kappa:
    -1 : 2nd order fully upwind,
    0  : 2nd order upwind biased
    1  : 両側のセル平均値の代数平均となり、分布はセル境界で連続となる。
    1/3: 3rd order upwind biased
    :param limit: 制限関数。minmod, vanAlbada, superbeeなど
    minmodは安定性に長けているが、不連続は鈍る傾向にある。
    :return QL, QR: 内挿したQL_j+1/2,QR_j-1/2(基本変数)
    """
    # 基本変数で内挿する
    q = Qtoq(Q)

    qLx, qRx = q.copy(), q.copy()  # 1st order
    qLy, qRy = q.copy(), q.copy()  # 1st order

    if order == 2 or order == 3:
        # 2nd / 3rd order & minmod limitter
        dqx = np.zeros([jmax, kmax, 4])
        dqy = np.zeros([jmax, kmax, 4])
        dqx[:-1, :] = q[1:, :] - q[:-1, :]
        dqy[:, :-1] = q[:, 1:] - q[:, :-1]

        b = (3.0 - kappa) / (1.0 - kappa)  # 式(2.74)

        # 制限関数の導入によりTVD化
        for k in range(1, kmax - 1):
            for j in range(1, jmax - 1):
                # x方向
                Dpx = limit(dqx[j, k], b * dqx[j - 1, k])  # 式(2.73a)
                Dmx = limit(dqx[j - 1, k], b * dqx[j, k])     # 式(2.73b)

                qLx[j, k] += 0.25 * ((1.0 + kappa) * Dpx + (1.0 - kappa) * Dmx)  # QL_j+1/2, k
                qRx[j, k] -= 0.25 * ((1.0 - kappa) * Dpx + (1.0 + kappa) * Dmx)  # QR_j-1/2 , k

                # y方向
                Dpy = limit(dqy[j, k], b * dqy[j, k - 1])  # 式(2.73a)
                Dmy = limit(dqy[j, k - 1], b * dqy[j, k])     # 式(2.73b)

                qLy[j, k] += 0.25 * ((1.0 + kappa) * Dpy + (1.0 - kappa) * Dmy)  # QL_j, k+1/2
                qRy[j, k] -= 0.25 * ((1.0 - kappa) * Dpy + (1.0 + kappa) * Dmy)  # QR_j, k-1/2
        # 境界
        # qL[0] = qL[1]
        # qR[-1] = qR[-2]

    return qLx, qRx, qLy, qRy


def MUSCL_x(Q, order: int, kappa: float, limit=minmod):
    """
    空間一次精度で求められた数値流速を代入することにより高次精度の数値流束を得る。
    保存変数ではなく基本変数で内挿しあとから戻した方が安定になりやすい。
    :param Q: 保存変数
    :param order: 精度のオーダー
    :param kappa:
    -1 : 2nd order fully upwind,
    0  : 2nd order upwind biased
    1  : 両側のセル平均値の代数平均となり、分布はセル境界で連続となる。
    1/3: 3rd order upwind biased
    :param limit: 制限関数。minmod, vanAlbada, superbeeなど
    minmodは安定性に長けているが、不連続は鈍る傾向にある。
    :return QL, QR: 内挿したQL_j+1/2,QR_j-1/2(基本変数)
    """
    # 基本変数で内挿する
    q = Qtoq(Q)

    qLx, qRx = q.copy(), q.copy()  # 1st order

    if order == 2 or order == 3:
        # 2nd / 3rd order & minmod limitter
        dqx = np.zeros([jmax, kmax, 4])
        dqx[:-1, :] = q[1:, :] - q[:-1, :]

        b = 2  # (3.0 - kappa) / (1.0 - kappa)  # 式(2.74)

        # 制限関数の導入によりTVD化
        for k in range(1, kmax - 1):
            for j in range(1, jmax - 1):
                # x方向
                Dpx = limit(dqx[j, k], b * dqx[j - 1, k])  # 式(2.73a)
                Dmx = limit(dqx[j - 1, k], b * dqx[j, k])     # 式(2.73b)

                qLx[j, k] += 0.25 * ((1.0 + kappa) * Dpx + (1.0 - kappa) * Dmx)  # QL_j+1/2, k
                qRx[j, k] -= 0.25 * ((1.0 - kappa) * Dpx + (1.0 + kappa) * Dmx)  # QR_j-1/2 , k

        # 境界
        # qL[0] = qL[1]
        # qR[-1] = qR[-2]

    return qLx, qRx


def MUSCL_y(Q, order: int, kappa: float, limit=minmod):
    """
    空間一次精度で求められた数値流速を代入することにより高次精度の数値流束を得る。
    保存変数ではなく基本変数で内挿しあとから戻した方が安定になりやすい。
    :param Q: 保存変数
    :param order: 精度のオーダー
    :param kappa:
    -1 : 2nd order fully upwind,
    0  : 2nd order upwind biased
    1  : 両側のセル平均値の代数平均となり、分布はセル境界で連続となる。
    1/3: 3rd order upwind biased
    :param limit: 制限関数。minmod, vanAlbada, superbeeなど
    minmodは安定性に長けているが、不連続は鈍る傾向にある。
    :return QL, QR: 内挿したQL_j+1/2,QR_j-1/2(基本変数)
    """
    # 基本変数で内挿する
    q = Qtoq(Q)

    qLy, qRy = q.copy(), q.copy()  # 1st order

    if order == 2 or order == 3:
        # 2nd / 3rd order & minmod limitter
        dqy = np.zeros([jmax, kmax, 4])
        dqy[:, :-1] = q[:, 1:] - q[:, :-1]

        b = 2  # (3.0 - kappa) / (1.0 - kappa)  # 式(2.74)

        # 制限関数の導入によりTVD化
        for k in range(1, kmax - 1):
            for j in range(1, jmax - 1):
                # y方向
                Dpy = limit(dqy[j, k], b * dqy[j, k - 1])  # 式(2.73a)
                Dmy = limit(dqy[j, k - 1], b * dqy[j, k])     # 式(2.73b)

                qLy[j, k] += 0.25 * ((1.0 + kappa) * Dpy + (1.0 - kappa) * Dmy)  # QL_j, k+1/2
                qRy[j, k] -= 0.25 * ((1.0 - kappa) * Dpy + (1.0 + kappa) * Dmy)  # QR_j, k-1/2
        # 境界
        # qL[0] = qL[1]
        # qR[-1] = qR[-2]

    return qLy, qRy


def Roe_flux(qLx, qRx, qLy, qRy, E, F):
    """
    van LeerのMUSTL型FDS法(高次精度)
    Roe平均を用いて数値流束E_j+1/2を求める
    :param qL, qR: 内挿する基本変数、qLはj_+1/2, qRはj-1/2
    :param E: 更新するE_j+1/2
    """
    dQ, EL, ER, FL, FR, AQ, BQ = np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4])
    Lambda, R, Rinv = np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4])

    # x方向
    rhoL, uL, vL, pL = qLx[:-1, :-1, 0], qLx[:-1, :-1, 1], qLx[:-1, :-1, 2], qLx[:-1, :-1, 3]
    rhoR, uR, vR, pR = qRx[1:, :-1, 0], qRx[1:, :-1, 1], qRx[1:, :-1, 2], qRx[1:, :-1, 3]

    eL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL**2 + vL**2)
    eR = pR / (gamma - 1.0) + 0.5 * rhoR * (uR**2 + vR**2)

    HL = (eL + pL) / rhoL
    HR = (eR + pR) / rhoR

    # Roe平均 式(6.38)
    sqrhoL = np.sqrt(rhoL)
    sqrhoR = np.sqrt(rhoR)

    uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
    vAVE = (sqrhoL * vL + sqrhoR * vR) / (sqrhoL + sqrhoR)
    qAVE = np.sqrt(uAVE**2 + vAVE**2)
    HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
    cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * (uAVE**2 + vAVE**2)))

    # # デバッグ用
    # if np.sum(np.isnan(cAVE))>0:
    #     print('cAVEx has nan')
    #     np.savetxt('error/rhoLx.csv', rhoL, fmt='%.2e')
    #     np.savetxt('error/rhoRx.csv', rhoR, fmt='%.2e')
    #     np.savetxt('error/HAVEx.csv', HAVE, fmt='%.2e')
    #     np.savetxt('error/cAVEx.csv', cAVE, fmt='%.2e')

    #     fig = plt.figure(figsize=(10, 8), dpi=100)  # グラフのサイズ
    #     plt.rcParams["font.size"] = 10  # グラフの文字サイズ

    #     ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    #     ax1.plot_surface(X[:-1, :-1] * 1e3, Y[:-1, :-1] * 1e3, rhoL, color='red', linewidth=1.5, label='Numerical')
    #     ax1.set_xlabel(r'$z$ mm')
    #     ax1.set_ylabel(r'$r$ mm')
    #     ax1.set_zlabel(r'$rhoL$')

    #     plt.show()
    #     sys.exit()

    b1 = 0.5 * (gamma - 1.0) * qAVE ** 2 / cAVE ** 2
    b2 = (gamma - 1.0) / cAVE ** 2

    Lambda[:, :, 0, 0] = np.abs(uAVE - cAVE)
    Lambda[:, :, 1, 1] = np.abs(uAVE)
    Lambda[:, :, 2, 2] = np.abs(uAVE + cAVE)
    Lambda[:, :, 3, 3] = np.abs(uAVE)

    dQ[:, :, 0] = rhoR - rhoL
    dQ[:, :, 1] = rhoR * uR - rhoL * uL
    dQ[:, :, 2] = rhoR * vR - rhoL * vL
    dQ[:, :, 3] = eR - eL

    EL[:, :, 0] = rhoL * uL
    EL[:, :, 1] = pL + rhoL * uL**2
    EL[:, :, 2] = rhoL * uL * vL
    EL[:, :, 3] = (eL + pL) * uL
    ER[:, :, 0] = rhoR * uR
    ER[:, :, 1] = pR + rhoR * uR**2
    ER[:, :, 2] = rhoR * uR * vR
    ER[:, :, 3] = (eR + pR) * uR

    R[:, :, 0, 0] = 1.0
    R[:, :, 0, 1] = 1.0
    R[:, :, 0, 2] = 1.0
    R[:, :, 1, 0] = uAVE - cAVE
    R[:, :, 1, 1] = uAVE
    R[:, :, 1, 2] = uAVE + cAVE
    R[:, :, 2, 0] = vAVE
    R[:, :, 2, 1] = vAVE
    R[:, :, 2, 2] = vAVE
    R[:, :, 2, 3] = 1.0
    R[:, :, 3, 0] = HAVE - uAVE * cAVE
    R[:, :, 3, 1] = 0.5 * qAVE**2
    R[:, :, 3, 2] = HAVE + uAVE * cAVE
    R[:, :, 3, 3] = vAVE

    Rinv[:, :, 0, 0] = 0.5 * (b1 + uAVE / cAVE)
    Rinv[:, :, 0, 1] = -0.5 * (b2 * uAVE + 1 / cAVE)
    Rinv[:, :, 0, 2] = -0.5 * b2 * vAVE
    Rinv[:, :, 0, 3] = 0.5 * b2
    Rinv[:, :, 1, 0] = 1.0 - b1
    Rinv[:, :, 1, 1] = b2 * uAVE
    Rinv[:, :, 1, 2] = b2 * vAVE
    Rinv[:, :, 1, 3] = -b2
    Rinv[:, :, 2, 0] = 0.5 * (b1 - uAVE / cAVE)
    Rinv[:, :, 2, 1] = -0.5 * (b2 * uAVE - 1 / cAVE)
    Rinv[:, :, 2, 2] = -0.5 * b2 * vAVE
    Rinv[:, :, 2, 3] = 0.5 * b2
    Rinv[:, :, 3, 0] = -vAVE
    Rinv[:, :, 3, 2] = 1.0

    for j in range(jmax - 1):
        for k in range(kmax - 1):
            AQ[j, k] = R[j, k] @ Lambda[j, k] @ Rinv[j, k] @ dQ[j, k]

    E[:-1, :-1] = 0.5 * (ER + EL - AQ)  # 式(6.43)

    # y方向
    rhoL, uL, vL, pL = qLy[-1, :-1, 0], qLy[-1, :-1, 1], qLy[-1, :-1, 2], qLy[-1, :-1, 3]
    rhoR, uR, vR, pR = qRy[:-1, 1:, 0], qRy[:-1, 1:, 1], qRy[:-1, 1:, 2], qRy[:-1, 1:, 3]

    eL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL**2 + vL**2)
    eR = pR / (gamma - 1.0) + 0.5 * rhoR * (uR**2 + vR**2)

    HL = (eL + pL) / rhoL
    HR = (eR + pR) / rhoR

    # Roe平均 式(6.38)
    sqrhoL = np.sqrt(rhoL)
    sqrhoR = np.sqrt(rhoR)

    uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
    vAVE = (sqrhoL * vL + sqrhoR * vR) / (sqrhoL + sqrhoR)
    qAVE = np.sqrt(uAVE**2 + vAVE**2)
    HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
    cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * (uAVE**2 + vAVE**2)))
    # デバッグ用
    # if np.sum(np.isnan(cAVE)) > 0:
    #     print('cAVEy has nan')
    #     np.savetxt('error/rhoLy.csv', rhoL, fmt='%.2e')
    #     np.savetxt('error/rhoRy.csv', rhoR, fmt='%.2e')

    #     fig = plt.figure(figsize=(10, 8), dpi=100)  # グラフのサイズ
    #     plt.rcParams["font.size"] = 10  # グラフの文字サイズ

    #     ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    #     ax1.plot_surface(X[:-1, :-1] * 1e3, Y[:-1, :-1] * 1e3, rhoL, color='red', linewidth=1.5, label='Numerical')
    #     ax1.set_xlabel(r'$z$ mm')
    #     ax1.set_ylabel(r'$r$ mm')
    #     ax1.set_zlabel(r'$rhoL$')

    #     plt.show()
    #     sys.exit()
    b1 = 0.5 * (gamma - 1.0) * qAVE ** 2 / cAVE ** 2
    b2 = (gamma - 1.0) / cAVE ** 2

    Lambda[:, :, 0, 0] = np.abs(vAVE - cAVE)
    Lambda[:, :, 1, 1] = np.abs(vAVE)
    Lambda[:, :, 2, 2] = np.abs(vAVE + cAVE)
    Lambda[:, :, 3, 3] = np.abs(vAVE)

    dQ[:, :, 0] = rhoR - rhoL
    dQ[:, :, 1] = rhoR * uR - rhoL * uL
    dQ[:, :, 2] = rhoR * vR - rhoL * vL
    dQ[:, :, 3] = eR - eL

    FL[:, :, 0] = rhoL * vL
    FL[:, :, 1] = rhoL * uL * vL
    FL[:, :, 2] = pL + rhoL * vL**2
    FL[:, :, 3] = (eL + pL) * vL
    FR[:, :, 0] = rhoR * vR
    FR[:, :, 1] = rhoR * uR * vR
    FR[:, :, 2] = pR + rhoR * vR**2
    FR[:, :, 3] = (eR + pR) * vR

    R[:, :, 0, 0] = 1.0
    R[:, :, 0, 1] = 1.0
    R[:, :, 0, 2] = 1.0
    R[:, :, 1, 0] = uAVE
    R[:, :, 1, 1] = uAVE
    R[:, :, 1, 2] = uAVE
    R[:, :, 1, 3] = 1.0
    R[:, :, 2, 0] = vAVE - cAVE
    R[:, :, 2, 1] = vAVE
    R[:, :, 2, 2] = vAVE + cAVE
    R[:, :, 3, 0] = HAVE - vAVE * cAVE
    R[:, :, 3, 1] = 0.5 * qAVE**2
    R[:, :, 3, 2] = HAVE + vAVE * cAVE
    R[:, :, 3, 3] = uAVE

    Rinv[:, :, 0, 0] = 0.5 * (b1 + vAVE / cAVE)
    Rinv[:, :, 0, 1] = -0.5 * b2 * uAVE
    Rinv[:, :, 0, 2] = -0.5 * (b2 * vAVE + 1 / cAVE)
    Rinv[:, :, 0, 3] = 0.5 * b2
    Rinv[:, :, 1, 0] = 1.0 - b1
    Rinv[:, :, 1, 1] = b2 * uAVE
    Rinv[:, :, 1, 2] = b2 * vAVE
    Rinv[:, :, 1, 3] = -b2
    Rinv[:, :, 2, 0] = 0.5 * (b1 - vAVE / cAVE)
    Rinv[:, :, 2, 1] = -0.5 * b2 * uAVE
    Rinv[:, :, 2, 2] = -0.5 * (b2 * vAVE - 1 / cAVE)
    Rinv[:, :, 2, 3] = 0.5 * b2
    Rinv[:, :, 3, 0] = -uAVE
    Rinv[:, :, 3, 1] = 1.0

    for j in range(jmax - 1):
        for k in range(kmax - 1):
            BQ[j, k] = R[j, k] @ Lambda[j, k] @ Rinv[j, k] @ dQ[j, k]

    F[:-1, :-1] = 0.5 * (FR + FL - BQ)  # 式(6.43)


def Roe_flux_x(qLx, qRx, E):
    """
    van LeerのMUSTL型FDS法(高次精度)
    Roe平均を用いて数値流束E_j+1/2を求める
    :param qL, qR: 内挿する基本変数、qLはj_+1/2, qRはj-1/2
    :param E: 更新するE_j+1/2
    """
    dQ, EL, ER, AQ = np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4])
    Lambda, R, Rinv = np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4])

    # x方向
    rhoL, uL, vL, pL = qLx[:-1, :-1, 0], qLx[:-1, :-1, 1], qLx[:-1, :-1, 2], qLx[:-1, :-1, 3]
    rhoR, uR, vR, pR = qRx[1:, :-1, 0], qRx[1:, :-1, 1], qRx[1:, :-1, 2], qRx[1:, :-1, 3]

    eL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL**2 + vL**2)
    eR = pR / (gamma - 1.0) + 0.5 * rhoR * (uR**2 + vR**2)

    HL = (eL + pL) / rhoL
    HR = (eR + pR) / rhoR

    # Roe平均 式(6.38)
    sqrhoL = np.sqrt(rhoL)
    sqrhoR = np.sqrt(rhoR)

    uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
    vAVE = (sqrhoL * vL + sqrhoR * vR) / (sqrhoL + sqrhoR)
    qAVE = np.sqrt(uAVE**2 + vAVE**2)
    HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
    cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * (uAVE**2 + vAVE**2)))

    b1 = 0.5 * (gamma - 1.0) * qAVE ** 2 / cAVE ** 2
    b2 = (gamma - 1.0) / cAVE ** 2

    Lambda[:, :, 0, 0] = np.abs(uAVE - cAVE)
    Lambda[:, :, 1, 1] = np.abs(uAVE)
    Lambda[:, :, 2, 2] = np.abs(uAVE + cAVE)
    Lambda[:, :, 3, 3] = np.abs(uAVE)

    dQ[:, :, 0] = rhoR - rhoL
    dQ[:, :, 1] = rhoR * uR - rhoL * uL
    dQ[:, :, 2] = rhoR * vR - rhoL * vL
    dQ[:, :, 3] = eR - eL

    EL[:, :, 0] = rhoL * uL
    EL[:, :, 1] = pL + rhoL * uL**2
    EL[:, :, 2] = rhoL * uL * vL
    EL[:, :, 3] = (eL + pL) * uL
    ER[:, :, 0] = rhoR * uR
    ER[:, :, 1] = pR + rhoR * uR**2
    ER[:, :, 2] = rhoR * uR * vR
    ER[:, :, 3] = (eR + pR) * uR

    R[:, :, 0, 0] = 1.0
    R[:, :, 0, 1] = 1.0
    R[:, :, 0, 2] = 1.0
    R[:, :, 1, 0] = uAVE - cAVE
    R[:, :, 1, 1] = uAVE
    R[:, :, 1, 2] = uAVE + cAVE
    R[:, :, 2, 0] = vAVE
    R[:, :, 2, 1] = vAVE
    R[:, :, 2, 2] = vAVE
    R[:, :, 2, 3] = 1.0
    R[:, :, 3, 0] = HAVE - uAVE * cAVE
    R[:, :, 3, 1] = 0.5 * qAVE**2
    R[:, :, 3, 2] = HAVE + uAVE * cAVE
    R[:, :, 3, 3] = vAVE

    Rinv[:, :, 0, 0] = 0.5 * (b1 + uAVE / cAVE)
    Rinv[:, :, 0, 1] = -0.5 * (b2 * uAVE + 1 / cAVE)
    Rinv[:, :, 0, 2] = -0.5 * b2 * vAVE
    Rinv[:, :, 0, 3] = 0.5 * b2
    Rinv[:, :, 1, 0] = 1.0 - b1
    Rinv[:, :, 1, 1] = b2 * uAVE
    Rinv[:, :, 1, 2] = b2 * vAVE
    Rinv[:, :, 1, 3] = -b2
    Rinv[:, :, 2, 0] = 0.5 * (b1 - uAVE / cAVE)
    Rinv[:, :, 2, 1] = -0.5 * (b2 * uAVE - 1 / cAVE)
    Rinv[:, :, 2, 2] = -0.5 * b2 * vAVE
    Rinv[:, :, 2, 3] = 0.5 * b2
    Rinv[:, :, 3, 0] = -vAVE
    Rinv[:, :, 3, 2] = 1.0

    for j in range(jmax - 1):
        for k in range(kmax - 1):
            AQ[j, k] = R[j, k] @ Lambda[j, k] @ Rinv[j, k] @ dQ[j, k]

    E[:-1, :-1] = 0.5 * (ER + EL - AQ)  # 式(6.43)


def Roe_flux_y(qLy, qRy, F):
    """
    van LeerのMUSTL型FDS法(高次精度)
    Roe平均を用いて数値流束E_j+1/2を求める
    :param qL, qR: 内挿する基本変数、qLはj_+1/2, qRはj-1/2
    :param E: 更新するE_j+1/2
    """
    dQ, FL, FR, BQ = np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4]), np.zeros([jmax - 1, kmax - 1, 4])
    Lambda, R, Rinv = np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4]), np.zeros([jmax - 1, kmax - 1, 4, 4])

    # y方向
    rhoL, uL, vL, pL = qLy[-1, :-1, 0], qLy[-1, :-1, 1], qLy[-1, :-1, 2], qLy[-1, :-1, 3]
    rhoR, uR, vR, pR = qRy[:-1, 1:, 0], qRy[:-1, 1:, 1], qRy[:-1, 1:, 2], qRy[:-1, 1:, 3]

    eL = pL / (gamma - 1.0) + 0.5 * rhoL * (uL**2 + vL**2)
    eR = pR / (gamma - 1.0) + 0.5 * rhoR * (uR**2 + vR**2)

    HL = (eL + pL) / rhoL
    HR = (eR + pR) / rhoR

    # Roe平均 式(6.38)
    sqrhoL = np.sqrt(rhoL)
    sqrhoR = np.sqrt(rhoR)

    uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
    vAVE = (sqrhoL * vL + sqrhoR * vR) / (sqrhoL + sqrhoR)
    qAVE = np.sqrt(uAVE**2 + vAVE**2)
    HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
    cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * (uAVE**2 + vAVE**2)))

    b1 = 0.5 * (gamma - 1.0) * qAVE ** 2 / cAVE ** 2
    b2 = (gamma - 1.0) / cAVE ** 2

    Lambda[:, :, 0, 0] = np.abs(vAVE - cAVE)
    Lambda[:, :, 1, 1] = np.abs(vAVE)
    Lambda[:, :, 2, 2] = np.abs(vAVE + cAVE)
    Lambda[:, :, 3, 3] = np.abs(vAVE)

    dQ[:, :, 0] = rhoR - rhoL
    dQ[:, :, 1] = rhoR * uR - rhoL * uL
    dQ[:, :, 2] = rhoR * vR - rhoL * vL
    dQ[:, :, 3] = eR - eL

    FL[:, :, 0] = rhoL * vL
    FL[:, :, 1] = rhoL * uL * vL
    FL[:, :, 2] = pL + rhoL * vL**2
    FL[:, :, 3] = (eL + pL) * vL
    FR[:, :, 0] = rhoR * vR
    FR[:, :, 1] = rhoR * uR * vR
    FR[:, :, 2] = pR + rhoR * vR**2
    FR[:, :, 3] = (eR + pR) * vR

    R[:, :, 0, 0] = 1.0
    R[:, :, 0, 1] = 1.0
    R[:, :, 0, 2] = 1.0
    R[:, :, 1, 0] = uAVE
    R[:, :, 1, 1] = uAVE
    R[:, :, 1, 2] = uAVE
    R[:, :, 1, 3] = 1.0
    R[:, :, 2, 0] = vAVE - cAVE
    R[:, :, 2, 1] = vAVE
    R[:, :, 2, 2] = vAVE + cAVE
    R[:, :, 3, 0] = HAVE - vAVE * cAVE
    R[:, :, 3, 1] = 0.5 * qAVE**2
    R[:, :, 3, 2] = HAVE + vAVE * cAVE
    R[:, :, 3, 3] = uAVE

    Rinv[:, :, 0, 0] = 0.5 * (b1 + vAVE / cAVE)
    Rinv[:, :, 0, 1] = -0.5 * b2 * uAVE
    Rinv[:, :, 0, 2] = -0.5 * (b2 * vAVE + 1 / cAVE)
    Rinv[:, :, 0, 3] = 0.5 * b2
    Rinv[:, :, 1, 0] = 1.0 - b1
    Rinv[:, :, 1, 1] = b2 * uAVE
    Rinv[:, :, 1, 2] = b2 * vAVE
    Rinv[:, :, 1, 3] = -b2
    Rinv[:, :, 2, 0] = 0.5 * (b1 - vAVE / cAVE)
    Rinv[:, :, 2, 1] = -0.5 * b2 * uAVE
    Rinv[:, :, 2, 2] = -0.5 * (b2 * vAVE - 1 / cAVE)
    Rinv[:, :, 2, 3] = 0.5 * b2
    Rinv[:, :, 3, 0] = -uAVE
    Rinv[:, :, 3, 1] = 1.0

    for j in range(jmax - 1):
        for k in range(kmax - 1):
            BQ[j, k] = R[j, k] @ Lambda[j, k] @ Rinv[j, k] @ dQ[j, k]

    F[:-1, :-1] = 0.5 * (FR + FL - BQ)  # 式(6.43)


# 陽解法
def explicit_fractional_timestep(Q, order, kappa, nmax, print_interval=8):
    """
    多次元問題の時間分割法
    """
    E = np.zeros([jmax, kmax, 4])
    F = np.zeros([jmax, kmax, 4])
    H = np.zeros([jmax, kmax, 4])
    Plateau = np.zeros([nmax, kmax])
    response = laser_intensity(nnum=2 * nmax, deltat=dt/2)
    x_laser = response['x_laser']
    S_laser = response['S_laser']
    for n in range(nmax):
        S = np.zeros([jmax, kmax, 4])
        S_ = np.zeros([jmax, kmax, 4])
        for k in range(kmax):
            S[(x > x_laser[2 * n, k] - l) & (x < x_laser[2 * n, k]), k, 3] = eta * S_laser[2 * n, k] / l * 1e9  # レーザー熱入力w
            S_[(x > x_laser[2 * n + 1, k] - l) & (x < x_laser[2 * n + 1, k]), k, 3] = eta * S_laser[2 * n + 1, k] / l * 1e9  # レーザー熱入力w

        qLx, qRx = MUSCL_x(Q, order, kappa)  # 保存変数→基本変数
        Roe_flux_x(qLx, qRx, E)  # 保存変数に戻り、Eをアップデート
        Q[1:, 1:] -= 0.5 * dtdx * (E[1:, 1:] - E[:-1, 1:]) + 0.5 * dt * (H[1:, 1:] / Y[1:, 1:, np.newaxis] - S[1:, 1:])  # 前進差分
        boundary_condition(Q)

        qLy, qRy = MUSCL_y(Q, order, kappa)  # 保存変数→基本変数
        Roe_flux_y(qLy, qRy, F)  # 保存変数に戻り、Eをアップデート
        Q[1:, 1:] -= dtdy * (F[1:, 1:] - F[1:, :-1])  # 前進差分
        boundary_condition(Q)

        # qLx, qRx, qLy, qRy = MUSCL(Q, order, kappa)  # 保存変数→基本変数
        # Roe_flux(qLx, qRx, qLy, qRy, E, F)  # 保存変数に戻り、Eをアップデート
        # Q[1:, 1:] -= 0.5 * dtdy * (F[1:, 1:] - F[1:, :-1])  # 前進差分
        # boundary_condition(Q)

        qLx, qRx = MUSCL_x(Q, order, kappa)  # 保存変数→基本変数
        Roe_flux_x(qLx, qRx, E)  # 保存変数に戻り、Eをアップデート
        Q[1:, 1:] -= 0.5 * dtdx * (E[1:, 1:] - E[:-1, 1:]) + 0.5 * dt * (H[1:, 1:] / Y[1:, 1:, np.newaxis] - S_[1:, 1:])  # 前進差分
        boundary_condition(Q)

        rho = Q[:, :, 0]
        u = Q[:, :, 1] / Q[:, :, 0]
        v = Q[:, :, 2] / Q[:, :, 0]
        e = Q[:, :, 3]
        p = (gamma - 1.0) * (e - 0.5 * rho * (u**2 + v**2))  # / 1e5
        Plateau[n, :] = p[0, :]

        if n % print_interval == 0:
            print(f'n= {n}', round((n + 1) / nmax, 2) * 100, '%')
            print('CFL:', calc_CFL(Q))
            np.savetxt(f'2Daxi/rho/{n}.csv', rho, fmt='%.2e')
            np.savetxt(f'2Daxi/u/{n}.csv', u, fmt='%.2e')
            np.savetxt(f'2Daxi/v/{n}.csv', v, fmt='%.2e')
            np.savetxt(f'2Daxi/p/{n}.csv', p, fmt='%.2e')
            np.savetxt(f'2Daxi/e/{n}.csv', e, fmt='%.2e')
            output_graphs(rho, u, v, p, n)
    response['Plateau'] = Plateau
    return response


def explicit_runge(Q, order, kappa, nmax, print_interval=2):
    """
    4段階ルンゲクッタ
    """
    RHS = np.zeros([jmax, kmax, 4])
    E = np.zeros([jmax, kmax, 4])
    F = np.zeros([jmax, kmax, 4])
    H = np.zeros([jmax, kmax, 4])
    response = laser_intensity()
    x_laser = response['x_laser']
    S_laser = response['S_laser']

    for n in range(nmax):
        if n % print_interval == 0:
            print(round((n + 1) / nmax, 2) * 100, '%')
            print('CFL:', calc_CFL(Q))

        Qk = Q.copy()
        Qn = Q.copy()
        S = np.zeros([jmax, kmax, 4])
        for k in range(kmax):
            S[(x > x_laser[n, k] - l) & (x < x_laser[n, k]), k, 3] = eta * S_laser[n, k] / l * 1e9  # レーザー熱入力w

        coefs = [[0.5, 1/6], [0.5, 1/3], [1.0, 1/3], [0, 1/6]]
        for coef in coefs:
            qLx, qRx, qLy, qRy = MUSCL(Qk, order, kappa)  # 保存変数→基本変数
            Roe_flux(qLx, qRx, qLy, qRy, E, F)  # 保存変数に戻り、Eをアップデート
            RHS[1:, 1:] = dtdx * (E[1:, 1:] - E[:-1, 1:]) + dtdy * (F[1:, 1:] - F[1:, :-1]) + dt * (H[1:, 1:] / Y[1:, 1:, np.newaxis] - S[1:, 1:])  # 前進差分
            Qk[1:, 1:] = Qn[1:, 1:] - coef[0] * RHS[1:, 1:]
            Q[1:, 1:] -= dt * coef[1] * RHS[1:, 1:]

            boundary_condition(Qk)
        boundary_condition(Q)
    return response


def calc_jacobian_matrix(Q):
    """
    A = dE/dQ, B = dF/dQ
    近似LDU分解
    """
    # 基本変数の用意
    rho, rhou, rhov, e = Q[:, :, 0], Q[:, :, 1], Q[:, :, 2], Q[:, :, 3]
    u = rhou / rho
    v = rhov / rho
    q = np.sqrt(u**2 + v**2)
    p = (gamma - 1) * (e - 0.5 * rho * q**2)
    # H = (e + p) / rho
    c = np.sqrt(gamma * p / rho)
    phi = (gamma - 1.0) * 0.5 * (u**2 + v**2)
    omega = gamma * e / rho

    A = np.zeros([jmax, kmax, 4, 4])

    A[:, :, 0, 1] = 1
    A[:, :, 1, 0] = -(3 - gamma) / 2 * u**2 + (gamma - 1) / 2 * v**2
    A[:, :, 1, 1] = (3 - gamma) * u
    A[:, :, 1, 2] = -(gamma - 1) * v
    A[:, :, 1, 3] = gamma - 1
    A[:, :, 2, 0] = -u * v
    A[:, :, 2, 1] = v
    A[:, :, 2, 3] = u
    A[:, :, 3, 0] = -gamma * u * e / rho + (gamma - 1) * u * q**2
    A[:, :, 3, 1] = gamma * e / rho + (gamma - 1) / 2 * (2 * u**2 + q**2)
    A[:, :, 3, 2] = -(gamma - 1) * u * v
    A[:, :, 3, 3] = gamma * u

    B = np.zeros([jmax, kmax, 4, 4])

    B[:, :, 0, 2] = 1
    B[:, :, 1, 0] = -u * v
    B[:, :, 1, 1] = v
    B[:, :, 1, 2] = u
    B[:, :, 2, 0] = -(3 - gamma) / 2 * v**2 + (gamma - 1) / 2 * u**2
    B[:, :, 2, 1] = -(gamma - 1) * u
    B[:, :, 2, 2] = (3 - gamma) * v
    B[:, :, 2, 3] = gamma - 1
    B[:, :, 3, 0] = -gamma * v * e / rho + (gamma - 1) * v * q**2
    B[:, :, 3, 1] = -(gamma - 1) * u * v
    B[:, :, 3, 2] = gamma * e / rho + (gamma - 1) / 2 * (2 * v**2 + q**2)
    B[:, :, 3, 3] = gamma * v

    C = np.zeros([jmax, kmax, 4, 4])

    C[:, :, 0, 2] = 1
    C[:, :, 1, 0] = - u * v
    C[:, :, 1, 1] = v
    C[:, :, 1, 2] = u
    C[:, :, 2, 0] = -v**2
    C[:, :, 2, 2] = 2 * v
    C[:, :, 3, 0] = v * (2 * phi**2 - omega)
    C[:, :, 3, 1] = -(gamma - 1) * u * v
    C[:, :, 3, 2] = omega - phi**2 - (gamma - 1) * v**2
    C[:, :, 3, 3] = gamma * v

    H = np.zeros([jmax, kmax, 4])

    H[:, :, 0] = rho * v
    H[:, :, 1] = rho * u * v
    H[:, :, 2] = rho * v**2
    H[:, :, 3] = (e + p) * v

    sigma_x = abs(u) + c  # 一番大きい値で近似する
    sigma_y = abs(v) + c  # 一番大きい値で近似する

    # A+
    Ap = A.copy()
    Ap[:, :, 0, 0] += sigma_x
    Ap[:, :, 1, 1] += sigma_x
    Ap[:, :, 2, 2] += sigma_x
    Ap[:, :, 3, 3] += sigma_x
    Ap = Ap / 2

    # A-
    Am = A.copy()
    Am[:, :, 0, 0] -= sigma_x
    Am[:, :, 1, 1] -= sigma_x
    Am[:, :, 2, 2] -= sigma_x
    Am[:, :, 3, 3] -= sigma_x
    Am = Am / 2

    # B+
    Bp = B.copy()
    Bp[:, :, 0, 0] += sigma_y
    Bp[:, :, 1, 1] += sigma_y
    Bp[:, :, 2, 2] += sigma_y
    Bp[:, :, 3, 3] += sigma_y
    Bp = Bp / 2

    # B-
    Bm = B.copy()
    Bm[:, :, 0, 0] -= sigma_y
    Bm[:, :, 1, 1] -= sigma_y
    Bm[:, :, 2, 2] -= sigma_y
    Bm[:, :, 3, 3] -= sigma_y
    Bm = Bm / 2

    return Ap, Am, sigma_x, Bp, Bm, sigma_y, C, H


def boundary_condition(Q):
    """
    境界条件を反映する
    :param Q:
    """
    # アルミ板
    Q[0, :] = Q[1, :]
    Q[0, :, 1] = 0
    Q[0, :, 2] = 0
    # 回転軸
    Q[:, 0] = Q[:, 1]
    Q[:, 0, 2] = 0  # 半径方向流出は0のはず
    # 半径方向流出
    Q[-1, :] = Q[-2, :]
    # 光軸方向流出
    Q[:, -1] = Q[:, -2]

    # 角
    Q[-1, 0] = (Q[-1, 1] + Q[-2, 0]) / 2
    Q[-1, -1] = (Q[-1, -2] + Q[-2, -1]) / 2


# 陰解法
# @jit(cache=True)
def implicit_solution(Q, order, kappa, nmax: int, norm_limit: float = 1e-6, iimax: int = 10):
    """
    陰解法では、かなりの近似を導入しているため、時間精度を高めるためには工夫が必要。
    1. 時間微分項を3点の後退差分で近似
    2. 内部反復を利用
    :param Q: 初期条件を反映済みのQ
    :param order: 精度の次数
    :param kappa: muscl法の精度と手法に影響
    :param nmax: time step num
    :param norm_limit: ノルムの上限。内部反復の収束判定に使う
    :param iimax: 内部反復を最大行う回数
    """
    RHS = np.zeros([jmax, kmax, 4])
    E = np.zeros([jmax, kmax, 4])
    F = np.zeros([jmax, kmax, 4])
    H = np.zeros([jmax, kmax, 4])
    dQ = np.zeros([jmax, kmax, 4])
    Qn = np.zeros([jmax, kmax, 4])
    response = laser_intensity()
    x_laser = response['x_laser']
    S_laser = response['S_laser']

    for n in range(nmax):
        if (n + 1) % round(nmax / 2) == 0:
            print(round((n + 1) / nmax, 2) * 100, '%')
            print('CFL:', calc_CFL(Q))
        dQold = Q - Qn  # Qn - Qn-1
        Qn = Q.copy()  # QがQmになる
        S = np.zeros([jmax, kmax, 4])
        for k in range(kmax):
            S[(x > x_laser[n, k] - l) & (x < x_laser[n, k]), k, 3] = eta * S_laser[n, k] / l * 1e9  # レーザー熱入力w
        # 内部反復(inner iteration)
        for ttt in range(iimax):
            Ap, Am, sigma_x, Bp, Bm, sigma_y, C, H = calc_jacobian_matrix(Q)

            # 右辺の計算準備
            qLx, qRx, qLy, qRy = MUSCL(Q, order, kappa)  # 保存変数→基本変数

            Roe_flux(qLx, qRx, qLy, qRy, E, F)  # 保存変数に戻り、Eをアップデート
            RHS[1:, 1:] = dtdx * (E[1:, 1:] - E[:-1, 1:]) + dtdy * (F[1:, 1:] - F[1:, :-1]) + dt * (H[1:, 1:] / Y[1:, 1:, np.newaxis] - S[1:, 1:])  # 前進差分
            # RHS[1:-1, 1:-1] = 0.5 * dtdx * (E[2:, 1:-1] - E[:-2, 1:-1]) + 0.5 * dtdy * (F[1:-1, 2:] - F[1:-1, :-2]) + dt * (H[1:-1, 1:-1] / Y[1:-1, 1:-1, np.newaxis] - S[1:-1, 1:-1])  # 中心差分
            if np.sum(np.isnan(F)) > 0:
                print('F has nan')
                sys.exit()
            dQ = np.zeros([jmax, kmax, 4])
            D = np.ones([jmax, kmax]) + (dtdx * sigma_x + dtdy * sigma_y) * theta_1 / (1 + theta_2)

            # スイープスタート(順次dQをアップデートするためスライスが使えない)
            # 第一スイープ
            for j in range(1, jmax - 1):
                for k in range(1, kmax - 1):
                    dQ[j, k] = (-(Q[j, k] - Qn[j, k])
                                - RHS[j, k] / (1 + theta_2)
                                + dtdx * theta_1 / (1 + theta_2) * Ap[j - 1, k] @ dQ[j - 1, k]
                                + dtdy * theta_1 / (1 + theta_2) * Bp[j, k - 1] @ dQ[j, k - 1]
                                + theta_2 / (1 + theta_2) * dQold[j, k]) / D[j, k]

            if np.sum(np.isnan(dQ))>0:
                print('dQ_rho_1 has nan')
                np.savetxt('error/dQ_rho.csv', dQ[:, :, 0], fmt='%.2e')
                np.savetxt('error/rho.csv', Q[:, :, 0], fmt='%.2e')
                np.savetxt('error/u.csv', Q[:, :, 1], fmt='%.2e')
                np.savetxt('error/v.csv', Q[:, :, 2], fmt='%.2e')
                np.savetxt('error/e.csv', Q[:, :, 3], fmt='%.2e')

                fig = plt.figure(figsize=(10,8), dpi=100)  # グラフのサイズ
                plt.rcParams["font.size"] = 10  # グラフの文字サイズ

                ax1 = fig.add_subplot(1, 1, 1, projection='3d')
                ax1.plot_surface(X * 1e3, Y * 1e3, dQ[:, :, 0], color='red', linewidth=1.5, label='Numerical')
                ax1.set_xlabel(r'$z$ mm')
                ax1.set_ylabel(r'$r$ mm')
                ax1.set_zlabel(r'$dQ_rho$')

                plt.show()
                sys.exit()

            # 第二, 三スイープ
            for j in range(jmax - 2, 0, -1):
                for k in range(kmax - 2, 0, -1):
                    dQ[j, k] = np.linalg.inv(np.eye(4) + dt * theta_1 / (1 + theta_2) * C[j, k] / Y[j, k]) @ (dQ[j, k] - (dtdx * theta_1 / (1 + theta_2) * Am[j + 1, k] @ dQ[j + 1, k] + dtdy * theta_1 / (1 + theta_2) * Bm[j, k + 1] @ dQ[j, k + 1])) / D[j, k]
                    # dQ[j, k] = dQ[j, k] - (dtdx * theta_1 / (1 + theta_2) * Am[j + 1, k] @ dQ[j + 1, k] + dtdy * theta_1 / (1 + theta_2) * Bm[j, k + 1] @ dQ[j, k + 1]) / D[j, k]
            if np.sum(np.isnan(dQ))>0:
                print('dQ_rho_2 has nan')
                np.savetxt('error/dQ_rho.csv', dQ[:, :, 0], fmt='%.2e')
                np.savetxt('error/rho.csv', Q[:, :, 0], fmt='%.2e')
                np.savetxt('error/u.csv', Q[:, :, 1], fmt='%.2e')
                np.savetxt('error/v.csv', Q[:, :, 2], fmt='%.2e')
                np.savetxt('error/e.csv', Q[:, :, 3], fmt='%.2e')

                fig = plt.figure(figsize=(10,8), dpi=100)  # グラフのサイズ
                plt.rcParams["font.size"] = 10  # グラフの文字サイズ

                ax1 = fig.add_subplot(1, 1, 1, projection='3d')
                ax1.plot_surface(X * 1e3, Y * 1e3, dQ[:, :, 0], color='red', linewidth=1.5, label='Numerical')
                ax1.set_xlabel(r'$z$ mm')
                ax1.set_ylabel(r'$r$ mm')
                ax1.set_zlabel(r'$dQ_rho$')

                plt.show()
                sys.exit()
            # 収束判定
            norm = np.zeros(4)
            for i in range(4):
                norm[i] = np.linalg.norm(dQ[1:-1, 1:-1, i], 1)  # / np.linalg.norm(RHS[1:-1, 1:-1, i], 1)
            if norm[0] < norm_limit and norm[1] < norm_limit and norm[2] < norm_limit and norm[3] < norm_limit:
                break

            Q += dQ
        boundary_condition(Q)
    return response


def output_graphs(rho, u, v, p, n, response=None, Q0=None, show=False):
    """
    結果の可視化
    """
    aspect = kmax * dy / jmax / dx
    fig = plt.figure(figsize=(8,8), dpi=100)  # グラフのサイズ
    plt.rcParams["font.size"] = 10  # グラフの文字サイズ

    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.plot_surface(X * 1e3, Y * 1e3, rho, color='red', linewidth=1.5, label='Numerical')
    ax1.set_xlabel(r'$z$ mm')
    ax1.set_ylabel(r'$r$ mm')
    ax1.set_zlabel(r'$\rho$')
    ax1.set_zlim([0, 5])
    ax1.get_proj = lambda: np.dot(Axes3D.get_proj(ax1), np.diag([1, aspect, 1, 1]))

    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    ax2.plot_surface(X * 1e3, Y * 1e3, u / 1e3, color='red', linewidth=1.5, label='Numerical')
    # plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax2.set_xlabel(r'$z$ mm')
    ax2.set_ylabel(r'$r$ mm')
    ax2.set_zlabel('$u$ km/s')
    ax2.set_zlim([-1, 10])
    ax2.get_proj = lambda: np.dot(Axes3D.get_proj(ax2), np.diag([1, aspect, 1, 1]))
    # plt.legend()

    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3.plot_surface(X * 1e3, Y * 1e3, v / 1e3, color='red', linewidth=1.5, label='Numerical')
    # plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax3.set_xlabel(r'$z$ mm')
    ax3.set_ylabel(r'$r$ mm')
    ax3.set_zlabel('$v$ km/s')
    ax3.set_zlim([0, 5])
    ax3.get_proj = lambda: np.dot(Axes3D.get_proj(ax3), np.diag([1, aspect, 1, 1]))

    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    ax4.plot_surface(X * 1e3, Y * 1e3, p / 1e5, color='red', linewidth=1.5, label='Numerical')
    # plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax4.set_xlabel(r'$z$ mm')
    ax4.set_ylabel(r'$r$ mm')
    ax4.set_zlabel('$p$ atm')
    ax4.set_zlim([0, 250])
    ax4.get_proj = lambda: np.dot(Axes3D.get_proj(ax4), np.diag([1, aspect, 1, 1]))
    # plt.legend()

    plt.tight_layout()
    fig.savefig(f'2Daxi/fig/img{n}.png')
    if show:
        plt.show()
    plt.clf()
    plt.close()

    if response is not None:
        x_laser_last = response['x_laser'][-1, :]
        u_ionz = response['u_ionz'][1:, 0]
        S_laser = response['S_laser'][1:, 0]
        Plateau = response['Plateau']
        t = np.arange(0, nmax) * dt
        T, Y_ = np.meshgrid(t, y, indexing='ij')

        fig = plt.figure(figsize=(12,4), dpi=100)  # グラフのサイズ
        plt.rcParams["font.size"] = 10  # グラフの文字サイズ
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.plot(S_laser, u_ionz / 1e3, color='green', linewidth=1.5)
        ax1.grid(color='black', linestyle='dotted', linewidth=0.5)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r'Peak Laser Intensity $S_{peak}$ [GW/m^2]')
        ax1.set_ylabel(r'Propagation velocity of ionization front $V$ [km/s]')
        ax1.set_xticks(np.array([100, 300, 500, 700, 1000, 2000]))
        ax1.set_xticklabels(np.array([100, 300, 500, 700, 1000, 2000]))
        ax1.set_yticks(np.array([1, 3, 5, 7, 10]))
        ax1.set_yticklabels(np.array([1, 3, 5, 7, 10]))

        ax2 = fig.add_subplot(1, 3, 3)
        ax2.plot(y * 1e3, x_laser_last * 1e3, color='green', linewidth=1.5, label='Numerical')
        ax2.grid(color='black', linestyle='dotted', linewidth=0.5)
        ax2.set_xlabel(r'$r$ mm')
        ax2.set_ylabel(r'Propagation velocity of ionization front $z$ mm')
        ax2.set_xlim([0, 2.5])
        ax2.set_ylim([0, 5])

        fig.savefig('2Daxi/summary/laser.png')
        fig = plt.figure(figsize=(8, 4), dpi=100)  # グラフのサイズ
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax1.plot_surface(T * 1e6, Y_ * 1e3, Plateau / 1e5, color='green', linewidth=1.5, label='Numerical')
        ax1.grid(color='black', linestyle='dotted', linewidth=0.5)
        ax1.set_xlabel(r'$t$ μs')
        ax1.set_ylabel(r'$r$ mm')
        ax1.set_zlabel(r'plateau pressure $P$ atm')

        PY = Plateau * Y_
        plateau = dy * (2 * np.sum(PY, axis=1) - PY[:, 0] - PY[:, -1]) / y[-1]**2
        np.savetxt('2Daxi/csv/plateau_ave.csv', plateau, fmt='%.2e')
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(t * 1e6, plateau / 1e5, color='green', linewidth=1.5, label='Numerical')
        ax2.grid(color='black', linestyle='dotted', linewidth=0.5)
        ax2.set_xlabel(r'$t$ μm')
        ax2.set_ylabel(r'Plateau pressure $p$ atm')
        # ax2.set_xlim([0, 2.5])
        # ax2.set_ylim([0, 5])

        fig.savefig('2Daxi/summary/plateau.png')

        if show:
            plt.show()
        plt.clf()
        plt.close()


def make_animation(print_interval, n=nmax):
    images = []
    for i in range(0, n, print_interval):
        im = Image.open('2Daxi/fig/img' + str(i) + '.png')
        images.append(im)

    images[0].save('2Daxi/summary/Qt.gif', save_all=True, append_images=images[1:n], duration=40, loop=0)


def make_plateau(print_interval, n=nmax):
    divn = int(n / print_interval)
    Plateau = np.zeros([divn, kmax])
    plateau_fit1 = np.zeros(divn)
    plateau_fit2 = np.zeros(divn)
    a = 22.28578216
    b = 0.440551007
    c = 0
    d = 4.608337327
    e = 0.892930883
    f = 5.033035056
    for i in range(0, n, print_interval):
        p = np.loadtxt(f'2Daxi/p/{i}.csv')
        Plateau[int(i / print_interval), :] = p[0, :]
        plateau_fit1[int(i / print_interval)] = a * np.exp(-b * i * dt * 1e6) + c + d * np.exp(-e * i * dt * 1e6) * np.cos(f * i * dt * 1e6)  # 単位atm
        plateau_fit2[int(i / print_interval)] = a * np.exp(-b * i * dt * 1e6) + c  # 単位atm
    t = np.arange(0, divn) * dt * print_interval
    T, Y_ = np.meshgrid(t, y, indexing='ij')

    fig = plt.figure(figsize=(8, 8), dpi=100)  # グラフのサイズ
    ax1 = fig.add_subplot(1, 1, 1, projection='3d')
    ax1.plot_surface(T * 1e6, Y_ * 1e3, Plateau / 1e5, color='red', linewidth=1.5, label='Numerical')
    ax1.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax1.set_xlabel(r'$t$ μs')
    ax1.set_ylabel(r'$r$ mm')
    ax1.set_zlabel(r'plateau pressure $P$ atm')
    fig.savefig('2Daxi/summary/plateau_2d.png')
    plt.show()

    fig = plt.figure(figsize=(8, 8), dpi=100)  # グラフのサイズ
    PY = Plateau * Y_
    plateau = dy * (2 * np.sum(PY, axis=1) - PY[:, 0] - PY[:, -1]) / y[-1]**2
    np.savetxt('2Daxi/csv/plateau_ave.csv', plateau, fmt='%.2e')
    ax2 = fig.add_subplot(1, 1, 1)
    ax2.plot(t * 1e6, plateau / 1e5, color='red', linewidth=1.5, label='Numerical')
    ax2.plot(t * 1e6, plateau_fit1, color='green', linewidth=1.5, label='experiment_fit1')
    ax2.plot(t * 1e6, plateau_fit2, color='grey', linewidth=1.5, label='experiment_fit2')
    ax2.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax2.set_xlabel(r'$t$ μm')
    ax2.set_ylabel(r'Plateau pressure $p$ atm')
    plt.legend()

    fig.savefig('2Daxi/summary/plateau_ave.png')
    plt.show()

    # ax2.set_xlim([0, 2.5])
    # ax2.set_ylim([0, 5])


if __name__ == '__main__':
    update_data = True
    show = False
    # shocktube = True
    print_interval = 40

    order = 2
    kappa = 1
    Q0 = init()
    # Qext = strict_answer()
    os.makedirs('2Daxi/rho', exist_ok=True)
    os.makedirs('2Daxi/u', exist_ok=True)
    os.makedirs('2Daxi/v', exist_ok=True)
    os.makedirs('2Daxi/p', exist_ok=True)
    os.makedirs('2Daxi/e', exist_ok=True)
    os.makedirs('2Daxi/csv', exist_ok=True)
    os.makedirs('2Daxi/fig', exist_ok=True)
    # if update_data is True:
    #     start = time()
    #     print(f'dt = {dt*1e9} ns, tmax = {dt*nmax*1e6} μs\ndx = {dx*1e3} mm, x = {dx*jmax*1e3} mm\ndy = {dy*1e3} mm, y = {dy*kmax*1e3} mm')
    #     Q = Q0.copy()
    #     # response = implicit_solution(Q, order, kappa, nmax, iimax=50)
    #     # response = explicit_runge(Q, order, kappa, nmax)
    #     response = explicit_fractional_timestep(Q, order, kappa, nmax, print_interval=print_interval)
    #     print(time() - start)
    #     # Roe_FDS(Q, order, kappa, nmax)
    #     rho = Q[:, :, 0]
    #     u = Q[:, :, 1] / Q[:, :, 0]
    #     v = Q[:, :, 2] / Q[:, :, 0]
    #     p = (gamma - 1.0) * (Q[:, :, 3] - 0.5 * (Q[:, :, 1]**2 + Q[:, :, 2]**2) / Q[:, :, 0])  # / 1e5
    #     np.savetxt('2Daxi/csv/rho.csv', Q[:, :, 0], fmt='%.2e')
    #     np.savetxt('2Daxi/csv/u.csv', u, fmt='%.2e')
    #     np.savetxt('2Daxi/csv/v.csv', v, fmt='%.2e')
    #     np.savetxt('2Daxi/csv/p.csv', p, fmt='%.2e')
    #     np.savetxt('2Daxi/csv/x_laser.csv', response['x_laser'], fmt='%.2e')
    #     np.savetxt('2Daxi/csv/S_laser.csv', response['S_laser'], fmt='%.2e')
    #     np.savetxt('2Daxi/csv/u_ionz.csv', response['u_ionz'], fmt='%.2e')
    #     np.savetxt('2Daxi/csv/plateau.csv', response['Plateau'], fmt='%.2e')
    # else:
    #     rho = np.loadtxt('2Daxi/csv/rho.csv')
    #     u = np.loadtxt('2Daxi/csv/u.csv')
    #     v = np.loadtxt('2Daxi/csv/v.csv')
    #     p = np.loadtxt('2Daxi/csv/p.csv')
    #     x_laser = np.loadtxt('2Daxi/csv/x_laser.csv')
    #     S_laser = np.loadtxt('2Daxi/csv/S_laser.csv')
    #     u_ionz = np.loadtxt('2Daxi/csv/u_ionz.csv')
    #     response = {}
    #     response['x_laser'] = x_laser
    #     response['S_laser'] = S_laser
    #     response['u_ionz'] = u_ionz

    # output_graphs(rho, u, v, p, response=response, n=nmax, show=show)
    # make_animation(print_interval, n=10840)
    make_plateau(print_interval, n=10560)
    # laser_propagation()
