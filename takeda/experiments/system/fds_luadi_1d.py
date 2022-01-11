import numpy as np
import matplotlib.pyplot as plt
from time import time


dt = 1e-7

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

PI = 1  # .013e5
RHOI = 1  # PI / R / T_0
UI = 0.0

PE = 0.1  # .013e4
RHOE = 0.1  # PE / R / T_0
UE = 0.0


jmax = 2000
dx = 0.01e-3
xmin, xmax = 0.0, dx * (jmax-1)
xmid = np.average([xmin, xmax])
x = np.linspace(xmin, xmax, jmax)
# dx = (xmax - xmin) / (jmax - 1)

dtdx = dt / dx

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
eta = 0.1


# ### Roeスキームによる計算


def init(shocktube=False):
    """
    Q の初期条件を決定
    :return Q:
    """
    Q = np.zeros([jmax, 3])
    if shocktube:
        Q[x <= xmid, 0] = RHOI
        Q[x <= xmid, 1] = RHOI * UI
        Q[x <= xmid, 2] = (PI / (gamma - 1.0) + 0.5 * RHOI * UI ** 2)

        Q[x > xmid, 0] = RHOE
        Q[x > xmid, 1] = RHOE * UE
        Q[x > xmid, 2] = (PE / (gamma - 1.0) + 0.5 * RHOE * UE ** 2)
    else:
        Q[:, 0] = rho_0
        Q[:, 1] = 0
        Q[:, 2] = p_0 / (gamma - 1.0)
    return Q


def laser_intensity(nmax):
    x_laser = np.zeros(nmax)
    u_laser = np.zeros(nmax)
    t = np.arange(0, dt * nmax, dt)
    Power_laser = 8.15 * np.exp(-0.866 * t * 1e6)

    W_G = W_G0 * 1e-3
    W_T = W_T0 * 1e-3
    S_laser0 = R_peak * Power_laser / 4 / W_G / W_T * 1e-3  # 波頭のレーザー強度

    for j in range(1, nmax):
        u_ionz_line1 = slope * (S_laser0[j] * 1e9 / rho_0 / ss**3) ** intercept * ss  # m/s 波頭の速度
        u_ionz_line3 = slope_low * (S_laser0[j] * 1e9 / rho_0 / ss**3) ** intercept_low * ss
        if (u_ionz_line1 > u_ionz_line3):
            u_ionz = u_ionz_line1  # Line1, 松井さん博論
            b = intercept
        else:
            u_ionz = u_ionz_line3  # Line3, 松井さん博論
            b = intercept_low
        u_laser[j] = u_ionz
        x_laser[j] = x_laser[j - 1] + u_ionz * dt  # m 電離波面の波頭進展位置=Laserによって加熱される最初の位置
    S_laser = R_peak * Power_laser / 4 / W_G / W_T * 1e-3  # 局所レーザー強度, GW/m2
    # print(S_laser, x_laser, sep='\n')
    return S_laser, x_laser, u_laser


def calc_CFL(Q):
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    u = rhou / rho
    p = (gamma - 1.0) * (e - 0.5 * rho * u ** 2)

    c = np.sqrt(gamma * p / rho)
    sp = c + np.abs(u)
    return max(sp) * dtdx


def Qtoq(Q):
    """
    保存量ベクトルQから基本量ベクトルqを形成する
    :param Q: 保存量ベクトルrho, rhou, e
    :return q: 基本量ベクトルrho, u, p
    """
    q = np.zeros([jmax, 3])
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]
    u = rhou / rho

    q[:, 0] = rho
    q[:, 1] = u
    q[:, 2] = (gamma - 1.0) * (e - 0.5 * rho * u**2)  # p

    return q


def qtoQ(q):
    """
    基本量ベクトルqから保存量ベクトルQを形成する
    :param q: 基本量ベクトル
    :return Q: 保存量ベクトル
    """
    Qc = q.copy()
    rho, u, p = q[:, 0], q[:, 1], q[:, 2]
    e = p / (gamma - 1.0) + 0.5 * rho * u**2

    Qc[:, 0] = rho
    Qc[:, 1] = rho * u
    Qc[:, 2] = e

    return Qc


def Roe_flux(qL, qR, E):
    """
    van LeerのMUSTL型FDS法(高次精度)
    Roe平均を用いて数値流束E_j+1/2を求める
    :param qL, qR: 内挿する基本変数、qLはj_+1/2, qRはj-1/2
    :param E: 更新するE_j+1/2
    """
    dQ, EL, ER, AQ = np.zeros([jmax-1, 3]), np.zeros([jmax-1, 3]), np.zeros([jmax-1, 3]), np.zeros([jmax-1, 3])
    Lambda, R, Rinv = np.zeros([jmax-1, 3, 3]), np.ones([jmax-1, 3, 3]), np.zeros([jmax-1, 3, 3])

    rhoL, uL, pL = qL[:-1, 0], qL[:-1, 1], qL[:-1, 2]
    rhoR, uR, pR = qR[1:, 0], qR[1:, 1], qR[1:, 2]

    eL = pL / (gamma - 1.0) + 0.5 * rhoL * uL**2
    eR = pR / (gamma - 1.0) + 0.5 * rhoR * uR**2

    HL = (eL + pL) / rhoL
    HR = (eR + pR) / rhoR

    # Roe平均 式(6.38)
    sqrhoL = np.sqrt(rhoL)
    sqrhoR = np.sqrt(rhoR)

    uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
    HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
    cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * uAVE**2))

    b1 = 0.5 * (gamma - 1.0) * uAVE**2 / cAVE**2
    b2 = (gamma - 1.0) / cAVE**2

    Lambda[:, 0, 0] = np.abs(uAVE - cAVE)
    Lambda[:, 1, 1] = np.abs(uAVE)
    Lambda[:, 2, 2] = np.abs(uAVE + cAVE)

    dQ[:, 0] = rhoR - rhoL
    dQ[:, 1] = rhoR * uR - rhoL * uL
    dQ[:, 2] = eR - eL

    EL[:, 0] = rhoL * uL
    EL[:, 1] = pL + rhoL * uL**2
    EL[:, 2] = (eL + pL) * uL
    ER[:, 0] = rhoR * uR
    ER[:, 1] = pR + rhoR * uR**2
    ER[:, 2] = (eR + pR) * uR

    R[:, 1, 0] = uAVE - cAVE
    R[:, 1, 1] = uAVE
    R[:, 1, 2] = uAVE + cAVE
    R[:, 2, 0] = HAVE - uAVE * cAVE
    R[:, 2, 1] = 0.5 * uAVE**2
    R[:, 2, 2] = HAVE + uAVE * cAVE

    Rinv[:, 0, 0] = 0.5 * (b1 + uAVE / cAVE)
    Rinv[:, 0, 1] = -0.5 * (b2 * uAVE + 1 / cAVE)
    Rinv[:, 0, 2] = 0.5 * b2
    Rinv[:, 1, 0] = 1.0 - b1
    Rinv[:, 1, 1] = b2 * uAVE
    Rinv[:, 1, 2] = -b2
    Rinv[:, 2, 0] = 0.5 * (b1 - uAVE / cAVE)
    Rinv[:, 2, 1] = -0.5 * (b2 * uAVE - 1 / cAVE)
    Rinv[:, 2, 2] = 0.5 * b2

    for j in range(jmax - 1):
        AQ[j] = R[j] @ Lambda[j] @ Rinv[j] @ dQ[j]

    E[:-1] = 0.5 * (ER + EL - AQ)  # 式(6.43)


def minmod(x, y):
    """
    制限関数
    :param x,y: 制限される前のdelta二つ。どちらもsize(3,1)のndarray
    :return delta: 制限されたdelta
    """
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(np.abs(x), sgn * y), 0.0)


def MUSCL(Q, order: int = 2, kappa: float = 0, limit=minmod):
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

    qL = q.copy()  # 1st order
    qR = q.copy()  # 1st order

    if order == 2 or order == 3:
        # 2nd / 3rd order & minmod limitter
        dq = np.zeros([jmax, 3])
        dq[:-1] = q[1:] - q[:-1]

        b = (3.0 - kappa) / (1.0 - kappa)  # 式(2.74)

        Dp = np.zeros([jmax, 3])
        Dm = np.zeros([jmax, 3])

        # 制限関数の導入によりTVD化
        for j in range(1, jmax - 1):
            Dp[j] = limit(dq[j], b * dq[j - 1])  # 式(2.73a)
            Dm[j] = limit(dq[j - 1], b * dq[j])     # 式(2.73b)

            qL[j] += 0.25 * ((1.0 + kappa) * Dp[j] + (1.0 - kappa) * Dm[j])  # QL_j+1/2 式(6.28a)
            qR[j] -= 0.25 * ((1.0 - kappa) * Dp[j] + (1.0 + kappa) * Dm[j])  # QR_j-1/2 式(6.28b)
        # 境界
        # qL[0] = qL[1]
        # qR[-1] = qR[-2]

    return qL, qR


# 陽解法
def Roe_FDS(Q, order, kappa, nmax, print_interval=2):
    """
    二段階ルンゲクッタ
    """
    E = np.zeros([jmax, 3])

    for n in range(nmax):
        if n % print_interval == 0:
            print(f'n = {n : 4d} : CFL = {calc_CFL(Q) : .4f}')

        Qold = Q.copy()

        coefs = [0.5, 1.0]
        for coef in coefs:
            QL, QR = MUSCL(Qold, order, kappa)

            Roe_flux(QL, QR, E)
            for j in range(1, jmax - 1):
                Qold[j] = Q[j] - coef * dtdx * (E[j] - E[j-1])

            Qold[0] = Q[0]
            Qold[-1] = Q[-1]

        Q[:] = Qold[:]


def calc_A(Q):
    """
    A = dE/dQ
    近似LDU分解
    """
    # 基本変数の用意
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]
    u = rhou / rho
    p = (gamma - 1) * (e - 0.5 * rho * u**2)
    H = (e + p) / rho
    c = np.sqrt(gamma * p / rho)

    A = np.zeros([jmax, 3, 3])

    A[:, 0, 1] = 1
    A[:, 1, 0] = -(3 - gamma) / 2 * u**2
    A[:, 1, 1] = (3 - gamma) * u
    A[:, 1, 2] = gamma - 1
    A[:, 2, 0] = ((gamma - 1) / 2 * u**2 - H) * u
    A[:, 2, 1] = H - (gamma - 1) * u**2
    A[:, 2, 2] = gamma * u

    sigma_x = abs(u) + c  # 一番大きい値で近似する

    # A+
    Ap = A.copy()
    Ap[:, 0, 0] += sigma_x
    Ap[:, 1, 1] += sigma_x
    Ap[:, 2, 2] += sigma_x
    Ap = Ap / 2

    # A-
    Am = A.copy()
    Am[:, 0, 0] -= sigma_x
    Am[:, 1, 1] -= sigma_x
    Am[:, 2, 2] -= sigma_x
    Am = Am / 2

    return Ap, Am, sigma_x


# 陰解法
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
    RHS = np.zeros([jmax, 3])
    E = np.zeros([jmax, 3])
    response = {}

    dQ = np.zeros([jmax, 3])
    S_laser, x_laser, u_laser = laser_intensity(nmax)
    response['S_laser'] = S_laser
    response['x_laser'] = x_laser
    response['u_laser'] = u_laser

    for n in range(nmax):
        if (n + 1) % round(nmax / 4) == 0:
            print(round(n / nmax, 2) * 100, '%')
        Qn = Q.copy()  # QがQmになる
        S = np.zeros([jmax, 3])
        S[(x > x_laser[n] - l) & (x <= x_laser[n]), 2] = eta * S_laser[n] / l * 1e9
        # print(S)

        # 内部反復(inner iteration)
        for ttt in range(iimax):
            Ap, Am, sigma_x = calc_A(Q)

            # 右辺の計算準備
            qL, qR = MUSCL(Q, order, kappa)  # 保存変数→基本変数
            Roe_flux(qL, qR, E)  # 保存変数に戻り、Eをアップデート
            RHS[1:] = dtdx * (E[1:] - E[:-1]) - dt * S[1:]

            dQ = np.zeros([jmax, 3])

            # 第一スイープ
            for j in range(1, jmax - 1):
                dQ[j] = (-(Q[j] - Qn[j]) - RHS[j] + dtdx * Ap[j - 1] @ dQ[j - 1]) / (1 + dtdx * sigma_x[j])

            # 第二, 三スイープ
            for j in range(jmax - 2, 0, -1):
                dQ[j] = dQ[j] - (dtdx * Am[j + 1] @ dQ[j + 1]) / (1 + dtdx * sigma_x[j])

            # 収束判定
            norm = np.zeros(3)
            for i in range(3):
                norm[i] = np.linalg.norm(dQ[1:-1, i], 1)   # / np.linalg.norm(RHS[1:-1, i], 1)
            if norm[0] < norm_limit and norm[1] < norm_limit and norm[2] < norm_limit:
                break
            Q += dQ
            Q[0] = Q[1]
            Q[-1] = Q[-2]

    return response


# 陽解法
def MacCormack(Q, eps_c, nmax, interval=5):
    E = np.zeros([jmax, 3])

    for n in range(nmax):
        if n % interval == 0:
            print(f'n = {n : 4d} : CFL = {calc_CFL(Q) : .4f}')

        Qs = Q.copy()

        E_flux(Q, E)
        for j in range(1, jmax - 1):
            Qs[j] = Q[j] - dtdx * (E[j] - E[j-1]) # 式(6.10)

        E_flux(Qs, E)
        for j in range(1, jmax - 2):
            Q[j] = 0.5 * (Q[j] + Qs[j]) - 0.5 * dtdx * (E[j + 1] - E[j]) # 式(6.10)

        Qb = Q.copy()
        for j in range(1, jmax - 1):
            D1 = Qb[j - 1] - 2.0 * Qb[j] + Qb[j + 1]
            D2 = Qb[j - 1] + 2.0 * Qb[j] + Qb[j + 1]
            k = eps_c * np.linalg.norm(D1) / np.linalg.norm(D2) # 式(6.12)
            Q[j] += k * D1 # 式(6.11)


# 厳密解の計算
def strict_answer():
    """
    :return Qext:
    """
    Pext = np.zeros([jmax, 3])
    Qext = np.zeros([jmax, 3])

    GUESS = 1.0
    FINC = 0.01
    itemax1 = 5000
    itemax2 = 500

    CI = np.sqrt(gamma * PI / RHOI)
    CE = np.sqrt(gamma * PE / RHOE)
    P1P5 = PI / PE

    GAMI = 1.0 / gamma
    GAMF = (gamma - 1.0) / (2.0 * gamma)
    GAMF2 = (gamma + 1.0) / (gamma - 1.0)
    GAMFI = 1.0 / GAMF

    for it1 in range(itemax1):
        for it2 in range(itemax2):
            SQRT1 = (gamma - 1.0) * (CE / CI) * (GUESS - 1.0)
            SQRT2 = np.sqrt(2.0 * gamma * (2.0 * gamma + (gamma + 1.0) * (GUESS - 1.0)))
            FUN = GUESS * (1.0 - (SQRT1 / SQRT2)) ** (-GAMFI)
            DIF = P1P5 - FUN

            if np.abs(DIF) <= 0.000002:
                break

            if DIF >= 0.0:
                GUESS += FINC
            else:
                GUESS -= FINC
                FINC = 0.5 * FINC
        else:
            continue

        break

    P4P5 = GUESS
    P4 = PE * P4P5
    P3P1 = P4P5 / P1P5
    P3 = P3P1 * PI

    R4R5 = (1.0 + GAMF2 * P4P5) / (GAMF2 + P4P5)
    RHO4 = RHOE * R4R5
    U4 = CE * (P4P5 - 1.0) * np.sqrt(2.0 * GAMI / ((gamma + 1.0) * P4P5 + (gamma - 1.0)))
    C4 = np.sqrt(gamma * P4 / RHO4)

    R3R1 = P3P1 ** GAMI
    RHO3 = RHOI * R3R1
    U3 = 2.0 * CI / (gamma - 1.0) * (1.0 - P3P1 ** GAMF)
    C3 = np.sqrt(gamma * P3 / RHO3)
    CS =  CE * np.sqrt(0.5 * ((gamma - 1.0) * GAMI + (gamma + 1.0) * GAMI * P4 / PE))

    TOT = 0.0
    EPST = 1.0e-14
    for n in range(nmax):
        TOT = TOT + dt
        rad = dt / dx

        x1 = xmid - CI * TOT
        x2 = xmid - (CI - 0.5 * (gamma + 1.0) * U3) * TOT
        x3 = xmid + U3 * TOT
        x4 = xmid + CS * TOT

        for j in range(jmax):
            xx = x[j]
            if xx <= x1:
                Qext[j, 0] = RHOI
                Qext[j, 1] = RHOI * UI
                Qext[j, 2] = PI / (gamma - 1.0) + 0.5 * UI * Qext[j, 1]
                Pext[j] = PI
            elif xx <= x2:
                UT = UI + (U3 - UI) / ((x2 - x1) + EPST) * ((xx - x1) + EPST)
                RTRI = (1.0 - 0.5 * (gamma - 1.0) * UT / CI) ** (2.0 / (gamma - 1.0))
                RT = RHOI * RTRI
                PT = RTRI ** gamma * PI
                Qext[j, 0] = RT
                Qext[j, 1] = RT * UT
                Qext[j, 2] = PT / (gamma - 1.0) + 0.5 * UT * Qext[j, 1]
                Pext[j] = PT
            elif xx <= x3:
                Qext[j, 0] = RHO3
                Qext[j, 1] = RHO3 * U3
                Qext[j, 2] = P3 / (gamma - 1.0) + 0.5 * U3 * Qext[j, 1]
                Pext[j] = P3
            elif xx <= x4:
                Qext[j, 0] = RHO4
                Qext[j, 1] = RHO4 * U4
                Qext[j, 2] = P4 / (gamma - 1.0) + 0.5 * U4 * Qext[j, 1]
                Pext[j] = P4
            else:
                Qext[j, 0] = RHOE
                Qext[j, 1] = RHOE * UE
                Qext[j, 2] = PE / (gamma - 1.0) + 0.5 * UE * Qext[j, 1]
                Pext[j] = PE
    return Qext


def output_graphs(Q0, rho, u, p, response):
    x_laser_last = response['x_laser'][-1]
    u_ionz = response['u_laser']
    S_laser = response['S_laser']

    # ### 結果の可視化
    fig = plt.figure(figsize=(10,7), dpi=100)  # グラフのサイズ
    plt.rcParams["font.size"] = 10  # グラフの文字サイズ
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(x * 1e3, Q0[:,0], color='green', linewidth=1.5, label='Numerical')
    ax1.plot(x * 1e3, rho, color='red', linewidth=1.5, label='Numerical')
    # ax1.plot(x, Qext[:,0], color='black', linewidth = 1.0, linestyle = 'dashed', label = 'Analytical')
    ax1.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax1.set_xlabel('$x$ mm')
    ax1.set_ylabel(r'$\rho$')
    # ax1.set_xlim([0,1e-5])
    ax1.axvline(x_laser_last * 1e3, ls="--", color="navy")
    # plt.legend()

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(x * 1e3, Q0[:,1]/Q0[:,0], color='green', linewidth=1.5, label='Numerical')
    ax2.plot(x * 1e3, u, color='red', linewidth=1.5, label='Numerical')
    # ax2.plot(x, Qext[:,1]/Qext[:,0], color='black', linewidth = 1.0, linestyle = 'dashed', label = 'Analitical')
    ax2.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax2.set_xlabel('$x$ mm')
    ax2.set_ylabel('$u$ km/s')
    # ax2.set_xlim([0,1e-3])
    ax2.axvline(x_laser_last * 1e3, ls="--", color="navy")
    # plt.legend()

    ax3 = fig.add_subplot(2, 2, 3)
    p0 = (gamma - 1.0) * (Q0[:,2] - 0.5 * Q0[:,1] ** 2 / Q0[:,0])  # / 1e5
    yext = (gamma - 1.0) * (Qext[:,2] - 0.5 * Qext[:,1] ** 2 / Qext[:,0])
    ax3.plot(x * 1e3, p0 / 1e5, color='green', linewidth=1.5, label='Numerical')
    ax3.plot(x * 1e3, p / 1e5, color='red', linewidth=1.5, label='Numerical')
    # ax3.plot(x, yext, color='black', linewidth = 1.0, linestyle = 'dashed',label = 'Analytical')
    ax3.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax3.set_xlabel('$x$ mm')
    ax3.set_ylabel('$p$ atm')
    # ax3.set_xlim([0, 1e-3])
    # ax3.set_ylim([0, 5])
    ax3.axvline(x_laser_last * 1e3, ls="--", color="navy")

    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(S_laser[1:], u_ionz[1:] / 1e3, color='green', linewidth=1.5, label='Numerical')
    ax4.grid(color='black', linestyle='dotted', linewidth=0.5)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('$S$ GW/m2')
    ax4.set_ylabel('$u_{ionz}$ km/s')
    # ax3.set_xlim([0, 1e-3])
    # ax3.set_ylim([0, 5])

    plt.legend()
    plt.show()

if __name__ == '__main__':
    start = time()
    nmax = 20
    print_interval = 4

    order = 2
    kappa = 0
    Q0 = init()
    Q = Q0.copy()
    response = implicit_solution(Q, order, kappa, nmax, iimax=200)
    print(round(time()-start, 2), 's')
    # Roe_FDS(Q, order, kappa, nmax)
    rho = Q[:, 0]
    u = Q[:, 1] / Q[:, 0]
    p = (gamma - 1.0) * (Q[:, 2] - 0.5 * rho * u**2)  # / 1e5

    Qext = strict_answer()
    output_graphs(Q0, rho, u, p, response=response)
