import numpy as np
import matplotlib.pyplot as plt


jmax = 101
dt = 0.002


gamma = 1.4

m = 28.966
R_0 = 8314.0
R = R_0 / m
T_0 = 293

PI = 1.013e5
RHOI = PI / R / T_0
UI = 0.0

PE = 1.013e4
RHOE = PE / R / T_0
UE = 0.0


xmin, xmid, xmax = 0.0, 0.5, 1.0
x = np.linspace(xmin, xmax, jmax)

dx = (xmax - xmin) / (jmax - 1)

dtdx = dt / dx


# ### Roeスキームによる計算


def init():
    """
    Q の初期条件を決定
    :return Q:
    """
    Q = np.zeros([jmax, 3])

    Q[x <= xmid, 0] = RHOI
    Q[x <= xmid, 1] = RHOI * UI
    Q[x <= xmid, 2] = (PI / (gamma - 1.0) + 0.5 * RHOI * UI ** 2)

    Q[x > xmid, 0] = RHOE
    Q[x > xmid, 1] = RHOE * UE
    Q[x > xmid, 2] = (PE / (gamma - 1.0) + 0.5 * RHOE * UE ** 2)

    return Q


def calc_CFL(Q):
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    u = rhou / rho
    p = (gamma - 1.0) * (e - 0.5 * rho * u ** 2)

    c = np.sqrt(gamma * p / rho)
    sp = c + np.abs(u)
    return max(sp) * dtdx


def Roe_flux(QL, QR, E):
    """
    van LeerのMUSTL型FDS法(高次精度)
    E_j+1/2を求める
    :param QL, QR: 内挿する基本変数、共にj_+1/2を採用すること
    :param E: 更新するE_j+1/2
    """
    for j in range(jmax - 1):
        rhoL, uL, pL = QL[j, 0], QL[j, 1], QL[j, 2]
        rhoR, uR, pR = QR[j + 1, 0], QR[j + 1, 1], QR[j + 1, 2]

        eL = pL / (gamma - 1.0) + 0.5 * rhoL * uL**2
        eR = pR / (gamma - 1.0) + 0.5 * rhoR * uR**2

        HL = (eL + pL) / rhoL
        HR = (eR + pR) / rhoR

        # cL = np.sqrt((gamma - 1.0) * (HL - 0.5 * uL ** 2))
        # cR = np.sqrt((gamma - 1.0) * (HR - 0.5 * uR ** 2))

        # Roe平均 式(6.38)
        sqrhoL = np.sqrt(rhoL)
        sqrhoR = np.sqrt(rhoR)

        # rhoAVE = sqrhoL * sqrhoR
        uAVE = (sqrhoL * uL + sqrhoR * uR) / (sqrhoL + sqrhoR)
        HAVE = (sqrhoL * HL + sqrhoR * HR) / (sqrhoL + sqrhoR)
        cAVE = np.sqrt((gamma - 1.0) * (HAVE - 0.5 * uAVE**2))
        # eAVE = rhoAVE * (HAVE - cAVE ** 2 / gamma)

        dQ = np.array([rhoR - rhoL, rhoR * uR - rhoL * uL, eR - eL])

        Lambda = np.diag([np.abs(uAVE - cAVE), np.abs(uAVE), np.abs(uAVE + cAVE)])

        b1 = 0.5 * (gamma - 1.0) * uAVE**2 / cAVE**2
        b2 = (gamma - 1.0) / cAVE**2

        R = np.array([[               1.0,           1.0,                1.0],
                      [       uAVE - cAVE,          uAVE,        uAVE + cAVE],
                      [HAVE - uAVE * cAVE, 0.5 * uAVE**2, HAVE + uAVE * cAVE]])

        Rinv = np.array([[0.5 * (b1 + uAVE / cAVE), -0.5 * (b2 * uAVE + 1 / cAVE), 0.5 * b2],
                         [                1.0 - b1,                     b2 * uAVE,      -b2],
                         [0.5 * (b1 - uAVE / cAVE), -0.5 * (b2 * uAVE - 1 / cAVE), 0.5 * b2]])

        AQ = R @ Lambda @ Rinv @ dQ

        EL = np.array([rhoL * uL, pL + rhoL * uL**2, (eL + pL) * uL])
        ER = np.array([rhoR * uR, pR + rhoR * uR**2, (eR + pR) * uR])

        E[j] = 0.5 * (ER + EL - AQ)  # 式(6.43)


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
    :return QL, QR: 内挿したQL,QR(基本変数)
    """
    # 基本変数で内挿する
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    Q[:, 1] = rhou / rho  # u
    Q[:, 2] = (gamma - 1.0) * (e - 0.5 * rho * Q[:, 1] ** 2)  # p

    QL = Q.copy()  # 1st order
    QR = Q.copy()  # 1st order

    if order == 2 or order == 3:
        # 2nd / 3rd order & minmod limitter
        dQ = np.zeros([jmax, 3])
        for j in range(jmax - 1):
            dQ[j] = Q[j + 1] - Q[j]

        b = (3.0 - kappa) / (1.0 - kappa)  # 式(2.74)

        Dp = np.zeros([jmax, 3])
        Dm = np.zeros([jmax, 3])

        # 制限関数の導入によりTVD化
        for j in range(1, jmax - 1):
            Dp[j] = limit(dQ[j], b * dQ[j - 1])  # 式(2.73a)
            Dm[j] = limit(dQ[j - 1], b * dQ[j])     # 式(2.73b)

            QL[j] += 0.25 * ((1.0 + kappa) * Dp[j] + (1.0 - kappa) * Dm[j])  # QL_j+1/2 式(6.28a)
            QR[j] -= 0.25 * ((1.0 - kappa) * Dp[j] + (1.0 + kappa) * Dm[j])  # QR_j-1/2 式(6.28b)
        # 境界
        QL[0] = QL[1]
        QR[-1] = QR[-2]

    return QL, QR


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
    dQ = np.zeros([jmax, 3])

    for n in range(1, nmax + 1):
        Qm = Q.copy()  # QがQnになる

        # 内部反復(inner iteration)
        for i in range(iimax):
            Ap, Am, sigma_x = calc_A(Qm)

            # 右辺の計算準備
            QL, QR = MUSCL(Qm, order, kappa)  # 保存変数→基本変数
            Roe_flux(QL, QR, E)  # 保存変数に戻る

            # 第一スイープ
            for j in range(1, jmax - 1):
                RHS = E[j] - E[j - 1]
                dQ[j] = (-(Qm - Q) - dt / dx * RHS + dt / dx * Ap[j - 1] @ dQ[j - 1]) / (1 + dt / dx * sigma_x[j])

            # 第二スイープ
            for j in range(1, jmax - 1):
                dQ[j] = (1 + dt / dx * sigma_x[j]) * dQ[j]

            # 第三スイープ
            for j in range(jmax - 2, 0, -1):
                dQ[j] = dQ[j] - dt / dx * Am[j + 1] @ dQ[j + 1] / (1 + dt / dx * sigma_x[j])

            # 収束判定
            norm_1, norm_2, norm_3 = np.linalg.norm(dQ[0], 1), np.linalg.norm(dQ[1], 1), np.linalg.norm(dQ[2], 1)
            if min(np.abs([norm_1, norm_2, norm_3])) < norm_limit:
                break
            else:
                Qm += dQ
        Q[1:-1] = Qm[1:-1]


# 陽解法
def MacCormack(Q, eps_c, nmax, interval=2):
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


if __name__ == '__main__':
    nmax = 1
    print_interval = 4

    order = 2
    kappa = 0

    Q = init()
    implicit_solution(Q, order, kappa, nmax)


    # ### 結果の可視化
    plt.figure(figsize=(8,8), dpi=100)  # グラフのサイズ
    plt.rcParams["font.size"] = 22  # グラフの文字サイズ
    plt.plot(x, Q[:,0], color='red', linewidth=1.5, label='Numerical')
    plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    # plt.legend()
    plt.show()


    plt.figure(figsize=(7,7), dpi=100)  # グラフのサイズ
    plt.rcParams["font.size"] = 22  # グラフの文字サイズ
    plt.plot(x, Q[:,1]/Q[:,0]/1e3, color='red', linewidth=1.5, label='Numerical')
    plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('$u$ km/s')
    # plt.legend()
    plt.show()


    plt.figure(figsize=(7,7), dpi=100)  # グラフのサイズ
    plt.rcParams["font.size"] = 22  # グラフの文字サイズ
    y = (gamma - 1.0) * (Q[:,2] - 0.5 * Q[:,1] ** 2 / Q[:,0]/1e5)
    plt.plot(x, y, color='red', linewidth=1.5,  label='Numerical')
    plt.grid(color='black', linestyle='dotted', linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('$p$ atm')
    plt.legend()
    plt.show()
