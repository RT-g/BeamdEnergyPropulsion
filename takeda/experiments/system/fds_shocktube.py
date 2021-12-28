import numpy as np
import matplotlib.pyplot as plt

jmax = 101
dt = 0.002

gamma = 1.4

PI = 1.0
RHOI = 1.0
UI = 0.0

PE = 0.1
RHOE = 0.1
UE = 0.0

xmin, xmid, xmax = 0.0, 0.5, 1.0
x = np.linspace(xmin, xmax, jmax)

dx = (xmax - xmin) / (jmax - 1)

dtdx = dt / dx


def init():
    Q = np.zeros([jmax, 3])

    Q[x <= xmid, 0] = RHOI
    Q[x <= xmid, 1] = RHOI * UI
    Q[x <= xmid, 2] = (PI / (gamma - 1) + 0.5 * RHOI * UI**2)

    Q[x > xmid, 0] = RHOE
    Q[x > xmid, 1] = RHOE * UE
    Q[x > xmid, 2] = (PE / (gamma - 1) + 0.5 * RHOE * UE**2)
    return Q


def calc_CFL(Q):
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    u = rhou / rho
    p = (gamma - 1) * (e - 0.5 * rho * u**2)
    c = np.sqrt(gamma * p / rho)
    sp = c + np.abs(u)

    return max(sp) * dtdx


def E_flux(Q, E):
    rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    u = rhou / rho
    p = (gamma - 1) * (e - 0.5 * rho * u**2)
    E[:, 0] = rhou
    E[:, 1] = p * rhou * u
    E[:, 2] = (e+rho) * u


def MacCormack(Q, nmax, eps_c=0.2, interval=2):
    E = np.zeros([jmax, 3])
    plt.figure(figsize=(7, 7), dpi=100)
    plt.rcParams["font.size"] = 22

    plt.plot(x, Q[:, 0], marker='o', lw=2, label='n=0')
    for n in range(nmax):
        Qs = Q.copy()

        E_flux(Q, E)
        for j in range(1, jmax-1):
            Qs[j] = Q[j] - dtdx * (E[j] - E[j-1])

        E_flux(Qs, E)
        for j in range(1, jmax-2):
            Qs[j] = 0.5 * (Q[j] + Qs[j]) - 0.5 * dtdx * (E[j+1] - E[j])

        Qb = Q.copy()

        for j in range(1, jmax-1):
            D1 = Qb[j-1] - 2.0 * Qb[j] + Qb[j+1]
            D2 = Qb[j-1] + 2.0 * Qb[j] + Qb[j+1]
            k = eps_c * np.linalg.norm(D1) / np.linalg.norm(D2)
            Q[j] += k * D1

        # if n % interval == 0:
        #     print(f'n = {n : 4d} : CFL = {calc_CFL(Q) : .4f}')
    plt.plot(x, Q[:, 0], marker='o', lw=2, label=f'n={n}')

    plt.grid(color='black', linestyle='dashed', linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('rho')
    # plt.legend()
    plt.show()

def minmod(x, y):
    sgn = np.sign(x)
    return sgn * np.maximum(np.minimum(np.abs(x), sgn * y), 0.0)


def roe_flux(QL, QR, E):
    for j in range(jmax - 1):
        rhoL, uL, eL = QL[j  , 0], QL[j  , 1] / QL[j  , 0], QL[j  , 2]
        rhoR, uR, eR = QR[j+1, 0], QR[j+1, 1] / QR[j+1, 0], QR[j+1, 2]

        pL = (gamma - 1.0) * (eL - 0.5 * rhoL * uL**2)
        pR = (gamma - 1.0) * (eR - 0.5 * rhoR * uR**2)

        HL = (eL + pL) / rhoL
        HR = (eR + pR) / rhoR

        # cL = np.sqrt((gamma - 1.0) * (HL - 0.5 * uL**2))
        # cR = np.sqrt((gamma - 1.0) * (HR - 0.5 * uR**2))

        # rho_ave = np.sqrt(rhoL * rhoR)
        u_ave = (np.sqrt(rhoL) * uL + np.sqrt(rhoR) * uR) / (np.sqrt(rhoL) + np.sqrt(rhoR))
        H_ave = (np.sqrt(rhoL) * HL + np.sqrt(rhoR) * HR) / (np.sqrt(rhoL) + np.sqrt(rhoR))
        c_ave = np.sqrt((gamma - 1.0) * (H_ave - 0.5 * u_ave**2))
        # e_ave = rho_ave * (H_ave - c_ave**2 / gamma)

        dQ = np.array([rhoR - rhoL, rhoR * uR - rhoL * uL, eR - eL])

        Lambda = np.abs(np.diag([u_ave-c_ave, u_ave, u_ave+c_ave]))

        b1 = 0.5 * (gamma - 1.0) * u_ave**2 / c_ave**2
        b2 = (gamma - 1) / c_ave**2

        R = np.array([[1, 1, 1],
                      [u_ave - c_ave, u_ave, u_ave + c_ave],
                      [H_ave - u_ave * c_ave, 0.5 * u_ave**2, H_ave + u_ave * c_ave]])
        Rinv = np.array([[0.5 * (b1 + u_ave / c_ave), -0.5 * (b2 * u_ave + 1 / c_ave), 0.5 * b2],
                         [1 - b1                    , b2 * u_ave                     , -b2],
                         [0.5 * (b1 - u_ave / c_ave), -0.5 * (b2 * u_ave - 1 / c_ave), 0.5 * b2]])

        AQ = R @ Lambda @ Rinv @ dQ

        EL = np.array([rhoL * uL, pL + rhoL * uL, (eL + pL) * uL])
        ER = np.array([rhoR * uR, pR + rhoR * uR, (eR + pR) * uR])

        E[j] = 0.5 * (ER + EL - AQ)
        return E


def MUSCL(Q, order, kappa):
    # rho, rhou, e = Q[:, 0], Q[:, 1], Q[:, 2]

    if order == 2 or order == 3:
        dQ = np.zeros([jmax, 3])
        for j in range(jmax - 1):
            dQ[j] = Q[j+1] - Q[j]

        b = (3-kappa) / (1-kappa)

        Dp = np.zeros([jmax, 3])
        Dm = np.zeros([jmax, 3])
        for j in range(1, jmax - 1):
            Dp[j] = minmod(dQ[ j ], b*dQ[j-1])
            Dm[j] = minmod(dQ[j-1], b*dQ[ j ])

        Dp[0] = Dp[1]
        Dm[0] = Dm[1]

        QL = Q.copy()
        QR = Q.copy()

        for j in range(1, jmax - 1):
            QL[j] += 0.25 * ((1-kappa) * Dp[j] + (1+kappa) * Dm[j])
            QR[j] -= 0.25 * ((1+kappa) * Dp[j] + (1-kappa) * Dm[j])

    else:
        QL = Q.copy()
        QR = Q.copy()

    return QL, QR


def roe_fds(Q, order, kappa, nmax, print_interval=2):
    E = np.zeros([jmax, 3])
    plt.figure(figsize=(7, 7), dpi=100)
    plt.rcParams["font.size"] = 22

    plt.plot(x, Q[:, 0], marker='o', lw=2, label='n=0')

    for n in range(nmax):
        Qold = Q.copy()
        # 2段階ルンゲクッタ
        coefs = [0.5, 1.0]
        for coef in coefs:
            QL, QR = MUSCL(Qold, order, kappa)
            E = roe_flux(QL, QR, E)

            for j in range(1, jmax - 1):
                Qold[j] = Q[j] - coef * dtdx * (E[j] - E[j-1])

            Qold[0] = Q[0]
            Qold[-1] = Q[-1]
        Q[:] = Qold[:]

        if n % print_interval == 0:
            print(f'n = {n : 4d} : CFL = {calc_CFL(Q) : .4f}')
    plt.plot(x, Q[:, 0], marker='o', lw=2, label=f'n={n}')

    plt.grid(color='black', linestyle='dashed', linewidth=0.5)
    plt.xlabel('x')
    plt.ylabel('rho')
    # plt.legend()
    plt.show()


nmax = 1000
print_interval = 4

order = 2

# -1 = 2nd order fully upwind
# 0 = 2nd order upwind biased
# 1/3 = 3rd order upwind biased
kappa = 0


if __name__ == '__main__':
    Q = init()
    MacCormack(Q, nmax, eps_c=0.2)
    # roe_fds(Q, order, kappa, nmax, print_interval)
