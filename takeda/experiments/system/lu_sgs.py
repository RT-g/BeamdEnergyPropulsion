import numpy as np
import os
import copy
import matplotlib.pyplot as plt

# --------------------------
# -- initial value        --
# --------------------------
nstep = 300                  # 時間ステップ数
nx0 = 100                    # 空間ステップ数
dt = 0.002                   # 時間刻み幅
dx = 0.01                    # 空間刻み幅

lbound = 1                     # 仮想境界セル数
nx = nx0+2*lbound              # 総空間セル数

# -- slope limiter --

k_muscl=1/3                        # muscl精度
b_muscl=(3-k_muscl)/(1-k_muscl)

# --定数--
gamma=1.4                    # 比熱比

norm_ok=1.0e-6

# -- 出力--
dir_name="sod_lusgs_0.002s_py" # 出力フォルダ名
out_name_front="time"                 # 出力ファイル名（先頭）
out_name_back="d-3"

# --------------------------
# -- function setup       --
# --------------------------
def setup():
    """
    基本量と保存量の設定
    初期値入力
    """
    global x, bol, bor, qf, Qc
    u = [0.0]*nx               # 速度
    rho = [0.0]*nx             # 密度
    p = [0.0]*nx               # 圧力
    e = [0.0]*nx               # エネルギー
    x = [0.0]*nx               # 位置

    """
    e=p/(r-1)+rho*u^2/2
    """
    for i in range(nx):
        u[i] = 0.0

        if i <= nx*0.5:
            rho[i] = 1.0
            p[i] = 1.0

        else:
            rho[i] = 0.125
            p[i] = 0.1

        e[i] = p[i]/(gamma-1)+rho[i]*(u[i]**2)/2
        x[i] = i*dx-dx/2


    bol = [0.0]*3             # 左端仮想セル
    bor = [0.0]*3             # 右端仮想セル
    for j in range(3):
        if j == 0:
            bol[j] = rho[0]
            bor[j] = rho[nx-1]
        elif j == 1:
            bol[j] = u[0]*rho[0]
            bor[j] = u[nx-1]*rho[nx-1]
        elif j == 2:
            bol[j] = e[0]
            bor[j] = e[nx-1]

    qf = [[0.0] * 3 for i in [1] * nx]  # 基本量
    Qc = [[0.0] * 3 for i in [1] * nx]  # 保存量
    for i in range(nx):
        for j in range(3):
            if j == 0:
                qf[i][j] = u[i]
                Qc[i][j] = rho[i]
            elif j == 1:
                qf[i][j] = rho[i]
                Qc[i][j] = u[i]*rho[i]
            elif j == 2:
                qf[i][j] = p[i]
                Qc[i][j] = e[i]

# --------------------------
# -- function cal_Q    --
# --------------------------
def cal_Q():
    """
    時間ステップを進める
    Qc=inner_ite(Qc)
    """
    global Qc

    Qc=inner_ite(Qc)

# --------------------------
# -- function bound       --
# --------------------------
def bound(lQc):  # 境界の計算
    """
    境界条件の設定
    """
    for i in range(3):
        lQc[0][i] = 2*bol[i]-lQc[1][i]  # 左端境界の計算
        lQc[nx-1][i] = lQc[nx-2][i]     # 右端の計算

    return lQc

# ----------------------------
# -- function cal_for_lusgs --
# ----------------------------
def cal_for_lusgs(lQc, lqf):
    """
    LU-SGSで行うLDU分解に向けて、ヤコビアン行列の近似(A+とA-)とbetaの計算を行う
    :param lQc: 保存変数
    :param lqf: 基本変数
    """
    global Amatrix_plus,Amatrix_minus,beta_sigma
    Amatrix_plus = np.empty((nx,3,3))
    Amatrix_minus = np.empty((nx,3,3))
    beta_sigma = np.empty(nx)
    beta = 1.1

    I = np.eye(3)

    for i in range(0, nx):
        H = (lQc[i][2]+lqf[i][2])/lQc[i][0]
        u = lqf[i][0]
        c = np.sqrt((gamma-1)*(H-0.5*u**2))

        sigma = abs(u) + c
        # iセルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(lQc[i],lqf[i])
        A_matrix = R @ Gam @ R_inv

        temp = I*beta*sigma

        Amatrix_plus[i] = 0.5*(A_matrix+temp)
        Amatrix_minus[i] = 0.5*(A_matrix-temp)
        beta_sigma[i] = beta*sigma

# --------------------------
# -- function inner_ite   --
# --------------------------
def inner_ite(lQc):
    """
    内部反復
    """
    global RHS

    delta_Q = np.zeros([nx, 3])
    delta_Q_temp = np.zeros([nx, 3])

    Qcn = np.array(copy.deepcopy(lQc))
    Qcm = np.array(copy.deepcopy(lQc))

    # inner iteration
    for ttt in range(10):

        qfm=Qctoqf(Qcm)

        cal_RHS(Qcm)  # Eの差分を算出

        cal_for_lusgs(Qcm,qfm)

        delta_Q = np.zeros([nx, 3])
        delta_Q_temp = np.zeros([nx, 3])

        lo_R = np.array(RHS)  # E_j+1/2 - E_j-1/2のarray

        sum_b = np.zeros(3)
        for i in range(lbound, nx - lbound):
            sum_b += abs(lo_R[i])  # E_j+1/2- E_j-1/2の絶対値の総和, 3*1
        print(sum_b)
        ite = 0
        con = 0  # 計算を続行するか否か
        # lusgs loop
        while con == 0:
            delta_Q_temp = copy.deepcopy(delta_Q)

            L = np.zeros([nx, 3])
            D = np.zeros(nx)
            U = np.zeros([nx, 3])
            for i in range(nx):
                L[i] = dt * (Amatrix_plus[i] @ delta_Q[i]) / 2
                D[i] = dx + dx + dt * beta_sigma[i]
                U[i] = dt * (Amatrix_minus[i] @ delta_Q[i]) / 2


            RHS = np.zeros([nx, 3])
            # RHS
            for i in range(lbound,nx-lbound):
                RHS[i] = -(Qcm[i] - Qcn[i]) * dx - dt * lo_R[i]

            # (D+L)Q=Q sweep1
            for i in range(lbound,nx-lbound):
                delta_Q[i] = (L[i-1] + RHS[i]) / D[i]

            # (D+U)Q=DQ sweep2,3
            for i in range(nx-lbound-1, 0, -1):
                delta_Q[i] = delta_Q[i] - U[i+1]/ D[i]

            # 収束判定
            if (ite + 1) % 100 == 0:
                sum_b_Ax = np.zeros(3)

                for i in range(lbound,nx-lbound):
                    for j in range(3):
                        sum_b_Ax += abs(delta_Q[i]-delta_Q_temp[i])

                norm2d = [0.0] * 3

                for i in range(3):
                    norm2d[i]=sum_b_Ax[i]/sum_b[i]  # (delta_Q-delta_Q_temp)/(E_j+1/2- E_j-1/2)

                if norm2d[0] < norm_ok and norm2d[1] < norm_ok and norm2d[2] < norm_ok:
                    con=1

            ite += 1

        delta_Q.tolist()

        for i in range(lbound,nx-lbound):
            for j in range(3):
                Qcm[i][j] += delta_Q[i][j]

        Qcm=bound(Qcm)

    return Qcm

# --------------------------
# -- function inv_matrix  --
# --------------------------
def inv_matrix(D):
    """
    ガウスの消去法により逆行列を求める
    :param D: 逆行列にしたい行列
    :return D: 逆行列
    """
    for l in range(3):
        aa       = 1.0 / D[l][l]
        D[l][l] = 1.0

        for m in range(3):
            D[l][m] = D[l][m] * aa

        for m in range(3):
            if m != l:
                bb       = D[m][l]
                D[m][l] = 0.0
                for n in range(3):
                    D[m][n] = D[m][n] - bb*D[l][n]

    return D

# --------------------------
# -- function cal_RHS     --
# --------------------------
def cal_RHS(lQc):  # 境界フラックスの計算
    """
    右辺の計算
    fvs法によってEを定義し、E_j+1/2-E_j-1/2を計算する
    :param lQc:
    """
    global RHS

    RHS = np.zeros([nx, 3])
    fvs(lQc)                             # FVS法によるフラックスの作成

    for i in range(1, nx-1):
        RHS[i] = Fplus[i] - Fplus[i-1]

    RHS.tolist()

# --------------------------
# -- function fvs         --
# --------------------------
def fvs(lQc):
    """
    FVS法によるフラックスの計算
    Fplus = E_j+1/2を定義、(セル1と2の境界をFplus[1]に格納)
    :param lQc:
    """
    global Fplus

    Fplus = np.zeros([nx+1,3])
    muscl(lQc)

    for i in range(0, nx-1):
        # i+1/2セルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(QcL[i],qfL[i])
        Ap = R @ (Gam + Gam_abs) @ R_inv  # 固有値が正のものを計算

        # i+1/2セルにおけるR,R^-1,Λ,|Λ|
        R, R_inv, Gam, Gam_abs = A_pm(QcR[i],qfR[i])
        Am = R @ (Gam - Gam_abs) @ R_inv  # 固有値が負のものを計算

        Fplus[i] = 0.5*(Ap @ QcL[i] + Am @ QcR[i])   # フラックスを計算

# --------------------------
# -- function A_pm       --
# --------------------------
def A_pm(lQc,lqf):
    """
    ヤコビアン行列の固有行列、固有値etc計算
    :return R, Rinv, lambda, abs(lambda):
    """
    H = (lQc[2]+lqf[2])/lQc[0]  # エンタルピー
    u = lqf[0]
    c = np.sqrt((gamma-1)*(H-0.5*u**2))
    b_para = (gamma-1)/c**2
    a_para = 0.5*b_para*u**2

    R = np.array([[1.0, 1.0, 1.0],
                  [u-c, u, u+c],
                  [H-u*c, 0.5*u**2, H+u*c]])
    R_inv = np.array([[0.5*(a_para+u/c), 0.5*(-b_para*u-1/c), 0.5*b_para], [
                        1-a_para, b_para*u, -b_para], [0.5*(a_para-u/c), 0.5*(-b_para*u+1/c), 0.5*b_para]])
    Gam = np.array([[(u-c), 0.0, 0.0], [0.0, u, 0.0], [0.0, 0.0, (u+c)]])
    Gam_abs = np.array(
        [[abs(u-c), 0.0, 0.0], [0.0, abs(u), 0.0], [0.0, 0.0, abs(u+c)]])

    return R, R_inv, Gam, Gam_abs

# --------------------------
# -- function yacobi_A    --
# --------------------------
def yacobi_A(lQc,lqf):
    """
    ヤコビ行列を計算し+と-の二種を返す
    :return A+, A-:
    """

    yacobiAp=[0.0] * nx
    yacobiAm=[0.0] * nx
    for i in range(nx):
        R, R_inv, Gam, Gam_abs = A_pm(lQc[i],lqf[i])
        yacobiAp[i] = np.dot((np.dot(R, Gam+Gam_abs)), R_inv)
        yacobiAm[i] = np.dot((np.dot(R, Gam-Gam_abs)), R_inv)

    return yacobiAp,yacobiAm

# -------------------------------
# -- function muscl            --
# -------------------------------
def muscl(lQc):
    """
    muscl法による補間
    QcLやQcRを更新する
    :param lQc:
    """
    global qf,qfL,qfR,QcL,QcR
    # 1と2の間を1に収納

    lqf=Qctoqf(lQc)

    qfL=[[0.0] * 3 for i in [1] * (nx+1)]
    qfR=[[0.0] * 3 for i in [1] * (nx+1)]

    for i in range(1,nx-2):
        for j in range(3):
            dplus_j=lqf[i+1][j]-lqf[i][j]
            dminus_j=lqf[i][j]-lqf[i-1][j]
            dplus_jp=lqf[i+2][j]-lqf[i+1][j]
            dminus_jp=lqf[i+1][j]-lqf[i][j]  # =dplus_j

            qfL[i][j]=lqf[i][j]  +1/4*((1-k_muscl)*minmod(dminus_j,dplus_j,b_muscl)+(1+k_muscl)*minmod(dplus_j,dminus_j,b_muscl))
            qfR[i][j]=lqf[i+1][j]-1/4*((1-k_muscl)*minmod(dplus_jp,dminus_jp,b_muscl)+(1+k_muscl)*minmod(dminus_jp,dplus_jp,b_muscl))

    # 境界内側用
    for j in range(3):
        dplus_jp=lqf[2][j]-lqf[1][j]
        dminus_jp=lqf[1][j]-lqf[0][j]
        qfR[0][j]=lqf[1][j]-1/4*((1-k_muscl)*minmod(dplus_jp,dminus_jp,b_muscl)+(1+k_muscl)*minmod(dminus_jp,dplus_jp,b_muscl))

        dplus_j=lqf[nx-1][j]-lqf[nx-2][j]
        dminus_j=lqf[nx-2][j]-lqf[nx-3][j]
        qfL[nx-2][j]=lqf[nx-2][j]+1/4*((1-k_muscl)*minmod(dminus_j,dplus_j,b_muscl)+(1+k_muscl)*minmod(dplus_j,dminus_j,b_muscl))


    QcL=qftoQc(qfL)
    QcR=qftoQc(qfR)

    # 境界外側用(境界は風上)
    qfL[0]=lqf[0][:]
    QcL[0]=lQc[0][:]
    qfR[nx-2]=lqf[nx-1][:]
    QcR[nx-2]=lQc[nx-1][:]

# -------------------------------
# -- function minmod           --
# -------------------------------
def minmod(x,y,b):
    """
    流速制限関数minmod
    """
    ans=np.sign(x)*max(0,min(abs(x),np.sign(x)*y*b))

    return ans

# -------------------------------
# -- function qftoQc           --
# -------------------------------
def qftoQc(qf):
    """
    基本量から保存量に変換
    """
    lo_Qc=[[0.0] * 3 for i in [1] * nx]
    for i in range(nx):
        for j in range(3):
            if j ==0:
                lo_Qc[i][j]=qf[i][1]
            elif j ==1:
                lo_Qc[i][j]=qf[i][1]*qf[i][0]
            elif j ==2:
                lo_Qc[i][j]=(qf[i][2]/(gamma-1)+1.0/2.0*qf[i][1]*(qf[i][0]**2))

    return lo_Qc

# -------------------------------
# -- function Qctoqf           --
# -------------------------------
def Qctoqf(Qc):
    """
    保存量から基本量に変換
    """
    lo_qf=[[0.0] * 3 for i in [1] * nx]
    for i in range(nx):
        for j in range(3):
            if j ==0:
                lo_qf[i][j]=Qc[i][1]/Qc[i][0]
            elif j ==1:
                lo_qf[i][j]=Qc[i][0]
            elif j ==2:
                lo_qf[i][j]=(gamma-1)*(Qc[i][2]-1.0/2.0*Qc[i][0]*((Qc[i][1]/Qc[i][0])**2))
    return lo_qf


# -------------------------------
# -- function output_q         --
# -------------------------------
def output_q(f_name):  # テキスト形式で出力
    """
    x,rho,u,pの出力
    """
    outlist=["x[m] u[m/s] rho[kg/m3] p[Pa]"]  # 出力するものの名前
    for i in range(len(qf)):
        outlist.append(str(x[i])+" "+str(qf[i][0])+" "+str(qf[i][1])+" "+str(qf[i][2]))
    outlist='\n'.join(outlist)

    with open(dir_name+"/"+f_name,'wt') as f:
        f.write(outlist)

if __name__ == '__main__':
    # preparetion
    os.makedirs(dir_name, exist_ok=True)
    setup()

    for k in range(nstep):
        print(k)
        cal_Q()
        qf=Qctoqf(Qc)

        #output_q(out_name_front+'{:0=4}'.format(int(k*dt*1000))+out_name_back)
        output_q(out_name_front+str(int(k*dt*1000))+out_name_back)
