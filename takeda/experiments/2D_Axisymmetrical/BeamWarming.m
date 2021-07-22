function [dq] = BeamWarming(q,dq_ex,dt,dxg,drg,M,N)
    %   A genuine 2d ChakaravarthyOsher Riemnan solver for Euler Equations using a Monotonic
    %   Upstreat Centered Scheme for Conservation Laws (MUSCL).
    %
    %   e.g. where: limiter='MC'; fluxMethod='HLLE1d';
    %
    %   Flux at j+1/2
    %
    %     j+1/2         Cell's grid:
    %   | wL|   |
    %   |  /|wR |           1   2   3   4        N-2 N-1  N
    %   | / |\  |   {x=0} |-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-| {x=L}
    %   |/  | \ |             1   2   3   4        N-2 N-1
    %   |   |  \|
    %   |   |   |       NC: Here cells 1 and N are ghost cells
    %     j  j+1            faces 1 and N-1, are the real boundary faces.
    %
    %   q = cat(3, rho, rho.*u, rho.*v, rho.*E);
    %   F = cat(3, rho.*u, rho.*u.^2+p, rho.*u.*v, u.*(rho.*E+p));
    %   G = cat(3, rho.*v, rho.*u.*v, rho.*v.^2+p, v.*(rho.*E+p));
    %
    %   Written by Ryota Takeda, JPN, 06.30.2021.
    global Param
    dq = zeros(M,N,4);
    gm = Param.GC.gamma;
    %%%%%%%%%%%%%%%%%
    % set constants %
    %%%%%%%%%%%%%%%%%
    theta_1 = Param.hp.theta_1;
    theta_2 = Param.hp.theta_2;

    % Build cells and Compute Jacobian and x,r
    % メモリ節約のため何回も参照する or 差分に使用するもの以外は一時変数にする
    cell(M,N).all = M*N;
    x=zeros(M,N); r=zeros(M,N); % cartesian grid
    for j = 1:N
        rg = j*drg; % 一般化されたr. cartesianでやる場合はrと一致
        for i = 1:N
            xg = i*dxg; %一般化されたx. cartesianでやる場合はxと一致
            [x,r] = getCartesian(xg,rg);
            r(i,j) = r; x(i,j) = x; % cartesian grid
            cell(i,j).q   = [q(i,j,1);q(i,j,2);q(i,j,3);q(i,j,4)];
            % メモリ節約のためrhs~dQは全てまとめて更新することにする。
            cell(i,j).res = [dq_ex(i,j,1);dq_ex(i,j,2);dq_ex(i,j,3);dq_ex(i,j,4)] * theta_2/(1+theta_2); %最初にdQ^n-1を入れておく
            J = Param.hp.Jacobian([x r]);
            cell(i,j).detJ   = det(J);
            cell(i,j).J  = J; % xy座標から一般化座標に変換するヤコビアン
        end
    end

    % Build Faces(cellの中間値)
    face(M-1,N-1).all = (M-1)*(N-1);
    for i = 1:M-1
        for j = 1:N-1
            face(i,j).G1d_E = zeros(4,1); %E_{i+1/2,j}
            face(i,j).lambda_E = zeros(4,1);
            face(i,j).T_E = zeros(4,4);
            face(i,j).G1d_F = zeros(4,1); %F_{i,j+1/2}
            face(i,j).lambda_F = zeros(4,1);
            face(i,j).T_F = zeros(4,4);
            face(i,j).dQ2 = zeros(4,1);
            face(i,j).dQ4 = zeros(4,1);
        end
    end

    % Compute fluxes across cells
    for i = 1:M-1     % all internal faces
        for j = 1:N-1 % all internal faces
            qSW = cell( i , j ).q;
            qSE = cell( i ,j+1).q;
            qNW = cell(i+1, j ).q;
            % compute HLLE1d flux
            [face(i,j).lambda_E,face(i,j).T_E,face(i,j).E_] = CalG1d(qSW,qSE,cell(i,j).J(1,:),'E');   % G1d_{i+1/2,  j  }
            [face(i,j).lambda_F,face(i,j).T_F,face(i,j).G1d_F] = CalG1d(qSW,qNW,cell(i,j).J(2,:),'F');   % G1d_{  i  ,j+1/2}
        end
    end

    %%%%%%%%%%%%%
    % Residuals %
    %%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate reft hand side %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 3:M-2     % internal cells
        for j = 3:N-2 % internal faces
            J = cell(i,j).J;

            % calculate flux E
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            % セル境界の物理量は空間一次精度で良い
            qLL = cell(i-2,j).q;
            qL  = cell(i-1,j).q;
            q   = cell( i,j ).q;
            qR  = cell(i+1,j).q;
            qRR = cell(i+2,j).q;
            % compute flux at  i+1/2 using
            dE = ChakravarthyOsher(qLL,qL,q,qR,qRR,'xg',cell(i,j).Jn,cell(i-1,j).Jn); % E_{i+1/2,j}-E_{i-1/2,j}

            % calculate flux F
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            qLL = cell(i,j-2).q;
            qL  = cell(i,j-1).q;
            qR  = cell(i,j+1).q;
            qRR = cell(i,j+2).q;
            % compute flux at j+1/2 using
            dF = ChakravarthyOsher(qLL,qL,q,qR,qRR,'rg',cell(i,j).Jn,cell(i,j-1).Jn); % F_{i,j+1/2}-F_{i,j-1/2}

            % compute H
            % Heat Input
            w = getHeatingSource(xg, cosine);
            rho = q(1);
            u   = q(2);
            v   = q(3);
            E   = q(4)/q(1);
            p = (gm-1)*(q(4) - rho*(u^2+v^2)/2);
            H = [rho*v; rho*u*v; rho*v^2; w + (rho*E+p)*v];
            % contributions to the residual of cell (i,j) and cells around it
            cell(i,j).res = cell(i,j).res - dE*dt/dxg/(1+theta_1)...
                                          - dF*dt/drg/(1+theta_1)...
                                          -  H*dt/cell(i,j).r/(1+theta_1);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % Inverting a matrix %
    %%%%%%%%%%%%%%%%%%%%%%
    % step1
    for i = 3:M-2     % internal cells
        for j = 3:N-2 % internal faces
            cell(i,j).res = face(i,j).T_E * cell(i,j).res;
        end
    end

    % step2
    for i = 3:M-2     % internal cells
        for j = 3:N-2 % internal faces
            dQ2L  = [cell(i-1,j).dQ2];
            dQ2R  = [cell(i+1,j).dQ2];
            cell(i,j).dQ2 = (eye(4)+theta_1/(1+theta_2)*dt/2/dx*diag(abs(lambdax)))\...
                            (cell(i,j).res+theta_1/(1+theta_2)*dt/2/dx*(lambdaxR*dQ2L-lambdaxL*dQ2R));

            % step3
            T_rT_x = [[1,              0,       0,                0];...
                      [0,        J/2/g22, g12/g22,           -J/g22];...
                      [0, G0/2-g12/2/g22,   J/g22,   G0/2+g12/2/g22];...
                      [0, G0/2+g12/2/g22,  -J/g22,   G0/2-g12/2/g22]];
            cell(i,j).res = T_rT_x * cell(i,j).dQ2;

            % step4
            dQ4L  = [cell(i,j-1).dQ4];
            dQ4R  = [cell(i,j+1).dQ4];
            cell(i,j).dQ4 = (eye(4)+theta_1/(1+theta_2)*dt/2/dr*diag(abs(lambdar)))\...
                            (cell(i,j).res+theta_1/(1+theta_2)*dt/2/dr*(lambdarR*dQ4L-lambdarL*dQ4R));

            % step5
            cell(i,j).res = T_r * cell(i,j).dQ4;

            % step6
            % ソースタームのヤコビ行列
            omega = gm*E(i,j);
            phi2 = (gm-1)/2*(u(i,j)^2+v(i,j)^2);
            C = [[                    0,                     0,                          1,         0];...
                 [       -u(i,j)*v(i,j),                v(i,j),                     u(i,j),         0];...
                 [            -v(i,j)^2,                     0,                   2*v(i,j),         0];...
                 [v(i,j)*(2*phi2-omega), -(gm-1)*u(i,j)*v(i,j), omega-phi2-(gm-1)*v(i,j)^2, gm*v(i,j)]];
            cell(i,j).res = (eye(4)+theta_1/(1+theta_2)*dt/cell(i,j).r*C) \ cell(i,j).res;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions %
%%%%%%%%%%%%%%%%%%%%%%%

function CHOS = ChakravarthyOsher(qLL,qL,q,qR,qRR,G1dL,G1dR)
    %{
    Obtain higher order accuracy numerical fluxes
    by interpolating the numerical fluxes calculated with spatial first order accuracy
    %}
    global Param
    % compute flux term
    bcomp   = Param.hp.bcomp;
    kappa   = Param.hp.kappa;
    limiter = Param.hp.limiter;
    function G2d = CalG2d(qL,q,qR,qRR,G1d)
        lambdaM = (lambda - abs(lambda))/2;
        lambdaP = (lambda + abs(lambda))/2;
        dWL = T\(q-qL);   % dW_{k-1/2}
        dWC = T\(qR-qL);  % dW_{k+1/2}
        dWR = T\(qRR-qR); % dW_{k+3/2}
        dG = zeros(4,1);
        for i = 1: 4
            dG(i) = (1-kappa)/4*(Limiter([lambdaP(i)*dWL(i),lambdaP(i)*dWC(i)],bcomp,limiter)*T(:,i))...
                  + (1+kappa)/4*(Limiter([lambdaP(i)*dWC(i),lambdaP(i)*dWL(i)],bcomp,limiter)*T(:,i))...
                  - (1-kappa)/4*(Limiter([lambdaM(i)*dWR(i),lambdaM(i)*dWC(i)],bcomp,limiter)*T(:,i))...
                  - (1+kappa)/4*(Limiter([lambdaM(i)*dWC(i),lambda<(i)*dWR(i)],bcomp,limiter)*T(:,i))...
                  + dG;
        end
        G2d = G1d + dG;
    end
    G2dR = CalG2d(qL,q,qR,qRR,G1dR);
    G2dL = CalG2d(qLL,qL,q,qR,G1dL);
    CHOS = G2dR-G2dL;
end



function [lambda,T,flux1d] = getRoeAverage(qL,qR,normal,J,axi)
    %{
    calculate Roe Average and Numerical flux calculated with spatial first-order accuracy.
    param:
    return: eigenvalue_matrix, eigenmatrix, flux
    %}
    global Param
    % normalにはcase Eの場合xi_x,xi_r
    nx = normal(1);
    nr = normal(2);

    g = nx^2+nr^2;

    gm1 = gm-1;
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*gx+vL*gr;
    pL = (gm1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    aL = sqrt(gm*pL/rL);
    HL = ( qL(4) + pL ) / rL;

    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*nr;
    pR = (gm1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gm*pR/rR);
    HR = (qR(4) + pR) / rR;

    % First compute the Roe Averages
    RT = sqrt(rR/rL); % r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt(gm1*(H-(u^2+v^2)/2));
    vn = u*nx+v*nr; % 反変速度

    % Riemann solver
    N = [[1 0 0 0];...
         [-u 1 0 0];...
         [-v 0 1 0];...
         [(u^2+v^2)/2 -u -v 1]*gm1];
    C = zeros(4,4);
    switch axi
        case 'E'
            lambda = [vn; vn+a*sqrt(g); vn; vn-a*sqrt(g)];
            C = [[1   0   0    -1/a^2]...
                 [0  nx  nr sqrt(g)/a]...
                 [0 -nr  nx         0]...
                 [0 -nx -nr sqrt(g)/a]];
        case 'F'
            lambda = [vn; vn; vn+a*sqrt(g); vn-a*sqrt(g)];
            C = [[1   0   0    -1/a^2]...
                 [0  nr -nx sqrt(g)/a]...
                 [0  nx  nr         0]...
                 [0 -nx -nr sqrt(g)/a]];
    end
    lambdaL = (lambda - abs(lambda))/2;
    lambdaR = (lambda + abs(lambda))/2;

    T = C*N;

    % % Left and Right fluxes
    flux_L=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*nr; rL*vnL*HL];
    flux_R=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*nr; rR*vnR*HR];

    dflux_L = T*diag(lambdaL)/T*(qR-qL);
    dflux_R = T*diag(lambdaR)/T*(qR-qL);

    flux1d = 1/2*(flux_L+flux_R - (dflux_R-dflux_L));
end

function HLLE = HLLE1Dflux(qL,qR,flux,normal)
    % Compute HLLE flux
    % normal vectors
    nx = normal(1);
    ny = normal(2);

    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    pL = (gm-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    aL = sqrt(gm*pL/rL);
    HL = ( qL(4) + pL ) / rL;

    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    pR = (gm-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gm*pR/rR);
    HR = ( qR(4) + pR ) / rR;

    % First compute the Roe Averages
    RT = sqrt(rR/rL); % r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gm-1)*(H-(u^2+v^2)/2) );
    vn = u*nx+v*ny;

    % % Wave speed estimates
    SLm = min([ vnL-aL, vn-a, 0]);
    SRp = max([ vnR+aR, vn+a, 0]);

    % % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];

    % Compute the HLL flux.
    HLLE = ( SRp*FL - SLm*FR + SLm*SRp*(qR-qL) )/(SRp-SLm);
end

function [fOO,gOO] = HLLE2Dflux(qSW,qSE,qNW,qNE)
    % Compute HLLE flux
    global Param

    % West state
    rSW = qSW(1);
    uSW = qSW(2)/rSW;
    vSW = qSW(3)/rSW;
    pSW = (gm-1)*( qSW(4) - rSW*(uSW^2+vSW^2)/2 );
    aSW = sqrt(gm*pSW/rSW);
    HSW = ( qSW(4) + pSW ) / rSW;

    % East state
    rSE = qSE(1);
    uSE = qSE(2)/rSE;
    vSE = qSE(3)/rSE;
    pSE = (gm-1)*( qSE(4) - rSE*(uSE^2+vSE^2)/2 );
    aSE = sqrt(gm*pSE/rSE);
    HSE = ( qSE(4) + pSE ) / rSE;

    % South state
    rNW = qNW(1);
    uNW = qNW(2)/rNW;
    vNW = qNW(3)/rNW;
    pNW = (gm-1)*( qNW(4) - rNW*(uNW^2+vNW^2)/2 );
    aNW = sqrt(gm*pNW/rNW);
    HNW = ( qNW(4) + pNW ) / rNW;

    % North state
    rNE = qNE(1);
    uNE = qNE(2)/rNE;
    vNE = qNE(3)/rNE;
    pNE = (gm-1)*( qNE(4) - rNE*(uNE^2+vNE^2)/2 );
    aNE = sqrt(gm*pNE/rNE);
    HNE = ( qNE(4) + pNE ) / rNE;




    % Compute Roe Averages - SW to SE
    rSroe = sqrt(rSE/rSW);
    uSroe = (uSW+rSroe*uSE)/(1+rSroe);
    vSroe = (vSW+rSroe*vSE)/(1+rSroe);
    HSroe = (HSW+rSroe*HSE)/(1+rSroe);
    aSroe = sqrt( (gm-1)*(HSroe-0.5*(uSroe^2+vSroe^2)) );

    % Compute Roe Averages - NW to NE
    rNroe = sqrt(rNE/rNW);
    uNroe = (uNW+rNroe*uNE)/(1+rNroe);
    vNroe = (vNW+rNroe*vNE)/(1+rNroe);
    HNroe = (HNW+rNroe*HNE)/(1+rNroe);
    aNroe = sqrt( (gm-1)*(HNroe-0.5*(uNroe^2+vNroe^2)) );

    % Compute Roe Averages - SW to NW
    rWroe = sqrt(rNW/rSW);
    uWroe = (uSW+rWroe*uNW)/(1+rWroe);
    vWroe = (vSW+rWroe*vNW)/(1+rWroe);
    HWroe = (HSW+rWroe*HNW)/(1+rWroe);
    aWroe = sqrt( (gm-1)*(HWroe-0.5*(uWroe^2+vWroe^2)) );

    % Compute Roe Averages - SE to NE
    rEroe = sqrt(rNE/rSE);
    uEroe = (uSE+rEroe*uNE)/(1+rEroe);
    vEroe = (vSE+rEroe*vNE)/(1+rEroe);
    HEroe = (HSE+rEroe*HNE)/(1+rEroe);
    aEroe = sqrt( (gm-1)*(HEroe-0.5*(uEroe^2+vEroe^2)) );




    % Wave speed estimates in the S
    sSW = min([ uSW-aSW, uSW+aSW, uSroe-aSroe, uSroe+aSroe ]);
    sSE = max([ uSE-aSE, uSE+aSE, uSroe-aSroe, uSroe+aSroe ]);

    % Wave speed estimates in the N
    sNW = min([ uNW-aNW, uNW+aNW, uNroe-aNroe, uNroe+aNroe ]);
    sNE = max([ uNE-aNE, uNE+aNE, uNroe-aNroe, uNroe+aNroe ]);

    % Wave speed estimates in the W
    sWS = min([ vSW-aSW, vSW+aSW, vWroe-aWroe, vWroe+aWroe ]);
    sWN = max([ vNW-aNW, vNW+aNW, vWroe-aWroe, vWroe+aWroe ]);

    % Wave speed estimates in the E
    sES = min([ vSE-aSE, vSE+aSE, vEroe-aEroe, vEroe+aEroe ]);
    sEN = max([ vNE-aNE, vNE+aNE, vEroe-aEroe, vEroe+aEroe ]);




    % Compute fluxes
    fSW = [rSW*uSW; rSW*uSW*uSW + pSW; rSW*vSW*uSW; rSW*uSW*HSW];
    fSE = [rSE*uSE; rSE*uSE*uSE + pSE; rSE*vSE*uSE; rSE*uSE*HSE];
    fNW = [rNW*uNW; rNW*uNW*uNW + pNW; rNW*vNW*uNW; rNW*uNW*HNW];
    fNE = [rNE*uNE; rNE*uNE*uNE + pNE; rNE*vNE*uNE; rNE*uNE*HNE];

    gSW = [rSW*vSW; rSW*vSW*uSW; rSW*vSW*vSW + pSW; rSW*vSW*HSW];
    gSE = [rSE*vSE; rSE*vSE*uSE; rSE*vSE*vSE + pSE; rSE*vSE*HSE];
    gNW = [rNW*vNW; rNW*vNW*uNW; rNW*vNW*vNW + pNW; rNW*vNW*HNW];
    gNE = [rNE*vNE; rNE*vNE*uNE; rNE*vNE*vNE + pNE; rNE*vNE*HNE];

    % Compute the intermediate states
    qSO = ( sSE*qSE - sSW*qSW + fSW-fSE )/(sSE-sSW);
    qNO = ( sNE*qNE - sNW*qNW + fNW-fNE )/(sNE-sNW);
    qOW = ( sWN*qNW - sWS*qSW + gSW-gNW )/(sWN-sWS);
    qOE = ( sEN*qNE - sES*qSE + gSE-gNE )/(sEN-sES);

    % Compute the intermediate states fluxes (normal HLLE 1d fluxes)
    fSO = ( sSE*fSW - sSW*fSE + sSW*sSE*(qSE-qSW) )/(sSE-sSW);
    fNO = ( sNE*fNW - sNW*fNE + sNW*sNE*(qNE-qNW) )/(sNE-sNW);
    gOW = ( sWN*gSW - sWS*gNW + sWS*sWN*(qNW-qSW) )/(sWN-sWS);
    gOE = ( sEN*gSE - sES*gNE + sES*sEN*(qNE-qSE) )/(sEN-sES);

    % Compute the transverse intermediate fluxes (Balsara's solution)
    fOW = [qOW(2);gOW(3)+(qOW(2)^2-qOW(3)^2)/qOW(1);qOW(3)*qOW(2)/qOW(1);qOW(2)*gOW(4)/qOW(3)];
    fOE = [qOE(2);gOE(3)+(qOE(2)^2-qOE(3)^2)/qOE(1);qOE(3)*qOE(2)/qOE(1);qOE(2)*gOE(4)/qOE(3)];
    gSO = [qSO(3);qSO(2)*qSO(3)/qSO(1);fSO(2)+(qSO(3)^2-qSO(2)^2)/qSO(1);qSO(3)*fSO(4)/qSO(2)];
    gNO = [qNO(3);qNO(2)*qNO(3)/qNO(1);fNO(2)+(qNO(3)^2-qNO(2)^2)/qNO(1);qNO(3)*fNO(4)/qNO(2)];


    % Strongly Interacting state q**
    qOO = 1/((sNE-sSW)*(sWN-sES)+(sEN-sWS)*(sSE-sNW))*( ...
            (sWN*sNE+sSE*sEN)*qNE - (sEN*sNW+sSW*sWN)*qNW + ...
            (sES*sSW+sNW*sWN)*qSW - (sWS*sSE+sNE*sES)*qSE ...
        - sWN*fNE+sEN*fNW - sES*fSW+sWS*fSE - (sEN-sES)*fOE+(sWN-sWS)*fOW ...
        - sSE*gNE+sSW*gNW - sNW*gSW+sNE*gSE - (sNE-sNW)*gNO+(sSE-sSW)*gSO );

    % Compute fluxes of the strongly interacting state:
    % Precompute deltas
    dq1 = sNW*sEN-sWN*sNE; df1 = sWN-sEN; dg1 = sNE-sNW;
    dq2 = sSW*sWN-sWS*sNW; df2 = sWS-sWN; dg2 = sNW-sSW;
    dq3 = sSE*sWS-sES*sSW; df3 = sES-sWS; dg3 = sSW-sSE;
    dq4 = sNE*sES-sEN*sSE; df4 = sEN-sES; dg4 = sSE-sNE;

    % Using LSQ
    b1 = dq1*(qNO-qOO) + df1*fNO + dg1*gNO;
    b2 = dq2*(qOW-qOO) + df2*fOW + dg2*gOW;
    b3 = dq3*(qSO-qOO) + df3*fSO + dg3*gSO;
    b4 = dq4*(qOE-qOO) + df4*fOE + dg4*gOE;

    % k-weights
    k11 = df1*(dg2^2+dg3^2+dg4^2) - dg1*(df2*dg2+df3*dg3+df4*dg4);
    k12 = df2*(dg1^2+dg3^2+dg4^2) - dg2*(df1*dg1+df3*dg3+df4*dg4);
    k13 = df3*(dg1^2+dg2^2+dg4^2) - dg3*(df1*dg1+df2*dg2+df4*dg4);
    k14 = df4*(dg1^2+dg2^2+dg3^2) - dg4*(df1*dg1+df2*dg2+df3*dg3);
    k21 = dg1*(df2^2+df3^2+df4^2) - df1*(df2*dg2+df3*dg3+df4*dg4);
    k22 = dg2*(df1^2+df3^2+df4^2) - df2*(df1*dg1+df3*dg3+df4*dg4);
    k23 = dg3*(df1^2+df2^2+df4^2) - df3*(df1*dg1+df2*dg2+df4*dg4);
    k24 = dg4*(df1^2+df2^2+df3^2) - df4*(df1*dg1+df2*dg2+df3*dg3);

    %A = [df1,dg1;df2,dg2;df3,dg3;df4,dg4]; M=A'*A; detM=det(M);
    detM = (df1*dg2-df2*dg1)^2 + (df1*dg3-df3*dg1)^2 + (df2*dg4-df4*dg2)^2 + ...
        (df3*dg2-df2*dg3)^2 + (df4*dg1-df1*dg4)^2 + (df4*dg3-df3*dg4)^2 ; % verified!

    % compute fluxes of Strongly Interacting state f** and g**
    fOO = (k11*b1 + k12*b2 + k13*b3 + k14*b4)/detM;
    gOO = (k21*b1 + k22*b2 + k23*b3 + k24*b4)/detM;

end

function l = Limiter(v,b,limiter)
    x=v(1); y=v(2);
    switch limiter
        case 'mm'
            l = minmod([x b*y]);
        case 'sb'
            if y==0; l=0; else,l = superbee(sign(y)*x/abs(y))*abs(y); end
    end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3;
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    s = sum(sign(v))/numel(v);
    mm = s*min(abs(v(:)));
    %m=size(v,1); mm=zeros(size(v,1),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
end

function sb = superbee(v)
    % Using Roe's definition
    x=v(1); y=v(2);
    s = sum(sign(v))/numel(v);
    sb = s*max(min(abs([2*x y])), min(abs([x 2*y])));
end

function vl = vanLeer(da,db)
    % Van Leer Slope Limiter Function
    vl = 0; if bd~=0, r=da/db; vl=(r+abs(r))/(1+abs(r)); end
end
