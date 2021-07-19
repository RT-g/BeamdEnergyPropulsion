function [dq] = ChakaravarthyOsher(q,dq_ex,dt,dxg,drg,M,N)
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
    %%%%%%%%%%%%%%%%%
    % set constants %
    %%%%%%%%%%%%%%%%%
    theta_1 = Param.hp.theta_1;
    theta_2 = Param.hp.theta_2;


    % Normal unitary face vectors: (nx,ny)
    % normals = {[0,1], [1,0], [0,-1], [-1,0]}; % i.e.: [N, E, S, W]

    % set variables
    rho = q(:,:,1);
    u   = q(:,:,2)./rho;
    v   = q(:,:,3)./rho;
    E   = q(:,:,4)./rho;
    p   = (Param.GC.gamma-1)*rho.*(E-0.5*(u.^2+v.^2));

    % Build cells
    cell(M,N).all = M*N;
    for j = 1:N
        rg = j*drg; % 一般化されたr. cartesianでやる場合はrと一致
        [~, ~, B, r, cosine] = getSWVelocity(rg);
        for i = 1:N
            xg = i*dxg; %一般化されたx. cartesianでやる場合はxと一致
            if r < 1/sqrt(Param.LV.a_s)
                x = xg - sqrt((1-Param.LV.a_s*r^2)^(-2*B)-1);
            else
                disp(x)
                break
            end

            % Heat Input
            w = getHeatingSource(xg, cosine);

            cell(i,j).q   = [q(i,j,1);q(i,j,2);q(i,j,3);q(i,j,4)];
            cell(i,j).H   = [rho(i,j)*v(i,j); rho(i,j)*u(i,j)*v(i,j); rho(i,j)*v(i,j)^2; w + (rho(i,j)*E(i,j)+p(i,j))*v(i,j)];
            % メモリ節約のためrhs~dQは差分を取るdQ2,dQ4以外は全てまとめて更新することにする
            cell(i,j).res = [dq_ex(i,j,1);dq_ex(i,j,2);dq_ex(i,j,3);dq_ex(i,j,4)] * theta_2 / (1+theta_2); %最初にdQ_n-1を入れておく
            cell(i,j).dQ2 = zeros(4,1);
            cell(i,j).dQ4 = zeros(4,1);
            cell(i,j).r   = r; % cartesian grid
            cell(i,j).x   = x; % cartesian grid
            J = Param.hp.Jacobian();
            cell(i,j).J   = det(J);
            cell(i,j).Jn  = J; % xy座標から一般化座標に変換するヤコビアン
            cell(i,j).Jg  = inv(J); % 一般化座標からxy座標に変換するヤコビアン
        end
    end

    %%%%%%%%%%%%%
    % Residuals %
    %%%%%%%%%%%%%

    % Compute residuals x-direction
    for i = 3:M-2     % internal cells
        for j = 3:N-2 % internal faces
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculate reft hand side %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            J = cell(i,j).J;
            g11 = cell(i,j).Jg(:,1)' * cell(i,j).Jg(:,1);
            g22 = cell(i,j).Jg(:,2)' * cell(i,j).Jg(:,2);
            g12 = cell(i,j).Jg(:,1)' * cell(i,j).Jg(:,2);
            G0 = sqrt(g11/g22);

            % calculate flux E
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            % セル境界の物理量は空間一次精度で良い
            qLL = [cell(i-2,j).q];
            qL  = [cell(i-1,j).q];
            q   = [cell( i,j ).q];
            qR  = [cell(i+1,j).q];
            qRR = [cell(i+2,j).q];
            % compute flux at j+1/2 using
            [flux_E,lambdax,T_x] = ChakravarthyOsher(qLL,qL,q,qR,qRR,'xg',cell(i,j).Jn,cell(i-1,j).Jn); % E_{i+1/2,j}-E_{i-1/2,j}
            lambdaxR = diag((lambdax + abs(lambdax))/2);
            lambdaxL = diag((lambdax - abs(lambdax))/2);

            % calculate flux F
            % Left (inside) and Right (outside) extrapolated q-values at the boundaries
            qLL = [cell(i,j-2).q];
            qL  = [cell(i,j-1).q];
            q   = [cell( i,j ).q];
            qR  = [cell(i,j+1).q];
            qRR = [cell(i,j+2).q];
            % compute flux at i+1/2 using
            [flux_F,lambdar,T_r] = ChakravarthyOsher(qLL,qL,q,qR,qRR,'rg',cell(i,j).Jn,cell(i,j-1).Jn); % F_{i,j+1/2}-F_{i,j-1/2}
            lambdarR = diag((lambdar + abs(lambdar))/2);
            lambdarL = diag((lambdar - abs(lambdar))/2);

            % contributions to the residual of cell (i,j) and cells around it
            cell(i,j).res = cell(i,j).res - flux_E*dt/dx/(1+theta_1)...
                                          - flux_F*dt/dr/(1+theta_1)...
                                          - cell(i,j).H*dt/cell(i,j).r/(1+theta_1);

            %%%%%%%%%%%%%%%%%%%%%%
            % Inverting a matrix %
            %%%%%%%%%%%%%%%%%%%%%%

            % step1
            cell(i,j).res = T_x \ cell(i,j).res;

            % step2 自信なし
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
            dQ4L  = [cell(i-1,j).dQ4];
            dQ4R  = [cell(i+1,j).dQ4];
            cell(i,j).dQ4 = (eye(4)+theta_1/(1+theta_2)*dt/2/dx*diag(abs(lambdar)))\...
                            (cell(i,j).res+theta_1/(1+theta_2)*dt/2/dx*(lambdarR*dQ4L-lambdarL*dQ4R));

            % step5
            cell(i,j).res = T_r * cell(i,j).dQ4;

            % step6
            % ソースタームのヤコビ行列
            omega = Param.GC.gamma*E(i,j);
            phi2 = (Param.GC.gamma-1)/2*(u(i,j)^2+v(i,j)^2);
            C = [[0,0,1,0];...
                    [-u(i,j)*v(i,j),v(i,j),u(i,j),0];...
                    [-v(i,j)^2,0,2*v(i,j),0];...
                    [v(i,j)*(2*phi2-omega),-(Param.GC.gamma-1)*u(i,j)*v(i,j),omega-phi2-(Param.GC.gamma-1)*v(i,j)^2,Param.GC.gamma*v(i,j)]];
            cell(i,j).res = (eye(4)+theta_1/(1+theta_2)*dt/cell(i,j).r*C) \ cell(i,j).res;
        end
    end
end

    %%%%%%%%%%%%%%%%%%%%%%%
    % Auxiliary Functions %
    %%%%%%%%%%%%%%%%%%%%%%%

    function [CHOS,lambda,T] = ChakravarthyOsher(qLL,qL,q,qR,qRR,axi,JnR,JnL)
        global Param
        % compute flux term
        bcomp           = Param.hp.bcomp;
        kappa           = Param.hp.kappa;
        limiter         = Param.hp.limiter;
        [G2dR,lambda,T] = CalG2d(qL,q,qR,qRR,bcomp,kappa,limiter,axi,JnR);
        G2dL            = CalG2d(qLL,qL,q,qR,bcomp,kappa,limiter,axi,JnL);
        CHOS            = G2dR-G2dL;
    end

    function [G2d,lambda,T]= CalG2d(qL,q,qR,qRR,bcomp,kappa,limiter,axi,Jn)
        [lambda,T,G1d] = CalG1d(q,qR,axi,Jn);
        lambdaL = (lambda - abs(lambda))/2;
        lambdaR = (lambda + abs(lambda))/2;
        dWL = T\(q-qL);
        dWC = T\(qR-qL);
        dWR = T\(qRR-qR);
        dG = zeros(4,1);
        for i = 1:3
            dG = (1-kappa)/4*(Limiter([lambdaR(i)*dWL(i),lambdaR(i)*dWC(i)],bcomp,limiter)*T(:,i))...
               + (1+kappa)/4*(Limiter([lambdaR(i)*dWC(i),lambdaR(i)*dWL(i)],bcomp,limiter)*T(:,i))...
               - (1-kappa)/4*(Limiter([lambdaL(i)*dWR(i),lambdaL(i)*dWC(i)],bcomp,limiter)*T(:,i))...
               - (1+kappa)/4*(Limiter([lambdaL(i)*dWC(i),lambdaL(i)*dWR(i)],bcomp,limiter)*T(:,i))...
               + dG;
        end
        G2d = G1d + dG;
    end

    function [lambda,T,G1d] = CalG1d(qL,qR,axi,Jn)
        global Param
        J = det(Jn);
        r_xg = -Jn(2,1)/J;
        r_rg =  Jn(1,1)/J;
        x_xg =  Jn(2,2)/J;
        x_rg = -Jn(1,2)/J;

        switch axi
        case 'xg'
            % normal vectors
            gx = Jn(1,1);
            gr = Jn(1,2);
            g22 = gx^2+gr^2;
        case 'rg'
            % normal vectors
            gx = Jn(2,1);
            gr = Jn(2,2);
            g11 = gx^2+gr^2;
        end

        gm1 = Param.GC.gamma-1;
        % Left state
        rL = qL(1);
        uL = qL(2)/rL;
        vL = qL(3)/rL;
        vnL = uL*gx+vL*gr;
        pL = (gm1)*( qL(4) - rL*(uL^2+vL^2)/2 );
        aL = sqrt(Param.GC.gamma*pL/rL);
        HL = ( qL(4) + pL ) / rL;

        % Right state
        rR = qR(1);
        uR = qR(2)/rR;
        vR = qR(3)/rR;
        vnR = uR*gx+vR*gr;
        pR = (gm1)*( qR(4) - rR*(uR^2+vR^2)/2 );
        aR = sqrt(Param.GC.gamma*pR/rR);
        HR = ( qR(4) + pR ) / rR;

        % First compute the Roe Averages
        RT = sqrt(rR/rL); % r = RT*rL;
        u = (uL+RT*uR)/(1+RT);
        v = (vL+RT*vR)/(1+RT);
        H = ( HL+RT* HR)/(1+RT);
        a = sqrt( (gm1)*(H-(u^2+v^2)/2) );
        vn = u*gx+v*gr;

        % Riemann solver
        lambda = [vn; vn; vn+a*sqrt(gx^2+gr^2); vn-a*sqrt(gx^2+gr^2)];
        lambdaL = (lambda - abs(lambda))/2;
        lambdaR = (lambda + abs(lambda))/2;

        switch axi
        case 'xg'
            C = zeros(4,4);
            C(1,1) = 1;
            C(1,4) = -1/a^2;
            C(2,2) = r_rg;
            C(2,3) = -x_rg;
            C(2,4) = sqrt(g22)/a;
            C(3,2) = x_rg;
            C(3,3) = r_rg;
            C(4,2) = -r_rg;
            C(4,3) = x_rg;
            C(4,4) = sqrt(g22)/a;
        case 'rg'
            C = zeros(4,4);
            C(1,1) = 1;
            C(1,4) = -1/a^2;
            C(2,2) = x_xg;
            C(2,3) = r_xg;
            C(3,2) = -r_xg;
            C(3,3) = x_xg;
            C(3,4) = sqrt(g11)/a;
            C(4,2) = r_xg;
            C(4,3) = -x_xg;
            C(4,4) = sqrt(g11)/a;
        end

        N = [[1 0 0 0];...
             [-u 1 0 0];...
             [-v 0 1 0];...
             [(u^2+v^2)/2 -u -v 1]*gm1];



        T = C*N;

        % % Left and Right fluxes
        FL=[rL*vnL; rL*vnL*uL + pL*gx; rL*vnL*vL + pL*gr; rL*vnL*HL];
        FR=[rR*vnR; rR*vnR*uR + pR*gx; rR*vnR*vR + pR*gr; rR*vnR*HR];

        dFL = T*diag(lambdaL)/T*(qR-qL);
        dFR = T*diag(lambdaR)/T*(qR-qL);

        G1d = 1/2*(FL+FR - (dFR-dFL));
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
        pL = (Param.GC.gamma-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
        aL = sqrt(Param.GC.gamma*pL/rL);
        HL = ( qL(4) + pL ) / rL;

        % Right state
        rR = qR(1);
        uR = qR(2)/rR;
        vR = qR(3)/rR;
        vnR = uR*nx+vR*ny;
        pR = (Param.GC.gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
        aR = sqrt(Param.GC.gamma*pR/rR);
        HR = ( qR(4) + pR ) / rR;

        % First compute the Roe Averages
        RT = sqrt(rR/rL); % r = RT*rL;
        u = (uL+RT*uR)/(1+RT);
        v = (vL+RT*vR)/(1+RT);
        H = ( HL+RT* HR)/(1+RT);
        a = sqrt( (Param.GC.gamma-1)*(H-(u^2+v^2)/2) );
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
        pSW = (Param.GC.gamma-1)*( qSW(4) - rSW*(uSW^2+vSW^2)/2 );
        aSW = sqrt(Param.GC.gamma*pSW/rSW);
        HSW = ( qSW(4) + pSW ) / rSW;

        % East state
        rSE = qSE(1);
        uSE = qSE(2)/rSE;
        vSE = qSE(3)/rSE;
        pSE = (Param.GC.gamma-1)*( qSE(4) - rSE*(uSE^2+vSE^2)/2 );
        aSE = sqrt(Param.GC.gamma*pSE/rSE);
        HSE = ( qSE(4) + pSE ) / rSE;

        % South state
        rNW = qNW(1);
        uNW = qNW(2)/rNW;
        vNW = qNW(3)/rNW;
        pNW = (Param.GC.gamma-1)*( qNW(4) - rNW*(uNW^2+vNW^2)/2 );
        aNW = sqrt(Param.GC.gamma*pNW/rNW);
        HNW = ( qNW(4) + pNW ) / rNW;

        % North state
        rNE = qNE(1);
        uNE = qNE(2)/rNE;
        vNE = qNE(3)/rNE;
        pNE = (Param.GC.gamma-1)*( qNE(4) - rNE*(uNE^2+vNE^2)/2 );
        aNE = sqrt(Param.GC.gamma*pNE/rNE);
        HNE = ( qNE(4) + pNE ) / rNE;




        % Compute Roe Averages - SW to SE
        rSroe = sqrt(rSE/rSW);
        uSroe = (uSW+rSroe*uSE)/(1+rSroe);
        vSroe = (vSW+rSroe*vSE)/(1+rSroe);
        HSroe = (HSW+rSroe*HSE)/(1+rSroe);
        aSroe = sqrt( (Param.GC.gamma-1)*(HSroe-0.5*(uSroe^2+vSroe^2)) );

        % Compute Roe Averages - NW to NE
        rNroe = sqrt(rNE/rNW);
        uNroe = (uNW+rNroe*uNE)/(1+rNroe);
        vNroe = (vNW+rNroe*vNE)/(1+rNroe);
        HNroe = (HNW+rNroe*HNE)/(1+rNroe);
        aNroe = sqrt( (Param.GC.gamma-1)*(HNroe-0.5*(uNroe^2+vNroe^2)) );

        % Compute Roe Averages - SW to NW
        rWroe = sqrt(rNW/rSW);
        uWroe = (uSW+rWroe*uNW)/(1+rWroe);
        vWroe = (vSW+rWroe*vNW)/(1+rWroe);
        HWroe = (HSW+rWroe*HNW)/(1+rWroe);
        aWroe = sqrt( (Param.GC.gamma-1)*(HWroe-0.5*(uWroe^2+vWroe^2)) );

        % Compute Roe Averages - SE to NE
        rEroe = sqrt(rNE/rSE);
        uEroe = (uSE+rEroe*uNE)/(1+rEroe);
        vEroe = (vSE+rEroe*vNE)/(1+rEroe);
        HEroe = (HSE+rEroe*HNE)/(1+rEroe);
        aEroe = sqrt( (Param.GC.gamma-1)*(HEroe-0.5*(uEroe^2+vEroe^2)) );




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
