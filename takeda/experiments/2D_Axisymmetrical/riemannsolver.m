function solver = riemannsolver(gamma,rho,u,v,H,P,U1,U2,U3,U4,matrix)
    D = sqrt(rho(2)/rho(1));
    if ~isreal(D) % たまにrhoがマイナスになってしまうのでその時は大きいほうをとる。
        D = max([rho(1) rho(2)]);
    end
    u_av = (u(1)+D*u(2))/(1+D); 
    v_av = (v(1)+D*v(2))/(1+D); 
    H_av = (H(1)+D*H(2))/(1+D);
    % 非線形なので両隣よりもすごく大きい値が出るかもしれないので、計算開始直後はできるだけ変な値が出ないように制限
    a_av = sqrt((gamma-1)*(H_av-1/2*(u_av.^2 + v_av.^2)));
    % if ~isreal(a_av)
    %     a_av = sqrt(max([gamma*P(1)/rho(1) gamma*P(2)/rho(2)])); % Square root of the derivative of P br rho, used in the eigenvalues of the Jacobian matrix
    % end
    if ~isreal(a_av)
        disp(a_av)
        disp(rho)
        disp(u)
        disp(v)
        disp(H)
        disp(u_av)
        disp(v_av)
        disp(H_av)
        disp((gamma-1)*(H_av-1/2*(u_av.^2 + v_av.^2)))
        disp(P)
        disp(rho)
        disp(gamma)
        disp(gamma*P(1)/rho(1))
        disp(gamma*P(2)/rho(2))
        disp(matrix)
        disp('a_avが崩壊')
        return 
    end

    % R行列の生成 from eigenvectors in w-space
    R = zeros(4,4);
    switch matrix
        case 'F'
            % ヤコビアン行列の固有値計算 
            lambda(1) = u_av;
            lambda(2) = u_av;
            lambda(3) = u_av+a_av;
            lambda(4) = u_av-a_av;
            R(1,:) = 1;
            R(1,2) = 0;
            R(2,:) = lambda;
            R(2,2) = 0;
            R(3,:) = v_av;
            R(3,2) = 1;
            R(4,1) = 1/2*(u_av.^2 + v_av.^2);
            R(4,2) = v_av;
            R(4,3) = H_av+u_av.*a_av;
            R(4,4) = H_av-u_av.*a_av;
        case 'G'
            % ヤコビアン行列の固有値計算 
            lambda(1) = v_av;
            lambda(2) = v_av;
            lambda(3) = v_av+a_av;
            lambda(4) = v_av-a_av;
            R(1,:) = 1;
            R(1,2) = 0;
            R(2,:) = u_av;
            R(2,2) = -1;
            R(3,:) = lambda;
            R(3,2) = 0;
            R(4,1) = 1/2*(u_av.^2 + v_av.^2);
            R(4,2) = -u_av;
            R(4,3) = H_av+v_av.*a_av;
            R(4,4) = H_av-v_av.*a_av;
    end
    % A =  diag(abs(lambda));
    % R_av = R*A/R; % Roe平均
    
    dU = zeros(4,1);
    dU(1) = U1(2)-U1(1);
    dU(2) = U2(2)-U2(1);
    dU(3) = U3(2)-U3(1);
    dU(4) = U4(2)-U4(1);

    % Solving the linear problem
    % 2次精度風上TVD
    alpha = R\dU; % 特性曲線を横切った時の変化量
    solver = {alpha, R, lambda};
end