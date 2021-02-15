function solver = riamansolver(gamma,rho,u,v,H,P,U1,U2,U3,U4,matrix)
    % rho
    % u
    % v
    % H
    % P
    % U1
    % U2
    % U3
    % U4    
    rho_av = sqrt(rho(1)*rho(2));
    u_av = (sqrt(rho(1))*u(1)+sqrt(rho(2))*u(2))/(sqrt(rho(1))+sqrt(rho(2))); 
    v_av = (sqrt(rho(1))*v(1)+sqrt(rho(2))*v(2))/(sqrt(rho(1))+sqrt(rho(2))); 
    H_av = (sqrt(rho(1))*H(1)+sqrt(rho(2))*H(2))/(sqrt(rho(1))+sqrt(rho(2)));
    % 非線形なので両隣よりもすごく大きい値が出るかもしれないので、計算開始直後はできるだけ変な値が出ないように制限
    a_av = sqrt(max((gamma-1)*(H_av-1/2*(u_av.^2 + v_av.^2)), min(gamma*P(1)/rho(1),gamma*P(2)/rho(2)))); % Square root of the derivative of P br rho, used in the eigenvalues of the Jacobian matrix

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
    dU = [U1(2)-U1(1) U2(2)-U2(1) U3(2)-U3(1) U4(2)-U4(1)]';
    % Solving the linear problem
    % 2次精度風上TVD
    alpha = inv(R)*dU;
    solver = {alpha, R, lambda};
end