function [u_ionz, intercept, b_s, r, cosine] = DecideSWVelocity(S_laser0, slope, slope_low, intercept, intercept_low, a, rho_0, a_s, h, sigma_G1, sigma_G2, sigma_T1, sigma_T2, A_G, A_T, B_G, B_T)
    % 時刻に合わせて電離波面伝播速度を計算する。電離波面速度が早い方を採用するという仮定
    % TODO もう少し厳密な定義を考えたい
    [u_ionz_line1, intercept_line1, b_s_line1, r_line1, cosine_line1] = SWvelocity(S_laser0, slope, intercept, rho_0, a, a_s, h, sigma_G1, sigma_G2, sigma_T1, sigma_T2, A_G, A_T, B_G, B_T)
    [u_ionz_line3, intercept_line3, b_s_line3, r_line3, cosine_line3] = SWvelocity(S_laser0, slope_low, intercept_low, rho_0, a, a_s, h, sigma_G1, sigma_G2, sigma_T1, sigma_T2, A_G, A_T, B_G, B_T)
    if (u_ionz_line1 > u_ionz_line3)
        [u_ionz, intercept, b_s, r, cosine] = [u_ionz_line1, intercept_line1, b_s_line1, r_line1, cosine_line1]
    else
        [u_ionz, intercept, b_s, r, cosine] = [u_ionz_line3, intercept_line3, b_s_line3, r_line3, cosine_line3]
    end
end

function [u_ionz, intercept, b_s, r, cosine] = SWvelocity(S_laser0, slope, intercept, rho_0, a, a_s, h, sigma_G1, sigma_G2, sigma_T1, sigma_T2, A_G, A_T, B_G, B_T)
    b_s = intercept/(1-intercept);
    r = -(2/(3*a_s^2*b_s^2*h + sqrt(9*a_s^4*b_s^4*h^2 + 4*a_s^3*b_s^3)))^(1/3) + ((3*a_s^2*b_s^2*h + sqrt(9*a_s^4*b_s^4*h^2 + 4*a_s^3*b_s^3))/2)^(1/3)/a_s/b_s; %m, Cartesian coordinates
    G_r = A_G*exp(-2*(r*1e3/sqrt(2)/sigma_G1)^4)+B_G*exp(-2*(r*1e3/sqrt(2)/sigma_G2)^2);
    T_r = A_T*exp(-2*(r*1e3/sqrt(2)/sigma_T1)^4)+B_T*exp(-2*(r*1e3/sqrt(2)/sigma_T2)^2);
    cosine = (G_r*T_r)^b_s;
    u_ionz = slope * (S_laser0 * cosine * 1e9 / rho_0 / a ^ 3) ^ intercept * a; %m/s 波頭の速度
end
