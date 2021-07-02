function [u_ionz, intercept, B, r, cosine] = getSWVelocity(y)
    % 時刻に合わせて電離波面伝播速度を計算する。電離波面速度が早い方を採用するという仮定
    % TODO もう少し厳密な定義を考えたい
    global Param
    [u_ionz_line1, intercept_line1, b_s_line1, r_line1, cosine_line1] = SWvelocity(Param.GC.slope, Param.GC.intercept, y);
    [u_ionz_line3, intercept_line3, b_s_line3, r_line3, cosine_line3] = SWvelocity(Param.GC.slope_low, Param.GC.intercept_low, y);

    if (u_ionz_line1 > u_ionz_line3)
        u_ionz    = u_ionz_line1;
        intercept = intercept_line1;
        B         = b_s_line1;
        r         = r_line1;
        cosine    = cosine_line1;
    else
        u_ionz    = u_ionz_line3;
        intercept = intercept_line3;
        B         = b_s_line3;
        r         = r_line3;
        cosine    = cosine_line3;
    end
end

function [u_ionz, intercept, B, r, cosine] = SWvelocity(slope, intercept, y)
    global Param
    A = Param.LV.a_s;
    B = intercept/(1-intercept);
    S_laser0 = Param.LV.S_laser0;
    a = Param.LV.S_laser0;
    rho_0 = Param.GC.rho_0;

    switch Param.hp.coodinate
        case 'general'
            r = -(2/(3*A^2*B^2*y + sqrt(9*A^4*B^4*y^2 + 4*A^3*B^3)))^(1/3) + ((3*A^2*B^2*y + sqrt(9*A^4*B^4*y^2 + 4*A^3*B^3))/2)^(1/3)/A/B; %m, Cartesian coordinates
        case 'cartesian'
            r = y;
    end
    G_r = Param.LV.A_G*exp(-2*(r*1e3/sqrt(2)/Param.LV.sigma_G1)^4)+Param.LV.B_G*exp(-2*(r*1e3/sqrt(2)/Param.LV.sigma_G2)^2);
    T_r = Param.LV.A_T*exp(-2*(r*1e3/sqrt(2)/Param.LV.sigma_T1)^4)+Param.LV.B_T*exp(-2*(r*1e3/sqrt(2)/Param.LV.sigma_T2)^2);
    cosine = (G_r*T_r)^B;
    u_ionz = slope * (S_laser0 * cosine * 1e9 / rho_0 / a^3)^intercept * a; %m/s 波頭の速度
end