function setLaserConstants
    global Param
    Param.LV.A_G = 0.224517656;
    Param.LV.B_G = 0.77548;
    Param.LV.sigma_G1 = 0.84473; %mm
    Param.LV.sigma_G2 = 1.75302; %mm
    Param.LV.A_T = 0.764525962;
    Param.LV.B_T = 0.235474467;
    Param.LV.sigma_T1 = 1.557528296;
    Param.LV.sigma_T2 = 4.050355336;
    Param.LV.laser_lambda = 10.6; %単位はum
    Param.LV.M2_G = 15;
    Param.LV.M2_T = 21;
    Param.LV.W_G0 = 1.7;%単位はmm
    Param.LV.W_T0 = 2.0;
    Param.LV.W_G = Param.LV.W_G0 * 1e-3;
    Param.LV.W_T = Param.LV.W_T0 * 1e-3;
    Param.LV.R_peak = 2.13;
    Param.LV.a_s = (Param.LV.B_G*1e6/Param.LV.sigma_G2^2 + Param.LV.B_T*1e6/Param.LV.sigma_T2^2);
    Param.LV.l = 0.2e-3; %m heating length 加熱長さ レーザーは0.2 mm
end
