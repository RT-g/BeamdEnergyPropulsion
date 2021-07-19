function Settings
    global Param
    % Computation Value
    Param.hp.CFL     = 0.50;     % CFL number. 1以下にする 0.5~0.99くらい;
    Param.hp.tEnd    = 0.05;     % Final time;
    Param.hp.nx      = 50;      % Number of cells/Elements in x;
    Param.hp.nr      = 50;      % Number of cells/Elements in r;
    Param.hp.dim     = '2dAxisymmetric';       % dimention
    Param.hp.IC      = 'constant';       % 19 IC cases are available;
    Param.hp.gas     = 'air';    % selected gas, air/argon/helium;
    Param.hp.method  = 4;	    % 1:Dim by Dim, 2:HLLE2d 1st-order, 3:HLLE2d 2nd-order, 4 2D Axisymmetrical, 5;
    Param.hp.plotFig = true;     % true:visualize evolution;
    Param.hp.Jacobian = @Jacobian; %coordinates, general/cartesian
    Param.hp.fluxMth ='ChakravarthyOsher';  % HLLE1d, HLLE2d;
    Param.hp.tineMth ='ADI';  % HLLE1d, HLLE2d;
    Param.hp.limiter = 'mm';
    Param.hp.theta_1 = 1;
    Param.hp.theta_2 = 0.5;
    Param.hp.kappa = 1/3;
    Param.hp.bcomp = 2; %圧縮パラメータ, (1,4]
    Param.hp.eta = 0.1; %圧縮パラメータ, (1,4]
    switch Param.hp.dim
    case '2dAxysymmetric'
        Param.hp.dnum = 2;
    case '2d'
        Param.hp.dnum = 2;
    end

    % Boundary condition
    Param.BC.left = 'closed'
    Param.BC.bottom = 'rotational'
    Param.BC.right = 'open'
    Param.BC.top = 'open'

    % Laser Constants
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

    % gas constants
    switch Param.hp.gas
        case 'air'
            fprintf('gas: air \n');
            Param.GC.m = 28;
            Param.GC.gamma = 1.4; %% Ratio of specific heats for ideal di-atomic gas.
            Param.GC.eta_trans = 0.48;
            Param.GC.slope = 0.16;
            Param.GC.intercept = 0.46;
            Param.GC.slope_low = 0.68;
            Param.GC.intercept_low = 0.3;
            Param.GC.a = 347;

        case 'argon'
            fprintf('gas: air \n');
            Param.GC.m = 40;
            Param.GC.gamma = 1.67; %% Ratio of specific heats for ideal di-atomic gas.
            Param.GC.eta_trans = 0.42;
            Param.GC.slope = 1.18;
            Param.GC.intercept = 0.21;
            % slope_low = 0.68;
            % intercept_low = 0.3;
            % a = 347

        case 'helium'
            fprintf('gas: helium \n');
            Param.GC.m = 4;
            Param.GC.gamma = 1.67; %% Ratio of specific heats for ideal di-atomic gas.
            Param.GC.eta_trans = 0.54;
            Param.GC.slope = 0.004519362;
            Param.GC.intercept = 1.10603162;
            % slope_low = 0.68;
            % intercept_low = 0.3;
            % a = 347
        otherwise
            error('only 3 gas species are available')
    end
    R_0 = 8314.0;
    R = R_0/Param.GC.m;
    Param.GC.R = R;
    P_0 = 1.013e5*1;
    T_0 = 293;
    Param.GC.rho_0 = P_0/R/T_0;

end

function j = Jacobian(vec)
% 自分で必要に応じて書き換えること
    arguments
        vec = none;
    end
    j = ones(Param.hp.dnum);
end
