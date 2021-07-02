function [r_0,u_0,v_0,p_0] = setEulerInitialCondition2d(x,y,IC,gas)
    % Load the IC of a 1D Riemann classical schok tube problem configuration.
    % In the notation we take advantage of the matlab array notation as follows
    %
    %   1.0 +-----------+-----------+
    %       |           |           |
    %       |   reg 2   |   reg 1   |
    %       |           |           |
    %   0.5 +-----------+-----------+
    %       |           |           |
    %       |   reg 3   |   reg 4   |
    %       |           |           |
    %   0.0 +-----------+-----------+
    %      0.0         0.5         1.0
    %
    % prop = [prop_reg1 , prop_reg2 , prop_reg3 , prop_reg4]
    %
    %   r = rho/density
    %   u = velocity in x direction
    %   v = velocity in y direction
    %   p = Pressure
    %
    % Manuel Diaz, NTU, 2014.06.27

    %% Initial Physical Properties per case:
    % Values specific to gas
    global Param
    switch gas
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

    switch IC
        case{1} % Configuration 1
            fprintf('Configuration 1 \n');
            p = [1.0  1.0  1.0  1.0];
            r = [1.0  1.0  1.0  1.0];
            u = [0.0  0.0  0.0  0.0];
            v = [0.0  0.0  0.0  0.0];

        case 'Sod_x'
            fprintf('Sods Shocktube in the x-direction (2-d test) \n');
            p = [0.1   1 1 0.1  ];
            r = [0.125 1 1 0.125];
            u = [0     0 0 0    ];
            v = [0     0 0 0    ];

        case 'constant'
            fprintf('constant state (2-d test) \n');
            p = [1.0  1.0  1.0  1.0];
            r = [Param.GC.rho_0  Param.GC.rho_0  Param.GC.rho_0  Param.GC.rho_0];
            u = [0.0  0.0  0.0  0.0];
            v = [0.0  0.0  0.0  0.0];

        otherwise
            error('only 18 cases are available');
    end
    %% Print configuration of selected IC
    fprintf('\n');
    fprintf('          reg 1 reg 2  reg 3  reg 4\n');
    fprintf('density : %2.4f %2.4f %2.4f %2.4f \n',r);
    fprintf('  x-vel : %2.4f %2.4f %2.4f %2.4f \n',u);
    fprintf('  y-vel : %2.4f %2.4f %2.4f %2.4f \n',v);
    fprintf('Presure : %2.4f %2.4f %2.4f %2.4f \n',p);
    fprintf('\n');

    %% Load Selected case Initial condition:

    % Parameters of regions dimensions
    reg1 = (x>=0.5 & y>=0.5); % region 1
    reg2 = (x <0.5 & y>=0.5); % region 2
    reg3 = (x <0.5 & y <0.5); % region 3
    reg4 = (x>=0.5 & y <0.5); % region 4

    % Initial Condition for our 2D domain
    r_0 = r(1)*reg1 + r(2)*reg2 + r(3)*reg3 + r(4)*reg4; % Density, rho
    u_0 = u(1)*reg1 + u(2)*reg2 + u(3)*reg3 + u(4)*reg4; % velocity in x
    v_0 = v(1)*reg1 + v(2)*reg2 + v(3)*reg3 + v(4)*reg4; % velocity in y
    p_0 = p(1)*reg1 + p(2)*reg2 + p(3)*reg3 + p(4)*reg4; % temperature.

    end
