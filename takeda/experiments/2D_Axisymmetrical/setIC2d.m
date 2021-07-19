function [Q0,a0] = setEulerInitialCondition2d(x,y,IC)
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

    switch IC
        case{1} % Configuration 1
            fprintf('Configuration 1 \n');
            p = [1.0  1.0  1.0  1.0]*1e5;
            r = [1.0  1.0  1.0  1.0];
            u = [0.0  0.0  0.0  0.0];
            v = [0.0  0.0  0.0  0.0];

        case 'Sod_x'
            fprintf('Sods Shocktube in the x-direction (2-d test) \n');
            p = [0.1   1 1 0.1  ]*1e5;
            r = [0.125 1 1 0.125];
            u = [0     0 0 0    ];
            v = [0     0 0 0    ];

        case 'constant'
            fprintf('constant state\n');
            p = [1.0  1.0  1.0  1.0]*1e5;
            r = [1.0  1.0  1.0  1.0]. * Param.GC.rho_0;
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
    r0 = r(1)*reg1 + r(2)*reg2 + r(3)*reg3 + r(4)*reg4; % Density, rho
    u0 = u(1)*reg1 + u(2)*reg2 + u(3)*reg3 + u(4)*reg4; % velocity in x
    v0 = v(1)*reg1 + v(2)*reg2 + v(3)*reg3 + v(4)*reg4; % velocity in y
    p0 = p(1)*reg1 + p(2)*reg2 + p(3)*reg3 + p(4)*reg4; % temperature.

    E0 = p0./((Param.GC.gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
    c0 = sqrt(Param.GC.gamma*p0./r0);                    % Speed of sound
    Q0 = cat(3, r0, r0.*u0, r0.*v0, r0.*E0);    % initial state
    vn = sqrt(u0.^2+v0.^2); lambda1=vn+c0; lambda2=vn-c0;
    a0 = max(abs([lambda1(:);lambda2(:)]));
    end
