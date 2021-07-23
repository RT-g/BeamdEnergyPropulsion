function [nx,nr,in,jn] = setGhostCells(nx,nr,O)
    global Param
    % plate
    switch Param.BC.left
    case 'closed'
        in_s = 1;
        in_e = nx;
    case 'open'
        in_s = 1 + O;
        in_e = nx + O;
        nx = nx + O;
    end

    % axis
    switch Param.BC.bottom
    case 'rotational'
        jn_s = 1 + O;
        jn_e = nr + O;
        nr = nr + O;
    end

    % outflow
    switch Param.BC.right
    case 'open'
        nx = nx + O;
    otherwise
        error('Boundary Conditiion right is only open');
    end

    % outflow or tube
    switch Param.BC.top
    case 'open'
        nr = nr + O;
    end

    in = in_s:in_e;
    jn = jn_s:jn_e;
end
