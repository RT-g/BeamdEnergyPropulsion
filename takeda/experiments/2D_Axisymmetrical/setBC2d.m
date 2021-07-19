function  q = setBC2d(q,nx,nr,O)
    % axis
    switch Param.BC.bottom
        case 'rotational'
            q(:,1,1) = q(:,1+O,1);
            q(:,1,2) = q(:,1+O,2);
            q(:,1,3) = 0; %中心ではur=0のはず
            q(:,1,4) = q(:,1+O,4);
            if O = 2
                q(:,2,1) = q(:,1+O,1);
                q(:,2,2) = q(:,1+O,2);
                q(:,2,3) = 0; %中心ではur=0のはず
                q(:,2,4) = q(:,1+O,4);
            end
    end

    % outflow or tube
    switch Param.BC.top
        case 'open'
            q(:,nr,1) = q(:,nr-O,1);
            q(:,nr,2) = q(:,nr-O,2);
            q(:,nr,3) = q(:,nr-O,3);
            q(:,nr,4) = q(:,nr-O,4);
            if O = 2
                q(:,nr-1,1) = q(:,nr-O,1);
                q(:,nr-1,2) = q(:,nr-O,2);
                q(:,nr-1,3) = q(:,nr-O,3);
                q(:,nr-1,4) = q(:,nr-O,4);
            end
        case 'closed'
            q(:,nr,1) = q(:,nr-O,1);
            q(:,nr,2) = 0;
            q(:,nr,3) = 0;
            q(:,nr,4) = q(:,nr-O,4);
            if O = 2
                q(:,nr-1,1) = q(:,nr-O,1);
                q(:,nr-1,2) = 0;
                q(:,nr-1,3) = 0;
                q(:,nr-1,4) = q(:,nr-O,4);
            end
    end

    % plate
    switch Param.BC.left
        case 'open'
            q(1,:,1) = q(1+O,:,1);
            q(1,:,2) = q(1+O,:,2);
            q(1,:,3) = q(1+O,:,3);
            q(1,:,4) = q(1+O,:,4);
            if O = 2
                q(2,:,1) = q(1+O,:,1);
                q(2,:,2) = q(1+O,:,2);
                q(2,:,3) = q(1+O,:,3);
                q(2,:,4) = q(1+O,:,4);
            end
        case 'closed'
            q(1,:,1) = q(1+O,:,1);
            q(1,:,2) = 0;
            q(1,:,3) = 0;
            q(1,:,4) = q(1+O,:,4);
            if O = 2
                q(2,:,1) = q(1+O,:,1);
                q(2,:,2) = 0;
                q(2,:,3) = 0;
                q(2,:,4) = q(1+O,:,4);
            end
    end

    % outflow
    switch Param.BC.right
        case 'open'
            q(nx,:,1) = q(nx-O,:,1);
            q(nx,:,2) = q(nx-O,:,2);
            q(nx,:,3) = q(nx-O,:,3);
            q(nx,:,4) = q(nx-O,:,4);
            if O = 2
                q(nx-1,:,1) = q(nx-O,:,1);
                q(nx-1,:,2) = q(nx-O,:,2);
                q(nx-1,:,3) = q(nx-O,:,3);
                q(nx-1,:,4) = q(nx-O,:,4);
            end
        otherwise
            error('Boundary Conditiion right is only open');
    end

end
