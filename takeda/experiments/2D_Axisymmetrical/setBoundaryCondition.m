function [Q1_cal, Q2_cal, Q3_cal, Q4_cal] = setBoundaryCondition(Q1_cal, Q2_cal, Q3_cal, Q4_cal, nx, nh, closed)
    arguments
        closed = false;
    end

    % Left(回転軸)
    for i = 2:nx-1
        Q1_cal(i,1) = Q1_cal(i,2);
        Q2_cal(i,1) = Q2_cal(i,2);
        Q3_cal(i,1) = 0; %中心ではur=0のはず
        Q4_cal(i,1) = Q4_cal(i,2);
    end

    % Right
    if closed
        % 閉じ込めあり(tube)の場合
        for i = 2:nx-1
            Q1_cal(i,nh) = Q1_cal(i,nh-1);
            Q2_cal(i,nh) = 0;
            Q3_cal(i,nh) = 0;
            Q4_cal(i,nh) = Q4_cal(i,nh-1);
        end
    else
        % 閉じ込めなしの場合
        for i = 2:nx-1
            Q1_cal(i,nh) = Q1_cal(i,nh-1);
            Q2_cal(i,nh) = Q2_cal(i,nh-1);
            Q3_cal(i,nh) = Q3_cal(i,nh-1);
            Q4_cal(i,nh) = Q4_cal(i,nh-1);
        end
    end

    % Top(光軸方向)流出
    for i = 2:nh-1
        Q1_cal(nx,i) = Q1_cal(nx-1,i);
        Q2_cal(nx,i) = Q2_cal(nx-1,i);
        Q3_cal(nx,i) = Q3_cal(nx-1,i);
        Q4_cal(nx,i) = Q4_cal(nx-1,i);
    end

    % Plate(アルミ板)
    Q1_cal(1,:) = Q1_cal(2,:);
    Q2_cal(1,:) = 0;
    Q3_cal(1,:) = 0;
    Q4_cal(1,:) = Q4_cal(2,:);

    % Corners
    Q1_cal(nx,1) = (Q1_cal(nx,2) + Q1_cal(nx-1,1))/2;
    Q2_cal(nx,1) = (Q2_cal(nx,2) + Q2_cal(nx-1,1))/2;
    Q3_cal(nx,1) = (Q3_cal(nx,2) + Q3_cal(nx-1,1))/2;
    Q4_cal(nx,1) = (Q4_cal(nx,2) + Q4_cal(nx-1,1))/2;

    Q1_cal(nx,nh) = (Q1_cal(nx-1,nh) + Q1_cal(nx,nh-1))/2;
    Q2_cal(nx,nh) = (Q2_cal(nx-1,nh) + Q2_cal(nx,nh-1))/2;
    Q3_cal(nx,nh) = (Q3_cal(nx-1,nh) + Q3_cal(nx,nh-1))/2;
    Q4_cal(nx,nh) = (Q4_cal(nx-1,nh) + Q4_cal(nx,nh-1))/2;
end
