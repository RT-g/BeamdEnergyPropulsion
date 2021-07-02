function [cell] = setBoundaryCondition(cell, nx, nr, closed)
    % Left(回転軸)
    cell(:,1).res(1) = cell(:,3).res(1);
    cell(:,2).res(1) = cell(:,3).res(1);
    cell(:,1).res(2) = cell(:,3).res(2);
    cell(:,2).res(2) = cell(:,3).res(2);
    cell(:,1).res(3) = 0; %中心ではur=0のはず
    cell(:,2).res(3) = cell(:,3).res(3);
    cell(:,1).res(4) = cell(:,3).res(4);
    cell(:,2).res(4) = cell(:,3).res(4);

    % Right
    if closed
        % 閉じ込めあり(tube)の場合
        cell(:,nr).res(1)   = cell(:,nr-2).res(1);
        cell(:,nr-1).res(1) = cell(:,nr-2).res(1);
        cell(:,nr).res(2)   = 0;
        cell(:,nr-1).res(2) = cell(:,nr-2).res(2); %0?
        cell(:,nr).res(3)   = 0;
        cell(:,nr-1).res(3) = cell(:,nr-2).res(3); %0?
        cell(:,nr).res(4)   = cell(:,nr-2).res(4);
        cell(:,nr-1).res(4) = cell(:,nr-2).res(4);
    else
        % 閉じ込めなしの場合
        cell(:,nr).res = cell(:,nr-2).res;
        cell(:,nr-1).res = cell(:,nr-2).res;
    end

    % Top(光軸方向)流出
    cell(nx,:).res = cell(nx-2,:).res;
    cell(nx-1,:).res = cell(nx-2,:).res;

    % Plate(アルミ板)
    cell(1,:).res(1) = cell(3,:).res(1);
    cell(2,:).res(1) = cell(3,:).res(1);
    cell(1,:).res(2) = 0;
    cell(2,:).res(2) = cell(3,:).res(2);%0?
    cell(1,:).res(3) = 0;
    cell(2,:).res(3) = cell(3,:).res(3);%0?
    cell(1,:).res(4) = cell(3,:).res(4);
    cell(2,:).res(4) = cell(3,:).res(4);

%     % Corners
%     Q1_cal(nx,1) = (Q1_cal(nx,2) + Q1_cal(nx-1,1))/2;
%     Q2_cal(nx,1) = (Q2_cal(nx,2) + Q2_cal(nx-1,1))/2;
%     Q3_cal(nx,1) = (Q3_cal(nx,2) + Q3_cal(nx-1,1))/2;
%     Q4_cal(nx,1) = (Q4_cal(nx,2) + Q4_cal(nx-1,1))/2;
% 
%     Q1_cal(nx,nr) = (Q1_cal(nx-1,nr) + Q1_cal(nx,nr-1))/2;
%     Q2_cal(nx,nr) = (Q2_cal(nx-1,nr) + Q2_cal(nx,nr-1))/2;
%     Q3_cal(nx,nr) = (Q3_cal(nx-1,nr) + Q3_cal(nx,nr-1))/2;
%     Q4_cal(nx,nr) = (Q4_cal(nx-1,nr) + Q4_cal(nx,nr-1))/2;
end
