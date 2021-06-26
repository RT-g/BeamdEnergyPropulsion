function [u_ionz0] = DecideSWVelocity(S_laser0, slope, slope_low, intercept, intercept_low, a, rho_0)
    % 時刻に合わせて電離波面伝播速度を計算する。電離波面速度が早い方を採用するという仮定
    % TODO もう少し厳密な定義を考えたい
    u_ionz_line1 = slope * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept * a; %m/s 波頭の速度
    u_ionz_line3 = slope_low * (S_laser0 * 1e9 / rho_0 / a ^ 3) ^ intercept_low * a;
    if (u_ionz_line1 > u_ionz_line3)
        u_ionz0 = u_ionz_line1; % Line1, 松井さん博論
    else
        u_ionz0 = u_ionz_line3; % Line3, 松井さん博論
    end
end
