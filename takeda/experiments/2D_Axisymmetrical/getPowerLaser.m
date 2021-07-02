function Power_laser = getPowerLaser(t, precise)
    % 時刻に合わせて電離波面伝播速度を計算する。電離波面速度が早い方を採用するという仮定
    arguments % デフォルト値を規定
        t = 0;
        precise = false;
    end

    if precise % 正確に再現すると温度などが高く出すぎる場合がある
        if (t*1e6<0.085)
            Power_laser = 287.9222823*t*1e6+0.0005175756469;
        elseif (t*1e6<0.125)
            Power_laser =-476.3277981*t*1e6+66.1043297;
        else
            Power_laser = 8.15*exp(-0.866*t*1e6);
        end
    else
        Power_laser = 8.15*exp(-0.866*t*1e6);
    end
end
