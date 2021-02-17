function pusai = enthalpy_fix(z,delta)
    % エンタルピー補正量を計算する
    if abs(z) > delta
        pusai = abs(z);
    else
        pusai = (z^2+delta^2)/2/delta;
    end
end