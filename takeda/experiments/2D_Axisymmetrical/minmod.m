function g = minmod(x,y)
    % 制限関数を計算する
    g = sign(x) * max([0 min([abs(x) sign(x)*y])]);
end