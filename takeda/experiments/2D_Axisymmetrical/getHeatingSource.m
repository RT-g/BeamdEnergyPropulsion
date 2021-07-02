function w = getHeatingSource(x, cosine)
    global Param
    % 加熱領域をx_laser0よりも前にしてしまうと計算が壊れるので実際に考えられる値よりも少し進めたほうがいい?
    if x > Param.LV.x_laser0-Param.LV.l && x < Param.LV.x_laser0 % ξが現実の長さを表しているわけではないので実際は少し補正が必要だが、ここでは無視
        w = Param.rp.eta * Param.LV.S_laser0 * cosine / Param.LV.l * 1e9; %W/m3
    else
        w = 0;
    end
end
