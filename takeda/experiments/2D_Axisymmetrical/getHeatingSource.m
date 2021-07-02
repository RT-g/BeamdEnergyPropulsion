function W = getHeatingSource(W, nx, dx x_laser0, l, S_laser0)

    for n2 = 1:nx
        x = n2*dx; %m, general coordinate
        % 加熱領域をx_laser0よりも前にしてしまうと計算が壊れるので実際に考えられる値よりも少し進めたほうがいい?
        if x > x_laser0-l && x < x_laser0 % ξが現実の長さを表しているわけではないので実際は少し補正が必要だが、ここでは無視
            % W(n2,n3) = eta * S_laser(n3) / l * 1e9; %W/m3
            W(n2,n3) = eta * S_laser0 / l * 1e9; %W/m3
        else
            W(n2,n3) = 0;å
        end
    end
end
