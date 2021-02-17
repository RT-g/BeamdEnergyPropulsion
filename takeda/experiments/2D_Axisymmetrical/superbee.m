function g = superbee(x,y)
    s = sign(y);
    g = s * max([0 min([2*abs(y) s*x]) min([abs(y) 2*s*x])]);
end