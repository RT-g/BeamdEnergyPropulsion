function s = tvd_sigma(z,lambda,delta)
    s = 1/2 * (enthalpy_fix(z,delta)-lambda*z^2);
end
