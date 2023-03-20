function peak = bessel_peak(r_max, w0, P)
    peak = P * 4.81 / (w0 * r_max * exp(0.5)) * 2 * pi * integral(@(r) r.*besselj(0, 2.405/w0 * r).^2, 0, r_max) / P0;
end