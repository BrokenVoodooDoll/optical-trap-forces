function peak = gauss_peak(r_max, w0)
    A = (1 - exp(-2*r_max.^2 / w0^2));
    peak = 2*A / (pi * w0^2);
end