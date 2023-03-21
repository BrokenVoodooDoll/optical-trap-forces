% Bessel beam
function b = bessel(r, r_max, w0, P)
    b = besselj(0, 2.405/w0 * r).^2;
end