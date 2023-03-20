% Gaussian TEM00 beam
function g = gauss(r, r_max, w0)
    g = exp(-2 * r.^2 / w0^2);
end