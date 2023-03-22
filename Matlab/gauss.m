% Gaussian TEM00 beam
function g = gauss(r, w0, r_max)
    A = (1 - exp(-2*r_max.^2 / w0^2)); % the fraction of power that falls on the pupil of the micro lens
    g = 2 * A * exp(-2 * r.^2 / w0^2);
end