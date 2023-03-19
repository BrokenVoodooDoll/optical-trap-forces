load_constants

function g = gauss(r)
    A = (1-exp(-2*r_max.^2 / w0^2));
    I0 = 2*P / (pi * w0^2 * A);
    g = I0 * exp(-2 * r.^2 / w0^2); % Gaussian TEM00 beam
end