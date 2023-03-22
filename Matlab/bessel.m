% Bessel beam
function b = bessel(r, w0, r_max)
    ring_radius = 2.405; % radisu of the first ring of the besselj_0
    A = 2 * ring_radius / (w0 * r_max * exp(0.5)) * 2 * pi * integral(@(r) r.*besselj(0, ring_radius/w0 * r).^2, 0, r_max);
    b = A * besselj(0, 2.405/w0 * r).^2;
end