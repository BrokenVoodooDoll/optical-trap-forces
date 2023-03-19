load_constants

function angle = th_r(theta)
    angle = asin(n1 / n2 * sin(theta));
end

function angle = ph_i(r)
    angle = atan(r / f);
end

function angle = th_i(r, z)
    angle = asin(z / Rsp .* sin(ph_i(r)));
end

function angle = gamma(beta, r)
    angle = acos(cos(pi/2 - ph_i(r)) .* cos(beta));
end

function angle = theta(beta, r, y)
    angle = asin(y / Rsp .* sin(gamma(beta, r)));
end