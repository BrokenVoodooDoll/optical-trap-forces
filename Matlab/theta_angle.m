function angle = theta_angle(beta, r, y, Rsp, focus)
    angle = asin(y / Rsp .* sin(gamma(beta, r, focus)));
end