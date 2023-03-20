function angle = th_i_y(beta, r, y, Rsp, focus)
    angle = asin(y / Rsp .* sin(gamma_angle(beta, r, focus)));
end