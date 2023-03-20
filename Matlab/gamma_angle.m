function angle = gamma_angle(beta, r, focus)
    angle = acos(cos(pi/2 - phi_i(r, focus)) .* cos(beta));
end