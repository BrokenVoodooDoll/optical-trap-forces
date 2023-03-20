function factor = qg_y_factor(beta, r, y, n1, n2, Rsp, focus)
    factor = qg_avg_factor(th_i_y(beta, r, y, Rsp, focus), n1, n2).*sin(gamma_angle(beta, r, focus));
end