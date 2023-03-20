function factor = qs_y_factor(beta, r, y, n1, n2, Rsp, focus)
   factor = qs_avg_factor(th_i_y(beta, r, y, Rsp, focus), n1, n2).*cos(phi_i(r, focus));
end