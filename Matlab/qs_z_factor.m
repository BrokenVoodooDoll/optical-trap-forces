function factor = qs_z_factor(r, z, n1, n2, Rsp, focus)
    factor = qs_avg_factor(th_i_z(r, z, Rsp, focus), n1, n2) .* cos(phi_i(r, focus));
end