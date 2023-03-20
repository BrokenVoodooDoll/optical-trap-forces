function factor = qg_z_factor(r, z, n1, n2, Rsp, focus)
    factor = -qg_avg_factor(th_i_z(r, z, Rsp, focus), n1, n2) .* sin(phi_i(r, focus));
end