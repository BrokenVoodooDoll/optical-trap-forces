function angle = th_i_z(r, z, Rsp, focus)
    angle = asin(z / Rsp .* sin(phi_i(r, focus)));
end