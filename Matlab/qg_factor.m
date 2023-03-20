function factor = qg_factor(th, psi, n1, n2)
    R = reflectivity(th, psi, n1, n2);
    T = transmittance(th, psi, n1, n2);
    theta_refl = th_r(th, n1, n2);
    factor = R .* sin(2*th) - T.^2 .* ...
        (sin(2*th - 2*theta_refl) + R .* sin(2*th)) ./ ...
        (1 + R.^2 + 2*R .* cos(2*theta_refl));
end