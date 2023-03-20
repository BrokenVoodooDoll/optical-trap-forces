function factor = qs_factor(th, psi, n1, n2)
    R = reflectivity(th, psi, n1, n2);
    T = transmittance(th, psi, n1, n2);
    theta_refl = th_r(th, n1, n2);

    factor = 1 + R .* cos(2*th) - T.^2.*...
    (cos(2*th - 2*theta_refl) + R .* cos(2*th)) ./ ...
    (1 + R.^2 + 2*R .* cos(2*theta_refl));
end