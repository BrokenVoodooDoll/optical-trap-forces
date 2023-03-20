function factor = qmag_factor(th, psi, n1, n2)
    factor = sqrt(qs_factor(th, psi, n1, n2).^2 + qg_factor(th, psi, n1, n2).^2);
end