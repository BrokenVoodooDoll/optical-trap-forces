function T = transmittance(th, psi, n1, n2)
    T = 1 - reflectivity(th, psi, n1, n2);
end