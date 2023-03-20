function factor = qg_avg_factor(th, n1, n2)
    factor = 0.5*(qg_factor(th, 0, n1, n2) + qg_factor(th, pi/2, n1, n2));
end