function factor = qmag_avg_factor(th, n1, n2)
    factor = sqrt(qs_avg_factor(th, n1, n2).^2 + qg_avg_factor(th, n1, n2).^2);
end