% Factors
function factor = Qs(th, psi)
    factor = 1 + R(th, psi) .* cos(2*th) - T(th, psi).^2.*...
    (cos(2*th - 2*th_r(th)) + R(th, psi) .* cos(2*th)) ./ ...
    (1 + R(th, psi).^2 + 2*R(th, psi) .* cos(2*th_r(th)));
end

function factor = Qg(th, psi)
    factor = R(th, psi) .* sin(2*th) - T(th, psi).^2 .* ...
        (sin(2*th - 2*th_r(th)) + R(th, psi) .* sin(2*th)) ./ ...
        (1 + R(th, psi).^2 + 2*R(th, psi) .* cos(2*th_r(th)));
end

function factor = Qmag(th, psi)
    factor = sqrt(Qs(th, psi).^2 + Qg(th, psi).^2);
end

% Average factors
function factor = Qs_avg(th)
    factor = 0.5*(Qs(th, 0) + Qs(th, pi/2));
end

function factor = Qg_avg(th)
    factor = 0.5*(Qg(th, 0) + Qg(th, pi/2));
end

function factor = Qmag_avg(th)
    factor = sqrt(Qs_avg(th).^2 + Qg_avg(th).^2);
end

function factor = Qgz(r, z)
    factor = -Qg_avg(th_i(r, z)) .* sin(ph_i(r));
end

function factor = Qsz(r, z)
    factor = Qs_avg(th_i(r, z)) .* cos(ph_i(r));
end