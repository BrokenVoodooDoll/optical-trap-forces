function reflectivity = R(th, psi)
    reflectivity = (tan(th - th_r(th)).^2 ./ tan(th + th_r(th)).^2) .* ...
        cos(psi).^2 + ...
        (sin(th - th_r(th)).^2 ./ sin(th + th_r(th)).^2) .* sin(psi).^2;
end
   
function transparency = T(th, psi)
    transparency = 1 - R(th, psi);
end