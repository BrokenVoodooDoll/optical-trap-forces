%all distances in m
close all
clear
clc
format compact

% These calculations are based on Ashkin's article "Forces of a single-beam
% gradient laser trap on a dielectric sphere in the ray optics regime".
% There are transverse forces only

load_constants

% Intensity profile graphics
rho = linspace(-r_max, r_max, 500);
I = gauss(rho, w0, r_max);
I0 = max(I);
figure
plot(rho, I/I0, 'k')
grid
xlabel('r, m')
ylabel('I(r)')

% Integration
Qres_g = @(y) 1 / (pi * w0^2) * integral2(@(beta, r) r .* gauss(r, w0, r_max) .* ...
    iscomplex(qg_y_factor(beta, r, y, n1, n2, Rsp, f)), 0, 2*pi, 0, r_max, ...
    'Method', 'iterated', 'AbsTol', 1e-8, 'RelTol', 1e-6);

Qres_s = @(y) 1 / (pi * w0^2) * integral2(@(beta, r) r .* gauss(r, w0, r_max) .* ...
    iscomplex(qs_y_factor(beta, r, y, n1, n2, Rsp, f)), 0, 2*pi, 0, r_max, ...
    'Method', 'iterated', 'AbsTol', 1e-8, 'RelTol', 1e-6);

% Calulation
N = 150;
y = linspace(-2*Rsp, 2*Rsp, N);
Transverse_g = zeros(1, N);
Transverse_s = zeros(1, N);

wb = waitbar(0, 'Calculating...');
for ii = 1:N
    Transverse_g(ii) = abs(Qres_g(y(ii)));
    Transverse_s(ii) = Qres_s(y(ii));
    waitbar(ii / N, wb, 'Calculating...');
end
close(wb);

Transverse = Transverse_g + Transverse_s;

%Graphics
figure
plot(y, F0 * Transverse_g,'r--', ...
    y, F0 * Transverse_s, 'b-.', ...
    y, F0*Transverse, 'k')
legend('F_{g}','F_{s}','F_{t}')
xlabel('r, �')
ylabel('F, �')
grid
