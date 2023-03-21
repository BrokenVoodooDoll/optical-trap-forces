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
I = gauss(rho, r_max, w0);
I0 = max(I);
figure
plot(rho, I/I0, 'k')
grid
xlabel('r, m')
ylabel('I(r)')

% Integration
G0 = gauss_peak(r_max, w0);
Qres_g = @(y) G0 * integral2(@(beta, r) r .* gauss(r, r_max, w0) .* ...
    iscomplex(qg_y_factor(beta, r, y, n1, n2, Rsp, f)), 0, 2*pi, 0, r_max, ...
    'Method', 'iterated', 'AbsTol', 1e-8, 'RelTol', 1e-6);

Qres_s = @(y) G0 * integral2(@(beta, r) r .* gauss(r, r_max, w0) .* ...
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
xlabel('r, ì')
ylabel('F, Í')
grid
