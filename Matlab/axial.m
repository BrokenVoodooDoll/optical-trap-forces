clear
close all
clc
format compact

% These calculations are based on Ashkin's article "Forces of a single-beam
% gradient laser trap on a dielectric sphere in the ray optics regime".
% There are axial forces only

load_constants

% Intensity profile plots
rho = linspace(-r_max, r_max, 500);
figure
I = gauss(rho, w0, r_max);
I0 = max(I);
plot(rho, I/I0, 'k')
grid
xlabel('r, м')
ylabel('I(r)')

% Integration
Qres_g = @(z) 1 / (pi * w0 ^ 2) * integral2(@(beta, r) r .* gauss(r, w0, r_max) .* ...
    iscomplex(qg_z_factor(r, z, n1, n2, Rsp, f)), 0, 2*pi, 0, r_max, ...
    'Method', 'iterated', 'AbsTol', 1e-12, 'RelTol', 1e-6);

Qres_s = @(z) 1 / (pi * w0 ^ 2) * integral2(@(beta, r) r .* gauss(r, w0, r_max) .* ...
    iscomplex(qs_z_factor(r, z, n1, n2, Rsp, f)), 0, 2*pi, 0, r_max, ...
    'Method', 'iterated', 'AbsTol', 1e-12, 'RelTol', 1e-6);

% Calulation
N = 200;
z = linspace(-2*Rsp, 2*Rsp, N);
Axial_g = zeros(1, N);
Axial_s = zeros(1, N);

wb = waitbar(0, 'Calculating...');
for ii = 1:N
    Axial_g(ii) = Qres_g(z(ii));
    Axial_s(ii) = Qres_s(z(ii));
    waitbar(ii / N, wb, 'Calculating...');
end
close(wb);

Axial_g = fliplr(Axial_g);
Axial_s = fliplr(Axial_s);
Axial = Axial_g + Axial_s;
z = -fliplr(z);

% Plots
figure
plot(z, F0 * Axial_g, 'b-.', ...
    z, F0 * Axial_s, 'r--', ...
    z, F0 * Axial, 'k')
legend('F_{g}','F_{s}','F_{t}')
xlabel('r, м')
ylabel('F, Н')
grid