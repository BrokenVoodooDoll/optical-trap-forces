close all
clear
format compact
clc

% These calculations are based on Ashkin's article "Forces of a single-beam
% gradient laser trap on a dielectric sphere in the ray optics regime

%all distances in mm
a = 1.0e-6; % radius of the bead
n1 = 1.0; % index of rafraction of the medium
n = 1.4607; % n2/n1
n2 = n*n1; % index of refraction of the fused silica
c0 = 3e8; % speed of light

%reflectivity
R = @(th,psi) (tan(th-asin(n1/n2*sin(th))).^2./...
    tan(th+asin(n1/n2*sin(th))).^2).*cos(psi).^2+...
    (sin(th-asin(n1/n2*sin(th))).^2./...
    sin(th+asin(n1/n2*sin(th))).^2).*sin(psi).^2;

%transparency
T = @(th,psi) 1-R(th,psi);

r = @(th) asin(n1/n2*sin(th));

% Factors
Qs = @(th, psi) 1 + R(th, psi) .* cos(2*th) - T(th,psi).^2 .* (cos(2*th -...
    2*r(th)) + R(th, psi) .* cos(2*th)) ./ (1 + R(th,psi).^2 +...
    2*R(th,psi) .* cos(2*r(th)));

Qg = @(th, psi) R(th, psi) .* sin(2*th) - T(th,psi).^2 .* (sin(2*th -...
    2*r(th)) + R(th, psi) .* sin(2*th)) ./ (1 + R(th,psi).^2 +...
    2*R(th,psi) .* cos(2*r(th)));

Qmag = @(th, psi) sqrt(Qs(th, psi).^2 + Qg(th, psi).^2);

t = linspace(0, pi/2, 1000);
t_deg = t*180/pi;
pol = pi/4;

figure
plot(t_deg, Qs(t, pol),'r--', t_deg, -Qg(t, pol),'b-.', t_deg, Qmag(t, pol),'k');
grid
xlabel('\theta, deg')
ylabel('Q')
legend('Q_s','Q_g','Q_t','location','northwest')
sdf('my')