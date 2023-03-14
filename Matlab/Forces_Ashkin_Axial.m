%all distances in m
clear
close all
clc
format compact

% These calculations are based on Ashkin's article "Forces of a single-beam
% gradient laser trap on a dielectric sphere in the ray optics regime".
% There are axial forces only

n1 = 1.33; % index of refraction of the immersion medium
n2 = 1.6; % index of refraction of the fused silica at wavelength 523 nm
n = n2/n1; % n2/n1
c0 = 3e8; % speed of light
NA = 1.25; % numerical aperture
th_max = asin(NA/n1); % maximum angle of incidence
f = 100.0e-3; % objective lens focus or WD
r_max = f*tan(th_max); % radius of a Gaussian beam (1:1 with input aperture condition)
Rsp = 1.0e-6; % sphere radius
P = 20.0e-3; % power of the laser

thr = @(th) asin(n1/n2*sin(th)); % refraction angle

%reflectivity
R = @(th,psi) (tan(th-thr(th)).^2./tan(th+thr(th)).^2).*cos(psi).^2+...
    (sin(th-thr(th)).^2./sin(th+thr(th)).^2).*sin(psi).^2;

%transparency
T = @(th,psi) 1-R(th,psi);

% Factors
Qs = @(th, psi) 1 + R(th, psi) .* cos(2*th) - T(th,psi).^2 .* (cos(2*th -...
    2*thr(th)) + R(th, psi) .* cos(2*th)) ./ (1 + R(th,psi).^2 +...
    2*R(th,psi) .* cos(2*thr(th)));

Qg = @(th, psi) R(th, psi) .* sin(2*th) - T(th,psi).^2 .* (sin(2*th -...
    2*thr(th)) + R(th, psi) .* sin(2*th)) ./ (1 + R(th,psi).^2 +...
    2*R(th,psi) .* cos(2*thr(th)));

Qmag = @(th, psi) sqrt(Qs(th, psi).^2 + Qg(th, psi).^2);

% Average factors (circular polarization
Qs_avg = @(th) 0.5*(Qs(th, 0) + Qs(th, pi/2));
Qg_avg = @(th) 0.5*(Qg(th, 0) + Qg(th, pi/2));
Qmag_avg = @(th) sqrt(Qs_avg(th).^2 + Qg_avg(th).^2);

% Angles
phi = @(r) atan(r/f); 
thi = @(r,z) asin(z/Rsp.*sin(phi(r))); 

Qgz = @(r,z) -Qg_avg(thi(r,z)).*sin(phi(r));
Qsz = @(r,z) Qs_avg(thi(r,z)).*cos(phi(r));

% Intensity profile
a = 1.0;
w0 = a*r_max;

%I = @(r) P/(pi*r_max^2); % uniform distribution

A = (1-exp(-2*r_max.^2/w0^2));
I0 = P*2/(pi*w0^2)/A;
I = @(r) I0*exp(-2*r.^2/w0^2); % Gaussian TEM00 beam

%A = 2*pi*integral(@(r) r.*besselj(0,2.405/w0*r).^2,0,r_max)/P0;
%w0_bb = 0.5*r_max;
%I0 = P*2/(pi*w0^2);
%I = @(r) I0*besselj(0,2.405/w0_bb*r).^2.*exp(-2*r.^2/w0^2); % Bessel beam

% Intensity profile graphics
rho = linspace(-r_max, r_max, 500);
figure
plot(rho, I(rho)/max(I(rho)),'k')
grid
xlabel('r, ì')
ylabel('I(r)')
sdf('my')

% Integration
Qres_g = @(z) 1/(pi*r_max^2)*2/(A*pi*w0^2)*integral2(@(beta,r) r.*I(r).*...
    iscomplex(Qgz(r,z)),0,2*pi,0,r_max,...
    'Method','iterated','AbsTol',1e-12,'RelTol',1e-6);

Qres_s = @(z) 1/(pi*r_max^2)*2/(A*pi*w0^2)*integral2(@(beta,r) r.*I(r).*...
    iscomplex(Qsz(r,z)),0,2*pi,0,r_max,...
    'Method','iterated','AbsTol',1e-12,'RelTol',1e-6);

% Calulation
N = 200;
z = linspace(-2*Rsp,2*Rsp,N);
Axial_g = zeros(1,N);
Axial_s = zeros(1,N);

for ii = 1:N
    Axial_g(ii) = Qres_g(z(ii));
    Axial_s(ii) = Qres_s(z(ii));
end

F0 = n1*P/c0; % net force

Axial_g = fliplr(Axial_g);
Axial_s = fliplr(Axial_s);
Axial = Axial_g + Axial_s;
z = -fliplr(z);

%Graphics
figure
plot(z,F0*Axial_g,'b-.',z,F0*Axial_s,'r--',z,F0*Axial,'k')
legend('F_{g}','F_{s}','F_{t}')
xlabel('r, ì')
ylabel('F, Í')
grid
sdf('my')