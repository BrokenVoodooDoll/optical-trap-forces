Rsp = 1.03e-6; % sphere radius [m]
n1 = 1.3337; % index of refraction of the immersion medium
n2 = 1.4607; % index of refraction of the fused silica at wavelength 523 nm
n = n2/n1; % relative refraction index
c0 = 3e8; % speed of light [m/s]
NA = 1.25; % numerical apertur  e
f = 2.0e-3; % objective lens focus or WD [m]
P = 51.7e-3; % power of the laser [W]
F0 = n1*P/c0; % resulting force [N]
th_max = asin(NA/n1); % maximum incidence angle

ratio = 1.0; % the ratio of the beam radius to the aperture radius
r_max = f * tan(th_max); % radius of a Gaussian beam (1:1 with input aperture condition)
w0 = ratio * r_max; % Gaussian beam waist radius [m]