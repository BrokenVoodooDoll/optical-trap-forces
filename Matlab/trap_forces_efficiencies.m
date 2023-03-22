close all
clear
format compact
clc

% These calculations are based on Ashkin's article "Forces of a single-beam
% gradient laser trap on a dielectric sphere in the ray optics regime

load_constants

t = linspace(0, pi/2, 1000);
t_deg = t * 180/pi;
psi = pi / 4;

figure
plot(t_deg, qs_factor(t, psi, n1, n2), ...
    t_deg, -qg_factor(t, psi, n1, n2), ...
    t_deg, qmag_factor(t, psi, n1, n2));
grid
title(['\psi = ', num2str(psi * 180/pi), '^\circ'])
xlabel('\theta, ^\circ')
ylabel('Q')
legend('Q_s','Q_g','Q_t','location','northwest')