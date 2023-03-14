# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime".
# There are axial forces only

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import integrate
from scipy import special
from scipy import constants

n1 = 1.3337  # index of refraction of the immersion medium
n2 = 1.4607 # index of refraction of the fused silica at wavelength 523 nm
n = n2 / n1 # n2/n1
NA = 1.25  # numerical aperture
th_max = np.arcsin(NA / n1)  # maximum angle of incidence
f = 2.0e-3  # objective lens focus or WD
r_max = f * np.tan(th_max)  # radius of a Gaussian beam (1:1 with input aperture condition)
Rsp = 1.03e-6  # sphere radius
P = 4.4e-3  # power of the laser


# angle of refraction
def r(th):
    return np.arcsin(n1 / n2 * np.sin(th))


# Fresnel reflectivity
def r_f(th, psi):
    return (np.tan(th - r(th)) ** 2 / np.tan(th + r(th)) ** 2) * np.cos(psi) ** 2 + \
           (np.sin(th - r(th)) ** 2 / np.sin(th + r(th)) ** 2) * np.sin(psi) ** 2


# Fresnel transparency
def t_f(th, psi):
    return 1 - r_f(th, psi)


# force factors
def q_s(th, psi):
    return 1 + r_f(th, psi) * np.cos(2 * th) - t_f(th, psi) ** 2 * \
           (np.cos(2 * th - 2 * r(th)) + r_f(th, psi) * np.cos(2*th)) / \
           (1 + r_f(th, psi) ** 2 + 2 * r_f(th, psi) * np.cos(2*r(th)))


def q_g(th, psi):
    return r_f(th, psi) * np.sin(2 * th) - t_f(th, psi) ** 2 * \
           (np.sin(2 * th - 2 * r(th)) + r_f(th, psi) * np.sin(2 * th)) / \
           (1 + r_f(th, psi) ** 2 + 2 * r_f(th, psi) * np.cos(2 * r(th)))


def q_mag(th, psi):
    return np.sqrt(q_s(th, psi) ** 2 + q_g(th, psi) ** 2)


# Average factors (circular polarization
def q_s_avg(th):
    return 0.5 * (q_s(th, 0) + q_s(th, np.pi/2))


def q_g_avg(th):
    return 0.5 * (q_g(th, 0) + q_g(th, np.pi/2))


def q_mag_avg(th):
    return np.sqrt(q_s_avg(th) ** 2 + q_g_avg(th) ** 2)


# Angles
def phi(dr):
    return np.arctan(dr / f)


def thi(dr, dz):
    return np.arcsin(dz / Rsp * np.sin(phi(dr)), dtype=np.cfloat)


def q_g_z(dr, dz):
    return -q_g_avg(thi(dr, dz)) * np.sin(phi(dr))


def q_s_z(dr, dz):
    return q_s_avg(thi(dr, dz)) * np.cos(phi(dr))


# Intensity profile
a = 1.0
w0 = a * r_max


def intensity_uniform():
    return P / (np.pi * r_max ** 2)


def intensity_gaussian_tem00(dr):
    i_0 = P * 2 / (np.pi*w0 ** 2)
    return i_0 * np.exp(-2 * dr ** 2 / w0 ** 2)


def intensity_bessel(dr):
    w0_bb = 0.5 * r_max
    i_0 = P * 2 / (np.pi * w0 ** 2)
    return i_0 * special.jv(0, 2.405 / w0_bb * dr) ** 2 * np.exp(- 2 * dr ** 2 / w0 ** 2)


# Intensity profile graphics
sns.set()
sns.set_style("darkgrid")

rho = np.linspace(-r_max, r_max, 500)
fig1 = plt.figure(1, figsize=(10, 6))
plt.plot(rho, intensity_gaussian_tem00(rho), 'k')
plt.xlabel('r, m', fontsize=18)
plt.ylabel('I(r)', fontsize=18)


# Integration
def q_res_g(dz, func):
    ans = integrate.quad(lambda x: x * func(x) * q_g_z(x, dz) * (~np.iscomplex(q_g_z(x, dz))).astype(float),
                               0, r_max, epsabs=1e-12, epsrel=1e-6)
    return 2 * np.pi * ans[0]


def q_res_s(dz, func):
    ans = integrate.quad(lambda x: x * func(x) * q_s_z(x, dz) * (~np.iscomplex(q_s_z(x, dz))).astype(float),
                               0, r_max, epsabs=1e-12, epsrel=1e-6)
    return 2 * np.pi * ans[0]


# Calculation
n = 200
z = np.linspace(-2 * Rsp, 2 * Rsp, n)

axial_g = [q_res_g(x, intensity_gaussian_tem00) for x in z]
axial_s = [q_res_s(x, intensity_gaussian_tem00) for x in z]

f_0 = n1 * P / constants.c  # net force

axial_g = np.array(axial_g[::-1])
axial_s = np.array(axial_s[::-1])
axial = axial_g + axial_s
z = -z[::-1]

# Graphics
fig2 = plt.figure(2, figsize=(10, 6))
plt.plot(z, f_0*axial_g, 'b-.', lw=1, label='$F_{g}$')
plt.plot(z, f_0*axial_s, 'r--', lw=1, label='$F_{s}$')
plt.plot(z, f_0*axial, 'k', lw=1, label='$F_{t}$')
plt.xlabel('r, m', fontsize=18)
plt.ylabel('F, N', fontsize=18)
plt.legend(fontsize=18)
plt.show()
