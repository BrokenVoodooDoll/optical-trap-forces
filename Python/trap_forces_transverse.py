# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime".
# There are transverse forces only

import numpy as np
import matplotlib.pyplot as plt
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
P = 14e-3  # power of the laser


# Angles
def r(th):
    return np.arcsin(n1 / n2 * np.sin(th, dtype=np.cfloat))


def phi(dr):
    return np.arctan(dr / f, dtype=np.cfloat)


def gamma(db, dr):
    return np.arccos(np.cos(np.pi / 2 - phi(dr)) * np.cos(db), dtype=np.cfloat)


def thi(db, dr, dy):
    return np.arcsin(dy / Rsp * np.sin(gamma(db, dr)), dtype=np.cfloat)


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


# Average factors (circular polarization)
def q_s_avg(th):
    return 0.5 * (q_s(th, 0) + q_s(th, np.pi/2))


def q_g_avg(th):
    return 0.5 * (q_g(th, 0) + q_g(th, np.pi/2))


def q_mag_avg(th):
    return np.sqrt(q_s_avg(th) ** 2 + q_g_avg(th) ** 2)


def q_g_z(db, dr, dy):
    return q_g_avg(thi(db, dr, dy)) * np.cos(phi(dr), dtype=np.cfloat)


def q_s_z(db, dr, dy):
    return q_s_avg(thi(db, dr, dy)) * np.sin(gamma(db, dr), dtype=np.cfloat)


# Intensity profile
a = 1.0
w0 = a * r_max


def intensity_uniform(dr):
    return 1 / (np.pi * r_max ** 2)


def intensity_gaussian_tem00(dr):
    amp = (1 - np.exp(-2 * r_max ** 2 / w0 ** 2))
    p_0 = np.pi * w0 ** 2 / 2
    return 1 / (amp * p_0) * np.exp(-2 * dr ** 2 / w0 ** 2)


def intensity_bessel(dr):
    p_0 = np.exp(0.5) * w0 * 0.0025 / 4.81

    def temp_f(w):
        return w * special.jv(0, 2.405 / w0 * w) ** 2
    amp = integrate.quad(temp_f, 0, r_max)
    return 1 / (2 * np.pi * amp[0]) * special.jv(0, 2.405 / w0 * dr) ** 2


# Intensity profile graphics
rho = np.linspace(-r_max, r_max, 500)
fig1 = plt.figure(1, figsize=(10, 6))
plt.plot(rho, intensity_gaussian_tem00(rho), 'k')
plt.grid()
plt.xlabel('r, m', fontsize=18)
plt.ylabel('I(r)', fontsize=18)
plt.show()


# Integration
def q_res_g(dy, func):
    ans = integrate.dblquad(lambda dr, db, dy: dr * func(dr) * q_g_z(db, dr, dy) * (~np.iscomplex(q_g_z(db, dr, P))).astype(float),
                            0, 2 * np.pi, 0, r_max, args=(dy, ),
                            epsabs=1e-4, epsrel=1e-6)
    return ans[0]


def q_res_s(dy, func):
    ans = integrate.dblquad(lambda dr, db: dr * func(dr) * q_s_z(db, dr, dy) * (~np.iscomplex(q_g_z(db, dr, P))).astype(float),
                            0, 2 * np.pi, 0, r_max,
                            epsabs=1e-4, epsrel=1e-6)
    return ans[0]


# Calculation
n = 150
y = np.linspace(-2 * Rsp, 2 * Rsp, n)

transverse_g = np.abs(np.array([q_res_g(x, intensity_gaussian_tem00) for x in y]))
transverse_s = np.array([q_res_s(x, intensity_gaussian_tem00) for x in y])

f_0 = n1 * P / constants.c  # net force

transverse = transverse_g + transverse_s


# Graphics
fig2 = plt.figure(2, figsize=(10, 6))
plt.plot(y, transverse_g, 'b-.', lw=1, label='$F_{g}$')
plt.plot(y, transverse_s, 'r--', lw=1, label='$F_{s}$')
plt.plot(y, transverse, 'k', lw=1, label='$F_{t}$')
plt.xlabel('r, m', fontsize=18)
plt.ylabel('F, m', fontsize=18)
plt.legend()
plt.grid()
plt.show()