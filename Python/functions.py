import numpy as np
from scipy import special
from scipy import integrate
from scipy import constants

n1 = 1.3337  # index of refraction of the immersion medium
n2 = 1.4607  # index of refraction of the fused silica at wavelength 523 nm
NA = 1.25  # numerical aperture
f = 2.0e-3  # objective lens focus or WD
Rsp = 1.03e-6  # sphere radius
P = 4.4e-3  # power of the laser
ratio = 1.0  # the ratio of the beam radius to the aperture radius

th_max = np.arcsin(NA / n1)  # maximum angle of incidence
# radius of a Gaussian beam (1:1 with input aperture condition)
r_max = f * np.tan(th_max)
n = n2 / n1  # n2/n1
w0 = ratio * r_max  # beam radius
F0 = n1 * P / constants.speed_of_light; # resulting force


# Angle of refraction
def theta_r(th):
    return np.arcsin(n1 / n2 * np.sin(th))


# Fresnel reflectivity
def reflectivity(th, psi):
    theta = theta_r(th)
    return (np.tan(th - theta) ** 2 / np.tan(th + theta) ** 2) * np.cos(psi) ** 2 + \
           (np.sin(th - theta) ** 2 / np.sin(th + theta) ** 2) * np.sin(psi) ** 2


# Fresnel transparency
def transparency(th, psi):
    return 1 - reflectivity(th, psi)


# force factors
def q_s(th, psi):
    R = reflectivity(th, psi)
    T = transparency(th, psi)
    theta = theta_r(th)
    return 1 + R * np.cos(2 * th) - T ** 2 * \
        (np.cos(2 * th - 2 * theta) + R * np.cos(2*th)) / \
        (1 + R ** 2 + 2 * R * np.cos(2*theta))


def q_g(th, psi):
    R = reflectivity(th, psi)
    T = transparency(th, psi)
    theta = theta_r(th)
    return R * np.sin(2 * th) - T ** 2 * \
        (np.sin(2 * th - 2 * theta) + R * np.sin(2 * th)) / \
        (1 + R ** 2 + 2 * R * np.cos(2 * theta))


def q_mag(th, psi):
    return np.sqrt(q_s(th, psi) ** 2 + q_g(th, psi) ** 2)


# Average factors (circular polarization)
def q_s_avg(th):
    return 0.5 * (q_s(th, 0) + q_s(th, np.pi/2))


def q_g_avg(th):
    return 0.5 * (q_g(th, 0) + q_g(th, np.pi/2))


def q_mag_avg(th):
    return np.sqrt(q_s_avg(th) ** 2 + q_g_avg(th) ** 2)


# Angles
def phi_i(r):
    return np.arctan(r / f)


def gamma(beta, r):
    return np.arccos(np.cos(np.pi / 2 - phi_i(r)) * np.cos(beta), dtype=np.cfloat)


def th_i_z(r, z):
    return np.arcsin(z / Rsp * np.sin(phi_i(r)), dtype=np.cfloat)


def th_i_y(beta, r, y):
    return np.arcsin(y / Rsp * np.sin(gamma(beta, r)), dtype=np.cfloat)


def q_g_z(r, z):
    return -q_g_avg(th_i_z(r, z)) * np.sin(phi_i(r))


def q_s_z(r, z):
    return q_s_avg(th_i_z(r, z)) * np.cos(phi_i(r))


def q_g_y(beta, r, y):
    return q_g_avg(th_i_y(beta, r, y)) * np.cos(phi_i(r), dtype=np.cfloat)


def q_s_y(beta, r, y):
    return q_s_avg(th_i_y(beta, r, y)) * np.sin(gamma(beta, r), dtype=np.cfloat)


# Intensity distributions
def gauss_peak():
    A = (1 - np.exp(-2*r_max ** 2 / w0 ** 2))
    return 2*A / (np.pi * w0 ** 2)


def gauss(r):
    return np.exp(-2 * r ** 2 / w0 ** 2)


def bessel(r):
    return special.jv(0, 2.405 / w0 * r) ** 2


def bessel_peak():
    return 4.81 / (w0 * r_max * np.exp(0.5)) * 2 * np.pi * integrate.quad(lambda r: r * special.jv(0, 2.405 / w0 * r) ** 2,
                                                                          0, r_max, epsabs=1e-12, epsrel=1e-6)[0]


def uniform(r):
    return P / (np.pi * r_max ** 2) * np.ones(r.shape)
