# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime

# all distances in mm

import numpy as np
import matplotlib.pyplot as plt


a = 1.0 * 1e-6  # radius of the bead
n1 = 1.33  # index of refraction of the medium
n = 1.8  # n2/n1
n2 = n * n1  # index of refraction of the fused silica
c0 = 3 * 1e8  # speed of light


# Fresnel reflectivity
def r_f(th, psi):
    return (np.tan(th - np.arcsin(n1 / n2 * np.sin(th))) ** 2 /
            np.tan(th + np.arcsin(n1 / n2 * np.sin(th))) ** 2) * np.cos(psi) ** 2 + \
           (np.sin(th - np.arcsin(n1 / n2 * np.sin(th))) ** 2 /
            np.sin(th + np.arcsin(n1 / n2 * np.sin(th))) ** 2) * np.sin(psi) ** 2


# Fresnel transparency
def t_f(th, psi):
    return 1 - r_f(th, psi)


# angle of refraction
def r(th):
    return np.arcsin(n1 / n2 * np.sin(th))


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


t = np.linspace(0, np.pi / 2, 1000)
t_deg = t * 180 / np.pi
pol = np.pi / 4

plt.figure(figsize=(13, 8))
plt.plot(t_deg, q_s(t, pol), 'r--', label='Q_s')
plt.plot(t_deg, -q_g(t, pol), 'b-.', label='Q_g')
plt.plot(t_deg, q_mag(t, pol), 'k', label='Q_t')
plt.grid()
plt.xlabel(r'$\theta$, deg', fontsize=18)
plt.ylabel('Q', fontsize=18)
plt.legend(fontsize=18)
plt.title('Beam efficiencies', fontsize=20)
plt.show()
