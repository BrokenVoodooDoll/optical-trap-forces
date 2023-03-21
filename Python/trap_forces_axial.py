# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime".
# There are axial forces only

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functions import *
from scipy import integrate


sns.set()


# Intensity profile graphics
rho = np.linspace(-r_max, r_max, 500)
fig1 = plt.figure(1, figsize=(10, 6))
I = gauss(rho)
I0 = np.max(I)
plt.plot(rho, I / I0)
plt.fill_between(rho, I / I0, 0, alpha=0.3)
plt.xlabel('r, m', fontsize=18)
plt.ylabel('I(r)', fontsize=18)

# Integration
def q_res_g(z, func):
    ans = integrate.quad(lambda x: x * func(x) * q_g_z(x, z) * (~np.iscomplex(q_g_z(x, z))).astype(float),
                         0, r_max, epsabs=1e-12, epsrel=1e-6)
    return ans[0]


def q_res_s(z, func):
    ans = integrate.quad(lambda x: x * func(x) * q_s_z(x, z) * (~np.iscomplex(q_s_z(x, z))).astype(float),
                         0, r_max, epsabs=1e-12, epsrel=1e-6)
    return ans[0]


# Calculation
n = 200
z = np.linspace(-2 * Rsp, 2 * Rsp, n)

axial_g = gauss_peak() * [q_res_g(x, gauss) for x in z]
axial_s = gauss_peak() * [q_res_s(x, gauss) for x in z]

f_0 = n1 * P / constants.c  # net force

axial_g = np.array(axial_g[::-1])
axial_s = np.array(axial_s[::-1])
axial = axial_g + axial_s
z = -z[::-1]

# Graphics
fig2 = plt.figure(2, figsize=(10, 6))
plt.plot(z, f_0*axial_g, '-.', lw=1, label='$F_{g}$')
plt.plot(z, f_0*axial_s, '--', lw=1, label='$F_{s}$')
plt.plot(z, f_0*axial, lw=1, label='$F_{t}$')
plt.xlabel('r, m', fontsize=18)
plt.ylabel('F, N', fontsize=18)
plt.legend(fontsize=18)
plt.show()
