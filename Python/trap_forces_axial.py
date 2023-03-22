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
    ans = integrate.quad(lambda dr: dr * func(dr) * q_g_z(dr, z) * (~np.iscomplex(q_g_z(dr, z))).astype(float),
                         0, r_max, epsabs=1e-12, epsrel=1e-6)
    return 2 * np.pi * ans[0] / (np.pi * w0 ** 2)

def q_res_s(z, func):
    ans = integrate.quad(lambda dr: dr * func(dr) * q_s_z(dr, z) * (~np.iscomplex(q_s_z(dr, z))).astype(float),
                         0, r_max, epsabs=1e-12, epsrel=1e-6)
    return 2 * np.pi * ans[0] / (np.pi * w0 ** 2)


# Calculation
n = 200
z = np.linspace(-2 * Rsp, 2 * Rsp, n)

axial_g = [q_res_g(dz, gauss) for dz in z]
axial_s = [q_res_s(dz, gauss) for dz in z]

axial_g = F0 * np.array(axial_g[::-1])
axial_s = F0 * np.array(axial_s[::-1])
axial = axial_g + axial_s
z = -z[::-1]

# Graphics
fig2 = plt.figure(2, figsize=(10, 6))
plt.plot(z, axial_g, '-.', lw=1, label='$F_{g}$')
plt.plot(z, axial_s, '--', lw=1, label='$F_{s}$')
plt.plot(z, axial, lw=1, label='$F_{t}$')
plt.xlabel('r, m', fontsize=18)
plt.ylabel('F, N', fontsize=18)
plt.legend(fontsize=18)
plt.show()
