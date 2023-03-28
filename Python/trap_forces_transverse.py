# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime".
# There are transverse forces only

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functions import *


sns.set()


# Intensity profile graphics
rho = np.linspace(-r_max, r_max, 500)
plt.figure(1, figsize=(10, 6))
I = gauss(rho)
I0 = np.max(I)
plt.plot(rho, I / I0)
plt.fill_between(rho, I / I0, 0, alpha=0.3)
plt.title("Forces along Z axis")
plt.xlabel('$r, m$', fontsize=18)
plt.ylabel('$I(r)$', fontsize=18)


# Integration
def q_res_g(dy, func):
    ans = integrate.dblquad(lambda dr, db, dy: dr * func(dr) * q_g_y(db, dr, dy) * (~np.iscomplex(q_g_y(db, dr, dy))).astype(float),
                            0, 2 * np.pi, 0, r_max, args=(dy, ),
                            epsabs=1e-6, epsrel=1e-6)
    return ans[0] / (np.pi * w0 ** 2)


def q_res_s(dy, func):
    ans = integrate.dblquad(lambda dr, db: dr * func(dr) * q_s_y(db, dr, dy) * (~np.iscomplex(q_g_y(db, dr, dy))).astype(float),
                            0, 2 * np.pi, 0, r_max,
                            epsabs=1e-6, epsrel=1e-6)
    return ans[0] / (np.pi * w0 ** 2)

# Calculation
n = 150
y = np.linspace(-2 * Rsp, 2 * Rsp, n)

transverse_g = F0 * np.abs(np.array([q_res_g(dy, gauss) for dy in y]))
transverse_s = F0 * np.array([q_res_s(dy, gauss) for dy in y])
transverse = transverse_g + transverse_s


# Graphics
plt.figure(2, figsize=(10, 6))
plt.plot(y, transverse_g, '-.', lw=1, label='$F_{g}$')
plt.plot(y, transverse_s, '--', lw=1, label='$F_{s}$')
plt.plot(y, transverse, lw=1, label='$F_{t}$')
plt.xlabel('$r, m$', fontsize=18)
plt.ylabel('$F, N$', fontsize=18)
plt.legend()

plt.show()