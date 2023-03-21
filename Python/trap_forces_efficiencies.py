# These calculations are based on Ashkin's article "Forces of a single-beam
# gradient laser trap on a dielectric sphere in the ray optics regime

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functions import *


sns.set()

t = np.linspace(0, np.pi / 2, 1000)
t_deg = t * 180 / np.pi
pol = np.pi / 4

plt.figure(figsize=(13, 8))
plt.plot(t_deg, q_s(t, pol), '--', label='Q_s')
plt.plot(t_deg, -q_g(t, pol), '-.', label='Q_g')
plt.plot(t_deg, q_mag(t, pol), label='Q_t')
plt.grid()
plt.xlabel(r'$\theta$, deg', fontsize=18)
plt.ylabel('Q', fontsize=18)
plt.legend(fontsize=18)
plt.title('Beam efficiencies', fontsize=20)
plt.show()
