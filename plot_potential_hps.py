#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

def hps(x, sigma, eps, lam):
    sig_x_6 = (sigma / x) ** 6
    lj = 4.0 * eps * (sig_x_6**2 - sig_x_6)
    if x < 2**(1/6) * sigma:
        return lj + (1-lam)*eps
    else:
        return lam * lj

sigma = 5.0
eps = 2.0
fig, ax = plt.subplots()

x = np.arange(0.0, 10.0, 0.1)

# Plot for five differnet lambda values
for lam in [0.0, 0.2, 0.5, 0.8, 1.0]:
    ax.plot(x, [hps(r, sigma, eps, lam) for r in x], label=f'$\lambda = ${lam}')

# Vertical line for the conversion point
ax.axvline(x=2**(1/6)*sigma, ymin=0, ymax=1, ls=':', c='k')

ax.set_xlim(3.5, 10)
ax.set_ylim(-3, 8)

ax.set_ylabel(r'$\Phi(r)$ / kcal/mol')
ax.set_xlabel(r'$r$ / angstrom')
ax.set_title(r'$\sigma = $' + f'{sigma}' + r', $\varepsilon = $' + f'{eps}')

ax.legend()
plt.show()
