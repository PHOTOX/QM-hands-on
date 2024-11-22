"""Picture for theoretical background of the imaginary-time quantum dynamics."""

import matplotlib.pyplot as plt
import numpy as np

### input ###
tplot = 5
tgrid = 1000
omega = 1
nstates = 10

colors = plt.cm.Blues(np.linspace(1, 0.2, nstates))

### code ###
tmax = tplot
t = np.linspace(0, tmax, tgrid)

ho_energies = np.array([0.5*omega + n for n in range(nstates)])
c0 = (np.exp(-ho_energies/8))
c0 /= np.sqrt(np.sum(c0**2))
c = np.zeros(shape=(nstates, tgrid))

for i in range(nstates):
    c[i] = c0[i]*np.exp(-ho_energies[i]*t)

norm = np.sqrt(np.sum(c**2, axis=0))

# plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.))

for state in range(nstates):
    axs[0].plot(t, c[state], color=colors[state], ls='-', label=f"$c_{state:d}$")
    axs[1].plot(t, c[state]/norm, color=colors[state], ls='-', label=r"$\tilde{c}$" + f"$_{state:d}$")

axs[0].set_title(r"Coefficients")
axs[0].set_xlabel(r"$t$")
axs[0].set_ylabel(r"$c_i(t)$")

axs[0].set_xlim(0, tplot)
axs[0].set_ylim(0)

axs[0].tick_params('both', direction='in', which='both', top=True, right=True)
axs[0].minorticks_on()
axs[0].legend(frameon=False, labelspacing=0, loc='upper right')

axs[1].set_title(r"Normalized coefficients")
axs[1].set_xlabel(r"$t$")
axs[1].set_ylabel(r"$\tilde{c}_i(t)$")

axs[1].set_xlim(0, tplot)
axs[1].set_ylim(0)

axs[1].tick_params('both', direction='in', which='both', top=True, right=True)
axs[1].minorticks_on()
axs[1].legend(frameon=False, labelspacing=0, loc='upper right')

plt.tight_layout()
plt.savefig('itqd1', dpi=300)
plt.show()
