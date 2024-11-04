"""Picture for theoretical background of the autocorrelation function."""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

### input ###
tplot = 15
tgrid = 100000
omega = 1
nstates = 10

colors = ['#701C1C', '#0F599F']

### code ###
tmax = tplot
tmin = -tmax
t = np.linspace(tmin, tmax, tgrid)

ho_energies = np.array([0.5*omega + n for n in range(nstates)])
ho_coeffs = ho_energies**3*1/(np.exp(ho_energies) - 1)
ho_coeffs /= np.sqrt(np.sum(ho_coeffs**2))

signal = np.zeros(shape=tgrid, dtype=complex)
for state in range(nstates):
    signal += ho_coeffs[state]**2*np.exp(-1j*ho_energies[state]*t)

# plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.))

axs[0].plot(t, np.real(signal), color=colors[1], ls='-', label=r"$\mathrm{Re}[S(t)]$")
axs[0].plot(t, np.imag(signal), color=colors[0], ls='--', label=r"$\mathrm{Im}[S(t)]$")

axs[0].set_title(r"Autocorrelation function")
axs[0].set_xlabel(r"$t$")
axs[0].set_ylabel(r"$S(t)$")

axs[0].set_xticks(np.arange(-int(tplot), int(tplot) + 10, 5))
axs[0].set_xlim(-tplot, tplot)
axs[0].set_ylim(-1.1, 1.55)

axs[0].tick_params('both', direction='in', which='both', top=True, right=True)
axs[0].minorticks_on()
axs[0].legend(frameon=False, labelspacing=0, loc='upper right')
axs[0].axhline(y=0, color='black', linewidth=0.5)
axs[0].axvline(x=0, color='black', linewidth=0.5)

for state in range(nstates):
    axs[1].plot([ho_energies[state]]*2, [0, ho_coeffs[state]**2], color=colors[1], ls='-')
    axs[1].text(x=ho_energies[state]+omega*0.2, y=ho_coeffs[state]**2+0.02, s=f"$|c_{state:d}|^2$", ha='center')

axs[1].scatter(ho_energies, ho_coeffs**2, color='black', ls='-', marker='x', zorder=2)

axs[1].set_title(r"Spectrum")
axs[1].set_xlabel(r"$E$")
axs[1].set_ylabel(r"$\sigma(E)$")

axs[1].set_xticks([0] + list(ho_energies), ['0'] + [f"$E_{e:d}$" for e in range(nstates)])
axs[1].set_xlim(0, omega*(nstates+1/2))
axs[1].set_ylim(0, 1.2*np.max(ho_coeffs**2))

axs[1].tick_params('both', direction='in', which='both', top=True, right=True)
axs[1].minorticks_on()
axs[1].xaxis.set_minor_locator(ticker.MultipleLocator(omega/4))  # minor ticks separator

plt.tight_layout()
plt.savefig('autocorr1', dpi=300)
plt.show()
