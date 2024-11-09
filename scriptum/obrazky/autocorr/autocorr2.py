"""Picture for theoretical background of the autocorrelation function."""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

### input ###
tmax = 1000
tplot = 18
tgrid = 100000
omega = 1
nstates = 20
sigma = 5

colors = ['#701C1C', '#0F599F']

### code ###
tmin = -tmax
t = np.linspace(tmin, tmax, tgrid)

ho_energies = np.array([0.5*omega + n for n in range(nstates)])
ho_coeffs = ho_energies**3*1/(np.exp(ho_energies) - 1)
ho_coeffs /= np.sqrt(np.sum(ho_coeffs**2))

signal = np.zeros(shape=tgrid, dtype=complex)
for state in range(nstates):
    signal += ho_coeffs[state]**2*np.exp(-1j*ho_energies[state]*t)

signal_damp =signal* np.exp(-0.5*t**2/sigma**2)

frequency = np.fft.ifftshift(2*np.pi*np.fft.fftfreq(tgrid, t[1] - t[0]))
spectrum = np.fft.ifftshift(np.fft.ifft(signal))
spectrum /= np.max(np.abs(spectrum))/np.max(ho_coeffs**2)

spectrum_damp = np.fft.ifftshift(np.fft.ifft(signal_damp))
spectrum_damp /= np.max(np.abs(spectrum_damp))/np.max(ho_coeffs**2)

# plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.))

axs[0].plot(t, np.real(signal), color=colors[1], ls='--', alpha=0.5, label=r"$S(t)$")
axs[0].plot(t, np.real(signal_damp), color=colors[1], ls='-', label=r"$\xi(t)S(t)$")
axs[0].fill_between(t, np.exp(-0.5*t**2/sigma**2), 0, alpha=0.1, color=colors[1])

axs[0].set_title(r"Autocorrelation function")
axs[0].set_xlabel(r"$t$")
axs[0].set_ylabel(r"$\mathrm{Re}[S(t)]$")

axs[0].set_xticks(np.arange(-25, 25, 5))
axs[0].set_xlim(-tplot, tplot)
axs[0].set_ylim(-1.1, 1.1)

axs[0].tick_params('both', direction='in', which='both', top=True, right=True)
axs[0].minorticks_on()
axs[0].legend(frameon=True, edgecolor='white',loc='lower right', handlelength=1.8, labelspacing=0)
axs[0].axhline(y=0, color='black', linewidth=0.5)
axs[0].axvline(x=0, color='black', linewidth=0.5)

# axs[1].plot(frequency, np.abs(spectrum), color=colors[1], ls='--')
axs[1].plot(frequency, np.abs(spectrum_damp), color=colors[1], ls='-')
axs[1].fill_between(frequency, np.abs(spectrum_damp), 0, alpha=0.1, color=colors[1])
for state in range(nstates):
    axs[1].plot([ho_energies[state]]*2, [0, ho_coeffs[state]**2], color=colors[1], ls='--', alpha=0.5)
# axs[1].scatter(ho_energies, ho_coeffs**2, color='black', ls='-', marker='x', zorder=2)

axs[1].set_title(r"Spectrum")
axs[1].set_xlabel(r"$E$")
axs[1].set_ylabel(r"$\sigma(E)$")

axs[1].set_xticks([0] + list(ho_energies), ['0'] + [f"$E_{e:d}$" for e in range(nstates)])
axs[1].set_xlim(0, 0.5*omega*nstates)
axs[1].set_ylim(0)

axs[1].tick_params('both', direction='in', which='both', top=True, right=True)
axs[1].minorticks_on()
axs[1].xaxis.set_minor_locator(ticker.MultipleLocator(omega/4))  # minor ticks separator

plt.tight_layout()
plt.savefig('autocorr2', dpi=300)
plt.show()
