"""Picture for theoretical background of the autocorrelation function."""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

### input ###
tmax = 50
tplot = 50
tgrid = 100000
omega = 1
nstates = 20
sigma = 5

colors = ['#701C1C', '#0F599F', '#B0DCE0']

### code ###
tmin = -tmax
t = np.linspace(0, tmax, tgrid)
t2 = np.linspace(tmin, tmax, 2*tgrid)

ho_energies = np.array([0.5*omega + n for n in range(nstates)])
ho_coeffs = ho_energies**3*1/(np.exp(ho_energies) - 1)
ho_coeffs /= np.sqrt(np.sum(ho_coeffs**2))

signal = np.zeros(shape=tgrid, dtype=complex)
signal2 = np.zeros(shape=2*tgrid, dtype=complex)
for state in range(nstates):
    signal += ho_coeffs[state]**2*np.exp(-1j*ho_energies[state]*t)
    signal2 += ho_coeffs[state]**2*np.exp(-1j*ho_energies[state]*t2)

gamma = 0.25
signal_damp = signal*np.exp(-gamma*t**2/sigma**2)
signal2 = signal2*np.exp(-gamma*t2**2/sigma**2)

frequency = np.fft.ifftshift(2*np.pi*np.fft.fftfreq(tgrid, t[1] - t[0]))
spectrum = np.fft.ifftshift(np.fft.ifft(signal))
spectrum /= np.max(np.abs(spectrum))/np.max(ho_coeffs**2)

spectrum_damp = np.fft.ifftshift(np.fft.ifft(signal_damp))
spectrum_damp /= np.max(np.abs(spectrum_damp))/np.max(ho_coeffs**2)

frequency2 = np.fft.ifftshift(2*np.pi*np.fft.fftfreq(2*tgrid, t2[1] - t2[0]))
spectrum2 = np.fft.ifftshift(np.fft.ifft(signal2))
spectrum2 /= np.max(np.abs(spectrum2))/np.max(ho_coeffs**2)

# plotting
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.))

axs[0].plot(t, np.real(signal), color=colors[1], ls='-', label=r"$S_0^{50}(t)$", zorder=-1)
axs[0].plot(t, np.real(signal_damp), color=colors[2], ls='-', label=r"$\xi(t)S_0^{50}(t)$", zorder=1)
axs[0].plot(t2, np.real(signal2), color=colors[0], ls='-', label=r"$\xi(t)S_{-50}^{50}(t)$", zorder=0)


axs[0].set_title(r"Autocorrelation function")
axs[0].set_xlabel(r"$t$")
axs[0].set_ylabel(r"$\mathrm{Re}[S(t)]$")

axs[0].set_xticks(np.arange(-100, 100, 10))
axs[0].set_xlim(-tplot, tplot)
axs[0].set_ylim(-1.1, 1.1)

axs[0].tick_params('both', direction='in', which='both', top=True, right=True)
axs[0].minorticks_on()
axs[0].legend(frameon=False, loc='lower left', handlelength=1.8, labelspacing=0.5)
axs[0].axhline(y=0, color='black', linewidth=0.5)
axs[0].axvline(x=0, color='black', linewidth=0.5)

axs[1].plot(frequency, np.abs(spectrum), color=colors[1], ls='-')
axs[1].plot(frequency, np.abs(spectrum_damp), color=colors[2], ls='-')
axs[1].plot(frequency2, np.abs(spectrum2), color=colors[0], ls='-')

axs[1].fill_between(frequency2, np.abs(spectrum2), 0, alpha=0.1, color=colors[0])
for state in range(nstates):
    axs[1].axvline(ho_energies[state], color='black', linewidth=0.5, alpha=0.2)

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
plt.savefig('autocorr3', dpi=300)
plt.show()
