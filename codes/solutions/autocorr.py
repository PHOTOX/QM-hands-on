"""Solution for chapter Energy spectrum and autocorrelation function."""

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -25, 25  # minimum and maximum of x
simtime = 10000.0  # simulation time in atomic time units
dt = 0.2  # time step in atomic time units
m = 1.0  # mass in atomic units
gamma = 1/200  # damping factor for autocorrelation function [exp(-kappa*t)]

# physical constants in atomic units
hbar = 1.0

# generate equidistant x grid from xmin to xmax
x = np.linspace(xmin, xmax, ngrid)

# generate momentum grid for the discrete Fourier transform
dx = (xmax - xmin)/(ngrid - 1)
k = 2*np.pi*np.fft.fftfreq(ngrid, d=dx)

# generate potential V
V = 0.005*x**2

# generate kinetic energy T in momentum space
T = hbar**2*k**2/2/m

# generate V propagator with time step dt/2
expV = np.exp(-1j*V*dt/2)

# generate T propagator in momentum space with time step dt
expT = np.exp(-1j*T*dt)

# initiate wave function; be careful, it must be a complex numpy array
alpha, x0, p0 = 0.3, 5.7, 0.0
psi = (2*alpha/np.pi)**0.25*np.exp(-alpha*(x - x0)**2 + 1j*p0*(x - x0), dtype=complex)
psi /= np.sqrt(np.real(np.trapz(y=np.conjugate(psi)*psi, x=x)))

# save the initial wave function for calculating the autocorrelation function
psi0 = psi  # initial wave function
time, autocorr = [], []  # empty lists for appending values of time and autocorrelation function

# propagation
t = 0  # set initial time to 0
print("\nLaunching quantum dynamics.\n---------------------------\n")
while t < simtime:  # loop until simulation time is reached
    # propagate half-step in V
    psi *= expV

    # propagate full step in T
    # first the inverse Fourier transform to momentum space
    psi_k = np.fft.ifft(psi, norm="ortho")
    # apply expT in momentum space
    psi_k *= expT
    # finally the Fourier transform back to coordinate space
    psi = np.fft.fft(psi_k, norm="ortho")

    # propagate half-step in V for the second time
    psi *= expV

    # calculate new time after propagation
    t += dt

    # calculate norm
    norm = np.real(np.trapz(y=np.conjugate(psi)*psi, x=x))

    # check that norm is conserved
    if np.abs(norm - 1) > 1e-5:
        print(f"ERROR: Norm ({norm:.9f}) is not conserved!")
        exit(1)

    # calculate expectation value of energy
    # potential energy <V>
    energyV = np.real(np.trapz(y=np.conjugate(psi)*V*psi, x=x))/norm
    # kinetic energy <T>, T operator will be again applied in the momentum space
    psi_k = np.fft.ifft(psi, norm="ortho")
    psi_t = np.fft.fft(T*psi_k, norm="ortho")
    energyT = np.real(np.trapz(y=np.conjugate(psi)*psi_t, x=x))/norm
    # total energy <E>
    energy = energyV + energyT

    # calculate the autocorrelation function S(t) = <psi(0)|psi(t)>
    overlap = np.trapz(y=np.conjugate(psi0)*psi, x=x)

    autocorr.append(overlap)  # appending the overlap to our autocorrelation function list
    time.append(t)  # appending t to our time list

    # print simulation data
    print(f"--Time: {t:.2f} a.t.u.")
    print(f"  <psi|psi> = {norm:.8f}, <E> = {energy:.6f}, <V> = {energyV:.6f}, <T> = {energyT:.6f}")

### autocorrelation function section ###
autocorr = np.array(autocorr) # converting the autocorrelation function to a numpy array
time = np.array(time) # converting the time to a numpy array

# apply the damping to the autocorrelation function
autocorr *= np.exp(-gamma*time**2)

# extend the autocorrelation function to negative times assuming that S(t) = S^*(-t)
time = np.concatenate([-time[::-1], time])  # new time array in range [-t_max, t_max]
autocorr = np.concatenate([np.conjugate(autocorr[::-1]), autocorr])  # new symmetric autocorr in range [-t_max, t_max]

# calculate spectrum from autocorrelation function and the frequency axis corresponding to it
spectrum = np.fft.ifft(autocorr)
freq = 2*np.pi*np.fft.fftfreq(len(time), d=dt)

# plot results
fig, axs = plt.subplots(1, 2, figsize=(8, 3), tight_layout=True)

# autocorrelation function
axs[0].plot(time, np.real(autocorr), label=r'$\mathcal{Re}[S(t)]$')
axs[0].plot(time, np.imag(autocorr), label=r'$\mathcal{Im}[S(t)]$')
axs[0].set_xlabel('Time (a.u.)')
axs[0].set_ylabel(r'$S(t)$')
axs[0].set_title('Autocorrelation Function')
axs[0].legend(frameon=False, labelspacing=0)

# spectrum
axs[1].plot(hbar*freq, np.abs(spectrum))
axs[1].set_xlim(0, np.max(hbar*freq[spectrum > np.max(spectrum)/1000]))
axs[1].set_ylim(0)
axs[1].set_xlabel('Energy (a.u.)')
axs[1].set_ylabel(r'$\mathcal{F}^{-1}[S(t)]$')
axs[1].set_title('Spectrum')

# searching for local maxima of the spectrum
print(f"\nMaxima of the spectrum:")
abs_spectrum = np.abs(spectrum)
loc_max_bool = (abs_spectrum[1:-1] > abs_spectrum[:-2]) & (abs_spectrum[1:-1] > abs_spectrum[2:])
loc_max_index = np.where(loc_max_bool)[0] + 1
loc_max_energies = hbar*freq[loc_max_index]
for index, en in enumerate(loc_max_energies):
    intensity = abs_spectrum[loc_max_index[index]]
    print(f" * State {index}: E = {en:.5f} a.u.; I = {intensity:.5e}")
    axs[1].axvline(en, lw=1, color='black', alpha=0.1)
    axs[1].scatter(en, intensity, marker='x', color='black', s=20)

plt.show()
