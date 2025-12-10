"""Solution for chapter Quantum dynamics: real-time quantum dynamics."""

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -1, 1  # minimum and maximum of x
simtime = 100.0  # simulation time in atomic time units
dt = 0.01  # time step in atomic time units
m = 20.0  # mass in atomic units

# physical constants in atomic units`````
hbar = 1.0

# generate equidistant x grid from xmin to xmax
x = np.linspace(xmin, xmax, ngrid)

# generate momentum grid for the discrete Fourier transform
dx = (xmax - xmin)/(ngrid - 1)
k = 2*np.pi*np.fft.fftfreq(ngrid, d=dx)

# generate potential V
omega = 2
V = 0.5*m*omega**2*x**2

# generate kinetic energy T in momentum space
T = hbar**2*k**2/2/m

# generate V propagator with time step dt/2
expV = np.exp(-1j*V*dt/2)

# generate T propagator in momentum space with time step dt
expT = np.exp(-1j*T*dt)

# initiate wave function; be careful, it must be a complex numpy array
alpha, x0, p0 = m*omega/2, 0.0, 20.0
psi = (2*alpha/np.pi)**0.25*np.exp(-alpha*(x - x0)**2 + 1j*p0*(x - x0), dtype=complex)

# initiate online plotting to follow the wave function on the fly
plt.ion()

# propagation
t = 0  # set initial time to 0
print("\nLaunching quantum dynamics.\n---------------------------\n")
while t < simtime:  # loop until simulation time is reached
    # propagate half-step in V
    psi *= expV

    # propagate full step in T
    # first Fourier transform to momentum space
    psi_k = np.fft.fft(psi, norm="ortho")
    # apply expT in momentum space
    psi_k *= expT
    # finally inverse Fourier transform back to coordinate space
    psi = np.fft.ifft(psi_k, norm="ortho")

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
    psi_k = np.fft.fft(psi, norm="ortho")
    psi_t = np.fft.ifft(T*psi_k, norm="ortho")
    energyT = np.real(np.trapz(y=np.conjugate(psi)*psi_t, x=x))/norm
    # total energy <E>
    energy = energyV + energyT

    # print simulation data
    print(f"--Time: {t:.2f} a.t.u.")
    print(f"  <psi|psi> = {norm:.4f}, <E> = {energy:.4f}, <V> = {energyV:.4f}, <T> = {energyT:.4f}")

    # plot
    plt.cla()  # clean figure
    plt.plot(x, np.conjugate(psi)*psi + energy, color='grey', label=r"$|\Psi|^2$")  # density
    plt.fill_between(x, np.real(psi) + energy, np.zeros(ngrid) + energy, alpha=0.2, label=r"$Re[\Psi]$")  # plot Re(wf)
    plt.fill_between(x, np.imag(psi) + energy, np.zeros(ngrid) + energy, alpha=0.2, label=r"$Im[\Psi]$")  # plot Im(wf)
    plt.plot(x, V, color='black')  # plot potential
    plt.title(f"$E = $ {energy:.4f} a.u.")
    plt.xlabel("$x$ (a.u.)")
    plt.ylabel("$\Psi$")
    plt.legend(frameon=False)
    plt.pause(interval=0.01)  # update plot and wait given interval

# close online plotting
plt.pause(2.0)  # wait 2s before closing the window
plt.ioff()
