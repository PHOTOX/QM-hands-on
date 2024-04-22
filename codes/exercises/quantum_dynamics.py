"""1D quantum dynamics code in real and imaginary time. Exercise for students"""

# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -15, 15  # minimum and maximum of x
simtime = 100.0  # simulation time in atomic time units
dt = 0.5  # time step in atomic time units
m = 1.0  # mass in atomic units
dynamics = ['rt', 'it'][1]  # select if real- or imaginary-time dynamics

# generate x grid
x = np.linspace(xmin, xmax, ngrid)

# generate momentum grid for the discrete Fourier transform
dx = (xmax - xmin)/(ngrid - 1)
p = 2*np.pi*np.fft.fftfreq(ngrid, d=dx)

# generate potential V
V = 0.005*x**2

# generate kinetic energy T in momentum space
T = p**2/2/m

# generate V propagator with time step dt/2
if dynamics == 'rt':
    expV =
elif dynamics == 'it':
    expV =

# generate T propagator in momentum space with time step dt
if dynamics == 'rt':
    expT =
elif dynamics == 'it':
    expT =

# initiate wave function; be careful, it must be a complex numpy array
psi =

# initiate online plotting to follow the wave function on the fly
plt.ion()

# propagation
t = 0  # set initial time to 0
print("\nLaunching quantum dynamics.\n---------------------------\n")
while t < simtime:  # loop until simulation time is reached
    # propagate half-step in V


    # propagate full step in T
    # first Fourier transform to momentum space
    psi_k = np.fft.fft(psi, norm="ortho")
    # apply expT in momentum space

    # finally inverse Fourier transform back to coordinate space
    psi = np.fft.ifft(psi_k, norm="ortho")

    # propagate half-step in V for the second time


    # calculate new time after propagation
    t += dt

    # if imaginary-time dynamics, renormalization is necessary
    if dynamics == 'it':


    # calculate norm
    norm =

    # calculate expectation value of energy
    # potential energy <V>
    energyV =
    # kinetic energy <T>, T operator will be again applied in the momentum space
    energyT =
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
