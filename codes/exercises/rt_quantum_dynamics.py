"""Exercise for chapter Quantum dynamics in real time."""

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -15, 15  # minimum and maximum of x
simtime = 100.0  # simulation time in atomic time units
dt = 0.2  # time step in atomic time units
m = 1.0  # mass in atomic units

# physical constants in atomic units
hbar =  # fill in

# generate equidistant x grid from xmin to xmax
x =  # fill in (we recommend using function linspace from numpy)

# generate momentum grid for the discrete Fourier transform
dx =  # fill in distance between two neighbouring (x_{i+1}-x_{i})
k = 2*np.pi*np.fft.fftfreq(ngrid, d=dx)

# generate potential V
V =  # fill in

# generate kinetic energy T in momentum space
T =  # fill in

# generate V propagator with time step dt/2
expV =  # fill in

# generate T propagator in momentum space with time step dt
expT =  # fill in

# initiate wave function; be careful, it must be a complex numpy array
psi =  # fill in

# initiate online plotting to follow the wave function on the fly
plt.ion()
ymin, ymax = np.min(V), np.max(V) # plotting range of the vertical axis

# propagation
t = 0  # set initial time to 0
print("\nLaunching quantum dynamics.\n---------------------------\n")
while t < simtime:  # loop until simulation time is reached
    # propagate half-step in V
    # fill in

    # propagate full step in T
    # first the inverse Fourier transform to momentum space
    psi_k = np.fft.ifft(psi, norm="ortho")
    # apply expT in momentum space
    # fill in
    # finally the Fourier transform back to coordinate space
    psi = np.fft.fft(psi_k, norm="ortho")

    # propagate half-step in V for the second time
    # fill in

    # calculate new time after propagation
    t += dt

    # calculate norm
    norm =  # fill in

    # check that norm is conserved
    # fill in, choose a reasonable threshold, e.g. 0.00001

    # calculate expectation value of energy
    # potential energy <V>
    energyV =  # fill in
    # kinetic energy <T>, T operator will be again applied in the momentum space
    energyT =  # fill in
    # total energy <E>
    energy =  # fill in

    # print simulation data
    print(f"--Time: {t:.2f} a.t.u.")
    print(f"  <psi|psi> = {norm:.8f}, <E> = {energy:.6f}, <V> = {energyV:.6f}, <T> = {energyT:.6f}")

    # plot
    plt.cla()  # clean figure
    plt.fill_between(x, np.abs(psi) + energy, energy, alpha=0.2, color='black', label=r"$|\Psi|$")  # plot |wf|
    plt.plot(x, np.real(psi) + energy, linestyle='--', label=r"$Re[\Psi]$")  # plot Re(wf)
    plt.plot(x, np.imag(psi) + energy, linestyle='--', label=r"$Im[\Psi]$")  # plot Im(wf)
    plt.plot(x, V, color='black')  # plot potential
    plt.title(f"Time: {t:.2f} a.t.u.\n" + r"$\langle \psi | \psi \rangle$" + f" = {norm:.8f}; $E = $ {energy:.6f} a.u.")
    plt.xlabel("$x$ (a.u.)")
    plt.ylabel("$\Psi$")
    # set new ymin, ymax
    max_wf = np.max(np.abs(psi))
    ymin, ymax = min([ymin, -max_wf+energy]), max([ymax, max_wf+energy])
    plt.ylim(ymin, ymax)
    plt.legend(frameon=False)
    plt.pause(interval=0.0001)  # update plot and wait given interval

# close online plotting
plt.pause(2.0)  # wait 2s before closing the window
plt.ioff()
plt.close()