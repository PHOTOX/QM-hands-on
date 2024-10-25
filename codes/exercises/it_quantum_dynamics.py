"""Solution for chapter Imaginary-time dynamics and stationary states."""

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -15, 15  # minimum and maximum of x
simtime = 200.0  # simulation time in atomic time units
dt = 0.2  # time step in atomic time units
m = 1.0  # mass in atomic units
nstates = 10  # number of states to be calculated

# physical constants in atomic units
hbar = 1.0

# generate equidistant x grid from xmin to xmax
x = np.linspace(xmin, xmax, ngrid)

# generate momentum grid for the discrete Fourier transform
dx = (xmax - xmin)/(ngrid - 1)
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
psi0 =  # fill in

# create and empty array for final optimized wave functions and energies
psi_optimized = np.zeros(shape=(nstates, len(x)))
energy_optimized = np.zeros(shape=(nstates))

# initiate online plotting to follow the wave function on the fly
plt.ion()
ymin, ymax = np.min(V), np.max(V)  # plotting range of the vertical axis

print("\nImaginary-time quantum dynamics.\n------------------------------------------\n")
# loop through states
for state in range(nstates):
    # creating wave function for propagation
    psi = psi0

    # first, we need to calculate energy of the initial wave function and store it into energy_old variable, so that we
    # can compare at the end of the cycle the new and old energies
    energy_old =  # fill in

    # propagation
    t = 0  # set initial time to 0
    print(f"\nState:  {state + 1:d}")
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

        # remove contributions from previous states, i.e. remove c_i*psi_optimized, where c_i = <psi_optimized|psi>
        for s in range(state):
        # fill in

        # calculate norm
        norm =  # fill in

        # renormalize the wave function so that the norm is one again
        # fill in

        # calculate new norm and check that it is equal to one
        norm =  # fill in
        if np.abs(norm - 1) > 1e-5:
            print(f"ERROR: Renormalization procedure failed, <psi|psi> = {norm:.9f}!")
            exit(1)

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
        # plot current wave function
        plt.fill_between(x, psi**2 + energy, energy, alpha=0.2, color='C1')  # plot |wf|
        plt.text(np.min(x), energy, f"{energy:.4f}")
        # plot previously optimized wave functions
        for s in range(state):
            plt.fill_between(x, psi_optimized[s]**2 + energy_optimized[s], energy_optimized[s], alpha=0.1,
                color='black')
            plt.text(np.min(x), energy_optimized[s], f"{energy_optimized[s]:.4f}")

        # plot potential
        plt.plot(x, V, color='black')
        # plot title and labels
        plt.title(
            f"State: {state + 1}; Time: {t:.2f} a.t.u.\n" + r"$\langle \psi | \psi \rangle$" + f" = {norm:.8f}; $E = $ {energy:.6f} a.u.")
        plt.xlabel("$x$ (a.u.)")
        plt.ylabel("$|\Psi|^2$")

        # set new ymin, ymax
        max_wf = np.max(psi**2)
        ymax = max([ymax, max_wf + energy])
        plt.ylim(ymin, ymax)
        plt.pause(interval=0.0001)  # update plot and wait given interval

        # check energy change between the current and the previous step
        # if |energy - energy_old| < threshold, then the wave function is optimized and break the cycle
        if : # fill in
            break
        else:
            energy_old = energy

    # save the final wave function and energy
    psi_optimized[state] = psi
    energy_optimized[state] = energy

# close online plotting
plt.pause(2.0)  # wait 2s before closing the window
plt.ioff()
plt.close()
