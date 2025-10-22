"""Solution for chapter Bohmian dynamics."""

# todo: note for scripta: examine the case of stationary wf -> there is no velocity field and trajs are stationary

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the simulation
ngrid = 500  # number of grid points
xmin, xmax = -20, 20  # minimum and maximum of x (necessary to set for each simulation)
simtime = 100.0  # simulation time in atomic time units
dt = 0.2  # time step in atomic time units
m = 1.0  # mass in atomic units
ntrajs = 10000  # number of trajectories

# physical constants in atomic units
hbar = 1.0

# todo: I don't need it as a function since I use it just once: put it at the spot of the calculation
# functions
def velocity_field(psi, hbar, m, k):
    # calculate of the gradient of the wave function using FFT
    # first the inverse Fourier transform to momentum space
    psi_k = np.fft.ifft(psi, norm="ortho")
    # apply expT in momentum space
    psi_k *= -1j*k
    # finally the Fourier transform back to coordinate space
    grad_psi = np.fft.fft(psi_k, norm="ortho")  # todo: resolve the issue with the gradient

    # Note: alternative using second order finite differences (shows better performance for edges, compare plots)
    grad_psi = np.gradient(psi, x)

    # calculate current
    v = hbar/m*np.imag(grad_psi/psi)
    return v


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
alpha, x0, p0 = 0.1, 2.5, 0.0
psi = (2*alpha/np.pi)**0.25*np.exp(-alpha*(x - x0)**2 + 1j*p0*(x - x0), dtype=complex)
# todo: remove this testing
psi = (x - x0)*(2*alpha/np.pi)**0.25*np.exp(-alpha*(x - x0)**2 + 1j*p0*(x - x0), dtype=complex)
psi /= np.trapz(y=np.conjugate(psi)*psi, x=x)**0.5  # normalize wave function

density = np.abs(psi)**2  # alternative to np.conjugate(psi)*psi but provides float while the other complex

# sample initial conditions for trajectories
traj = np.zeros((ntrajs), dtype=float)
for t in range(ntrajs):
    max_prob = np.max(density)
    while True:
        x0 = np.random.uniform(xmin, xmax)
        rnd = np.random.uniform(0, max_prob)
        prob = np.interp(x=x0, xp=x, fp=density)
        if rnd < prob:
            traj[t] = x0
            break

# initiate online plotting to follow the wave function on the fly
plt.ion()
ymin, ymax = np.min(V), np.min(V)  # plotting range of the vertical axis

# propagation
t = 0  # set initial time to 0
print("\nLaunching quantum dynamics.\n---------------------------\n")
while t < simtime:  # loop until simulation time is reached
    # first we will propagate the trajectories
    # calculate the velocity field
    v_field = velocity_field(psi, hbar, m, k)

    # propagate the trajectories
    traj = traj + np.interp(x=traj, xp=x, fp=v_field)*dt

    # now propagate the wave function
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

    # print simulation data
    print(f"--Time: {t:.2f} a.t.u.")
    print(f"  <psi|psi> = {norm:.8f}, <E> = {energy:.6f}, <V> = {energyV:.6f}, <T> = {energyT:.6f}")

    # plot
    plt.cla()  # clean figure
    # wave function
    plt.fill_between(x, np.conjugate(psi)*psi, 0, alpha=0.2, color='black', label=r"$|\Psi|$")  # plot |wf|
    plt.plot(x, V - energy, color='black')  # plot potential
    # trajectories
    plt.scatter(traj, np.full_like(traj, 0), color='red', s=1, label="Trajectories")  # plot trajectories
    plt.hist(traj, bins=100, range=(xmin, xmax), color='red', alpha=0.2, density=True,
        label="Trajectory histogram")  # plot histogram of trajectories

    plt.title(f"Time: {t:.2f} a.t.u.\n" + r"$\langle \psi | \psi \rangle$" + f" = {norm:.8f}; $E = $ {energy:.6f} a.u.")
    plt.xlabel("$x$ (a.u.)")
    plt.ylabel("$\Psi$")
    # set new ymin, ymax
    max_wf = np.max(np.abs(psi))
    ymin, ymax = min([ymin, -max_wf]), max([ymax, max_wf])
    plt.ylim(ymin, ymax)
    plt.legend(frameon=False)
    plt.pause(interval=0.0001)  # update plot and wait given interval

# close online plotting
plt.pause(2.0)  # wait 2s before closing the window
plt.ioff()
plt.close()
