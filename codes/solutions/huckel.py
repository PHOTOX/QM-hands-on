"""Solution for chapter Hückel theory."""

# import necessary libraries
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerical python

# parameters of the molecule
natoms = 4  # number of atoms in the molecule
neighbors = [[1], [0, 2], [1, 3], [2]]  # neighbors of each atom

# parameters of the Hamiltonian matrix
alpha = -11.4  # alpha parameter (eV)
beta = -3.48  # beta parameter (eV)

# generate the Hamiltonian matrix full of zeros
H = np.zeros(shape=(natoms, natoms), dtype=int)

# fill the Hamiltonian matrix with the ones for the neighbors
for atom in range(natoms):
    for neighbor in neighbors[atom]:
        H[atom, neighbor] = 1

# print the Hamiltonian matrix
print(f"Hamiltonian matrix H: \n {np.array2string(H, separator='  ')}\n")

# diagonalize the Hamiltonian matrix
eigenvalues, eigenvectors = np.linalg.eigh(H)
eigenvectors = eigenvectors.T  # transposition because the eigenvectors are stored in columns and we want rows

# NOTE: the eigenvectors are sorted in ascending order, so we need to reverse them
eigenvectors = eigenvectors[::-1,:]

print(np.linalg.eigvalsh(H))

# calculate the energies
energies = alpha - beta*eigenvalues

# print the results for each orbital
for state in range(natoms):
    print(f"* Orbital {state + 1:d}")
    print(f"  Energy = alpha - {eigenvalues[state]:5.2f} * beta = {energies[state]:.2f} eV")
    print("  Psi = " + " + ".join(f"{eigenvectors[state, i]:.4f} * phi_{i + 1:d}" for i in range(natoms)))

# plot orbitals, this works onlu for polyenes
fig, ax = plt.subplots(natoms, 1, figsize=(4.0, 2.0*natoms))
for state in range(natoms):
    eigvec = eigenvectors[state]
    ax[state].plot(range(natoms), np.zeros(natoms), color='black', marker='o', ls='-', lw=1)
    ax[state].scatter(range(0, natoms), np.sign(eigvec)*np.ones(natoms)/2, marker='o', color='blue', s=4000*abs(eigvec),
        alpha=abs(eigvec))
    ax[state].scatter(range(0, natoms), -np.sign(eigvec)*np.ones(natoms)/2, marker='o', color='red', s=4000*abs(eigvec),
        alpha=abs(eigvec))
    ax[state].set_ylabel(f"$\psi_{state + 1:d}$")
    ax[state].set_xticks(range(0, natoms), [f"C$_{i + 1:d}$" for i in range(0, natoms)])
    ax[state].set_yticks([])
    ax[state].set_xlim(-0.5, natoms - 0.5)
    ax[state].set_ylim(-1, 1)
    ax[state].tick_params('both', direction='in', which='both', top=True, right=True)

plt.tight_layout()
plt.show()
