import h5py
import numpy as np

# Open the HDF5 file
h5file = h5py.File("seperate_mscapst2_1.rasscf.h5", "r")

# Explore the groups
def print_groups(name, obj):
    print(name, type(obj))
    
h5file.visititems(print_groups)

# Usually the effective Hamiltonian is under: /CASPT2/Hamiltonian
# or /CAS/Hamiltonian
if "/CASPT2/Hamiltonian" in h5file:
    H_eff = np.array(h5file["/CASPT2/Hamiltonian"])
    print("Effective Hamiltonian (full precision):")
    print(H_eff)
else:
    print("Check the groups in the HDF5 file; look for Hamiltonian or coupling matrices.")

h5file.close()

for i in range(H_eff.shape[0]):
    for j in range(H_eff.shape[1]):
        print(f"{H_eff[i,j]:.20e}", end=" ")
    print()
