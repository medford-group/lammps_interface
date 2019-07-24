from lammps_interface.tools import make_rdf_based_descriptors, make_params_file, gaussian_fit_descriptors
from ase.build import molecule
from ase.io import read
import numpy as np

elements = ['H','O']

atoms = read('traj.traj', index = ':')

etas, rs_s = make_rdf_based_descriptors(atoms, plot = True)

etas = np.append(np.linspace(0.0001, 4, 10), np.logspace(-4,1,10))
rs_s = [0] * 20

make_params_file(elements, etas, rs_s)

#gaussian_fit_descriptors(atoms, n_gaussians = 20, plot = True, nbins = 200)
