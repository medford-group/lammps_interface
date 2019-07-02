from pymatgen.analysis.wulff import WulffShape
from pymatgen.io.ase import AseAtomsAdaptor as ad
from ase.spacegroup import crystal
import numpy as np
from lammps_interface.tools import make_wulffish_nanoparticle, surround_with_water 
from ase.visualize import view
from ase.neighborlist import neighbor_list
from ase.geometry.analysis import Analysis
from pymatgen.analysis.local_env import CrystalNN, NearNeighbors


# some parameters
rmax = 20 # max radius
surface_energies = [0.6, 0.47, 0.95] # surface energies
hkl_family = [(1,0,0),(1,1,0),(0,1,1)] # corresponding surfaces
# make a rutile lattice
a = 4.6
c = 2.95
rutile =crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

atoms = make_wulffish_nanoparticle(rutile, millers = hkl_family,
                                    surface_energies = surface_energies,
                                    rmax = 24, prune = True)



view(atoms)
#wet_rutile = surround_with_water(atoms)

#wet_rutile.write('nanoparticle.traj')

