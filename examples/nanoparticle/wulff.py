from pymatgen.analysis.wulff import WulffShape
from pymatgen.io.ase import AseAtomsAdaptor as ad
from ase.spacegroup import crystal
import numpy as np
from lammps_interface.convience_tools import make_wulffish_nanoparticle, surround_with_water 


# some parameters
rmax = 20 # max radius
surface_energies = [0.6, 0.47, 0.95] # surface energies
hkl_family = [(1,0,0),(1,1,0),(0,1,1)] # corresponding surfaces
# make a rutile lattice
a = 4.6
c = 2.95
rutile =crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

rutile = make_wulffish_nanoparticle(rutile)

wet_rutile = surround_with_water(rutile)

wet_rutile.write('nanoparticle.traj')

