from pymatgen.io.ase import AseAtomsAdaptor as adaptor
from ase.spacegroup import crystal
import numpy as np
from lammps_interface.convience_tools import make_wulffish_nanoparticle, surround_with_water,put_water_on_slab
from ase.visualize import view
from pymatgen.core.surface import generate_all_slabs, SlabGenerator


# some parameters
rmax = 20 # max radius
surface_energies = [0.6, 0.47, 0.95] # surface energies
hkl_family = [(1,0,0),(1,1,0),(0,1,1)] # corresponding surfaces
# make a rutile lattice
a = 4.6
c = 2.95
rutile =crystal(['Ti', 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                spacegroup=136, cellpar=[a, a, c, 90, 90, 90])

structure = adaptor.get_structure(rutile)

gen = SlabGenerator(structure, miller_index = (1,1,0), lll_reduce = True,
                    min_slab_size = 10, min_vacuum_size = 20)

slab = gen.get_slabs()

#slabs = generate_all_slabs(structure, max_index = 1,
#                           min_slab_size = 10, min_vacuum_size = 10)

#print(slabs[3].miller_index)
slab = adaptor.get_atoms(slab[1].get_orthogonal_c_slab()) * (3,3,1)
slab.center()
slab = put_water_on_slab(slab)
view(slab)

