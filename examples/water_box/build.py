from ase.build import molecule
from tools import make_box_of_molecules, write_lammps_inputs_moltemplate
import os

mol = molecule('H2O')

# make a box of water
atoms = make_box_of_molecules(mol, 700, [[40,0,0],[0,40,0],[0,0,40]])

# make LAMMPS input files
write_lammps_inputs_moltemplate(atoms, 'tip3p_2004', 700)

# run lammps
os.system('lmp < mt.in > lmp.log')
