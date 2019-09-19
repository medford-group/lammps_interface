from lammps_interface.tools import single_point_lammps
from ase.build import molecule

atoms = molecule('H2O')
atoms.set_cell([10,10,10])
atoms.set_pbc([True] * 3)

atoms = single_point_lammps(atoms, method='simple_nn_single_point',
                            ff_file='ffield.reax.water_2017')

print(atoms.get_potential_energy())
