from ase.build import bulk, surface
from lammps_interface.tools import put_molecules_on_slab, make_standard_input, write_lammps_data
from ase.visualize import view

iron_surface = surface('Fe', (1,0,0), layers = 4, vacuum = 10)
iron_surface *= (2,2,1)

hydrated_iron = put_molecules_on_slab(iron_surface, 
                                      fluid_molecule = 'H2O',
                                      molar_density = 55.5556,
                                      offset = 1)
view(hydrated_iron)

make_standard_input(calculation_type = 'reaxff')

write_lammps_data(hydrated_iron, filename = 'lmp.data', bonding = False)
