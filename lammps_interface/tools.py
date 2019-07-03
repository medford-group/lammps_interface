# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 12:31:50 2019

@author: benjamin
"""
import importlib

import numpy as np
import os
from ase import io
from pymatgen.io.ase import AseAtomsAdaptor as adaptor
from ase.spacegroup import crystal
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from shutil import copyfile

force_fields = {'tip3p_2004':'TIP3P_2004','tip3p_1983_charmm_hybrid':'tip3p_1983_charmm_hybrid'}

def make_box_of_molecules(molecules, num_molecules, box,
                          tolerance=2.0, outputfile='box.pdb',
                          radius = None):
    """
    This function makes a box of molecules suitable for use in
    classical simulations using the packmol package. The default settings
    assume you want to make a box of water
    
    inputs:
        molecules (list of ASE atoms objects):
            a list of molecules you want in the box
        num_molecules (list):
            a list of the number of each molecule desired, the numnbers
            should be in the same order as the molecules arguement.
        box (3x3 array or empty atoms object with specified cell):
            This is a way to specify the size of the box
        tolerance (float):
            the closest two molecules are allows to be to one another
        outputfile (str):
            a file to place the output in. must end in either .pdb or 
            .xyz, the two formats packmol supports
        radius (list):
            The radius of various molecules, used to ensure molecules do
            not overlap
    
    returns:
        atoms (ASE atoms object):
            an ASE atoms object which has the box created by packmol
    """
    
    if hasattr(box, 'cell'):
        box = box.cell
    if outputfile.split('.')[-1] not in ['pdb','xyz']:
        raise ValueError('file must end with .pdb or .xyz')
    if type(molecules) != list:
        molecules = [molecules]
    if type(num_molecules) != list:
        num_molecules = [num_molecules]
     
    f = open('pk.inp','w')
    f.write('tolerance ' + str(tolerance) + '\n')
    f.write('output ' + outputfile+ '\n')
    f.write('filetype ' + outputfile.split('.')[-1] + '\n')
    
    for i, molecule, num in zip(range(len(molecules)), molecules, num_molecules):
        filename = molecule.get_chemical_formula() + '_' + str(i) + '.pdb'
        molecule.write(filename)
        f.write('structure ' + filename + '\n')
        f.write('  number ' + str(num) + '\n')
        l = [np.linalg.norm(a) for a in box]
        f.write('  inside box 0. 0. 0. '  # Write the box size
                + str(l[0]) + ' ' 
                + str(l[1]) + ' '
                + str(l[2]) + '\n')
        if radius is not None:
            f.write('   radius ' + str(radius) + '\n')
        f.write('end structure' + '\n\n')
    f.close()
    os.system('packmol < pk.inp > pk.log')
    atoms = io.read(outputfile)
    atoms.cell = box    
    os.system('rm pk* *.pdb')
    return atoms

def write_lammps_inputs_moltemplate(atoms, force_field, num_molecules):
    """
    A function to wrap the moltemplate package, to make it easier to call
    from python. This writes the input files for lammps under the following
    names: mt.in, mt.in.init, mt.in, settings, mt.data. The output of 
    moltemplate is saved to mt.log.

    inputs:
        atoms (ASE atoms object):
            The MD simulation box, complete with boundary conditions and
            cell size
        force_field (str):
            The name of the force field you'd like to use. These are stored
            in moltemplate and may not be 100% accurate, make sure you check.
        num_molecules (int):
            The number of molecules in the atoms object. Currently this only
            supports a single molecule type

    returns:
        None
    """
    if force_field not in force_fields:
        raise Exception('the force field you selected, {}, is not listed as being included in moltemplate'.format(force_field))
    
    atoms.write('mt.pdb')
    f = open('mt.lt','w')
    f.write('import "' + force_field + '.lt"'+ '\n\n')
    f.write('write_once("Data Boundary") {\n')
    f.write('0.0 ' + str(np.linalg.norm(atoms.cell[0])) + ' xlo xhi\n')
    f.write('0.0 ' + str(np.linalg.norm(atoms.cell[1])) + ' ylo yhi\n')
    f.write('0.0 ' + str(np.linalg.norm(atoms.cell[2])) + ' zlo zhi\n')
    f.write('}\n\n')
    
    f.write('atm = new '+ force_field.upper() +
            ' [' + str(num_molecules) + ']')
    f.close()
    os.system('moltemplate.sh -pdb mt.pdb -atomstyle full mt.lt &> mt.log')

def write_lammps_data(atoms, filename = 'lmp.data', 
                      atoms_style = 'full', bonding = False):
    """
    a function for writing LAMMPS data files. These files are the "atomic positions"
    input file. Currently only the "full" atom type is supported. This is probably all
    you need.

    inputs:
        atoms:
            The atoms object
        filename:
            whatever you'd like the filename to be
        atoms_style:
            the atoms style to be used in LAMMPS
        bonding:
            if you want bonding included

    returns:
        None
    """
    from ase.data import atomic_masses, chemical_symbols

    supported_styles = ['full']

    # check inputs
    if atoms_style not in supported_styles:
        raise Exception('atoms style {} is not supported currently, only syles: {} are supported'.format(atoms_style, supported_styles))

    if bonding == True:
        raise Exception('This function has not implemented bonding')

    # sort out the numbering of the atomic symbols
    number_dict = {}
    numbers = list(set(atoms.get_atomic_numbers()))
    # order lowest to highest
    numbers.sort()
    for i, number in enumerate(numbers):
        number_dict[number] = i + 1 # maps atomic numbers to LAMMPS atom type numbering
    
    # sort out charges
    charges = atoms.get_initial_charges()

    f = open(filename, 'w')

    # header stuff
    f.write('Made with Medford Group LAMMPS Interface\n\n\n')
    f.write('{}  atoms\n'.format(len(atoms)))
    f.write('{}  atom types\n'.format(len(number_dict))) 
    f.write('0 {} xlo xhi\n'.format(atoms.cell[0][0]))
    f.write('0 {} ylo yhi\n'.format(atoms.cell[1][1]))
    f.write('0 {} zlo zhi\n\n'.format(atoms.cell[2][2]))

    # atomic masses
    f.write('Masses\n\n')
    for atomic_number, atom_type_number in number_dict.items():
        f.write('{} {} # {}\n'.format(atom_type_number,
                               atomic_masses[atomic_number],
                               chemical_symbols[atomic_number]))
    
    # The actual atomic inputs
    f.write('\n\n')
    f.write('Atoms\n\n')

    for i, atom in enumerate(atoms):
        x, y, z = atom.position
        f.write('{} {} {} {} {} {} {}\n'.format(
                                            str(i + 1), # index
                                            '0', # bonding index
                                            number_dict[atom.number], # atom type numbers
                                            atom.charge, # charges
                                            x, y, z))
    f.close() 

    
def fix_xyz_files(fil, data_file):
    """
    replaces the atom numbers with the actual atom names in the xyz file
    LAMMPS dumps out.

    inputs:
        fil (str):
            The xyz file you'd like to fix

    returns:
        data_file (str):
            The LAMMPS input data file
    """
    element_key = elements_from_datafile(data_file)
    f = open(fil,'r')
    xyz = f.read()
    for key, value in element_key.items():
        xyz = xyz.replace('\n{} '.format(key),'\n{} '.format(value))
    f.close()
    f = open(fil,'w')
    f.write(xyz)
    f.close()

def parse_custom_dump(dump, datafile, label = 'atoms',
                      energyfile = None, write_traj = False):
    """
    This function parses the output of a LAMMPS custom dump file. Currently
    this function assumes you've put in this for the custom dump: 
    "element x y z fx fy fz". There are plans to expand this to be more general
    if the need arises. This assumes units are set to real.

    inputs:
        dump (str):
            The filename of the dump from LAMMPS
        datafile (str):
            The atomic positions/bonds datafile from lammps
        
    returns:
        atoms_list (list):
            A list contianing atoms objects for each recorded timestep.
    """
    from ase.atoms import Atoms    
    f = open(dump,'r')
    text = f.read()
    if 'type x y z fx fy fz' not in text:
        raise Exception('the dump file is not in the correct format. Please make the dump "type x y z fx fy fz". Further functionality is planned to be added.')
    if datafile is not None:
        element_key = elements_from_datafile(datafile)
    steps = text.split('ITEM: TIMESTEP')[1:]
    atoms_list = []
    if energyfile is not None:
        g = open('energy.txt')
        pe = g.readlines()[1:]
    for j, step in enumerate(steps):
        step = step.split('ITEM') 
        # parse the body
        dump_format = step[3].split('\n')[0].split('ATOMS')[1].strip()
        dump_format = dump_format.split()
        data = np.empty([len(step[3].split('\n')[1:-1]),len(dump_format)])
        for i, line in enumerate(step[3].split('\n')[1:-1]):
            data[i] = [float(a) for a in line.split()]
        # parse the elements
        if 'type' in dump_format:
            typ = dump_format.index('type')
            elements = [element_key[a] for a in data[:,typ]]
            atoms = Atoms(elements)
            del elements
        else:
            raise Exception('The atom type must be in the dump file to parse it into an atoms object')
        # parse the atomic positions
        if 'x' in dump_format and 'y' in dump_format and 'z' in dump_format:
            x = dump_format.index('x')
            y = dump_format.index('y')
            z = dump_format.index('z')
            pos = data[:,[x, y, z]]
            atoms.positions = pos
            del pos
        # parse the forces
        if 'fx' in dump_format and 'fy' in dump_format and 'fz' in dump_format:
            from ase.units import kcal, mol
            fx = dump_format.index('fx')
            fy = dump_format.index('fy')
            fz = dump_format.index('fz')
            forces = data[:, [fx, fy, fz]]
            forces *= (kcal / mol)  # convert to eV/A
        # parse the potential energy
        if 'c_energy' in dump_format:
            eng = dump_format.index('c_energy')
            per_atom_energy = data[:,eng] * (kcal / mol)  # convert to eV/A
            atoms.set_initial_charges(per_atom_energy)
            energy = sum(per_atom_energy)
        # recover the cell
        cell = step[2].strip().split('\n')
        if cell[0].split()[-3:] == ['pp', 'pp', 'pp']:
            atoms.set_pbc([True, True, True])
        cell_mat = np.zeros([3,3])
        for i, direction in enumerate(cell[-3:]):
            bot, top = [float(a) for a in direction.split()]
            cell_mat[i][i] = top - bot
            atoms.positions[:,i] -= bot
        atoms.cell = cell_mat
        # build calculator, it is critical that this is done last
        from ase.calculators.singlepoint import SinglePointCalculator as SP
        calc = SP(atoms, forces = forces, energy = energy)
        atoms.set_calculator(calc)
        atoms_list.append(atoms)
        del data
    if len(atoms_list) == 1:
        return atoms_list[0]
    if write_traj == True:
        from ase.io.trajectory import TrajectoryWriter 
        tw = TrajectoryWriter('{}.traj'.format(label))
        for atoms in atoms_list:
            tw.write(atoms)
    return atoms_list

def elements_from_datafile(data_file):
    """
    This function reads a LAMMPS datafile and returns a dictionary
    mapping the atom type numbers from the file to the chemical symbols
    based on the masses.
    
    inputs:
        data_file (str):
            The path of the LAMMPS datafile

    returns:
        element_key (dict):
            A dictionary mapping the atom type numbers to their chemical
            symbols
    """
    from ase.data import atomic_masses_iupac2016, chemical_symbols

    atomic_masses_iupac2016 = [np.round(a,decimals = 1) for a in atomic_masses_iupac2016]
    atomic_masses_iupac2016[0] = 0  # just makes everything easier
    f = open(data_file,'r')
    data = f.read()
    f.close()
    data = data.split('Masses')[1]
    data = data.split('Atoms')[0].strip()
    data = data.split('\n')
    element_key = {}
    for atom in data:
        atm = atom.strip().split()
        mass = float(atm[1])
        atm_num = atomic_masses_iupac2016.index(np.round(mass,decimals=1))
        element_key[int(atm[0])] = chemical_symbols[atm_num]
    return element_key

def strip_bonding_information(filename):
    """
    function to remove the bonding information from a lammps datafile
    to be used as a reaxff input
    """

    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    
    for i, line in enumerate(lines):
        if 'bond' in line or 'angle' in line:
            lines[i] = ''
        if 'Bond' in line:
            j = i
    if 'j' in locals():
        del lines[j:]
    f = open(filename,'w')
    f.writelines(lines)
    f.close()
    
def make_standard_input(calculation_type = 'reaxff',
                        ff_file = 'TiO2_water.reax',
                        timestep = 0.2, 
                        ensemble = 'nvt',
                        steps = 20000,
                        temp = 300):
    """
    function to just generate a generic input file to run a basic
    simulation. No frills, it just works.
    """
    from .standard_inputs import input_files

    if calculation_type not in input_files.keys():
        raise Exception('the calculation type you selected, {}, is not implemented yet'.format(calculation_type))
     
    from .standard_inputs import input_files

    with open('lmp.input','w') as f:
        f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     timestep, 
                                                     temp, 
                                                     steps))
    if calculation_type == 'reaxff':
        import inspect
        cwd = os.getcwd()
        package_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        ff = os.path.join(package_directory, 'special_input_files/' + ff_file)
        control = os.path.join(package_directory, 'special_input_files/control.reaxff')
        copyfile(ff, cwd)
        copyfile(control, cwd)

    
    

def surround_with_molecules(atoms, particle_spacing = 8, 
                            fluid_molecule = 'H2O', molar_density = 55.5556,
                            metal = 'Ti', shape = 'spherical',
                            shape_factor = 0.9):
    """
    takes in an object (usually a nano-particle) and surrounds it with water 
    molecules at approximately the right for room temperature density. It assumes 
    your object is roughly spherical.

    inputs:
        atoms:
            the atoms object you'd like to surround with water
        particle_spacing:
            how much spacing you'd like around the particle
        fluid_molecule:
            The molecule you want to cover your slab with, can be either an atoms object
            or a string (must be in the g2 set if it is a string)

        metal:
            If you have a particular atom you'd like centered, in the middle
            of the cell, input it here
        shape:
            the approximate shape of your particle, right now only spherical
            and rectangular are implemented
        shape_factor:
            the factor (between 0 and 1) by which you think the spherical/
            rectangular volume is off. This is multiplied by the calculated
            volume

    returns:
        atoms:
            a periodic atoms object surrounded with water

    """
    from ase.build import molecule
    from ase.atoms import Atoms
    from ase.units import mol, m
    if type(fluid_molecule) == str:
        fluid_molecule = molecule(fluid_molecule)

    # Figure out what molecule is provided
    if type(fluid_molecule) == str:
        fluid_molecule = molecule(fluid_molecule)
    # Calculate number density
    number_density = molar_density * 1000 * mol / m ** 3
    d = []
    for dimension in range(3):
        smallest_value = min(atoms.positions[:, dimension])
        largest_value = max(atoms.positions[:,dimension])
        max_dist = largest_value - smallest_value
        d.append(max_dist)
    atoms.set_cell([a + particle_spacing for a in d])
    atoms.center()
    if shape == 'spherical':
        # 4/3pi*r**3
        volume = 4 / 3 * np.pi * (max(d) / 2) ** 3
    elif shape == 'rectangular':
        volume = np.product(d)
    volume *= shape_factor

    number_of_molecules = int(np.floor((atoms.get_volume() - (d) ** 3) * 0.0333679)) # approximate density of water

    # this next part is such a mess, I'm so sorry
    # all this is doing is trying to make the particle more centered in the unit cell
    atoms = make_box_of_molecules([fluid_molecule,atoms], [number_of_molecules,1], atoms.cell)
    metal_s = [a for a in atoms if a.symbol == metal]
    metal_object = Atoms(cell = atoms.cell)
    for atom in metal_s:
        metal_object += atom
    initial_metal_positions = metal_object.positions.copy()
    metal_object.center()
    shift = initial_metal_positions[0] - metal_object.positions[0]
    atoms.positions -= shift
    # </mess>
    atoms.wrap(pbc = [True] * 3)
    return atoms

def make_wulffish_nanoparticle(atoms, millers, surface_energies,
                               rmax, prune = False):
    """
    wraps the nanoparticle generation functionality of MPInterface (link below).
    It's not perfect for complicated structures, but it gets you pretty close.

    inputs:
        atoms:
            the bulk structure you want to make a nanoparticle of
        millers:
            a list of miller indicies you want to include i.e. [(1,1,1),(1,0,0)]
        surface_energies:
            a list of energies that corresponds to the miller indicies

    returns:
        particle:
            the atoms object of the wulff construction


    https://github.com/henniggroup/MPInterfaces
    """
    from mpinterfaces.nanoparticle import Nanoparticle
    structure = adaptor.get_structure(atoms)

    sa = SpacegroupAnalyzer(structure)
    structure_conventional = sa.get_conventional_standard_structure()

    nanoparticle = Nanoparticle(structure_conventional, rmax=rmax,
                            hkl_family=millers,
                            surface_energies=surface_energies)

    nanoparticle.create()
    #particle = adaptor.get_atoms(nanoparticle)
    particle = adaptor.get_atoms(nanoparticle.get_boxed_structure(10**6,10**6,10**6))
    particle.set_cell([0] * 3)
    if prune == True:
        particle = prune_oxygens(particle, metal = 'Ti')
    return particle

def prune_oxygens(atoms, metal = None):
    """
    a function that removes all oxygens that have a coordination less
    that 2 based on a 2.2 A maximum bond distance. This is meant for
    nanoparticles

    inputs:
        atoms:
            the atoms object of the nanoparticle
        metal:
            If not set to none, a selected metal atom will also be 
            checked

    returns:
        atoms:
            the new (pruned) atoms object
    """
    inds = []
    for i in range(len(atoms)):
        if atoms[i].symbol == 'O':
            bond_distance = 2.2
            minumum_coordination = 2
        elif metal is not None:
            if atoms[i].symbol == metal:
                bond_distance = 3
                minumum_coordination = 4 # arbitrary
        search = list(range(len(atoms)))
        search.remove(i)
        distances = atoms.get_distances(i, search)
        mask = distances < bond_distance
        cn = list(mask).count(True)
        if cn < minumum_coordination:
            inds.append(i)
    for index in inds[::-1]:
        del atoms[index]
    return atoms



def put_molecules_on_slab(atoms, offset = 1.5, 
                         molar_density = 55.55556,
                         fluid_molecule = 'H2O'):
    """
    overlays a water layer on top of a slab. This does not care about what the physics
    of water on this slab should be, it just puts some waters above it from packmol. The
    cell of the slab should be a orthogonal

    inputs:
        atoms:
            the atoms object of the slab
        offset:
            how far the water layer should be offset
        molar_density:
            the molar density (mol/L) of the water you'd like to put on the slab
        fluid_molecule:
            The molecule you want to cover your slab with, can be either an atoms object
            or a string (must be in the g2 set if it is a string)

    returns:
        slab:
            the water covered slab
    """
    from ase.build import molecule
    from ase.units import mol, m
    
    # Figure out what molecule is provided
    if type(fluid_molecule) == str:
        fluid_molecule = molecule(fluid_molecule)
    # Calculate number density
    number_density = molar_density * 1000 * mol / m ** 3
    # Find dimensions of the box
    highest_atom = max(atoms.positions[:,2])
    length, width, top_of_cell = atoms.cell[0,0], atoms.cell[1,1],atoms.cell[2,2]
    height = top_of_cell - highest_atom
    volume = length * width * height
    number_of_molecules = int(np.ceil(volume * number_density)) # rounding up by default
    # Make the box
    water_layer = make_box_of_molecules([fluid_molecule], [number_of_molecules],
                                        box = [length, width, height - offset * 2])
    water_layer.positions += np.array([0, 0, highest_atom + offset])
    slab = atoms + water_layer
    return slab


def make_rdf_based_descriptors(images, n_descriptors = 20,
                               index = 0, cutoff = 6.5, 
                               fall_off_percent = 0.01,
                               descriptor_type = 'simplenn',
                               plot = False):
    """
    generates reasonable values for eta and rs for a given image in a trajectory
    based on the radial distribution function.
    """
    from scipy.integrate import trapz
    from ase.ga.utilities import get_rdf
    #fall_off_percent = 0.2 # arbitrarily chosen
    #localization_distance = 0.2 # arbitrarily chosen

    if type(images) == list:
        atoms = images[index]
    else:
        atoms = images

    rdf = get_rdf(atoms, rmax = cutoff, nbins = 200)
    rdf_1, distances = rdf
    # find 75% index
    cut = int(-1 * np.ceil(len(rdf) * 0.75))
    rdf_mean = np.mean(rdf_1[cut:])
    # subtract off the mean and make all positive
    #rdf_1 -= rdf_mean
    rdf_1 = abs(rdf_1)
    integrals = []
    for i in range(len(rdf_1)):
        integrals.append(trapz(rdf_1[:i+1]))
    # divide the integral evenly
    increment = integrals[-1] / (n_descriptors + 1)
    n = 1
    descriptor_distances = []
    for distance, integral in zip(distances, integrals):
        if integral > (n * increment):
            descriptor_distances.append(distance)
            n += 1
            if len(descriptor_distances) == n_descriptors:
                break
    etas = []
    rs_s = [0] * n_descriptors
    for distance in descriptor_distances:
        etas.append(-1 * np.log(0.1) / distance ** 2)
    for i, distance in enumerate(descriptor_distances):
        if i == 0:
            localization_distance = abs(descriptor_distances[i+1] - distance)
        elif i == len(descriptor_distances) - 1:
            localization_distance = (abs(descriptor_distances[i-1] -\
                                        distance) + abs(cutoff - distance)) / 2
        else:
            localization_distance = (abs(descriptor_distances[i+1] - distance) +\
                                    abs(descriptor_distances[i-1] - distance)) / 2
            #localization_distance = abs(descriptor_distances[i+1] - distance)
        etas.append(-1 * np.log(fall_off_percent) / localization_distance ** 2)
        #rs_s.append(distance - 0.1)
        rs_s.append(distance)
    if plot == True:
        import matplotlib.pyplot as plt
        plt.plot(distances, rdf_1)
        for rs, eta in zip(rs_s, etas):
            x = np.linspace(0,cutoff,1000)
            y = [np.exp(-1 * eta * (a - rs) ** 2) for a in x]
            plt.plot(x,y)
        plt.show()
    return etas, rs_s

def make_params_file(files, etas, rs_s, dist_dict,cutoff = 6.0):
    """
    make a params file for simple_NN
    """
    for file in files:
        with open('params_{}'.format(file),'w') as f:
            # G2
            for species in range(1,len(files)+1):
                for eta, Rs in zip(etas, rs_s):
                    f.write('2 {} 0 {} {} {} 0.0\n'.format(species, cutoff,
                                                           np.round(eta, 6), Rs))
            # G4 
            for i in range(1,len(files)+1):
                n = i
                while True:
                    for eta in np.logspace(-4,-2,4):
                        for lamda in [1.0,-1.0]:
                            for zeta in [1.0,4.0]:
                                f.write('4 {} {} {} {} {} {}\n'.format(i, n, cutoff,
                                                                         np.round(eta, 6),
                                                                         zeta, lamda))
                    n += 1
                    if n > len(files):
                        break


def convert_to_csv_file(atoms, filename = 'atoms.csv'):
    """
    takes an atoms object and converts it to a csv format: 'x, y, z, atom'

    inputs:
        atoms:
            the atoms object
        filename:
            the name of the file you want to write to

    returns:
        None
    """
    with open(filename,'w') as f:
        for atom in atoms:
            x, y, z = atom.position
            f.write('{}, {}, {}, {}\n'.format(x,y,z,atom.symbol))

