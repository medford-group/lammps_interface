# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 12:31:50 2019
author: Ben Comer (Georgia Tech)
"""
import importlib
from collections import namedtuple, defaultdict
import itertools
import shutil

import numpy as np
import os
from ase import io
import json
from io import StringIO
import pickle

from ase.geometry.analysis import Analysis
from pymatgen.io.ase import AseAtomsAdaptor as adaptor
from ase.spacegroup import crystal
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from shutil import copyfile

force_fields = {'tip3p_2004':'TIP3P_2004','tip3p_1983_charmm_hybrid':'tip3p_1983_charmm_hybrid'}

def make_box_of_molecules(molecules, num_molecules, box,
                          tolerance = 2.0, outputfile = 'box.pdb',
                          clean_folder = True,
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
    #os.sysitem('packmol < pk.inp > pk.log')
    run_packmol()
    atoms = io.read(outputfile)
    atoms.cell = box
    if clean_folder == True:
        os.system('rm pk* *.pdb')
    return atoms


def run_packmol():
    """
    executes packmol and makes sure to let the user know to cite
    packmol (giving credit is important!)
    """
    print('you are using a function that utilizes packmol, please'
          ' ensure you cite packmol in your paper.\n L. Martínez, '
          'R. Andrade, E. G. Birgin, J. M. Martínez. Packmol: A '
          'package for building initial configurations for '
          'molecular dynamics simulations. Journal of Computational '
          'Chemistry, 30(13):2157-2164, 2009.')
    os.system('packmol < pk.inp > pk.log')

def equilibrate_tip3p(atoms):
    """
    THIS FUCNTION DOES NOT WORK YET
    takes a box of water and equilibrates it using TIP3P. Bonds
    are infered. Much of this is taken from the ASE tutorial on
    equilibrating TIPnP boxes:
    https://wiki.fysik.dtu.dk/ase/tutorials/tipnp_equil/tipnp_equil.html

    """
    from ase.constraints import FixBondLengths
    from ase.md import Langevin
    from ase.units import fs, kB
    from ase.calculators.tip3p import TIP3P
    analysis = Analysis(atoms)
    bonds = analysis.get_bonds('H','O')
    bonds = [a for a in bonds if len(a) == 2]

    atoms.set_constraint(FixBondLengths(bonds))
    atoms.set_calculator(TIP3P(rc = 6))
    md = Langevin(atoms, 1 * fs, temperature=300 * kB,
              friction=0.01)
    from ase.io.trajectory import Trajectory
    traj = Trajectory('test.traj', 'w', atoms)
    md.attach(traj.write, interval=1)
    md.run(4000)
    return atoms

def write_lammps_inputs_moltemplate(atoms, force_field, num_molecules,
                                    clean_folder = True):
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
    import subprocess
    f = open('mt.log','w')
    p = subprocess.call('moltemplate.sh -pdb mt.pdb -atomstyle full mt.lt',
                        shell = True, stdout = f, stderr = f)
    if clean_folder:
        os.system('rm -r output_ttree')
    #os.system('moltemplate.sh -pdb mt.pdb -atomstyle full mt.lt &> mt.log')

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
    LAMMPS dumps out. This is needed to translate back from lammps to ase

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

def center_traj(traj, output_name = 'fixed.traj'):
    """
    center the atoms in a trajectory

    Parameters:
        traj (list):
            a list of ase atoms objects that you'd like to center
        output_name (str):
            the name you'd like to give the output file

    returns:
        None
    """

    from ase.io import read, write
    from ase.calculators.singlepoint import SinglePointCalculator as sp

    new_traj = []
    for image in traj:
        energy = image.get_potential_energy()
        forces = image.get_forces()
        if sum(sum(image.cell)) == 0:  # Crudely check if it has a cell
            image.set_cell([10,10,10])
        image.center()
        image.set_pbc([False] * 3)
        image.set_calculator(sp(atoms = image,
                                energy = energy,
                                forces = forces))
        new_traj.append(image)

    write(ouput_name, new_traj)


def parse_custom_dump(dump, datafile, label='atoms',
                      energyfile=None, write_traj=False,
                      units='real'):
    """
    This function parses the output of a LAMMPS custom dump file. Currently
    this function assumes you've put in this for the custom dump: 
    "element x y z fx fy fz". There are plans to expand this to be more general
    if the need arises. This assumes units are set to real. Only works with
    orthogonal unit cells currently

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
    from ase.units import kcal, mol, fs
    f = open(dump,'r')
    text = f.read()
    #if 'type x y z fx fy fz' not in text:
    #    raise Exception('the dump file is not in the correct format.'
    #                    ' Please make the dump "type x y z fx fy fz".'
    #                    ' Further functionality is planned to be added.')
    if datafile is not None:
        element_key = elements_from_datafile(datafile)
    if units == 'metal':
        ps = fs * 1000 # I /think/ this is right
        F_conversion = 1
        E_conversion = 1
        V_conversion =(1 / ps)
    elif units == 'real':
        F_conversion = (kcal / mol)
        E_conversion = (kcal / mol)
        V_conversion = (1 / fs)


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
            fx = dump_format.index('fx')
            fy = dump_format.index('fy')
            fz = dump_format.index('fz')
            forces = data[:, [fx, fy, fz]]
            forces *= F_conversion  # convert to eV/A
        # parse the potential energy
        if 'c_energy' in dump_format:
            eng = dump_format.index('c_energy')
            per_atom_energy = data[:,eng] * E_conversion  # convert to eV/A
            atoms.set_initial_charges(per_atom_energy)
            energy = sum(per_atom_energy)
        # parse the velocities
        if 'vx' in dump_format and 'vy' in dump_format and 'vz' in dump_format:
            vx = dump_format.index('fx')
            vy = dump_format.index('fy')
            vz = dump_format.index('fz')
            velocities = data[:, [vx, vy, vz]]
            velocities *= V_conversion  # convert to A/[AU time]
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
        if 'velocities' in locals():
            atoms.set_velocities(velocities)
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
    to be used as a reaxff input. This is useful if you've generated
    a lammps input file with moltemplate and need to remove the bonding
    information.

    Parameters:
        filename (str):
            The name of the file

    returns
        None
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
    
def make_standard_input(calculation_type='reaxff',
                        ff_file='TiO2_water.reax',
                        timestep=0.2, 
                        ensemble='nvt',
                        steps=20000,
                        temp=300,
                        elements=None):
    """
    TODO: make sure that the atoms are placed in order in the .in file
    function to just generate a generic input file to run a basic
    simulation. No frills, it just works.
    """
    from .standard_inputs import input_files

    if calculation_type not in input_files.keys():
        raise Exception('the calculation type you selected, {}, is not implemented yet'.format(calculation_type))
     
    from .standard_inputs import input_files

    atoms_in_order = ' '.join(elements)

    with open('lmp.input','w') as f:
        if 'single_point_reaxff' in calculation_type:
            atoms_in_order = ' '.join(elements)
            f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     atoms_in_order))
        if 'single_point_simple_nn' in calculation_type:
            atoms_in_order = ' '.join(elements)
            f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     atoms_in_order))
        if 'single_point_eam' in calculation_type:
            atoms_in_order = ' '.join(elements)
            f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     atoms_in_order))
        if 'eam' in calculation_type:
            atoms_in_order = ' '.join(elements)
            f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     atoms_in_order,
                                                     temp,
                                                     steps))

        else:
            f.write(input_files[calculation_type].format(
                                                     ff_file,
                                                     timestep, 
                                                     temp, 
                                                     steps))
    if 'reaxff' in calculation_type:
        import inspect
        cwd = os.getcwd()
        package_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        ff = os.path.join(package_directory, 'special_input_files/' + ff_file)
        control = os.path.join(package_directory, 'special_input_files/control.reaxff')
        copyfile(ff, os.path.join(cwd, ff_file))
        copyfile(control, os.path.join(cwd, 'control.reaxff'))
    elif 'eam' in  calculation_type:
        import inspect
        cwd = os.getcwd()
        package_directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        ff = os.path.join(package_directory, 'special_input_files/' + ff_file)
        copyfile(ff, os.path.join(cwd, ff_file))


    
    

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

    number_of_molecules = int(np.floor((atoms.get_volume() - volume) * 0.0333679)) # approximate density of water

    atoms = make_box_of_molecules([fluid_molecule,atoms], [number_of_molecules,1], atoms.cell)
    # this next part is such a mess, I'm so sorry
    # all this is doing is trying to make the particle more centered in the unit cell
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

def prune_oxygens(atoms, metal=None, bond_distance=2.2):
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
    TODO: make this code work for non-orthogonal slabs
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
    if number_of_molecules < 0:
        raise ValueError('A negative number of fluid molecules was calculated.'
                         ' Make sure your slab has a vacuum layer')
    # Make the box
    water_layer = make_box_of_molecules([fluid_molecule], [number_of_molecules],
                                        box = [length, width, height - offset * 2])
    water_layer.positions += np.array([0, 0, highest_atom + offset])
    slab = atoms + water_layer
    return slab


def make_rdf_based_descriptors(images, n_descriptors=20,
                               cutoff=6.5, 
                               fall_off_percent=0.2,
                               descriptor_type='simple-nn',
                               nbins=200,
                               plot=False):
    """
    generates reasonable values for eta and rs for a given image in a trajectory
    based on the radial distribution function.

    Parameters:
        images (list):
            a list of ASE atoms objects to be used in the radial
            distribution function calculation
        n_descriptors (int):
            the number of g2 descriptors you want to use
        fall_off_percent (float):
            the value of desriptor_n/descriptor_n+1 at the center of
            descriptor_n+1. Higher values are more smeared gaussians
            lower values are tighter gaussians
        descriptor_type (str):
            the program into which you will be putting these eta
            values. simplenn and amp have different conventions for
            what eta is. possible values are `amp` and `simplenn`
        nbins (int):
            how many bins you want the radial distribution function
            to be calculated with
        plot (bool):
            If True, the descriptors will be plotted over top of the
            radial distribution function
    """
    from scipy.integrate import trapz
    from ase.ga.utilities import get_rdf

    #if type(images) == list:
    #    atoms = images[index]
    #else:
    #    atoms = images

    analysis = Analysis(images)
    rdf = analysis.get_rdf(rmax = cutoff, nbins = nbins,
                           return_dists = True
                           )

    total_rdf = np.empty((nbins,), dtype = 'float64')
    for image, distances in rdf:
        total_rdf += image
    total_rdf /= len(rdf)
    rdf = total_rdf

    #rdf = get_rdf(atoms, rmax = cutoff, nbins = 200)
    #rdf_1, distances = rdf
    # ignore the stupid code in the next few lines
    """
    # find 75% index
    cut = int(-1 * np.ceil(len(rdf) * 0.75))
    rdf_mean = np.mean(rdf_1[cut:])
    # subtract off the mean and make all positive
    #rdf_1 -= rdf_mean
    rdf_1 = abs(rdf_1)
    """
    integrals = []
    for i in range(len(rdf)):
        integrals.append(trapz(rdf[:i+1]))
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
    descriptor_distances[0] -= 0.08 # offset the first one slightly
    etas = []
    rs_s = [0] * n_descriptors
    for distance in descriptor_distances:
        etas.append(-1 * np.log(fall_off_percent) / distance ** 2)
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
        plt.plot(distances, np.array(rdf))
        for rs, eta in zip(rs_s, etas):
            x = np.linspace(0, cutoff, 1000)
            y = [np.exp(-1 * eta * (a - rs) ** 2) for a in x]
            plt.plot(x,np.array(y) * max(rdf), alpha = 0.4)
        plt.xlabel('Distance (Angstrom)')
        plt.ylabel('Total RDF')
        plt.title('Radial Distribution Function and Derived Descriptors')
        plt.savefig('rdf.png')
        plt.show()

    
    if descriptor_type.lower() in ['amp']:
        etas = [a / cutoff for a in etas]
    return etas, rs_s


def gaussian_basis(x, a, xk, sigma):
    """
    helper function for gaussian_fit_descriptors

    a: the aplitude
    xk: the x offset
    sigma: the standard deviation
    """
    return a * np.exp( -1 * ((x - xk) ** 2 / (2 * sigma ** 2)))


def n_sized_gaussian(x, *a):
    """
    helper function for gaussian_fit_descriptors
    """

    n = int(len(a) / 3)
    a = np.array(a)
    a = np.reshape(a, (n, 3))
    s = 0
    for row in a:
        s += gaussian_basis(x, *list(row))
    return s


def gaussian_fit_descriptors(traj, n_gaussians=5, cutoff=6.5, 
                             nbins=10, plot=False):
    """
    approximates the RDF using a sum of gaussian functions then 
    uses the centers and standard deviations of those gaussians to
    generate descriptors. This is still under development

    Parameters:
        traj (list):
            a list of ASE atoms objects with which you'd like to 
            calculate the radial distribution function
        n_gaussians (int):
            the number of gaussians you'd like to fit the radial
            distribution function to
        cutoff (float):
            the distance in angstroms at which you'd like to cut 
            off the descriptors 
        nbins (int):
            the number of bins you'd like to calculate the radial
            distribution function with.
        plot (bool):
            If set to True, a plot will be displayed showing the
            radial distribution function with the gaussians 
            overlayed.

    returns:
        None
    """
    from ase.ga.utilities import get_rdf
    from scipy.optimize import curve_fit
    from sklearn.neighbors import KernelDensity


    n = n_gaussians

    atoms = traj[0]
    
    full_dists = []
    for i in range(len(atoms)):
        indices = list(range(len(atoms)))
        indices.remove(i)
        dists = atoms.get_distances(i, indices)
        dists = [a for a in dists if a < cutoff]
        full_dists += dists

    full_dists = np.array(full_dists)
    full_dists = full_dists.reshape((len(full_dists), 1))

    kde = KernelDensity(kernel='gaussian',bandwidth=0.2).fit(full_dists)
    dens = kde.score_samples(np.reshape(np.linspace(0.8,cutoff,1000), (1000, 1)))

    analysis = Analysis(traj)
    rdf = analysis.get_rdf(rmax = cutoff, nbins = nbins, 
                           return_dists = True
                           )

    total_rdf = np.empty((nbins,), dtype = 'float64')
    for image, distances in rdf:
        total_rdf += image
    total_rdf /= len(rdf)


    # cut off anything smaller than 0.8 A for stability
    distances_2 = [a for a in distances if a > 0.8]
    index = list(distances).index(distances_2[0])
    total_rdf = total_rdf[index:]
    distances = distances_2
    del distances_2

    tries = 0
    while tries < 1000:

        a0 = np.array([7,1.,0.2]) + np.random.random(3) * 0.01 # a good first guess
        if n > 1:
            a0 = np.append(a0, np.array([2,1.6,0.5]) + np.random.random(3) * 0.01)
        if n > 2:
            a0 = np.append(a0, np.random.random((n - 2) * 3) * np.array([2, 2, 0.1] * (n - 2)) + \
                                   np.array([1.2, 2, 0.2] * (n - 2)))

        a0 = list(a0)
        try:
            opt_params, covariance = curve_fit(n_sized_gaussian, distances, total_rdf,
                                       p0 = a0,
                                       maxfev = 10 ** 3 # Allow tons of evaluations
                                       )
        except RuntimeError:
            tries += 1
    sigmas = opt_params[1::3]
    rs_s = opt_params[0::3]

    # this this the relationship between sigma and eta
    etas = [1 / (2 * a ** 2) for a in sigmas]
    print(etas, rs_s)

    if plot == True:
        from matplotlib import pyplot as plt
        plt.plot(distances, total_rdf)
        #plt.plot(np.linspace(0.8,cutoff,1000), dens)
        plt.plot(distances, n_sized_gaussian(distances, *opt_params))
        for rs, eta in zip(rs_s, etas):
            x = np.linspace(0,cutoff,1000)
            y = [np.exp(-1 * eta * (a - rs) ** 2) for a in x]
            plt.plot(x,y)
        #plt.plot(distances, n_sized_gaussian(distances, *a0))
        plt.show()


def make_params_file(elements, etas, rs_s, g4_eta = 4, cutoff = 6.5,
                     g4_zeta=[1.0, 4.0], g4_gamma=[1, -1],
                     convert_from_amp=False):
    """
    makes a params file for simple_NN. This is the file containing
    the descriptors. This function makes g2 descriptos for the eta
    and rs values that are input, and g4 descriptors that are log
    spaced between 10 ** -5 and 10 ** -1. The number of these
    that are made is controlled by the `n_g4_eta` variable

    Parameters:
        elements (list):
            a list of elements for which you'd like to make params
            files for
        etas (list):
            the eta values you'd like to use for the descriptors
        rs_s (list):
            a list corresponding to `etas` that contains the rs
            values for each descriptor
        g4_eta (int or list):
            the number of g4 descriptors you'd like to use. if a
            list is passed in the values of the list will be used
            as eta values
        cutoff (float):
            the distance in angstroms at which you'd like to cut 
            off the descriptors
    
    returns:
        None

    
    """
    if len(etas) != len(rs_s):
        raise ValueError('the length of etas and rs_s must be the same')
    if type(g4_eta) == int:
        g4_eta = np.logspace(-4, -1, num = g4_eta)
    if convert_from_amp:
        etas = [a * cutoff for a in etas]
        g4_eta = [a * cutoff for a in g4_eta]
    for element in elements:
        with open('params_{}'.format(element),'w') as f:
            # G2
            for species in range(1, len(elements) + 1):
                for eta, Rs in zip(etas, rs_s):
                    f.write('2 {} 0 {} {} {} 0.0\n'.format(species, cutoff,
                                                           np.round(eta, 6), Rs))
            # G4 
            for i in range(1, len(elements) + 1):
                n = i
                while True:
                    for eta in g4_eta:
                        for gamma in g4_gamma:
                            for zeta in g4_zeta:
                                f.write('4 {} {} {} {} {} {}\n'.format(i, n, cutoff,
                                                                       np.round(eta, 6),
                                                                       zeta, gamma))
                    n += 1
                    if n > len(elements):
                        break

def reorganize_simple_nn_derivative(image, dx_dict):
    """
    reorganizes the fingerprint derivatives from simplen_nn into
    amp format

    Parameters:
        image (ASE atoms object):
            the atoms object used to make the finerprint
        dx_dict (dict):
            a dictionary of the fingerprint derivatives from simple_nn

    """
    # TODO check for bugs
    d = defaultdict(list)
    sym_dict = defaultdict(list)
    syms = image.get_chemical_symbols()
    for i, sym in enumerate(syms):
        sym_dict[sym].append(i)
    # the structure is:
    # [elements][atom i][symetry function #][atom j][derivitive in direction]
    for element, full_arr in dx_dict.items():
        for i, arr_t in enumerate(full_arr):
            true_i = sym_dict[element][i]
            for sf in arr_t:
                for j, dir_arr in enumerate(sf):
                    for k, derivative in enumerate(dir_arr):
                        d[(true_i,element,j,syms[j],k)].append(derivative)
    zero_keys = []
    for key, derivatives in d.items():
        zero_check = [a == 0 for a in derivatives]
        if zero_check == [True] * len(derivatives):
            zero_keys.append(key)
    for key in zero_keys:
        del d[key]
    d = dict(d)
    return d

def reorganize_simple_nn_fp(image, x_dict):
    """
    reorganizes the fingerprints from simplen_nn into
    amp format

    Parameters:
        image (ASE atoms object):
            the atoms object used to make the finerprint
        x_dict (dict):
            a dictionary of the fingerprints from simple_nn

    """
    # TODO check for bugs
    # the structure is:
    # [elements][atom i][symetry function #][fp]
    fp_l = []
    sym_dict = defaultdict(list)
    syms = image.get_chemical_symbols()
    for i, sym in enumerate(syms):
        sym_dict[sym].append(i)
    for element, full_arr in x_dict.items():
        for i, fp  in enumerate(full_arr):
            true_i = sym_dict[i]
            fp_l.append((element,list(fp)))
    return fp_l

def get_hash(atoms):
    import hashlib
    """Creates a unique signature for a particular ASE atoms object.

    This is used to check whether an image has been seen before. This is just
    an md5 hash of a string representation of the atoms object.

    Parameters
    ----------
    atoms : ASE dict
        ASE atoms object.

    Returns
    -------
        Hash string key of 'atoms'.
    """
    string = str(atoms.pbc)
    try:
        flattened_cell = atoms.cell.array.flatten()
    except AttributeError:  # older ASE
        flattened_cell = atoms.cell.flatten()
    for number in flattened_cell:
        string += '%.15f' % number
    for number in atoms.get_atomic_numbers():
        string += '%3d' % number
    for number in atoms.get_positions().flatten():
        string += '%.15f' % number

    md5 = hashlib.md5(string.encode('utf-8'))
    hash = md5.hexdigest()
    return hash


def convert_simple_nn_fps(traj, delete_old=True):
    from multiprocessing import Pool
    # make the directories
    if not os.path.isdir('./amp-fingerprints.ampdb'):
        os.mkdir('./amp-fingerprints.ampdb')
    if not os.path.isdir('./amp-fingerprints.ampdb/loose'):
        os.mkdir('./amp-fingerprints.ampdb/loose')
    if not os.path.isdir('./amp-fingerprint-primes.ampdb'):
        os.mkdir('./amp-fingerprint-primes.ampdb')
    if not os.path.isdir('./amp-fingerprint-primes.ampdb/loose'):
        os.mkdir('amp-fingerprint-primes.ampdb/loose')
    # perform the reorganization
    """
    for i, image in enumerate(traj):
        pic = pickle.load(open('./data/data{}.pickle'.format(i + 1), 'rb'))
        im_hash = get_hash(image)
        x_list = reorganize_simple_nn_fp(image, pic['x'])
        pickle.dump(x_list, open('./amp-fingerprints.ampdb/loose/' + im_hash, 'wb'))
        del x_list  # free up memory just in case
        x_der_dict = reorganize_simple_nn_derivative(image, pic['dx'])
        pickle.dump(x_der_dict, open('./amp-fingerprint-primes.ampdb/loose/' + im_hash, 'wb'))
        del x_der_dict  # free up memory just in case
        if delete_old:  # in case disk space is an issue
            os.remove('./data/data{}.pickle'.format(i + 1))
    """
    with Pool(10) as p:
        l_trajs = list(enumerate(traj))
        p.map(reorganize, l_trajs)
    if delete_old:
        os.rmdir('./data')

def reorganize(inp, delete_old=True):
    i, image = inp
    pic = pickle.load(open('./data/data{}.pickle'.format(i + 1), 'rb'))
    im_hash = get_hash(image)
    x_list = reorganize_simple_nn_fp(image, pic['x'])
    pickle.dump(x_list, open('./amp-fingerprints.ampdb/loose/' + im_hash, 'wb'))
    del x_list  # free up memory just in case
    x_der_dict = reorganize_simple_nn_derivative(image, pic['dx'])
    pickle.dump(x_der_dict, open('./amp-fingerprint-primes.ampdb/loose/' + im_hash, 'wb'))
    del x_der_dict  # free up memory just in case
    if delete_old:  # in case disk space is an issue
        os.remove('./data/data{}.pickle'.format(i + 1))


class DummySimple_nn(object):
    """
    a dummy class to fool the simple_nn descriptor class into
    thinking it's attached to a simple_nn instance
    """
    def __init__(self, atom_types):
        self.inputs = {
            'generate_features': True,
            'preprocess': False,
            'train_model': True,
            'atom_types': atom_types}
        self.logfile = open('simple_nn_log', 'w')


def make_simple_nn_fps(traj, descriptors, clean_up_directory=False,
                       elements='all'):
    """
    generates descriptors using simple_nn. The files are stored in the
    ./data folder. These descriptors will be in the simple_nn form and
    not immediately useful for other programs

    Parameters:
        traj (list of ASE atoms objects):
            a list of the atoms you'd like to make descriptors for
        descriptors (tuple):
            a tuple containing (g2_etas, g2_rs_s, g4_etas, cutoff, g4_zetas, g4_gammas)
        clean_up_directory (bool):
            if set to True, the input files made by simple_nn will
            be deleted

    returns:
        None
    """
    from simple_nn.features.symmetry_function import Symmetry_function

    # handle inputs
    if type(traj) != list:
        traj = [traj]

    # clean up any previous runs
    if os.path.isdir('./data'):
        shutil.rmtree('./data')

    # set up the input files
    io.write('simple_nn_input_traj.traj', traj)
    with open('str_list', 'w') as f:
        f.write('simple_nn_input_traj.traj :') # simple_nn requires this file

    if elements == 'all':
        atom_types = []
        # TODO rewrite this
        for image in traj:
            atom_types += image.get_chemical_symbols()
            atom_types = list(set(atom_types))
    else:
        atom_types = elements

    make_params_file(atom_types, *descriptors, convert_from_amp=True)

    # build the descriptor object
    descriptor = Symmetry_function(fp_dir='.')
    params = {a:'params_{}'.format(a) for a in atom_types}

    descriptor.inputs = {'params': params, 
                         'refdata_format': 'traj', 
                         'compress_outcar': False,
                         'data_per_tfrecord': 150, 
                         'valid_rate': 0.1, 
                         'remain_pickle': False, 
                         'continue': False, 
                         'add_atom_idx': True, 
                         'num_parallel_calls': 5, 
                         'atomic_weights': {'type': None, 'params': {}}, 
                         'weight_modifier': {'type': None, 'params': {}}, 
                         'scale_type': 'minmax', 
                         'scale_scale': 1.0, 
                         'scale_rho': None}
    dummy_class = DummySimple_nn(atom_types=atom_types)
    descriptor.parent = dummy_class

    # generate the descriptors
    descriptor.generate(label='')
    
    if False:
        # clean the folder of all the junk
        files = ['simple_nn_input_traj.traj', 'str_list',
                 'pickle_list', 'simple_nn_log']
        files += list(params.values())
        for file in files:
            os.remove(file)

def make_amp_descriptors_simple_nn(traj, g2_etas, g2_rs_s, g4_etas, g4_zetas, g4_gammas, cutoff,
                                   elements='all',descriptor_type='simple_nn'):
    """
    uses simple_nn to make descriptors in the amp format.
    Only creates the same symmetry functions for each element
    for now.
    """
    if descriptor_type == 'amp':
        g2_etas = [a / cutoff ** 2 for a in g2_etas]
        g4_etas = [a / cutoff ** 2 for a in g4_etas]

    c = cutoff
    make_simple_nn_fps(traj,
                       (g2_etas, g2_rs_s, g4_etas, 
                        cutoff, 
                        g4_zetas, g4_gammas),
                        clean_up_directory=True,
                        elements=elements)
    convert_simple_nn_fps(traj, delete_old=True)

def make_fingerprint_matrix(traj, descriptors, clean_up_directory=True,
                            elements='all', return_image_inds=False):
    """
    make a massive array of the fingerprints of a trajectory

    inputs:
        traj: (list of ASE atoms objects)
            The list of atoms objects you'd like to include in the fingerprint
            matrix.
        descriptors: (tuple)
            a tuple containing (g2_etas, g2_rs_s, g4_etas, cutoff, g4_zetas, g4_gammas)
        clean_up_directory: (bool)
            If set to True, the input files used will be deleted after calculations
            are complete
        elements: (list or str)
            if a list of elements is passed in, fingerprints will only be calculated
            for those elements (other atoms are just ignored, as if the did not exist)
            if 'all' is passed in all elements found in the trajectory are included
        return_image_inds: (bool)
            If set to True, the indicies of the original image for each fingerprint
            are also returned

    returns:
        arrays_dict: (dict)    
            a dictionary containing the chemical elements in the images as keys
            and the fingerprints of those elements as values
        image_inds: (dict)
            a dictionary with the same keys as `arrays_dict` but whose values
            contain the index of the input variable `traj` where that finerprint
            is found

    """
    arrays_dict = {}
    image_inds = defaultdict(list)
    if elements == 'all':
        elements = []
        for image in traj:
            elements += image.get_chemical_symbols()
            elements = list(set(elements))


    # calculate how many different combinations of atoms there are
    combinations = itertools.combinations_with_replacement(elements, 2)
    num_combinations = len([a for a in combinations])

    fp_len = len(descriptors[0]) * len(elements) + \
             len(descriptors[2]) * len(descriptors[4]) * len(descriptors[5]) * num_combinations

    for element in elements:
        arrays_dict[element] = np.array([], dtype='float64').reshape(0, fp_len)
    make_simple_nn_fps(traj, descriptors, clean_up_directory=clean_up_directory, 
                       elements=elements)
    for i, image in enumerate(os.listdir('./data')):
        dat = pickle.load(open('./data/' + image, 'rb'))
        for element in elements:
            if element not in dat['x'].keys():
                continue
            array = np.array(dat['x'][element], dtype='float64')
            current = arrays_dict[element]
            new = np.concatenate((current, array), axis=0)
            image_inds[element] += [i] * len(array)
            arrays_dict[element] = new
    
    arrays_dict['total'] = np.empty((0, fp_len))

    for array in arrays_dict.values():
        arrays_dict['total'] = np.concatenate((arrays_dict['total'], array), axis=0)
    if return_image_inds:
        return arrays_dict, image_inds
    else:
        return arrays_dict


def extract_rdf(filename, plot = False):
    """
    pulls the rdf data out of the structured LAMMPS output and
    converts it to a pandas dataframe, optionally plotting it.

    inputs:
        filename:
            the name of the rdf file
        plot:
            if set the True, it will plot the first RDF

    returns:
        None
    """
    import pandas as pd
    with open(filename, 'r') as f:
        text = f.read()
    rdfs = text.split('# TimeStep Number-of-rows\n')[1:]
    rdf_dfs = []
    for rdf in rdfs:
        rdf = rdf.replace('# ', '')
        rdf = rdf.split('\n')
        del rdf[1]
        if '#' in rdf[-1]:
            del rdf[-1]
        rdf = '\n'.join(rdf)
        data = StringIO(rdf)
        df = pd.read_csv(data, sep = ' ')
        rdf_dfs.append(df)

    if plot == True:
        import matplotlib.pyplot as plt
        plt.plot(df['c_myRDF[1]'],df['c_myRDF[2]'])
        plt.title('Radial Distribution Function Water Reaxff')
        plt.ylabel('g(r)')
        plt.xlabel('r ($\AA$)')
        plt.show()

    return df



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

def kernel_density_radial_distribution_function(traj, bandwidth = 0.2, 
                                                cutoff = 6.5, 
                                                nbins = 100,
                                                plot = False):
    """
    approximates the RDF using kernel density estimation

    Parameters:
        traj (list):
            a list of ASE atoms objects with which you'd like to 
            calculate the radial distribution function
        bandwidth (float):
            the bandwidth to be used for the kernel density
            appoximation
        cutoff (float):
            the distance in angstroms at which you'd like to cut 
            off the descriptors
        nbins (int):
            the number of bins you'd like to calculate your radial
            distribution function with
        plot (bool):
            if set to True, a plot will be returned of the 
            radial distribution function and the kernel density
            approximation

    returns:
        None
    """
    from ase.ga.utilities import get_rdf
    from scipy.optimize import curve_fit
    from sklearn.neighbors import KernelDensity

    if type(traj) != list:
        traj = [traj]

    full_dists = []
    for atoms in traj:
        for i in range(len(atoms)):
            indices = list(range(len(atoms)))
            indices.remove(i)
            dists = atoms.get_distances(i, indices)
            dists = [a for a in dists if a < cutoff]
            full_dists += dists

    full_dists = np.array(full_dists)
    full_dists = full_dists.reshape((len(full_dists), 1))


    analysis = Analysis(traj)
    rdf = analysis.get_rdf(rmax = cutoff, nbins = nbins,
                           return_dists = True
                           )

    total_rdf = np.empty((nbins,), dtype = 'float64')
    for image, distances in rdf:
        total_rdf += image
    total_rdf /= len(rdf)


    distances_kde = np.linspace(0.7, cutoff, 1000)

    kde = KernelDensity(kernel='epanechnikov', bandwidth = bandwidth).fit(full_dists)
    dens = kde.score_samples(np.reshape(distances_kde, (1000, 1)))

    dens -= min(dens)

    #for i, d in enumerate(dens):
    #    dens[i] = d / distances_kde[i] ** 2 / np.pi

    if plot == True:
        from matplotlib import pyplot as plt
        plt.plot(distances, total_rdf)
        plt.plot(np.linspace(0.7,cutoff,1000), dens)
        plt.show()

def atoms_to_json(atoms):
    """
    convert an atoms object to a json form
    """
    with StringIO() as f:
        io.write(f, atoms, format = 'json')
        json_str_atoms = f.getvalue()

    with StringIO(json_str_atoms) as g:
        json_atoms = json.load(g)
    json_atoms = json_atoms['1']
    superfluous_entries = ['ctime', 'mtime', 'unique_id', 'user']
    for entry in superfluous_entries:
        json_atoms.pop(entry, None)
    return json_atoms

def json_to_atoms(json_atoms):
    """
    convert back from json to an atoms object
    """
    json_atoms = {'1': json_atoms,
                  'ids': [1],
                  'nextid': 2}
    json_str_atoms = json.dumps(json_atoms)
    with StringIO(json_str_atoms) as f:
        atoms = io.read(f, format='json')
    return atoms

def rereference_traj(traj, reference_energy):
    """
    subtracts off the reference energy provided from all the other
    energies and returns the result

    Parameters:
        traj (list of ASE atoms objects):
            a list of images you'd like to re-reference
        reference_energy (float):
            the reference energy you'd like to subtract from the images

    returns:
        traj (list of ASE atoms objects):
            the re-referenced trajectory

    """
    from ase.calculators.singlepoint import SinglePointCalculator as sp
    from ase.io.trajectory import Trajectory

    # in case they didn't read the documenation well
    if type(traj) != list:
        traj = [traj]

    for image in traj:
        frc = image.get_forces()
        eng = image.get_potential_energy()
        image.set_calculator(sp(image, energy = eng - reference_energy,
                            forces = frc))
    return traj

def restart_simple_nn(num):
    e1 = os.system('cp SAVER_epoch{}.meta SAVER.meta'.format(num))
    e2 = os.system('cp SAVER_epoch{}.data-00000-of-00001 SAVER.data-00000-of-00001'.format(num))
    e3 = os.system('cp SAVER_epoch{}.index SAVER.index'.format(num))
    if e1 or e2 or e3:
        raise RuntimeError('could not restart simple_nn')

def fix_pbc(traj):
    from ase.calculators.singlepoint import SinglePointCalculator as sp
    if type(traj) != list:
        traj = [traj]
    for image in traj:
        energy = image.get_potential_energy()
        forces = image.get_forces()
        image.set_pbc([True] * 3)
        image.wrap(pbc=[True] * 3)
        image.set_calculator(sp(image, energy = energy,
                            forces = forces))
    return traj

def clean_traj(traj):
    from ase.calculators.singlepoint import SinglePointCalculator as sp
    if type(traj) != list:
        traj = [traj]
    energies = np.array([a.get_potential_energy() for a in traj])
    forces = [a.get_forces() for a in traj]
    new_images = []
    std_energy = np.std(energies)
    median_energy = np.median(energies)
    for image, energy, force in zip(traj, energies, forces):
        if abs(energy - median_energy) > 3 * std_energy:
            continue
        image.set_calculator(sp(image, energy = energy,
                            forces = force))
        new_images.append(image)
    return new_images


def parse_simple_nn_log(directory='.', plot=False):
    import matplotlib.pyplot as plt
    with open('LOG', 'r') as f:
        txt = f.read()

    epochs = txt.split('epoch')[1:]
    train_engs, test_engs, train_frcs, test_frcs = [], [], [], []

    for epoch in epochs:
        eng = epoch.split('E RMSE(T V) = ')[1]
        eng, frc = eng.split('F RMSE(T V) =')
        frc = frc.split(' learning_rate')[0]

        train_eng, test_eng = tuple([float(a) for a in eng.split()])
        train_frc, test_frc = tuple([float(a) for a in frc.split()])
        train_engs.append(train_eng)
        test_engs.append(test_eng)
        train_frcs.append(train_frc)
        test_frcs.append(test_frc)

    x = range(1, len(train_engs) + 1)

    # generate a moving average
    avg_test_engs = []
    avg_test_frcs = []
    for i, eng in enumerate(test_engs):
        if i < 19:
            avg_test_engs.append(np.mean(test_engs[:i+1]))
            avg_test_frcs.append(np.mean(test_frcs[:i+1]))
        else:
            avg_test_engs.append(np.mean(test_engs[i-19:i+1]))
            avg_test_frcs.append(np.mean(test_frcs[i-19:i+1]))

    if plot:
        plt.rcParams["figure.figsize"] = (15,5)
        fig, _axs = plt.subplots(nrows = 1, ncols = 2)
        axs = _axs.flatten()
        axs[0].plot(x, test_frcs, label='learning rate')
        axs[0].plot(x, avg_test_frcs, label='moving average')
        axs[0].set_title('Learning Curve')
        axs[0].set_ylabel('Test Set Force RMSE (eV/A)')
        axs[0].set_xlabel('Checkpoint')
        axs[0].set_ylim([0,1])
        axs[1].plot(x, test_engs, label='learning rate')
        axs[1].plot(x, avg_test_engs,label='moving average')
        axs[1].set_title('Learning Curve')
        axs[1].set_ylabel('Test Set Energy RMSE (eV)')
        axs[1].set_xlabel('Checkpoint')
        #axs[1].set_ylim([0,0.1])
        axs[1].legend()
        plt.savefig('learning.png')
        plt.show()

        fig = plt.figure()
        plt.plot(x, test_engs)
        #plt.show()
    Results = namedtuple('Results', 'train_engs test_engs train_frcs test_frcs')
    return Results(train_engs, test_engs, train_frcs, test_frcs)

def change_yaml_line(line_text, filename):
    """
    helper function to modify a single line of a .yaml
    file

    Parameters:
        line_text:
            the text you'd like the line to be, this is
            used to find the relevant line in the .yaml
            file
        filename:
            the name of the file you'd like to modify
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    new_lines = []
    line_identifier = line_text.split(':')[0].strip()
    for line in lines:
        if line_identifier in line:
            new_lines.append(line_text + '\n')
        else:
            new_lines.append(line)
    with open(filename, 'w') as f:
        f.writelines(new_lines)

def change_num_epoch_simple_nn(num_epoch, filename='input.yaml'):
    """
    a function to change the number of epochs being run in a
    simple_nn input file.
    """
    change_yaml_line('  total_epoch: {}'.format(num_epoch), filename)

def toggle_continue_simple_nn(cont, filename='input.yaml'):
    """
    allows you to toggle if simple_nn re-imports the last
    model
    """
    if cont == 'yes':
        change_yaml_line('  continue: true', filename)
    elif cont == 'no':
        change_yaml_line('  continue: false', filename)

def write_fp_code_input(atoms):
    os.mkdir('fp_calc')
    os.chdir('./fp_calc')

    with open('cell_martix', 'w') as f:
        for line in atoms.get_cell():
            for entry in line:
                f.write(str(entry) + ' ')
            f.write('\n')

    with open('positions', 'w') as f:
        for atom_pos in atoms.get_scaled_positions():
            for entry in atom_pos:
                f.write(str(entry) + ' ')
            f.write('\n')

    with open('elements', 'w') as f:
        for element in atoms.get_chemical_symbols():
            f.write(str(element) + '\n')

    with open('element_alist', 'w') as f:
        elements = list(set(atoms.get_chemical_symbols()))
        elements.sort()
        for element in elements:
            f.write(element + '\n')

    os.chdir('../')


def atomic_parity_plot(traj1, traj2):
    """
    makes a parity plot between the energies of two trajectories
    which are assumed to be the same elementwise
    """
    from matplotlib import pyplot as plt
    E1 = []
    E2 = []
    for i1, i2, in zip(traj1, traj2):
        E1.append(i1.get_potential_energy())
        E2.append(i2.get_potential_energy())
    plt.scatter(E1, E2)
    plt.plot([min(E1 + E2)] * 2, [max(E1 + E2)] * 2)
    plt.show()

def amp_parity_plot(amp_calc='amp-checkpoint.amp', images=None,
                    training_images=None, per_atom=True):
    from amp import Amp
    from matplotlib import pyplot as plt
    true_engs = []
    amp_engs = []
    for run_set in [images, training_images]:
        if run_set is None:
            continue
        for image in run_set:
            calc = Amp.load(amp_calc)
            true_engs.append(image.get_potential_energy())
            image.set_calculator(calc)
            amp_engs.append(image.get_potential_energy())
            print(true_engs[-1], amp_engs[-1])
            del calc
            #if per_atom:
            #    true_engs = [a / len(image) for a in true_engs]
            #    amp_engs = [a / len(image) for a in amp_engs]
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.plot(true_engs[:len(images)], amp_engs[:len(images)], 'ob',
            label='test set')
    if training_images is not None:
        ax.plot(true_engs[len(images):], amp_engs[len(images):],
                'ok', label='training set')
    print([a for a in zip(true_engs[:len(images)], amp_engs[:len(images)])])
    ax.plot([min(true_engs),max(true_engs)], [min(true_engs), max(true_engs)])
    ax.set_title('Parity Plot With AMP')
    ax.set_ylabel('AMP energy')
    ax.set_xlabel('Method Energy')
    ax.legend()
    plt.show()

def calc_rmse(array1, array2):
    sq_errs = []
    for i, j in zip(array1, array2):
        sq_errs.append((i - j) ** 2)
    return np.sqrt(np.mean(sq_errs))

def single_point_lammps(atoms, method='simple_nn_single_point',
                        ff_file='ffield.reax.water_2017', atoms_order=['H', 'O']):
    make_standard_input(calculation_type=method,
                        ff_file=ff_file,elements=atoms_order)
    write_lammps_data(atoms)
    os.system('lmp < lmp.input > lmp.log')
    atoms = parse_custom_dump('atoms.atm', 'lmp.data', label='atoms',
                              energyfile=None, write_traj=False,
                              units='real')
    return atoms



def run_schnetpack(db_file):
        import os
        import torch.nn.functional as F

        import logging
        from torch.optim import Adam
        import schnetpack as spk
        from schnetpack.train import Trainer, CSVHook, ReduceLROnPlateauHook, TensorboardHook
        from schnetpack.train.metrics import MeanAbsoluteError, RootMeanSquaredError
        #from schnetpack.metrics import build_mse_loss
        from schnetpack.train.metrics import MeanSquaredError
        from schnetpack.data import AtomsData
        import schnetpack.atomistic as atm
        import schnetpack.representation as rep
        import numpy as np

        #mse_loss = MeanSquaredError()
        mse_loss = spk.train.loss.build_mse_loss

        logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


        # basic settings
        model_name = db_file.split('.')[0]
        model_dir = model_name + '_model'  # directory that will be created for storing model
        if not os.path.isdir(model_name):
            os.mkdir(model_dir)
        #os.makedirs(model_dir)
        properties = ["energy", "forces"]  # properties used for training

        # data preparation
        logging.info("get dataset")
        dataset = AtomsData(db_file, 
                            available_properties=properties,
                            #required_properties=properties,
                            collect_triples=True)

        train, val, test = spk.train_test_split(
            data=dataset,
            num_train=400,
            num_val=20,
            split_file=os.path.join(model_dir, "split.npz"),
        )
        train_loader = spk.AtomsLoader(train, batch_size=32, shuffle=True)
        val_loader = spk.AtomsLoader(val, batch_size=5)

        # get statistics
        #atomrefs = dataset.get_atomrefs(properties)

        atomrefs={'energy':{'H':-4.578106,
                  'O':2052.084023}}
        energy_array = np.array([0,-6.879145831,0,0,0,0,0,0,-2055.665708] + [0] * 91)
        energy_array = energy_array.reshape((100,1))
        atomrefs={'energy':energy_array}
        per_atom = dict(energy=True, forces=False)
        means, stddevs = train_loader.get_statistics(
            ['energy'], 
            #atomrefs=atomrefs, 
            single_atom_ref=atomrefs,
            #per_atom=per_atom
            divide_by_atoms = True,
        )

        # model build
        logging.info("build model")
        #representation = spk.SchNet(n_interactions=6)
        #representation = spk.representation.SymmetryFunctions(elements={6})
        """
        output_modules = [
            spk.Atomwise(
                property="energy",
                derivative="forces",
                mean=means["energy"],
                stddev=stddevs["energy"],
                negative_dr=True,
            )
        ]
        """
        """
        reps = rep.BehlerSFBlock(elements={13},
                                 n_radial = 40,
                                 n_angular = 40,
                                 cutoff_radius = 8,
                                 crossterms = False,
                                 zetas = {1,4},
                                 
                                 )
        """
        reps = rep.SymmetryFunctions(
                n_radial=40,
                n_angular=40,
                zetas={1,4},
                cutoff=spk.nn.CosineCutoff,
                cutoff_radius=4.0,
                centered=False,
                crossterms=False,
                #elements=frozenset((1, 6, 7, 8, 9)),
                sharez=True,
                trainz=False,
                initz="weighted",
                len_embedding=5,
                pairwise_elements=False,)

        std_reps = rep.StandardizeSF(reps, train_loader)

        output = spk.atomistic.output_modules.ElementalAtomwise(reps.n_symfuncs,
                                                      n_layers=4,
                                                      atomref=atomrefs['energy'],
                                                      n_hidden=80,
                                                      property = 'energy',
                                                      derivative = 'forces',
                                                      mean=means["energy"],
                                                      stddev=stddevs["energy"],
                                                      negative_dr = True
                                                      )


        #output = spk.output_modules.Atomwise(reps.n_symfuncs)
        model = atm.AtomisticModel(std_reps, output)
        #model = spk.AtomisticModel(representation, output_modules)

        # build optimizer
        optimizer = Adam(params=model.parameters(), lr=1e-4)

        # hooks
        logging.info("build trainer")
        metrics = [MeanAbsoluteError(p, p) for p in properties] + \
                  [RootMeanSquaredError(p, p) for p in properties]
        hooks = [CSVHook(log_path=model_dir, metrics=metrics), 
                ReduceLROnPlateauHook(optimizer), 
                TensorboardHook(log_path=model_dir,
                metrics=metrics)]

        # trainer
        #loss = lambda b, p: F.mse_loss(p["y"], b['energy','force'])
        loss = mse_loss(properties)


        trainer = Trainer(
            model_dir,
            model=model,
            hooks=hooks,
            loss_fn=loss,
            optimizer=optimizer,
            train_loader=train_loader,
            validation_loader=val_loader,
            #device = 'cuda'
        )

        # run training

        logging.info("training")
        trainer.train(device="cuda", n_epochs=1000000)

