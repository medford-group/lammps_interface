#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(name='lammps_interface',
      version='0.1',
      description='Simple Python Tools to Get LAMMPS Working',
      author='Ben Comer',
      author_email='ben.comer@gatech.edu',
      url='https://github.com/medford-group/lammps_interface',
      packages=find_packages(),
      install_requires=['spglib', 'numpy','ase','scipy',
                        'MPInterfaces_Latest','pymatgen',
                        'moltemplate'],
     )

