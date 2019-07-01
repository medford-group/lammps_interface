# Medford Group LAMMPS Interface Tools
This is a collection of python tools used to interface with [LAMMPS](https://lammps.sandia.gov/). The interface is Linux only for the moment, if you're interested in Widnows support, open an issue and I can have a look at it. This package relies on `packmol`, please ensure that the [packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml) executable is in your PATH somewhere. 

This package utilizes the Atomistic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/index.html)) fairly extensively, but does not use the calculator system within it. This enables the package to be very modular.

## Installation
Install with:
`pip install git+https://github.com/medford-group/lammps_interface`

ensure that packmol is installed for full functionality, but some stuff works without it

## Usage
Since this package consists of python functions you'll need to import them from the `tools.py` file. The functions have long names so that they are as descriptive and simple as possible.

```
from lammps_interface.lammps_interface.tools import make_box_of_molecules
from ase.build import molecule

atoms = make_box_of_molecules(molecule[H2O],700,[40,40,40])
```

There are lots of goodies in there that you may find valuble.

If you have feature requests, open an issue and I can have a look.
