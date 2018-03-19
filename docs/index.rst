.. _index:

=======================
Introducing OMMProtocol
=======================

.. image:: https://travis-ci.org/insilichem/ommprotocol.svg?branch=master
    :target: https://travis-ci.org/insilichem/ommprotocol

.. image:: https://ci.appveyor.com/api/projects/status/3sobexd0dobfha09?svg=true
    :target: https://ci.appveyor.com/project/jaimergp/ommprotocol

.. image:: https://anaconda.org/insilichem/ommprotocol/badges/downloads.svg
    :target: https://anaconda.org/InsiliChem/ommprotocol

A command line application to launch molecular dynamics simulations with OpenMM


+ No coding required - just a YAML input file!
+ Smart support for different input file formats:
    + **Topology**: PDB/PDBx, Mol2, Amber's PRMTOP, Charmm's PSF, Gromacs' TOP, Desmond's DMS
    + **Positions**: PDB, COOR, INPCRD, CRD, GRO
    + **Velocities**: PDB, VEL
    + **Box vectors**: PDB, XSC, CSV, INPCRD, GRO
    + A fallback method is implemented and will attempt to load everything else that might be supported by `ParmEd <http://parmed.github.io/ParmEd/html/index.html>`_.
+ Choose your preferred **trajectory** format (PDB, PDBx, DCD, HDF5, NETCDF, MDCRD) and **checkpoints** (Amber's, CHARMM's, OpenMM XML states).
+ Live report of simulation progress, with estimated ETA and speed.
+ Checkpoint every *n* steps. Also, emergency rescue files are created if an error occurs.
+ Autochunk the trajectories for easy handling.


.. toctree::
    :maxdepth: 1
    :caption: For users

    install.rst
    usage.rst
    input.rst
    forcefields.rst

.. toctree::
    :maxdepth: 1
    :caption: For developers

    development.rst

.. toctree::
    :maxdepth: 1
    :caption: Support

    support.rst
