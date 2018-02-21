===========
OMMProtocol
===========

.. image:: https://travis-ci.org/insilichem/ommprotocol.svg?branch=master
    :target: https://travis-ci.org/insilichem/ommprotocol

.. image:: https://ci.appveyor.com/api/projects/status/3sobexd0dobfha09?svg=true
    :target: https://ci.appveyor.com/project/jaimergp/ommprotocol

.. image:: https://anaconda.org/insilichem/ommprotocol/badges/downloads.svg
    :target: https://anaconda.org/InsiliChem/ommprotocol

.. image:: https://readthedocs.org/projects/ommprotocol/badge/?version=latest
    :target: http://ommprotocol.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


A command line application to launch molecular dynamics simulations with OpenMM

.. image:: https://raw.githubusercontent.com/insilichem/ommprotocol/master/docs/img/ommprotocol.gif


Some cool features
------------------

+ No coding required - just a YAML input file!
+ Smart support for different input file formats:
    + **Topology**: PDB/PDBx, Mol2, Amber's PRMTOP, Charmm's PSF, Gromacs' TOP, Desmond's DMS
    + **Positions**: PDB, COOR, INPCRD, CRD, GRO
    + **Velocities**: PDB, VEL
    + **Box vectors**: PDB, XSC, CSV, INPCRD, GRO
    + A fallback method is implemented and will attempt to load verything else that might be supported by `ParmEd <http://parmed.github.io/ParmEd/html/index.html>`_.
+ Choose your preferred **trajectory** format (PDB, PDBx, DCD, HDF5, NETCDF, MDCRD) and **checkpoints** (Amber's, CHARMM's, OpenMM XML states).
+ Live report of simulation progress, with estimated ETA and speed.
+ Checkpoint every *n* steps. Also, emergency rescue files are created if an error occurs.
+ Autochunk the trajectories for easy handling.


Installation & usage
--------------------
Download the `latest installer <https://github.com/insilichem/ommprotocol/releases/latest>`_ or use ``conda install -c omnia -c insilichem ommprotocol`` if you already have Anaconda/Miniconda installed. Further details `here <http://ommprotocol.readthedocs.io/en/latest/install.html>`_.

When installed, you should be able to run:

::

    ommprotocol <inputfile.yaml>

Check the `documentation <http://ommprotocol.readthedocs.io/en/latest/input.html>`_ to read more on how to create input files for OMMProtocol.


Get help
--------

.. image:: https://readthedocs.org/projects/ommprotocol/badge/?version=latest
    :target: http://ommprotocol.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Documentation is always available at `ReadTheDocs <http://ommprotocol.readthedocs.io>`_. If you have problems running ``ommprotocol``, feel free to `create an issue <https://github.com/insilichem/ommprotocol/issues>`_! Also, make sure to visit our main webpage at `insilichem.com <http://www.insilichem.com>`_.


Citation
--------

OMMProtocol is scientific software, funded by public research grants (Spanish MINECO's project ``CTQ2014-54071-P``, Generalitat de Catalunya's project ``2014SGR989`` and research grant ``2015FI_B00768``, COST Action ``CM1306``). If you make use of Ommprotocol in scientific publications, please cite it. It will help measure the impact of our research and future funding! A manuscript is in progress. In the meantime, please cite this repository URL.

.. code-block:: latex

    @misc{ommprotocol,
    author       = {Jaime Rodríguez-Guerra Pedregal and
                    Lur Alonso-Cotchico and
                    Lorea Velasco and
                    Jean-Didier Maréchal},
    title        = {OMMProtocol: A command line application to launch molecular dynamics simulations with OpenMM},
    url          = {https://github.com/insilichem/ommprotocol}
    }
