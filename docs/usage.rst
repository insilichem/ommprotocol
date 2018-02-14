======================
How to use OMMProtocol
======================

Once `installed <install.rst>`_, first thing to do is `creating the input file <input.rst>`_ for your simulation. This task usually involve two different steps:

1. Getting the structural data (topology, coordinates)
2. Specifying the simulation details:
    - Forcefield parameters
    - Solvation conditions: explicit, implicit, no solvent?
    - Simulation conditions: temperature, pressure, NPT, NVT?
    - Technical details: how to compute non-bonded interactions, whether to use periodic boundary conditions, whether to constrain some specific types of bonds, the integration method and timestep...

These details are probably out of the scope of this documentation, and the reader is encouraged to read specific tutorials about this, such as:

- `An Introduction to Molecular Dynamics Simulations using AMBER <http://ambermd.org/tutorials/basic/tutorial0/index.htm>`_
- `NAMD tutorials <http://www.ks.uiuc.edu/Training/Tutorials/namd-index.html>`_
- `GROMACS tutorials <http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/>`_

With a correctly formed YAML input file named, for example, _simulation.yaml_, the user can now run:

::

    ommprotocol simulation.yaml

If the structure is correctly formed and the forcefield parameters are well defined, the screen will now display a status like this:

__GIF animation pending__


The generated files will be written to the directory specified in the ``outputpath`` key (or, if omitted, to the same directory _simulation.yaml_ is in), with the following name format: ``[globalname]_[stagename].[extension]``, where ``globalname`` is the value of the global ``name`` key in the input file, and ``stagename`` is the value of the ``stage`` key in each stage.

Most of this files can be opened during the simulation. That way you can check the progress of the trajectory in viewers like VMD, PyMol or UCSF Chimera. Since .log files are created by default with some metadata about the simulation (temperature, potential energy, volume...), they are a convenient way of checking if everything is working OK. For example, the energies and temperatures should be more or less constant. To do that, a helper utility called ``ommanalyze`` is included, which is able to produce interactive plots of such properties:

.. image:: img/ommanalyze.png