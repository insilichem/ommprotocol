.. _forcefields:

==================
OpenMM forcefields
==================

Since OMMProtocol is compatible with a multiple formats thanks to OpenMM itself and other excellent packages (MDTraj, ParmEd...), you can use the forcefield formats defined in other MD suites. Namely, PRMTOP for Amber, PSF+PAR+STR in CHARMM or TOP in Gromacs. If you are already using some of those packages, you don't need to do anything else: just provide the paths in the topology section (CHARMM parameters must be specified in ``charmm_parameters``).

However, OpenMM does provide its own set of forcefields, converted from the original formats (Amber, Charmm and others) to its FFXML format. The following section lists all the built in forcefields in OpenMM as of v7.2. The updated list will be available at the `OpenMM repo <https://github.com/pandegroup/openmm/tree/master/wrappers/python/simtk/openmm/app/data>`_.

::

    amber14/DNA.OL15.xml
    amber14/DNA.bsc1.xml
    amber14/RNA.OL3.xml
    amber14/lipid17.xml
    amber14/protein.ff14SB.xml
    amber14/protein.ff15ipq.xml
    amber14/spce.xml
    amber14/tip3p.xml
    amber14/tip3pfb.xml
    amber14/tip4pew.xml
    amber14/tip4pfb.xml
    charmm36/spce.xml
    charmm36/tip3p-pme-b.xml
    charmm36/tip3p-pme-f.xml
    charmm36/tip4p2005.xml
    charmm36/tip4pew.xml
    charmm36/tip5p.xml
    charmm36/tip5pew.xml
    charmm36/water.xml
    absinth.xml
    amber03.xml
    amber03_obc.xml
    amber10.xml
    amber10_obc.xml
    amber14-all
    amber96.xml
    amber96_obc.xml
    amber99Test.xml
    amber99_obc.xml
    amber99sb.xml
    amber99sbildn.xml
    amber99sbnmr.xml
    amberfb15.xml
    amoeba2009.xml
    amoeba2009_gk.xml
    amoeba2013.xml
    amoeba2013_gk.xml
    charmm36.xml
    charmm_polar_2013.xml
    hydrogens.xml
    iamoeba.xml
    pdbNames.xml
    residues.xml
    spce.xml
    swm4ndp.xml
    tip3p.xml
    tip3pfb.xml
    tip4pew.xml
    tip4pfb.xml
    tip5p.xml


To use them with a PDB file, just specify them in a list for the ``forcefield`` key, like:

::

    topology: some.pdb
    forcefield: [amber99sbildn.xml, tip3p.xml]

or, if you prefer this other syntax:

::

    topology: some.pdb
    forcefield:
    - amber99sbildn.xml
    - tip3p.xml


More forcefields
----------------

The OpenMM team is doing a tremendous effort towards the next release, which will include even more forcefields. You can check the progress `here <https://github.com/choderalab/openmm-forcefields/projects/1>`_. This will include more builtin forcefields and also a separate package called ``openmm-forcefields``, developed `here <https://github.com/choderalab/openmm-forcefields>`_. When this is available, it will be shipped with OMMProtocol.


Custom forcefields
------------------

While the best option to generate custom parameters is to use something like AmberTools to create a PRMTOP topology and use that, there are options to develop custom parameters with OpenMM. Check these links for further information:

- `Creating and Customizing Force Fields in OpenMM <https://www.youtube.com/watch?v=xap418xVjNI>`_ (YouTube video).
- ``openmm-forcefields`` also features Python converters for Amber & CHARMM forcefields. As a result, automated tooling for those forcefields can be used and then converted to OpenMM, like ``antechamber`` or `cgenff``.
- `openmoltools <https://github.com/choderalab/openmoltools>`_ (included with OMMProtocol) provides some functions to process and convert forcefields. Specifically, ``openmoltools.amber.run_antechamber`` for parameterizing small molecules through AmberTools' ``antechamber``, and ``openmoltools.utils.create:_ffxml_file`` to convert the result to OpenMM XML forcefield format.