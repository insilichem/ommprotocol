#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

"""
ommprotocol.io
--------------

Handle IO stuff
"""
# Python stdlib
from __future__ import print_function, division, absolute_import
import os
import sys
from collections import namedtuple
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
# OpenMM and 3rd party helpers
try:
    import ruamel_yaml as yaml
except ImportError:
    import yaml
from simtk import unit as u
from simtk.openmm.app import (PDBFile, ForceField, AmberPrmtopFile, PDBReporter,
                              AmberInpcrdFile, CharmmPsfFile, CharmmParameterSet)
from simtk.openmm import XmlSerializer
from parmed.namd import NamdBinCoor, NamdBinVel
from parmed import load_file as parmed_load_file
from openmoltools.utils import create_ffxml_file
from mdtraj.reporters import DCDReporter
# 3rd party
from simtk.openmm import app
from parmed.openmm import RestartReporter, NetCDFReporter, MdcrdReporter
from mdtraj.reporters import HDF5Reporter
# Own
from .utils import sanitize_path_for_file
from .md import NONBONDEDMETHODS, CONSTRAINTS, FORCEFIELDS

if sys.version_info.major == 3:
    basestring = str

class YamlLoader(yaml.Loader):

    """
    YAML Loader with `!include` constructor. Straight from
    https://gist.github.com/joshbode/569627ced3076931b02f
    """

    def __init__(self, stream, version):
        """Initialise Loader."""

        try:
            self._root = os.path.dirname(stream.name)
        except AttributeError:
            self._root = os.path.curdir

        yaml.Loader.__init__(self, stream, version)

    def construct_include(self, node):
        """Include file referenced at node."""

        filename = os.path.join(self._root, self.construct_scalar(node))
        filename = os.path.abspath(filename)
        extension = os.path.splitext(filename)[1].lstrip('.')

        with open(filename, 'r') as f:
            if extension in ('yaml', 'yml'):
                return yaml.load(f, self)
            else:
                return ''.join(f.readlines())

YamlLoader.add_constructor('!include', YamlLoader.construct_include)


class MultiFormatLoader(object):

    """
    A base class to load different formats of the same type of file with a
    single method. It is meant to be inherited by handlers that have to deal
    with this situation.

    A basic `load` classmethod is provided to handle the delegation of the
    parsing to the proper loader. To do so, it depends on another classmethod
    `_loaders` that acts as a dict that maps file extensions to handler methods.
    """

    @classmethod
    def load(cls, path, *args, **kwargs):
        name, ext = os.path.splitext(path)
        try:
            return cls._loaders(ext.lstrip('.'))(path, *args, **kwargs)
        except KeyError:
            raise NotImplementedError('Unknown loader for format {}'.format(ext))
        except IOError:
            raise IOError('Could not access file {}'.format(path))

    @classmethod
    def _loaders(cls, ext):
        raise NotImplementedError('Override this method')


class InputContainer(object):

    """
    A base class to storage system parameters in an easy way, such as positions
    or velocities, with optional validation of input data
    """

    def __init__(self, topology=None, positions=None, velocities=None, box=None, **kwargs):
        self._topology = None
        self._positions = None
        self._velocities = None
        self._box = None
        self.topology = topology
        self.positions = positions
        self.velocities = velocities
        self.box = box

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, obj):
        self._topology = obj  # assertinstance(obj, (Topology, None))

    @property
    def positions(self):
        return self._positions

    @positions.setter
    def positions(self, obj):
        self._positions = obj  # assertinstance(obj, (u.Quantity, None))

    @property
    def velocities(self):
        return self._velocities

    @velocities.setter
    def velocities(self, obj):
        self._velocities = obj  # assertinstance(obj, (u.Quantity, None))

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, obj):
        self._box = obj  # assertinstance(obj, (u.Quantity, None))

    @property
    def has_topology(self):
        return self.topology is not None

    @property
    def has_positions(self):
        return self.positions is not None

    @property
    def has_velocities(self):
        return self.velocities is not None

    @property
    def has_box(self):
        return self.box is not None


class SystemHandler(MultiFormatLoader, InputContainer):

    """
    Loads an OpenMM topology from `path`

    Parameters
    ----------
    path : str
        Path to desired topology file. Supports pdb, prmtop, psf.
    """

    @classmethod
    def _loaders(cls, ext):
        return {'pdb': cls.from_pdb,
                'prmtop': cls.from_prmtop,
                'top': cls.from_prmtop,
                'psf': cls.from_psf}[ext]

    @classmethod
    def from_pdb(cls, path, forcefields=None, **kwargs):
        """
        Loads topology, positions and, potentially, velocities and vectors,
        from a PDB file

        Parameters
        ----------
        path : str
            Path to PDB file
        forcefields : list of str
            Paths to FFXML and/or FRCMOD forcefields. REQUIRED.

        Returns
        -------
        pdb : SystemHandler
            SystemHandler with topology, positions, and, potentially, velocities and
            box vectors. Forcefields are embedded in the `master` attribute.
        """
        pdb = PDBFile(path)
        box = kwargs.pop('box', pdb.topology.getPeriodicBoxVectors())
        positions = kwargs.pop('positions', pdb.positions)
        velocities = kwargs.pop('velocities', getattr(pdb, 'velocities', None))

        if not forcefields:
            forcefields = FORCEFIELDS
            print('INFO: Forcefields for PDB not specified. Using default:\n ',
                  ', '.join(forcefields))
        pdb.forcefield = ForceField(*process_forcefield(*forcefields))

        return cls(master=pdb, topology=pdb.topology, positions=positions,
                   velocities=velocities, box=box, path=path, **kwargs)

    @classmethod
    def from_prmtop(cls, path, positions=None, **kwargs):
        """
        Loads Amber Parm7 parameters and topology file

        Parameters
        ----------
        path : str
            Path to *.prmtop or *.top file
        positions : simtk.unit.Quantity
            Atomic positions

        Returns
        -------
        prmtop : SystemHandler
            SystemHandler with topology
        """
        if positions is None:
            raise ValueError('TOP/PRMTOP files require initial positions.')
        prmtop = AmberPrmtopFile(path)
        box = kwargs.pop('box', prmtop.topology.getPeriodicBoxVectors())
        return cls(master=prmtop, topology=prmtop.topology, positions=positions, box=box,
                   path=path, **kwargs)

    @classmethod
    def from_psf(cls, path, positions=None, charmm_parameters=None, **kwargs):
        """
        Loads PSF Charmm structure from `path`. Requires `charmm_parameters`.

        Parameters
        ----------
        path : str
            Path to PSF file
        charmm_parameters : list of str
            Paths to Charmm parameters files, such as *.par or *.str. REQUIRED

        Returns
        -------
        psf : SystemHandler
            SystemHandler with topology. Charmm parameters are embedded in
            the `master` attribute.
        """
        psf = CharmmPsfFile(path)
        if charmm_parameters is None:
            raise ValueError('PSF files require key `charmm_parameters`.')
        if positions is None:
            raise ValueError('PSF files require key `positions`.')
        psf.parmset = CharmmParameterSet(*charmm_parameters)
        psf.loadParameters(psf.parmset)
        return cls(master=psf, topology=psf.topology, positions=positions, path=path,
                   **kwargs)

    def __init__(self, master=None, **kwargs):
        InputContainer.__init__(self, **kwargs)
        if isinstance(master, str):
            raise ValueError('To instantiate from file, use .load() or '
                             'one of the .from_*() methods.')
        self.master = master
        self._path = kwargs.get('path')

    def create_system(self, **system_options):
        """
        Create an OpenMM system for every supported topology file with given system options
        """
        if self.master is None:
            raise ValueError('This instance is not able to create systems.')

        if isinstance(self.master, PDBFile):
            if not hasattr(self.master, 'forcefield'):
                raise ValueError('PDB topology files must be instanciated with forcefield paths.')
            system = self.master.forcefield.createSystem(self.topology, **system_options)

        elif isinstance(self.master, AmberPrmtopFile):
            system = self.master.createSystem(**system_options)

        elif isinstance(self.master, CharmmPsfFile):
            if not hasattr(self.master, 'parmset'):
                raise ValueError('PSF topology files must be instanciated with Charmm parameters.')
            system = self.master.createSystem(self.master.parmset, **system_options)
        
        else:
            raise NotImplementedError('This handler is not able to create systems.')

        if self.has_box:
            system.setDefaultPeriodicBoxVectors(*self.box)
        return system


class Positions(MultiFormatLoader):

    """
    Set of loaders to get position coordinates from file `path`

    Parameters
    ----------
    path : str
        Path to desired coordinates file. Supports pdb, coor, inpcrd.

    Returns
    -------
    positions : simtk.unit.Quantity([atoms,3])
    """

    @classmethod
    def _loaders(cls, ext):
        return {'pdb': cls.from_pdb,
                'coor': cls.from_coor,
                'inpcrd': cls.from_inpcrd,
                'crd': cls.from_inpcrd}[ext]

    @classmethod
    def from_pdb(cls, path):
        pdb = PDBFile(path)
        return pdb.positions

    @classmethod
    def from_coor(cls, path):
        coor = NamdBinCoor.read(path)
        positions = u.Quantity(coor.coordinates[0], unit=u.angstroms)
        return positions

    @classmethod
    def from_inpcrd(cls, path):
        inpcrd = AmberInpcrdFile(path)
        return inpcrd.positions


class Velocities(MultiFormatLoader):

    """
    Set of loaders to get velocities for a given topology from file `path`

    Parameters
    ----------
    path : str
        Path to desired velocities file. Supports vel.

    Returns
    -------
    velocities : simtk.unit.Quantity([atoms,3])
    """

    @classmethod
    def _loaders(cls, ext):
        return {'vel': cls.from_vel}[ext]

    @classmethod
    def from_vel(cls, path):
        vel = NamdBinVel.read(path)
        velocities = u.Quantity(vel.velocities[0], unit=u.angstroms/u.picosecond)
        return velocities


class BoxVectors(MultiFormatLoader):

    """
    Set of loaders to get vectors from file `path`

    Parameters
    ----------
    path : str
        Path to desired velocities file. Supports vel.

    Returns
    -------
    velocities : simtk.unit.Quantity([atoms,3])
    """

    @classmethod
    def _loaders(cls, ext):
        return {'xsc': cls.from_xsc}[ext]

    @classmethod
    def from_xsc(cls, path):
        """ Returns u.Quantity with box vectors from XSC file """

        def parse(path):
            """
            Open and parses an XSC file into its fields

            Parameters
            ----------
            path : str
                Path to XSC file

            Returns
            -------
            namedxsc : namedtuple
                A namedtuple with XSC fields as names
            """
            with open(path) as f:
                lines = f.readlines()
            NamedXsc = namedtuple('NamedXsc', lines[1].split()[1:])
            return NamedXsc(*map(float, lines[2].split()))

        xsc = parse(path)
        return u.Quantity([[xsc.a_x, xsc.a_y, xsc.a_z], 
                           [xsc.b_x, xsc.b_y, xsc.b_z],
                           [xsc.c_x, xsc.c_y, xsc.c_z]], unit=u.angstroms)

    @classmethod
    def from_csv(cls, path):
        """
        Get box vectors from comma-separated values in file `path`.

        The csv file must containt only one line, which in turn can contain
        three values (orthogonal vectors) or nine values (triclinic box).

        The values should be in nanometers.

        Parameters
        ----------
        path : str
            Path to CSV file

        Returns
        -------
        vectors : simtk.unit.Quantity([3, 3], unit=nanometers
        """
        with open(path) as f:
            fields = map(float, next(f).split(','))
        if len(fields) == 3:
            return u.Quantity([[fields[0], 0, 0], 
                               [0, fields[1], 0],
                               [0, 0, fields[2]]], unit=u.nanometers)
        if len(fields) == 9:
            return u.Quantity([fields[0:3], 
                               fields[3:6],
                               fields[6:9]], unit=u.nanometers)


class Restart(MultiFormatLoader, InputContainer):

    """
    Loads a restart file that can contain positions, velocities and box vectors.

    Parameters
    ----------
    path : str
        Restart file

    Returns
    -------
    positions : simtk.unit.Quantity([atoms,3])
    velocities : simtk.unit.Quantity([atoms,3])
    vectors : simtk.unit.Quantity([1,3])
    """

    @classmethod
    def from_xml(cls, path):
        with open(path) as f:
            xml = XmlSerializer.deserialize(f.read())
        positions = xml.getPositions()
        velocities = xml.getVelocities()
        box = xml.getPeriodicBoxVectors()
        return cls(positions=positions, velocities=velocities, box=box)

    @classmethod
    def from_rst(cls, path):
        positions, velocities, box = None, None, None
        rst = parmed_load_file(path)
        positions = u.Quantity(rst.coordinates[0], unit=u.angstrom)
        if rst.hasvels:
            velocities = u.Quantity(rst.velocities[0], unit=u.angstrom/u.picosecond)
        if rst.hasbox:
            vectors = [[rst.cell_lengths[0], 0, 0],
                       [0, rst.cell_lengths[1], 0],
                       [0, 0, rst.cell_lengths[2]]]
            box = u.Quantity(vectors, unit=u.angstrom)
        return cls(positions=positions, velocities=velocities, box=box)

    @classmethod
    def _loaders(cls, ext):
        return {'xml': cls.from_xml,
                'xmlstate': cls.from_xml,
                'state': cls.from_xml,
                'rs': cls.from_rst,
                'rst': cls.from_rst,
                'restart': cls.from_rst}[ext]


class ProgressBarReporter(object):
    """
    A simple progress bar reporter for stdout.
    """
    def __init__(self, file, interval, total_steps=None, margin=4):
        """Create a ProgressBarReporter.

        Parameters
        ----------
        file : string or open file object
            The file to write to. Normally stdout. Will overwrite!
        interval : int
            The interval (in time steps) at which to write checkpoints.
        margin : int
            Blank padding to write before bar.
        """
        if isinstance(file, str):
            self._own_handle = True
            self._out = open(file, 'w', 0)
        else:
            self._out = file
            self._own_handle = False

        self.interval = interval
        self.margin = margin
        self.total_steps = total_steps
        self._initialized = False
        

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self.interval - simulation.currentStep%self.interval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._initialized:
            self._initial_clock_time  = datetime.now()
            self._initial_simulation_time = state.getTime()
            self._initial_steps = simulation.currentStep
            self._initialized = True

        steps = simulation.currentStep
        time = datetime.now() - self._initial_clock_time
        days = time.total_seconds()/86400.0
        ns = (state.getTime()-self._initial_simulation_time).value_in_unit(u.nanosecond)

        margin = ' ' * self.margin
        ns_day = ns/days
        delta = (self.total_steps-steps)*time/steps
        # remove microseconds to have cleaner output
        remaining = timedelta(days=delta.days, seconds=delta.seconds)
        percentage = 100.0*steps/self.total_steps
        if ns_day:
            template = '{}{}/{} steps ({:.1f}%) - {} left @ {:.1f} ns/day                    \r'
        else:
            template = '{}{}/{} steps ({:.1f}%)                                              \r'
        report = template.format(margin, steps, self.total_steps, percentage, remaining, ns_day)
        self._out.write(report)
        if hasattr(self._out, 'flush'):
            self._out.flush()

    def __del__(self):
        self._out.write('\n')
        if self._own_handle:
            self._out.close()

class SegmentedDCDReporter(DCDReporter):

    """
    A DCD reporter based on mdtraj's that creates a new file every n `new_every`
    """

    def __init__(self, file, reportInterval, atomSubset=None, new_every=0):
        self._original_filename = file if isinstance(file, basestring) else None
        self.new_every = new_every
        super(SegmentedDCDReporter, self).__init__(file, reportInterval, atomSubset=atomSubset)

    def report(self, simulation, state):
        if self.new_every > 0:
            self._check_size(simulation.currentStep)
        super(SegmentedDCDReporter, self).report(simulation, state)

    def _check_size(self, current_step):
        if current_step and not current_step % self.new_every and self._original_filename:
            path, ext = os.path.splitext(self._original_filename)
            filename = '{}.{}{}'.format(path, current_step, ext)
            self._traj_file = self.backend(filename, 'w')


###########################
# Input preparation
###########################

def prepare_input(argv=None):
    """
    Get, parse and prepare input file.
    """
    p = ArgumentParser(description='insiliChem.bio OpenMM launcher: '
                       'easy to deploy MD protocols for OpenMM')
    p.add_argument('input', metavar='INPUT FILE', type=str,
                   help='YAML input file')
    args = p.parse_args(argv if argv else sys.argv[1:])

    # Load config file
    with open(args.input) as f:
        cfg = yaml.load(f, YamlLoader)
    # Paths and dirs
    cfg['_path'] = os.path.abspath(args.input)
    cfg['system_options'] = prepare_system_options(cfg)
    cfg['outputpath'] = sanitize_path_for_file(cfg.get('outputpath', '.'), args.input)
    try: 
        os.makedirs(cfg['outputpath'])
    except OSError:
        pass

    handler = prepare_handler(cfg)

    return handler, cfg


def prepare_handler(cfg):
    """
    Load all files into single object.
    """
    positions, velocities, boxvectors = None, None, None

    if 'checkpoint' in cfg:
        restart_path = sanitize_path_for_file(cfg['checkpoint'], cfg['_path'])
        restart = Restart.load(restart_path)
        positions = restart.positions
        velocities = restart.velocities
        boxvectors = restart.box

    if 'positions' in cfg:
        positions_path = sanitize_path_for_file(cfg.pop('positions'), cfg['_path'])
        positions = Positions.load(positions_path)

    if 'velocities' in cfg:
        velocities_path = sanitize_path_for_file(cfg.pop('velocities'), cfg['_path'])
        velocities = Velocities.load(velocities_path)

    if 'box' in cfg:
        boxvectors_path = sanitize_path_for_file(cfg.pop('box', None), cfg['_path'])
        boxvectors = BoxVectors.load(boxvectors_path)

    return SystemHandler.load(cfg.pop('topology'), positions=positions, velocities=velocities,
                              box=boxvectors, forcefield=cfg.pop('forcefield', None),
                              charmm_parameters=cfg.pop('charmm_parameters', None))


def prepare_system_options(cfg, fill_not_found=True):
    """
    Retrieve and delete (pop) system options from input configuration.
    """
    d = {}
    if fill_not_found:
        d['nonbondedMethod'] = NONBONDEDMETHODS.get(cfg.pop('nonbondedMethod', None))
        d['nonbondedCutoff'] = cfg.pop('nonbondedCutoff', None) * u.nanometers
        d['constraints'] = CONSTRAINTS.get(cfg.pop('constraints', None))
        for key in ['rigidWater', 'ewaldErrorTolerance']:
            d[key] = cfg.pop(key, None)
    else:
        if 'nonbondedMethod' in cfg:
            d['nonbondedMethod'] = NONBONDEDMETHODS.get(cfg.pop('nonbondedMethod'))
        if 'nonbondedCutoff' in cfg:
            d['nonbondedCutoff'] = cfg.pop('nonbondedCutoff') * u.nanometers
        if 'constraints' in cfg:
            d['constraints'] = CONSTRAINTS.get(cfg.pop('constraints'))
        for key in ['rigidWater', 'ewaldErrorTolerance']:
            if key in 'cfg':
                d[key] = cfg.pop(key)
    return d


###########################
# Helpers
###########################

def process_forcefield(*forcefields):
    """
    Given a list of filenames, check which ones are `frcmods`. If so,
    convert them to ffxml. Else, just return them.
    """
    for forcefield in forcefields:
        if forcefield.endswith('.frcmod'):
            yield create_ffxml_file(forcefield)
        else:
            yield forcefield

###########################
# Defaults
###########################

REPORTERS = {
    'PDB': PDBReporter,
    'DCD': SegmentedDCDReporter,
    'CHK': app.CheckpointReporter,
    'HDF5':  HDF5Reporter,
    'NETCDF': NetCDFReporter,
    'MDCRD': MdcrdReporter,
    'RS': partial(RestartReporter, write_multiple=True, netcdf=True),
}