#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

"""
ommprotocol.io
--------------

IO stuff:

- Parser logic for YAML input files
- Handlers for source topologies, coordinates, velocities,
box vectors and restart checkpoints, and all the business
logic to deal with the precedence of those files.
- Output reporters for stdout and segmented DCD trajectories.
"""

# Python stdlib
from __future__ import print_function, division, absolute_import
import os
import pickle
import sys
from collections import namedtuple
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
import logging
# OpenMM and 3rd party helpers
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml
import jinja2
from simtk import unit as u
from simtk.openmm.app import (PDBFile, PDBxFile, ForceField,
                              PDBReporter, PDBxReporter,
                              AmberPrmtopFile, AmberInpcrdFile,
                              CharmmPsfFile, CharmmCrdFile, CharmmParameterSet,
                              GromacsTopFile, GromacsGroFile,
                              DesmondDMSFile, CheckpointReporter)
from simtk.openmm import XmlSerializer, app as openmm_app
import mdtraj
from mdtraj.reporters import DCDReporter, HDF5Reporter
import parmed
from parmed.namd import NamdBinCoor, NamdBinVel
from parmed.openmm import RestartReporter, NetCDFReporter, MdcrdReporter
from openmoltools.utils import create_ffxml_file
# Own
from .utils import sanitize_path_for_file, ignored_exceptions, warned_getattr, extant_file
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

if sys.version_info.major == 3:
    basestring = str

logger = logging.getLogger(__name__)


class YamlLoader(yaml.Loader):

    """
    YAML Loader with `!include` constructor. Straight from
    https://gist.github.com/joshbode/569627ced3076931b02f
    """

    def __init__(self, stream, version=None, **kwargs):
        """Initialise Loader."""

        try:
            self._root = os.path.dirname(stream.name)
        except AttributeError:
            self._root = os.path.curdir

        yaml.Loader.__init__(self, stream, version=None, **kwargs)

    def construct_include(self, node):
        """Include file referenced at node."""

        filename = os.path.join(self._root, self.construct_scalar(node))
        filename = os.path.abspath(filename)
        extension = os.path.splitext(filename)[1].lstrip('.')

        with open(filename, 'r') as f:
            if extension in ('yaml', 'yml'):
                return yaml.load(f, Loader=self)
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
            logger.error('! Unknown loader for format %s. '
                         'Trying with ParmEd as fallback', ext)
            return cls.from_parmed(path, *args, **kwargs)
        except IOError:
            raise IOError('Could not access file {}'.format(path))

    @classmethod
    def _loaders(cls, ext):
        raise NotImplementedError('Override this method')

    @classmethod
    def from_parmed(cls, *args, **kwargs):
        raise NotImplementedError('ParmEd fallback strategy not available here.')

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
        Path to desired topology file. Supports pdb, prmtop, psf, top, dms,
        pickled OpenMM Topology objects.
    """

    @classmethod
    def _loaders(cls, ext):
        return {'pdb': cls.from_pdb,
                'pdbx': cls.from_pdbx,
                'cif': cls.from_pdbx,
                'prmtop': cls.from_amber,
                'psf': cls.from_charmm,
                'dms': cls.from_desmond,
                'top': cls.from_gromacs,
                'pickle': cls.from_pickle,
                'pickle2': cls.from_pickle,
                'pickle3': cls.from_pickle}[ext]

    @classmethod
    def from_pdb(cls, path, forcefield=None, loader=PDBFile, **kwargs):
        """
        Loads topology, positions and, potentially, velocities and vectors,
        from a PDB or PDBx file

        Parameters
        ----------
        path : str
            Path to PDB/PDBx file
        forcefields : list of str
            Paths to FFXML and/or FRCMOD forcefields. REQUIRED.

        Returns
        -------
        pdb : SystemHandler
            SystemHandler with topology, positions, and, potentially, velocities and
            box vectors. Forcefields are embedded in the `master` attribute.
        """
        pdb = loader(path)
        box = kwargs.pop('box', pdb.topology.getPeriodicBoxVectors())
        positions = kwargs.pop('positions', pdb.positions)
        velocities = kwargs.pop('velocities', getattr(pdb, 'velocities', None))

        if not forcefield:
            from .md import FORCEFIELDS as forcefield
            logger.info('! Forcefields for PDB not specified. Using default: %s',
                        ', '.join(forcefield))
        pdb.forcefield = ForceField(*list(process_forcefield(*forcefield)))

        return cls(master=pdb.forcefield, topology=pdb.topology, positions=positions,
                   velocities=velocities, box=box, path=path, **kwargs)

    @classmethod
    def from_pdbx(cls, *args, **kwargs):
        __doc__ = cls.from_pdb.__doc__
        return cls.from_pdb(loader=PDBxFile, *args, **kwargs)

    @classmethod
    def from_amber(cls, path, positions=None, **kwargs):
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
            raise ValueError('Amber TOP/PRMTOP files require initial positions.')
        prmtop = AmberPrmtopFile(path)
        box = kwargs.pop('box', prmtop.topology.getPeriodicBoxVectors())
        return cls(master=prmtop, topology=prmtop.topology, positions=positions, box=box,
                   path=path, **kwargs)

    @classmethod
    def from_charmm(cls, path, positions=None, forcefield=None, **kwargs):
        """
        Loads PSF Charmm structure from `path`. Requires `charmm_parameters`.

        Parameters
        ----------
        path : str
            Path to PSF file
        forcefield : list of str
            Paths to Charmm parameters files, such as *.par or *.str. REQUIRED

        Returns
        -------
        psf : SystemHandler
            SystemHandler with topology. Charmm parameters are embedded in
            the `master` attribute.
        """
        psf = CharmmPsfFile(path)
        if forcefield is None:
            raise ValueError('PSF files require key `forcefield`.')
        if positions is None:
            raise ValueError('PSF files require key `positions`.')
        psf.parmset = CharmmParameterSet(*forcefield)
        psf.loadParameters(psf.parmset)
        return cls(master=psf, topology=psf.topology, positions=positions, path=path,
                   **kwargs)

    @classmethod
    def from_desmond(cls, path, **kwargs):
        """
        Loads a topology from a Desmond DMS file located at `path`.

        Arguments
        ---------
        path : str
            Path to a Desmond DMS file
        """
        dms = DesmondDMSFile(path)
        pos = kwargs.pop('positions', dms.getPositions())
        return cls(master=dms, topology=dms.getTopology(), positions=pos, path=path,
                   **kwargs)

    @classmethod
    def from_gromacs(cls, path, positions=None, forcefield=None, **kwargs):
        """
        Loads a topology from a Gromacs TOP file located at `path`.

        Additional root directory for parameters can be specified with `forcefield`.

        Arguments
        ---------
        path : str
            Path to a Gromacs TOP file
        positions : simtk.unit.Quantity
            Atomic positions
        forcefield : str, optional
            Root directory for parameter files
        """
        if positions is None:
            raise ValueError('Gromacs TOP files require initial positions.')
        box = kwargs.pop('box', None)
        top = GromacsTopFile(path, includeDir=forcefield, periodicBoxVectors=box)
        return cls(master=top, topology=top.topology, positions=positions, box=box,
                   path=path, **kwargs)

    @classmethod
    def from_parmed(cls, path, *args, **kwargs):
        """
        Try to load a file automatically with ParmEd. Not guaranteed to work, but
        might be useful if it succeeds.

        Arguments
        ---------
        path : str
            Path to file that ParmEd can load
        """
        st = parmed.load_file(path, structure=True, *args, **kwargs)
        box = kwargs.pop('box', getattr(st, 'box', None))
        velocities = kwargs.pop('velocities', getattr(st, 'velocities', None))
        positions = kwargs.pop('positions', getattr(st, 'positions', None))
        return cls(master=st, topology=st.topology, positions=positions, box=box,
                   velocities=velocities, path=path, **kwargs)

    @classmethod
    def from_pickle(cls, path, positions=None, forcefield=None, **kwargs):
        if positions is None:
            raise ValueError('Pickled topology files require initial positions.')
        if forcefield is None:
            raise ValueError('Pickled topology files require XML/FRCMOD forcefields.')

        topology = cls._pickle_load(path)
        forcefield = ForceField(*list(process_forcefield(*forcefield)))
        return cls(master=forcefield, topology=topology, positions=positions, path=path,
                   **kwargs)

    @staticmethod
    def _pickle_load(path):
        """
        Loads pickled topology. Careful with Python versions though!
        """
        _, ext = os.path.splitext(path)
        topology = None
        if sys.version_info.major == 2:
            if ext == '.pickle2':
                with open(path, 'rb') as f:
                    topology = pickle.load(f)
            elif ext in ('.pickle3', '.pickle'):
                with open(path, 'rb') as f:
                    topology = pickle.load(f, protocol=3)
        elif sys.version_info.major == 3:
            if ext == '.pickle2':
                with open(path, 'rb') as f:
                    topology = pickle.load(f)
            elif ext in ('.pickle3', '.pickle'):
                with open(path, 'rb') as f:
                    topology = pickle.load(f)
        if topology is None:
            raise ValueError('File {} is not compatible with this version'.format(path))
        return topology

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
            raise ValueError('Handler {} is not able to create systems.'.format(self))

        if isinstance(self.master, ForceField):
            system = self.master.createSystem(self.topology, **system_options)
        elif isinstance(self.master, (AmberPrmtopFile, GromacsTopFile, DesmondDMSFile)):
            system = self.master.createSystem(**system_options)
        elif isinstance(self.master, CharmmPsfFile):
            if not hasattr(self.master, 'parmset'):
                raise ValueError('PSF topology files require Charmm parameters.')
            system = self.master.createSystem(self.master.parmset, **system_options)
        else:
            raise NotImplementedError('Handler {} is not able to create systems.'.format(self))

        if self.has_box:
            system.setDefaultPeriodicBoxVectors(*self.box)
        return system

    def write_pdb(self, path):
        """
        Outputs a PDB file with the current contents of the system
        """
        if self.master is None and self.positions is None:
            raise ValueError('Topology and positions are needed to write output files.')
        with open(path, 'w') as f:
            PDBFile.writeFile(self.topology, self.positions, f)


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
                'coor': cls.from_namd,
                'gro': cls.from_gromacs,
                'inpcrd': cls.from_amber,
                'crd': cls.from_charmm}[ext]

    @classmethod
    def from_pdb(cls, path):
        return PDBFile(path).positions

    @classmethod
    def from_namd(cls, path):
        coor = NamdBinCoor.read(path)
        positions = u.Quantity(coor.coordinates[0], unit=u.angstroms)
        return positions

    @classmethod
    def from_amber(cls, path):
        return AmberInpcrdFile(path).positions

    @classmethod
    def from_gromacs(cls, path):
        return GromacsGroFile(path).getPositions()

    @classmethod
    def from_charmm(cls, path):
        return CharmmCrdFile(path).positions

    @classmethod
    def from_parmed(cls, path):
        return parmed.load_file(path, structure=True).positions


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
        return {'vel': cls.from_namd}[ext]

    @classmethod
    def from_namd(cls, path):
        vel = NamdBinVel.read(path)
        velocities = u.Quantity(vel.velocities[0], unit=u.angstroms/u.picosecond)
        return velocities

    @classmethod
    def from_parmed(cls, path):
        return parmed.load_file(path, structure=True).velocities


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
        return {'xsc': cls.from_xsc,
                'csv': cls.from_csv,
                'pdb': cls.from_pdb,
                'gro': cls.from_gromacs,
                'inpcrd': cls.from_amber}[ext]

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
        elif len(fields) == 9:
            return u.Quantity([fields[0:3],
                               fields[3:6],
                               fields[6:9]], unit=u.nanometers)
        else:
            raise ValueError('This type of CSV is not supported. Please '
                             'provide a comma-separated list of three or nine '
                             'floats in a single-line file.')

    @classmethod
    def from_pdb(cls, path):
        return PDBFile(path).topology.getPeriodicBoxVectors()

    @classmethod
    def from_amber(cls, path):
        return AmberInpcrdFile(path).boxVectors

    @classmethod
    def from_gromacs(cls, path):
        return GromacsGroFile(path).getPeriodicBoxVectors()

    @classmethod
    def from_parmed(cls, path):
        return parmed.load_file(path, structure=True, hasbox=True).box


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
    def _loaders(cls, ext):
        return {'xml': cls.from_xml,
                'xmlstate': cls.from_xml,
                'state': cls.from_xml,
                'rs': cls.from_rst,
                'rst': cls.from_rst,
                'restart': cls.from_rst}[ext]

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
        rst = parmed.load_file(path)
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
    def from_parmed(cls, path):
        st = parmed.load_file(path)
        return cls(positions=st.positions, velocities=st.velocities, box=st.box)


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
        steps = self.interval - simulation.currentStep % self.interval
        return steps, False, False, False, False

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
            self._initial_clock_time = datetime.now()
            self._initial_simulation_time = state.getTime()
            self._initial_steps = simulation.currentStep
            self._initialized = True

        steps = simulation.currentStep
        time = datetime.now() - self._initial_clock_time
        days = time.total_seconds()/86400.0
        ns = (state.getTime()-self._initial_simulation_time).value_in_unit(u.nanosecond)

        margin = ' ' * self.margin
        ns_day = ns/days
        delta = ((self.total_steps-steps)*time.total_seconds())/steps
        # remove microseconds to have cleaner output
        remaining = timedelta(seconds=int(delta))
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


class SerializedReporter(object):

    """
    Report progress and state in a serialized format,
    enabling inter-process reporting.

    Parameters
    ----------
    file : string or open file object
        The file to write to. Normally stdout. Will overwrite!
    interval : int
        The interval (in time steps) at which to write checkpoints.
    """

    _steps = [0]

    def __init__(self, file, interval):
        if isinstance(file, str):
            self._own_handle = True
            self._out = open(file, 'w', 0)
        else:
            self._out = file
            self._own_handle = False

        self.interval = interval
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
        steps = self.interval - simulation.currentStep % self.interval
        return steps, True, False, False, False

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
            self._initialized = True

        self._steps[0] += self.interval
        positions = state.getPositions()

        # Serialize
        self._out.write(b''.join([b'\nSTARTOFCHUNK\n',
                                  pickle.dumps([self._steps[0], positions._value]),
                                  b'\nENDOFCHUNK\n']))
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
    p = ArgumentParser(description='InsiliChem Ommprotocol: '
                       'easy to deploy MD protocols for OpenMM')
    p.add_argument('input', metavar='INPUT FILE', type=extant_file,
                   help='YAML input file')
    p.add_argument('--version', action='version', version='%(prog)s v{}'.format(__version__))
    p.add_argument('-c', '--check', action='store_true',
                   help='Validate input file only')
    args = p.parse_args(argv if argv else sys.argv[1:])

    jinja_env = jinja2.Environment(trim_blocks=True, lstrip_blocks=True)
    # Load config file
    with open(args.input) as f:
        rendered = jinja_env.from_string(f.read()).render()
        cfg = yaml.load(rendered, Loader=YamlLoader)
    # Paths and dirs
    from .md import SYSTEM_OPTIONS
    cfg['_path'] = os.path.abspath(args.input)
    cfg['system_options'] = prepare_system_options(cfg, defaults=SYSTEM_OPTIONS)
    cfg['outputpath'] = sanitize_path_for_file(cfg.get('outputpath', '.'), args.input)

    if not args.check:
        with ignored_exceptions(OSError):
            os.makedirs(cfg['outputpath'])

    handler = prepare_handler(cfg)

    return handler, cfg, args


def prepare_handler(cfg):
    """
    Load all files into single object.
    """
    positions, velocities, box = None, None, None
    _path = cfg['_path']
    forcefield = cfg.pop('forcefield', None)
    topology = sanitize_path_for_file(cfg.pop('topology'), _path)

    if 'checkpoint' in cfg:
        restart_path = sanitize_path_for_file(cfg['checkpoint'], _path)
        restart = Restart.load(restart_path)
        positions = restart.positions
        velocities = restart.velocities
        box = restart.box

    if 'positions' in cfg:
        positions_path = sanitize_path_for_file(cfg.pop('positions'), _path)
        positions = Positions.load(positions_path)
        box = BoxVectors.load(positions_path)

    if 'velocities' in cfg:
        velocities_path = sanitize_path_for_file(cfg.pop('velocities'), _path)
        velocities = Velocities.load(velocities_path)

    if 'box' in cfg:
        box_path = sanitize_path_for_file(cfg.pop('box'), _path)
        box = BoxVectors.load(box_path)

    options = {}
    for key in 'positions velocities box forcefield'.split():
        value = locals()[key]
        if value is not None:
            options[key] = value

    return SystemHandler.load(topology, **options)


def prepare_system_options(cfg, defaults=None):
    """
    Retrieve and delete (pop) system options from input configuration.
    """
    d = {} if defaults is None else defaults.copy()
    if 'nonbondedMethod' in cfg:
        d['nonbondedMethod'] = warned_getattr(openmm_app, cfg.pop('nonbondedMethod'), None)
    if 'nonbondedCutoff' in cfg:
        d['nonbondedCutoff'] = cfg.pop('nonbondedCutoff') * u.nanometers
    if 'constraints' in cfg:
        d['constraints'] = warned_getattr(openmm_app, cfg.pop('constraints'), None)
    for key in ['rigidWater', 'ewaldErrorTolerance']:
        if key in cfg:
            d[key] = cfg.pop(key)
    if 'extra_system_options' in cfg:
        if 'implicitSolvent' in cfg['extra_system_options']:
            implicit_solvent = warned_getattr(
                openmm_app, cfg['extra_system_options']['implicitSolvent'], None)
            cfg['extra_system_options']['implicitSolvent'] = implicit_solvent
        d.update(cfg.pop('extra_system_options'))
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
            gaffmol2 = os.path.splitext(forcefield)[0] + '.gaff.mol2'
            yield create_ffxml_file([gaffmol2], [forcefield])
        else:
            yield forcefield


def statexml2pdb(topology, state, output=None):
    """
    Given an OpenMM xml file containing the state of the simulation,
    generate a PDB snapshot for easy visualization.
    """
    state = Restart.from_xml(state)
    system = SystemHandler.load(topology, positions=state.positions)
    if output is None:
        output = topology + '.pdb'
    system.write_pdb(output)


def export_frame_coordinates(topology, trajectory, nframe, output=None):
    """
    Extract a single frame structure from a trajectory.
    """
    if output is None:
        basename, ext = os.path.splitext(trajectory)
        output = '{}.frame{}.inpcrd'.format(basename, nframe)

    # ParmEd sometimes struggles with certain PRMTOP files
    if os.path.splitext(topology)[1] in ('.top', '.prmtop'):
        top = AmberPrmtopFile(topology)
        mdtop = mdtraj.Topology.from_openmm(top.topology)
        traj = mdtraj.load_frame(trajectory, int(nframe), top=mdtop)
        structure = parmed.openmm.load_topology(top.topology, system=top.createSystem())
        structure.box_vectors = top.topology.getPeriodicBoxVectors()

    else:  # standard protocol (the topology is loaded twice, though)
        traj = mdtraj.load_frame(trajectory, int(nframe), top=topology)
        structure = parmed.load_file(topology)

    structure.positions = traj.openmm_positions(0)

    if traj.unitcell_vectors is not None:  # if frame provides box vectors, use those
        structure.box_vectors = traj.openmm_boxes(0)

    structure.save(output, overwrite=True)


###########################
# Defaults
###########################

REPORTERS = {
    'PDB': PDBReporter,
    'PDBX': PDBxReporter,
    'DCD': SegmentedDCDReporter,
    'CHK': CheckpointReporter,
    'HDF5':  HDF5Reporter,
    'NETCDF': NetCDFReporter,
    'MDCRD': MdcrdReporter,
    'RS': partial(RestartReporter, write_multiple=True, netcdf=True),
}
