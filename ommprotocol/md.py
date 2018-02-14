#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ommprotocol: A command line application to launch
#              MD protocols with OpenMM
# By Jaime RGP <@jaimergp>

"""
ommprotocol.md
--------------

All the logic that governs the simulation of a MD protocol,
stage by stage.

A protocol is a chained list of stages, which are Stage instances.
"""

# Stdlib
from __future__ import print_function, division
import os
import sys
from contextlib import contextmanager
import logging
# 3rd party
from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
from mdtraj import Topology as MDTrajTopology
# Own
from .io import REPORTERS, ProgressBarReporter, SerializedReporter, prepare_system_options
from .utils import (random_string, assert_not_exists, timed_input,
                    available_platforms, warned_getattr)

logger = logging.getLogger(__name__)

###########################
# Defaults
###########################

FORCEFIELDS = ['amber99sbildn.xml', 'tip3p.xml']
SELECTORS = {
    None: 'none',
    'protein_no_H': 'protein and element != H',
    'calpha': 'name == CA'
}
PRECISION = {
    'CPU': None,
    'CUDA': 'CudaPrecision',
    'OpenCL': 'OpenCLPrecision'
}
SYSTEM_OPTIONS = {
    'nonbondedMethod': app.NoCutoff,
}
DEFAULT_OPTIONS = {
    'system_options': SYSTEM_OPTIONS,
}


def protocol(handler, cfg):
    """
    Run all the stages in protocol

    Parameters
    ----------
    handler : SystemHandler
        Container of initial conditions of simulation

    cfg : dict
        Imported YAML file.
    """
    # Stages
    if 'stages' not in cfg:
        raise ValueError('Protocol must include stages of simulation')

    pos, vel, box = handler.positions, handler.velocities, handler.box
    stages = cfg.pop('stages')
    for stage_options in stages:
        options = DEFAULT_OPTIONS.copy()
        options.update(cfg)
        stage_system_options = prepare_system_options(stage_options)
        options.update(stage_options)
        options['system_options'].update(stage_system_options)
        stage = Stage(handler, positions=pos, velocities=vel, box=box,
                      total_stages=len(stages), **options)
        pos, vel, box = stage.run()
        del stage


class Stage(object):

    """
    Controls a simulation stage from a SystemHandler instance. It will handle
    the actual OpenMM system and then the Simulation object. Integrators,
    barostat and restraints are all easily handled too.

    Using it is easy: instantiate with a SystemHandler object and then call
    `run()`. However, you can also use it as an OpenMM high-level controller.

    Parameters
    ----------
    handler : simtk.openmm.Topology
        The topology input file (PRMTOP, PDB)
    positions : simtk.Quantity, optional
        The starting coordinates of this stage. Only needed if
        handler is a PRMTOP file.
    steps : int, optional
        Number of MD steps to simulate. If 0, no MD will take place
    timestep : float, optional
        Integration timestep, in fs. Defaults to 1.0.
    forcefields : list of str or file-like, optional
        Forcefields to apply in PDB inputs.
    velocities : simtk.unit.Quantity, optional
        The initial velocities of this stage. If None, they will be set
        to the requested temperature
    box_vectors : simtk.unit.Quantity, optional
        Replacement periodic box vectors, instead of handler's.
    barostat : bool, optional
        True for NPT @ 1 atmosphere. False for NVT
    restrained_atoms, constrained_atoms : str or None, optional
        Parts of the system that should remain restrained or constrained
        during the stage. Available values in SELECTORS dict.
        If None, no atoms will be fixed.
    minimization : bool, optional
        If True, minimize before MD
    minimization_tolerance : float, optional, default=10 kJ/mol
        Threshold value minimization should converge to
    minimization_max_iterations : int, optional, default=10000
        Limit minimization iterations up to this value. If zero, don't limit.
    temperature : float, optional
        Target temperature of system in Kelvin, defaults to 300K
    trajectory : 'PDB' or 'DCD', optional
        Output format of trajectory file, if desired.
    trajectory_every : int, optional
        Frequency of trajectory write, in number of simulation steps
    trajectory_new_every : int, optional
        Create a new file for trajectory (only DCD) every n steps.
    trajectory_atom_subset : int, optional
        Save information for just these atoms (only DCD).
    restart_every : int, optional
        Frequencty of restart file creation. Defaults to 1E6 steps (1ns)
    report_every : int, optional
        Frequency of stdout print, in number of simulation steps
    verbose : bool, optional
        Whether to report information to stdout or not
    project_name : str, optional
        Name of the essay (common for several stages). If not set,
        five random characters will be used.
    name : str, optional
        Name of the stage, used as a suffix for the output files generated
        by this stage. If not supplied, a random string will be used.
    output : str, optional
        Location of output files. Working directory by default.
    platform : str, optional
        Which platform to use ('CPU', 'CUDA', 'OpenCL'). If not set,
        OpenMM will choose the fastest available.
    platform_properties : dict, optional
        Additional options to be passed to the platform constructor.
    system_options : dict, optional
        Set of options to configure the system. See SYSTEM_OPTIONS dict
        for defaults.
    restraint_strength : float, optional
        If restraints are in use, the strength of the applied force in
        kJ/mol. Defaults to 5.0.
    pressure : float, optional
        Barostat pressure, in bar. Defaults to 1.01325.
    integrator : simtk.openmm.Integrator, optional
        Which integrator to use. Defaults to LangevinIntegrator.
    friction : float, optional
        Friction coefficient for LangevinIntegrator, in 1/ps. Defaults to 1.0.
    barostat_interval : float, optional
        Interval of steps at which barostat updates. Defaults to 25 steps.
    save_state_at_end : bool, optional
        Whether to create a state.xml file at the end of the stage or not.
    attempt_rescue : bool, optional
        Whether to try to generate an emergency state file if an exception
        is raised.
    total_stages : int, optional
    """
    _PROJECTNAME = random_string(length=5)
    _stage_number = [0]

    def __init__(self, handler, positions=None, velocities=None, box=None,
                 steps=0, minimization=False, barostat=False, temperature=300,
                 timestep=1.0, pressure=1.01325, integrator='LangevinIntegrator',
                 barostat_interval=25, system_options=None, platform=None,
                 platform_properties=None, trajectory=None, trajectory_every=2000,
                 outputpath='.', trajectory_atom_subset=None, trajectory_new_every=0,
                 restart=None, restart_every=1000000, report=True, report_every=1000,
                 project_name=None, name=None, restrained_atoms=None,
                 restraint_strength=5, constrained_atoms=None, friction=1.0,
                 minimization_tolerance=10, minimization_max_iterations=10000,
                 save_state_at_end=True, attempt_rescue=True,
                 total_stages=None, verbose=True,
                 **kwargs):
        for k in kwargs:
            if not k.startswith('_'):
                logger.warning('Option %s not recognized!', k)

        # System properties
        self.handler = handler
        self.positions = positions
        self.velocities = velocities
        self.box = box
        self.system_options = system_options if system_options else {}
        self.restrained_atoms = restrained_atoms
        self.restraint_strength = restraint_strength
        self.constrained_atoms = constrained_atoms
        # Simulation conditions
        self.steps = int(steps)
        self.minimization = minimization
        self.minimization_tolerance = minimization_tolerance
        self.minimization_max_iterations = minimization_max_iterations
        self.barostat = barostat
        self.temperature = temperature
        self.timestep = timestep
        self.pressure = pressure
        self._integrator_name = integrator
        self.friction = friction
        self.barostat_interval = int(barostat_interval)
        # Hardware
        self._platform = platform
        self.platform_properties = platform_properties
        # Output parameters
        self.project_name = project_name if project_name is not None else self._PROJECTNAME
        self.name = name if name is not None else random_string(length=5)
        self.outputpath = outputpath
        self.verbose = verbose
        self.trajectory = trajectory
        self.trajectory_every = int(trajectory_every)
        self.trajectory_new_every = int(trajectory_new_every)
        self.trajectory_atom_subset = self._mask_selection(trajectory_atom_subset)
        self.restart = restart
        self.restart_every = int(restart_every)
        self.report = report
        self.report_every = int(report_every)
        self.save_state_at_end = save_state_at_end
        self.attempt_rescue = attempt_rescue
        self.total_stages = total_stages
        # Private attributes
        self._system = None
        self._simulation = None
        self._integrator = None
        self._progress_reporter = None
        self._log_reporter = None
        self._trajectory_reporter = None
        self._restart_reporter = None
        self._mass_options = {}
        self._stage_number[0] += 1

    def run(self):
        """
        Launch MD simulation, which may consist of:
        1. Optional minimization
        2. Actual MD simulation, with n steps.

        This method also handles reporters.

        Returns
        -------
        positions, velocities : unit.Quantity([natoms, 3])
            Position, velocity of each atom in the system
        box: unit.Quantity([1, 3])
            Periodic conditions box vectors
        """
        if self.verbose:
            status = '#{}'.format(self._stage_number[0])
            if self.total_stages is not None:
                status += '/{}'.format(self.total_stages)
            status += ': {}'.format(self.name)
            if self.restrained_atoms:
                status += ' [Restrained {}]'.format(self.restrained_atoms)
            elif self.constrained_atoms:
                status += ' [Constrained {}]'.format(self.constrained_atoms)
            logger.info(status)

        # Add forces
        if self.restrained_atoms:
            self.apply_restraints()
        if self.constrained_atoms:
            self.apply_constraints()
        if self.barostat:
            self.apply_barostat()

        if self.minimization:
            if self.verbose:
                logger.info('  Minimizing...')
            self.minimize()

        uses_pbc = self.system.usesPeriodicBoundaryConditions()
        if self.steps:
            # Stdout progress
            if self.report and self.progress_reporter not in self.simulation.reporters:
                self.simulation.reporters.append(self.progress_reporter)

            # Log report
            if self.report and self.log_reporter not in self.simulation.reporters:
                self.simulation.reporters.append(self.log_reporter)

            # Trajectory / movie files
            if self.trajectory and self.trajectory_reporter not in self.simulation.reporters:
                self.simulation.reporters.append(self.trajectory_reporter)

            # Checkpoint or restart files
            if self.restart and self.restart_reporter not in self.simulation.reporters:
                self.simulation.reporters.append(self.restart_reporter)

            # MD simulation
            if self.verbose:
                pbc = 'PBC ' if uses_pbc else ''
                conditions = 'NPT' if self.barostat else 'NVT'
                logger.info('  Running {}MD for {} steps @ {}K, {}'.format(pbc, self.steps,
                                                                     self.temperature,
                                                                     conditions))

            with self.handle_exceptions():
                self.simulate()

        if self.save_state_at_end:
            path = self.new_filename(suffix='.state')
            self.simulation.saveState(path)

        # Save and return state
        state = self.simulation.context.getState(getPositions=True, getVelocities=True,
                                                 enforcePeriodicBox=uses_pbc)

        return state.getPositions(), state.getVelocities(), state.getPeriodicBoxVectors()

    def minimize(self, tolerance=None, max_iterations=None):
        """
        Minimize energy of the system until meeting `tolerance` or
        performing `max_iterations`.
        """
        if tolerance is None:
            tolerance = self.minimization_tolerance
        if max_iterations is None:
            max_iterations = self.minimization_max_iterations
        self.simulation.minimizeEnergy(tolerance * u.kilojoules_per_mole, max_iterations)

    def simulate(self, steps=None):
        """
        Advance simulation n steps
        """
        if steps is None:
            steps = self.steps
        self.simulation.step(steps)

    @property
    def system(self):
        if self._system is None:
            if self.constrained_atoms and self.system_options.pop('constraints', None):
                logger.warning('  Warning: `constraints` and `constrained_atoms` are incompatible. '
                               'Removing `constraints` option for this stage.')
            self._system = self.handler.create_system(**self.system_options)
        return self._system

    @system.deleter
    def system(self):
        del self._system
        self._system = None

    @property
    def simulation(self):
        if self._simulation is None:
            platform = self.platform
            try:
                sim = self._simulation = app.Simulation(self.handler.topology, self.system,
                                                        self.integrator, *platform)
            except Exception as e:
                template = '{}. Try with: {}.'
                if 'Illegal property name' in str(e):
                    msg = template.format(e, ', '.join(platform[0].getPropertyNames()))
                    raise ValueError(msg)
                elif 'There is no registered Platform' in str(e):
                    msg = template.format(e, ', '.join(available_platforms()))
                    raise ValueError(msg)
                raise e
            # Box vectors
            box = self.box if self.box is not None else self.handler.box
            if box is not None:
                sim.context.setPeriodicBoxVectors(*box)

            # Positions
            pos = self.positions if self.positions is not None else self.handler.positions
            if pos is None:
                raise ValueError('Positions must be set to start a simulation.')
            sim.context.setPositions(pos)

            # Velocities
            vel = self.velocities if self.velocities is not None else self.handler.velocities
            if vel is not None:
                sim.context.setVelocities(vel)
            else:
                sim.context.setVelocitiesToTemperature(self.temperature*u.kelvin)
        return self._simulation

    @simulation.deleter
    def simulation(self):
        del self._simulation
        self._simulation = None

    @property
    def integrator(self):
        if self._integrator is None:
            try:
                i = getattr(mm, self._integrator_name)
            except (TypeError, AttributeError):
                raise NotImplementedError('Integrator {} not found'
                                          .format(self._integrator_name))
            self._integrator = i(self.temperature * u.kelvin,
                                 self.friction / u.picoseconds,
                                 self.timestep * u.femtoseconds)
        return self._integrator

    @integrator.deleter
    def integrator(self):
        del self._integrator
        self._integrator = None

    @property
    def platform(self):
        if self._platform is None:
            return None,
        platform = mm.Platform.getPlatformByName(self._platform)
        if self.platform_properties is None:
            return platform,

        # Patch to allow env-defined GPUs
        device = self.platform_properties.get('DeviceIndex', '')
        if str(device).startswith('ENV_'):
            envvar = os.environ.get(device[4:], None)
            if envvar is not None:
                logger.warning('Setting DeviceIndex from env var %s to %s', device[4:], envvar)
                self.platform_properties['DeviceIndex'] = envvar

        return platform, self.platform_properties

    def reporter(self, name):
        try:
            return REPORTERS[name.upper()]
        except KeyError:
            raise NotImplementedError('Reporter {} not found'.format(name))

    def apply_barostat(self):
        if not self.system.usesPeriodicBoundaryConditions():
            raise ValueError('Barostat can only be used with PBC conditions.')
        self.system.addForce(mm.MonteCarloBarostat(self.pressure*u.bar,
                                                   self.temperature*u.kelvin,
                                                   self.barostat_interval))

    def apply_constraints(self):
        indices = self.subset(self.constrained_atoms)
        system = self.system
        for i in indices:
            system.setParticleMass(int(i), 0.0)

    def apply_restraints(self):
        force = self.restraint_force(self.restraint_strength)
        indices = self.subset(self.restrained_atoms)
        self.apply_force(force, indices)

    @property
    def progress_reporter(self):
        if self._progress_reporter is None:
            if os.environ.get('OMMPROTOCOL_SLAVE'):
                rep = SerializedReporter(sys.stdout, self.report_every)
            else:
                rep = ProgressBarReporter(sys.stdout, self.report_every, total_steps=self.steps)
            self._progress_reporter = rep
        return self._progress_reporter

    @progress_reporter.deleter
    def progress_reporter(self):
        try:
            self.simulation.reporters.remove(self._progress_reporter)
        except ValueError:
            pass
        self._progress_reporter = None

    @property
    def log_reporter(self):
        if self._log_reporter is None:
            mass = {'systemMass': self.system_mass} if self.constrained_atoms else {}
            path = self.new_filename(suffix='.log')
            rep = app.StateDataReporter(path, self.report_every, step=True,
                                        potentialEnergy=True, kineticEnergy=True,
                                        temperature=True, volume=True, progress=True,
                                        remainingTime=True, speed=True,
                                        totalSteps=self.steps, separator='\t', **mass)
            self._log_reporter = rep
        return self._log_reporter

    @log_reporter.deleter
    def log_reporter(self):
        try:
            self.simulation.reporters.remove(self._log_reporter)
        except ValueError:
            pass
        self._log_reporter = None

    @property
    def trajectory_reporter(self):
        if self._trajectory_reporter is None:
            suffix = '.{}'.format(self.trajectory.lower())
            path = self.new_filename(suffix=suffix)
            options = {}
            if self.trajectory == 'DCD':
                options.update({'new_every': self.trajectory_new_every,
                                'atomSubset': self.trajectory_atom_subset})
            rep = self.reporter(self.trajectory)(path, self.trajectory_every, **options)
            self._trajectory_reporter = rep
        return self._trajectory_reporter

    @trajectory_reporter.deleter
    def trajectory_reporter(self):
        try:
            self.simulation.reporters.remove(self._trajectory_reporter)
        except ValueError:
            pass
        self._trajectory_reporter = None

    @property
    def restart_reporter(self):
        if self._restart_reporter is None:
            suffix = '.{}.{}'.format(self.restart.lower(), self.restart_every)
            path, ext_or_int, n = self.new_filename(suffix=suffix).rsplit('.', 2)
            try:
                ext_or_int = int(ext_or_int)  # Is ext an integer?
            except ValueError:  # Ext is the actual file extension
                path = '{}.{}'.format(path, ext_or_int)
            else:  # Ext is an int! Reformat
                name, ext = os.path.splitext(path)
                path = '{}.{}{}'.format(name, ext_or_int, ext)
            rep = self.reporter(self.restart)(path, self.restart_every)
            self._restart_reporter = rep
        return self._restart_reporter

    @restart_reporter.deleter
    def restart_reporter(self):
        try:
            self.simulation.reporters.remove(self._restart_reporter)
        except ValueError:
            pass
        self._restart_reporter = None

    def subset(self, selector):
        """
        Returns a list of atom indices corresponding to a MDTraj DSL
        query. Also will accept list of numbers, which will be coerced
        to int and returned.
        """
        if isinstance(selector, (list, tuple)):
            return map(int, selector)
        selector = SELECTORS.get(selector, selector)
        mdtop = MDTrajTopology.from_openmm(self.handler.topology)
        return mdtop.select(selector)

    @property
    def system_mass(self):
        system_mass = sum(a.element.mass._value for a in self.handler.topology.atoms())
        return system_mass * u.dalton

    def apply_force(self, force, indices=None):
        positions = self.positions if self.positions is not None else self.handler.positions
        if indices is None:
            indices = range(self.handler.topology.getNumAtoms())
        for i, index in enumerate(indices):
            force.addParticle(i, positions[index].value_in_unit(u.nanometers))

    @staticmethod
    def restraint_force(strength=5.0):
        """
        Force that restrains atoms to fix their positions, while allowing
        tiny movement to resolve severe clashes and so on.

        Returns
        -------
        force : simtk.openmm.CustomExternalForce
            A custom force to restrain the selected atoms
        """
        force = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        force.addGlobalParameter('k', strength*u.kilocalories_per_mole/u.angstroms**2)
        force.addPerParticleParameter('x0')
        force.addPerParticleParameter('y0')
        force.addPerParticleParameter('z0')
        return force

    def new_filename(self, suffix='', prefix='', avoid_overwrite=True):
        filename = '{}{}_{}{}'.format(prefix, self.project_name, self.name, suffix)
        path = os.path.join(self.outputpath, filename)
        if avoid_overwrite:
            path = assert_not_exists(path)
        return path

    @contextmanager
    def handle_exceptions(self, verbose=True):
        """
        Handle Ctrl+C and accidental exceptions and attempt to save
        the current state of the simulation
        """
        try:
            yield
        except (KeyboardInterrupt, Exception) as ex:
            if not self.attempt_rescue:
                raise ex
            if isinstance(ex, KeyboardInterrupt):
                reraise = False
                answer = timed_input('\n\nDo you want to save current state? (y/N): ')
                if answer and answer.lower() not in ('y', 'yes'):
                    if verbose:
                        sys.exit('Ok, bye!')
            else:
                reraise = True
                logger.error('\n\nAn error occurred: %s', ex)
            if verbose:
                logger.info('Saving state...')
            try:
                self.backup_simulation()
            except Exception:
                if verbose:
                    logger.error('FAILED :(')
            else:
                if verbose:
                    logger.info('SUCCESS!')
            finally:
                if reraise:
                    raise ex
                sys.exit()

    def backup_simulation(self):
        """
        Creates an emergency report run, .state included
        """
        path = self.new_filename(suffix='_emergency.state')
        self.simulation.saveState(path)
        uses_pbc = self.system.usesPeriodicBoundaryConditions()
        state_kw = dict(getPositions=True, getVelocities=True,
                        getForces=True, enforcePeriodicBox=uses_pbc,
                        getParameters=True, getEnergy=True)
        state = self.simulation.context.getState(**state_kw)
        for reporter in self.simulation.reporters:
            if not isinstance(reporter, app.StateDataReporter):
                reporter.report(self.simulation, state)

    def _mask_selection(self, expression):
        pass
