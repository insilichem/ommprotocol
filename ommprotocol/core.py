#!/usr/bin/env python
# -*- coding: utf-8 -*-

#################################################
#           insiliChem OpenMM launcher          #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.com> @ 2016    #
#################################################

import os
import sys
from contextlib import contextmanager
from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
from mdtraj.reporters import HDF5Reporter
from ommprotocol.io import assert_not_exists, random_string

class Protocol(object):
    pass


class Stage(object):

    """
    Controls a simulation stage from a SystemHandler instance. It will handle
    the actual OpenMM system and then the Simulation object. Integrators,
    barostat and restraints are all easily handled too.

    Using it is easy: instantiate with a SystemHandler object and then call
    `run()`. However, you can also use it as an OpenMM high-level controller.


    Parameters
    ----------
    input_top : simtk.openmm.Topology
        The topology input file (PRMTOP, PDB)
    positions : simtk.Quantity, optional
        The starting coordinates of this stage. Only needed if
        input_top is a PRMTOP file.
    md_steps : int, optional
        Number of MD steps to simulate. If 0, no MD will take place
    timestep : float, optional
        Integration timestep, in fs. Defaults to 1.0.
    forcefields : list of str or file-like, optional
        Forcefields to apply in PDB inputs.
    velocities : simtk.unit.Quantity, optional
        The initial velocities of this stage. If None, they will be set
        to the requested temperature
    box_vectors : simtk.unit.Quantity, optional
        Replacement periodic box vectors, instead of input_top's.
    barostat : bool, optional
        True for NPT @ 1 atmosphere. False for NVT
    restrained_atoms, constrained_atoms : str or None, optional
        Parts of the system that should remain restrained or constrained
        during the stage. Available values in SELECTORS dict.
        If None, no atoms will be fixed.
    minimize : bool, optional
        If True, minimize before MD
    minimization_tolerance : float, optional, default=10kJ/mol
        Threshold value minimization should converge to
    minimization_max_iterations : int, optional, default=0
        Limit minimization iterations up to this value. If zero, don't limit.
    temperature : float, optional
        Target temperature of system in Kelvin, defaults to 300K
    trajectory : 'PDB' or 'DCD', optional
        Output format of trajectory file, if desired.
    trajectory_step : int, optional
        Frequency of trajectory write, in number of simulation steps
    restart_step : int, optional
        Frequencty of restart file creation. Defaults to 1E6 steps (1ns)
    stdout_step : int, optional
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
    precision : str, optional
        Precision model to use: single, double or mixed.
    system_kwargs : dict, optional
        Set of options to configure the system. See SYSTEM_KWARGS dict
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
        Wether to create a state.xml file at the end of the stage or not.
    """

    SELECTORS = {
        'all': lambda a: True, 
        None: lambda a: False,
        'protein': lambda a: a.residue.name not in ('WAT', 'HOH', 'TIP3') \
                   and a.name not in ('Cl-', 'Na+', 'SOD', 'CLA'),
        'protein_no_H': lambda a: a.residue.name not in ('WAT', 'HOH', 'TIP3') \
                        and a.name not in ('Cl-', 'Na+', 'SOD', 'CLA') and a.element.symbol != 'H',
        'backbone': lambda a: a.name in ('CA', 'C', 'N')
    }
    NONBONDEDMETHODS = {
        'NoCutoff': app.NoCutoff, '': app.NoCutoff, None: app.NoCutoff, 'None': app.NoCutoff,
        'CutoffNonPeriodic': app.CutoffNonPeriodic,
        'CutoffPeriodic': app.CutoffPeriodic,
        'Ewald': app.Ewald,
        'PME': app.PME
    }
    CONSTRAINTS = {
        '': None, None: None, 'None': None,
        'HBonds': app.HBonds,
        'AllBonds': app.AllBonds,
        'HAngles': app.HAngles
    }
    REPORTERS = {
        'PDB': app.PDBReporter,
        'DCD': app.DCDReporter,
        'HDF5':  HDF5Reporter
    }
    INTEGRATORS = {
        'BrownianIntegrator': mm.BrownianIntegrator,
        'LangevinIntegrator': mm.LangevinIntegrator,
    }
    PRECISION = {
        'CUDA': 'CudaPrecision',
        'OpenCL': 'OpenCLPrecision'
    }

    _PROJECTNAME = random_string(length=5)

    def __init__(self, handler, positions=None, velocities=None, box=None,
                 steps=0, minimization=True, barostat=True, temperature=300,
                 timestep=1.0, pressure=1.01325, integrator=None, friction=1.0,
                 barostat_interval=25, system_options=None, platform=None, precision=None,
                 trajectory=None, trajectory_every=10000, output='.', verbose=True,
                 restart=None, restart_every=1000000, report=True, report_every=1000,
                 project_name=None, name=None, restrained_atoms=None,
                 restraint_strength=5, constrained_atoms=None,
                 minimization_tolerance=10, minimization_max_iterations=0,
                 save_state_at_end=True):
        # System properties
        self.handler = handler
        self.positions = positions
        self.velocities = velocities
        self.box = box
        self.system_options = system_options
        self.restrained_atoms = restrained_atoms
        self.restraint_strength = restraint_strength
        self.constrained_atoms = constrained_atoms
        # Simulation conditions
        self.steps = steps
        self.minimization = minimization
        self.minimization_tolerance = minimization_tolerance
        self.minimization_max_iterations = minimization_max_iterations
        self._barostat = barostat
        self.temperature = temperature
        self.timestep = timestep
        self.pressure = pressure
        self._integrator = integrator
        self.friction = friction
        self.barostat_interval = barostat_interval
        # Hardware
        self._platform = platform
        self.precision = precision
        # Output parameters
        self.project_name = project_name if project_name is not None else self._PROJECTNAME
        self.name = name if name is not None else random_string(length=5)
        self.output = output
        self.verbose = verbose
        self.trajectory = trajectory
        self.trajectory_every = trajectory_every
        self.restart = restart
        self.restart_every = restart_every
        self.report = report
        self.report_every = report_every
        self.save_state_at_end = save_state_at_end
        # Private attributes
        self._system = None
        self._simulation = None
        self._mass_options = {}

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
            status = 'Running {} @ {}K'.format(self.name, self.temperature)
            status += ', NPT' if self.barostat else ', NVT'
            if self.restrained_atoms:
                status += ' [Restrained {}]'.format(self.restrained_atoms)
            if self.constrained_atoms:
                status += ' [Constrained {}]'.format(self.constrained_atoms)
            print(status)

        if self.minimization:
            if self.verbose:
                print('  Minimizing...')
            self.minimize()
        if not self.steps:
            return

        # Stdout report
        if self.report:
            mass = {'systemMass': self.system_mass} if self.constrained_atoms else {}
            rep = app.StateDataReporter(sys.stdout, self.report_every, step=True, 
                                        potentialEnergy=True, kineticEnergy=True, 
                                        temperature=True, volume=True, progress=True, 
                                        remainingTime=True, speed=True, 
                                        totalSteps=self.steps, separator='\t', **mass)
            self.simulation.reporters.append(rep)
        # Trajectory / movie files
        if self.trajectory:
            suffix = '.{}'.format(self.trajectory.lower())
            path = self.new_filename(suffix=suffix)
            rep = self.reporter(self.trajectory)(path, self.trajectory_every)
            self.simulation.reporters.append(rep)
        # Checkpoint or restart files
        if self.restart:
            suffix = '.{}'.format(self.restart.lower())
            path = self.new_filename(suffix=suffix)
            rep = self.reporter(self.restart)(path, self.restart_every)
            self.simulation.reporters.append(rep)
        
        if self.verbose:
            pbc = 'PBC' if self.system.usesPeriodicBoundaryConditions() else ''
            print('  Running {} MD for {} steps'.format(pbc, self.steps))
        with self.attempt_rescue(self.simulation):
            self.simulate()
       
        uses_pbc = self.system.usesPeriodicBoundaryConditions()
        state = self.simulation.context.getState(getPositions=True, getVelocities=True,
                                                 enforcePeriodicBox=uses_pbc)
        if self.save_state_at_end:
            path = self.new_filename(suffix='.state.xml')
            self.simulation.saveState(path)
        
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
            self._system = self.handler.create_system(**self.system_options)
        return self._system

    @system.deleter
    def system(self):
        self._system = None

    @property
    def simulation(self):
        if self._simulation is None:
            self._simulation = app.Simulation(self.handler.topology, self.system, 
                                                 self.integrator, self.platform)
            sim = self._simulation

            # Box vectors
            box = self.box if self.box is not None else self.handler.box
            if box is not None:
                sim.context.setPeriodicBoxVectors(*box)

            # Positions
            pos = self.positions if self.positions is not None else self.handler.positions
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
        self._simulation = None


    @property
    def integrator(self):
        try:
            return self.INTEGRATORS[self._integrator]
        except KeyError:
            raise NotImplementedError('Integrator {} not found'.format(self._integrator))

    @property
    def platform(self):
        if self._platform is not None:
            return {'platform': mm.Platform.getPlatformByName(self._platform),
                    'properties': {self.PRECISION[self._platform]: self.precision}}
    
    def reporter(self, name):
        try:
            return self.REPORTERS[name]
        except KeyError:
            raise NotImplementedError('Integrator {} not found'.format(name))

    def apply_barostat(self):
        if not self.system.usesPeriodicBoundaryConditions():
            raise ValueError('Barostat should not be used without PBC conditions.')
        self.system.addForce(mm.MonteCarloBarostat(self.pressure*u.bar,
                                                   self.temperature*u.kelvin,
                                                   self.barostat_interval))

    def apply_constraints(self):
        subset = self.subset(self.constrained_atoms)
        for i, atom in enumerate(self.handler.topology.atoms()):
            if subset(atom):
                self.system.setParticleMass(i, 0)

    def apply_restraints(self):
        force = self.restrain_force(self.restraint_strength)
        subset = self.subset(self.restrained_atoms)
        self.apply_force(self.handler, self.positions, force, subset)

    def subset(self, selector):
        try:
            return self.SELECTORS[selector]
        except KeyError:
            raise NotImplementedError('Selector {} not found'.format(selector))

    @property
    def system_mass(self):
        system_mass = sum(a.element.mass._value for a in self.handler.topology.atoms())
        return system_mass * u.dalton

    @staticmethod
    def apply_force(topology, positions, force, subset=None):
        if subset is None:
            subset = lambda x: True
        for i, (atom, position) in enumerate(zip(topology.atoms(), positions)):
            if subset(atom):
                force.addParticle(i, position.value_in_unit(u.nanometers))

    @staticmethod
    def restrain_force(strength=5.0):
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
        path = os.path.join(self.output, filename)
        if avoid_overwrite:
            path = assert_not_exists(path)
        return path


    @contextmanager
    def attempt_rescue(self, simulation, prompt=True, verbose=True, reraise=True):
        try:
            yield simulation
        except (Exception, KeyboardInterrupt) as ex:
            try:
                if prompt:
                    answer = input('Interruption detected. Attempt rescue? (y/N): ')
                    if answer not in ('Y', 'y', 'yes', 'YES'):
                        if verbose:
                            sys.exit('Ok, bye!')
                    elif verbose:
                        print('Ok!', end=' ')
                elif verbose:
                    print ('Interruption detected. Attempting rescue... ', end=' ')
                path = self.new_filename(suffix='.state.xml')
                self.simulation.saveState(path)
                state = simulation.context.getState(getPositions=True, getVelocities=True,
                                                    getForces=True, getEnergy=True,
                                                    getParameters=True)
                for reporter in simulation.reporters:
                    if not isinstance(reporter, app.StateDataReporter):
                        reporter.report(simulation, state)
            except Exception:
                if verbose:
                    print('FAILED :(')
            else:
                if verbose:
                    print('SUCCESS!')
            finally:
                if reraise:
                    raise ex
