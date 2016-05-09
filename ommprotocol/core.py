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

    SELECTORS = {
        'all': lambda a: True, 
        None: lambda a: False,
        'protein': lambda a: a.residue.name not in ('WAT', 'HOH', 'TIP3') and a.name not in ('Cl-', 'Na+', 'SOD', 'CLA'),
        'protein_no_H': lambda a: a.residue.name not in ('WAT', 'HOH', 'TIP3') and a.name not in ('Cl-', 'Na+', 'SOD', 'CLA') and a.element.symbol != 'H',
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
                 restraint_strength=5, constrained_atoms=None):
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
        self.project_name = project_name if project_name else self._PROJECTNAME
        self.name = name if name else random_string(length=5)
        self.output = output
        self.verbose = verbose
        self.trajectory = trajectory
        self.trajectory_every = trajectory_every
        self.restart = restart
        self.restart_every = restart_every
        self.report = report
        self.report_every = report_every

        # Private attributes
        self._system = None
        self._simulation = None
        self._mass_options = {}

    def run(self):
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

        if self.report:
            mass = {'systemMass': self.system_mass} if self.constrained_atoms else {}
            rep = app.StateDataReporter(sys.stdout, self.report_every, step=True, 
                                        potentialEnergy=True, kineticEnergy=True, 
                                        temperature=True, volume=True, progress=True, 
                                        remainingTime=True, speed=True, 
                                        totalSteps=self.steps, separator='\t', **mass)
            self.simulation.reporters.append(rep)
        if self.trajectory:
            suffix = '.{}'.format(self.trajectory.lower())
            path = self.new_filename(suffix=suffix)
            rep = self.reporter(self.trajectory)(path, self.trajectory_every)
            self.simulation.reporters.append(rep)
        if self.restart:
            suffix = '.{}'.format(self.restart.lower())
            path = self.new_filename(suffix=suffix)
            rep = self.reporter(self.restart)(path, self.restart_every)
            self.simulation.reporters.append(rep)
        
        if self.verbose:
            print('  Running MD with {} steps'.format(self.steps))
        with self.attempt_rescue(self.simulation):
            self.simulate()
       

    def minimize(self):
        self.simulation.minimizeEnergy()

    def simulate(self, steps=None):
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
