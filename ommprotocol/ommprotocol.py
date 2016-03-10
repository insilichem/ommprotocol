#!/usr/bin/env python -W ignore::DeprecationWarning
# -*- coding: utf-8 -*-

#################################################
#         insiliChem.bio OpenMM launcher        #
# --------------------------------------------- #
# By Jaime RGP <jaime@insilichem.bio> @ 2016    #
#################################################

"""
ommprotocol
===========

insiliChem.bio OpenMM launcher for a standard MD workflow.
"""

# Python imports
from __future__ import division, print_function
from argparse import ArgumentParser
from collections import namedtuple
from datetime import datetime
import os
import random
import string
import sys
import warnings
# prepatch deprecation warnings (tmp fix)
warnings.warn = lambda *args, **kwargs: None
# OpenMM imports
try:
    from simtk.openmm import app
    import simtk.openmm as mm
    from simtk import unit
    from parmed.namd import NamdBinCoor, NamdBinVel
    from parmed.openmm import RestartReporter
    from parmed import load_file as parmed_load_file
    from openmoltools.utils import create_ffxml_file
    from mdtraj.reporters import HDF5Reporter
    import yaml
except ImportError:
    sys.exit('ERROR: Some dependencies are missing. Make sure you have '
             'installed openmm, parmed and openmoltools. Using Anaconda '
             'Python distribution is recommended. Install miniconda from '
             '(http://conda.pydata.org/miniconda.html) and run:\n'
             '    conda install -c omnia openmm parmed openmoltools pyyaml')
# Python2vs3
if sys.version_info.major == 2:
    input = raw_input

#=========================================================================
# Constants, defaults and mappings
#=========================================================================

PROJECT_NAME = ''.join(random.choice(string.ascii_letters) for _ in range(5))
FORCEFIELDS = ['amber99sbildn.xml', 'tip3p.xml']
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
SELECTORS = {
    'protein': lambda a: a.residue.name not in ('WAT', 'HOH') and a.name not in ('Cl-', 'Na+'),
    'protein_no_H': lambda a: a.residue.name not in ('WAT', 'HOH') and a.name not in ('Cl-', 'Na+') and a.element.symbol != 'H',
    'backbone': lambda a: a.name in ('CA', 'C', 'N')
}
INTEGRATORS = {
    'BrownianIntegrator': mm.BrownianIntegrator,
    'CustomIntegrator': mm.CustomIntegrator,
    'DrudeLangevinIntegrator': mm.DrudeLangevinIntegrator,
    'DrudeSCFIntegrator': mm.DrudeSCFIntegrator,
    'LangevinIntegrator': mm.LangevinIntegrator,
    'RPMDIntegrator': mm.RPMDIntegrator,
    'VariableLangevinIntegrator': mm.VariableLangevinIntegrator,
    'VariableVerletIntegrator': mm.VariableVerletIntegrator,
    'VerletIntegrator': mm.VerletIntegrator,
}
#=========================================================================
# Main functions
#=========================================================================


def stage(input_top, positions=None, forcefields=None, velocities=None, box_vectors=None,
          md_steps=0, minimize=True, barostat=True, temperature=300,
          restrained_atoms=None, constrained_atoms=None, timestep=1.0,
          trajectory=None, trajectory_step=10000, stdout_step=1000, restart_step=1000000,
          output='.', verbose=True, project_name=PROJECT_NAME, name=None,
          _system_kwargs=None, _platform=None, _restraint_strength=5,
          _pressure=1.01325, _integrator=None, _friction=1.0, **kwargs):
    """
    Create a new stage of the MD protocol, before production

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
    _platform: str, optional
        Which platform to use ('CPU', 'CUDA', 'OpenCL'). If not set,
        OpenMM will choose the fastest available.
    _system_kwargs: dict, optional
        Set of options to configure the system. See SYSTEM_KWARGS dict
        for defaults.
    _restraint_strength: float, optional
        If restraints are in use, the strength of the applied force in
        kJ/mol. Defaults to 5.0.
    _pressure: float, optional
        Barostat pressure, in bar. Defaults to 1.01325.
    _integrator: simtk.openmm.Integrator, optional
        Which integrator to use. Defaults to LangevinIntegrator.
    _friction: float, optional
        Friction coefficient for LangevinIntegrator, in 1/ps. Defaults to 1.0.

    Returns
    -------
    positions, velocities
    """
    # Defaults checking
    if kwargs:
        print('WARNING: Some options were not recognized: {}'.format(', '.join(kwargs)))
    if name is None:
        name = ''.join(random.choice(string.ascii_letters) for _ in range(5))
    if _system_kwargs is None:
        _system_kwargs = {}
        print('INFO: System options not specified. Using OpenMM default config:\n'
              '      nonbondedMethod=NoCutoff, nonbondedCutoff=1nm, constraints=None,\n'
              '      rigidWater=True, removeCMMotion=True')

    if verbose:
        status = 'Running {} @ {}K'.format(name, temperature)
        status += ', NPT' if barostat else ', NVT'
        if restrained_atoms:
            status += ' [Restrained {}]'.format(restrained_atoms)
        if constrained_atoms:
            status += ' [Constrained {}]'.format(constrained_atoms)
        print(status)

    # Create the system (universe)
    if isinstance(input_top, app.PDBFile):  # PDB in use
        if not forcefields:
            forcefields = FORCEFIELDS
            print('INFO: Forcefields for PDB not specified. Using default:\n      ',
                  ', '.join(forcefields))
        processed_forcefields = [_process_forcefield(ff) for ff in forcefields]
        forcefield = app.ForceField(*processed_forcefields)
        system = forcefield.createSystem(input_top.topology, **_system_kwargs)
    elif isinstance(input_top, app.AmberPrmtopFile):  # prmtop was used
        system = input_top.createSystem(**_system_kwargs)
    elif isinstance(input_top, app.CharmmPsfFile):  # psf was used
        system = input_top.createSystem(input_top.parmset, **_system_kwargs)
    else:
        sys.exit('ERROR: Input file format not recognized.')

    # Apply restraints or constraints, if requested
    mass_kwargs = {}
    if restrained_atoms in SELECTORS:
        force = apply_force(input_top.topology, positions, subset=SELECTORS[restrained_atoms],
                            strength=_restraint_strength)
        system.addForce(force)
    elif constrained_atoms in SELECTORS:
        system_mass = sum(a.element.mass._value for a in input_top.topology.atoms())
        mass_kwargs['systemMass'] = system_mass * unit.dalton
        apply_constraint(input_top.topology, system, subset=SELECTORS[constrained_atoms])

    if barostat:
        system.addForce(mm.MonteCarloBarostat(_pressure*unit.bar, temperature*unit.kelvin, 25))

    platform_kwargs = {}
    if _platform in ('Reference', 'CPU', 'CUDA', 'OpenCL'):
        platform = mm.Platform.getPlatformByName(_platform)
        platform_kwargs['platform'] = platform

    if _integrator is None:
        _integrator = mm.LangevinIntegrator(temperature*unit.kelvin,
                                            _friction/unit.picoseconds,
                                            timestep*unit.femtoseconds)
    simulation = app.Simulation(input_top.topology, system, _integrator, **platform_kwargs)
    simulation.context.setPositions(positions)

    if velocities is not None:
        simulation.context.setVelocities(velocities)
    else:
        simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)

    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)

    if minimize:
        if verbose:
            print('  Minimizing...')
        simulation.minimizeEnergy()

    if md_steps:
        if md_steps > stdout_step:
            simulation.reporters.append(
                app.StateDataReporter(sys.stdout, stdout_step, step=True, potentialEnergy=True,
                                      kineticEnergy=True, temperature=True,
                                      volume=True, progress=True, remainingTime=True, speed=True,
                                      totalSteps=md_steps, separator='\t', **mass_kwargs))

        if trajectory in REPORTERS and md_steps > trajectory_step:
            trajectory_out = new_filename_from(os.path.join(
                output, '{}_{}.{}'.format(project_name, name, trajectory.lower())))
            restart_out = new_filename_from(os.path.join(
                output, '{}_{}.restart'.format(project_name, name)))
            simulation.reporters.append(REPORTERS[trajectory](trajectory_out, trajectory_step))
            simulation.reporters.append(RestartReporter(restart_out, restart_step, True, True))

        if verbose:
            print('  Running MD with {} steps'.format(md_steps))
        try:
            simulation.step(md_steps)
        except (Exception, KeyboardInterrupt) as ex:
            if isinstance(ex, KeyboardInterrupt):
                report = input('\nCtr+C detected. Save current status? (y/N): ')
                if report not in ('Y', 'y', 'yes', 'YES'):
                    sys.exit('Ok, bye!')
            try:
                print('-'*70)
                print('Interruption detected. Attempting emergency report...', end=' ')
                state_xml = '{}_{}.state.xml'.format(project_name, name)
                simulation.saveState(new_filename_from(os.path.join(output, state_xml)))
                state = simulation.context.getState(getPositions=True, getVelocities=True,
                                                    getForces=True, getEnergy=True,
                                                    getParameters=True)
                for reporter in simulation.reporters:
                    if not isinstance(reporter, app.StateDataReporter):
                        reporter.report(simulation, state)
            except Exception:
                print('FAILED!')
            else:
                print('SUCCESS!')
            finally:
                print('-'*70)
                raise ex

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    state_xml = '{}_{}.state.xml'.format(project_name, name)
    simulation.saveState(new_filename_from(os.path.join(output, state_xml)))
    return state.getPositions(), state.getVelocities(), state.getPeriodicBoxVectors()


def run_protocol(argv=None):
    """
    Run a batch of stages with no additional boilerplate.

    Parameters
    ----------
    protocol: list of dict, optional
        A protocol is implemented as a list of dictionaries, each containing the
        desired arguments to call ``stage()``. See section *Default stages* for
        examples, as well as the list ``PROTOCOL``.
    default_user_options: dict, optional
        Set of options that should be shared across all stages of the protocol.
        Use it to avoid inserting the same values in all stages manually.
    """
    start_time = datetime.now()

    loaded_input, pos, vel, filename, out, box, protocol = args_parse(argv)
    basedir, name_with_ext = os.path.split(filename)
    output = out if out else basedir
    try:
        os.mkdir(output)
    except OSError:
        pass
    project_name, ext = os.path.splitext(name_with_ext)

    # Check Charmm parameters
    if isinstance(loaded_input, app.CharmmPsfFile):
        try:
            loaded_input.parmset = app.CharmmParameterSet(*protocol['charmm_parameters'])
        except KeyError:
            sys.exit('ERROR: For PSF input files you must specify CHARMM parameters in '
                     'protocol with key "charmm_parameters"')
        else:
            loaded_input.loadParameters(loaded_input.parmset)

    # Parse system options
    system_options = protocol.get('system', {})
    system_options['nonbondedMethod'] = NONBONDEDMETHODS.get(system_options.get('nonbondedMethod'))
    system_options['nonbondedCutoff'] = system_options.get('nonbondedCutoff', 1.0)*unit.nanometers
    system_options['constraints'] = CONSTRAINTS.get(system_options.get('nonbondedMethod'))
    # Prepare the options dict
    default_options = {'forcefields': protocol.get('forcefields'),
                       'project_name': project_name,
                       'output': output,
                       '_system_kwargs': system_options}

    for stage_options in protocol['stages']:
        options = default_options.copy()
        options.update(stage_options)
        pos, vel, box = stage(
            loaded_input, positions=pos, velocities=vel, box_vectors=box, **options)

    print('Done in', datetime.now()-start_time)


def args_parse(argv=None):
    parser = ArgumentParser(
        description='insiliChem.bio OpenMM launcher: easy to deploy MD protocols for OpenMM')
    parser.add_argument('input', metavar='INPUT FILE', type=str,
                        help='Topology file, in PDB, PSF, or PRMTOP format')
    parser.add_argument('-c', '--coordinates', type=str,
                        help='Initial coordinates (.inpcrd, .coor, .pdb), required for PRMTOP/PSF files')
    parser.add_argument('-v', '--velocities', type=str,
                        help='Initial velocities (.vel)')
    parser.add_argument('-b', '--box_vectors', type=str,
                        help='Box vectors in .xsc format (for namd restarts)')
    parser.add_argument('-r', '--restart', type=str,
                        help='Restart file of a previous run with this topology '
                        '(ignores velocities and positions)')
    parser.add_argument('-o', '--output', type=str,
                        help='Location of output files. If not set, same folder as INPUT.')
    parser.add_argument('-p', '--protocol', type=str,
                        help='MD protocol to follow, in YAML format.')
    args = parser.parse_args(argv if argv else sys.argv[1:])

    positions, velocities, box_vectors = None, None, None
    # Read restart file
    if args.restart:
        print('INFO: Loading positions, vectors and velocities from restart file')
        if args.restart.endswith('.xml'):
            with open(args.restart) as f:
                xml = mm.XmlSerializer.deserialize(f.read())
                positions = xml.getPositions()
                velocities = xml.getVelocities()
                box_vectors = xml.getPeriodicBoxVectors()
        else:
            rst = parmed_load_file(args.restart)
            positions = unit.Quantity(rst.coordinates[0], unit=unit.angstrom)
            if rst.hasvels:
                velocities = unit.Quantity(rst.velocities[0], unit=unit.angstrom/unit.picosecond)
            if rst.hasbox:
                box_vectors = [[v if i == j else 0 for j in range(3)]
                               for i, v in enumerate(rst.cell_lengths)] * unit.angstrom
    # no restart file, try with separate positions and velocities files
    else:
        # Read initial coordinates
        if args.coordinates:
            try:
                positions = app.AmberInpcrdFile(args.coordinates).positions
            except TypeError:  # coordinates are not .inpcrd, maybe Namd .coor?
                try:
                    namdbincoor = NamdBinCoor.read(args.coordinates)
                    positions_array = namdbincoor.coordinates[0]
                    positions = unit.Quantity(positions_array, unit=unit.angstroms)
                except:  # maybe pdb?
                    pdb = app.PDBFile(args.coordinates)
                    positions = pdb.positions
                    box_vectors = pdb.topology.getUnitCellDimensions()
        # Read initial velocities
        if args.velocities:
            vel = NamdBinVel.read(args.velocities)
            velocities_array = vel.velocities[0] / vel.SCALE_FACTOR
            velocities = unit.Quantity(velocities_array, unit=unit.angstroms/unit.picosecond)

        if args.box_vectors and args.box_vectors.endswith('.xsc'):
            xsc = parse_xsc(args.box_vectors)
            box_vectors = unit.Quantity([[xsc.a_x, 0, 0], [0, xsc.b_y, 0], [0, 0, xsc.c_z]],
                                        unit=unit.angstroms).in_units_of(unit.nanometers)

    # Topology is a PRMTOP file
    if args.input.endswith('.top') or args.input.endswith('.prmtop'):
        loaded_input = app.AmberPrmtopFile(args.input)
    # Topology is a PDB file
    elif args.input.endswith('.pdb'):
        loaded_input = app.PDBFile(args.input)
        if positions is None:
            positions = loaded_input.positions
        if velocities is None:
            velocities = getattr(loaded_input, 'velocities', None)
    elif args.input.endswith('.psf'):
        loaded_input = app.CharmmPsfFile(args.input)
        if positions is None:
            sys.exit('ERROR: PSF files require coordinates (with -c).')
        if box_vectors is None:
            print('    Info: Attempting automatic calculation of box vectors...')
            x, y, z = tuple(max((pos[i] for pos in positions))
                            - min((pos[i] for pos in positions)) for i in range(3))
            box_vectors = unit.Quantity([[x, 0, 0], [0, y, 0], [0, 0, z]],
                                        unit=unit.nanometers)
        loaded_input.setBox(*[v[i]._value for i, v in enumerate(box_vectors)])
    else:
        sys.exit('ERROR: Input file {} not recognized'.format(args.input))

    if not args.protocol:
        args.protocol = 'standard.yaml'
    try:
        yamlfile = open(args.protocol)
    except IOError:  # protocol not found in working dir, try builtin
        fname = args.protocol + '.yaml' if not args.protocol.endswith('.yaml') else args.protocol
        builtin_protocols = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
        try:
            yamlfile = open(os.path.join(builtin_protocols, fname))
        except IOError:
            sys.exit('ERROR: Protocol {} not recognized. Bad path?\n'
                     '       Also, you can try with the built-in protocols:'
                     ' {}'.format(args.protocol, ', '.join(sorted(os.listdir(builtin_protocols)))))

    protocol = yaml.load(yamlfile)
    yamlfile.close()

    return loaded_input, positions, velocities, args.input, args.output, box_vectors, protocol


#=========================================================================
# Helpers
#=========================================================================
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
    force.addGlobalParameter('k', strength*unit.kilocalories_per_mole/unit.angstroms**2)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    return force


def apply_force(topology, coordinates, subset=None, force=None, strength=5.0):
    """
    Apply a given force to a set of atoms
    """
    if force is None:
        force = restrain_force(strength=strength)
    if subset is None:
        subset = lambda x: True

    for i, (atom, coord) in enumerate(zip(topology.atoms(), coordinates)):
        if subset(atom):
            force.addParticle(i, coord.value_in_unit(unit.nanometers))
    return force


def apply_constraint(topology, system, subset=None):
    """
    Another style of fixing atom positions. Instead of a strong force on their
    initial coordinates, set their masses to zero. This is the OpenMM convention
    to apply constraints.
    """
    if subset is None:
        subset = lambda x: True
    for i, atom in enumerate(topology.atoms()):
        if subset(atom):
            system.setParticleMass(i, 0*unit.dalton)


def parse_xsc(path):
    with open(path) as f:
        lines = f.readlines()
        NamedXsc = namedtuple('NamedXsc', lines[1].split()[1:])
        return NamedXsc(*map(float, lines[2].split()))


def new_filename_from(path):
    name, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(path):
        path = '{}.{}{}'.format(name, i, ext)
        i += 1
    return path


def _process_forcefield(forcefield):
    if forcefield.endswith('.frcmod'):
        return create_ffxml_file(forcefield)
    return forcefield


if __name__ == '__main__':
    run_protocol()
