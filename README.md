ommprotocol
===========

Easy to deploy MD protocols for OpenMM

By Jaime Rodríguez-Guerra ([@jaimergp](https://github.com/jaimergp))

Contributors: Lur Alonso, Patricia Saura

Installation
------------
Really easy with Conda. If you don't know what Conda is, check [its webpage](http://conda.pydata.org/docs/). Quick steps:

1. Download and install Miniconda.

        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3*.sh

2. Restart terminal and add `omnia` and `insilichem` to the channels list.

        conda config --add channels omnia
        conda config --add channels insilichem

3. Install ommprotocol in an environment called `openmm`:

        conda create -n openmm ommprotocol

4. Activate the environment

        source activate ommprotocol

5. Test

        ommprotocol -h


Usage
-----
Quick help:

    ommprotocol -h

Input file (last argument, no flags) can be:

 - Standard PDB (topology and coordinates). Protocol MUST include 
   [OpenMM forcefields](https://github.com/pandegroup/openmm/tree/master/wrappers/python/simtk/openmm/app/data) to use.
 - Amber's PRMTOP (topology and forcefields).
 - Charmm PSF with accompanying topology. Protocol MUST include parameters files in `charmm_parameters`.


Additional inputs:

 - Amber's inpcrd or Namd's coor for input coordinates (`-c`). Required for PRMTOP and PSF.
 - Namd's vel for velocities (`-v`).
 - Charmm .restart (used by this script) to restore positions and velocities (`-r`).
 - Protocol file (`-p`). This file can specify forcefields, Charmm parameters, and finally,
   the MD system options and stages.


How to implement a new protocol
-------------------------------
*Some example protocols ready to use are present in folder `data/`. Use them as examples.*

Protocols use a YAML file with at least one entry called `stages`.
This entry is a list of dictionaries, each specifying the options for that stage.
See example.yaml for completeness.

The YAML file can also contain another top-level entry called `system`,
specifying the global options for the system creation (non bonded methods, etc).

Additionally, the protocol can include different forcefields than default ones (amber99sbildn, tip3p)
with key `forcefields`. For example:
    
    forcefields:
        - amber10.xml
        - tip3p.xml

If you are using a PSF file as input, then it MUST contain a key `charmm_parameters` 
with the list of charmm parameters to load. It can be empty (`[]`), but it must be present. For example:

    charmm_parameters:
        - top.par
        - system.str

System defaults
---------------
These are the same as OpenMM default values. Just for quick reference:

| Variable                     | System key        | Default  | Other values                                  | Comments                                                                                                                                                  |
|------------------------------|-------------------|----------|-----------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Non bonded method            | `nonbondedMethod` | NoCutoff | CutoffNonPeriodic, CutoffPeriodic, Ewald, PME | [OpenMM docs](http://docs.openmm.org/6.2.0/userguide/application.html#nonbonded-interactions)                                                             |
| Non bonded cutoff (nm)       | `nonbondedCutoff` | 1.0      |                                               |                                                                                                                                                           |
| Constraints                  | `constraints`     | None     | HBonds, AllBonds, HAngles                     | [OpenMM docs](http://docs.openmm.org/6.2.0/userguide/application.html#constraints)                                                                        |
| Water rigidity               | `rigidWater`      | True     | False                                         |                                                                                                                                                           |
| Remove center of mass motion | `removeCMMotion`  | True     | False                                         |                                                                                                                                                           |
| Mass of hydrogen atoms       | `hydrogenMass`    | None     |                                               | The mass to use for hydrogen atoms bound to heavy atoms. Any mass added to a hydrogen is subtracted from the heavy atom to keep their total mass the same |


Stage defaults
--------------

| Variable                    | Stage key             | Default              | Other values                                                                                                                                                                                  | Comments                                                                                                           |
|-----------------------------|-----------------------|----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------|
| Forcefields                 | `forcefield`          | amber99sbildn, tip3p | Check [OpenMM forcefields](https://github.com/pandegroup/openmm/tree/master/wrappers/python/simtk/openmm/app/data)                                                                            | Only applies if the input is a PDB file                                                                            |
| Timestep (fs)               | `timestep`            | 1.0                  |                                                                                                                                                                                               |                                                                                                                    |
| Temperature (K)             | `temperature`         | 300                  |                                                                                                                                                                                               |                                                                                                                    |
| Barostat                    | `barostat`            | True                 | False                                                                                                                                                                                         |                                                                                                                    |
| Barostat pressure (bar)     | `_pressure`           | 1.01325              |                                                                                                                                                                                               |                                                                                                                    |
| Simulation steps            | `md_steps`            | 0                    |                                                                                                                                                                                               |                                                                                                                    |
| Restraints on atoms         | `restrained_atoms`    | None                 | backbone, protein, protein_no_H                                                                                                                                                               | A force is applied on selected atoms to largely reduce movement.                                                   |
| Restraint strength (kJ/mol) | `_restraint_strength` | 5.0                  |                                                                                                                                                                                               |                                                                                                                    |
| Constraints on atoms        | `constrained_atoms`   | None                 | backbone, protein, protein_no_H                                                                                                                                                               | Selected atoms won't move at all                                                                                   |
| Minimize structure          | `minimize`            | True                 | False                                                                                                                                                                                         |                                                                                                                    |
| Trajectory output format    | `trajectory`          | None                 | DCD, PDB, HDF5                                                                                                                                                                                |                                                                                                                    |
| Trajectory output frequency | `trajectory_step`     | 10000                |                                                                                                                                                                                               |                                                                                                                    |
| Console output frequency    | `stdout_step`         | 1000                 |                                                                                                                                                                                               | 0 to disable                                                                                                       |
| Restart output frequency    | `restart_step`        | 1000000              |                                                                                                                                                                                               | 0 to disable                                                                                                       |
| Output folder               | `output`              | .                    |                                                                                                                                                                                               | Relative to working directory                                                                                      |
| Protocol name output files  | `project_name`        | random string        |                                                                                                                                                                                               | The same for all stages                                                                                            |
| Stage name                  | `name`                | random string        |                                                                                                                                                                                               | Stage dependent                                                                                                    |
| Integrator                  | `_integrator`         | LangevinIntegrator   | BrownianIntegrator, CustomIntegrator, DrudeLangevinIntegrator, DrudeSCFIntegrator, LangevinIntegrator, RPMDIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator, VerletIntegrator | [OpenMM docs](http://docs.openmm.org/6.2.0/userguide/application.html#integrators)                                 |
| Integrator friction (1/ps)  | `_friction`           | 1.0                  |                                                                                                                                                                                               |                                                                                                                    |
| Override system defaults    | `_system_kwargs`      | Global ones          |                                                                                                                                                                                               | Check *System defaults* section and [1]                                                                            |

[1] System options can be also specified for each stage independently with `_system_kwargs`,
updating the global options. This means that a per-stage system options will take
precedence over the global system options in case of key conflict. If you already
defined these global options:

    system: {a: 0, b: 1}

and then in a stage you specify:

    _system_kwargs: {a: 1, c: 2}

The resulting options would be:

    system: {a:1, b: 1, c:2}