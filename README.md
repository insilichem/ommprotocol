ommprotocol
===========

Easy to deploy MD protocols for OpenMM

By Jaime Rodríguez-Guerra ([@jaimergp](https://github.com/jaimergp))

Contributors: Lur Alonso

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
Quick launch:

    ommprotocol <inputfile.yaml>

All input configuration is done through YAML files. Some examples are included in `examples` folder; parameters should be self-explaining.

There's two main parts in these files: the global parameters, which are common for all stages, and the `stages` section, which lists all the stages to be simulated in the requested order. Each stage can override one or more global parameter, if needed.


Available parameters
--------------------
The following keys are available for the input file. They are listed in different categories, but you can mix them up as you want. If an optional key is not provided, it will be set to the default value.

### Input options
- `topology`: Main input file. Should contain, at least, the topology, but it can also contain positions, velocities, PBC vectors, forcefields... Supported file types: PDB, PRMTOP, PSF.
- `positions`: File with the initial coordinates of the system. Overrides those in topology, if needed. Supported file types: PDB, INPCRD, COOR.
- `velocities`: File containing the initial velocities of this stage. If not set, they will be set to the requested temperature. Supported file types: PDB, VEL.
- `box_vectors`: File with replacement periodic box vectors, instead of those in the topology or positions file. Supported file types: XSC.
- `forcefield`: Which forcefields should be used, if not provided in topology. Needed for PDB topologies.
- `charmm_parameters`: Parameters for PSF files. Like: *something.par, something.str.*
- `checkpoint`: Restart simulation from this file. It can provide one or more of the options above. Supported file types: xmlstate, rs.

### Output options
- `project_name`: Name for this simulation.
- `outputpath`: Path to output folder. If relative, it'll be relative to input file. 
- `report`: True for live report of progress.
- `report_every`: Update interval of live progress reports.
- `trajectory`: Output format of trajectory file, if desired.
- `trajectory_every`: Write trajectory every n steps.
- `trajectory_new_every`: Create a new file for trajectory every n steps.
- `restart`: Output format for restart formats, if desired
- `restart_every`: Write restart format every n steps.
- `save_state_at_end`: Whether to save the state of the simulation at the end of every stage.

### General conditions of simulation
- `minimize`: If *True*, minimize before MD.
- `steps`: Number of MD steps to simulate. If 0, no MD will take place.
- `timestep`: Integration timestep, in fs. Defaults to 1.0.
- `temperature`: In Kelvin.
- `barostat`: *True* for NPT, *False* for NVT.
- `pressure`: In bar. Only used if barostat is *True*.
- `barostat_every`: Update interval of barostat, in steps.
- `restrained_atoms`, `constrained_atoms`: Parts of the system that should remain restrained (a force is applied to minimize movement) or constrained (no movement at all) during the simulation. Available values: *all*, *protein*, *protein_no_H*, *backbone*, *calpha*. If *null*, no atoms will be fixed.
- `restraint_strength`: If restraints are in use, the strength of the applied force in kJ/mol. Defaults to 5.0.
- `integrator`: Which integrator should be used. Langevin by default.
- `friction`: Friction coefficient for integrator, if needed. In 1/ps. Defaults to 1.0.
- `minimization_tolerance` : Threshold value minimization should converge to.
- `minimization_max_iterations` : Limit minimization iterations up to this value. If zero, don't limit.

### OpenMM simulation parameters
- `nonbondedMethod`: Choose between *NoCutoff*, *CutoffNonPeriodic*, *CutoffPeriodic*, *Ewald*, *PME*.
- `nonbondedCutoff`: In nm. Defaults to 1.0.   
- `constraints`: Choose between *null*, *HBonds*, *AllBonds*, *HAngles*.
- `rigidWater`: *True* or *False*.
- `removeCMMotion`: Whether to remove center of mass motion during simulation. Defaults to *True*.  
- `hydrogenMass`: The mass to use for hydrogen atoms bound to heavy atoms. Any mass added to a hydrogen is subtracted from the heavy atom to keep their total mass the same.

### Hardware options
- `platform`: Which platform to use: *CPU*, *CUDA*, *OpenCL*. If not set, OpenMM will choose the fastest available.
- `precision`: Precision model to use: *single*, *double* or *mixed*.

### Stage-only parameters
- `name`: A name for this stage. 