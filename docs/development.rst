.. _development:

Software architecture
=====================

OMMProtocol is a glue application, which means that the main business logic resides within the third-party modules it depends on. Nonetheless, this should not necessarily imply a disorganized architecture. The main codebase is clearly divided in two categories: input and output handling (*io* module) and MD settings (*md* module). A third module, utils, collects miscellaneous functions that do not fall within the previous scopes. Finally, code concerning *ommanalyze* is stored in the analyze module.

Module *io*
-----------

This module hosts the input file handling logic, such as the precedence of format files (function *io.prepare_handler*), and the main container class (*io.SystemHandler*) that gives access to the components needed to create an OpenMM System object. Each of those components (*io.Topology*, *io.Positions*, *io.Velocities*, *io.BoxVectors*, and *io.Restart* objects) inherit from *io.MultiFormatLoader*, which supports the automated load of different formats based on the file extension, and io.InputContainer, a simple class that supports different attributes with light validation of the proper data structure.

The custom reporters provided by OMMProtocol are also contained here: *SegmentedDCDReporter* and *ProgressBarReporter*. The first allows the generation of DCD trajectories in chunked files to prevent huge file sizes, and the second converts OpenMM’s *StateDataReporter* in a more interactive console reporter (only one line dynamically updated for each protocol stage).

Module *md*
-----------

The goal of this module is to thread together the different stages of the protocol and run the corresponding simulations one after another. The main actor in this module is the *md.Stage* class, which contains all the needed logic to run a simulation in OpenMM: creation of the *System* object, application of restraints or constraints, preparation of the universe conditions such as temperature or pressure, configuration of the platform properties, construction of the Simulation object, setup of the output reporters… Each of these components is encapsulated in cached properties for maximum performance and ease of use in interactive sessions.

A helper function, *md.run_protocol*, takes the options for each stage specified in the input file and builds the needed *Stage* objects to execute them one after the other, passing the final state of each stage as the initial state of the next one. Since each stage must be named uniquely in the input file, the generated output files are meaningfully titled, leading to easy identification during the analysis.

Module *analyze*
----------------

The *ommanalyze* executable provides commands to perform routinary plots in trajectory analysis, like RMSD or potential energy plots. Currently, it only provides two subcommands: ``ommanalyze rmsd``, which requires the topology and one or more trajectory files, and outputs an interactive plot with matplotlib and ``ommanalyze log``, which simply plots the contents of the .log files generated during the trajectory. This module is only a stub that, if successful, could be further extended with more common analysis procedures thanks to the MDTraj library.