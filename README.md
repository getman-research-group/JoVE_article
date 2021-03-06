# Files and scripts for JoVE article.

<h2> MCPliQ </h2>

See https://github.com/getman-research-group/MCPliQ for more information about MCPliQ.

<h3> Pre-compiled executable </h3>

mcpliq is a pre-compiled executable for the MCPliQ code used to add water molecules to the simulation cell. The provided mcpliq executable was compiled with GCC version 7.1.0. 

<h3> Source code </h3>

The MCPliQ.c file contains the source code for this for this executable. If changes to the code are desired, the appropriate changes can be made and the code can be recompiled using the gcc compiler and calling the command 'gcc -lm MCPliQ.c -o mcpliq'. The -lm option is required to include the proper math libraries during compilation.


<h2> Python Scripts </h2>

All provided Python scripts were developed and tested using Python version 3.6.


<h2> Example Files </h2>

Example LAMMPS files for a CH2OH* adsorbate on a Pt(111) surface with 24 water molecules are supplied in the Example_Files directory. These include a data file, an input.equil file for an NPT simulation, and an input.prod file for an NVT simulation.

Note: Since the CH2OH* adsorbate is held fixed during the LAMMPS simulations, the bonds and angles within the adsorbate and the corresponding bond and angle force field parameters are not defined in these files.
