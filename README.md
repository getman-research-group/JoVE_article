# Supporting files and scripts for JoVE article.

mcpliq is a pre-compiled executable for the MCPliQ code used to add water molecules to the simulation cell. The provided mcpliq executable was compiled with GCC version 7.1.0. 

The MCPliQ.c file contains the source code for this for this executable. If changes to the code are desired, the appropriate changes can be made and the code can be recompiled using the gcc compiler and calling the command 'gcc -lm MCPliQ.c -o mcpliq'. The -lm option is required to include the proper math libraries during compilation.
