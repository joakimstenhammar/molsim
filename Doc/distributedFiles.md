Distributed files
=================
This chapter contains a complete list of the distributed source files.

| Filename         | Description                                            |
| :------:         | :----------------------------------------------------: |
| `molsim.F90`     | Routines specific for  MOLSIM                          |
| `particle.F90`   | Routines for particles                                 |
| `potential.F90`  | Routines for potential and force tables                |
| `coordinate.F90` | Routines for coordinates including start configuration |
| `md.F90`         | Routines for MD integration                            |
| `mc.F90`         | Routines for MC integration                            |
| `bd.F90`         | Routines for BD integration                            |
| `nlist.F90`      | Routines for neighbour lists                           |
| `energy.F90`     | Routines calculating energy and forces                 |
| `dump.F90`       | Routines for dumping and reading variables             |
| `group.F90`      | Routines for particle group division                   |
| `static.F90`     | Routines for static analysis                           |
| `dynamic.F90`    | Routines for dynamic analysis                          |
| `image.F90`      | Routines for preparing image data files                |
| `mixed.F90`      | Routines specific for mixed tasks                      |
| `molaux.F90`     | Auxiliary routines                                     |
| `mesh.F90`       | Routines for making meshes                             |
| `statistics.F90` | Routines for statistical analyses                      |
| `mollib.F90`     | General library routines                               |
| `parallel.F90`   | Parallel interface between  MOLSIM and MPI             |
| `moluser.F90`    | User-provided routines                                 |
| `molsim.lib`     | Potential library                                      |

| Command    | Description                                               |
| :--------: | :-------------------------------------------------------: |
| molsim_ser | Command file for execution of serial version of  MOLSIM   |
| molsim_par | Command file for execution of parallel version of  MOLSIM |
| makefile   | Makefile for compiling and linking                        |
| archive    | Command file for filing                                   |
