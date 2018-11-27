# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [4.6.1] - 2018-11-27
### Changed
- The recursive subroutine `UndoPBC` was replaced by an iterative version.
### Fixed
- Due to an overflow of the stack, the recursive subroutine `UndoPBC` crushed for systems with large numbers of particles per chain. The newly implemented version of `UndoPBC` chooses an iterative instead of a recursive strategy and does not suffer from the limited size of the stack.

## [4.6.0] - 2018-11-09
### Added
- The input variable `lreadbondcl` in `nmlParticle` is introduced. It is a flag for reading cross-linking information when `txstart=zero`.
### Fixed
- The problem that the crosslinks are not set in non-hierarchical structures if `txstart='zero'`- is fixed. The informations about crosslinking in the .cnf-file are read in the variable `bondcl` if `lreadbondcl=.true.`.

## [4.5.4] - 2018-11-08
### Fixed
- Changed the maximum number of cross-links of chain particles in a network to `nbondcl = 2` to allow for the simulation of networks with one particle per chain.

## [4.5.3] - 2018-11-07
### Fixed
- A check has been implemented in subroutine `SetNetwork` to prevent networks to interpenetrate.

## [4.5.2] - 2018-11-06
### Changed
- Subroutine `FileOpen` was modernized with regards to its usage of the `open` function
- Before, an `append` state of files was achieved by "manually" playing forward
- Now, the `append` state is explicitly achieved by using the `position` keyword of the `open` function

## [4.5.1] - 2018-10-30
### Fixed
- Prevented the dumping of the charged state for systems without weak charges

## [4.5.0] - 2018-09-24
### Added
- A manual generating by Doxygen was added. For that reason the code was newly commented. Some files concerning the manual were added in the directory doc to the  repository.
### Changed
- The calculation of the radius of gyration of the hierarchical structures was improved. In detail:
 - Now the root mean squared radius of gyration is calculated, before it was the mean absolute radius of gyration
 - The possibility to calculate the radius of gyration in the case of multiple hierarchical structures is added
- Read `bondcl` only when using continue, not for zero
- The Repository was moved to github

## [4.4.0] - 2018-01-03
### Fixed
- the variable `nvar` is set correctly for `chaindf`. This was not the case before (when `lhierarchical = .true.`), resulting in division by zero for the uninitialized positions in the `var` array in `DistFuncSamp ►le`.
- Compatibility issues with gfortran 4.8 were fixed
- The way the configure script determines which compiler to use is improved.
**NOTE** You have rerun the configure script when upgrading to this version
- The calculation of the cellsize in the case that a cell is as wide as the box is fixed.
- Some makefile dependencies are fixed
### Added
- a framework to give explicitly the random number seeds was added
### Changed
- the random number seeds are given in the output file
- `ntrydef`, the number of tries to set a configuration, is now an input parameter
- The makefile now has the binary executable as the ultimate target, so that make only runs when it is really of use

## [4.3.2] - 2017-08-23
### Fixed
- the speed of creating the celllist was increased, to have a more efficient method of carrying out npt calculations
- the periodic boundary conditions of `PBC` work also if the coordinates are more than one cell away
- when setting chains on a lattice, and the creation of the configuration for one chain fails, it tries again on the same lattice point, and not on the next one
- if the distance is smaller than `r2umin` and also smaller than `radatat` then it is treated as an hard core overlap
### Added
- added a `gprof` mode to the `gfortran` compilation to allow for profiling with `gprof`

## [4.3.1] - 2017-08-23
### Fixed
- fixed reading of variables from ucnf file, as the parameter have changed with version 2.2.0

## [4.3.0] - 2017-08-08
### Added
- Added an new and improved version of the cell list

## [4.2.2] - 2017-07-26
### Fixed
- Minor bug fixed in ImageVTF: Use if-construct instead of merge function to satisfy gfortran compiler

## [4.2.1] - 2017-07-26
### Fixed
- Adding stop for system with empty chains/ empty particle types

## [4.2.0] - 2017-07-25
### Changed
- Grouping scheme for systems containing weak charges (ionizable species)

## [4.1.0] - 2017-07-25
### Changed
- Molsim has been extended by the possibility to simulate multiple ionizable species with explicit counterions.

## [4.0.0] - 2017-07-17
### Changed
- The image generation via ImageVTF has been revised
- Coloring according to group assignment is now possible
- Splitting of VTF file is now possible

## [3.2.2] - 2017-07-11
### Fixed
- Fixed bugs in the parallel routines

## [3.2.1] - 2017-07-10
### Fixed
- Fixed a bug in the Testin Makefile

## [3.2.0] - 2017-07-05
### Changed
- Allow for parallel execution of Molsim of systems containing weak charges
- Enable the possibility to use Ewald summation together with weak charges
### Fixed
- Minor bugs related to weak charges

## [3.1.0] - 2017-07-05
### Changed
- Extended features for the usage of finite networks
- Change in network topology
- Network Statistics
- Network related grouping scheme

## [3.0.0] - 2017-07-05 - Myrapla
### Fixes
- fixed some minor bugs

### Changed
- changed the way the testin directory works
- molsim now uses the command line argumends in fortran


## [2.4.5] - 2017-04-24 - Quapsel
### Fixed
- fixes a bug in the sso routine, which prevented continuation of parallel simulations

## [2.4.4] - 2017-04-12 - Quapsel
### Fixed
- fixes a bug in the sso routine, which lead to "Not a Number" error.
- using the system wide fftw version, instead of a locally installed. One needs to run the configure script once when updating

## [2.4.3] - 2017-04-12 - Quapsel
### Fixed
- fixes a bug where sso calculated no results when using txstart=`continue`

## [2.4.2] - 2017-04-06 - Quapsel
### Fixed
- Fixed #33: the interactions of systems with weak charges and polyatomic particles was calculated wrongly

## [2.4.1] - 2017-03-22 - Quapsel
### Fixed
- removed all warnings occuding during the compilation of the code

## [2.4.0] - 2017-03-21 - Quapsel
### Changed
- the distribution of the complexation rate as well as the complexation rate of segments can be measured

## [2.3.0] - 2017-02-16 - Klikk
### Changed
- changed output about polymer chains and added txcopolymer == 'sequence'

## [2.2.0] - 2017-02-14 - Kokuna
### Changed
- changed the random number generator to one recomended by the numerical recipies (Press)

## [2.1.2] - 2016-01-30 - Tangoloss
### Fixed
- fixed wrong passing of arrays in some functions

### Added
- gfortran compile flags


## [2.1.1] - 2016-12-12 - Tangoloss
### Fixed
- when using `txpivot` as `short` or `lower`, all beads connected to the moved particles are also moved

## [2.1.0] - 2016-11-29 - Tangoloss
### Changed
- new network properties are being calculated: theta, rg2x, rg2y and rg2z
- please confer the [old Wiki](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/networks) for further information

## [2.0.0] - 2016-11-11 - Tangela
- Changes described in [Netzwerke (old wiki)](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/networks)

## [1.3.4] - 2016-11-04 Porygon
### Fixed
- fixed some errors occuring when compiling with `mode=warn`

## [1.3.3] - 2016-11-01 Porygon
### Changed
- allow `Hierarchical Move` to be used in par

## [1.3.2] - 2016-11-01 Porygon
### Fixed
- updated faulty `Testin` file

## [1.3.1] - 2016-10-26 - Porygon
### Fixed
- added support for `gfortran`

## [1.3.0] - 2016-10-17 - Porygon
- Changes described in [complexation-analysis (old wiki)](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/complexation-analysis)

## [1.2.1] - 2016-10-17 - Ditto
### Changed
- `molsim_par` also passes on exit code of Molsim

## [1.2.0] - 2016-10-12 - Ditto
- Changes described in [SSO (old wiki)](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/sso)

## [1.1.3] - 2016-09-07 - Smettbo
### Fixed
- Added requirement for group devision of `SubstructureDF`

## [1.1.2] - 2016-09-06 - Smettbo
- Changes described in [Bugfixes (old wiki)](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/bugfixes)


## [1.0.0] - 2016-08-15 - Raupy
### Changed
- The Makefile was adapted to work on the RWTH Cluster
- The README was updated
- Molsim exists with exit code `1` when an error occurs
- Giving `hpc` as a final argument saves the dump files on the `hpcwork` partition of the RWTH-Cluster

### Added
- Added `configure.sh` which is to be run before `make`. This file
    -  checks dependencies
    -  installs FFTW if neccesary
    -  manages naming the executable
- Added the Testin-directory to perform unit tests

## [0.0.0] - 2016-08-11 - Simsala
### Added
- The raw version from Per Linse was added.
- A README was added
- The version number was incoprorated

# Changed Input
 -  [4.0.0] in `nmlVtfImage` the second argument was removed from `tximage`, and the forner third and fourth argument are moved to position 2 and 3. See [
Image generation for VMD (old wiki)](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/imagevtf)
