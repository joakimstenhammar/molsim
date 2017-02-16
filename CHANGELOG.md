# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

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
- please confer the [Wiki](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/networks) for further information

## [2.0.0] - 2016-11-11 - Tangela
- Changes described in [Netzwerke](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/networks)

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
- Changes described in [complexation-analysis](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/complexation-analysis)

## [1.2.1] - 2016-10-17 - Ditto
### Changed
- `molsim_par` also passes on exit code of Molsim

## [1.2.0] - 2016-10-12 - Ditto
- Changes described in [SSO](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/sso)

## [1.1.3] - 2016-09-07 - Smettbo
### Fixed
- Added requirement for group devision of `SubstructureDF`

## [1.1.2] - 2016-09-06 - Smettbo
- Changes described in [Bugfixes](https://git.rwth-aachen.de/pascal.hebbeker/Molsim/wikis/bugfixes)


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
