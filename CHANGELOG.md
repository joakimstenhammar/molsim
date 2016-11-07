# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

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