Developing
============
When developing Molsim, please adhere to the contributing guide, which is shipped in the MOLSIM repository.

# Testin {#Testin}
In the Molsim directory you will find a directory called `Testin`. The Testin offers an efficient way to check whether changes of the Molsim source code lead to any malfunction of the program.

### General idea
The general idea of the Testin-directory is to account for the sustainment of the functionality of Molsim. Whenever the source code of Molsim is changed, one has to ensure that none of the functionalities of Molsim have been broken. Despite of a few aspects (e.g. format changes) the program execution should yield consistent results. By comparing the results of test simulations of the current with the previous version it is easy to detect whatever kind of malfunctions.

The Testin-directory provides a method to compare the output of the current version to the last stable version.

### Usage
* **First you need to generate the output from the stable version.** The stable version should be the latest `master` branch from which your branch diverged. If you are unsure, check when the file `Testin/stable.md5sum` was last changed. The repository at this stage is the latest stable version. When you have checked out the stable version you want to compare with, run `make stable` in the `Testin` directory. This will compile the current source code with the `mode = test` flag and create the output in the `Testin/out_stable` directory.
* at any time you can run `make out` to run the same input files and store the output in the `Testin/out`.
* to compare the output of the `out` and `out_stable` directories, run `make diff`. Note that if the out files are not present, `make diff` will create them accordingly. The collected differences are written into the file `diff.out`

The `diff` should, in general, yield no differences. If differences were found, one needs to carefully inspect each of these instances for whether one deals with expected changes or unexpected changes.
* *Unexpected changes: Investigate where the difference is coming from and fix the bug!*
* Expected changes: When you have only expected changes, add comments to the `diff.out` file, which describe the changes. This commented version of the `diff.out` files should be attached to the pull request when you want to merge into the stable version. Note that the `diff.out` file will be overwritten when running `make` so you might want to add you comments to a copy of `diff.out`
* After your merge request has been positively reviewed, run `make declarestable` to declare your current version of the code as the stable one. The files included in the check for the stable version are all `.F90` files, the `makefile` and the `molsim_ser` file in the `Src` directory. Additionally also the `configure.sh` file in the molsim root is checked.

* to delete all output of the `out` runs (and the `diff.out` file!) run `make clean`, to delete the `out_stable` dir run `make cleanstable`. To delete both run `make cleanall`.

### Further notes:
* You can speed up the process, by running make in parallel (e.g. use `make -j 4`, to run on 4 cores).
* When running `make` all needed files are generated automatically, so you do not need to run `make out` before running `make diff`. A single call of `make` is enough.
* The `out_stable` files can not be created automatically, when you are in the directory of a modified version of molsim. Checkout the stable version, run `make stable`, and go back to your branch. The `out_stable` dir should stay, and can be used to compare with.

### List of elements contained in the Testin-directory
| element                    | description                                        |
|----------------------------|----------------------------------------------------|
| `in`                       | input files covering all functionalities of Molsim |
| `todiff.txt`               | all output files to diff                           |
| `Makefile`                 | Makefile to run all tests                          |
| `scripts/diff1.sh`         | script to compare the output files                 |
| `scripts/molsim_stable.sh` | script to generate outputs from the stable version |
| `stable.md5sum`            | file defining the current stable version           |

After executing `make` you will additionally find:

| element      | description                                                                                   |
|--------------|-----------------------------------------------------------------------------------------------|
| `out`        | output files resulting from the Molsim version under examination                              |
| `out_stable` | output files corresponding to the input file in `in`, created with the latest stable version. |
| `diff.out`   | Output of the diff between the current and the stable versions                                |

# Random Numbers {#RandomNumbers}

The random number generator used in Molsim is `rand` function as described by Press et al. in Numerical Recipes in Fortran 77, *The Art of Scientific Computing*, Second Edition, 1997. It has a a period of \f$\approx 3.1 \times 10^{18}\f$ and returns a uniform random number in the range (0.0, 1.0).

The routine uses a pair of of integers (`ixseed` and `iyseed`) to generate random numbers. The `iseed` you give in the `nmlSystem` is used to generate the first integer number pair. When `iseed` is negative the current time is used to set `iseed`. In most cases it is sufficient to set only `iseed`. The case where you also want to set `ixseed` and `iyseed` are when you want to use a specific random number sequence from a previous simulation.

* You want to set `ixseed` and `iyseed` to some specific values:
 * In `nmlSystem` set `iseed`, `iyseed` and `ixseed` to the desired values and set `luseXYseed` to `.t.`.
 * Note that if you set `iseed` to a negative value, then the seed will not be set by the cpu time. A negative seed will reinitialize the random number generator and overwrite `ixseed` and `iyseed`.
* Given an output file where you want to recreate the random numbers you have three possibilities:
 * If you want to generate a new simulation with the same initial seed: set `iseed` to the value of the seed given at the beginning of the input file
 * If you want to generate a new simulation with the same initial seed of the previous simulation which was started by setting `iyseed` and `ixseed` explicitly: set `iseed`, `iyseed` and `ixseed` to the ones given in the beginning of the output file set `luseXYseed` in `nmlSystem` to `.t.`.
 * If you want to run a simulation, continuing with the same chain of random numbers where a the previous simulation left of: set `iseed`, `iyseed` and `ixseed` to the ones given in the end of the output file set `luseXYseed` in `nmlSystem` to `.t.`.

# Compilation Modes {#Compilation}

To make the development easier, several compilation modes are provided. They are used as following:
* when making Molsim a second argument can be massed to make: `mode=<mode>`. Example:

```shell
make ser mode=test
```
to compile Molsim in the test mode. The following modes are supported:
 * `normal`: compile with normal options to achieve high computational efficiency.
 * `test`: compile in such a way, that the results are reproducible. Use this mode when performing tests in the Testin directory.
 * `debug`: During run-time the compiler checks for out-of-bounds arrays and uninitialized variables. Only use for debug purposes.
 * `quick`: Fast compilation without optimization. Use only for testing purposes.
 * `warn`: The compiler gives all warnings during compile time.
 * `gprof`: For profiling the execution with gprof (only for the `gfortran` compiler)

## Gprof
To profile you program with [`gprof`](https://sourceware.org/binutils/docs/gprof/), do the following (only for the `gfortran` compiler):
* Compile with `mode=gprof`:

```shell
make ser mode=gprof
```
* run Molsim on you project in the directory of you input file. A file `gmon.out` should be generated

```shell
molsim_ser <Projectname>
```
* run `gprof`

```shell
gprof <path to the molsim exe> gmon.out
```
