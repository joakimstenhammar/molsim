Input
=====
All input data read from the file FIN are done with the namelist procedure. Each namelist contains input related variables. Further information on the namelist procedure is given in Appendix A.

The following namelists are available. Compulsory nameslist are checked with X.

| file           | namelist                         |      | variables describing                                    |
| :------------- | :-------------------             | :--- | :-----------------------------------------------------  |
| molsim.F90     | \subpage nmlSystem               | X    | general system variables                                |
|                | \subpage nmlScale                | X    | scaling units                                           |
|                | \subpage nmlThermoInteg          |      | thermodynamic integration                               |
|                | \subpage nmlDist                 |      | distribution functions                                  |
| particle.F90   | \subpage nmlParticle             | X    | particle variables                                      |
|                | \subpage nmlCopolymerSequence    |      | describe the sequence of copolymers                     |
|                | \subpage nmlNetworkConfiguration |      | describe network configuration                          |
|                | \subpage nmlRepeating            |      | define the repeating block structure of copolymers      |
| potential.F90  | \subpage nmlPotential            | X    | potentials and forces                                   |
|                | \subpage nmlPotentialChain       |      | bond and angle potential                                |
|                | \subpage nmlPotentialExternal    |      | external potential                                      |
|                | \subpage nmlPolarizationIter     |      | many-body polarization calculation                      |
| coordinate.F90 | \subpage nmlSetConfiguration     | X    | initial start configuration                             |
| md.F90         | \subpage nmlMD                   |      | molecular dynamics simulation                           |
| mc.F90         | \subpage nmlMC                   |      | Monte Carlo simulation                                  |
|                | \subpage nmlMCAll                |      | Monte Carlo (all) simulation                            |
|                | \subpage nmlMCWeight             |      | applying weights for pmf calculations                   |
|                | \subpage nmlUmbrella             |      | umbrella potential sampling                             |
|                | \subpage nmlMCPmf                |      | direct pmf calculation using updated weights            |
| sso.F90        | \subpage nmlSPartSSO             |      | optimization of  the single particle move               |
| bd.F90         | \subpage nmlBD                   |      | Brownian dynamics simulation                            |
| nlist.F90      | \subpage nmlIntList              |      | calculation of list of nonbonded pairs                  |
| dump.F90       | \subpage nmlDump                 |      | dumping/reading of data                                 |
| group.F90      | \subpage nmlGroup                |      | dividing particles into groups                          |
| static.F90     | \subpage nmlStatic               |      | call of static analysis routines                        |
|                | \subpage nmlSPDF                 |      | single particle distribution functions                  |
|                | \subpage nmlRDF                  |      | radial distribution functions                           |
|                | \subpage nmlRDFChain             |      | radial distribution functions for chains                |
|                | \subpage nmlRDFSph               |      | projected radial distribution functions                 |
|                | \subpage nmlRDFCond              |      | conditional radial distribution functions               |
|                | \subpage nmlG3Dist               |      | normalized triplet correlation functions                |
|                | \subpage nmlSF                   |      | structure factors                                       |
|                | \subpage nmlScatIntens           |      | scattering intensities                                  |
|                | \subpage nmlAngDF                |      | angular distribution functions                          |
|                | \subpage nmlAngExtDF             |      | angular distribution functions wrt external frame       |
|                | \subpage nmlOriDipDF             |      | orientation/dipole distribution functions               |
|                | \subpage nmlRadAngDF             |      | radial-angular 2d distribution functions                |
|                | \subpage nmlKirkwoodgk           |      | Kirkwood's Gk factor                                    |
|                | \subpage nmlOriPolDF             |      | orientaitional polarization distribution functions      |
|                | \subpage nmlNNHB                 |      | neighbors and hydrogen bonds                            |
|                | \subpage nmlNNDF                 |      | nearest neighbor distribution functions                 |
|                | \subpage nmlChainDF              |      | chain distribution functions                            |
|                | \subpage nmlChainTypeDF          |      | chain type distribution functions                       |
|                | \subpage nmlChainTypeExtDF       |      | chain type distribution functions, external frame       |
|                | \subpage nmlCBPC                 |      | probability of bead-particle contact                    |
|                | \subpage nmlLoopTailTrain        |      | loop, tail, and train characteristics                   |
|                | \subpage nmlCluster              |      | cluster size distribution functions                     |
|                | \subpage nmlMultipoleDF          |      | electrostatic multipole moment distribution functions   |
|                | \subpage nmlEnergyDF             |      | energy distribution functions                           |
|                | \subpage nmlWidom1               |      | chemical potentials using Widom's method                |
|                | \subpage nmlWidom2               |      | chemical potentials using Widom's method                |
|                | \subpage nmlMeanForce1           |      | mean force between two particles                        |
|                | \subpage nmlMeanForce2           |      | mean force between two particles                        |
|                | \subpage nmlPotMeanForce         |      | potential of mean force between two particles           |
|                | \subpage nmlSurfaceArea          |      | surface area exposed by all particles of one type       |
|                | \subpage nmlTrajectory           |      | write trajectory on unit FLIST                          |
|                | \subpage nmlSubStructureDF       |      | distribution function of substructures                  |
|                | \subpage nmlNetworkDF            |      | network distribution functions                          |
|                | \subpage nmlNetworkRadialDF      |      | radial network distribution functions                   |
| dynamic.F90    | \subpage nmlDynamic              |      | call of static analysis routines                        |
|                | \subpage nmlMSD                  |      | mean square displacement                                |
|                | \subpage nmlOriXTCF              |      | orientational tcf of particle ``x'``/``y'``/``z'``-axis |
|                | \subpage nmlLinVelTCF            |      | linear / angular velocity tcf                           |
|                | \subpage nmlForTCF               |      | force  and torque tcf                                   |
|                | \subpage nmlIDMTCF               |      | induced dipole moment tcf                               |
|                | \subpage nmlUtotTCF              |      | energy tcf                                              |
| image.F90      | \subpage nmlImage                |      | call of image data writing routines                     |
|                | \subpage nmlVRML                 |      | generation of VRML coordinate files                     |
|                | \subpage nmlVTF                  |      | generation of VTF coordinate files                      |
| mixed.F90      | \subpage nmlMixed                |      | general Molmix variables                                |
|                | \subpage nmlMixed1               |      | generation of potential energy curves                   |
|                | \subpage nmlMixed2               |      | calculation of potential energies on a lattice          |
|                | \subpage nmlMixed3               |      | calculation of global potential energy minimum          |
|                | \subpage nmlMixed4               |      | generation of random coordinates                        |
|                | \subpage nmlMixed5               |      | calculation of second virial coefficients               |
|                | \subpage nmlMixed6               |      | calculation of orientational averaged potential energy  |
| moluser.F90    | \subpage nmlComplexation         |      | analysis of interparticle complexation                  |
|                | \subpage nmlComplexDist          |      | calculation of complexation distribution functions      |

The following subsections contain all input variables that are used for reading data from file FIN. The
variables are grouped together and listed below their namelist name. The first line contains the name,
the type, the dimension, and, if any, the default value of the variable. The following lines briefly
explain the use of the variable. If the variable only can attain a limited number of values, these are
listed. If description about writing data is given, the output is made on file FOUT (if nothing else is
stated).

A practical point. Often there is a need of reading a variable conditionally, i.e., if some condition is
true the variable is assigned a value and used later. If the condition is false, the presence of the
variable in the input file does no harm; thus it has not to be deleted. In the latter case the input file
contains redundant data.
