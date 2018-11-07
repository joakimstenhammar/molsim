! ... 'version 6.4.7, Sep 18, 2015'

!************************************************************************
!************************************************************************
!**                                                                    **
!**  Copyright 1990 - 2014                                             **
!**                                                                    **
!**  Per Linse                                                         **
!**  Physical Chemistry                                                **
!**  Department of Chemistry                                           **
!**  Lund University                                                   **
!**  Sweden                                                            **
!**                                                                    **
!**  All rights reserved. The code may not be modified or              **
!**  redistributed without the written conscent of the copyright       **
!**  owner. The copyright owner does not take any responsibility       **
!**  for any error in the code or its documentation.                   **
!**                                                                    **
!************************************************************************
!************************************************************************

!************************************************************************
!> \page mol mol.F90
!!  **MolModule**
!! *main module for the Molsim software*
!************************************************************************

module MolModule

   use StatisticsModule
   use MeshModule
   use Random_Module, only: k4b

!define standard output unit
   use, intrinsic :: iso_fortran_env, only : ustdout=>output_unit !, &
                                          ! ustdin=>input_unit !, &
                                          ! ustderr=>error_unit


! ... parameters

   integer(4), parameter :: mngen  = 5                ! max number of generations
   integer(4), parameter :: mnapt = 64                ! max number of atoms belonging to one particle
   integer(4), parameter :: mnct  = 10                ! max number of chain types (used for input variables)
   integer(4), parameter :: mnpt  = 10                ! max number of particle types (used for input variables)
   integer(4), parameter :: mnat  = 10                ! max number of atom types (used for input variables)

   real(8),    parameter :: Zero        = 0.0d0
   real(8),    parameter :: One         = 1.0d0
   real(8),    parameter :: ThreeHalf   = 1.5d0
   real(8),    parameter :: Two         = 2.0d0
   real(8),    parameter :: Three       = 3.0d0
   real(8),    parameter :: Four        = 4.0d0
   real(8),    parameter :: Five        = 5.0d0
   real(8),    parameter :: Six         = 6.0d0
   real(8),    parameter :: Half        = One/Two
   real(8),    parameter :: Third       = One/Three
   real(8),    parameter :: Fourth      = One/Four
   real(8),    parameter :: Sixth       = One/Six
   real(8),    parameter :: SqTwo       = sqrt(Two)
   real(8),    parameter :: SqThree     = sqrt(Three)
   real(8),    parameter :: expmax      = 87.0d0     ! exp(expmax) < overflow is assumed (expmax is ca 700 for real(8))
   real(8),    parameter :: Pi          = 3.1415926535897932d0
   real(8),    parameter :: TwoPi       = Two*Pi
   real(8),    parameter :: FourPi      = Four*Pi
   real(8),    parameter :: FourPiThird = Four*Pi/Three
   real(8),    parameter :: RadToDeg    = 180.0d0/Pi
   real(8),    parameter :: DegToRad    = Pi/180.0d0
   real(8),    parameter :: Bohr        = 0.52917720859d0    ! Å
   real(8),    parameter :: AvNo        = 6.02214179d23      ! 1/mol
   real(8),    parameter :: Boltz       = 1.3806504d-23      ! J/K
   real(8),    parameter :: ECh         = 1.602176487d-19    ! C
   real(8),    parameter :: GasConstant = 8.314472d0         ! J/(mol K)
   real(8),    parameter :: Epsi0       = 8.854187817d-12    ! F/m
   real(8),    parameter :: CalToJ      = 3600.d0/860.d0
   real(8),    parameter :: AuTokJ      = ECh**2/(Four*Pi*Epsi0)/(Bohr*1d-10)*AvNo*1d-3

   character(1), parameter :: tab = char(9)

! ... enumeration constants

   integer(4), parameter :: iReadInput        = 1 ! initiate and read input data
   integer(4), parameter :: iWriteInput       = 2 ! prepare and write input data
   integer(4), parameter :: iBeforeSimulation = 3 ! initiation before simulation
   integer(4), parameter :: iBeforeMacrostep  = 4 ! initiation before a macrostep
   integer(4), parameter :: iSimulationStep   = 5 ! one simulation time step/pass
   integer(4), parameter :: iAfterMacrostep   = 6 ! averaging etc after a macrostep
   integer(4), parameter :: iAfterSimulation  = 7 ! averaging etc after the simulation

! ... data structure for potential energy
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type potenergy_var
      real(8), allocatable :: twob(:)       ! two-body contribution (excluding ewald contribution)
      real(8), allocatable :: oneb(:)       ! one-body contribution (dielectric discontinuity)
      real(8)              :: tot           ! total
      real(8)              :: rec           ! reciprocal electrostatic space contribution (UEwald)
      real(8)              :: stat          ! static electrostatic contribution (umbodyp)
      real(8)              :: pol           ! polarization electrostatic contribution (umbodyp)
!> \page bond
!! `bond_var(real, integer, real)`(1:\ref nct )
!! **default:** \ref nct*`bond_var(0.0, 2, 0.0)`
!! * Force constant, power, and equilibrium separation of bond potential.
      real(8)              :: bond
!> \page angle
!! `bond_var(real, integer, real)`(1:\ref nct )
!! **default:** \ref nct `*bond_var(0.0, 2, 0.0)`
!! * Force constant, power, and equilibrium separation of angular potential.
      real(8)              :: angle
!> \page clink
!! `bond_var(real, integer, real)`
!! **default:** `bond_var(0.0, 2, 0.0)`
!! * Force constant, power, and equilibrium separation of crosslinks.
      real(8)              :: crosslink
      real(8)              :: external      ! external contribution
   end type potenergy_var

! ... data structure for bonds etc: These are documented in the manual in Chapter 7 (file datastructures.md)

   type bond_var
      real(8)      :: k                     ! force constant
      integer(4)   :: p                     ! power
      real(8)      :: eq                    ! equilibrium value
   end type bond_var

! ... data structure for asorption conditions: These are documented in the manual in Chapter 7 (file datastructures.md)

   type adscond_var
      character(8) :: txobject              ! "plane_xy": lower and upper xy-plane of a parallelepiped
                                            ! "part=xxx": surface of a particle xxx consitutes the adsorbing surface
      character(3) :: txsurface             ! "-", "+", "-|+" for lower, upper, and both surfaces
      real(8)      :: dist                  ! threshold distance of adsorption for the com of a particle from plane or center of particle
   end type adscond_var

! ... data structure for chain properties: These are documented in the manual in Chapter 7 (file datastructures.md)

   type chainprop_var
      real(8)    :: ro(3)                   ! center of mass
      real(8)    :: rbb2                    ! bead-to-bead distance squared
      real(8)    :: angle                   ! angle between consecutive beads
      real(8)    :: cos                     ! cos(180 - angle between consecutive beads)
      real(8)    :: ree(3)                  ! end-to-end vector
      real(8)    :: ree2                    ! end-to-end distance squared
      real(8)    :: rg2                     ! radius of gyration squared
      real(8)    :: rg2s                    ! square extention along principal axes (smallest)
      real(8)    :: rg2m                    ! square extention along principal axes (middle)
      real(8)    :: rg2l                    ! square extention along principal axes (largest)
      real(8)    :: rg2z                    ! radius of gyration squared projected on the z-axis
      real(8)    :: rg2xy                   ! radius of gyration squared projected on the xy-plane
      real(8)    :: lpree                   ! persistence length based on end-to-end separation
      real(8)    :: lprg                    ! persistence length based on radius of gyration
      real(8)    :: shape                   ! square end-to-end distance / square of radius of gyration
      real(8)    :: asph                    ! asphericity (JCP 100, 636 (1994))
      real(8)    :: torp                    ! toroidicity
   end type chainprop_var

! ... data structure for network (finite) properties: These are documented in the manual in Chapter 7 (file datastructures.md)

   type networkprop_var
      real(8)    :: ro(3)                   ! center of mass
      real(8)    :: rg2                     ! radius of gyration squared
      real(8)    :: rg2x                    ! radius of gyration squared projected on the x-axis
      real(8)    :: rg2y                    ! radius of gyration squared projected on the y-axis
      real(8)    :: rg2z                    ! radius of gyration squared projected on the z-axis
      real(8)    :: rg2s                    ! square extention along principal axes (smallest)
      real(8)    :: rg2m                    ! square extention along principal axes (middle)
      real(8)    :: rg2l                    ! square extention along principal axes (largest)
      real(8)    :: eivr(3,3)               ! normalized eigenvectors of the principal frame
      real(8)    :: theta(3)                ! angles of axes of largest extension and x-, y-, and z-axes of main frame
      real(8)    :: asph                    ! asphericity (JCP 100, 636 (1994))
      real(8)    :: alpha                   ! degree of ionization (for titrating systems)
   end type networkprop_var

! ... data structure for simple averaging: These are documented in the manual in Chapter 7 (file datastructures.md)

   type aver_var
      real(8)      :: s2                    ! summation/averaging over steps
      real(8)      :: s1                    ! summation/averaging over macrosteps
   end type aver_var

! ... data structure for 1D static variables: These are documented in the manual in Chapter 7 (file datastructures.md)

   type static1D_var
      logical      :: l                     ! logical flag for engeagement
      real(8)      :: min                   ! minimum value
      real(8)      :: max                   ! maximum value
      integer(4)   :: nbin                  ! number of bins
      logical      :: lnorm                 ! logical flag for normalization
      character(40):: label                 ! title
      integer(4)   :: nvar                  ! expanded into nvar variables
   end type static1D_var

! ... data structure for 2D static variables: These are documented in the manual in Chapter 7 (file datastructures.md)

   type static2D_var
      logical      :: l                     ! logical flag for engeagement
      real(8)      :: min(2)                ! minimum value
      real(8)      :: max(2)                ! maximum value
      integer(4)   :: nbin(2)               ! number of bins
      logical      :: lnorm                 ! logical flag for normalization
      character(40):: label                 ! title
      integer(4)   :: nvar                  ! expanded into nvar variables
   end type static2D_var

! ... version, date and author

   character(29) :: txVersionDate = 'version 6.4.7, v4.5.3'
   character(9)  :: txAuthor      = 'Per Linse'

! ... external units

   character(240) :: project ! name of the project (15 characters shorter than the filenames to have enough space for suffixes
   ! filenames
   character(255) :: fin     ! input data
   character(255) :: fout    ! output data
   character(255) :: fcnf    ! configuration data

#define STRINGIFY(x) x
   character(255) :: flib = trim(adjustl(STRINGIFY(FLIBMACRO)))

   character(255) :: flist   ! list data
   character(255) :: fuser   ! user-provided data
   character(255) :: fwrl    ! image data in wrl format
   character(255) :: fvtf    ! VMD: vtf data
   character(255) :: ftcl    ! VMD: tcl script
   character(255) :: fgroup  ! group data
   character(255) :: fpos    ! position data
   character(255) :: fori    ! orientation data
   character(255) :: fliv    ! linear velocity data
   character(255) :: fanv    ! angular velocity data
   character(255) :: ffor    ! force data
   character(255) :: ftor    ! torque data
   character(255) :: fidm    ! induced dipole moment data
   character(255) :: flaz    ! atom charge status
   character(255) :: futot   ! potential energy data

   integer(4),    parameter :: uin   = 1   ! input data
   integer(4),    parameter :: uout  = 2   ! output data
   integer(4),    parameter :: ucnf  = 3   ! configuration data
   integer(4),    parameter :: ulib  = 4
   integer(4),    parameter :: ulist = 8   ! list data
   integer(4),    parameter :: uuser = 9   ! user-provided data
   integer(4),    parameter :: uwrl  = 10  ! image data
   integer(4),    parameter :: uvtf  = 11  ! VMD: vtf data
   integer(4),    parameter :: utcl  = 12  ! VMD: tcl script
   integer(4),    parameter :: ugroup= 15  ! group data
   integer(4),    parameter :: upos  = 20  ! position data
   integer(4),    parameter :: uori  = 21  ! orientation data
   integer(4),    parameter :: uliv  = 22  ! linear velocity data
   integer(4),    parameter :: uanv  = 23  ! angular velocity data
   integer(4),    parameter :: ufor  = 24  ! force data
   integer(4),    parameter :: utor  = 25  ! torque data
   integer(4),    parameter :: uidm  = 26  ! induced dipole moment data
   integer(4),    parameter :: ulaz  = 27  ! atom charge status
   integer(4),    parameter :: uutot = 28  ! potential energy data

! ... terms
!> \page txtitle
!! `character(90)`
!! * User-provided title.
   character(90) :: txtitle

!> \page txmode
!! `character(10)`
!!  **default:** '`simulation`'
!!  * '`simulation`': Simulation and analyses. Further specification of method, ensemble, boundary conditions, and initial
!!  configuration are specified by variables \ref txmethod, \ref txensemb, \ref txbc, and \ref txstart. On a top level, analyses are controlled
!!  by \ref lgroup , \ref lstatic , and \ref ldynamic.
!!  * '`analysis`': Reading of dump data from DUMP files and analyses. On a top level, analyses are controlled by \ref lgroup ,
!!  \ref lstatic , and \ref ldynamic.
!!  * '`mixed`': Mixed activities controlled by namelist \ref nmlMixed.
   character(10) :: txmode

!> \page txmethod
!! `character(5)`
!! * '`md`': Molecular dynamic simulation. Further specification is given in namelist \ref nmlMD.
!! * '`mc`': Monte Carlo simulation. Further specification is given in namelist \ref nmlMC.
!! * '`mcall`': Monte Carlo with simultaneous movement of all particles. Further specification is given in namelist \ref nmlMCAll.
!! * '`bd`': Brownian dynamic simulation (configuration space). Further specificati on is given in namelist \ref nmlBD.
   character(5)  :: txmethod

!> \page txensemb
!! `character(3)`
!! * `nve`: Microcanonical ensemble (only MD). This option should also be used for MD with temperature and/or volume scaling.
!! Further specification is given in namelist \ref nmlMD.
!! * '`ntv`': Canonical ensemble (only MC or BD).
!! * '`nvt`': Canonical ensemble (only MC or BD).
!! * '\ref npt': Isothermal-isobaric ensemble (only MC).
!! * '`ntp`': Isothermal-isobaric ensemble (only MC).
!! * '`mvt`': Grand canonical ensemble (only MC). Still not fully tested.
!! * '`mtv`': Grand canonical ensemble (only MC). Still not fully tested.
   character(3)  :: txensemb
!> \page txstart
!! `character(8)`
!! * '`setconf`' : Generation of a start configuration and accumulation variables are set to zero. This option is used to obtain a
!!                 start configuration that should be equilibrated. Further specification is given in namelist \ref nmlSetConfiguration.
!! * '`readfin`' : Read start configuration from file FIN and accumulation variables are set to zero. The format
!!                 is(ro(1:3,ip),ori(1:3,1:3,ip),ip = 1,np), i.e., the same as the output of particle coordinates on FOUT.
!! * '`zero`' : Start configuration is read from FCNF and accumulation variables are set to zero. This option is used to start of a production run composed of \ref nstep1*\ref nstep2 steps/passes.
!! * '`continue`' : Start configuration and accumulated averages are read from FCNF. This option is used to continue an
!!                  equilibration or production run if the execution was aborted due to exceeded time limit, system malfunction etc. To continue,
!!                  change start to 'continue' and resubmit the job. This option may also be used to extend a completed production run. Then increase
!!                  \ref nstep1 to the new total number of macrosteps, change start to 'continue', and resubmit the job. The simulation will continue from
!!                  the last run to a new total \ref nstep1* \ref nstep2 steps/passes.
   character(8)  :: txstart
!> \page lcont
!! `logical`
!! **default:** `.false.`
!! * Quantities which use is primary to check that the simulation is properly advancing may be calculated (by driver ControlAver).
!! The quantities are averaged over particle types and they are average square force and torque (per particle), linear and angular
!! moments (per particle), as well as translational and rotational temperatures (only MD). Fraction accepted and rejected attempts to
!! move (only MC). Mean square displacement per step/pass. Orientational order. Defined as the scalar product of the direction of a
!! molecular axis and its direction at start. The value is one for a perfect alignment and approximately zero for a fluid. Useful for
!! monitoring the relaxation of an initial lattice start. Position and orientational means. The position of the center of mass denoted
!! \f$ \langle r0\rangle \f$, and the projection of the molecular axes on the box axes denoted ``<x'>``, ``<y'>``, and ``<z'>``.
!! * `.true.`: Control quantities are calculated.
!! * `.false.`: No calculation of control quantities.
   logical       :: lcont
!> \page laver
!! `logical`
!! **default:** `.false.`
!! * Thermodynamic averages and their precision as well as fluctuations and their precision may be calculated (by routine MainAver).
!! Consider the quantity Q. The precision its average \f$ \langle Q\rangle\f$ and fluctuation \f$ \sqrt{\langle Q^2\rangle -
!! \langle Q\rangle ^2} \f$, both given as one standard deviation are evaluated by block averaging with variable block length and extrapolation to infinite block length. The quantities considered are:
!!
!! 1. Total energy (MD)
!! 2. Kinetic energy (MD)
!! 3. Potential energy, total
!! 4. Potential energy, total two-body contribution (only \ref txelec='dip' or 'pol' the electrostatic interaction is excluded)
!! Potential energy, partial two-body contributions (only \ref txelec='dip' or 'pol' the electrostatic interaction is excluded)
!! 5. Potential energy from the reciprocal space (only \ref txelec='charge' .and. \ref lewald)
!! 6. Electrostatic energy (only \ref txelec='dip' or 'pol'), including the reciprocal space (only \ref lewald)
!! 7. Polarization energy (only \ref txelec='pol'), including the reciprocal space (only \ref lewald)
!! 8. Enthalpy
!! 9. Heat capacity
!! 10. Temperature
!! 11. Pressure
!! 12. Volume
!! 13. Induced dipole moment: total (only \ref txelec='pol')
!! 14. Induced dipole moment: particle (only \ref txelec='pol')
!! 15. Induced dipole moment: atom (only \ref txelec='pol')
!! * `.true.`: Thermodynamic averages are calculated.
!! * `.false.`: No calculation of thermodynamic averages.
   logical       :: laver
!> \page lti
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Handle thermodynamic integration (by routine ThermoInteg). Further specification is given in namelist \ref nmlThermoInteg.
!! * `.false.`: No thermodynamic integration.
   logical       :: lti
!> \page ldist
!! `logical`
!! **default:** `.false.`
!! * `.true.` Distribution functions are calculated (by routine DistFunc). Further specification is given in namelist \ref nmlDist.
!! * `.false.` No calculation of distribution functions.
   logical       :: ldist
!> \page ldump
!! `logical`
!! **default:** `.false.`
!! * `.true.` Variables are dumped or read. Further specification is given in namelist \ref nmlDump.
!! * `.false.` No dumping/readiing
   logical       :: ldump
!> \page lgroup
!! `logical`
!! **default:** `.false.`
!! * `.true.` Particles are divided into groups. This division is used in routine static and optionally by moldyn. Further specification is given in namelist \ref nmlGroup.
!! * `.false.` No division into groups
   logical       :: lgroup
!> \page lstatic
!! `logical`
!! **default:** `.false.`
!! * `.true.` Makes it possible to call several routines calculating static properties. Further specification is given in namelist \ref nmlStatic (require \ref lgroup=.true.).
!! * `.false.` No static analysis.
   logical       :: lstatic
!> \page ldynamic
!! `logical`
!! **default:** `.false.`
!! * `.true.` Makes it possible to call several routines calculating dynamic properties. Further specification is given in namelist \ref nmlDynamic (require \ref lgroup=.true.).
!! * `.false.` No dynamic analysis.
   logical       :: ldynamic
!> \page limage
!! `logical`
!! **default:** `.false.`
!! * `.true.` Makes it possible to call several routines preparing files for generating images. Further specification is given in
!! namelist \ref nmlImage.
!! * `.false.` No image analysis.
   logical       :: limage
!> \page nmlSystem_itest itest
!! `integer`
!! **default:** `0`
!! * Flag for test output. This possibility is for maintenance purposes.
!! * `0`: Nothing. The normal option.
!! * '1': Energy, pressure, particle and atom coordinates, forces, torques, dipole moments, polarizability tensor, linear and
!!        angular accelerations, linear and angular velocities, and particle distance matrix are written after each step/configuration.
!! * `3`: Intermediate energies (energy.F90).
!! * `4`: Neighbour list (nlist.F90).
   integer(4)    :: itest
!> \page ltrace
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Information upon entering and exiting subroutines on three levels are written on unit 40 for the master and unit 41 for slaves.
!! * `.false.`: Nothing.
   logical       :: ltrace
!> \page lblockaver
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Intermediate data concerning the block averaging and its extrapolation to infinite block length are written on file blockaver.data.
!! * `.false`: Nothing.
   logical       :: lblockaver
!> \page txuser
!! `character(80)`
!! * Character string for engaging user-provided code in the MOLSIM core.
   character(80) :: txuser
!> \page ipart
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Initial and final particle coordinates and velocities are written for every \ref ipart particle.
   integer(4)    :: ipart
!> \page iatom
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Initial atomic coordinates are written for every \ref iatom atom.
   integer(4)    :: iatom
!> \page iaver
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Cumulative (for the macrostep) thermodynamic averages are written for every \ref iaver steps/passes.
   integer(4)    :: iaver
!> \page ishow
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Calculated functions are listed at every \ref ishow bin. Compact list.
   integer(4)    :: ishow
!> \page iplot
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Calculated functions are plotted.
   integer(4)    :: iplot
!> \page ilist
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `>0`: Calculated functions are listed at every \ref ilist bin on file FLIST. Extended list.
   integer(4)    :: ilist
!> \page idump
!! `integer`
!! **default:** `10`
!! * Dumping/reading is performed every \ref idump step/pass. The types of quantities dumped are controlled by the remaining variables.
!! If \ref txstart='continue' the dump files are properly positioned according to the step/configuration saved in file FCNF.
   integer(4)    :: idump                  ! if idump > 0, dump every dump:th step/pass
!> \page istatic
!! `integer`
!! **default:** `1`
!! * Conformation sampling interval. Holds both for \ref txmode = 'simulation' and 'analysis'.
   integer(4)    :: istatic                ! if istatic > 0, call static analysis routines every istatic step/pass:

   logical       :: lsim                   !=.true. simulation
   logical       :: lana                   !=.true. analysis
   logical       :: lmix                   !=.true. mixed tasks
!> \page ltime
!! `logical`
!! **default:** `.true.`
!! * `.true.`: Timing statistics are carried out.
!! * `.false.`: No timing statistics.
   logical       :: ltime

   logical       :: lnve                   !=.true. if NVE ensemble
   logical       :: lnvt                   !=.true. if NVT ensemble
   logical       :: lntp                   !=.true. if NTP ensemble
   logical       :: lmvt                   !=.true. if MVT ensemble
   real(8)       :: facmvt = 4.0           ! factorial increase of memory allocation, applied on np and na (lmvt)

! ... cell
!> \page txbc
!! `character(3)`
!! * '`xyz`': Periodical boundary condition, x, y, and z-directions.
!! * '`xy`': Periodical boundary condition, x and y-directions.
!! * '`z`': Periodical boundary condition, z-direction.
!! * '`rd`': Periodical boundary condition, rhombic dodecahadral.
!! * '`to`': Periodical boundary condition, truncated octahedral.
!! * '`sph`': Hard sphere boundary (only MC).
!! * '`cyl`': Hard cylinder boundary (only MC).
!! * '`ell`': Hard ellipsoidal boundary (only MC).
   character(3)  :: txbc                   ! control boundary conditions
   logical       :: lbcbox                 ! box-like cell (rätblock)
   logical       :: lbcrd                  ! rhombic dodecahedral cell
   logical       :: lbcto                  ! truncated octahedral cell
   logical       :: lbcsph                 ! spherical cell
   logical       :: lbccyl                 ! cylindrical cell
   logical       :: lbcell                 ! ellipsoidal cell
!> \page boxlen
!! `real`(1:3)
!! * Box length in x-, y-, and z-directions.
   real(8)       :: boxlen(3)
   real(8)       :: boxlen2(3)             ! boxlen/2
   real(8)       :: boxleni(3)             ! 1/boxlen
   real(8)       :: boxlen2i(3)            ! 1/(2*boxlen)
   real(8)       :: boxlenshort            ! minval(boxlen(1:3))
   real(8)       :: TwoPiBoxi(3)           ! 2*pi/boxlen
!> \page cellside
!! `real`
!! * Side length of a rhobic dodecaheron (only \ref txbc='rd') and truncated octahedron cell (only \ref txbc='to').
   real(8)       :: cellside
!> \page sphrad
!! `real`
!! * Spherical cell radius (only \ref txbc='sph').
   real(8)       :: sphrad
   real(8)       :: sphrad2                ! sphrad**2
!> \page ellrad
!! `real`(1:3)
!! * Ellipsoidal cell radii (only \ref txbc='ell').
   real(8)       :: ellrad(3)
   real(8)       :: ellradi(3)             ! 1/ellrad
!> \page cylrad
!! `real`
!! * Cylindrical cell radius (only \ref txbc='cyl').
   real(8)       :: cylrad
   real(8)       :: cylrad2                ! cylrad**2
!> \page cyllen
!! `real`
!! * Cylindrical cell length (only \ref txbc='cyl').
   real(8)       :: cyllen = 0.0d0
   real(8)       :: cyllen2                ! length of cylindrical cell/2
!> \page lenscl
!! `real`
!! **default:** `1.0`
!! * Scaling constant which scales box, and particle coordinates (only \ref txstart='zero'). Useful for creating a system similar to a previous one, but with different box lengths.
   real(8)       :: lenscl
   logical       :: lPBC                   ! periodic boundary conditions
   real(8)       :: dpbc(3)                ! =boxlen for some pbc, otherwise zero

! ... scaling
!> \page scllen
!! `real`
!! **default:** `1.0d-10` m
!! * Scaling factor for length.
    real(8)      :: scllen
!> \page sclmas
!! `real`
!! **default:** `1.0d-3` kg/mol
!! * Scaling factor for mass.
    real(8)      :: sclmas
!> \page scltem
!! `real`
!! **default:** `1.0` K
!! * Scaling factor for temperature
    real(8)      :: scltem
!> \page scltim
!! `real`
!! **default:** `1.0d-12` s
!! * Scaling factor for time.
    real(8)      :: scltim
!> \page sclene
!! `real`
!! **default:** `1.0d+3`J/mol
!! * Scaling factor for energy.
    real(8)      :: sclene
!> \page sclhca
!! `real`
!! **default:** `1.0` J/(mol*K)
!! * Scaling factor for heat capacity.
    real(8)      :: sclhca
!> \page sclpre
!! `real`
!! **default:** `1.0d+6` Pa
!! * Scaling factor for pressure.
    real(8)      :: sclpre
!> \page scldif
!! `real`
!! **default:** `1.0d-9` \f$m^2\f$/s
!! * Scaling factor for diffusion coefficient.
    real(8)      :: scldif
    real(8)      :: sclvol                 ! volume
    real(8)      :: sclvel                 ! linear velocity
    real(8)      :: sclacc                 ! linear acceleration
    real(8)      :: sclfor                 ! force
    real(8)      :: sclmin                 ! moment of intertia
!> \page sclang
!! `real`
!! **default:** `PI/180` 1/rad
!! * Scaling factor for angle.
    real(8)      :: sclang
    real(8)      :: beta                   ! 1/(kT)
    real(8)      :: Epsi0FourPi            ! ECh**2/(Four*pi*Epsi0*scllen)*AvNo/sclene
    real(8)      :: EpsiFourPi             ! Epsi0FourPi/relpermitt

! ... particle
!> \page txelec
!! `character(11)`
!! **default:** '`charge`'
!! * Describes the complexity of the electrostatic interactions. Available sources of electrostatic interactions are: (permanent)
!! charges, weak charges, dipoles, and induced dipoles. Charges, dipole moments, and polarizability tensors employed are specified in
!! namelist \ref nmlParticle.
!! * '`charge`': Atoms possessing charges only.
!! * '`weakcharge`': Atoms possessing charges and weak charges only. Limited to hard-sphere monoatomic particles.
!! * '`dip`': Atoms possessing charges and static dipoles only. MD for all boundary conditions, and MC for all except Ewald summation.
!! * '`pol`': Atoms possessing charges, static dipoles, and induced dipoles only. MD or MCALL. Ref. Chem. Phys. 191, 195 (1995).
!!            Charges, dipole moments, and polarizability tensors employed are specified in namelist \ref nmlParticle.
!! * '`dipsph`': Atoms possessing charges and/or dipoles in a spherical geometry, no neighbour list; special energy routines.
!! * '`dieldis`': Atoms possessing charges with the presence of a planar or spherical dielectric discontinuity.
   character(11) :: txelec                 !*"charge", "weakcharge", "dip", "pol", "dipsph", "dipsphimage", "dieldis"
   logical       :: lcharge                ! enable atomic charge interactions
   logical       :: lweakcharge            ! enable atomic charge interactions, weak charges
                                           ! 1) no explicit counterions: single particle type, several titrable atom types
                                           ! 2) explicit counterions: several particle types, at most one titrable atom type per particle
   logical       :: ldipole                ! enable atomic charge and dipole interactions
   logical       :: lpolarization          ! enable atomic charge, dipole, and polarization interactions
   logical       :: ldipolesph             ! enable atomic charge and dipole interactions, spherical geometry, explicit energy evaluation
   logical       :: laimage                ! enable atomic image charge and dipole interactions
   logical       :: ldieldis               ! enable atomic charge and dielectric discontinuities
   integer(4)    :: lq2sum                 ! =0 system contains charges
                                           ! =1 system contains dipoles and no charges
                                           ! =-1 else
   real(8)                 :: q2sum        ! sum(q**2), where q is either atomic charge or dipole according to lq2sum
!> \page nnwt
!! `integer`
!! **default:** `0`
!! * Number of network types (only \ref txsetconf = 'network').
   integer(4)              :: nnwt         !*number of network types
!> \page nct
!! `integer`
!! **default:** `0`
!! * Number of chain types.
   integer(4)              :: nct
!> \page npt
!! `integer`
!! * Number of particle types.
   integer(4)              :: npt
   integer(4)              :: nat          ! number of atom types
   integer(4)              :: nnw          ! number of networks
   integer(4)              :: nc           ! number of chains
   integer(4)              :: np           ! number of particles
   integer(4)              :: na           ! number of atoms
   integer(4)              :: na_alloc     ! number of atoms (for memory allocation)
   integer(4)              :: np_alloc     ! number of particles (for memory allocation)
!> \page nnwnwt
!! `integer`(1:\ref nnwt)
!! **default:** \ref nnwt*`0`
!! * Number of networks of network type inwt (only \ref txsetconf = 'network').
   integer(4), allocatable :: nnwnwt(:)
!> \page ncct
!! `integer`(1:\ref nct)
!! * Number of chains of each chain type.
   integer(4)              :: ncct(mnct)
!> \page nppt
!! `integer`(1:\ref npt)
!! * Number of particles of each particle type.
   integer(4)              :: nppt(mnpt)
   integer(4), allocatable :: nctnwt(:)    ! number of chain types belonging to one network type [allocate with nnwt]
   integer(4)              :: nptct(mnct)  ! number of particle types belonging to one chain type
   integer(4), allocatable :: ncnwt(:)     ! number of chains belonging to a network type [allocate with nnwt]
   integer(4), allocatable :: npnwt(:)     ! number of particles belonging to a network type [allocate with nwwt]
   integer(4), allocatable :: nclnwt(:)    ! number of cross-links belonging to a network type [allocate with nnwt]
   integer(4)              :: npct(mnct)   ! number of particles belonging to a chain type
!> \page natpt
!! `integer`(1:\ref npt)
!! * Number of atom types of each particle type.
   integer(4)              :: natpt(mnpt)  !*number of atom types belonging to a particle type
   integer(4), allocatable :: napt(:)      ! number of atoms belonging to a particle type
   integer(4), allocatable :: naat(:)      ! number of atoms of an atom type in one particle
!> \page ncctnwt
!! `integer`(1:\ref nct,1:\ref nnwt)
!! **default:** \ref nct *\ref nnwt *`0`
!! * Number of chains of chain type ict of network type inwt (only \ref txsetconf = ‘network’).
   integer(4), allocatable :: ncctnwt(:,:)
!> \page npptct
!! `integer`(1:\ref npt ,1:\ref nct )
!! * \ref npptct (ipt,ict) is the number of particles of type ipt belonging to chain of type ict.  Note, either non or all particles of a
!!   given type has to belong to chains, in the latter case sum(\ref ncct (1:\ref nct )*\ref npptct (ipt,1:\ref nct )=\ref nppt (ipt) is required.
   integer(4)    :: npptct(mnpt,mnct)      !*number of particles of a particle type of chain type
!> \page naatpt
!! `integer`(1:nat,1:\ref npt )
!! * \ref naatpt (ialoc,ipt) denotes the no of atoms of type ialoc (local list) on a particle of type ipt. ialoc runs from 1 to \ref natpt (ipt), the no of atom types of particle type ipt.
   integer(4)    :: naatpt(mnat,mnpt)
   integer(4)    :: nnwtnwt                ! number of different network type pairs
   integer(4)    :: nctct                  ! number of different chain type pairs
   integer(4)    :: nptpt                  ! number of different particle type pairs
   integer(4)    :: natat                  ! number of different atom type pairs
!> \page txtoponwt
!! `character(30)`(1:\ref nnwt)
!! **default:** \ref nnwt *'`default`'
!! * Controls the network topology to be set (currently only ‘default’ possible and only \ref txsetconf = ‘network’).
   character(30), allocatable :: txtoponwt(:)
!> \page txcopolymer
!! `character(30)`(1:\ref nct)
!! **default:** \ref nct *'`block`'
!! * '`block`': Block copolymer.
!! * '`regular`': Regular copolymer (alternating as possible).
!! * '`sequence`': Copolymer with highly specific monomer distribution (the control is given in \ref nmlCopolymerSequence).
!! * '`repeating`': Copolymers with a defined repeating block structure (the control is given in \ref nmlRepeating).
!! * '`random`': Random Copolymer.
   character(30) :: txcopolymer(mnct)
   logical       :: lspma                  !*control tacticity
!> \page txnwt
!! `character(10)`(1:\ref nnwt)
!! **default:** \ref nnwt *`'network'`
!! * Text label for each network type inwt (only \ref txsetconf = ‘network’).
   character(10), allocatable :: txnwt(:)
!> \page txct
!! `character(10)`(1:\ref nct)
!! * Text label for each chain type.
   character(10) :: txct(mnct)
!> \page txpt
!! `character(10)`(1:\ref npt)
!! * Text label for each particle type.
   character(10) :: txpt(mnpt)
!> \page txat
!! `character(10)`(1:nat)
!! * Text label for each atom type (global atom type list).
   character(10) :: txat(mnat)
   character(21), allocatable :: txnwtnwt(:) ! name of network type-network type pair
   character(21), allocatable :: txctct(:) ! name of chain type-chain type pair
   character(21), allocatable :: txptpt(:) ! name of particle type-particle type pair
   character(21), allocatable :: txatat(:) ! name of atom type-atom type pair
!> \page lmonoatom
!! `logical`
!! **default:** `.true.`
!! * `.true.`: The program checks to see if all particles have only one atom each. If so, 1) sections involving orientational
!!            integration, orientational movement, and some orientational output are omitted and 2) simpler and faster potential and force
!!            routines are used.
!! * `.false.`: Forces the program to treat the particles as polyatomic. This option, which may make the program slower, is for maintaining and checking purposes.
   logical                    :: lmonoatom
   logical                    :: lpolyatom     !
!> \page massat
!! `real`(1:nat)
!! * Mass of atom type.
   real(8)                    :: massat(mnat)
   real(8), allocatable       :: masspt(:)     ! mass of particle type
   real(8), allocatable       :: massnwt(:)    ! mass of network type
   real(8), allocatable       :: massipt(:)    ! inverse mass of particle type
   real(8), allocatable       :: massinwt(:)   ! inverse mass of network type
   real(8), allocatable       :: mompt(:,:)    ! moments of inertia of particle type
   real(8), allocatable       :: momipt(:,:)   ! inverse moments of inertia of particle type
   real(8), allocatable       :: massp(:)      ! mass of particle
   real(8), allocatable       :: massip(:)     ! inverse mass of particle
   real(8), allocatable       :: momp(:,:)     ! moments of inertia of particle
   real(8), allocatable       :: momip(:,:)    ! inverse moments of inertia of particle
!> \page radat
!! `real`(1:nat)
!! * Hard-sphere radius of atom type. It is used to avoid an unreasonable start configuration (see namelist \ref nmlSetConfiguration
!!  ) and to check that atoms do not come too close to each other due to potential or program error. The check is performed after each
!!  macrostep by routine  CheckHSOverlap . The value of \ref radat should be smaller than the van der Waals radius, approximately 75% of it.
!!  In the case of a hard-core potential and MC, \ref radat should be equal to the hard-core radius of the atom.
   real(8)                    :: radat(mnat)
   integer(4), allocatable    :: itradegfree(:)! number of translational degrees of freedom of a particle type
   integer(4), allocatable    :: irotdegfree(:)! number of rotational degrees of freedom of a particle type
   real(8), allocatable       :: racom(:)      !
   real(8), allocatable       :: r1atat(:)     !
   real(8), allocatable       :: r2atat(:)     !
!> \page zat
!! `real`(1:nat)
!! **default:** nat*`0.0`
!! * Charge of atom type. The charge should be given in number of elementary charges.
   real(8)                    :: zat(mnat)
!> \page zatalpha
!! `real`(1:nat)
!! **default:** `0.0`
!! * >`0.0`: Gaussian charge distribution with width 1/(sqrt(2)*\ref zatalpha)
!! * =`0.0`: Point charge
   real(8)                    :: zatalpha(mnat)
!> \page sigat
!! `real`(1:nat)
!! **default:** nat*`0.0`
!! * Lennard-Jones \f$ \sigma \f$ parameter of atom type;\f$ u(r) = 4\epsilon [(\sigma/r)^{12}-(\sigma/r)^6] \f$.
   real(8)                    :: sigat(mnat)
!> \page epsat
!! `real`(1:nat)
!! **default:** (1:nat)
!! * Lennard-Jones \f$ \epsilon \f$ parameter of atom type;\f$ u(r) = 4\epsilon [(\sigma/r)^{12}-(\sigma/r)^6] \f$.
   real(8)                    :: epsat(mnat)   !*LJ epsilon parameter for atomes of a given type
   real(8), allocatable       :: az(:)         ! atom charge
   real(8), allocatable       :: aztm(:)       ! atom charge (trial move)
!> \page latweakcharge
!! `logical`(1:nat)
!! **default:** nat*`.false.`
!! * `.true.`: Weak (titrating) charge of atom type (only \ref txelec ='weakcharge').
!! * `.false.`: No weak charge.
   logical                    :: latweakcharge(mnat)
!> \page jatweakcharge
!! `integer`(1:nat)
!! **default:** nat*`0`
!! * type of atom carrying counter charge to weak charge iat (0 means no counter charge)
   integer(4)                 :: jatweakcharge(mnat)
!> \page pH
!! `real`
!! **default:** `0.0`
!! * pH of the solution.
   real(8)                    :: pH
!> \page pK
!! `real`(1:nat)
!! **default:** nat*`0.0`
!! * \ref pK value of the weak charge.
   real(8)                    :: pK(mnat)
   real(8)                    :: pHmpK(mnat)   ! pH - pKa
   integer(4), allocatable    :: iananweakcharge(:) ! atom carrying weak charge -> its atom carrying its couterion charge
   logical, allocatable       :: laz(:)        ! .true. if atom is charged
   logical, allocatable       :: laztm(:)      ! .true. if atom is charged (trial move)
   logical, allocatable       :: lpset(:)      ! logical valiable being .true. for setted paricles
!> \page lradatbox
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Check that atoms including their hard-sphere radii are inside the box.
!! * `.false.`: No such check.
   logical                    :: lradatbox     ! use atom and atom radius when checking if particle inside box

                                           ! equivalences:
                                           ! nnw          = sum(nnwnwt(1:nnwt))
                                           ! nc           = sum(ncct(1:nct))
                                           ! np           = sum(nppt(1:npt))
                                           ! na           = sum(napt(1:npt)*nppt(1:npt))
                                           ! nat          = sum(natpt(1:npt))
                                           ! nctnwt(inwt) = count(ncctnwt(1:nct,inwt) > 0)
                                           ! nptct(ict)   = count(npptct(1:npt,ict) > 0)
                                           ! ncnwt(inwt)  = sum(ncctnwt(1:nct,inwt))
                                           ! npct(ict)    = sum(npptct(1:npt,ict))
                                           ! napt(ipt)    = sum(naatpt(1:natpt(ipt),ipt)))
                                           ! naat(iat)    = naatpt(ka,iptat(iat))

! ... chain bond

   logical       :: lchain                 ! flag for connecting particles into chains by bonds
   logical       :: lfixedori              ! fix the particle orientation relative to the bond vectors, mirrorplane is xy plane
   integer(4), allocatable :: bondnn(:,:)  ! bond and particle             -> bonded particle

! ... crosslinks between chains
!> \page lclink
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Enabling cross-links between particles and chains or between chains. Diamond- like network and bottle-brushes are currently supported.
!! * `.false.`: No cross-links.
   logical       :: lclink
!> \page maxnbondcl
!! `integer`(1:\ref npt)
!! **default:** \ref npt*`4`
!! * Maximum number of cross-links of a particle (only \ref txsetconf='hierarchical').
   integer(4)    :: maxnbondcl(mnpt)
   integer(4)    :: ncl                    ! number of crosslinks
   integer(4), allocatable :: nbondcl(:)   ! actual number of crosslinks to/from particle ip
   integer(4)    :: maxvalnbondcl          ! maxval(nbondcl(:))
   integer(4), allocatable :: bondcl(:,:)  ! crosslink and particle        -> crosslinked particle
!> \page lmultigraft
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Enabling multigrafted polymers made of chains.
!! * `.false.`: No multigrafting.
   logical       :: lmultigraft

! ... hierarchical structure
!> \page lhierarchical
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Enabling hierarchical structures made of chains.
!! * `.false.`: No hierarchical structures.
   logical       :: lhierarchical
!> \page ngen
!! `integer`
!! **default:** `1`
!! * Number of generations of hierarchical structure (only \ref lhierarchical=.true. and \ref txsetconf ='hierarchical').
   integer(4)    :: ngen
   integer(4)    :: nh                     ! number of hierarchical structures
   integer(4)    :: nch(0:mngen)           ! number of chains of a generation in a hierarchial structure
   integer(4)    :: nphn                   ! number of particles of a hierarchical structure
!> \page ictgen
!! `integer`(1:\ref ngen)
!! **default:** `1`
!! * Generation number -> chain type (only \ref txsetconf='hierarchical').
   integer(4)    :: ictgen(0:mngen)        !*generation number             -> chain type
   integer(4), allocatable :: genic(:)     ! chain number                  -> generation number
!> \page nbranch
!! `integer`(0:\ref ngen -1)
!! **default:** `0`
!! * Number of branches of a branching point of generation igen (only \ref lhierarchical=.true. and \ref txsetconf = 'hierarchical').
   integer(4)    :: nbranch(0:mngen-1)
!> \page ibranchpbeg
!! `integer`(0:\ref ngen -1)
!! **default:** `0`
!! * Chain segment of first branching point of generation igen (only \ref lhierarchical=.true.  and \ref txsetconf = 'hierarchical').
   integer(4)    :: ibranchpbeg(0:mngen-1)
!> \page ibranchpinc
!! `integer`(0:\ref ngen -1)
!! **default:** `1`
!! * Segment increment between branching point of generation igen (only \ref lhierarchical=.true. and \ref txsetconf = 'hierarchical').
   integer(4)    :: ibranchpinc(0:mngen-1)

!  note: all chains of a given generation have to be of same type &
!        chains of a given type may appear in more than one generation
!        no further restriction on the nature of chains
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!           |          |          |                  | <--- ictgen(igen+1)
!           |          |          |                  |
!           |          |          |                  |
!           |          |          |        ...       |
!     ---------------------------------------------------------------- <--- ictgen(igen)
!     <----->          <---------->
!   ibranchpbeg(igen)       ibranchpinc(igen)
!
!           1          2          3        ...     nbranch(igen)
!
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
! Example 1: bottle-brush polymer
!
! first generation brush polymer: backbone containing 100 segments, 10 side chains each containing 5 segments
!
!     npct(1:2)       = 100, 5,
!     ngen            = 1,
!     ictgen(0)       = 1, 2,
!     nbranch(0)      = 10,
!     ibranchpbeg(0)  = 5,
!     ibranchpinc(0)  = 10,
!
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
! Example 2: dendrimer
!
!  second order generation bifunctional dendrimer: first generation chains have 5 and second 10 segments
!
!     npct(1:3)       = 1, 5, 10,
!     ngen            = 2,
!     ictgen(0)       = 1, 2, 3,
!     nbranch(0)      = 2, 2,
!     ibranchpbeg(0)  = 1, 5,
!     ibranchpinc(0)  = 0, 0,

! ... finite network structures

   logical                 :: lnetwork           ! flag for networks
   logical   , allocatable :: lptnwt(:,:)        ! true, if particle type ipt is part of network type inwt [allocate with npt,nnwt]
   logical   , allocatable :: lpnnwn(:,:)        ! true if particle number ip is part of network number inw
   integer(4), allocatable :: npweakchargenwt(:) ! number of titratable beads in network type inwt [allocate with nnwt]
!> \page iptclnwt
!! `integer`(1:\ref nnwt)
!! **default:** \ref nnwt*`0`
!! * Type of particles forming cross-links of networks of network type inwt (only \ref txsetconf = 'network').
   integer(4), allocatable :: iptclnwt(:)

! ... pointers

   integer(4), allocatable :: inwtnwn(:)   ! network (1:nnw)               -> its network type (1:nnwt)
   integer(4), allocatable :: inwtct(:)    ! chain type (1:nct)            -> its network type (1:nnwt)
   integer(4), allocatable :: inwtcn(:)    ! chain (1:nc)                  -> its network type (1:nnwt)
   integer(4), allocatable :: ictcn(:)     ! chain (1:nc)                  -> its chain type (1:nct)
   integer(4), allocatable :: inwncn(:)    ! chain (1:nc)                  -> its network (1:nnw)
   integer(4), allocatable :: ictpt(:)     ! particle type (1:npt)         -> its chain type (1:nct)
   integer(4), allocatable :: ihnpn(:)     ! particle (1:np)               -> its hierarchical structure (1:nh)
   integer(4), allocatable :: ictpn(:)     ! particle (1:np)               -> its chain type (1:nct)
   integer(4), allocatable :: icnpn(:)     ! particle (1:np)               -> its chain (1:nc)
   integer(4), allocatable :: iptpn(:)     ! particle (1:np)               -> its particle type (1:npt)
   integer(4), allocatable :: iptat(:)     ! atom type (1:nat)             -> its particle type (1:npt)
   integer(4), allocatable :: iptan(:)     ! atom (1:na)                   -> its particle type (1:npt)
   integer(4), allocatable :: ipnan(:)     ! atom (1:na)                   -> its particle (1:np)
   integer(4), allocatable :: iatan(:)     ! atom (1:na)                   -> its atom type (1:na)

   integer(4), allocatable :: inwnnwt(:)   ! network type (1:nnwt)         -> its first network (1:nnw)
   integer(4)              :: ipnhn        ! hierarchical structure        -> its first particle
   integer(4), allocatable :: icnct(:)     ! chain type (1:nct)            -> its first chain (1:nc)
   integer(4), allocatable :: ipnpt(:)     ! particle type (1:npt)         -> its first particle (1:np)
   integer(4), allocatable :: iatpt(:)     ! particle type (1:npt)         -> its first atom type (1:nat)
   integer(4), allocatable :: ianpn(:)     ! particle (1:np)               -> its first atom (1:na)
   integer(4), allocatable :: ianat(:)     ! atom type (1:nat)             -> its first atom (1:na)
   integer(4), allocatable :: inwtnwt(:,:) ! two network types (1:nnwt)    -> network type pair (1:nnwt)
   integer(4), allocatable :: ictct(:,:)   ! two chain types (1:nct)       -> chain type pair (1:nctct)
   integer(4), allocatable :: iptpt(:,:)   ! two particle types (1:npt)    -> particle type pair (1:nptpt)
   integer(4), allocatable :: iatat(:,:)   ! two atom types (1:nat)        -> atom type pair (1:natat)
   integer(4), allocatable :: isegpn(:)    ! particle (1:np)               -> segment number
   integer(4), allocatable :: ipnsegcn(:,:)! seg (1:npct) and chain (1:nc) -> particle (1:np)
   integer(4), allocatable :: icihigen(:,:)! hier (1:nh) and gen (0:ngen)  -> its first chain
   integer(4), allocatable :: icnclocnwn(:,:)  ! local chain (1:ncnwt) and network (1:nnw)    -> its chain (1:nc)
   integer(4), allocatable :: ipncllocnwn(:,:) ! cross-link (1:nclnwt) and network (1:nnw)    -> its particle (1:np)
   integer(4), allocatable :: ipnplocnwn(:,:)  ! local particle (1:npnwt) and network (1:nnw) -> its particle (1:np)

! ... md

   logical       :: lmd                    ! flag for molecular dynamics simulation

! ... mc

   logical       :: lmc                    ! flag for monte carlo simulation
   logical       :: lmcall                 ! flag for monte carlo simulation with all particle moves
   logical       :: lptmdutwob             ! flag for calulating dutobdy among moving particles
   logical, allocatable :: lptm(:)         ! flag for moving particles
   integer(4)    :: nptm                   ! number of moving particles
   integer(4)    :: natm                   ! number of moving atoms
   integer(4), allocatable :: ipnptm(:)    ! particle (local) (1:np)       -> particle (1:np)
   integer(4), allocatable :: iptmpn(:)    ! particle (1:np)               -> particle (local) (1:np)
   integer(4), allocatable :: ianatm(:)    ! atom (local) (1:na)           -> atom (1:na)
!> \page lmcweight
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Applying a weighting function for potential of mean force MC simulations. Further specification is given in namelist \ref nmlMCWeight .
!! * `.false.`:  No weighting function applied.
   logical       :: lmcweight              ! flag for using weighting function

! ... bd

   logical       :: lbd                    ! flag for Brownian dynamics simulation
!> \page dcoeff
!! `real`
!! * Isotropic particle self-diffusion coefficient.
   real(8), allocatable :: dcoeff(:)

! ... averages

   real(8)       :: tempst                 ! temperature, start
!> \page temp
!! `real`
!! * System temperature (only \ref txensemb='\ref npt' or 'nvt'). Desired temperature (only \ref txensemb='nve' ensemble and temperature scaling).
!!   Temp is also used if the velocities are set according to a Maxwell distribution.
   real(8)       :: temp
   real(8)       :: temptra                ! temperature, translational
   real(8)       :: temprot                ! temperature, rotational
   real(8)       :: prsrst                 ! pressure, start
!> \page prsr
!! `real`
!! * Desired pressure (only \ref txensemb='\ref npt' or 'nve' and with volume scaling).
   real(8)       :: prsr
   real(8)       :: prsrreds3              ! pressure, final reduced pressure
   real(8)       :: prsrredsd              ! pressure, sd of final reduced pressure
   real(8)       :: virial                 ! virial
   real(8)       :: volst                  ! volume, start
   real(8)       :: vol                    ! volume

   real(8)       :: tempaver               ! average temperature (after a macrostep and simulation)
   real(8)       :: prsraver               ! average pressure (after a macrostep and simulation)
   real(8)       :: volaver                ! average volume (after a macrostep and simulation)
   real(8)       :: npartaver              ! average numer of particles (after a macrostep and simulation)

! ... random number
!> \page iseed
!! `integer`
!! * Seed of random number generator. \ref iseed > 0 is required.
   integer(k4b)  :: iseed
!> \page maxcpu
!! `integer`
!! **default:** `0`
!! * `0`: Nothing.
!! * `<0`: The program stops before starting next macrostep if the total cpu time including the next one exceeds \ref maxcpu (in
!!         seconds). The cpu time for the last macrostep is used as the estimate of the cpu time for the next macrostep. This is handy if
!!         batch queue installation is used. In such cases, if the job exceeds the maximum cpu time for the batch job, the job might stop
!!         abruptly (depending on system installation) without executing the remaining commands in the flink command file and possibly
!!         corrupting the file FCNF. If so, the entirely job is lost. The variable \ref maxcpu allows the program to stop itself and promptly send
!!         the FCNF file which then can be used to continue the simulation by using the option \ref txstart='continue'.
   integer(4)    :: maxcpu                 ! maximum number of cpu time (seconds)

! ... configuration

   integer(4)    :: nstep                  ! nstep1*nstep2
!> \page nstep1
!! `integer`
!! * Number of macrosteps. The total number of steps/passes are \ref nstep1*\ref nstep2. After the simulation grand averages are calculated and written for the \ref nstep1*\ref nstep2 steps/passes.
   integer(4)    :: nstep1
!> \page nstep2
!! `integer`
!! * Number of steps (only MD or BD), or number of passes (only MC or MCALL) per macrostep. One pass is one attempt to move each
!!   particle in average (only MC) or one attempt to simultaneously move all particles (only MCALL). Averages are calculated and written
!!   for each set of \ref nstep2 steps/passes. After each set the FCNF file is updated with new coordinates and accumulated averages.
   integer(4)    :: nstep2
   integer(4)    :: istep                  ! present step/pass
   integer(4)    :: istep1                 ! present macrostep
   integer(4)    :: istep2                 ! present step/pass of one macrostep
   integer(4)    :: ipass                  ! present mc step of one pass
   integer(4)    :: nstep1beg              ! first macrostep to be performed
   integer(4)    :: nstep1done             ! number of macrosteps performed

! ... coordinate, etc
!> \page lintsite
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Positions of interaction sites are given by \ref raintin.
!! * `.false.`: Positions of interaction sites are given by \ref rain.
   logical              :: lintsite
   real(8), allocatable :: rasite(:,:,:)   ! interaction-site coordinates in molecular frame
   real(8), allocatable :: ra(:,:,:)       ! atom coordinates in molecular frame
   real(8), allocatable :: ro(:,:)         ! particle position
   real(8), allocatable :: rotm(:,:)       ! particle position, for trial configuration
   real(8), allocatable :: drotm(:,:)      ! suggested particle move
   real(8), allocatable :: ori(:,:,:)      ! particle orientation
   real(8), allocatable :: oritm(:,:,:)    ! particle orientation, for trial configuration
   real(8), allocatable :: rod(:,:)        ! particle linear velocity
   real(8), allocatable :: qua(:,:)        ! particle orientation, quaterion
   real(8), allocatable :: quad(:,:)       ! quaterion velocity
   real(8), allocatable :: angvelo(:,:)    ! angular velocities of particles, particle frame
   real(8), allocatable :: forceo(:,:)     ! forces on particles
   real(8), allocatable :: torqueo(:,:)    ! torques on particles
   real(8), allocatable :: drostep(:,:)    ! particle displacement in one step/pass

   real(8), allocatable :: r(:,:)          ! atom coordinates in lab frame
   real(8), allocatable :: rtm(:,:)        ! atom coordiantes, for trial configuration
   real(8), allocatable :: force(:,:)      ! forces on atoms
   real(8), allocatable :: torque(:,:)     ! torques on atoms

! ... neighbour and linked cell lists

   logical                 :: lvlist       ! flag for neighbour lists
   integer(4), allocatable :: ipnploc(:)   ! particle (local) -> particle (1:np)
   integer(4), allocatable :: nneighpn(:)  ! particle (local) -> number of neighbours
   integer(2), allocatable :: jpnlist(:,:) ! ineigh (local list) and ip (global or local) -> neigbour particle (1:np)
   logical                 :: lllist       ! flag for linked lists
   logical   , allocatable :: lcellllist(:)! list of cells to be used for a position (local)
   integer(4), allocatable :: headllist(:) ! head of the linked list
   integer(4), allocatable :: jpllist(:)   ! linked list
   logical                 :: lCellList       ! flag for cell lists

! ... two-body potential and its lookup table

   integer(4)              :: nbuf         ! length of buffer
   real(8),    allocatable :: ubuf(:)      ! buffer for potential table

!> \page rcut
!! `real`
!! **default:** `0.0`
!!
!! Cutoff  distance  of  the  potential  and  forces  based  on  the  particle - particle  center  of mass distance. If \ref rcut=0.0, the cutoff distance sets to
!! * `(boxlen2(1)**2 + boxlen2(2)**2 + boxlen2(3)**2)**1/2` ( only \ref txbc= ' xyz '),
!! * `(boxlen2(1)**2 + boxlen2(2)**2 + boxlen(3)**2)**1/2` ( only \ref txbc= ' xy '),
!! * `(boxlen(1)**2 + boxlen(2)**2 + boxlen2(3)**2)**1/2` ( only \ref txbc= ' z '),
!! * `2*rsph` ( only \ref txbc= ' sph '), or
!! * `(2*rcyl)**2+lcyl**2)**1/2` ( only \ref txbc= ' cyl'),
   real(8)                 :: rcut
   real(8)                 :: rcut2        ! rcut**2
   real(8),    allocatable :: r2umin(:)    ! lower limit squared of potential table
   integer(4), allocatable :: iubuflow(:)  ! points on the first entry for iatjat

! ... ewald summation
!> \page lewald
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Ewald summation. Implemented for \ref txpot='(1,6,12)', \ref txpot=nemo... with \ref txelec='pol', and default potential. Periodic boundary conditions are required.
!! * `.false.`: No Ewald summation.
   logical       :: lewald
!> \page lsurf
!! `logical`
!! **default:** `.true.`
!! * `.true.`: Inclusion of the surface term of the Ewald summation (corresponding to \f$ \epsilon\f$ (surrounding) = 1).
!! * `.false.`: Exclusion of the surface term of the Ewald summation (corresponding to \f$ \epsilon\f$ (surrounding) -> infinity).
   logical       :: lsurf
!> \page lewald2dlc
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Correction for making a 3d-periodic system 2d-periodic according to Arnold et al.  JCP 2002. Should only be applied
!!             with considerable insight on, e.g, the need of making the simulation box longer than the simulated system in the z-direction and
!!             the convergence.
!! * `.false.`: No such correction.
   logical       :: lewald2dlc
!> \page iewaldopt
!! `integer`
!! * Controls the choice of Ewald truncation analysis. \ref ualphared = \ref ualpha*\ref rcut.
!! \n
!! * For \ref txewaldrec = 'std':
!! * `0`: \ref ualphared, \ref rcut, and \ref ncut are used.
!! * `1`: \ref ualphared and \ref rcut are used to calculate \ref uewaldtol and \ref ncut.\f$^1\f$
!! * `2`: \ref uewaldtol and \ref ualpha are used to calculate \ref rcut and \ref ncut.\f$^1\f$
!! * `3`: \ref uewaldtol and \ref rcut are used to calculate \ref ualpha and \ref ncut.\f$^1\f$
!! * `4`: \ref uewaldtol and \ref ncut are used to calculate \ref ualpha and \ref rcut.\f$^1\f$
!! \n \f$^1\f$ According to Kolafa and Perram (charges) and Wang and Holm (dipoles).
!! \n
!! * For \ref txewaldrec = 'spm':
!! * `0`: \ref ualphared, \ref rcut, order and \ref nmesh are used.
!! * `3`: \ref uewaldtol, \ref rcut, order and \ref nmesh are use.
   integer(4)    :: iewaldopt
!> \page txewaldrec
!! `character`
!! **default:** '`std`'
!! * '`std`': standard Ewald summation.
!! * '`spm`': smooth particle mesh Ewald summation (require installation of the FFTW library from [www.fftw.org](http://www.fftw.org/); see
!!     makefile). If only the standard Ewald summation is used, calls to subroutines fftw… in files energy.F90 and denergy.F90 can be
!!     comment away.
   character(3)  :: txewaldrec
!> \page uewaldtol
!! `real`
!! **default:** `0.0`
!! * Potential energy tolerance in Ewald summation (see \ref iewaldopt).
   real(8)       :: uewaldtol
!> \page ualpha
!! `real`
!! * Ewald parameter (see \ref iewaldopt).
   real(8)       :: ualpha
!> \page ualphared
!! `real`
!! **default:** `3.0`
!! * Reduced Ewald parameter used for the Ewald summation (see \ref iewaldopt) and the Spherical Ewald potential.
   real(8)       :: ualphared
!> \page ncut
!! `integer`
!! * Largest number of k-vectors in one direction in the reciprocal space (see \ref iewaldopt).
   integer(4)    :: ncut
   integer(4)    :: ncut2                  ! square of largest k-vector
!> \page ncutregion
!! `character(6)`
!! **default:** '`sphere`'
!! * '`sphere`': Spherical k-space region.
!! * '`cube`': Cubic k-space region.
   character(6)  :: ncutregion
   integer(4)    :: ncut2d                 ! largerst 1d k-vector (2d correction)
   integer(4)    :: ncut2d2                ! square of largest k-vector (2d correction)
!> \page nmesh
!! `integer`
!! **default:** `48`
!! * Number of meshes used in the reciprocal space (see \ref iewaldopt ).
   integer(4)    :: nmesh                  ! number of meshes (SPME)  TODO: non-cubic
!> \page order
!! `integer`
!! **default:** `5`
!! * Interpolation order in the reciprocal space (see \ref iewaldopt ).
   integer(4)    :: order                  ! order of cardinal B-splines (SPME)

! ... reaction field method
!> \page lrf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Reaction field method (only \ref txelec = 'pol').
!! * `.false.`: No reaction field method.
   logical       :: lrf
!> \page epsrf
!! `real`
!! **default:** `78.0`
!! * Relative dielectric permittivity of the surrounding beyond \ref rcut for reaction field method. Zero means infinity.
   real(8)       :: epsrf
   real(8)       :: rffac                  ! 2*(epsstat-1)/((2*epsstat+1)*rcut**3)
   real(8)       :: rffac2                 ! rffac/2

! ... cutoff and shift of lj potential
!> \page lljcut
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Apply a cutoff and shift of the Lennard-Jones potential (only \ref txpot (ipt) = '(1,6,12)').
!! * `.false.`: No cut and shift.
   logical       :: lljcut
!> \page ljrcut
!! `real`
!! **default:** `2.0**(1.0/6.0)`
!! * Cutoff distance in LJ-sigma units.
   real(8)       :: ljrcut
!> \page ljushift
!! `real`
!! **default:** `1.0`
!! * Shift of LJ potential in LJ-epsilon units.
   real(8)       :: ljushift               ! shift in epsilon units

! ... screened coloumb interaction
!> \page lscrc
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Screened Coulomb potential. See also \ref lscrhs.
!! * `.false.`: No screened Coulomb potential.
   logical       :: lscrc
!> \page scrlen
!! `real`
!! * Screening length for the screened Coulomb potential.
   real(8)       :: scrlen
!> \page lscrhs
!! `logical`
!! **default:** `.false.`
!! * `.true.`: U(r) = \ref zat (iat) * jat)/r*exp(-r/\ref scrlen )*FAC(iat)*FAC(jat) with FAC(iat) = exp(-r*\ref radat (iat))/(1+\ref radat (iat)/\ref scrlen )
!! * `.false.`: U(r) = \ref zat (iat)*\ref zat (jat)/r*exp(-r/\ref scrlen ). If \ref zat all are zero \ref ucoff (1,iatjat) is used.
   logical       :: lscrhs

! ... image charge interaction
!> \page epsimg
!! `real`
!! * Relative dielectric permittivity of the surrounding medium (only \ref txelec = 'dipsphimage').
   real(8)       :: epsimg
!> \page radimg
!! `real`
!! * Radius of the surrounding medium (only \ref txelec = 'dipsphimage').
   real(8)       :: radimg
   real(8), allocatable :: zimg(:)         ! charge of image atom
   real(8), allocatable :: rimg(:,:)       ! position of image atom
   real(8), allocatable :: dipimg(:,:)     ! image dipole moment vector

! ... dielectric discontinuities (for box or sphere)
!> \page epsi1
!! `real`
!! * Dielectric constant of the sphere (only \ref txelec = 'dipdiel').
   real(8)       :: epsi1
!> \page epsi2
!! `real`
!! * Dielectric constant outside the sphere (only \ref txelec = 'dipdiel').
   real(8)       :: epsi2
!> \page boundaryrad
!! `real`
!! * Radius of dielectric boundary (only \ref txelec = 'dipdiel').
   real(8)       :: boundaryrad
!> \page lmaxdiel
!! `integer`
!! * Truncation of l-sum (only \ref txelec = 'dipdiel').
   integer(4)    :: lmaxdiel
!> \page lbg
!! `logical`
!! * `true.`: Volume charge density of inside the sphere neutralizing the system (only \ref txelec = 'dipdiel').
   logical       :: lbg                    ! flag for neutralizing background of atome type 1 inside the dielectric sphere
   real(8)       :: ubgfac                 ! for neutralizing background of atom type 1 inside the dielectric sphere

! ... gravitation potential

   real(8)       :: gravitation_force      ! gravitation force

! ... ramp potential
!> \page lambda_ramp
!! `real`
!! **default:** `1.1`
!! * End of the ramp potential in units of hard-sphere diameter.
   real(8)       :: lambda_ramp
!> \page epsilon_ramp
!! `real`
!! **default:** `1.0`
!! * Depth of the ramp potential.
   real(8)       :: epsilon_ramp
!> \page alpha_ramp
!! `real`
!! **default:** `1.0d-3`
!! * Regulate the soft slope change at \ref lambda_ramp. Higher alpha implies softer change. u(r) =
!! slope(r-r_ramp)(0.5+(1/Pi)*tan(\ref alpha_ramp(r-r_ramp), slope = -\ref epsilon_ramp /r1atat(\ref lambda_ramp -1), r_ramp = r1atat*\ref lambda_ramp for
!! r1atat < r1atat(\ref lambda_ramp -1).
   real(8)       :: alpha_ramp

! ... square_well potential
!> \page lambda_sw
!! `real`
!! **default:** `1.1`
!! * End of the square-well potential in units of hard-sphere diameter.
   real(8)       :: lambda_sw
!> \page epsilon_sw
!! `real`
!! **default:** `1.0`
!! * Depth of the square-well potential.
   real(8)       :: epsilon_sw
!> \page alpha_sw
!! `real`
!! **default:** `1.0d-3`
!! * Regulate the soft slope change at \ref lambda_sw . Higher alpha implies softer change.  u(r) =
!! slope(r-r_sw)(0.5+(1/Pi)*tan(\ref alpha_sw (r-r_sw), slope = -\ref epsilon_sw /r1atat(\ref lambda_sw -1), r_sw = r1atat*\ref lambda_sw for r1atat <
!! r1atat(\ref lambda_sw -1).
   real(8)       :: alpha_sw               ! parameter dermining the smoothness of the potential at lambda_sw

! ... depletion potentials (asakura_oosawa and lekkerkerker_tuinier)
!> \page rad_dep
!! `real`
!! * Radius of penetrable hard sphere (Asakura-Oosawa model)
   real(8)       :: rad_dep                ! radius of penetrable hard sphere (asakura-oosawa model)
!> \page rho_dep
!! `real`
!! * number density of penetrable hard sphere (asakura-oosawa model)
   real(8)       :: rho_dep
!> \page factor_dep
!! `real`
!! **default:** `1.0`
!! * depletion-thickness factor
   real(8)       :: factor_dep

! ... external field

   real(8)       :: ofaxis(1:3)            ! orienting field axis
   real(8)       :: ofstrength             ! orienting field strength

! ... chain and crosslink potential

   type(bond_var), allocatable :: bond(:)  ! atom-atom bond
   type(bond_var), allocatable :: angle(:) ! atom-atom-atom angle
   type(bond_var):: clink                  ! atom-atom crosslink

! ... external potential
!> \page luext
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Application of external potentials on the particles.
!! * `.false.`: No external potential.
   logical       :: luext
!> \page txuext
!! `character(20)`(1:\ref npt)
!! **default:** \ref npt*' '
!! * Text label used for selecting external potential. The number of external potential is growing and a number of parameters specifying these parameters are available. See routine  IOPotExternal for details.
!! * '`wall_z`': hard walls at abs(z) = wall_z-ext
!! * '`sw_wall_zlow`': square-well wall at z = -wall_z-ext
!! * '`ramp_wall_z`': ramp walls at abs(z) = wall_z-ext
!! * '`lj_wall_z`': LJ walls at abs(z) = wall_z-ext
!! * '`lj_wall_z_ts`': truncated and shifted LJ walls at abs(z) = wall_z-ext
!! * '`lj_wall_z_mod`': heterogeneous LJ walls at abs(z) = wall_z-ext
!! * '`lj_wall_zlow`': LJ + ramp wall at z = -wall_z-ext
!! * '`lj_wall_desorb`': truncated LJ + ramp wall at z = -wall_z-ext (Niklas)
!! * '`estat_field_z`': homogeneous electrical field in z-direction
!! * '`hom_charged_walls`': interaction with a charged surface
!! * '`i_soft_sphere`': external, soft, and spherical wall
!! * '`Gunnar_soft_sphere`': external, soft, and spherical wall (Gunnar)
!! * '`out_hard_ellipsoid`': hard ellipsoidal wall
!! * '`capsid_shell`'; hard spherical capsid shell
!! * '`uniform_shell`': hard spherical capsid shell with a uniform surface charge density
!! * '`sphdielboundary_q`': dielectric sphere (multipole expansion)
!! * '`sphdielboundary_p`': dielectric sphere (pairwise interaction)
!! * '`sphdielboundary`': spherical dielectric boundary (pairwise interaction)
!! * '`core_shell`': hard inner and outer spherical walls (Steffi)
!! * '`insulating_sphere`': penetrable uniformly charged sphere (Steffi)
!! * '`hollow_sphere' hollow and charged sphere (Steffi)
!! * See code for further details.
   character(20), allocatable :: txuext(:)
   real(8)       :: wall_z_ext             !*parameter describing location of external wall
   real(8)       :: lambda_sw_ext          !*parameter describing external square-well potential
   real(8)       :: epsilon_sw_ext         !*parameter describing external square-well potential
   real(8)       :: alpha_sw_ext           !*parameter describing external square-well potential
   real(8)       :: lambda_ramp_ext        !*parameter describing external ramp potential
   real(8)       :: epsilon_ramp_ext       !*parameter describing external ramp potential
   real(8), allocatable :: sigma_ext(:)    !*parameter describing external potential
   real(8), allocatable :: epsilon_ext(:)  !*parameter describing external potential
   real(8), allocatable :: z3coeff_ext(:)  !*parameter describing external potential
   real(8), allocatable :: z9coeff_ext(:)  !*parameter describing external potential
   real(8)       :: u0_drift               ! coefficient of drift potential udrift = u0 + u1*z
   real(8)       :: u1_drift               !*coefficient of drift potential udrift = u0 + u1*z
   real(8)       :: c_ext                  !*parameter describing external potential
   real(8)       :: lx_ext, TwoPilxi_ext   !*parameter describing external potential
   real(8)       :: ly_ext, TwoPilyi_ext   !*parameter describing external potential
   real(8), allocatable :: zmin_ext(:)     ! parameter describing external potential
   real(8), allocatable :: delta_ext(:)    ! parameter describing external potential
   real(8)       :: efield_ext(3)          !*parameter describing external potential
   integer(4), parameter :: mninchden = 202
   real(8)       :: surfchargeden          !*surface charge density of walls
   logical       :: llongrangecontr        !*flag for adding long-range (outside-boxlen) contrbution
   real(8)       :: zdist(mninchden)       !*array holding the z-coordinate
   real(8)       :: chden(mninchden)       !*array holding the charge distribution

   real(8)       :: rInsSphere             !*radius of the charged insulating sphere
   real(8)       :: zInsSphere             !*charge of the charged insulating sphere
   real(8)       :: rChargeIn              !*inner radius of the shell/core wall
   real(8)       :: rChargeOut             !*outer radius of the shell/core wall
   real(8)       :: rInSphere              !*inner radius of the shell
   real(8)       :: rOutSphere             !*outer radius of the shell
   real(8)       :: rCylinder              !*radius of the hard cylinder
   real(8)       :: zCylinder              !*charge density of the cylinder, elementary charges per unit length

   real(8), allocatable :: ruext(:,:)      !*parameter describing external potential
   real(8), allocatable :: ruexti(:,:)     ! parameter describing external potential
   real(8)       :: auext                  !*parameter describing external potential
   integer(4)    :: nuext                  !*parameter describing external potential
   real(8)       :: rcap                   !*capside (shell) inner radius
   real(8)       :: dcap                   !*capside (shell) thickness

! ... energy

   real(8)       :: ekin = 0.0d0           ! kinetic energy
   type(potenergy_var) :: u                ! potential energy
   type(potenergy_var) :: du               ! difference in potential energy between two configurations
   real(8)       :: htot                   ! total enthalpy

! ... static dipole moment

   real(8), allocatable :: dipa(:,:,:)     ! static dipole moment vector in molecular frame
   real(8), allocatable :: dip(:,:)        ! static dipole moment vector
   real(8), allocatable :: diptm(:,:)      ! static dipole moment vector, for trial move

! ... polarization

   logical, allocatable :: lapolarization(:)! true if polarization /= 0
!> \page tpolit
!! `real`
!! **default:** `10**-4`
!! * Relative tolerance of the induced dipole moment in the iteration.
   real(8)              :: tpolit
!> \page mpolit
!! `integer`
!! **default:** `15`
!! * Maximum number of iterations.
   integer(4)           :: mpolit
!> \page npolit
!! `integer`
!! **default:** `5`
!! * Interval of iteration.
   integer(4)           :: npolit
   logical              :: lidmconv        ! flag indicating convergence of idm iteration
   real(8), allocatable :: poltensa(:,:,:) ! polarisability tensor (xx,yy,zz,xy,xz,yz) in molecular frame
   real(8), allocatable :: poltens(:,:)    ! polarisability tensor (xx,yy,zz,xy,xz,yz)
   real(8), allocatable :: idm(:,:)        ! induced dipole moment
   real(8), allocatable :: idm1(:,:)       ! induced dipole moment, previous step
   real(8), allocatable :: idm2(:,:)       ! induced dipole moment, two steps back
   real(8), allocatable :: idmo(:,:)       ! induced dipole moment, particle
   real(8)              :: idmsys(3)       ! induced dipole moment, system
   real(8), allocatable :: diptot(:,:)     ! static + induced dipole moment vector
!> \page ldamping
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Short-distance damping of electrostatic force when calculating induced dipole moment.
!! * `.false.`: Nothing.
   logical              :: ldamping

! ... ellipsoidal particles
!> \page lellipsoid
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Hard-core ellipsoidal (prolate) particles.
!! * `.false.`: No hard-core ellipoidal particles.
   logical       :: lellipsoid
!> \page radellipsoid
!! `real`
!! **default:** `1.0`
!! * Radius of degenerated axes.
   real(8)       :: radellipsoid
!> \page aellipsoid
!! `real`
!! **default:** `1.0`
!! * Aspect ratio: >1 prolate and <1 oblate.
   real(8)       :: aellipsoid            ! aspect ratio (>1 prolate, <1 oblate)
   real(8)       :: radellipsoid2         ! radellipsoid**2

! ... superball particles

!> \page lsuperball
!! `logical`
!! **default:** `.false.`
!! * `.true.`: particles are superballs
   logical       :: lsuperball
!> \page radsuperball
!! `real`
!! **default:** `One`
!! radius of superballs
   real(8)       :: radsuperball
   real(8)       :: radsuperball2         ! radsuperball**2
   real(8)       :: radsuperball2i        ! 1/radsuperball2
!> \page qsuperball
!! `real`
!! **default:** `One`
!! * q parameter of superballs
   real(8)       :: qsuperball
   real(8)       :: qsuperball2           ! 2q
   real(8)       :: qsuperballi           ! 1/qsuperball
!> \page txmethodsuperball
!! `character(4)`
!! **default:** `'nr'`
!! * 'nr', 'mesh', 'test'
   character(4)  :: txmethodsuperball
!> \page nitersuperball
!! `integer`
!! **default:** `10`
!! * maximal number of nr iterations
   integer(4)    :: nitersuperball
!> \page tolsuperball
!! `real`
!! **default:** `1.0d-4`
!! * tolerance of nr iterations
   real(8)       :: tolsuperball
!> \page dl_damp
!! `real`
!! **default:** `One`

!> \page dl_cut
!! `real`
!! **default:** `1d10`
   real(8)       :: dl_damp, dl_cut
!> \page dr_damp
!! `real`
!! **default:** `One`

!> \page dr_cut
!! `real`
!! **default:** `1d10`
   real(8)       :: dr_damp, dr_cut
!> \page meshdepthsuperball
!! `integer`
!! **default:** `4`
!! * depth of mesh
   integer(4)    :: meshdepthsuperball
   logical       :: lsuperballtest        ! for test output
   real(8)       :: rcut2superball(2)     ! separations squared where refined overlap check is made
   type(TriMesh) :: superBallMesh         ! triangle mesh with DOP-tree
!> \page lstatsuperball
!! `logical`
!! **default:** `.false.`
!! * `.true`: engage time statistics
   logical       :: lstatsuperball
   real(8), parameter :: qsuperball_max_nr = 1.7 ! maximal q parameter of superballs for nr algorithm

! ... group division

   integer(4)                :: ngr(2)     ! number of groups of two different types
   integer(4)                :: ngrgr      ! number of group pairs (ngr(1)*ngr(2))
   integer(4)                :: maxngr     ! maximal number of groups (max of ngr(1:2))
   integer(4)                :: ngrvar     ! number of group variables (2+ngr(1)+ngr(2))
   type(scalar_var), allocatable, save :: grvar(:) ! containing group number averages
   integer(4),   allocatable :: igrpn(:,:) ! particle (1:np)               -> its group number (1:igr)
   integer(4),   allocatable :: iptgr(:,:) ! gruop number (1:ngr)          -> particle type (1: npt)
   integer(4),   allocatable :: iatgr(:,:) ! group number (1:ngr)          -> its first atom type (1:nat)
   integer(4),   allocatable :: natgr(:,:) ! number of atom types of particle in group igr
   integer(4),   allocatable :: igrgr(:,:) ! two groups (1:ngr)            -> group pair (1:ngrgr)
   integer(4),   allocatable :: igrpnt(:,:)! group number and "m"          -> group variale number
   character(4), allocatable :: txgr(:)    ! number label of group igr
   character(6), allocatable :: txgrgr(:)  ! number label of group pair igr, jgr

! ... scratch variables

   logical                   :: laux
   integer(4)                :: iaux
   real(8)                   :: raux
   real(8), allocatable      :: vaux(:,:)
   integer(4), allocatable   :: ivaux(:,:)

! ... parallel parameters (also used in the serial code)

   integer(4), parameter  :: mnproc = 1024 ! maximum number of processors

   integer(4)    :: nproc                  ! number of processors
   integer(4)    :: myid                   ! identity of processor, 0 to nproc-1 (processor specific)
   logical       :: master                 ! =.true. if master, otherwise .false. (processor specific)
   logical       :: slave                  ! =.true. if slave, otherwise .false. (processor specific)
   integer(4)    :: icmyid(2)              ! chain boundary (processor specific)
   integer(4)    :: ipmyid(2)              ! particle boundary (processor specific)
   integer(4)    :: iamyid(2)              ! atom boundary (processor specific)
   integer(4)    :: kvecmyid(2)            ! k-vector boundary (processor specific)
   integer(4)    :: kvecoffmyid            ! k-vector offset (processor specific)

end module MolModule
