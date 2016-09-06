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
!*                                                                      *
!*     MolModule                                                        *
!*                                                                      *
!************************************************************************

! ... main module for the Molsim software

module MolModule

   use StatisticsModule
   use MeshModule

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

   type potenergy_var
      real(8), allocatable :: twob(:)       ! two-body contribution (excluding ewald contribution)
      real(8), allocatable :: oneb(:)       ! one-body contribution (dielectric discontinuity)
      real(8)              :: tot           ! total
      real(8)              :: rec           ! reciprocal electrostatic space contribution (UEwald)
      real(8)              :: stat          ! static electrostatic contribution (umbodyp)
      real(8)              :: pol           ! polarization electrostatic contribution (umbodyp)
      real(8)              :: bond          ! bond contribution
      real(8)              :: angle         ! angle contribution
      real(8)              :: crosslink     ! bond (chains-nodes) contribution in hierarcical structures and networks
      real(8)              :: external      ! external contribution
   end type potenergy_var

! ... data structure for bonds etc

   type bond_var
      real(8)      :: k                     ! force constant
      integer(4)   :: p                     ! power
      real(8)      :: eq                    ! equilibrium value
   end type bond_var

! ... data structure for asorption conditions

   type adscond_var
      character(8) :: txobject              ! "plane_xy": lower and upper xy-plane of a parallelepiped
                                            ! "part=xxx": surface of a particle xxx consitutes the adsorbing surface
      character(3) :: txsurface             ! "-", "+", "-|+" for lower, upper, and both surfaces
      real(8)      :: dist                  ! threshold distance of adsorption for the com of a particle from plane or center of particle
   end type adscond_var

! ... data structure for chain properties

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
      real(8)    :: asph                    ! asphericity (JCP 100, 636(1994))
      real(8)    :: torp                    ! toroidicity
   end type chainprop_var

! ... data structure for simple averaging

   type aver_var
      real(8)      :: s2                    ! summation/averaging over steps
      real(8)      :: s1                    ! summation/averaging over macrosteps
   end type aver_var

! ... data structure for 1D static variables

   type static1D_var
      logical      :: l                     ! logical flag for engeagement
      real(8)      :: min                   ! minimum value
      real(8)      :: max                   ! maximum value
      integer(4)   :: nbin                  ! number of bins
      logical      :: lnorm                 ! logical flag for normalization
      character(10):: label                 ! title
      integer(4)   :: nvar                  ! expanded into nvar variables
   end type static1D_var

! ... data structure for 2D static variables

   type static2D_var
      logical      :: l                     ! logical flag for engeagement
      real(8)      :: min(2)                ! minimum value
      real(8)      :: max(2)                ! maximum value
      integer(4)   :: nbin(2)               ! number of bins
      logical      :: lnorm                 ! logical flag for normalization
      character(10):: label                 ! title
      integer(4)   :: nvar                  ! expanded into nvar variables
   end type static2D_var

! ... version, date and author

   character(29) :: txVersionDate = 'version 6.4.7, Sep 18, 2015 v1.2.0'
   character(9)  :: txAuthor      = 'Per Linse'

! ... external units

   character(10), parameter :: fin   = 'FIN  '
   character(10), parameter :: fout  = 'FOUT '
   character(10), parameter :: fcnf  = 'FCNF '
   character(10), parameter :: flist = 'FLIST'
   character(10), parameter :: fuser = 'FUSER'
   character(10), parameter :: fimg  = 'FIMG '
   character(10), parameter :: fvtf  = 'FVTF '
   character(10), parameter :: ftcl  = 'FTCL '
   integer(4),    parameter :: uin   = 1   ! input data
   integer(4),    parameter :: uout  = 2   ! output data
   integer(4),    parameter :: ucnf  = 3   ! configuration data
   integer(4),    parameter :: ulist = 8   ! list data
   integer(4),    parameter :: uuser = 9   ! user-provided data
   integer(4),    parameter :: uimg  = 10  ! image data
   integer(4),    parameter :: uvtf  = 11  ! VMD: vtf data
   integer(4),    parameter :: utcl  = 12  ! VMD: tcl script

! ... terms

   character(90) :: txtitle                ! user-provided title
   character(10) :: txmode                 ! 'simulation' generate configuration/trajectory by simulation
                                           ! 'analysis', read configuration/trajectory
                                           ! 'mixed', mixed tasks
   character(5)  :: txmethod               ! 'md', 'mc', 'mcall', 'bd'
   character(3)  :: txensemb               ! 'nve', 'nvt', ntp', 'mvt'
   character(8)  :: txstart                ! 'setconf', 'zero', 'continue', 'readfin'
   logical       :: lcont                  ! write control data
   logical       :: laver                  ! write averages
   logical       :: lti                    ! make thermodynamic integration
   logical       :: ldist                  ! calculate distribution functions
   logical       :: ldump                  ! dump data
   logical       :: lgroup                 ! make group division
   logical       :: lstatic                ! perform static analysis
   logical       :: ldynamic               ! perform dynamic analysis
   logical       :: limage                 ! write coordinates for images
   integer(4)    :: itest                  ! =1, write various data after each step/trial move
                                           ! =3, write energies, forces, and virials
                                           ! =4, write neighbour lists
                                           ! =5, write detailed superball overlap-check data
                                           ! =7, examine truncation error of the Ewald summation
   logical       :: ltrace                 ! write tracing data on trace.master.data and trace.slave.data
   logical       :: lblockaver             ! write blockavering data on blockaver.data
   character(80) :: txuser                 ! enable user-provided code in the Molsim core, code is not verified
   integer(4)    :: ipart                  ! if ipart > 0, write every ipart:th particle
   integer(4)    :: iatom                  ! if iatom > 0, write every iatom:th atom
   integer(4)    :: iaver                  ! if iaver > 0, write cumulative averages every iaver:th step/pass
   integer(4)    :: ishow                  ! if ishow > 0, write every ishow:th value of functions
   integer(4)    :: iplot                  ! if iplot > 0, plot functions
   integer(4)    :: ilist                  ! if ilist > 0, write every ilist:th value of functions on unit flist
   integer(4)    :: idump                  ! if idump > 0, dump every dump:th step/pass
   integer(4)    :: istatic                ! if istatic > 0, call static analysis routines every istatic step/pass:

   logical       :: lsim                   !=.true. simulation
   logical       :: lana                   !=.true. analysis
   logical       :: lmix                   !=.true. mixed tasks

   logical       :: ltime                  !=.true. timing statistics

   logical       :: lnve                   !=.true. if NVE ensemble
   logical       :: lnvt                   !=.true. if NVT ensemble
   logical       :: lntp                   !=.true. if NTP ensemble
   logical       :: lmvt                   !=.true. if MVT ensemble
   real(8)       :: facmvt = 4.0           ! factorial increase of memory allocation, applied on np and na (lmvt)

! ... cell

   character(3)  :: txbc                   ! control boundary conditions
   logical       :: lbcbox                 ! box-like cell (rätblock)
   logical       :: lbcrd                  ! rhombic dodecahedral cell
   logical       :: lbcto                  ! truncated octahedral cell
   logical       :: lbcsph                 ! spherical cell
   logical       :: lbccyl                 ! cylindrical cell
   logical       :: lbcell                 ! ellipsoidal cell
   real(8)       :: boxlen(3)              ! box lengths
   real(8)       :: boxlen2(3)             ! boxlen/2
   real(8)       :: boxleni(3)             ! 1/boxlen
   real(8)       :: boxlen2i(3)            ! 1/(2*boxlen)
   real(8)       :: boxlenshort            ! minval(boxlen(1:3))
   real(8)       :: TwoPiBoxi(3)           ! 2*pi/boxlen
   real(8)       :: cellside               ! length of cell side (RD or TO cell)
   real(8)       :: sphrad                 ! radius of spherical cell
   real(8)       :: sphrad2                ! sphrad**2
   real(8)       :: ellrad(3)              ! radius of ellipsoidal cell
   real(8)       :: ellradi(3)             ! 1/ellrad
   real(8)       :: cylrad                 ! radius of cylindrical cell
   real(8)       :: cylrad2                ! cylrad**2
   real(8)       :: cyllen                 ! length of cylindrical cell
   real(8)       :: cyllen2                ! length of cylindrical cell/2
   real(8)       :: lenscl                 ! scaling factor of box length and particle positions
   logical       :: lPBC                   ! periodic boundary conditions
   real(8)       :: dpbc(3)                ! =boxlen for some pbc, otherwise zero

! ... scaling

    real(8)      :: scllen                 ! length
    real(8)      :: sclmas                 ! mass
    real(8)      :: scltem                 ! temperature
    real(8)      :: scltim                 ! time
    real(8)      :: sclene                 ! energy
    real(8)      :: sclhca                 ! heat capacity
    real(8)      :: sclpre                 ! pressure
    real(8)      :: scldif                 ! diffusion coefficient
    real(8)      :: sclvol                 ! volume
    real(8)      :: sclvel                 ! linear velocity
    real(8)      :: sclacc                 ! linear acceleration
    real(8)      :: sclfor                 ! force
    real(8)      :: sclmin                 ! moment of intertia
    real(8)      :: sclang                 ! angle
    real(8)      :: beta                   ! 1/(kT)
    real(8)      :: Epsi0FourPi            ! ECh**2/(Four*pi*Epsi0*scllen)*AvNo/sclene
    real(8)      :: EpsiFourPi             ! Epsi0FourPi/relpermitt

! ... particle

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
   integer(4)              :: nct          !*number of chain types
   integer(4)              :: npt          !*number of particle types
   integer(4)              :: nat          ! number of atom types
   integer(4)              :: nc           ! number of chains
   integer(4)              :: np           ! number of particles
   integer(4)              :: na           ! number of atoms
   integer(4)              :: na_alloc     ! number of atoms (for memory allocation)
   integer(4)              :: np_alloc     ! number of particles (for memory allocation)
   integer(4)              :: ncct(mnct)   !*number of chains of a chain type
   integer(4)              :: nppt(mnpt)   !*number of particles of a particle type
   integer(4)              :: nptct(mnct)  ! number of particle types belonging to one chain type
   integer(4)              :: npct(mnct)   ! number of particles belonging to a chain type
   integer(4)              :: natpt(mnpt)  !*number of atom types belonging to a particle type
   integer(4), allocatable :: napt(:)      ! number of atoms belonging to a particle type
   integer(4), allocatable :: naat(:)      ! number of atoms of an atom type in one particle
   integer(4)    :: npptct(mnpt,mnct)      !*number of particles of a particle type of chain type
   integer(4)    :: naatpt(mnat,mnpt)      !*number of atoms of an atom type of a particle type
   integer(4)    :: nctct                  ! number of different chain type pairs
   integer(4)    :: nptpt                  ! number of different particle type pairs
   integer(4)    :: natat                  ! number of different atom type pairs
   character(30) :: txcopolymer(mnct)      !*type of copolymer
   logical       :: lspma                  !*control tacticity
   character(10) :: txct(mnct)             !*name of chain type
   character(10) :: txpt(mnpt)             !*name of particle type
   character(10) :: txat(mnat)             !*name of atom type
   character(20), allocatable :: txctct(:) ! name of chain type-chain type pair
   character(20), allocatable :: txptpt(:) ! name of particle type-particle type pair
   character(20), allocatable :: txatat(:) ! name of atom type-atom type pair
   logical                :: lmonoatom     !
   logical                :: lpolyatom     !
   real(8)                :: massat(mnat)  !*mass of atom type
   real(8), allocatable   :: masspt(:)     ! mass of particle type
   real(8), allocatable   :: massipt(:)    ! inverse mass of particle type
   real(8), allocatable   :: mompt(:,:)    ! moments of inertia of particle type
   real(8), allocatable   :: momipt(:,:)   ! inverse moments of inertia of particle type
   real(8), allocatable   :: massp(:)      ! mass of particle
   real(8), allocatable   :: massip(:)     ! inverse mass of particle
   real(8), allocatable   :: momp(:,:)     ! moments of inertia of particle
   real(8), allocatable   :: momip(:,:)    ! inverse moments of inertia of particle
   real(8)                :: radat(mnat)   !*hard-core radius of atom type
   integer(4), allocatable:: itradegfree(:)! number of translational degrees of freedom of a particle type
   integer(4), allocatable:: irotdegfree(:)! number of rotational degrees of freedom of a particle type
   real(8), allocatable   :: racom(:)      !
   real(8), allocatable   :: r1atat(:)     !
   real(8), allocatable   :: r2atat(:)     !
   real(8)                :: zat(mnat)     !*charge of atoms of a given type
   real(8)                :: zatalpha(mnat)!*1/(sqrt(2) times the width of Gaussian charge distribution)
   real(8)                :: sigat(mnat)   !*LJ sigma parpameter for atomes of a given type
   real(8)                :: epsat(mnat)   !*LJ epsilon parameter for atomes of a given type
   real(8), allocatable   :: az(:)         ! atom charge
   logical                :: latweakcharge(mnat) !*.true. if weak charge among atoms of a given type
   real(8)                :: pK(mnat)      !*pKa
   real(8)                :: pH            ! pH
   real(8)                :: pHmpK(mnat)   ! pH - pKa
   integer(4)             :: iatweakcharge ! type of atom carrying weak charge
   integer(4)             :: jatweakcharge !*type of atom carrying counter charge to weak charge (0 means no counter charge)
   integer(4), allocatable:: iananweakcharge(:) ! atom carrying weak charge -> its atom carrying its couterion charge
   logical, allocatable   :: laz(:)        ! .true. if atom is charged
   logical, allocatable   :: laztm(:)      ! .true. if atom is charged (trial move)
   logical, allocatable   :: lpset(:)      ! logical valiable being .true. for setted paricles
   logical                :: lradatbox     ! use atom and atom radius when checking if particle inside box

                                           ! equivalences:
                                           ! nc         = sum(ncct(1:nct))
                                           ! np         = sum(nppt(1:npt))
                                           ! na         = sum(napt(1:npt)*nppt(1:npt))
                                           ! nat        = sum(natpt(1:npt))
                                           ! nptct(ict) = count(npptct(1:npt,ict) > 0)
                                           ! npct(ict)  = sum(npptct(1:npt,ict))
                                           ! napt(ipt)  = sum(naatpt(1:natpt(ipt),ipt)))
                                           ! naat(iat)  = naatpt(ka,iptat(iat))

! ... chain bond

   logical       :: lchain                 ! flag for connecting particles into chains by bonds
   logical       :: lfixedori              ! fix the particle orientation relative to the bond vectors, mirrorplane is xy plane
   integer(4), allocatable :: bondnn(:,:)  ! bond and particle             -> bonded particle

! ... crosslinks between chains

   logical       :: lclink                 ! flag for enabeling crosslinks between particles of diff. chains
   integer(4)    :: maxnbondcl(mnpt)       !*maximum number of crosslink to/from a particle allowed
   integer(4)    :: ncl                    ! number of crosslinks
   integer(4), allocatable :: nbondcl(:)   ! actual number of crosslinks to/from particle ip
   integer(4)    :: maxvalnbondcl          ! maxval(nbondcl(:))
   integer(4), allocatable :: bondcl(:,:)  ! crosslink and particle        -> crosslinked particle

   logical       :: lmultigraft            ! flag for multigrafted chains

! ... hierarchical structure

   logical       :: lhierarchical          ! flag for hierarchical structures
   integer(4)    :: ngen                   !*number of generations
   integer(4)    :: nh                     ! number of hierarchical structures
   integer(4)    :: nch(0:mngen)           ! number of chains of a generation in a hierarchial structure
   integer(4)    :: nphn                   ! number of particles of a hierarchical structure
   integer(4)    :: ictgen(0:mngen)        !*generation number             -> chain type
   integer(4), allocatable :: genic(:)     ! chain number                  -> generation number
   integer(4)    :: nbranch(0:mngen-1)     !*number of branches
   integer(4)    :: ibranchpbeg(0:mngen-1) !*particle number of the chain for first branch point
   integer(4)    :: ibranchpinc(0:mngen-1) !*particle increment of the chain between branch points

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

! ... pointers

   integer(4), allocatable :: ictcn(:)     ! chain (1:nc)                  -> its chain type (1:nct)
   integer(4), allocatable :: ictpt(:)     ! particle type (1:npt)         -> its chain type (1:nct)
   integer(4), allocatable :: ihnpn(:)     ! particle (1:np)               -> its hierarchical structure (1:nh)
   integer(4), allocatable :: ictpn(:)     ! particle (1:np)               -> its chain type (1:nct)
   integer(4), allocatable :: icnpn(:)     ! particle (1:np)               -> its chain (1:nc)
   integer(4), allocatable :: iptpn(:)     ! particle (1:np)               -> its particle type (1:npt)
   integer(4), allocatable :: iptat(:)     ! atom type (1:nat)             -> its particle type (1:npt)
   integer(4), allocatable :: iptan(:)     ! atom (1:na)                   -> its particle type (1:npt)
   integer(4), allocatable :: ipnan(:)     ! atom (1:na)                   -> its particle (1:np)
   integer(4), allocatable :: iatan(:)     ! atom (1:na)                   -> its atom type (1:na)
   integer(4)              :: ipnhn        ! hierarchical structure        -> its first particle
   integer(4), allocatable :: icnct(:)     ! chain type (1:nct)            -> its first chain (1:nc)
   integer(4), allocatable :: ipnpt(:)     ! particle type (1:npt)         -> its first particle (1:np)
   integer(4), allocatable :: iatpt(:)     ! particle type (1:npt)         -> its first atom type (1:nat)
   integer(4), allocatable :: ianpn(:)     ! particle (1:np)               -> its first atom (1:na)
   integer(4), allocatable :: ianat(:)     ! atom type (1:nat)             -> its first atom (1:na)
   integer(4), allocatable :: ictct(:,:)   ! two chain types (1:nct)       -> chain type pair (1:nctct)
   integer(4), allocatable :: iptpt(:,:)   ! two particle types (1:npt)    -> particle type pair (1:nptpt)
   integer(4), allocatable :: iatat(:,:)   ! two atom types (1:nat)        -> atom type pair (1:natat)
   integer(4), allocatable :: isegpn(:)    ! particle (1:np)               -> segment number
   integer(4), allocatable :: ipnsegcn(:,:)! seg (1:npct) and chain (1:nc) -> particle (1:np)
   integer(4), allocatable :: icihigen(:,:)! hier (1:nh) and gen (0:ngen)  -> its first chain

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
   logical       :: lmcweight              ! flag for using weighting function

! ... bd

   logical       :: lbd                    ! flag for Brownian dynamics simulation
   real(8), allocatable :: dcoeff(:)       !*self-diffusion coefficient

! ... averages

   real(8)       :: tempst                 ! temperature, start
   real(8)       :: temp                   ! temperature
   real(8)       :: temptra                ! temperature, translational
   real(8)       :: temprot                ! temperature, rotational
   real(8)       :: prsrst                 ! pressure, start
   real(8)       :: prsr                   ! pressure
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

   integer(4)    :: iseed                  ! seed of random number generator
   integer(4)    :: maxcpu                 ! maximum number of cpu time (seconds)

! ... configuration

   integer(4)    :: nstep                  ! nstep1*nstep2
   integer(4)    :: nstep1                 ! number of macrosteps
   integer(4)    :: nstep2                 ! number of steps/passes per macrostep
   integer(4)    :: istep                  ! present step/pass
   integer(4)    :: istep1                 ! present macrostep
   integer(4)    :: istep2                 ! present step/pass of one macrostep
   integer(4)    :: ipass                  ! present mc step of one pass
   integer(4)    :: nstep1beg              ! first macrostep to be performed
   integer(4)    :: nstep1done             ! number of macrosteps performed

! ... coordinate, etc

   logical              :: lintsite        ! true if special interaction sites are used instead of atomic sites
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

! ... two-body potential and its lookup table

   integer(4)              :: nbuf         ! length of buffer
   real(8),    allocatable :: ubuf(:)      ! buffer for potential table
   real(8)                 :: rcut         ! interaction cutoff
   real(8)                 :: rcut2        ! rcut**2
   real(8),    allocatable :: r2umin(:)    ! lower limit squared of potential table
   integer(4), allocatable :: iubuflow(:)  ! points on the first entry for iatjat

! ... ewald summation

   logical       :: lewald                 ! flag for ewald summation
   logical       :: lsurf                  ! flag for surface term
   logical       :: lewald2dlc             ! flag for 2d layer correction for effectively 2d periodic systems
                                           ! see arnold et al. JCP 117, 2496 (2002), 117, 2503 (2002)
   integer(4)    :: iewaldopt              ! selection of method to determine alpha, rcut, and ncut
   character(3)  :: txewaldrec             ! select 'std' (standard ewald) or 'spm' (smooth particle mesh)
   real(8)       :: uewaldtol              ! potential energy tolerance of ewald summation (kJ/mol particles)
   real(8)       :: ualpha                 ! damping factor
   real(8)       :: ualphared              ! damping factor / rcut
   integer(4)    :: ncut                   ! largerst 1d k-vector
   integer(4)    :: ncut2                  ! square of largest k-vector
   character(6)  :: ncutregion             ! 'cube' or 'sphere'
   integer(4)    :: ncut2d                 ! largerst 1d k-vector (2d correction)
   integer(4)    :: ncut2d2                ! square of largest k-vector (2d correction)
   integer(4)    :: nmesh                  ! number of meshes (SPME)  TODO: non-cubic
   integer(4)    :: order                  ! order of cardinal B-splines (SPME)

! ... reaction field method

   logical       :: lrf                    ! flag for reaction field method
   real(8)       :: epsrf                  ! static dielectric constant of surroundings
   real(8)       :: rffac                  ! 2*(epsstat-1)/((2*epsstat+1)*rcut**3)
   real(8)       :: rffac2                 ! rffac/2

! ... cutoff and shift of lj potential

   logical       :: lljcut                 ! flag for cutoff and shift of lj potential
   real(8)       :: ljrcut                 ! cutoff distance in sigma units
   real(8)       :: ljushift               ! shift in epsilon units

! ... screened coloumb interaction

   logical       :: lscrc                  ! flag for screened coloumb interaction
   real(8)       :: scrlen                 ! screening length
   logical       :: lscrhs                 ! flag for engaging incomplete screening due to excluded volume

! ... image charge interaction

   real(8)       :: epsimg                 ! dielectric constant of surrounding
   real(8)       :: radimg                 ! location of dielectric discontinuity
   real(8), allocatable :: zimg(:)         ! charge of image atom
   real(8), allocatable :: rimg(:,:)       ! position of image atom
   real(8), allocatable :: dipimg(:,:)     ! image dipole moment vector

! ... dielectric discontinuities (for box or sphere)

   real(8)       :: epsi1                  ! dielectric permittivity of left/inner region
   real(8)       :: epsi2                  ! dielectric permittivity of right/outer region
   real(8)       :: boundaryrad            ! radius of dielectric discontinuity (spherical bc)
   integer(4)    :: lmaxdiel               ! maximal number of terms (spherical bc)
   logical       :: lbg                    ! flag for neutralizing background of atome type 1 inside the dielectric sphere
   real(8)       :: ubgfac                 ! for neutralizing background of atom type 1 inside the dielectric sphere

! ... gravitation potential

   real(8)       :: gravitation_force      ! gravitation force

! ... ramp potential

   real(8)       :: lambda_ramp            ! r1atat*(lambda_ramp-One) distance where ramp potential becomes zero
   real(8)       :: epsilon_ramp           ! energy at the hards-sphere contact
   real(8)       :: alpha_ramp             ! parameter dermining the smoothness of the potential at lambda_ramp

! ... square_well potential

   real(8)       :: lambda_sw              ! r1atat*(lambda_sw-One) distance where sw potential becomes zero
   real(8)       :: epsilon_sw             ! energy at the hard-core contact
   real(8)       :: alpha_sw               ! parameter dermining the smoothness of the potential at lambda_sw

! ... depletion potentials (asakura_oosawa and lekkerkerker_tuinier)

   real(8)       :: rad_dep                ! radius of penetrable hard sphere (asakura-oosawa model)
   real(8)       :: rho_dep                ! number density of penetrable hard sphere (asakura-oosawa model)
   real(8)       :: factor_dep             ! depletion-thickess factor

! ... external field

   real(8)       :: ofaxis(1:3)            ! orienting field axis
   real(8)       :: ofstrength             ! orienting field strength

! ... chain and crosslink potential

   type(bond_var), allocatable :: bond(:)  ! atom-atom bond
   type(bond_var), allocatable :: angle(:) ! atom-atom-atom angle
   type(bond_var):: clink                  ! atom-atom crosslink

! ... external potential

   logical       :: luext                  ! flag for external potential
   character(20), allocatable :: txuext(:) !*selecting type of external potential
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

   real(8)       :: ekin                   ! kinetic energy
   type(potenergy_var) :: u                ! potential energy
   type(potenergy_var) :: du               ! difference in potential energy between two configurations
   real(8)       :: htot                   ! total enthalpy

! ... static dipole moment

   real(8), allocatable :: dipa(:,:,:)     ! static dipole moment vector in molecular frame
   real(8), allocatable :: dip(:,:)        ! static dipole moment vector
   real(8), allocatable :: diptm(:,:)      ! static dipole moment vector, for trial move

! ... polarization

   logical, allocatable :: lapolarization(:)! true if polarization /= 0
   real(8)              :: tpolit          ! relative tolerance of the idm for convergence
   integer(4)           :: mpolit          ! max number of iterations
   integer(4)           :: npolit          ! interval of interation
   logical              :: lidmconv        ! flag indicating convergence of idm iteration
   real(8), allocatable :: poltensa(:,:,:) ! polarisability tensor (xx,yy,zz,xy,xz,yz) in molecular frame
   real(8), allocatable :: poltens(:,:)    ! polarisability tensor (xx,yy,zz,xy,xz,yz)
   real(8), allocatable :: idm(:,:)        ! induced dipole moment
   real(8), allocatable :: idm1(:,:)       ! induced dipole moment, previous step
   real(8), allocatable :: idm2(:,:)       ! induced dipole moment, two steps back
   real(8), allocatable :: idmo(:,:)       ! induced dipole moment, particle
   real(8)              :: idmsys(3)       ! induced dipole moment, system
   real(8), allocatable :: diptot(:,:)     ! static + induced dipole moment vector
   logical              :: ldamping        ! damping of electrostatics when evaluating the polarization

! ... ellipsoidal particles

   logical       :: lellipsoid            ! .true.; particles are to be treated as ellipsoids !restrictions apply!
   real(8)       :: radellipsoid          ! radius of degenerated axes
   real(8)       :: aellipsoid            ! aspect ratio (>1 prolate, <1 oblate)
   real(8)       :: radellipsoid2         ! radellipsoid**2

! ... superball particles

   logical       :: lsuperball            ! .true.; particles are superballs
   real(8)       :: radsuperball          ! radius of superballs
   real(8)       :: radsuperball2         ! radsuperball**2
   real(8)       :: radsuperball2i        ! 1/radsuperball2
   real(8)       :: qsuperball            ! q parameter of superballs
   real(8)       :: qsuperball2           ! 2q
   real(8)       :: qsuperballi           ! 1/qsuperball
   character(4)  :: txmethodsuperball     ! 'nr', 'mesh', 'test'
   integer(4)    :: nitersuperball        ! maximal number of nr iterations
   real(8)       :: tolsuperball          ! tolerance of nr iterations
   real(8)       :: dl_damp, dl_cut
   real(8)       :: dr_damp, dr_cut
   integer(4)    :: meshdepthsuperball    ! deepth of mesh (4 is normally OK)
   logical       :: lsuperballtest        ! for test output
   real(8)       :: rcut2superball(2)     ! separations squared where refined overlap check is made
   type(TriMesh) :: superBallMesh         ! triangle mesh with DOP-tree
   logical       :: lstatsuperball        ! .true.; engage time statistics
   real(8), parameter :: qsuperball_max_nr = 1.7 ! maximal q parameter of superballs for nr algorithm

! ... group division

   integer(4)                :: ngr(2)     ! number of groups of two different types
   integer(4)                :: ngrgr      ! number of group pairs (ngr(1)*ngr(2))
   integer(4)                :: maxngr     ! maximal number of groups (max of ngr(1:2))
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
