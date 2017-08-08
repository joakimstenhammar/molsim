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
!*     MCModule                                                         *
!*                                                                      *
!************************************************************************

! ... module for mc

module MCModule

   use MolModule

   real(8), parameter :: ddelta = 1.0d-10  ! to ensure correct roundoff

! ... to be used in later developments

   type cluster2_tm_var
      logical      :: l                    ! flag for cluster1 tm
      real(8)      :: p                    ! relative probability of cluster2 move
      real(8)      :: dtran                ! maximal translational trial displacement
      real(8)      :: drot                 ! maximal rotational trial displacement
      integer(4)   :: mode                 ! =0 : search members only of type iptmove
                                           ! =1 : search members across all particle types
   end type cluster2_tm_var

   type cluster1_tm_var
      logical      :: l                    ! flag for cluster1 tm
      real(8)      :: p                    ! relative probability of cluster1 move
      real(8)      :: rad                  ! radius of region for cluster members
      real(8)      :: psel                 ! probability to select a particle within rad
      real(8)      :: dtran                ! maximal translational trial displacement
      real(8)      :: drot                 ! maximal rotational trial displacement
   end type cluster1_tm_var

   type trialmove_var
      logical      :: l                    ! flag for type of trial move
      real(8)      :: p                    ! relative probability of type of trial move
      real(8)      :: dtran                ! maximal translational trial displacement
      real(8)      :: drot                 ! maximal rotational trial displacement
      integer(4)   :: mode                 ! specific for type of trial displacment
      logical      :: lcl1                 ! flag for cluster1 trial move
      real(8)      :: radcl1               ! radius of region for cluster members
      real(8)      :: pselcl1              ! probability to select a particle within rad
   end type trialmove_var

! ... mc trial move variables

   logical, allocatable       :: lptmove(:)          ! flag of moving particle of given type
   integer(4)                 :: ipmove              ! particle to be moved
   integer(4)                 :: iptmove             ! type of particle to be moved
   integer(4)                 :: isamp               ! = 0 sequential selection of particles
                                                     ! = 1 random selection of particles

   logical                    :: lpspart             ! flag of single-particle move
   real(8), allocatable       :: pspart(:)           ! probability of single-particle move
   real(8), allocatable       :: dtran(:)            ! translation parameter of single-particle move
   real(8), allocatable       :: drot(:)             ! rotational parameter of single-particle move
   logical, allocatable       :: lcl1spart(:)        ! cluster1 move
   logical, allocatable       :: lfixzcoord(:)       ! no displacement in z-direction (applicable to spart and spartcl2 trial moves)
   logical, allocatable       :: lfixxycoord(:)      ! no displacement in xy-plane (applicable to spart and spartcl2 trial moves)
   logical, allocatable       :: lshiftzcom(:)       ! shift to get com at z = 0 if possible
   logical                    :: lfixchainstartspart ! if first particle of a chain, no move
   integer(4)                 :: ispart              ! =0 normal, >0 for development/special use
   logical                    :: lpspartcl2          ! flag of single-particle + cluster2 move
   real(8), allocatable       :: pspartcl2(:)        ! probability of single-particle + cluster2 move
   character(3), allocatable  :: txmembcl2(:)        ! = 'ipt' seach members only of same particle type as that displaced
                                                     ! = 'all' seach members among all partciles
   real(8), allocatable       :: radcl2(:)           ! max separation between clusters of type 1 for belonging to same cluster of type 2
   real(8), allocatable       :: dtrancl2(:)         ! translation parameter of single-particle + cluster2 move

   logical                    :: lppivot             ! flag of pivot rotation move
   real(8), allocatable       :: ppivot(:)           ! probability of pivot rotation move
   character(5), allocatable  :: txpivot(:)          ! = 'short' rotation of the shorter subchain
                                                     ! = 'lower' rotation of the lower subchain
                                                     ! = 'upper' rotation of the upper subchain
   character(6)               :: txdualpivot         ! = 'nodual' individual pivot  move
                                                     ! = 'linear' linear chain ict = 3 as reference (dual pivot = linear chain + branched diblock copolymer)
                                                     ! = 'combed' combed chain ict = 1 as reference (dual pivot = branched diblock +  copolymerlinear chain)
                                                     ! = 'mullin' idem as 'combed' and additionally works for icnct(ict) > 1
   real(8), allocatable       :: drotpivot(:)        ! rotation parameter of pivot rotation move
   real(8), allocatable       :: drotminpivot(:)     ! rotation parameter of pivot rotation move
   integer(4)                 :: ipivotrotmode       ! = 1 rotation around bond vector
                                                     ! = 2 rotation around a normal to a plane
                                                     ! = 3 rotation around a random direction
   logical, allocatable       :: lcl1pivot(:)        ! cluster1 move

   logical                    :: lpchain             ! flag of chain move
   real(8), allocatable       :: pchain(:)           ! probability of chain move
   real(8), allocatable       :: dtranchain(:)       ! translation parameter of chain move
   real(8), allocatable       :: drotchain(:)        ! rotational parameter of chain move
   logical, allocatable       :: lcl1chain(:)        ! cluster1 move

   logical                    :: lpslither           ! flag of chain slithering  move
   real(8), allocatable       :: pslither(:)         ! probability of chain slithering move

   logical                    :: lpbrush             ! flag of brush move
   real(8), allocatable       :: pbrush(:)           ! probability of brush move
   real(8), allocatable       :: dtranbrush(:)       ! translation parameter of brush move
   real(8), allocatable       :: drotbrush(:)        ! rotational parameter of brush move
   logical, allocatable       :: lcl1brush(:)        ! cluster1 move
   logical                    :: lpbrushcl2          ! flag of brush + cluster2 move
   real(8), allocatable       :: pbrushcl2(:)        ! probability of brush + cluster2 move
   real(8), allocatable       :: dtranbrushcl2(:)    ! translation parameter of brush + cluster2 move
   real(8), allocatable       :: drotbrushcl2(:)     ! rotational parameter of brush + cluster2 move

   logical                    :: lphierarchical      ! flag of hierarchical move
   real(8), allocatable       :: phierarchical(:)    ! probability of hierarchical move
   real(8), allocatable       :: dtranhierarchical(:)! translation parameter of hierarchical move

   logical                    :: lpnetwork           ! flag of network move
   real(8), allocatable       :: pnetwork(:)         ! probability of network move
   real(8), allocatable       :: dtrannetwork(:)     ! translation parameter of network move

   logical                    :: lpvol               ! flag of chain volume move
   real(8), allocatable       :: pvol(:)             ! probability of volume move
   real(8)                    :: dvol                ! volume change parameter

   logical                    :: lpnpart             ! flag of number of particle move
   real(8), allocatable       :: pnpart(:)           ! probability of number of particle move
   real(8), allocatable       :: chempot(:)          ! chemical potential

   real(8), allocatable       :: radcl1(:)           ! radius of cluster of type 1
   real(8), allocatable       :: pselectcl1(:)       ! probability of selecting a particle to belong to a cluster of type 1
   integer(4)                 :: npclnew             ! number of cluster particles of new conf
   integer(4)                 :: npclold             ! number of cluster particles of old conf
   integer(4)                 :: nprimpartcl         ! number of primary cluster particles of new & old conf

   logical                    :: lpcharge            ! flag of charge-change move
   real(8), allocatable       :: pcharge(:)          ! probability of charge-change move

   logical                    :: lpspartsso          ! flag for single particle move sso  ! Pascal Hebbeker
   logical, allocatable       :: lssopt(:)           ! flag for single particle move sso of particle types  ! Pascal Hebbeker
   logical                    :: lmcsep              ! flag for separating local from non-local moves
   real(8), allocatable :: pspartsso(:)              ! probability of single particle move sso
   real(8), allocatable :: plocal(:)
   real(8), allocatable       :: curdtranpt(:)            ! translation parameter of single-particle move



! ... mcall trial move variables

   real(8), allocatable        :: dtranall(:)
   real(8), allocatable        :: drotall(:)

! ... for mc statistics

   integer(4)    :: imovetype                        ! type of move
   integer(4), parameter :: nmovetype           = 13 ! number of different moves
   integer(4), parameter :: ispartmove          = 1  ! single-particle move
   integer(4), parameter :: ispartcl2move       = 2  ! single-particle + cluster2 move
   integer(4), parameter :: ipivotmove          = 3  ! pivot move
   integer(4), parameter :: ichainmove          = 4  ! chain move
   integer(4), parameter :: islithermove        = 5  ! slithering move
   integer(4), parameter :: ibrushmove          = 6  ! brush move
   integer(4), parameter :: ibrushcl2move       = 7  ! brush + cluster2 move
   integer(4), parameter :: ihierarchicalmove   = 8  ! hierarchical move
   integer(4), parameter :: inetworkmove        = 9  ! network move
   integer(4), parameter :: ivolumemove         = 10 ! volume change move
   integer(4), parameter :: inpartmove          = 11 ! number of particle change move
   integer(4), parameter :: ichargemove         = 12 ! charge change move
   integer(4), parameter :: ispartsso           = 13 ! single particle + SSO move   ! Pascal Hebbeker
   character(30), parameter :: txmovetype(1:nmovetype) = &
                                               [ 'single-particle move          ', &
                                                 'single-particle + cl2 move    ', &
                                                 'pivot move                    ', &
                                                 'chain move                    ', &
                                                 'slithering move               ', &
                                                 'brush move                    ', &
                                                 'brush + cl2                   ', &
                                                 'hierarchical move             ', &
                                                 'network move                  ', &
                                                 'volume change move            ', &
                                                 'number of particle change move', &
                                                 'charge-change move            ', &
                                                 'single-particle move + SSO    ' ]   ! Pascal Hebbeker

   integer(4)    :: ievent                   ! outcome of Metropolis test
   integer(4), parameter :: nevent       = 5 ! number of different Metropolis test outcome
   integer(4), parameter :: imcaccept    = 1 ! mc trial move accepted
   integer(4), parameter :: imcreject    = 2 ! mc trial move energy rejected
   integer(4), parameter :: imcboxreject = 3 ! mc trial move outside-box rejected
   integer(4), parameter :: imchsreject  = 4 ! mc trial move hard-core rejected
   integer(4), parameter :: imchepreject = 5 ! mc trial move hard-external-potential rejected
   character(35), parameter :: txevent(0:nevent) = &
                                               [ 'total number of conf.          = ', &
                                                 'accepted number of conf.       = ', &
                                                 'energy rejected conf.          = ', &
                                                 'outside-box rejected conf.     = ', &
                                                 'hard-core rejected conf.       = ', &
                                                 'hard-ext.-pot. rejected conf.  = ' ]

! ... unbrella sampling using a weighting function depending on the separation between two particles

   integer(4), parameter :: mnpolmcw = 10  ! maximum number of parameters of the weighting function
   integer(4)    :: ipmcw1                 ! particle number 1
   integer(4)    :: ipmcw2                 ! particle number 2
   character(11) :: txpotmcw               ! 'polynomial': acoeffmcw(0) + acoeffmcw(1)*r + ...
                                           ! 'exponential': acoeffmcw(0)*exp(-acoeffmcw(1)*(r-acoeffmcw(2)))
   integer(4)    :: npolmcw                ! degree of the polynomial
   real(8)       :: acoeffmcw(0:mnpolmcw)  ! coefficients of the weighting function

! ... automatic update of umbrella potential

   integer(4), parameter :: mumbgrid = 500 ! maximum grid for umbrella potential
   logical      :: lautumb                 ! true if automatic update of umbrella potential is used
   character(6) :: typeumb                 ! type of umbrella potential (e.g. particle-particle, particle-wall ...)
   real(8)      :: xumb(mumbgrid)          ! the potential of mean force
   real(8)      :: xumbmax, xumbmin        ! maximum(minimum) element of xumb
   integer(4)   :: iumb, iumbnew           ! index for bin in xumb
   integer(4)   :: ipumb1, ipumb2          ! for particle-particle or atom-atom ump: particle identifiers
   integer(4)   :: iaumb1, iaumb2          ! for atom-atom ump: atom identifiers
   real(8)      :: rminumbrella            ! for particle-particle or atom-atom ump: minimum particle-particle distance
   real(8)      :: delumb                  ! distance between two grid points in xumb
   integer(4)   :: numbgrid                ! number of grid points in xumb
   character(4) :: cupdate                 ! type of update of the weighting function
   logical      :: lradumb                 ! true if the umbrella potential is scaled to give g(r)
   character(1) :: umbcoord                ! if set, particles ipumb1 and ipumb2 can only move along the coordinate
                                           ! set by umbcoord with fixed orientation
   logical      :: lreadumb                ! true if the initial umbrella potential is read from file

! ... calculation of potential of mean force by updating weights

   integer(4), parameter :: mnbinmcpmf = 500 ! maximum number of bins of mc potential of mean force
   logical      :: lmcpmf                    ! =.true. calculation calculation of potential of mean force by updating weights
   integer(4)   :: iptmcpmf                  ! mc pmf between two particles of type iptmcpmf
   integer(4)   :: nbinmcpmf                 ! number of grid points of mcpmf
   real(8)      :: rlowmcpmf                 ! lower limit of mcpmf
   real(8)      :: ruppmcpmf                 ! upper limit of mcpmf
   real(8)      :: termmcpmf                 ! control the update of the mcpmf
   real(8)      :: binmcpmf                  ! bin of mc pmf
   real(8)      :: binimcpmf                 ! 1/binmcpmf
   real(8)      :: mcpmf(mnbinmcpmf)         ! potential of mean force
   real(8)      :: mcpmfmin                  ! maximum value of of mcpmf
   real(8)      :: mcpmfmax                  ! maximum value of of mcpmf
   integer(4)   :: ibinmcpmf, ibinnewmcpmf   ! index for bin in mcpmf

   integer(4)   :: idum

! ... test output

   integer(4)   :: itestmc                   ! = 1, call of TestIOMCProb
                                             ! = 2, call of TestMCMove
                                             ! = 3, call of TestVolChange from VolChange
                                             ! = 3, call of TestNPartChange1 from VolNParChange
                                             ! = 3, call of TestNPartChange2 from VolChange
                                             ! = 3, call of TestChargeChange1 from ChargeChange
                                             ! = 3, call of TestChargeChange2 from ChargeChange

   interface
      subroutine SPartMove(iStage, loptsso)
         integer(4), intent(in) :: iStage
         logical, optional, intent(in) :: loptsso
      end subroutine SPartMove
   end interface


   contains
      subroutine CallMove(prandom, iptmove, iStage)

         implicit none
         real(8), intent(in) :: prandom
         integer(4), intent(in) :: iptmove
         integer(4), intent(in) :: iStage

         !first local moves
         if (prandom < pspart(iptmove)) then
            call SPartMove(iStage)                             ! single-particle trial move
         else if (prandom < pspartsso(iptmove)) then
            call SPartMove(iStage, loptsso=.true.)

         !then non-local moves
         else if (prandom < pspartcl2(iptmove)) then
            call SPartCl2Move(iStage)                          ! single-particle + cluster2 trial move
         else if (prandom < ppivot(iptmove)) then
            call PivotMove(iStage)                             ! pivot rotation trial move
         else if (prandom < pchain(iptmove)) then
            call ChainMove(iStage)                             ! chain trial move
         else if (prandom < pslither(iptmove)) then
            call SlitherMove(iStage)                           ! chain slithering trial move
         else if (prandom < pbrush(iptmove)) then
            call BrushMove(iStage)                             ! brush trial move
         else if (prandom < pbrushcl2(iptmove)) then
            call BrushCl2Move(iStage)                          ! brush + cluster2 trial move
         else if (prandom < phierarchical(iptmove)) then
            call HierarchicalMove(iStage)                       ! hierarchical trial move
         else if (prandom < pnetwork(iptmove)) then
            call NetworkMove(iStage)                           ! network trial move
         else if (prandom < pvol(iptmove)) then
            call VolChange(iStage)                             ! volume change trial move
         else if (prandom < pnpart(iptmove)) then
            call NPartChange(iStage)                           ! number of particle change trial move
         !else if (prandom < pcharge(iptmove)) then
         else
            call ChargeChange(iStage)                          ! charge change trial move
         end if

      end subroutine


end module MCModule

!************************************************************************
!*                                                                      *
!*     MCDriver                                                         *
!*                                                                      *
!************************************************************************

! ... monte carlo driver

subroutine MCDriver(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IOMC(iStage)
      if (lmcweight) call MCWeightIO(iStage)
      if (lautumb) call UmbrellaIO(iStage)
      if (lmcpmf) call MCPmfIO(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iWriteInput)

      call IOMC(iStage)
      if (lmcweight) call MCWeightIO(iStage)
      if (lautumb) call UmbrellaIO(iStage)
      if (lmcpmf) call MCPmfIO(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iBeforeSimulation)

      if (lcont) call MCAver(iStage)
      if (lautumb) call UmbrellaIO(iStage)
      if (lmcpmf) call MCPmfIO(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iBeforeMacrostep)

      call IOMC(iStage)
      if (lcont) call MCAver(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iSimulationStep)

      call MCPass(iStage)                         ! MCAver(iSimulationStep) is called in MCPass
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iAfterMacrostep)

      if (lcont) call MCAver(iStage)
      if (lmcpmf) call MCPmfIO(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   case (iAfterSimulation)

      if (lcont) call MCAver(iStage)
      if (lautumb) call UmbrellaIO(iStage)
      if (lmcpmf) call MCPmfIO(iStage)
      if (lpspartsso) call SSODriver(iStage)   ! Pascal Hebbeker

   end select

end subroutine MCDriver

!************************************************************************
!*                                                                      *
!*     IOMC                                                             *
!*                                                                      *
!************************************************************************

! ... perform i/o on monte carlo variables

subroutine IOMC(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOMC'
   integer(4) :: ipt, ict, iat, iatloc
   real(8), allocatable :: psum(:)
   real(8) :: drnlist, Getdrnlist
   logical :: lzero

   namelist /nmlMC/ isamp,                                                                       &
                    pspart, dtran, drot,                                                         &
                    lcl1spart, lfixzcoord, lfixxycoord, lshiftzcom, lfixchainstartspart, ispart, &
                    pspartcl2, txmembcl2, radcl2, dtrancl2,                                      &
                    ppivot, txpivot, txdualpivot, drotpivot, drotminpivot, ipivotrotmode, lcl1pivot, &
                    pchain, dtranchain, drotchain, lcl1chain,                                    &
                    pslither,                                                                    &
                    pbrush, dtranbrush, drotbrush, lcl1brush,                                    &
                    pbrushcl2, dtranbrushcl2, drotbrushcl2,                                      &
                    phierarchical, dtranhierarchical,                                            &
                    pnetwork, dtrannetwork,                                                      &
                    pvol, dvol,                                                                  &
                    pnpart, chempot,                                                             &
                    pcharge,                                                                     &
                    radcl1, pselectcl1,                                                          &
                    lmcweight, lautumb, lmcpmf,                                                  &
                    itestmc,                                                                     &
                    pspartsso, lmcsep                                                               ! Pascal Hebbeker

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      if (.not. allocated(lptmove)) then
         allocate(lptmove(npt))
      end if
      if (.not. allocated(pspart)) then
         allocate(pspart(npt))
      end if
      if (.not. allocated(dtran)) then
         allocate(dtran(npt))
      end if
      if (.not. allocated(drot)) then
         allocate(drot(npt))
      end if
      if (.not. allocated(lcl1spart)) then
         allocate(lcl1spart(npt))
      end if
      if (.not. allocated(lfixzcoord)) then
         allocate(lfixzcoord(npt))
      end if
      if (.not. allocated(lfixxycoord)) then
         allocate(lfixxycoord(npt))
      end if
      if (.not. allocated(lshiftzcom)) then
         allocate(lshiftzcom(npt))
      end if
      if (.not. allocated(pspartcl2)) then
         allocate(pspartcl2(npt))
      end if
      if (.not. allocated(txmembcl2)) then
         allocate(txmembcl2(npt))
      end if
      if (.not. allocated(radcl2)) then
         allocate(radcl2(npt))
      end if
      if (.not. allocated(dtrancl2)) then
         allocate(dtrancl2(npt))
      end if
      if (.not. allocated(ppivot)) then
         allocate(ppivot(npt))
      end if
      if (.not. allocated(txpivot)) then
         allocate(txpivot(npt))
      end if
      if (.not. allocated(drotpivot)) then
         allocate(drotpivot(npt))
      end if
      if (.not. allocated(drotminpivot)) then
         allocate(drotminpivot(npt))
      end if
      if (.not. allocated(lcl1pivot)) then
         allocate(lcl1pivot(npt))
      end if
      if (.not. allocated(pchain)) then
         allocate(pchain(npt))
      end if
      if (.not. allocated(dtranchain)) then
         allocate(dtranchain(npt))
      end if
      if (.not. allocated(drotchain)) then
         allocate(drotchain(npt))
      end if
      if (.not. allocated(lcl1chain)) then
         allocate(lcl1chain(npt))
      end if
      if (.not. allocated(pslither)) then
         allocate(pslither(npt))
      end if
      if (.not. allocated(pbrush)) then
         allocate(pbrush(npt))
      end if
      if (.not. allocated(dtranbrush)) then
         allocate(dtranbrush(npt))
         dtranbrush = 0.0E+00
      end if
      if (.not. allocated(drotbrush)) then
         allocate(drotbrush(npt))
         drotbrush = 0.0E+00
      end if
      if (.not. allocated(lcl1brush)) then
         allocate(lcl1brush(npt))
         lcl1brush = .false.
      end if
      if (.not. allocated(pbrushcl2)) then
         allocate(pbrushcl2(npt))
      end if
      if (.not. allocated(dtranbrushcl2)) then
         allocate(dtranbrushcl2(npt))
         dtranbrushcl2 = 0.0E+00
      end if
      if (.not. allocated(drotbrushcl2)) then
         allocate(drotbrushcl2(npt))
         drotbrushcl2 = 0.0E+00
      end if
      if (.not. allocated(phierarchical)) then
         allocate(phierarchical(npt))
      end if
      if (.not. allocated(dtranhierarchical)) then
         allocate(dtranhierarchical(npt))
         dtranhierarchical = 0.0E+00
      end if
      if (.not. allocated(pnetwork)) then
         allocate(pnetwork(npt))
      end if
      if (.not. allocated(dtrannetwork)) then
         allocate(dtrannetwork(npt))
         dtrannetwork = 0.0E+00
      end if
      if (.not. allocated(pvol)) then
         allocate(pvol(npt))
      end if
      if (.not. allocated(pnpart)) then
         allocate(pnpart(npt))
      end if
      if (.not. allocated(chempot)) then
         allocate(chempot(npt))
         chempot = 0.0E+00
      end if
      if (.not. allocated(radcl1)) then
         allocate(radcl1(npt))
      end if
      if (.not. allocated(pselectcl1)) then
         allocate(pselectcl1(npt))
      end if
      if (.not. allocated(pcharge)) then
         allocate(pcharge(npt))
      end if
      if (.not. allocated(pspartsso)) then
         allocate(pspartsso(npt))
      end if
      if (.not. allocated(plocal)) then
         allocate(plocal(npt))
      end if

      isamp           = 1

      pspart          = One
      dtran           = Zero
      drot            = Zero
      lcl1spart       =.false.
      lfixzcoord      =.false.
      lfixxycoord     =.false.
      lshiftzcom      =.false.
      lfixchainstartspart =.false.
      ispart          = 0
      pspartcl2       = Zero
      txmembcl2       ='ipt'
      radcl2          = Zero
      dtrancl2        = Zero
      ppivot          = Zero
      txpivot         ='short'
      txdualpivot     ='nodual'
      drotpivot       = 360.d0
      drotminpivot    = Zero
      ipivotrotmode   = 1
      lcl1pivot       =.false.

      pchain          = Zero
      dtranchain      = Zero
      drotchain       = Zero
      lcl1chain       =.false.

      pbrush          = Zero
      pbrushcl2       = Zero

      phierarchical   = Zero

      pnetwork        = Zero

      pslither        = Zero

      pvol            = Zero
      dvol            = Zero

      pnpart          = Zero

      pcharge         = Zero

      pspartsso       = Zero

      radcl1          = Zero
      pselectcl1      = One

      lmcweight       = .false.
      lautumb         = .false.
      lmcpmf          = .false.

      lmcsep          = .false.

      itestmc         = 0

      rewind(uin)
      read(uin,nmlMC)

      do ipt = 1, npt
         call LowerCase(txpivot(ipt))
      end do
      call LowerCase(txdualpivot)


      ! already set lpspartsso: it is neede for the SSO-Driver----------
      if (count(pspartsso(1:npt) > Zero) > 0) then                      ! any prob > 0 ?
         lpspartsso =.true.                                                                             ! set lprob = .true.
      end if
      if(lpspartsso) then
         if(.not. allocated(lssopt)) then
            allocate(lssopt(npt))
         end if
         lssopt = .false.
         where ( pspartsso(1:npt) > Zero) lssopt = .true.
      end if


   case (iWriteInput)

! ... allocate memory for MC-specific pointers

      if (.not.allocated(ianatm)) then
         allocate(ianatm(na_alloc))
         ianatm = 0
      end if
      if (.not.allocated(iptmpn)) then
         allocate(iptmpn(np_alloc))
         iptmpn = 0
      end if
      if (.not.allocated(ipnptm)) then
         allocate(ipnptm(np_alloc))
         ipnptm = 0
      end if

! ... initiate lptmdutwob

      lptmdutwob =.false.

! ... initiate lptm

      if (.not.allocated(lptm)) then
         allocate(lptm(np_alloc))
         lptm = .false.
      end if
      lptm =.false.

! ... initiate lmcsep

      if (lmcsep .and. (.not. allocated(plocal))) then
         allocate(plocal(npt))
      end if

! ... check conditions

      if (.not.lntp) pvol = Zero
      if (.not.lmvt) pnpart = Zero
      if ((ipivotrotmode < 1) .or. (ipivotrotmode > 3)) call stop(txroutine, 'ipivotrotmode out of boundary', uout)
      if (.not.lclink) then
         if(maxval(phierarchical) > Zero) call stop(txroutine, 'hierarchical move is requested but link == .f.',uout)
         if(maxval(pnetwork) > Zero) call stop(txroutine, 'network move is requested but link == .f.',uout)
      end if

      if (lweakcharge) then  ! exclude call of ChargeCharge if no weak charge
         iat = 0
         do ipt = 1, npt
            lzero = .true.
            do iatloc = 1, natpt(ipt)
               iat = iat +1
               if (latweakcharge(iat) .and. (naatpt(iatloc,ipt) > 0)) lzero = .false.  ! at least one weak charge
            end do
            if (lzero) pcharge(ipt) = Zero
         end do
      end if

! ... change representation from ict to ipt       (simplify, provide one entry per particle type in input)

  !   call MCChangeRep(npt,nct,ictpt,ppivot)
  !   call MCChangeRep(npt,nct,ictpt,drotpivot)
  !   call MCChangeRep(npt,nct,ictpt,drotminpivot)
  !   call MCChangeRep(npt,nct,ictpt,pchain)
  !   call MCChangeRep(npt,nct,ictpt,dtranchain)
  !   call MCChangeRep(npt,nct,ictpt,drotchain)
  !   call MCChangeRep(npt,nct,ictpt,pslither)
      if (itestmc == 1) call TestIOMCProb('before normalization', uout)

! ... normalize probabilities of different trial moves

      if (.not.allocated(psum)) then
         allocate(psum(npt))
         psum = 0.0E+00
      end if

      psum(1:npt)  = pspart(1:npt)   + pspartcl2(1:npt) + &
                     ppivot(1:npt)   +                    &
                     pchain(1:npt)   +                    &
                     pslither(1:npt) +                    &
                     pbrush(1:npt)   + pbrushcl2(1:npt) + &
                     phierarchical(1:npt) +               &
                     pnetwork(1:npt)      +               &
                     pvol(1:npt)          +               &
                     pnpart(1:npt)        +               &
                     pcharge(1:npt)       +               &
                     pspartsso(1:npt)                     ! Pascal Hebbeker
      where(psum <= 1.0d-15) psum = 1.0d-15               ! to avoid division by Zero
      pspart(1:npt)        = pspart(1:npt)/psum(1:npt)
      pspartcl2(1:npt)     = pspartcl2(1:npt)/psum(1:npt)
      ppivot(1:npt)        = ppivot(1:npt)/psum(1:npt)
      pchain(1:npt)        = pchain(1:npt)/psum(1:npt)
      pslither(1:npt)      = pslither(1:npt)/psum(1:npt)
      pbrush(1:npt)        = pbrush(1:npt)/psum(1:npt)
      pbrushcl2(1:npt)     = pbrushcl2(1:npt)/psum(1:npt)
      phierarchical(1:npt) = phierarchical(1:npt)/psum(1:npt)
      pnetwork(1:npt)      = pnetwork(1:npt)/psum(1:npt)
      pvol(1:npt)          = pvol(1:npt)/psum(1:npt)
      pnpart(1:npt)        = pnpart(1:npt)/psum(1:npt)
      pcharge(1:npt)       = pcharge(1:npt)/psum(1:npt)
      pspartsso(1:npt)     = pspartsso(1:npt)/psum(1:npt)   ! Pascal Hebbeker
      if (itestmc == 1) call TestIOMCProb('after normalization', uout)

! ... set lptmove which controls if some particles should not be moved at all

      lptmove = .true.
      where(psum(1:npt) <= ddelta) lptmove(1:npt) = .false.

      if (allocated(psum)) deallocate(psum)

! ... angles in radians

      drot(1:npt)            = drot(1:npt)*sclang
      drotpivot(1:npt)       = drotpivot(1:npt)*sclang
      drotminpivot(1:npt)    = drotminpivot(1:npt)*sclang
      drotchain(1:npt)       = drotchain(1:npt)*sclang
      drotbrush(1:npt)       = drotbrush(1:npt)*sclang

! ... set logical flags and check some conditions

      call CheckMCProb('pspart',        pspart,       .false.,    .false., lpspart)
!     call CheckMCProb('pspartcl2',     pspartcl2,     lpolyatom,  lchain, lpspartcl2)
      call CheckMCProb('pspartcl2',     pspartcl2,    .false.,     lchain, lpspartcl2)
      call CheckMCProb('ppivot',        ppivot,       .false.,    .false., lppivot)
      call CheckMCProb('pchain',        pchain,       .false.,    .false., lpchain)
      call CheckMCProb('pslither',      pslither,     .false.,    .false., lpslither)
      call CheckMCProb('pbrush',        pbrush,       .false.,    .false., lpbrush)
      call CheckMCProb('pbrushcl2',     pbrushcl2,    .false.,    .false., lpbrushcl2)
      call CheckMCProb('phierarchical', phierarchical,.false.,    .false., lphierarchical)
      call CheckMCProb('pnetwork',      pnetwork,     .false.,    .false., lpnetwork)
      call CheckMCProb('pvol',          pvol,         .false.,    .false., lpvol)
      call CheckMCProb('pnpart',        pnpart,       .false.,    .false., lpnpart)
      call CheckMCProb('pcharge',       pcharge,      .false.,    .false., lpcharge)
      call CheckMCProb('pspartsso',     pspartsso,    .false.,    .false., lpspartsso)   ! Pascal Hebbeker

! ???
!     if (lpspart .and. lchain .and. count(lcl1spart(1:npt))>0) call Stop(txroutine, 'lpspart .and. lchain .and. lcl1part', uout)
      if (ldieldis .and. count(lcl1spart(1:npt))>0) call Stop(txroutine, 'ldieldis .and. lcl1part', uout)

! ... check chain contour length and box length if pivot rotation (potential pivot trouble)

      if (lppivot .and. txbc == 'xyz') then
         do ict = 1, nct
            if (bond(ict)%eq*1.1*(npct(ict)-1) > boxlenshort)  &
            call Warn(txroutine, 'chain contour length > box length, potential trouble with pbc & pivot', uout)
         end do
      end if

#if defined (_PAR_)
      if (count(pselectcl1(1:npt) > Zero .and. pselectcl1(1:npt) < One) > 0) call Stop(txroutine, 'invalid pselectcl1 for _PAR_ ', uout)
      if (count(pslither(1:npt) > Zero) > 0) call Stop(txroutine, 'pslither > 0 .and. _PAR_ not supported', uout)
      if (lpspartcl2) call Stop(txroutine, 'pspartcl2(ipt) > 0 .and. _PAR_ not supported', uout)
      if (lpbrush) call Stop(txroutine, 'pbrush(ipt) > 0 .and. _PAR_ not supported', uout)
      if (lpbrushcl2) call Stop(txroutine, 'pbrushcl2(ipt) > 0 .and. _PAR_ not supported', uout)
      !if (lphierarchical) call Stop(txroutine, 'phierarchical(ipt) > 0 .and. _PAR_ not supported', uout)
      if (lpnetwork) call Stop(txroutine, 'pnetwork(ipt) > 0 .and. _PAR_ not supported', uout)
      !if (lpcharge) call Stop(txroutine, 'pcharge(ipt) > 0 .and. _PAR_ not supported', uout)
#endif

      if (master) then

! ... check if drnlist is sufficiently large (not complete)

         drnlist = Getdrnlist()
         if (lvlist) then
            if (lpspart) then
               if (count(lcl1spart(1:npt)) == 0) then
                  if (maxval(dtran(1:npt)/2) > drnlist) call Warn(txroutine, 'maxval(dtran(1:npt)/2) > drnlist', uout)
               else
                  if (maxval(dtran(1:npt)/2+radcl1(1:npt)) > drnlist) &
                     call Warn(txroutine, 'maxval(dtran(1:npt)/2+radcl1(1:npt)) > drnlist', uout)
               end if       !!! above ?? !!
            end if
            if (lpspartcl2     .and. count(dtrancl2(1:npt)/2+radcl2(1:npt) > drnlist) > 0) &
                                         call Warn(txroutine, 'dtrancl2(ipt)/2+radcl2(ipt) > drnlist', uout)
            if (lpchain        .and. count(dtranchain(1:npt)/2 > drnlist) > 0) &
                                         call Warn(txroutine, 'dtranchain(ipt)/2 > drnlist', uout)
            if (lpbrush        .and. count(dtranbrush(1:npt)/2 > drnlist) > 0) &
                                         call Warn(txroutine, 'dtranbrush(ipt)/2 > drnlist', uout)
            if (lpbrushcl2     .and. count(dtranbrushcl2(1:npt)/2+radcl1(1:npt) > drnlist) > 0) &
                                         call Warn(txroutine, 'dtranbrushcl2(ipt)/2+radcl1(ipt) > drnlist', uout)
            if (lphierarchical .and. count(dtranhierarchical(1:npt)/2 > drnlist) > 0) &
                                         call Warn(txroutine, 'dtranhierarchical(ipt)/2 > drnlist', uout)
            if (lpnetwork      .and. count(dtrannetwork(1:npt)/2 > drnlist) > 0) &
                                         call Warn(txroutine, 'dtrannetwork(ipt)/2 > drnlist', uout)
         end if

         call WriteHead(2, 'mc data', uout)
         if (lmcsep) write(uout,'(a)') 'separating move types (lmcsep)'
         if (isamp == 0) write(uout,'(a)') 'uniform sampling, sequence'
         if (isamp == 1) write(uout,'(a)') 'uniform sampling, random'
         if (dtran(1) < 0) write(uout,'(a)') 'spherical displacement area'
         if (dtran(1) > 0) write(uout,'(a)') 'cubic displacement volume'
         write(uout,'()')
         if (lpspart) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ispartmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pspart(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtran(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drot(1:npt)/sclang
            if (count(lcl1spart(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' cluster1 move                 =', lcl1spart(1:npt)
            if (count(lfixzcoord(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' fixed z-coordinate           =', lfixzcoord(1:npt)
            if (count(lfixxycoord(1:npt)) > 0)write(uout,'(a,t40,6(l10,5x))') ' fixed xy-coordinates         =', lfixxycoord(1:npt)
            if (count(lshiftzcom(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' z_com = 0 if possible        =', lshiftzcom(1:npt)
            if (lfixchainstartspart) write(uout,'(a)') ' no displacement of first chain particle'
            if (ispart > 0) write(uout,'(a,t35,i10)') 'ispart                         = ', ispart
            write(uout,'()')
         end if
         if (lpspartcl2) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ispartcl2move), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pspartcl2(1:npt)
            write(uout,'(a,t40,6(a10  ,5x))') 'selction of cluster members    = ', txmembcl2(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') ' radius of secondary cluster   = ', radcl2(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtrancl2(1:npt)
            if (count(lfixzcoord(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' fixed z-coordinate           =', lfixzcoord(1:npt)
            if (count(lfixxycoord(1:npt)) > 0)write(uout,'(a,t40,6(l10,5x))') ' fixed xy-coordinates         =', lfixxycoord(1:npt)
            write(uout,'()')
         end if
         if (lppivot) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ipivotmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', ppivot(1:npt)
            if (txdualpivot /= 'nodual') write(uout,'(a,t40,6(a10  ,5x))') 'dual pivot move                = ', txdualpivot
            write(uout,'(a,t40,6(a10  ,5x))') 'subchain to be rotated         = ', txpivot(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drotpivot(1:npt)/sclang
            write(uout,'(a,t40,6(f10.3,5x))') 'auxillary displacement param.  = ', drotminpivot(1:npt)/sclang
            write(uout,'(a,t35,i6)')          'pivot rotation mode            = ', ipivotrotmode
            if (count(lcl1pivot(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' cluster1 move                 =', lcl1pivot(1:npt)
            write(uout,'()')
         end if
         if (lpchain) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ichainmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pchain(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtranchain(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drotchain(1:npt)/sclang
            if (count(lcl1chain(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' cluster1 move                 =', lcl1chain(1:npt)
            write(uout,'()')
         end if
         if (lpslither) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(islithermove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pslither(1:npt)
            write(uout,'()')
         end if
         if (lpbrush) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ibrushmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pbrush(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtranbrush(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drotbrush(1:npt)/sclang
            if (count(lcl1brush(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' cluster1 move                 =', lcl1brush(1:npt)
            if (count(lshiftzcom(1:npt)) > 0) write(uout,'(a,t40,6(l10,5x))') ' z_com = 0 if possible        =', lshiftzcom(1:npt)
            write(uout,'()')
         end if
         if (lpbrushcl2) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ibrushcl2move), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pbrushcl2(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtranbrushcl2(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drotbrushcl2(1:npt)*RadToDeg
            write(uout,'()')
         end if
         if (lphierarchical) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ihierarchicalmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', phierarchical(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtranhierarchical(1:npt)
            write(uout,'()')
         end if
         if (lpnetwork) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(inetworkmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pnetwork(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtrannetwork(1:npt)
            write(uout,'()')
         end if
         if (lpvol) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ivolumemove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pvol(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') 'maximum volume change          = ', dvol
            write(uout,'()')
         end if
         if (lpnpart) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(inpartmove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pnpart(1:npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'chemical potential             = ', chempot(1:npt)
            write(uout,'()')
         end if
         if (lpcharge) then
            write(uout,'(a,t35,6(5x,a))') txmovetype(ichargemove), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pcharge(1:npt)
            write(uout,'()')
         end if
         if (lpspartsso) then   ! Pascal Hebbeker
            write(uout,'(a,t35,6(5x,a))') txmovetype(ispartsso), txpt(1:npt)
            write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
            write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pspartsso(1:npt)
            write(uout,'()')
         end if
         if ((count(lcl1spart(1:npt))>0) .or. lpspartcl2 .or.                      &
            (count(lcl1pivot(1:npt))>0) .or.                                      &
            (count(lcl1chain(1:npt))>0) .or.                                      &
            (count(lcl1brush(1:npt))>0) .or. lpbrushcl2) then
            write(uout,'()')
            write(uout,'(a,t35,6(f10.3,5x))') ' radius of primary cluster     = ', radcl1(1:npt)
            write(uout,'(a,t40,6(f10.3,5x))') ' particle selection prob.      = ', pselectcl1(1:npt)
         end if

         if (lmcweight) write(uout, '(a)') 'apply weighting function based on the position of two particles'
         if (lautumb) write(uout,'(a)') 'generate umbrella potential with automatic update procedure'

      end if

#if defined (_PAR_)
      if (lautumb) call Stop(txroutine, 'generation of umbreall potential not adapted for _PAR_', uout)
#endif

! ... modify the representation of the probability for their use in MCPass

! ... first the local moves
      pspartsso(1:npt)      = pspartsso(1:npt)      + pspart(1:npt)
      if (lmcsep) then
         plocal(1:npt) = pspartsso(1:npt)       !the probability to perform a local move
      end if

! ... then the non-local moves
      pspartcl2(1:npt)      = pspartcl2(1:npt)      + pspartsso(1:npt)
      ppivot(1:npt)         = ppivot(1:npt)         + pspartcl2(1:npt)
      pchain(1:npt)         = pchain(1:npt)         + ppivot(1:npt)
      pslither(1:npt)       = pslither(1:npt)       + pchain(1:npt)
      pbrush(1:npt)         = pbrush(1:npt)         + pslither(1:npt)
      pbrushcl2(1:npt)      = pbrushcl2(1:npt)      + pbrush(1:npt)
      phierarchical(1:npt)  = phierarchical(1:npt)  + pbrushcl2(1:npt)
      pnetwork(1:npt)       = pnetwork(1:npt)       + phierarchical(1:npt)
      pvol(1:npt)           = pvol(1:npt)           + pnetwork(1:npt)
      pnpart(1:npt)         = pnpart(1:npt)         + pvol(1:npt)
      pcharge(1:npt)        = pcharge(1:npt)        + pnpart(1:npt)

      if (itestmc == 1) call TestIOMCProb('after modification', uout)

   case (iBeforeMacrostep)

      if (txuser == 'sim_annealing') then
         call WriteHead(2, 'simulated annealing', uout)
         dtran = dtran*0.99
         drot  = drot *0.99
         write(uout,'(a,t35,6(5x,a))') txmovetype(ispartmove), txpt(1:npt)
         write(uout,'(a,t35,6(5x,a))') '-------------------', ('----------',ipt = 1,npt)
         write(uout,'(a,t35,6(f10.3,5x))') 'probability                    = ', pspart(1:npt)
         write(uout,'(a,t40,6(f10.3,5x))') 'tran. displacement parameter   = ', dtran(1:npt)
         write(uout,'(a,t40,6(f10.3,5x))') 'rot. displacement parameter    = ', drot(1:npt)/sclang
      end if

   end select

contains

!........................................................................

subroutine MCChangeRep(npt, nct, ictpt, var)  ! changing representation from nct to npt   (obsolete)

   implicit none

   integer(4), intent(in)    :: npt
   integer(4), intent(in)    :: nct
   integer(4), intent(in)    :: ictpt(*)
   real(8),    intent(inout) :: var(*)

   integer(4) :: ipt, ict
   real(8)    :: varsave(nct)

   varsave(1:nct) = var(1:nct)
   var(1:npt) = 0.0d0
   do ipt = 1, npt
     ict = ictpt(ipt)
     if (ict > 0) var(ipt) = varsave(ict)
   end do

end subroutine MCChangeRep

!........................................................................

subroutine CheckMCProb(label, prob, lpolyatom, lchain, lprob)   ! set logical flag and check some conditions
   character(*), intent(in) :: label
   real(8), intent(in)      :: prob(*)
   logical, intent(in)      :: lpolyatom
   logical, intent(in)      :: lchain
   logical, intent(out)     :: lprob
   lprob =.false.
   if (count(prob(1:npt) > Zero) > 0) then                      ! any prob > 0 ?
      lprob =.true.                                                                             ! set lprob = .true.
      if (count(prob(1:npt) < Zero) > 0) call Stop(txroutine, trim(label)//'(ipt) < Zero', uout)! do not allow any prob < 0
      if (lpolyatom) call Stop(txroutine, trim(label)//' and polyatom is not supported', uout)  ! check if polyatom
      if (lchain)    call Stop(txroutine, trim(label)//' and lchain is not supported', uout)    ! check if lchain
   end if
end subroutine CheckMCProb

!........................................................................

subroutine TestIOMCProb(txheading,unit)   ! test output of probabilities of MC moves
   implicit none
   character(*) :: txheading
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//' '//txheading, unit)
   write(unit,'(a,6f10.3)') 'pspart(1:npt)        ', pspart(1:npt)
   write(unit,'(a,6f10.3)') 'pspartsso(1:npt)     ', pspartsso(1:npt)
   write(unit,'(a,6f10.3)') 'pspartcl2(1:npt)     ', pspartcl2(1:npt)
   write(unit,'(a,6f10.3)') 'ppivot(1:npt)        ', ppivot(1:npt)
   write(unit,'(a,6f10.3)') 'pchain(1:npt)        ', pchain(1:npt)
   write(unit,'(a,6f10.3)') 'pslither(1:npt)      ', pslither(1:npt)
   write(unit,'(a,6f10.3)') 'pbrush(1:npt)        ', pbrush(1:npt)
   write(unit,'(a,6f10.3)') 'pbrushcl2(1:npt)     ', pbrushcl2(1:npt)
   write(unit,'(a,6f10.3)') 'phierarchical(1:npt) ', phierarchical(1:npt)
   write(unit,'(a,6f10.3)') 'pnetwork(1:npt)      ', pnetwork(1:npt)
   write(unit,'(a,6f10.3)') 'pvol(1:npt)          ', pvol(1:npt)
   write(unit,'(a,6f10.3)') 'pnpart(1:npt)        ', pnpart(1:npt)
   write(unit,'(a,6f10.3)') 'pcharge(1:npt)       ', pcharge(1:npt)
end subroutine TestIOMCProb

end subroutine IOMC

!************************************************************************
!*                                                                      *
!*     MCPass                                                           *
!*                                                                      *
!************************************************************************

! ... perform one mc pass (np trial moves)

subroutine MCPass(iStage)

   use MCModule
   use NListModule, only : drnlist, drosum
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCPass'
   integer(4) :: ipt, ict
   real(8)    :: Random, prandom, drnold, rchain
   logical :: lnonloc

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   drostep= Zero

   if(lmcsep) then ! call only local or non-local moves, not both
      prandom = Random(iseed)
      lnonloc = .false.
      if (any( prandom > plocal(1:npt) )) then !only also non-local moves are present
         lnonloc = .true.
      end if

      if(lnonloc) then
         drosum = Zero
         drnold = drnlist
         rchain = Zero
         do ipt = 1, npt
            ict = ictpt(ipt)
            if ((prandom > plocal(ipt)) .and. (ict > 0)) then !non-local moves of this particle type and particle is in chain
               rchain=max(rchain, sum(npptct(1:npt,ict))*bond(ict)%eq)
            end if
         end do

         drnlist = max(4*rchain, nphn*maxval(dtranhierarchical(:))) + drnold !set drnlist to four times contour length + drnlist (old)
         ! four times as 1 particle can move at max 2 times contour length using pivot move, and therefore two particles can approach each at max 4 times the contour length
         ! added drnold to reflect any possible local moves

         !get neighbor list
         if (lvlist) then
            call SetVList
            call VListAver(iStage)
         end if
         if (lllist) then
            call SetLList(rcut+drnlist)
            call LListAver(iStage)
         end if
      end if
   end if

   do ipass = 1, np

      if (ltrace) call WriteTrace(3, txroutine//' (start of new ipass)', iStage)

! ... select a particle

      if (isamp == 0) ipmove = ipass
      if (isamp == 1) ipmove = 1+int(np*Random(iseed))
      ipmove = max(1,int(min(ipmove,np)))
      iptmove = iptpn(ipmove)


! ... check if particle should be moved

      if (.not.lptmove(iptmove)) cycle

! ... select a trial move method

      if (pspart(iptmove) > One - 1d-15) then
         call SPartMove(iStage)
      else
         if(.not. lmcsep) prandom = Random(iseed)
         call CallMove(prandom, iptmove, iStage)
      end if

      call Restorelptm(nptm, ipnptm, lptm)                     ! restore lptm

      if (itest == 1) call TestSimulation
      if (lcont) call MCAver(iSimulationStep)

   end do

   if(lmcsep) then
      ! restore neighbour list
      if(lnonloc) then
         drnlist = drnold
         drosum = Zero

         if (lvlist) then
            call SetVList
            call VListAver(iStage)
         end if
         if (lllist) then
            call SetLList(rcut+drnlist)
            call LListAver(iStage)
         end if
      end if
   end if

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine MCPass

!************************************************************************
!*                                                                      *
!*     SPartMove                                                        *
!*                                                                      *
!************************************************************************

! ... perform one single-particle  trial move

subroutine SPartMove(iStage, loptsso)

   use MolModule
   use MCModule, only: imovetype, ispart
   use MCModule, only: ispartmove, ipmove, iptmove
   use MCModule, only: ievent, imcaccept
   use MCModule, only: itestmc, lautumb, lfixchainstartspart, lmcpmf
   use MCModule, only: npclnew, npclold, radcl1
   use MCModule, only: lcl1spart, pselectcl1, dtran, curdtranpt
   use MCModule, only: lfixzcoord, lfixxycoord, drot, lshiftzcom

   implicit none

   integer(4), intent(in) :: iStage
   logical, optional, intent(in) :: loptsso

   character(40), parameter :: txroutine ='SPartMove'
   character(40), parameter :: txroutinesso ='SSO'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc, dnpcl
   real(8)    :: weight, MCWeight, UmbrellaWeight, MCPmfWeight

   integer(4) :: ipsurf
   integer(4) :: ihost

   logical  :: lsso

   real(8)  :: dtr

   if (ltrace) call WriteTrace(3, txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   !sneak in sso---------------------------------------------------------------
   if( present(loptsso)) then
      lsso = loptsso
   else
      lsso = .false.
   end if

   if(lsso) then
      if (ltrace) call WriteTrace(4, txroutinesso, iStage)
   end if
   !---------------------------------------------------------------------------


   imovetype = ispartmove

#if defined (_PAR_)
! ... the following is not (yet) allowed since iseed will not be syncronized
   if (pselectcl1(iptmove) < One) call Stop(txroutine, '_PAR_ .and. pselectcl1(iptmove) < One', uout)
   ipnptm(1:np) = 0                                 ! necessary for later allreduce
#endif

! .............. define particle(s) to be moved ..............

! ... consider first particle ipmove

   nptm          = 1
   iploc         = 1
   ipnptm(iploc) = ipmove
   iptmpn(ipmove) = iploc
   lptm(ipmove)  =.true.

! ... and then the remaining particles

   if (lcl1spart(iptmove)) then
      if (lvlist) call ClusterMember('old', .false., .false., radcl1, pselectcl1)
      if (lllist) call ClusterMemberLList('old', .false., .false., radcl1, pselectcl1)
   end if

! .. get displacement parameter
   dtr=dtran(iptmove)
   if(lsso) then
      dtr=-curdtranpt(iptmove) !sso always produces positiv values in curdtranpt; but requires sperical displacement volume. GetRandomTrialPos uses spherical displacement volume id dtr<0
   end if


! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtr, iseed, nptm, ipnptm, ro, rotm, drotm)
   if (lfixzcoord(iptmove)) call FixedZCoord
   if (lfixxycoord(iptmove)) call FixedXYCoord
   if (lfixchainstartspart) call FixedChainStart

!-------------------------------------------------------------------------------------

   if (ispart == 1) then                         ! for pmf calculations
      ipsurf = 1
      if (iptmove == 2) call SurfacePart(ipsurf) ! place particle ipmove on the surface of ip = 1
      ipsurf = 2
      if (iptmove == 3) call SurfacePart(ipsurf) ! place particle ipmove on the surface of ip = 2
   end if

   if (ispart == 2) then                         ! for macroions (ipt = 1) with mobile surface charges (ipt = 2)
      if (iptmove == 2) then
         ihost = 1 + (ipmove-ipnpt(iptmove))/(nppt(2)/nppt(1))
         call SurfacePart(ihost)                ! place particle ipmove on the surface of ip = ihost
      end if
   end if

   if (ispart == 3) call SymmMove                ! symmetric z-move of particle 1 and 2

!-------------------------------------------------------------------------------------

   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lcl1spart(iptmove)) then
      if (lpolyatom) call GetFixedTrialOri(nptm, ipnptm, ori, oritm)
   end if
   if (lpolyatom .or. lellipsoid .or. lsuperball) &
      call GetRandomTrialOri(drot(iptmove), iseed, ori(1,1,ipmove), oritm(1,1,iploc))

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (.not. lboxoverlap) then

! ............. evaluate energy difference ...............
      call DUTotal(lhsoverlap, lhepoverlap)

      if (.not. (lhsoverlap .or. lhepoverlap)) then
! ............. calculate nonenergetic weights .............
         weight = One
         dnpcl = Zero
         if (lcl1spart(iptmove)) then
            if (lvlist) call ClusterMember('new', .false., .false., radcl1, pselectcl1)       ! calculate npclnew
            if (lllist) call ClusterMemberLList('new', .false., .false., radcl1, pselectcl1)  ! calculate npclnew
            dnpcl = npclnew-npclold
            if (dnpcl /= 0) weight = weight*(One-pselectcl1(iptmove))**dnpcl
         end if

         if (lmcweight) weight = weight*MCWeight()
         if (lautumb) weight = weight*UmbrellaWeight(1)
         if (lmcpmf) weight = weight*MCPmfWeight(1)
      end if

   end if

! ............. decide new configuration .............

   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............
   if(lsso) then
      do iploc = 1, nptm
         call SSOUpdate((ievent == imcaccept), iptmove, drotm(1:3,iploc))
      end do
   end if
   if (ievent == imcaccept) call MCUpdate       ! update energies and coordinates

   if (lautumb) call UmbrellaUpdate              ! update weight function for umbrella potential
   if (lmcpmf) call MCPmfUpdate                  ! update weight function for mc pmf
   if (lshiftzcom(iptmove)) call ShiftZDirection
   if (ispart == 3) lptmdutwob = .false.

contains

!........................................................................

subroutine SurfacePart(iref)  ! place particle iploc on the surface of particle iref
   integer(4), intent(in) :: iref
   real(8) :: dx, dy, dz, r1, r2, norm

! ... get the separation vector between iref and iploc

   dx = ro(1,iref) - rotm(1,iploc)
   dy = ro(2,iref) - rotm(2,iploc)
   dz = ro(3,iref) - rotm(3,iploc)
   call PBCr2(dx,dy,dz,r2)
   r1 = sqrt(r2)

! ... get the fraction of dx, dy, dz to be used for placing iploc on the surface of iref

   norm = One -(One+1d-10)*(radat(iptmove)+radat(iptpn(iref))) / r1

! ... limit the displacement to 0.8 of the rcut

   if (norm*r1 > rcut) norm = 0.8d0*rcut/r1

! ... new trial move of iploc

   rotm(1,iploc) = rotm(1,iploc) + dx*norm
   rotm(2,iploc) = rotm(2,iploc) + dy*norm
   rotm(3,iploc) = rotm(3,iploc) + dz*norm

   call PBC(rotm(1,iploc),rotm(2,iploc),rotm(3,iploc))

   drotm(1,iploc) = rotm(1,iploc) - ro(1,ipmove)
   drotm(2,iploc) = rotm(2,iploc) - ro(2,ipmove)
   drotm(3,iploc) = rotm(3,iploc) - ro(3,ipmove)
end subroutine SurfacePart

!........................................................................

subroutine SymmMove
   character(40), parameter :: txroutine ='SymmMove'
   if (lcl1spart(iptmove)) call Stop(txroutine, 'lcl1spart', uout)
   nptm = 2                          ! two particle to move
   ipnptm(2) = 3 - ipnptm(1)         ! id of the second particle
   lptm(ipnptm(2)) = .true.          ! mark move
   rotm(1:3,2) = -rotm(1:3,1)        ! trial coordinate of the second particle
   lptmdutwob = .true.               ! engage energy evaluation "inside" the moving group

    write(*,*) 'nptm', nptm
    write(*,*) 'ipnptm(1:nptm)', ipnptm(1:nptm)
    write(*,*) 'lptm(1:np)', lptm(1:np)
    write(*,*) 'rotm(1:3,1)', rotm(1:3,1)
    write(*,*) 'rotm(1:3,2)', rotm(1:3,2)

end subroutine SymmMove

!........................................................................

end subroutine SPartMove

!************************************************************************
!*                                                                      *
!*     SPartCl2Move                                                     *
!*                                                                      *
!************************************************************************

! ... perform one single-particle + cluster2 trial move

subroutine SPartCl2Move(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SPartCl2Move'
   integer(4) :: mnpair = 10000
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc, iplow, ipupp, jp
   integer(4) :: npair, dnpcl, npcl
   integer(4), allocatable :: n1(:), n2(:), iobjcluster(:), icllis(:)
   real(8)    :: weight, r2, dx, dy, dz
   integer(4) :: ipshift

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

#if defined (_PAR_)
   if (.not.master) call Stop(txroutine, 'SPartCl2Move not addapted for _PAR_', uout)
#endif

   imovetype = ispartcl2move

   iplow = ipnpt(iptmove)
   ipupp = ipnpt(iptmove)+nppt(iptmove)-1

! .............. define particle(s) to be moved ..............

! ... consider first primary cluster particles

   if (.not.allocated(n1)) then
      allocate(n1(mnpair), n2(mnpair), iobjcluster(mnpair), icllis(mnpair))
      n1 = 0
      n2 = 0
      iobjcluster = 0
      icllis = 0
   end if

   if (txmembcl2(iptmove) == 'ipt') then   ! selection among particles of type iptmove including particle ipmove
      call CalcPartPairListIpt(iptmove,radcl2(iptmove)**2,mnpair,npair,n1,n2) ! list of particle pair for particles of type iptmove
      call MakeCluster(nppt(iptmove), npair, n1, n2, iobjcluster)
      call CalcClusterMember(nppt(iptmove), iobjcluster, ipmove-ipnpt(iptmove)+1, npcl, icllis)
      ipshift = ipnpt(iptmove) - 1
   else if (txmembcl2(iptmove) == 'all') then ! selection among all particles
      call CalcPartPairListAll(np, radcl2(iptmove)**2, mnpair, npair, n1, n2)
      call MakeCluster(np, npair, n1, n2, iobjcluster)
      call CalcClusterMember(np, iobjcluster, ipmove, npcl, icllis)
      ipshift = 0
   end if

   nptm = 0
   do iploc = 1, npcl
      nptm = nptm + 1
      ipnptm(nptm) = icllis(iploc) + ipshift
      lptm(ipnptm(nptm)) = .true.
   end do
   nprimpartcl = nptm

   deallocate(n1, n2, iobjcluster, icllis)

! ... and then the remaining particles (if any)

   if (radcl1(iptmove) > Zero) then
      if (lvlist) call ClusterMember('old', .false., .false., radcl1, pselectcl1)
      if (lllist) call ClusterMemberLList('old', .false., .false., radcl1, pselectcl1)
   end if

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtrancl2(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   if (lfixzcoord(iptmove)) call FixedZCoord
   if (lfixxycoord(iptmove)) call FixedXYCoord
   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lpolyatom) call GetFixedTrialOri(nptm, ipnptm, ori, oritm)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)
   if (lhsoverlap .or. lhepoverlap) goto 200

! ... discard moves where nprimpartcl would change

   do iploc = 1, nprimpartcl
      do jp = iplow, ipupp
         if (lptm(jp)) cycle                  ! consider only particles not already in the cluster
         dx = ro(1,jp)-rotm(1,iploc)
         dy = ro(2,jp)-rotm(2,iploc)
         dz = ro(3,jp)-rotm(3,iploc)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < radcl2(iptmove)**2) then
           ievent = imcreject
           goto 300
         end if
      end do
   end do

! ............. calculate nonenergetic weights .............

   weight = One
   dnpcl = Zero
   if (radcl1(iptmove) > Zero) then ! take into account the change of the number of remaining particles
      if (lvlist) call ClusterMember('new', .false., .false., radcl1, pselectcl1)       ! calculate npclnew
      if (lllist) call ClusterMemberLList('new', .false., .false., radcl1, pselectcl1)  ! calculate npclnew
      dnpcl = npclnew-npclold
      if (dnpcl/= 0) weight = (One-pselectcl1(iptmove))**dnpcl
   end if
! ............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

300 continue

! .............. update .............

   if (ievent == imcaccept) call MCUpdate        ! update energies and coordinates
!old   if (ievent == imcaccept .and. lllist) call IntListDriver(4)

end subroutine SPartCl2Move

!************************************************************************
!*                                                                      *
!*     PivotMove                                                        *
!*                                                                      *
!************************************************************************

! ... perform pivot rotation trial move

subroutine PivotMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='PivotMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iseg, isegpivot, iseg1, iseg2, ic, ict, ip, iploc, ip1, ip2, ip3, dnpcl
   integer(4) :: ip4, iseg3, icl, igen, nptmfirst, nptmlast
   real(8)    :: alpha, rotaxis(3), dx, dy, dz, dxpbc, dypbc, dzpbc, Random, weight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (icnpn(ipmove) == 0) call Stop(txroutine, 'pivot rotation based on a non-chain particle', uout)

   imovetype = ipivotmove

   isegpivot = isegpn(ipmove)    ! segment number isegpivot of pivot particle ipmove
   ic = icnpn(ipmove)            ! chain number ic where particle ipmove is residing
   ict = ictcn(ic)               ! chain type ict of chain ic

! ... early exit if first or last particle in the chain

   if (isegpivot == 1 .or. isegpivot == npct(ictpn(ipmove))) then
      ievent = 0
      return
   end if

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

! .............. define particle(s) to be moved ..............

! ... get the particle(s) defining the rotation axis and the chain particle(s) to be rotated

   ip1 = ipnsegcn(isegpivot,ic)                    ! pivot particle (ipmove)
   if (txpivot(iptmove) == 'short') then
      if (isegpivot <= npct(ictpn(ipmove))/2) then  ! 'lower' part of the chain
         ip2 = ipnsegcn(isegpivot+1,ic)               ! nonmoving neigbour (ipmove+1)
         ip3 = ipnsegcn(isegpivot-1,ic)               ! moving neigbour (ipmove-1)
         iseg1 = 1                                    ! first particle in rotating subchain
         iseg2 = isegpivot                            ! last particle in rotating subchain
      else                                         ! 'upper' part of the chain
         ip2 = ipnsegcn(isegpivot-1,ic)               ! nonmoving neigbour (ipmove+1)
         ip3 = ipnsegcn(isegpivot+1,ic)               ! moving neigbour (ipmove-1)
         iseg1 = isegpivot                            ! first particle in rotating subchain
         iseg2 = npct(ictpn(ipmove))                  ! last particle in rotating subchain
      end if
   else if (txpivot(iptmove) == 'lower') then
      ip2 = ipnsegcn(isegpivot+1,ic)               ! nonmoving neigbour (ipmove+1)
      ip3 = ipnsegcn(isegpivot-1,ic)               ! moving neigbour (ipmove-1)
      iseg1 = 1                                    ! first particle in rotating subchain
      iseg2 = isegpivot                            ! last particle in rotating subchain
   else if (txpivot(iptmove) == 'upper') then
      ip2 = ipnsegcn(isegpivot-1,ic)               ! nonmoving neigbour (ipmove+1)
      ip3 = ipnsegcn(isegpivot+1,ic)               ! moving neigbour (ipmove-1)
      iseg1 = isegpivot                            ! first particle in rotating subchain
      iseg2 = npct(ictpn(ipmove))                  ! last particle in rotating subchain
   else
      call Stop(txroutine, 'error in txpivot', uout)
   end if

! ... store chain particles to be rotated

   nptm = 0
   do iseg = iseg1, iseg2
      ip = ipnsegcn(iseg,ic)
      nptm = nptm+1
      ipnptm(nptm) = ip
      iptmpn(ip) = nptm
      lptm(ip) = .true.
   end do
   call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)   ! restore the chain with UndoPBCChain

!  ... if chain is crosslinked, move also particles in chains of higher generation connected

   if (lclink) then
      nptmfirst = nptm + 1                                                ! used for the next generation
      do iseg = iseg1, iseg2                                              ! loop over main chain segments
         ip = ipnsegcn(iseg,ic)                                           ! particle number of main chain segment
         do icl = 1, nbondcl(ip)                                          ! loop over number of crosslinks
            if ((genic(icnpn(bondcl(icl,ip))) < genic(ic)) .and. (txpivot(iptmove) == 'upper')) cycle           ! ignore crosslinks to lower generations if txpivot = upper
            call UndoPBCChain(vaux(1,ip),icnpn(bondcl(icl,ip)), 1, vaux)  ! restore the chain with UndoPBCChain
            do iseg3 = 1, sum(npptct(1:npt,ictpn(bondcl(icl,ip))))        ! loop over segments of relevant chain type
               ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ip)))                ! particle number in side chain
               nptm = nptm+1
               ipnptm(nptm) = ip4
               iptmpn(ip4) = nptm
               lptm(ip4) = .true.
            end do
         end do
      end do
      nptmlast = nptm                                                                  ! used for next generation
      do igen = merge(genic(ic),0,(txpivot(iptmove) == 'upper'))+1, ngen
         if (igen == genic(ic)) cycle
         do iploc = nptmfirst, nptmlast
            do icl=1, nbondcl(ipnptm(iploc))
               if (genic(icnpn(bondcl(icl,ipnptm(iploc)))) < igen) cycle               ! ignore crosslinks to lower generations
               call UndoPBCChain(vaux(1,ip),icnpn(bondcl(icl,ipnptm(iploc))), 1, vaux) ! restore the chain with UndoPBCChain
               do iseg3=1, sum(npptct(1:npt,ictpn(bondcl(icl,ipnptm(iploc)))))         ! loop over segments of relevant chain type
                  ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ipnptm(iploc))))
                  nptm = nptm+1
                  ipnptm(nptm) = ip4
                  iptmpn(ip4) = nptm
                  lptm(ip4) = .true.
               end do
            end do
         end do
         nptmfirst = nptmlast+1      ! update for next generation
         nptmlast = nptm             ! update for next generation
      end do
    end if

    if (txdualpivot /= 'nodual') call PivotDual

! ... remaining particles (if cluster1 move)

   if (lcl1pivot(iptmove)) then
      call ClusterMember('old', .false., .true., radcl1, pselectcl1)
      do iploc = nprimpartcl+1, nptm
         ip = ipnptm(iploc)
         vaux(1,ip) = ro(1,ip)
         vaux(2,ip) = ro(2,ip)
         vaux(3,ip) = ro(3,ip)
         dx = vaux(1,ip)-vaux(1,ip1)
         dy = vaux(2,ip)-vaux(2,ip1)
         dz = vaux(3,ip)-vaux(3,ip1)
         call PBC2(dx,dy,dz,dxpbc,dypbc,dzpbc)
         vaux(1,ip) = vaux(1,ip) - dxpbc
         vaux(2,ip) = vaux(2,ip) - dypbc
         vaux(3,ip) = vaux(3,ip) - dzpbc
      end do
   end if

! .............. calculate a trial configuration ...............

! ... get rotation axis

   if (ipivotrotmode == 1) then        ! rotation axis given by the vector joining ip1 and ip2 (bond angle is preserved)
      rotaxis(1) = vaux(1,ip2)-vaux(1,ip1)
      rotaxis(2) = vaux(2,ip2)-vaux(2,ip1)
      rotaxis(3) = vaux(3,ip2)-vaux(3,ip1)
   else if (ipivotrotmode == 2) then   ! rotation axis normal to the plane defined by ip1, ip2, and ip3 (bond angle is preserved)
      rotaxis(1) =+(vaux(2,ip2)-vaux(2,ip1))*(vaux(3,ip3)-vaux(3,ip1))-(vaux(2,ip3)-vaux(2,ip1))*(vaux(3,ip2)-vaux(1,ip1))
      rotaxis(2) =-(vaux(1,ip2)-vaux(1,ip1))*(vaux(3,ip3)-vaux(3,ip1))+(vaux(1,ip3)-vaux(1,ip1))*(vaux(3,ip2)-vaux(1,ip1))
      rotaxis(3) =+(vaux(1,ip2)-vaux(1,ip1))*(vaux(2,ip3)-vaux(2,ip1))-(vaux(1,ip3)-vaux(1,ip1))*(vaux(2,ip2)-vaux(2,ip1))
      if (sum(abs(rotaxis(1:3))) < 1.0d-10) call SphRandom(iseed, rotaxis(1), rotaxis(2), rotaxis(3))   ! exception: ip1, ip2, and ip3 on a straight line
   else if (ipivotrotmode == 3) then   ! random rotation axis (bond angle is not preserved)
      call SphRandom(iseed, rotaxis(1), rotaxis(2), rotaxis(3))
   end if

! ... get rotation angle

   alpha = drotpivot(iptmove)*(Random(iseed)-0.5)
   if (drotminpivot(iptmove) > Zero) then
      alpha = drotminpivot(iptmove)+(drotpivot(iptmove)-drotminpivot(iptmove))*Random(iseed)
      if (Random(iseed) > Half) alpha =-alpha
      alpha = Half*alpha
   end if

! ... rotate and get rotm and drotm

   call RotateSetPart(rotaxis, alpha, vaux(1,ip1), nptm, ipnptm, vaux, rotm, drotm)
   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lpolyatom) call GetRotatedTrialOri(rotaxis, alpha, nptm, ipnptm, ori, oritm)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! .............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)
   if (lhsoverlap .or. lhepoverlap) goto 200

! ............. calculate nonenergetic weights .............

   weight = One
   dnpcl = Zero
   if (lcl1pivot(iptmove)) then
      call ClusterMember('new', .false., .true., radcl1, pselectcl1)       !  calculate npclnew
      dnpcl = npclnew-npclold
      if (dnpcl /= 0) weight = (One-pselectcl1(iptmove))**dnpcl
   end if
! .............. decide new configuration .............

200 continue
!  if (abs(du%tot) < 1d-13) du%tot = Zero  ! Temporary for testing differnt compiler with different roundoff
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate       ! update energies and coordinates

contains

!........................................................................

subroutine PivotDual    ! dual pivot rotation

   real(8), allocatable ::  dtemp(:,:)
   integer(4) :: ict_temp, ic_temp, iseg_loc, iseg_min, jp1, jp3, jp_loc, jp_min, jp_max
   integer(4) :: ih, ict_loc, ic_loc, iseg_ini, iseg_fin, jseg, jp, direction, iseg_move, iseg_nomove , ip_nomove, iseg_lastmove , ip_lastmove
   real(8)    :: dx1, dy1, dz1, dx2, dy2, dz2, dtemp_min(2), r1, r2
   real(8)    :: dinf, r_dist

   if (.not.allocated(dtemp)) then
      allocate(dtemp(np_alloc,2))
      dtemp = 0.0E+00
   end if

   dtemp(1:np,2) = 0.0
   dtemp_min(1:2) = 10000
   iseg_min = 0
! ..............   linear .............

   if (txdualpivot == 'linear') then      ! dual pivot: reference: linear chain (ict=3)
                                          ! 2ns tructure involved: comb structures (ict =1 and ict = 2)
   if(ict==3) then
      do ih = 1, nh                       ! loop over number of hierarchic structures
         do igen = 0, ngen                                                       ! loop over generations
            ict_loc = ictgen(igen)                                               ! chain type
            do ic_loc = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1      ! loop over chains of the structure
               do iseg = 1, npct(ict_loc)                                        ! loop over segments
                  ip = ipnsegcn(iseg,ic_loc)
                  if (ict_loc == 1) then  ! UNDO backbone (first chain type) reference chain ic (linear chain)
                     call UndoPBCChain(vaux(1,ipnsegcn(1,ic)), ic_loc, -1, vaux)
                     dinf = boxlen(1)
                     do jseg = 1, npct(ict)
                        jp = ipnsegcn(jseg,ic)
                        dx = vaux(1,ip) - vaux(1,jp)
                        dy = vaux(2,ip) - vaux(2,jp)
                        dz = vaux(3,ip) - vaux(3,jp)
                        r_dist = dx**2+dy**2+dz**2
                        if (sqrt(r_dist) < dinf) dinf = sqrt(r_dist)
                     end do
                     direction = -1
                     if (dinf > boxlen(1)/3.0)   direction = +1
                     call UndoPBCChain(vaux(1,ipnsegcn(1,ic)), ic_loc, direction, vaux)
                  end if
                  if (nbondcl(ip) == 1.and.ict_loc > 1 ) call UndoPBCChain(vaux(1,bondcl(1,ip)), ic_loc, direction, vaux)  !linked beads belonging to the backbone as reference
               end do
            end do
         end do
      end do

      ict_temp = 1                        ! first chain type, i.e. backbone
      ic_temp = icnct(ict_temp)
      do iseg_loc = 1, npct(ict_temp)     ! estimate distance between all particles of 3rd chain ip1;
         dx1 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip1)
         dy1 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip1)
         dz1 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip1)
         dtemp(iseg_loc,1) = sqrt(dx1**2+dy1**2+dz1**2)
         if (dtemp(iseg_loc,1) < dtemp_min(1))  dtemp_min(1) = dtemp(iseg_loc,1)  ! find the closest segment
    end do
    do iseg_loc = 1, npct(ict_temp)
       if (dtemp_min(1) == dtemp(iseg_loc,1)) iseg_min = iseg_loc    ! store closest segment as iseg_min
    end do

    do iseg_loc = iseg_min-1, iseg_min +1, 2  !  compare distance form the two neighbours of the closest particle ipnsegcn(iseg_min,ic_temp) to particle ip3
       if( (iseg_loc < 1) .or. (iseg_loc > npct(ict_temp)) ) then !if neighbour does not exists
          cycle
       end if
       dx2 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip3)
       dy2 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip3)
       dz2 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip3)
       dtemp(iseg_loc,2) = sqrt(dx2**2+dy2**2+dz2**2)
       if (dtemp(iseg_loc,2) < dtemp_min(2))  dtemp_min(2) = dtemp(iseg_loc,2)
    end do
    do iseg_loc = iseg_min-1, iseg_min +1, 2
       if( (iseg_loc < 1) .or. (iseg_loc > npct(ict_temp)) ) then !if neighbour does not exists
          cycle
       end if
       if (dtemp_min(2) == dtemp(iseg_loc,2)) then
          jp1 =  ipnsegcn(iseg_min,ic_temp) ! store closest particle to ip1 as jp1
          jp3 = ipnsegcn(iseg_loc,ic_temp)  ! store closest particle to ip3 as jp3
       end if
    end do

    if (jp1 > jp3 ) then                    ! jp_min  and jp_max: first and last particle of the linear chain in pivot move
       jp_min = ipnsegcn(1,ic_temp)
       jp_max = jp1
       iseg_ini = 1
       iseg_fin = isegpn(jp1)
    else if (jp1 < jp3 )  then
       jp_min = jp1
       jp_max = ipnsegcn(npct(ict_temp),ic_temp)
       iseg_ini = isegpn(jp1)
       iseg_fin = npct(ictpn(jp_max))
    end if

    do jp_loc = jp_min, jp_max
       nptm = nptm+1
       ipnptm(nptm) = jp_loc
       iptmpn(jp_loc) = nptm
       lptm(jp_loc) = .true.
    end do

    if (lclink) then
       nptmfirst = nptm + 1                                               ! used for the next generation
       do iseg = iseg_ini, iseg_fin                                       ! loop over main chain segments
          ip = ipnsegcn(iseg,ic_temp)                                     ! particle number of main chain segment
          do icl = 1, nbondcl(ip)                                         ! loop over number of crosslinks
             if (genic(icnpn(bondcl(icl,ip))) < genic(ic)) cycle          ! ignore crosslinks to lower generations.
             do iseg3 = 1, sum(npptct(1:npt,ictpn(bondcl(icl,ip))))       ! loop over segments of relevant chain type
                ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ip)))               ! particle number in side chain
                nptm = nptm+1
                ipnptm(nptm) = ip4
                iptmpn(ip4) = nptm
                lptm(ip4) = .true.
             end do
          end do
       end do
       nptmlast = nptm                                                    ! used for next generation
       do igen = genic(ic)+1, ngen
          do iploc = nptmfirst, nptmlast
             do icl=1, nbondcl(ipnptm(iploc))
                if (genic(icnpn(bondcl(icl,ipnptm(iploc)))) < igen) cycle ! ignore crosslinks to lower generations.
                do iseg3=1, sum(npptct(1:npt,ictpn(bondcl(icl,ipnptm(iploc))))) ! loop over segments of relevant chain type
                   ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ipnptm(iploc))))
                   nptm = nptm+1
                   ipnptm(nptm) = ip4
                   iptmpn(ip4) = nptm
                   lptm(ip4) = .true.
                end do
             end do
          end do
          nptmfirst = nptmlast+1      ! update for next generation
          nptmlast = nptm             ! update for next generation
       end do
    end if
   end if
   end if

! ..............   combed   .............

   if (txdualpivot == 'combed') then       ! dual pivot: reference: comb structures (ict =1 and ict = 2)
                                           ! 2ns tructure involved: linear chain (ict=3)
      if(ict == 1) then
         do ic = 1, nc                     ! undo the linear chain
            ict = ictcn(ic)
            ip = ipnsegcn(1,ic)
            if(ihnpn(ip) /= 0) cycle       ! exclude hierarchical structures
            call UndoPBCChain(vaux(1,ipnsegcn(1,1)), ic, -1, vaux)! it assumes ic =1 is the backbone of the hierarchy; direction should be  -1
         end do

         ict_temp = 3                      ! third chain type
         ic_temp = icnct(ict_temp)
         do iseg_loc = 1, npct(ict_temp)   ! estimate distance between all particles of 3rd chain ip1;
            dx1 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip1)
            dy1 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip1)
            dz1 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip1)
            dtemp(iseg_loc,1) = sqrt(dx1**2+dy1**2+dz1**2)
            if (dtemp(iseg_loc,1) < dtemp_min(1))  dtemp_min(1) = dtemp(iseg_loc,1)   !    find the closest segment
         end do
         do iseg_loc = 1, npct(ict_temp)
            if (dtemp_min(1) == dtemp(iseg_loc,1)) iseg_min = iseg_loc  ! store closest segment as iseg_min
         end do
         do iseg_loc = iseg_min-1, iseg_min +1, 2        ! compare distance form the two neighbours of the closest particle ipnsegcn(iseg_min,ic_temp) to particle ip3
            dx2 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip3)
            dy2 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip3)
            dz2 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip3)
            dtemp(iseg_loc,2) = sqrt(dx2**2+dy2**2+dz2**2)
            if (dtemp(iseg_loc,2) < dtemp_min(2))  dtemp_min(2) = dtemp(iseg_loc,2)
         end do
         do iseg_loc = iseg_min-1, iseg_min +1, 2
            if (dtemp_min(2) == dtemp(iseg_loc,2)) then
               jp1 =  ipnsegcn(iseg_min,ic_temp)  ! store closest particle to ip1 as jp1
               jp3 = ipnsegcn(iseg_loc,ic_temp)   ! store closest particle to ip3 as jp3
            end if
         end do

         if (jp1 > jp3 ) then       ! determin first and last particle of the linear chain in pivot move
            jp_min = ipnsegcn(1,ic_temp)
            jp_max = jp1
         else if (jp1 < jp3 )  then
            jp_min = jp1
            jp_max = ipnsegcn(npct(ict),ic_temp)
         end if

         do jp_loc = jp_min, jp_max
            nptm = nptm+1
            ipnptm(nptm) = jp_loc
            iptmpn(jp_loc) = nptm
            lptm(jp_loc) = .true.
         end do
      end if
   end if

! ..............   mullin   .............

    if (txdualpivot == 'mullin') then     ! dual pivot: reference: linear chain (ict=3)
                                          ! 2ns tructure involved: comb structures (ict =1 and ict = 2)  icnct(ict) > 1 (more combs)
      if(ict==3) then
         do ih = 1, nh                    ! loop over number of hierarchic structures
            do igen = 0, ngen             ! loop over generations
               ict_loc = ictgen(igen)     ! chain type
               do ic_loc = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1   ! loop over chains of the structure
                  do iseg = 1, npct(ict_loc)                                     ! loop over segments
                     ip = ipnsegcn(iseg,ic_loc)
                     if (ict_loc == 1) then            ! UNDO backbone (first chain type) reference chain ic (linear chain)
                        call UndoPBCChain(vaux(1,ipnsegcn(1,ic)), ic_loc, -1, vaux)
                        dinf = boxlen(1)
                        do jseg = 1, npct(ict)!   ict ==3
                           jp = ipnsegcn(jseg,ic)
                           dx = vaux(1,ip) - vaux(1,jp)
                           dy = vaux(2,ip) - vaux(2,jp)
                           dz = vaux(3,ip) - vaux(3,jp)
                           r_dist = dx**2+dy**2+dz**2
                           if (sqrt(r_dist) < dinf) dinf = sqrt(r_dist)
                        end do
                        direction = -1
                        if (dinf > boxlen(1)/3.0)   direction = +1
                        call UndoPBCChain(vaux(1,ipnsegcn(1,ic)), ic_loc, direction, vaux)
                     end if
                  end do
               end do
            end do
         end do

         ict_temp = 1                     ! first chain type, i.e. backbone
         do ic_temp = icnct(ict_temp), icnct(ict_temp)-1+ncct(ict_temp)
            iseg_move = 0.0
            if(iseg1 == 1) then
               iseg_nomove = npct(ictpn(ipmove))
               iseg_lastmove = 1
               ip_nomove = ipnsegcn(iseg_nomove,ic)
               ip_lastmove = ipnsegcn(iseg_lastmove,ic)
            else if (iseg2 == npct(ictpn(ipmove)))  then
               iseg_nomove = 1
               iseg_lastmove = npct(ictpn(ipmove))
               ip_nomove = ipnsegcn(iseg_nomove,ic)
               ip_lastmove = ipnsegcn(iseg_lastmove,ic)
            end if
            do iseg_loc = 1, npct(ict_temp)! estimate distance between all particles of 3rd chain ip1;
               dx1 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip1)
               dy1 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip1)
               dz1 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip1)
               dtemp(iseg_loc,1) = sqrt(dx1**2+dy1**2+dz1**2)
               if (dtemp(iseg_loc,1) < dtemp_min(1))  dtemp_min(1) = dtemp(iseg_loc,1)   ! find the closest segment
            end do
            do iseg_loc = 1, npct(ict_temp)
               if (dtemp_min(1) == dtemp(iseg_loc,1)) iseg_min = iseg_loc                ! store closest segment as iseg_min
            end do

            jp1 =  ipnsegcn(iseg_min,ic_temp)                  ! store closest particle to ip1 as jp1
            if(iseg_min /= 1 .or. iseg_min /= npct(ict_temp))  then
               do iseg_loc = iseg_min-1, iseg_min +1, 2        ! compare distance form the two neighbours of the closest particle ipnsegcn(iseg_min,ic_temp) to particle ip3
                  dx2 = vaux(1,ipnsegcn(iseg_loc,ic_temp)) - vaux(1,ip3)
                  dy2 = vaux(2,ipnsegcn(iseg_loc,ic_temp)) - vaux(2,ip3)
                  dz2 = vaux(3,ipnsegcn(iseg_loc,ic_temp)) - vaux(3,ip3)
                  dtemp(iseg_loc,2) = sqrt(dx2**2+dy2**2+dz2**2)
                  if (dtemp(iseg_loc,2) < dtemp_min(2))  then
                     dtemp_min(2) = dtemp(iseg_loc,2)
                     jp3 = ipnsegcn(iseg_loc,ic_temp)          ! store closest particle to ip3 as jp3
                     iseg_move = iseg_loc                      ! store closest segment to ip3 as iseg_move
                  end if
               end do
            elseif(iseg_min == 1) then                           ! end chain
               jp3 = ipnsegcn(2,ic_temp)
            elseif (iseg_min == npct(ict_temp)) then            !  end chain
               jp3 = ipnsegcn(npct(ict_temp)-1,ic_temp)
            end if

            if(iseg_min == 1.or.iseg_min == npct(ict_temp).or.dtemp_min(1)>3*bond(1)%eq) then ! closest segment to ip1 is end-chain
               dx1 = vaux(1,ipnsegcn(iseg_min,ic_temp)) - vaux(1,ip3)
               dy1 = vaux(2,ipnsegcn(iseg_min,ic_temp)) - vaux(2,ip3)
               dz1 = vaux(3,ipnsegcn(iseg_min,ic_temp)) - vaux(3,ip3)
               r1 = sqrt(dx1**2+dy1**2+dz1**2)
               dx2 = vaux(1,ipnsegcn(iseg_min,ic_temp)) - vaux(1,ip2)
               dy2 = vaux(2,ipnsegcn(iseg_min,ic_temp)) - vaux(2,ip2)
               dz2 = vaux(3,ipnsegcn(iseg_min,ic_temp)) - vaux(3,ip2)
               r2 = sqrt(dx2**2+dy2**2+dz2**2)
               if (r1 < r2) then  ! end chain closer to the moving particle  ip3 ==> all particle of ic_temp to be moved
                  jp_min = ipnsegcn(1,ic_temp)
                  jp_max = ipnsegcn(npct(ict_temp),ic_temp)
                  iseg_ini = 1
                  iseg_fin = npct(ictpn(jp_max))
               else               ! end chain closer to the nonmoving particle  ip2 ==> all particle of ic_temp are not moved
                  jp_min = 0
                  jp_max = 0
                  iseg_ini = 0
                  iseg_fin = 0
                  cycle           ! go to the next ic_temp chain
               end if
            else
               if (jp1 > jp3) then ! jp_min and jp_max: first and last particle of the linear chain in pivot move
                  jp_min = ipnsegcn(1,ic_temp)
                  jp_max = jp1
                  iseg_ini = 1
                  iseg_fin = isegpn(jp1)
               else if (jp1 < jp3)  then
                  jp_min = jp1
                  jp_max = ipnsegcn(npct(ict_temp),ic_temp)
                  iseg_ini = isegpn(jp1)
                  iseg_fin = npct(ictpn(jp_max))
               end if
            end if

            do jp_loc = jp_min, jp_max
               nptm = nptm+1
               ipnptm(nptm) = jp_loc
               iptmpn(jp_loc) = nptm
               lptm(jp_loc) = .true.
            end do

         if (lclink) then
            nptmfirst = nptm + 1                                         ! used for the next generation
            do iseg = iseg_ini, iseg_fin                                 ! loop over main chain segments
               ip = ipnsegcn(iseg,ic_temp)                               ! particle number of main chain segment
               do icl = 1, nbondcl(ip)                                   ! loop over number of crosslinks
                  if (genic(icnpn(bondcl(icl,ip))) < genic(ic)) cycle    ! ignore crosslinks to lower generations.
                  do iseg3 = 1, sum(npptct(1:npt,ictpn(bondcl(icl,ip)))) ! loop over segments of relevant chain type
                     ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ip)))         ! particle number in side chain
                     nptm = nptm+1
                     ipnptm(nptm) = ip4
                     iptmpn(ip4) = nptm
                     lptm(ip4) = .true.
                  end do
               end do
            end do
            nptmlast = nptm                                              ! used for next generation
            do igen = genic(ic)+1, ngen
               do iploc = nptmfirst, nptmlast
                  do icl=1, nbondcl(ipnptm(iploc))
                     if (genic(icnpn(bondcl(icl,ipnptm(iploc)))) < igen) cycle      ! ignore crosslinks to lower generations.
                     do iseg3=1, sum(npptct(1:npt,ictpn(bondcl(icl,ipnptm(iploc)))))! loop over segments of relevant chain type
                        ip4 = ipnsegcn(iseg3,icnpn(bondcl(icl,ipnptm(iploc))))
                        nptm = nptm+1
                        ipnptm(nptm) = ip4
                        iptmpn(ip4) = nptm
                        lptm(ip4) = .true.
                     end do
                  end do
               end do
               nptmfirst = nptmlast+1      ! update for next generation
               nptmlast = nptm             ! update for next generation
            end do
         end if
      end do
   end if
   end if

   deallocate(dtemp)

end subroutine PivotDual

!........................................................................

end subroutine PivotMove

!************************************************************************
!*                                                                      *
!*     ChainMove                                                        *
!*                                                                      *
!************************************************************************

! ... perform one chain trial move

subroutine ChainMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iseg, ic, ict, ip, dnpcl
   real(8)    :: weight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (icnpn(ipmove) == 0) call Stop(txroutine, 'chain move based on a non-chain particle', uout)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ichainmove

   ic = icnpn(ipmove)     ! chain number ic where particle ipmove is residing
   ict = ictcn(ic)        ! chain type ict of chain ic

! .............. define particle(s) to be moved ..............

   nptm = 0                                 ! particles in chain ic
   do iseg = 1, npct(ict)
      ip = ipnsegcn(iseg,ic)
      nptm = nptm+1
      ipnptm(nptm) = ip
      iptmpn(ip) = nptm
      lptm(ip) =.true.
   end do

! ... and then the remaining particles

   if (lcl1chain(iptmove)) call ClusterMember('old', .false., .true., radcl1, pselectcl1)

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtranchain(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lpolyatom) call GetFixedTrialOri(nptm, ipnptm, ori, oritm)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

!  ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! ............. calculate nonenergetic weights .............

   weight = One
   dnpcl = Zero
   if (lcl1chain(iptmove)) then
      call ClusterMember('new', .false., .true., radcl1, pselectcl1)       !  calculate npclnew
      dnpcl = npclnew-npclold
      if (dnpcl/= 0) weight = (One-pselectcl1(iptmove))**dnpcl
   end if

!  ............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate       ! update energies and coordinates

end subroutine ChainMove

!************************************************************************
!*                                                                      *
!*     SlitherMove                                                      *
!*                                                                      *
!************************************************************************

! ... perform one chain slithering trial move
!     biased sampling of bond length and bond angle

subroutine SlitherMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SlitherMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iseg, ic, ict, ip, ip1, ip2, iploc, ih, igen, nbondloc, jp, jploc, jpt, jchain, jp_bonded, jseg
   integer(4) :: iseghead, idir, jc, jct, iplochead
   real(8)    :: rotaxis(3), theta, rotmat(3,3)
   real(8)    :: Random, dx, dy, dz, r1, r1i, r2, weight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (icnpn(ipmove) == 0) then
      call Stop(txroutine, 'slither move on non-chain particle', uout)
   else
      if (nh /= 0 .and. ngen /= 1) call Stop(txroutine, 'ngen /= 1', uout)
   end if

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = islithermove

   lptmdutwob = .true.   ! engage energy evaluation "inside" the moving group

! .............. define particle(s) to be moved ..............

   ic = icnpn(ipmove)                       ! chain number ic where particle ipmove is residing
   ict = ictcn(ic)                          ! chain type ict of chain ic
   ih = 0
   if (lhierarchical) ih = ihnpn(ipmove)

   nptm = 0                                 ! will be number of particles to be moved
   if (ih == 0) then                        ! nonhierarchical linear chain
      do iseg = 1, npct(ictpn(ipmove))
         ip = ipnsegcn(iseg,ic)
         nptm = nptm+1
         ipnptm(nptm) = ip
         iptmpn(ip) = nptm
         lptm(ip) =.true.
      end do
   else                                     ! hierarchical linear chain: slither along the backbone
      do igen = 0, ngen                                              ! loop over generations
         jct = ictgen(igen)                                          ! chain type
         do jc = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1 ! loop over chains of the structure
            do iseg = 1, npct(jct)                                   ! loop over segments
               ip = ipnsegcn(iseg,jc)
               nptm = nptm+1
               iptmpn(ip) = nptm
               ipnptm(nptm) = ip
               lptm(ip) = .true.
            end do
         end do
      end do
   end if

! .............. calculate a trial configuration ...............

   if (Random(iseed) < Half) then
      iseghead = 1
      iplochead = 1
      idir = +1
   else
      iseghead = npct(ict)
      iplochead = nptm
      idir = -1
   end if

! ... get the displacement in x- y- z- direction

   do iseg = 1, npct(ict)
      if (iseg == iseghead) then            ! head segment
         ip1 = ipnsegcn(iseg,ic)
         ip2 = ipnsegcn(iseg+idir,ic)
         call SlitherMoveSub                ! move head to a new position
      else                                  ! remaining segments
         ip1 = ipnsegcn(iseg-idir,ic)
         ip2 = ipnsegcn(iseg,ic)
         drotm(1,iseg) = ro(1,ip1)-ro(1,ip2)
         drotm(2,iseg) = ro(2,ip1)-ro(2,ip2)
         drotm(3,iseg) = ro(3,ip1)-ro(3,ip2)
      end if
   end do

! ... add the displacement to ro and put in rotm

   do iseg = 1, npct(ict)
      ip = ipnsegcn(iseg,ic)
      call PBC(drotm(1,iseg),drotm(2,iseg),drotm(3,iseg))
      rotm(1,iseg) = ro(1,ip)+drotm(1,iseg)
      rotm(2,iseg) = ro(2,ip)+drotm(2,iseg)
      rotm(3,iseg) = ro(3,ip)+drotm(3,iseg)
      if (lclink) then
         if (nbondcl(ip) /= 0) then                    ! check for side chains
            do nbondloc = 1, nbondcl(ip)
               jp = bondcl(nbondloc,ip)                ! particle bonded to ip
               jpt = ictpn(jp)                         ! chain type of jp particle
               jchain = icnpn(jp)                      ! chain of jp particle
               do jseg = 1, npct(jpt)
                  jp_bonded = ipnsegcn(jseg,jchain)    ! particles belonging to the chain jchain bound to particle ip
                  jploc = iptmpn(jp_bonded)
                  rotm(1,jploc) = ro(1,jp_bonded)+drotm(1,iseg)
                  rotm(2,jploc) = ro(2,jp_bonded)+drotm(2,iseg)
                  rotm(3,jploc) = ro(3,jp_bonded)+drotm(3,iseg)
               end do
            end do
         end if
      end if
   end do

   call CheckPartBCTM(nptm, rotm, lboxoverlap)

   if (lpolyatom) then
         do iploc = iplochead, nptm+1-iplochead,  idir
            ip = ipnptm(iploc)
            if (iploc == iplochead) then
               oritm(1:3,1:3,iploc) = ori(1:3,1:3,ip)
            else
               oritm(1:3,1:3,iploc) = ori(1:3,1:3,ip-idir)
            end if
         end do
   end if

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

!  ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! ............. calculate nonenergetic weights .............

   weight = One      ! note, the biased sampling is not necessarily correct for hetrochains

!  ............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, (du%tot-(du%bond+du%angle))*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate     ! update energies and coordinates

   lptmdutwob =.false.                        ! restore lptmdutwob

contains

!........................................................................

subroutine SlitherMoveSub

   real(8)    :: dot

   dx = (ro(1,ip1)-ro(1,ip2))
   dy = (ro(2,ip1)-ro(2,ip2))
   dz = (ro(3,ip1)-ro(3,ip2))
   call PBCr2(dx,dy,dz,r2)
   r1i = One/sqrt(r2)
   dx = dx*r1i
   dy = dy*r1i
   dz = dz*r1i

   rotaxis(1) = Random(iseed)-Half             ! rotation axis rotaxis perpendicular to (dx,dy,dz)
   rotaxis(2) = Random(iseed)-Half
   rotaxis(3) = Random(iseed)-Half
   dot = rotaxis(1)*dx+rotaxis(2)*dy+rotaxis(3)*dz
   rotaxis(1) = rotaxis(1)-dot*dx
   rotaxis(2) = rotaxis(2)-dot*dy
   rotaxis(3) = rotaxis(3)-dot*dz

   call BondAngleTab('value',ict,theta)        ! get bond angle between three consequtive beads
   theta = Pi-theta

   call AxisAngToOri(rotaxis,theta,rotmat)     ! get rotation matrix rotmat

   call BondLengthTab('value',ict,r1)          ! get bond length r1

   drotm(1,iseg) = r1*(rotmat(1,1)*dx+rotmat(1,2)*dy+rotmat(1,3)*dz)
   drotm(2,iseg) = r1*(rotmat(2,1)*dx+rotmat(2,2)*dy+rotmat(2,3)*dz)
   drotm(3,iseg) = r1*(rotmat(3,1)*dx+rotmat(3,2)*dy+rotmat(3,3)*dz)

end subroutine SlitherMoveSub

!........................................................................

end subroutine SlitherMove

!************************************************************************
!*                                                                      *
!*     BrushMove                                                        *
!*                                                                      *
!************************************************************************

! ... perform one brush trial move

subroutine BrushMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='BrushMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc, ic, iseg, ip, dnpcl
   real(8)    :: weight, MCWeight, UmbrellaWeight, MCPmfWeight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ibrushmove

! .............. define particle(s) to be moved ..............

! ... central particle

   nptm = 1
   iploc = 1
   ipnptm(iploc) = ipmove
   lptm(ipmove) = .true.

! ... particles in the brush

   do ic = ipmove, nc, nppt(iptpn(ipmove))     ! chains grafted to particle ipmove
      do iseg = 1, npct(ictcn(ic))             ! particles of chain ic
         nptm = nptm+1
         ip = ipnsegcn(iseg,ic)
         ipnptm(nptm) = ip
         lptm(ip) =.true.
      end do
   end do

! ... and then the remaining particles

   if (lcl1brush(iptmove)) call ClusterMember('old', .true., .false., radcl1, pselectcl1)

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtranbrush(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   call RotTranBrush
   call CheckPartBCTM(nptm, rotm, lboxoverlap)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! .............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! .............. calculate nonenergetic weights ...............

   weight = One
   dnpcl = Zero
   if (lcl1brush(iptmove)) then
      if (lvlist) call ClusterMember('new', .true., .false., radcl1, pselectcl1)       ! calculate npclnew
      dnpcl = npclnew-npclold
      if (dnpcl /= 0) weight = weight*(One-pselectcl1(iptmove))**dnpcl
   end if

   if (lmcweight) weight = weight*MCWeight()
   if (lautumb) weight = weight*UmbrellaWeight(1)
   if (lmcpmf) weight = weight*MCPmfWeight(1)

! .............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate        ! update energies and coordinates

   if (lshiftzcom(iptmove)) call ShiftZDirection

   if (lautumb) call UmbrellaUpdate              ! update weight function for umbrella potential
   if (lmcpmf) call MCPmfUpdate                  ! update weight function for mc pmf

contains

!........................................................................

subroutine RotTranBrush

   real(8) :: dx, dy, dz, dxpbc, dypbc, dzpbc, dxtran, dytran, dztran, alpha, rotaxis(3)
   integer(4) :: ic, ip
   real(8) :: Random

! ... get rotation displacement

   call SphRandom(iseed, rotaxis(1), rotaxis(2), rotaxis(3))
   alpha = drotbrush(iptmove)*(Random(iseed)-0.5)

! ... get translational displacement

   dxtran = rotm(1,iploc) - ro(1,ipmove)  ! x trial-displacement of partice ipmove
   dytran = rotm(2,iploc) - ro(2,ipmove)  ! y trial-displacement of partice ipmove
   dztran = rotm(3,iploc) - ro(3,ipmove)  ! z trial-displacement of partice ipmove
   if (lbccyl) then                       ! move only in the z-direction (for pmf calculations)
      dxtran = Zero
      dytran = Zero
   end if

! ... undo pbc for particles in the brush

   vaux(1:3,ipmove) = ro(1:3,ipmove)
   do ic = ipmove, nc, nppt(iptpn(ipmove))
     call UndoPBCChain(ro(1:3,ipmove), ic, 1, vaux)
   end do

! ... and remaining particles (if cluster1 move)

   if (lcl1brush(iptmove)) then
     do iploc = nprimpartcl+1, nptm
        ip = ipnptm(iploc)
        vaux(1,ip) = ro(1,ip)
        vaux(2,ip) = ro(2,ip)
        vaux(3,ip) = ro(3,ip)
        dx = vaux(1,ip)-ro(1,ipmove)
        dy = vaux(2,ip)-ro(2,ipmove)
        dz = vaux(3,ip)-ro(3,ipmove)
        call PBC2(dx,dy,dz,dxpbc,dypbc,dzpbc)
        vaux(1,ip) = vaux(1,ip) - dxpbc
        vaux(2,ip) = vaux(2,ip) - dypbc
        vaux(3,ip) = vaux(3,ip) - dzpbc
     end do
  end if

! ... rotate

   call RotateSetPart(rotaxis, alpha, ro(1:3,ipmove), nptm, ipnptm, vaux, rotm, drotm)

! ... translate

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      rotm(1,iploc) = rotm(1,iploc)+dxtran    ! add contribution from translation
      rotm(2,iploc) = rotm(2,iploc)+dytran
      rotm(3,iploc) = rotm(3,iploc)+dztran
      drotm(1,iploc) = drotm(1,iploc)+dxtran
      drotm(2,iploc) = drotm(2,iploc)+dytran
      drotm(3,iploc) = drotm(3,iploc)+dztran
   end do

end subroutine RotTranBrush

!........................................................................

end subroutine BrushMove

!************************************************************************
!*                                                                      *
!*     BrushCl2Move                                                     *
!*                                                                      *
!************************************************************************

! ... perform one brush + cluster2 trial move

subroutine BrushCl2Move(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='BrushiCl2Move'
   integer(4) :: mnpair = 10000
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc, iplow, ipupp, ic, iseg, jp, iprimpart, iprimpartloc
   integer(4) :: npair, npcl
   integer(4), allocatable :: n1(:), n2(:), iobjcluster(:), icllis(:)
   real(8)    :: weight, r2, dx, dy, dz

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ibrushcl2move

   iplow = ipnpt(iptmove)
   ipupp = ipnpt(iptmove)+nppt(iptmove)-1

! .............. define particle(s) to be moved ..............

! ... consider first primary cluster particles (particles of type iptmove including particle ipmove)

   if (.not.allocated(n1)) then
      allocate(n1(mnpair), n2(mnpair), iobjcluster(mnpair), icllis(mnpair))
      n1 = 0
      n2 = 0
      iobjcluster = 0
      icllis = 0
   end if

   call CalcPartPairListIpt(iptmove,radcl2(iptmove)**2,mnpair,npair,n1,n2) ! list of particle pair for particles of type iptmove
   call MakeCluster(nppt(iptmove), npair, n1, n2, iobjcluster)
   call CalcClusterMember(nppt(iptmove), iobjcluster, ipmove-ipnpt(iptmove)+1, npcl, icllis)

   nptm = 0
   do iploc = 1, npcl
      nptm = nptm + 1
      ipnptm(nptm) = icllis(iploc) + (ipnpt(iptmove) - 1)
      lptm(ipnptm(nptm)) = .true.
   end do
   nprimpartcl = nptm

   deallocate(n1, n2, iobjcluster, icllis)

! ... its brush (here grafted chains cyclically placed on the particles)

   do iprimpartloc = 1, nprimpartcl
      iprimpart = ipnptm(iprimpartloc)
      do ic = iprimpart, nc, nppt(iptpn(iprimpart))
         do iseg = 1, npct(ictcn(ic))
            nptm = nptm+1
            iploc = iploc+1
            jp = ipnsegcn(iseg,ic)
            ipnptm(nptm) = jp
            lptm(jp) =.true.
         end do
      end do
   end do

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtranbrush(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lpolyatom) call GetRandomTrialOri(drot(iptmove), iseed, ori(1,1,ipmove), oritm(1,1,iploc))

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! .............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)
   if (lhsoverlap .or. lhepoverlap) goto 200

! ... discard moves where nprimpartcl would change

! write(*,*) 'check if new cluster2 members'
   do iploc = 1, npcl
!     write(*,*) 'iploc = ', iploc, 'ip', ipnptm(iploc)
      do jp = iplow, ipupp
!        if (lptm(jp)) write(*,*) jp, 'is already a member'
         if (lptm(jp)) cycle                  ! consider only particles not already in the cluster
         dx = ro(1,jp)-rotm(1,iploc)
         dy = ro(2,jp)-rotm(2,iploc)
         dz = ro(3,jp)-rotm(3,iploc)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < radcl2(iptmove)**2) then
            ievent = imcreject
!           write(*,*) jp , 'is too close', sqrt(r2)
            goto 300
         end if
      end do
   end do

! ............. calculate nonenergetic weights .............

   weight = One

! .............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

300 continue

! .............. update .............

   if (ievent == imcaccept) call MCUpdate        ! update energies and coordinates

end subroutine BrushCl2Move

!************************************************************************
!*                                                                      *
!*     HierarchicalMove                                                 *
!*                                                                      *
!************************************************************************

! ... perform one hierarchical trial move

subroutine HierarchicalMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='HierarchicalMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: igen, iseg, ic, ict, ip, ih
   real(8)    :: weight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

!xx  if (icnpn(ipmove) == 0) call Stop(txroutine, 'hierarchical move based on a non-chain particle', uout)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ihierarchicalmove

! .............. define particle(s) to be moved ..............

   ih = ihnpn(ipmove)                                                      ! hierarchical strucure
   nptm = 0
   do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
            nptm = nptm+1
            iptmpn(ip) = nptm
            ipnptm(nptm) = ip
            lptm(ip) =.true.
         end do
      end do
   end do

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtranhierarchical(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   call CheckPartBCTM(nptm, rotm, lboxoverlap)
   if (lpolyatom) call GetFixedTrialOri(nptm, ipnptm, ori, oritm)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

!  ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! ............. calculate nonenergetic weights .............

   weight = One

!  ............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate       ! update energies and coordinates

end subroutine HierarchicalMove

!************************************************************************
!*                                                                      *
!*     NetworkMove                                                      *
!*                                                                      *
!************************************************************************

! ... perform one network trial move

subroutine NetworkMove(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NetworkMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc
   real(8)    :: weight

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = inetworkmove

   lptmdutwob = .true.   ! engage energy evaluation "inside" the moving group

! .............. define particle(s) to be moved ..............

! ... central particle

   nptm = 1
   iploc = 1
   ipnptm(iploc) = ipmove
   lptm(ipmove) = .true.

! ... and then the remaining particles

!   if (lcl1network(iptmove)) call ClusterMember('old', .true., .false., radcl1, pselectcl1)

! .............. calculate a trial configuration ...............

   call GetRandomTrialPos(dtrannetwork(iptmove), iseed, nptm, ipnptm, ro, rotm, drotm)
   call TranNetwork
   call CheckPartBCTM(nptm, rotm, lboxoverlap)

!   write(*,*) 'MC: network: lfixedori:',lfixedori
   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   call SetTrialAtomProp
   if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)

   if (lweakcharge) then
      laztm(1:natm) = laz(ianatm(1:natm))
      if (lewald) aztm(1:natm) = az(ianatm(1:natm))
   end if

   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! .............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! .............. calculate nonenergetic weights ...............

   weight = One

! .............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) call MCUpdate        ! update energies and coordinates

   lptmdutwob =.false.                           ! restore lptmdutwob

contains

!........................................................................

subroutine TranNetwork

   real(8) :: dx, dy, dz, dxtran, dytran, dztran, fac, term
   integer(4) :: ic, ip, ip1, ip2, icl, nseg

! ... get translational displacement of node

   dxtran = rotm(1,iploc) - ro(1,ipmove)  ! x trial-displacement of partice ipmove
   dytran = rotm(2,iploc) - ro(2,ipmove)  ! y trial-displacement of partice ipmove
   dztran = rotm(3,iploc) - ro(3,ipmove)  ! z trial-displacement of partice ipmove

! ... identify network particles to be moved and move them

   do icl = 1, nbondcl(ipmove)            ! loop over crosslinks to particle ipmove
      ip1 = bondcl(icl,ipmove)            ! id of crosslink end
      ic = icnpn(ip1)                     ! id chain containing the crosslink end
      nseg = npct(ictcn(ic))              ! number of segments of the chain
      ip2 = ip1 + nseg-1                  ! proposed other end of the chain
      if(ip2 > np) then
         ip2 = ip1 - (nseg-1)                        ! select the other posibility
      else
         if (icnpn(ip2) /= ic) ip2 = ip1 - (nseg-1)  ! select the other posibility
      end if
      fac = One
      term = -One/nseg

      do ip = ip1, ip2 - sign(1, int(ip2-ip1)), sign(1, int(ip2-ip1))
         fac = fac + term
         nptm = nptm + 1
         ipnptm(nptm) = ip
         iptmpn(ip) = nptm
         lptm(ip) = .true.
         dx = fac*dxtran
         dy = fac*dytran
         dz = fac*dztran

         rotm(1,nptm) = ro(1,ip)+dx
         rotm(2,nptm) = ro(2,ip)+dy
         rotm(3,nptm) = ro(3,ip)+dz
         drotm(1,nptm) = dx
         drotm(2,nptm) = dy
         drotm(3,nptm) = dz
      end do
   end do

   end subroutine TranNetwork

end subroutine NetworkMove

!************************************************************************
!*                                                                      *
!*     VolChange                                                        *
!*                                                                      *
!************************************************************************

! ... perform one volume change trial move

subroutine VolChange(iStage)

   use MCModule
   use CellListModule, only: InitCellList, SetCellList
   use NListModule, only: drnlist
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='VolChange'

   type(potenergy_var) :: uold
   real(8)    :: volold, boxlenratio, dpv, dlnvol, dvaux, Random
   integer(4) :: i

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ivolumemove

! .............. calculate a trial configuration ...............

! ... save old configurational energy and boxlen lengths

   uold = u
   volold = vol

! ... calculate new box lengths and coordinates

   boxlenratio = ((volold+dvol*(Random(iseed)-Half))/volold)**(Third)
   boxlen(1:3) = boxlen(1:3)*boxlenratio
   call SetBoxParam
   if (lewald) call EwaldSetup
   do i = 1, np
      ro(1,i) = ro(1,i)*boxlenratio
      ro(2,i) = ro(2,i)*boxlenratio
      ro(3,i) = ro(3,i)*boxlenratio
   end do
   call SetAtomPos(1, np, lintsite)
   if (lCellList) then
      call InitCellList(rcut + drnlist, iStage)
      call SetCellList()
   end if

!  ............. evaluate energy difference ...............

   call UTotal(iStage)
   du%tot = u%tot-uold%tot
   if (lintsite) call SetAtomPos(1, np, .false.)

! ............... get "enthalpy" difference .............

   dpv = prsrst*sclpre*(vol-volold)*sclvol*AvNo/sclene
   dlnvol =-np/beta*log(vol/volold)
   dvaux = du%tot+dpv+dlnvol

   if (itestmc == 3) call TestVolChange(uout)

!  ............. decide new configuration .............

   call Metropolis(.false., .false., .false., One, dvaux*beta)

! .............. update .............

   if (ievent == imcaccept) then       ! accepted volume change

   else                                ! rejected volume change

      u = uold
      boxlenratio = One/boxlenratio
      boxlen = boxlen*boxlenratio
      call SetBoxParam
      if (lewald) call EwaldSetup
      ro = ro*boxlenratio
      call SetAtomPos(1, np, .false.)
      if (lCellList) then
         call InitCellList(rcut + drnlist, iStage)
         call SetCellList()
      end if

   end if

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

!........................................................................

contains

subroutine TestVolChange(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,t20,3f15.5)') 'uold%tot, u%tot     ', uold%tot, u%tot
   write(unit,'(a,t20,3f15.5)') 'du%tot              ', du%tot
   write(unit,'(a,t20,3f15.5)') 'pressure            ', prsr
   write(unit,'(a,t20,3f15.5)') 'boxlenratio         ', boxlenratio
   write(unit,'(a,t20,3f15.5)') 'volold, vol         ', volold, vol
   write(unit,'(a,t20,3f15.5)') 'd(pv)               ', dpv
   write(unit,'(a,t20,3f15.5)') 'd(lnvol)            ', dlnvol
   write(unit,'(a,t20,3f15.5)') 'dvaux               ', dvaux
end subroutine TestVolChange

!........................................................................

end subroutine VolChange

!************************************************************************
!*                                                                      *
!*     NPartChange                                                      *
!*                                                                      *
!************************************************************************

! ... perform one number of particle change trial move

!     note, not check to work with weak charges

subroutine NPartChange(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NPartChange'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap, lDelete
   integer(4) :: jp, jpt, iptjpt
   real(8) :: prandom, uuu, fforce(3), weight, Random

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = inpartmove

! ... deletion or insertion?

   lDelete = .true.
   prandom = Random(iseed)
   if (prandom > half) lDelete = .false.

   du%twob = Zero
   lboxoverlap = .false.
   lhsoverlap = .false.
   lhepoverlap = .false.

! ... calculate energy difference

   if (lDelete) then    ! trial delete of particle ipmove

      do jp = 1, np
         jpt = iptpn(jp)
         iptjpt = iptpt(iptmove, jpt)
         if (jp == ipmove) cycle
         call UTwoBodyPair(ipmove, jp, uuu, fforce)
         du%twob(iptjpt) = du%twob(iptjpt) - uuu    ! note, minus sign
      end do
      du%twob(0) = sum(du%twob(1:nptpt))
      weight = nppt(iptmove)/vol * exp(-beta*chempot(iptmove))

   else                 ! trial insert a particle of type iptmove

      nppt(iptmove) = nppt(iptmove) + 1            ! add one particle of type iptmove

      if (sum(nppt(1:npt)) > np_alloc) then
          write(*,*) 'sum(nppt), np_alloc', sum(nppt(1:npt)), np_alloc
          call stop(txroutine, 'sum(nppt(1:npt)) > np_alloc', uout)
      end if

      call SetObjectParam1
      call SetObjectParam2
      ipmove = np                                  ! place the new particle last
      continue
      call SetPartPosRandomMC(ipmove)              ! set trial random particle position
!      write(*,'(a,3f10.3)') 'trial position', ro(1:3,ipmove)
!      if (abs(ro(3,ipmove)) > 0.8*boxlen2(3)) then
!         write(*,*) 'too close a surface'
!         goto 10
!      end if
      call SetPartOriRandom(iseed,ori(1,1,ipmove)) ! set trial random particle orientation
      call SetAtomProp(ipmove,ipmove,.false.)      ! trial atom properties

      do jp = 1, np-1
         jpt = iptpn(jp)
         iptjpt = iptpt(iptmove, jpt)
         call UTwoBodyPair(ipmove, jp, uuu, fforce)
         if (uuu >= 1d10) then
            lhsoverlap = .true.
            goto 400
         end if
         du%twob(iptjpt) = du%twob(iptjpt) + uuu    ! note, plus sign
      end do
      du%twob(0) = sum(du%twob(1:nptpt))
      weight = vol/nppt(iptmove)*exp(beta*chempot(iptmove))
400   continue

   end if

   du%tot = du%twob(0)                            ! could be generalized

   if (itestmc == 3) call TestNPartChange1(uout)

!  ............. decide new configuration .............

   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) then        ! accepted change of number of partciles

      if (lDelete) then
         nppt(iptmove) = nppt(iptmove) - 1
         ro(1:3,ipmove) = ro(1:3,np)
         ori(1:3,1:3,ipmove) = ori(1:3,1:3,np)
         call SetObjectParam1
         call SetObjectParam2
         call SetAtomProp(ipmove,ipmove,.false.)      ! trial atom properties
      else
      end if
      call SetVList                                   ! update neighbour list
      u%tot           = u%tot           + du%tot
      u%twob(0:nptpt) = u%twob(0:nptpt) + du%twob(0:nptpt)

   else                                ! rejected change of number of paricles

      if (lDelete) then
      else
         nppt(iptmove) = nppt(iptmove) - 1
         call SetObjectParam1
         call SetObjectParam2
      end if

   end if
   if (itestmc == 3) call TestNPartChange2(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine TestNPartChange1(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1', unit)
   write(unit,'(a,t20,i4    )') 'ipmove              ', ipmove
   write(unit,'(a,t20,i4    )') 'iptmove             ', iptmove
   write(unit,'(a,t20,l     )') 'lDelete             ', lDelete
   write(unit,'(a,t20,l     )') 'lhsoverlap          ', lhsoverlap
   write(unit,'(a,t20,3f15.5)') 'volume              ', vol
   write(unit,'(a,t20,i4    )') 'nppt(iptmove)       ', nppt(iptmove)
   write(unit,'(a,t20,3f15.5)') 'weight              ', weight
   write(unit,'(a,t20,3f15.5)') 'du%tot              ', du%tot
   write(unit,'(a,t20,3f15.5)') 'du%twob(0:npt)     ', du%twob(0:npt)
   if (.not. lDelete) then
   write(unit,'(a,t20,3f15.5)') 'trial coordinates   ', ro(1:3,ipmove)
   end if
end subroutine TestNPartChange1

!........................................................................

subroutine TestNPartChange2(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a,t20,i4    )') 'istep1              ', istep1
   write(unit,'(a,t20,i4    )') 'istep2              ', istep2
   write(unit,'(a,t20,i4    )') 'ipass               ', ipass
   write(unit,'(a,t20,i4    )') 'ievent              ', ievent
   write(unit,'(a,t20,i4    )') 'np                  ', np
   write(unit,'(a,t20,i4    )') 'nppt(iptmove)       ', nppt(iptmove)
   call WritePartAtomProp
end subroutine TestNPartChange2

!........................................................................

end subroutine NPartChange

!************************************************************************
!*                                                                      *
!*     ChargeChange                                                     *
!*                                                                      *
!************************************************************************

! ... perform a charge-change trial move

!  zat = -1:   HA   <->   A(-) + H(+);        HA weak acid and A(-) weak base
!  zat = +1:   BH(+)   <->   B + H(+);        BH(+) weak acid and B weak base

!  pH-pKa < 0 (low  pH) => reaction shifted left
!  pH-pKa > 0 (high pH) => reaction shifted right

!  devloped for jatweakcharge = 0:  wihtout counterions
!               jatweakcharge > 0:  atom type of counter ions (atoms)

subroutine ChargeChange(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChargeChange'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   real(8) ::  random,  weight, xsign
   integer(4) :: ip, ipt, ia, ialoc, ianmove, iatmove

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ichargemove

   lptmdutwob = .true.   ! engage energy evaluation "inside" the moving group

! .............. define charge(s) to be changed ..............

   nptm = 1
   ipnptm(1) = ipmove
   lptm(ipmove)  =.true.
   if (lmonoatom) then                 ! trial charge change of the single atom
      ia = ipmove
   else                                ! select a single atom for trail charge change ; e.g. lysozyme:
      ipt = iptpn(ipmove)
      do                                ! select one eligible atom of the particle
         ia = ianpn(ipmove) + int(random(iseed)*napt(ipt)) ! select one atom of the particle
         if (latweakcharge(iatan(ia))) exit                ! check if eligible for trial charge change
      end do
   end if
   ianatm(1) = ia               ! update ianatm
   if (jatweakcharge(iatan(ia)) > 0) then  ! get id of atom and particle carrying counter charge
      nptm = 2
      ia = iananweakcharge(ianatm(1))  ! id of atom
      ianatm(2) = ia
      ip = ipnan(ia)                   ! id of particle
      ipnptm(2) = ip
      lptm(ip) = .true.
   end if
   natm = nptm

   rotm(1:3,1:nptm)      = ro(1:3,ipnptm(1:nptm))
   oritm(1:3,1:3,1:nptm) = ori(1:3,1:3,ipnptm(1:nptm))
   rtm(1:3,1:natm) = r(1:3,ianatm(1:natm))
   call CheckPartBCTM(nptm, rotm, lboxoverlap)

   if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
       call AddNeighbours
       call UpdateOri
   end if

   laztm(1:natm) =.not.laz(ianatm(1:natm))  ! get trial charge states

   if (lewald) then
      do ialoc = 1, natm
         ia = ianatm(ialoc)
         if (laztm(ialoc)) then
            aztm(ialoc) = zat(iatan(ia))
         else
            aztm(ialoc) = Zero
         end if
      end do
   end if

if (itest == 90) then
   call writehead(3, txroutine, uout)
   write(uout,'(a,2i5)') 'ipnapm(iploc)', (ipnptm(1:nptm))
   write(uout,'(a,2i5)') 'ianatm(ialoc)', (ianatm(1:natm))
   write(uout,'(a,2l5)') 'laz(ialoc)', (laz(ianatm(1:natm)))
   write(uout,'(a,2l5)') 'laztm(ialoc)', (laztm(1:natm))
end if


!   call SetTrialAtomProp

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)

! ............. calculate nonenergetic weights .............

   ianmove = ianatm(1)                                ! id of titrating atom
   iatmove = iatan(ianmove)                           ! atome type of titrating atom
   xsign = sign(one,zat(iatmove))                     ! sign of charge of the charged state
   if (.not.laz(ianmove)) xsign = -xsign              ! sign to be used in weight
   weight = exp(xsign*log(10.0d0)*pHmpK(iatmove))

 !  write(*,*) 'iptmove', iptmove
 !  write(*,*) 'iatmove', iatmove
 !  write(*,*) 'ipnptm', ipnptm(1:nptm)
 !  write(*,*) 'ianmove', ianmove
 !  write(*,*) 'weight', weight

!  ............. decide new configuration .............

200 continue
   if (itestmc == 3) call TestChargeChange1(uout)
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) then
      call MCUpdate       ! update energies and coordinates
      do ialoc = 1, natm
         laz(ianatm(ialoc)) = laztm(ialoc)    ! update charge status
         if (lewald) az(ianatm(ialoc)) = aztm(ialoc) ! update atom charge
      end do
   end if

if (itest == 90) then
   call writehead(3, txroutine, uout)
   write(uout,'(a,e12.5)')  'weight', weight                           !cc
   write(uout,'(a,i5)')     'ievent', ievent                          !cc
   write(uout,'(a,10l5)')   'laztm', laztm(1:2)                          !cc
   write(uout,'(a,10l5)')   'laz', laz(ianatm(1:2))                          !cc
   write(uout,'(a,10f12.5)')'u%twob(0:nptpt)', real(u%twob(0:nptpt))  !cc
end if

   lptmdutwob =.false.                        ! restore lptmdutwob

   if (itestmc == 3) call TestChargeChange2(uout)

contains

!........................................................................

subroutine TestChargeChange1(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1', unit)
   write(unit,'(a,t20,i4    )') 'ipmove              ', ipmove
   write(unit,'(a,t20,i4    )') 'iptmove             ', iptmove
   write(unit,'(a,t20,2l    )') 'laz(ipnptm(1:nptm)) ', laz(ipnptm(1:nptm))
   write(unit,'(a,t20,2l    )') 'laztm(1:nptm)       ', laztm(1:nptm)
   write(unit,'(a,t20,3f15.5)') 'weight              ', weight
   write(unit,'(a,t20, l    )') 'lboxoverlap         ' ,lboxoverlap
   write(unit,'(a,t20, l    )') 'lhsoverlap          ', lhsoverlap
   write(unit,'(a,t20, l    )') 'lhepoverlap         ', lhepoverlap
   write(unit,'(a,t20,1f15.5)') 'du%tot              ', du%tot
   write(unit,'(a,t20,7f15.5)') 'du%twob(0:nptpt)    ', du%twob(0:nptpt)
   write(unit,'(a,t20,1f15.5)') 'u%tot       (old)   ', u%tot
   write(unit,'(a,t20,100l  )') 'laz(1:na) (old)     ' ,laz(1:na)
end subroutine TestChargeChange1

!........................................................................

subroutine TestChargeChange2(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a,t20, i5   )') 'ievent              ', ievent
   write(unit,'(a,t20,1f15.5)') 'u%tot     (new)     ', u%tot
   write(unit,'(a,t20,100l  )') 'laz(1:na) (new)     ', laz(1:na)
end subroutine TestChargeChange2

!........................................................................

end subroutine ChargeChange

!************************************************************************
!*                                                                      *
!*     SetPartPosRandomMC                                               *
!*                                                                      *
!************************************************************************

! ... generate a random position

subroutine SetPartPosRandomMC(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip           ! paricle

   real(8), parameter :: fac = One - 0.1D-1
   real(8) :: Random, dxpbc, dypbc, dzpbc

   if (lbcbox) then
      ro(1,ip) = -boxlen2(1) + boxlen(1)*Random(iseed)
      ro(2,ip) = -boxlen2(2) + boxlen(2)*Random(iseed)
      ro(3,ip) = -boxlen2(3) + boxlen(3)*Random(iseed)
   else if (lbcrd) then
      do
         ro(1,ip) = boxlen(1)*(Random(iseed) - Half)
         ro(2,ip) = boxlen(2)*(Random(iseed) - Half)
         ro(3,ip) = boxlen(3)*(Random(iseed) - Half)
         call PBC2(ro(1,ip),ro(2,ip),ro(3,ip),dxpbc,dypbc,dzpbc)
         if ((dxpbc == Zero) .and. (dypbc == Zero) .and. (dzpbc == Zero)) exit
      end do
   else if (lbcto) then
      do
         ro(1,ip) = boxlen(1)*(Random(iseed) - Half)
         ro(2,ip) = boxlen(2)*(Random(iseed) - Half)
         ro(3,ip) = boxlen(3)*(Random(iseed) - Half)
         call PBC2(ro(1,ip),ro(2,ip),ro(3,ip),dxpbc,dypbc,dzpbc)
         if ((dxpbc == Zero) .and. (dypbc == Zero) .and. (dzpbc == Zero)) exit
      end do
   else if (lbcsph) then
      do
         ro(1,ip) = Two*sphrad*(Random(iseed)-0.5)
         ro(2,ip) = Two*sphrad*(Random(iseed)-0.5)
         ro(3,ip) = Two*sphrad*(Random(iseed)-0.5)
         if (sum(ro(1:3,ip)**2) <= fac*sphrad2) exit
      end do
   else if (lbccyl) then
      do
         ro(1,ip) = Two*cylrad*(Random(iseed)-0.5)
         ro(2,ip) = Two*cylrad*(Random(iseed)-0.5)
         if (sum(ro(1:2,ip)**2) <= fac*cylrad2) exit
      end do
      ro(3,ip) = -boxlen2(3) + boxlen(3)*Random(iseed)
   else if (lbcell) then
      do
         ro(1,ip) = Two*ellrad(1)*(Random(iseed)-0.5)
         ro(2,ip) = Two*ellrad(2)*(Random(iseed)-0.5)
         ro(3,ip) = Two*ellrad(2)*(Random(iseed)-0.5)
         if (sum((ro(1:3,ip)*ellradi(1:3))**2) <= fac) exit
      end do
   end if

end subroutine SetPartPosRandomMC

!************************************************************************
!*                                                                      *
!*     MCAllDriver                                                      *
!*                                                                      *
!************************************************************************

! ... mcall driver

subroutine MCAllDriver(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCAllDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IOMCAll(iStage)
      if (lautumb)  call UmbrellaIO(iStage)

   case (iWriteInput)

      call IOMCAll(iStage)
      if (lautumb)  call UmbrellaIO(iStage)

   case (iBeforeSimulation)

      call MCAllAver(iStage)
      if (lautumb)  call UmbrellaIO(iStage)

   case (iBeforeMacrostep)

      call MCAllAver(iStage)

   case (iSimulationStep)

      call MCAllPass(iStage)
      call MCAllAver(iStage)

   case (iAfterMacrostep)

      call MCAllAver(iStage)

   case (iAfterSimulation)

      call MCAllAver(iStage)
      if (lautumb)  call UmbrellaIO(iStage)

   end select

end subroutine MCAllDriver

!************************************************************************
!*                                                                      *
!*     IOMCAll                                                          *
!*                                                                      *
!************************************************************************

! ... performing i/o on monte carlo all variables

subroutine IOMCAll(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOMCAll'

   namelist /nmlMCAll/ dtranall, drotall, lautumb

   select case (iStage)
   case (iReadInput)

      if (.not. allocated(dtranall)) then
         allocate(dtranall(npt))
         dtranall = 0.0E+00
      end if
      if (.not. allocated(drotall)) then
         allocate(drotall(npt))
         drotall = 0.0E+00
      end if

      lautumb = .false.

      rewind(uin)
      read(uin, nmlMCAll)

   case (iWriteInput)

      if (master) then
         call WriteHead(2, 'mcall data', uout)
         if (dtranall(1) < 0)  write(uout,'(a)') 'spherical displacement area'
         if (dtranall(1) > 0)  write(uout,'(a)') 'cubic displacement volume'
         write(uout,'(a,t35,4(f10.3,2x))') 'trans. displacement param.  = ', dtranall(1:npt)
         write(uout,'(a,t35,4(f10.3,2x))') 'rot. displacement param.    = ', drotall(1:npt)
         if (lautumb)          write(uout,'(a)') 'generate umbrella potential with automatic update procedure'

         drotall(1:npt)=drotall(1:npt)*sclang
      end if
#if defined (_PAR_)
      if (lautumb) call Stop(txroutine, 'generation of umbreall potential not adapted for _PAR_', uout)
#endif

   end select

end subroutine IOMCAll

!************************************************************************
!*                                                                      *
!*     MCAllPass                                                        *
!*                                                                      *
!************************************************************************

! ... perform one mc pass (simultaneous move of all particles)

!     one configuration involves

!             a) get a trial configuration by translating and rotation all particles
!             b) evaluate energy difference
!             c) decide new configuration

subroutine MCAllPass(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCallPass'
   logical    :: lboxoverlap, lmetropolis
   !logical    :: lhsoverlap, lhepoverlap
   integer(4) :: ip, ipt
   real(8)    :: weight, UmbrellaWeight,  Random
   type(potenergy_var) :: usave
   real(8)             :: prsrsave
   real(8), allocatable, save :: rosave(:,:)    ! temporary use (need not to be saved)
   real(8), allocatable, save :: orisave(:,:,:) ! temporary use (need not to be saved)
   real(8), allocatable, save :: boxlensave(:)  ! temporary use (need not to be saved)
   real(8), allocatable, save :: idmsyssave(:)  ! temporary use (need not to be saved)
   real(8), allocatable, save :: idmosave(:,:)  ! temporary use (need not to be saved)
   real(8), allocatable, save :: idmsave(:,:)   ! temporary use (need not to be saved)
   real(8), allocatable, save :: forcesave(:,:) ! temporary use (need not to be saved)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   if (ltime) call CpuAdd('start', trim(txroutine)//'_move', 1, uout)

   if(.not.allocated(rosave)) then
      allocate(rosave(3,np_alloc), orisave(3,3,np_alloc), boxlensave(1:3), &
           idmsyssave(1:3), idmosave(1:3,1:np), idmsave(1:3,1:na), forcesave(1:3,1:na))
      rosave = 0.0E+00
      orisave = 0.0E+00
      boxlensave = 0.0E+00
      idmsyssave = 0.0E+00
      idmosave = 0.0E+00
      idmsave = 0.0E+00
      forcesave = 0.0E+00
   end if

! .............. calculate a trial configuration ...............

! ... save old configuration energy and coordinates

   usave      = u
   prsrsave   = prsr
   rosave     = ro
   orisave    = ori
   boxlensave = boxlen
   idmsyssave = idmsys
   idmosave   = idmo
   idmsave    = idm
   forcesave  = force

! ... get trial coordinates

   do ip = 1, np
      ipt = iptpn(ip)
      if ((umbcoord/= 'r') .and. ((ip == ipumb1) .or. (ip == ipumb2))) then
         drostep(1:3,ip) = Zero
         if (umbcoord == 'x') drostep(1,ip) = (Random(iseed)-Half)*dtranall(ipt)
         if (umbcoord == 'y') drostep(2,ip) = (Random(iseed)-Half)*dtranall(ipt)
         if (umbcoord == 'z') drostep(3,ip) = (Random(iseed)-Half)*dtranall(ipt)
      else
         drostep(1,ip) = (Random(iseed)-Half)*dtranall(ipt)
         drostep(2,ip) = (Random(iseed)-Half)*dtranall(ipt)
         drostep(3,ip) = (Random(iseed)-Half)*dtranall(ipt)
      end if
      ro(1:3,ip) = ro(1:3,ip)+drostep(1:3,ip)
      call CheckPartBCTM(1, ro(1:3,ip), lboxoverlap)
      if (lboxoverlap) goto 200
      if ((umbcoord/= 'r') .and. ((ip == ipumb1) .or. (ip == ipumb2))) then
      else
         call GetRandomTrialOri(drotall(ipt), iseed, ori(1,1,ip), oritm)
         ori(1:3,1:3,ip) = oritm(1:3,1:3,1)
      end if
   end do

! ... scale box lengths

   if (lntp) call ScaleLength

! ... if we use nonatomic interaction sites, change to interaction sites

   call SetAtomProp(1, np, lintsite)
 ! if (lradatbox) call CheckAtomBCTM(na, r, lboxoverlap)  ! CheckAtomBCTM does not work here

   if (ltime) call CpuAdd('stop', trim(txroutine)//'_move', 1, uout)

! ............. calculate energy difference ...............

   call UTotal(iStage)
   du%tot = u%tot-usave%tot

! ... if nonatomic interaction sites, change back to geometrical sites

   if (lintsite) call SetAtomPos(1, np, .false.)

! ............. calculate nonenergetic weights .............

   weight = One
   if (lautumb) weight = weight*UmbrellaWeight(0)

! .............. decide new configuration .............

!.. not working if lpolarization =.false. / peo
!   call Metropolis(lboxoverlap,lhsoverlap .or. .not.lidmconv,weight,du%tot*beta)

!  to be simplified

! peo: do not work for UTwoBodyP  lmetropolis = lhsoverlap
!   if (lpolarization) lmetropolis = lhsoverlap .or. .not.lidmconv
   lmetropolis =.false.
   if (lpolarization) lmetropolis =.not.lidmconv

200 continue
   call metropolis(lboxoverlap, lmetropolis, .false., weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) then           ! accepted configuration

   else                                    ! rejected configuration

      u    = usave
      prsr = prsrsave
      ro   = rosave
      ori  = orisave
      call SetAtomProp(1, np, lintsite)
      drostep = Zero
      if (lntp) then
        boxlen = boxlensave
        call SetBoxParam
        if (lewald) call EwaldSetup
      end if
      idmsys = idmsyssave
      idmo   = idmosave
      idm    = idmsave
      force  = forcesave

   end if

   if (lautumb) call UmbrellaUpdate                 ! update weight function for umbrella potential

   if (itest == 1) call TestSimulation

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine MCAllPass

!************************************************************************
!*                                                                      *
!*     ClusterMember                                                    *
!*                                                                      *
!************************************************************************

! ... calculate secondary cluster members using Verlet neighbour list

subroutine ClusterMember(str, lonlyipmove, lnoselftype, radcl, pselectcl)

   use MCModule
   implicit none

   character(3), intent(in) :: str
   logical, intent(in) :: lonlyipmove      ! = .true. include only cluster members to ipmove
   logical, intent(in) :: lnoselftype      ! = .true. exclude particles of type ipt
   real(8), intent(in) :: radcl(*)         ! cluster radius
   real(8), intent(in) :: pselectcl(*)     ! probability of acceptance

   integer(4), save :: naux
   integer(4) :: ip, iploc, ipt, jp, jploc
   real(8)    :: dx, dy, dz, r2
   real(8)    :: Random

   if (str(1:3) == 'old') then    ! calculate nprimpartcl, npcold, nptm, ipnptm, and lptm

      nprimpartcl = nptm         ! set number of primary cluster members
      npclold  = nptm
      naux = nprimpartcl
      if (lonlyipmove) naux = 1
      do iploc = 1, naux
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = 1, nneighpn(ip)
            jp = jpnlist(jploc,ip)
            if (lnoselftype) then
               if (ipt == iptpn(jp)) cycle
            end if
            if (lptm(jp)) cycle
            dx = ro(1,jp)-ro(1,ip)
            dy = ro(2,jp)-ro(2,ip)
            dz = ro(3,jp)-ro(3,ip)
            call PBCr2(dx,dy,dz,r2)
            if (abs(dx) > radcl(ipt)) cycle
            if (abs(dy) > radcl(ipt)) cycle
            if (abs(dz) > radcl(ipt)) cycle
            if (r2 < radcl(ipt)**2) then
               npclold = npclold+1
               if (pselectcl(ipt) >= One) then
                  nptm = nptm+1
                  ipnptm(nptm) = jp
                  lptm(jp) =.true.
               else
                  if (Random(iseed) < pselectcl(ipt)) then
                     nptm = nptm+1
                     ipnptm(nptm) = jp
                     lptm(jp) =.true.
                  end if
               end if
            end if
         end do
      end do
#if defined (_PAR_)
      call all_reduce_nptm_ipnptm_lptm           ! get global nptm, ipnptm, lptm, and npclold
#endif

   else if (str == 'new') then

      npclnew = 0
      do iploc = 1, naux
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = 1, nneighpn(ip)
            jp = jpnlist(jploc,ip)
            if (lptm(jp)) cycle                    ! consider only particles not already in the cluster
            if (lnoselftype) then
               if (ipt == iptpn(jp)) cycle
            end if
            dx = ro(1,jp)-rotm(1,iploc)
            dy = ro(2,jp)-rotm(2,iploc)
            dz = ro(3,jp)-rotm(3,iploc)
            call PBCr2(dx,dy,dz,r2)
            if (abs(dx) > radcl(ipt)) cycle
            if (abs(dy) > radcl(ipt)) cycle
            if (abs(dz) > radcl(ipt)) cycle
            if (r2 < radcl(ipt)**2) then
               npclnew = npclnew+1
!!!            write(*,'(a,2i5,3g11.4)') 'new member, jpt, jp, ro', iptpn(jp), jp, ro(1:3,jp)
            end if
         end do
      end do
#if defined (_PAR_)
      call par_allreduce_int(npclnew, iaux)  ! allreduce npclnew
#endif
      npclnew = npclnew+nptm                     ! add the original partciles
   end if

#if defined (_PAR_)

contains

!........................................................................

subroutine all_reduce_nptm_ipnptm_lptm                ! get global nptm, ipnptm, lptm, and npclold

   integer(4) :: nptmaux(0:mnproc-1), nptmoff

   nptmaux(0:nproc-1) = 0                            ! exchange nptm-1 using nptmaux
   nptmaux(myid) = nptm-1

   call par_allreduce_ints(nptmaux(0),ivaux,nproc)

   nptmoff = 1+sum(nptmaux(0:myid-1))                ! shift nptm-1 data of ipmtlis
   do iploc = nptm-1, 1, -1                           ! starting from location 2
      ipnptm(nptmoff+iploc) = ipnptm(1+iploc)
   end do
   do iploc = 1, min(nptmoff-1,nptm-1)
      ipnptm(1+iploc) = 0
   end do

   nptm = 1+sum(nptmaux(0:nproc-1))                 ! get global nptm

   call par_allreduce_ints(ipnptm(2),ivaux,nptm-1)    ! get global ipnptm

   do iploc = 1, nptm                                ! get global lptm
      lptm(ipnptm(iploc)) =.true.
   end do

   npclold = nptm                                    ! get global npclold

! example:
!
! myid     nptm      ipnptm     lptm
! 0        2         1 5 0 ...   all .false.
! 1        3         1 8 3 0 ... all .false.
! 2        2         1 6 0 ...   all .false.
!
! myid     nptmaux(before) nptmaux(after)
! 0        1 0 0            1 2 1
! 1        0 2 0            1 2 1
! 2        0 0 1            1 2 1
!
! myid     nptoff
! 0        1
! 1        2
! 2        4
!
! myid     ipnptm
! 0        1 5 0 ...
! 1        1 0 8 3 0 ...
! 2        1 0 0 0 6 0 ...
!
! nptm = 1 + 4 = 5
!
! ipnptm = 1 5 8 3 6 0 ...

end subroutine all_reduce_nptm_ipnptm_lptm

#endif

!........................................................................

end subroutine ClusterMember

!************************************************************************
!*                                                                      *
!*     ClusterMemberLList                                               *
!*                                                                      *
!************************************************************************

! ... calculate secondary cluster members using linked lists

subroutine ClusterMemberLList(str, lonlyipmove, lnoselftype, radcl, pselectcl)

   use MCModule
   implicit none

   character(3), intent(in) :: str
   logical, intent(in) :: lonlyipmove      ! = .true. include only cluster members to ipmove
   logical, intent(in) :: lnoselftype      ! =.true. exclude particles of type ipt
   real(8), intent(in) :: radcl(*)         ! cluster radius
   real(8), intent(in) :: pselectcl(*)     ! probability of acceptance

   integer(4), save :: naux
   integer(4) :: ip, iploc, ipt, jp, icell, Getncellllist
   real(8)    :: dx, dy, dz, r2
   real(8)    :: Random

   if (str(1:3) == 'old') then  ! calculate nprimpartcl, npcold, nptm, ipnptm, and lptm

      nprimpartcl = nptm         ! set number of primary cluster members
      npclold  = nptm
      naux = nprimpartcl
      if (lonlyipmove) naux = 1
      do iploc = 1, naux
         ip = ipnptm(iploc)
         ipt = iptpn(ip)

         call SetNCell(ro(1:3,ip), rcut, lcellllist)
         do icell = 1, Getncellllist()
           if (.not.lcellllist(icell)) cycle
           jp = headllist(icell)

           do
            if (jp == 0) exit
            if (lnoselftype) then
               if (ipt == iptpn(jp)) goto 102
            end if
            if (lptm(jp)) goto 102

            dx = ro(1,jp)-ro(1,ip)
            dy = ro(2,jp)-ro(2,ip)
            dz = ro(3,jp)-ro(3,ip)
            call PBCr2(dx,dy,dz,r2)
            if (abs(dx) > radcl(ipt)) goto 102
            if (abs(dy) > radcl(ipt)) goto 102
            if (abs(dz) > radcl(ipt)) goto 102
            if (r2 < radcl(ipt)**2) then
               npclold = npclold+1
               if (pselectcl(ipt) >= One) then
                  nptm = nptm+1
                  ipnptm(nptm) = jp
                  lptm(jp) =.true.
               else
                  if (Random(iseed) < pselectcl(ipt)) then
                     nptm = nptm+1
                     ipnptm(nptm) = jp
                     lptm(jp) =.true.
                  end if
               end if
            end if

 102     continue
         jp = jpllist(jp)
        end do
      end do
    end do
#if defined (_PAR_)
      call all_reduce_nptm_ipnptm_lptm           ! get global nptm, ipnptm, lptm, and npclold
#endif

   else if (str == 'new') then
      npclnew = 0
      do iploc = 1, naux
         ip = ipnptm(iploc)
         ipt = iptpn(ip)

         call SetNCell(ro(1:3,ip), rcut, lcellllist)
         do icell = 1, Getncellllist()
            if (.not.lcellllist(icell)) cycle

            if (headllist(icell) == 0) cycle
            jp = headllist(icell)

            do
               if (lptm(jp) .or. ip == jp) goto 103
               if (lnoselftype) then
                  if (ipt == iptpn(jp)) goto 103
               end if
               dx = ro(1,jp)-rotm(1,iploc)
               dy = ro(2,jp)-rotm(2,iploc)
               dz = ro(3,jp)-rotm(3,iploc)
               call PBCr2(dx,dy,dz,r2)
               if (abs(dx) > radcl(ipt)) goto 103
               if (abs(dy) > radcl(ipt)) goto 103
               if (abs(dz) > radcl(ipt)) goto 103
               if (r2 < radcl(ipt)**2) then
                  npclnew = npclnew+1
!!!               write(*,'(a,2i5,3g11.4)') 'new member, jpt, jp, ro', iptpn(jp), jp, ro(1:3,jp)
               end if
 103        continue
            if (jpllist(jp) == 0) exit
            jp = jpllist(jp)
           end do
         end do

      end do
#if defined (_PAR_)
      call par_allreduce_int(npclnew,iaux)    ! allreduce npclnew
#endif
      npclnew = npclnew+nptm                      ! add the original partciles
   end if

#if defined (_PAR_)

contains

!........................................................................

subroutine all_reduce_nptm_ipnptm_lptm                ! get global nptm, ipnptm, lptm, and npclold

   integer(4) :: nptmaux(0:mnproc-1), nptmoff

   nptmaux(0:nproc-1) = 0                             ! exchange nptm-1 using nptmaux
   nptmaux(myid) = nptm-1

   call par_allreduce_ints(nptmaux(0),ivaux,nproc)

   nptmoff = 1+sum(nptmaux(0:myid-1))                 ! shift nptm-1 data of ipmtlis
   do iploc = nptm-1, 1, -1                           ! starting from location 2
      ipnptm(nptmoff+iploc) = ipnptm(1+iploc)
   end do
   do iploc = 1, min(nptmoff-1,nptm-1)
      ipnptm(1+iploc) = 0
   end do

   nptm = 1+sum(nptmaux(0:nproc-1))                   ! get global nptm

   call par_allreduce_ints(ipnptm(2),ivaux,nptm-1)     ! get global ipnptm

   do iploc = 1, nptm                                 ! get global lptm
      lptm(ipnptm(iploc)) =.true.
   end do

   npclold = nptm                                     ! get global npclold

! example:
!
! myid     nptm      ipnptm     lptm
! 0        2         1 5 0 ...   all .false.
! 1        3         1 8 3 0 ... all .false.
! 2        2         1 6 0 ...   all .false.
!
! myid     nptmaux(before) nptmaux(after)
! 0        1 0 0            1 2 1
! 1        0 2 0            1 2 1
! 2        0 0 1            1 2 1
!
! myid     nptoff
! 0        1
! 1        2
! 2        4
!
! myid     ipnptm
! 0        1 5 0 ...
! 1        1 0 8 3 0 ...
! 2        1 0 0 0 6 0 ...
!
! nptm = 1 + 4 = 5
!
! ipnptm = 1 5 8 3 6 0 ...

end subroutine all_reduce_nptm_ipnptm_lptm

#endif

!........................................................................

end subroutine ClusterMemberLList

!************************************************************************
!*                                                                      *
!*     GetRandomTrialPos                                                *
!*                                                                      *
!************************************************************************

! ... calculate translational displacement and trial position

subroutine GetRandomTrialPos(dtran, iseed, nptm, ipnptm, ro, rotm, drotm)

   implicit none
   real(8),    intent(in)  :: dtran       ! translational displacement
                                          ! > 0 square region, < 0 spherical region
   integer(4), intent(inout)  :: iseed       ! random number seed
   integer(4), intent(in)  :: nptm        ! number of moving particles
   integer(4), intent(in)  :: ipnptm(*)  ! moving particles
   real(8),    intent(in)  :: ro(3,*)     ! particle position
   real(8),    intent(out) :: rotm(3,*)   ! particle position, for trial configuration
   real(8),    intent(out) :: drotm(3,*)  ! suggested particle move

   real(8), parameter :: Zero = 0.0d0, Half = 0.5d0, One = 1.0d0, Two = 2.0d0
   integer(4) :: iploc, ip
   real(8) :: dx, dy, dz
   real(8) :: Random

! ... get translational displacement

   if (dtran >= Zero) then                     ! in a box
      dx = (Random(iseed)-Half)*dtran
      dy = (Random(iseed)-Half)*dtran
      dz = (Random(iseed)-Half)*dtran
   else                                       ! in a sphere of radius 0.5
      do
         dx = Random(iseed)-Half
         dy = Random(iseed)-Half
         dz = Random(iseed)-Half
         if (dx**2 + dy**2 + dz**2 < Half**2) exit
      end do
      dx = dx*abs(dtran)
      dy = dy*abs(dtran)
      dz = dz*abs(dtran)
   end if

! ... get trial coordinates and store translational displacemnt

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      rotm(1,iploc) = ro(1,ip)+dx
      rotm(2,iploc) = ro(2,ip)+dy
      rotm(3,iploc) = ro(3,ip)+dz
      drotm(1,iploc) = dx
      drotm(2,iploc) = dy
      drotm(3,iploc) = dz
   end do

end subroutine GetRandomTrialPos

!************************************************************************
!*                                                                      *
!*     GetRandomTrialOri                                                *
!*                                                                      *
!************************************************************************

! ... calculate a trail orientation

subroutine GetRandomTrialOri(drot, iseed, ori, oritm)

   implicit none

   real(8),    intent(in)  :: drot         ! rotational displacement
                                           ! > 0 rotation around a box axis, < 0 rotation around a particle axis
   integer(4), intent(inout)  :: iseed        ! random number seed
   real(8),    intent(in)  :: ori(3,3)     ! particle orientation
   real(8),    intent(out) :: oritm(3,3)   ! particle orientation, for trial configuration

   integer(4) :: iaxis, m1, m2, m3
   real(8)    :: drottemp, cr, sr
   real(8)    :: Random

   iaxis = int(1.0d0+3.0d0*Random(iseed))
   iaxis = max(1,int(min(int(iaxis),3)))
   m1 = 1+mod(int(iaxis  ),3)
   m2 = 1+mod(int(iaxis+1),3)
   m3 = 1+mod(int(iaxis-1),3)
   drottemp = (Random(iseed)-0.5d0)*drot
   cr = cos(drottemp)
   sr = sin(drottemp)
   if (drot >= 0) then
      oritm(m1,1:3) = cr*ori(m1,1:3) - sr*ori(m2,1:3)
      oritm(m2,1:3) = sr*ori(m1,1:3) + cr*ori(m2,1:3)
      oritm(m3,1:3) =                                     ori(m3,1:3)
   else if (drot < 0) then
      oritm(1:3,m1) = cr*ori(1:3,m1) - sr*ori(1:3,m2)
      oritm(1:3,m2) = sr*ori(1:3,m1) + cr*ori(1:3,m2)
      oritm(1:3,m3) =                                     ori(1:3,m3)
   end if

end subroutine GetRandomTrialOri

!************************************************************************
!*                                                                      *
!*     GetFixedTrialOri                                                 *
!*                                                                      *
!************************************************************************

! ... get fixed trial orientation

subroutine GetFixedTrialOri(nptm, ipnptm, ori, oritm)

   implicit none

   integer(4), intent(in)  :: nptm         ! number of moving particles
   integer(4), intent(in)  :: ipnptm(*)    ! moving particles
   real(8),    intent(in)  :: ori(3,3,*)   ! particle orientation
   real(8),    intent(out) :: oritm(3,3,*) ! particle orientation, for trial configuration

   integer(4) :: iploc, ip

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      oritm(1:3,1:3,iploc) = ori(1:3,1:3,ip)
   end do

end subroutine GetFixedTrialOri

!************************************************************************
!*                                                                      *
!*     GetRotatedTrialOri                                               *
!*                                                                      *
!************************************************************************

! ... get rotated trial orientation

subroutine GetRotatedTrialOri(rotaxis, alpha, nptm, ipnptm, ori, oritm)

   implicit none

   real(8),    intent(in)  :: rotaxis(3)      ! direction of rotation axis
   real(8),    intent(in)  :: alpha           ! rotation angle
   integer(4), intent(in)  :: nptm            ! number of rotating particles
   integer(4), intent(in)  :: ipnptm(*)       ! moving particles
   real(8),    intent(in)  :: ori(3,3,*)      ! particle orientation
   real(8),    intent(out) :: oritm(3,3,*)    ! particle orientation, for trial configuration

   integer(4) :: iploc, ip, m
   real(8)    :: rotmat(3,3)

! ... get rotation matrix rotmat

   call AxisAngToOri(rotaxis,alpha,rotmat)

! ... rotate and get oritm

   do iploc = 1,nptm
      ip = ipnptm(iploc)
      do m = 1,3
         oritm(1,m,iploc) = rotmat(1,1)*ori(1,m,ip)+rotmat(1,2)*ori(2,m,ip)+rotmat(1,3)*ori(3,m,ip)
         oritm(2,m,iploc) = rotmat(2,1)*ori(1,m,ip)+rotmat(2,2)*ori(2,m,ip)+rotmat(2,3)*ori(3,m,ip)
         oritm(3,m,iploc) = rotmat(3,1)*ori(1,m,ip)+rotmat(3,2)*ori(2,m,ip)+rotmat(3,3)*ori(3,m,ip)
      end do
   end do

end subroutine GetRotatedTrialOri

!************************************************************************
!*                                                                      *
!*     RotateSetPart                                                    *
!*                                                                      *
!************************************************************************

! ... rigid-body rotation of a set of particles

subroutine RotateSetPart(rotaxis, alpha, rotorigin, nptm, ipnptm, ro, rotm, drotm)

   implicit none

   real(8),    intent(in)  :: rotaxis(3)      ! direction of rotation axis
   real(8),    intent(in)  :: alpha           ! rotation angle
   real(8),    intent(in)  :: rotorigin(3)    ! origin of rotation
   integer(4), intent(in)  :: nptm            ! number of rotating particles
   integer(4), intent(in)  :: ipnptm(*)       ! moving particles
   real(8),    intent(in)  :: ro(3,*)         ! particle position
   real(8),    intent(out) :: rotm(3,*)       ! particle position, trial configuration
   real(8),    intent(out) :: drotm(3,*)      ! suggested particle move

   integer(4) :: iploc, ip
   real(8)    :: rotmat(3,3), dx, dy, dz, tempx, tempy, tempz

! ... get rotation matrix rotmat

   call AxisAngToOri(rotaxis, alpha, rotmat)

! ... rotate and get rotm and drotm

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      dx = ro(1,ip)-rotorigin(1)
      dy = ro(2,ip)-rotorigin(2)
      dz = ro(3,ip)-rotorigin(3)
      tempx = rotmat(1,1)*dx+rotmat(1,2)*dy+rotmat(1,3)*dz
      tempy = rotmat(2,1)*dx+rotmat(2,2)*dy+rotmat(2,3)*dz
      tempz = rotmat(3,1)*dx+rotmat(3,2)*dy+rotmat(3,3)*dz
      rotm(1,iploc) = rotorigin(1)+tempx
      rotm(2,iploc) = rotorigin(2)+tempy
      rotm(3,iploc) = rotorigin(3)+tempz
      drotm(1,iploc) = rotm(1,iploc)-ro(1,ip)
      drotm(2,iploc) = rotm(2,iploc)-ro(2,ip)
      drotm(3,iploc) = rotm(3,iploc)-ro(3,ip)
   end do

end subroutine RotateSetPart

!************************************************************************
!*                                                                      *
!*     FizedZCoord                                                      *
!*                                                                      *
!************************************************************************

! ... no displacement in z-direction

subroutine FixedZCoord
   use MCModule
   implicit none
   integer(4) :: iploc, ip
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      rotm(3,iploc) = ro(3,ip)
      drotm(3,iploc) = Zero
   end do
end subroutine FixedZCoord

!************************************************************************
!*                                                                      *
!*     FizedXYCoord                                                     *
!*                                                                      *
!************************************************************************

! ... no displacement in the xy-plane

subroutine FixedXYCoord
   use MCModule
   implicit none
   integer(4) :: iploc, ip
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      rotm(1:2,iploc) = ro(1:2,ip)
      drotm(1:2,iploc) = Zero
   end do
end subroutine FixedXYCoord

!************************************************************************
!*                                                                      *
!*     FixedChainStart                                                  *
!*                                                                      *
!************************************************************************

! ... put rotm = ro and drotm to zero if particle is the first segment in a chain

subroutine FixedChainStart
   use MCModule
   implicit none
   character(40), parameter :: txroutine ='FixedChainStart'
   if (nptm /= 1) call Stop(txroutine,'nptm /=1', uout)
   if (icnpn(ipmove) /= 0) then                      ! particle ipmove belonging to a chain ?
      if (isegpn(ipmove) == 1) then                  ! first particle in the chain ?
         rotm(1:3,1) = ro(1:3,ipmove)                ! no displacement of ipmove
         drotm(1:3,1) = Zero                         ! no displacement of ipmove
      end if
   end if
end subroutine FixedChainStart

!************************************************************************
!*                                                                      *
!*     ShiftZDirection                                                  *
!*                                                                      *
!************************************************************************

! ... center the system at z = 0 as much as possible (only for lbccyl)

subroutine ShiftZDirection
   use MCModule
   implicit none
   character(40), parameter :: txroutine ='ShiftZDirection'
   real(8) :: zcom, dist, dz
   if (.not.lbccyl) call Stop(txroutine, '.not.lbccyl''', uout)
   zcom = Half*(ro(3,1)+ro(3,2))                         ! center-of-mass in the z-direction of particle 1 and 2
   dist = Half*cyllen + minval(ro(3,1:np)*sign(One,zcom))! smallest distance between relevant surface and closest particle
   dz = - sign(min(abs(zcom),dist),zcom)                 ! actual displacement
   ro(3,1:np) = ro(3,1:np)+dz                            ! displace particles
   call SetAtomPos(1,np,.false.)                         ! get atom positions
end subroutine ShiftZDirection

!************************************************************************
!*                                                                      *
!*     CheckPartBCTM                                                    *
!*                                                                      *
!************************************************************************

! ... check boundary conditions of trial particle coordinate

subroutine CheckPartBCTM(nploc, rr, lboxoverlap)

   use MCModule
   implicit none

   integer(4)   , intent(in)    :: nploc        ! number of particles to check
   real(8)      , intent(inout) :: rr(3,*)      ! in: trial coordinate
                                                ! out: possibly multiple reflected trial coordinate
   logical      , intent(out)   :: lboxoverlap  !=.true. if outside the box
                                                !=.false. if inside the box
   integer(4) :: ip

   lboxoverlap =.false.
   do ip = 1, nploc
      if (lbcbox) then
         if (txbc == 'xyz') then
            call PBC(rr(1,ip),rr(2,ip),rr(3,ip))
         else if (txbc == 'xy') then
            if (abs(rr(3,ip)) > boxlen2(3)) lboxoverlap =.true.
            if (lboxoverlap) exit
            call PBC(rr(1,ip),rr(2,ip),rr(3,ip))
         else if (txbc == 'z') then
            if (abs(rr(1,ip)) > boxlen2(1)) lboxoverlap =.true.
            if (abs(rr(2,ip)) > boxlen2(2)) lboxoverlap =.true.
            if (lboxoverlap) exit
            call PBC(rr(1,ip),rr(2,ip),rr(3,ip))
         end if
      end if
      if (lbcrd .or. lbcto) then
         call PBC(rr(1,ip),rr(2,ip),rr(3,ip))
      else if (lbcsph) then
         if (rr(1,ip)**2+rr(2,ip)**2+rr(3,ip)**2 > sphrad2) lboxoverlap =.true.
         if (lboxoverlap) exit
      else if (lbccyl) then
         if (rr(1,ip)**2+rr(2,ip)**2 > cylrad2) lboxoverlap =.true.
         if (abs(rr(3,ip)) > Half*cyllen) lboxoverlap =.true.
         if (lboxoverlap) exit
      else if (lbcell) then
         if (sum((rr(1:3,ip)*ellradi(1:3))**2) > One) lboxoverlap =.true.
         if (lboxoverlap) exit
      end if
   end do

end subroutine CheckPartBCTM

!************************************************************************
!*                                                                      *
!*     AddNeighbours                                                    *
!*                                                                      *
!************************************************************************

subroutine AddNeighbours()

use MCModule
implicit none

integer(4)        :: iploc, ip, ipn, icn, ict, isn, nptmadd

nptmadd=nptm            !storage location for next particle

!write (*,*) 'start adding neighbours'
do iploc=1, nptm
!   write (*,*) 'checking ', iploc
   ip = ipnptm(iploc)
   ict = ictpn(ip)
   icn = icnpn(ip)
   isn = isegpn(ip)


   if (ict<1) exit        !particle not in a chain


   if (isn > 1) then
      !check for previous particle
      ipn = ipnsegcn(isn-1,icn)
      if (.not.lptm(ipn)) then    !particle not yet in trial move
         nptmadd=nptmadd+1
         rotm(1:3,nptmadd)=ro(1:3,ipn)
         ipnptm(nptmadd)=ipn
         iptmpn(ipn)=nptmadd
         lptm(ipn)=.true.
      end if
   end if

   if (isn<npct(ict)) then
      !check for next particle
      ipn = ipnsegcn(isn+1,icn)
      if (.not.lptm(ipn)) then    !particle not yet in trial move
         nptmadd=nptmadd+1
         rotm(1:3,nptmadd)=ro(1:3,ipn)
         ipnptm(nptmadd)=ipn
         iptmpn(ipn)=nptmadd
         lptm(ipn)=.true.
      end if
   end if
end do

!write(*,*) nptmadd-nptm, ' particles added to the trial move'
lptmdutwob = .true.
nptm = nptmadd

end subroutine AddNeighbours

!************************************************************************
!*                                                                      *
!*     UpdateOri                                                        *
!*                                                                      *
!************************************************************************

! ... get fixed trial orientation

subroutine UpdateOri()
   use MCModule
   implicit none

   integer(4)     :: iploc, ip, ipn, ict, icn, isn
   real(8)        :: r1(3), r2(3), r3(3), r12(3), r32(3)
   !character(20)  :: format = '(2I3, 4F7.2)'

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ict = ictpn(ip)
      icn = icnpn(ip)
      isn = isegpn(ip)

      oritm(1:3,1:3,iploc) = ori(1:3,1:3,ip)

      if (ict>0) then                                                                              ! Particle is part of a chain
!write(*,*) 'is =', isn                                       ! update trial orientation
         r2 = rotm(:,iploc)
!write(*,format) ip, ip, ro(:,ip)
!write(*,format) ip, ip, r2
         if ( isn == 1) then                                                                            ! if first segment
            !determine correct particle positions

            ipn = ipnsegcn(isn+1,icn)
            r3 = ro(:,ipn)
!write(*,format) ipn, ip, r3
            if (lptm(ipn)) r3 = rotm(:,iptmpn(ipn))
!write(*,format) ipn, ip, r3

            ! set particle orientation

            r32 = r3-r2
            call PBC(r32(1),r32(2),r32(3))

            call SetPartOriBond('end', r32, r32, oritm(1,1,iploc))

         else if (isn==npct(ict)) then                                                                  !  if last particle in a chain

            ipn = ipnsegcn(isn-1,icn)
            r1 = ro(:,ipn)
!write(*,format) ipn, ip, r1
            if (lptm(ipn)) r1 = rotm(:,iptmpn(ipn))
!write(*,format) ipn, ip, r1

            r12 = r1-r2
            call PBC(r12(1),r12(2),r12(3))

            call SetPartOriBond('end', r12, r12, oritm(1,1,iploc))
         else                                                                                          !  particle is midpoint orient to both vectors
            ipn = ipnsegcn(isn-1,icn)
            r1 = ro(:,ipn)
            if (lptm(ipn)) r1 = rotm(:,iptmpn(ipn))
!write(*,format) ipn, ip, r1

            ipn = ipnsegcn(isn+1,icn)
            r3 = ro(:,ipn)
            if (lptm(ipn)) r3 = rotm(:,iptmpn(ipn))
!write(*,format) ipn, ip, r3

            r12 = r1-r2
            call PBC(r12(1),r12(2),r12(3))

            r32 = r3-r2
            call PBC(r32(1),r32(2),r32(3))

            call SetPartOriBond('mid', r12, r32, oritm(1,1,iploc))
         end if
!        call UpdateOriTest()
!       write(*,*) 'orientation set'
      end if
   end do
end subroutine UpdateOri

subroutine UpdateOriTest(iploc, ip, ipn, r1, r2, r3)
use MCModule
implicit none
real(8), intent(in)        ::r1(3),r2(3),r3(3)
integer(4), intent(in)        :: iploc, ip, ipn
character(20)  :: format = '(2I3, 4F7.2)'

if (.not. all(oritm(1:3,1:3,iploc).eq.ori(1:3,1:3,ip))) then
  write(*,*) '------------------------------------------------------------------------'
  write(*,*) oritm(1:3,1:3,iploc)
  write(*,*)
  write(*,*) ori(1:3,1:3,ip)
  write(*,*) '------------------------------------------------------------------------'
end if

if (.not. all(rotm(1:3,iploc).eq.ro(1:3,ip))) then
  write(*,*) '------------------------------------------------------------------------'
  write(*,format) ipn, ip, r1
  write(*,format) ip, ip, r2
  write(*,format) ipn, ip, r3
  write(*,*) '------------------------------------------------------------------------'
end if

end subroutine UpdateOriTest

!************************************************************************
!*                                                                      *
!*     SetTrialAtomProp                                                 *
!*                                                                      *
!************************************************************************

! ... set trial atom properties from trial particles properties

subroutine SetTrialAtomProp

   use MCModule
   implicit none

   integer(4) :: ip, iploc, ipt, ia, ialoc

   if (lmonoatom) then             ! set rtm
      natm = nptm
      do iploc = 1, nptm
         ianatm(iploc) = ipnptm(iploc)
         rtm(1,iploc) = rotm(1,iploc)
         rtm(2,iploc) = rotm(2,iploc)
         rtm(3,iploc) = rotm(3,iploc)
      end do
   else                          ! set natm, ianatm, and rtm
      natm = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do ia = ianpn(ip), ianpn(ip)+napt(iptpn(ip))-1
            natm = natm+1
            ianatm(natm) = ia
            ialoc = ia-(ianpn(ip)-1)
            rtm(1,natm) = rotm(1,iploc)+oritm(1,1,iploc)*ra(1,ialoc,ipt)+oritm(1,2,iploc)*ra(2,ialoc,ipt)+oritm(1,3,iploc)*ra(3,ialoc,ipt)
            rtm(2,natm) = rotm(2,iploc)+oritm(2,1,iploc)*ra(1,ialoc,ipt)+oritm(2,2,iploc)*ra(2,ialoc,ipt)+oritm(2,3,iploc)*ra(3,ialoc,ipt)
            rtm(3,natm) = rotm(3,iploc)+oritm(3,1,iploc)*ra(1,ialoc,ipt)+oritm(3,2,iploc)*ra(2,ialoc,ipt)+oritm(3,3,iploc)*ra(3,ialoc,ipt)
         end do
      end do
   end if

   if (ldipole .or. ldipolesph) then                ! set diptm
      natm = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do ia = ianpn(ip), ianpn(ip)+napt(iptpn(ip))-1
            natm = natm+1
            ialoc = ia-(ianpn(ip)-1)
            diptm(1,natm) = oritm(1,1,iploc)*dipa(1,ialoc,ipt)+oritm(1,2,iploc)*dipa(2,ialoc,ipt)+oritm(1,3,iploc)*dipa(3,ialoc,ipt)
            diptm(2,natm) = oritm(2,1,iploc)*dipa(1,ialoc,ipt)+oritm(2,2,iploc)*dipa(2,ialoc,ipt)+oritm(2,3,iploc)*dipa(3,ialoc,ipt)
            diptm(3,natm) = oritm(3,1,iploc)*dipa(1,ialoc,ipt)+oritm(3,2,iploc)*dipa(2,ialoc,ipt)+oritm(3,3,iploc)*dipa(3,ialoc,ipt)
         end do
      end do
   end if

end subroutine SetTrialAtomProp

!************************************************************************
!*                                                                      *
!*     CheckAtomBCTM                                                    *
!*                                                                      *
!************************************************************************

! ... check boundary conditions of trial atom coordinate

subroutine CheckAtomBCTM(naloc, rr, lboxoverlap)

   use MCModule
   implicit none

   integer(4)   , intent(in)    :: naloc        ! number of atoms to check
   real(8)      , intent(in)    :: rr(3,*)      ! trial coordinate
   logical      , intent(out)   :: lboxoverlap  !=.true. if outside the box
                                                !=.false. if inside the box

   character(40), parameter :: txroutine ='CheckAtomMCTM'
   integer(4) :: ialoc, ia, iat

   lboxoverlap =.false.
   do ialoc = 1, naloc
      ia = ianatm(ialoc)
      iat = iatan(ia)

      if (lbcbox) then
         if (txbc == 'xyz') then
            continue
         else if (txbc == 'xy') then
            if (abs(rr(3,ialoc)) > boxlen2(3)-radat(iat)) lboxoverlap =.true.
            if (lboxoverlap) exit
         else if (txbc == 'z') then
            if (abs(rr(1,ialoc)) > boxlen2(1)-radat(iat)) lboxoverlap =.true.
            if (abs(rr(2,ialoc)) > boxlen2(2)-radat(iat)) lboxoverlap =.true.
            if (lboxoverlap) exit
         end if
      else if (lbcrd .or. lbcto) then
         continue
      else if (lbcsph) then
         if (rr(1,ialoc)**2+rr(2,ialoc)**2+rr(3,ialoc)**2 > (sphrad-radat(iat))**2) lboxoverlap =.true.
         if (lboxoverlap) exit
      else if (lbccyl) then
         if (rr(1,ialoc)**2+rr(2,ialoc)**2 > (cylrad-radat(iat))**2) lboxoverlap =.true.
         if (abs(rr(3,ialoc)) > Half*cyllen-radat(iat)) lboxoverlap =.true.
         if (lboxoverlap) exit
      else if (lbcell) then
         call Stop(txroutine, 'not implemented for lbcell', uout)
         if (sum((rr(1:3,ialoc)*ellradi(1:3))**2) > One) lboxoverlap =.true.
         if (lboxoverlap) exit
      end if
   end do

end subroutine CheckAtomBCTM

!************************************************************************
!*                                                                      *
!*     Metropolis                                                       *
!*                                                                      *
!************************************************************************

! ... Metropolis test

subroutine Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, dured)

   use MCModule
   implicit none

   logical, intent(in)  :: lboxoverlap ! =.true. if box overlap
   logical, intent(in)  :: lhsoverlap  ! =.true. if hard-core overlap
   logical, intent(in)  :: lhepoverlap ! =.true. if hard-external-potential overlap
   real(8), intent(in)  :: weight      ! nonenergetic weight
   real(8), intent(in)  :: dured       ! reduced energy difference

   character(40), parameter :: txroutine ='Metropolis'
   real(8) :: fac
   real(8) :: Random

   if (lboxoverlap) then
      ievent = imcboxreject
   else if (lhsoverlap) then
      ievent = imchsreject
   else if (lhepoverlap) then
      ievent = imchepreject
   else if (dured > expmax) then
      ievent = imcreject
   else if (dured <-expmax) then
      if (weight == Zero) then
         ievent = imcreject
      else
         ievent = imcaccept
      end if
   else
      fac = weight*exp(-dured)
      if (fac > One) then
         ievent = imcaccept
      else if (fac > Random(iseed)) then
         ievent = imcaccept
      else
         ievent = imcreject
      end if
   end if

!  call TestMetropolis(uout)

contains

!........................................................................

subroutine TestMetropolis(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,t20, l    )') 'lboxoverlap         ', lboxoverlap
   write(unit,'(a,t20, l    )') 'lhsoverlap          ', lhsoverlap
   write(unit,'(a,t20, l    )') 'lhepoverlap         ' ,lhepoverlap
   write(unit,'(a,t20,3f15.5)') 'weight              ', weight
   write(unit,'(a,t20,1f15.5)') 'dured               ', dured
   write(unit,'(a,t20, i5   )') 'ievent              ', ievent
end subroutine TestMetropolis

!........................................................................

end subroutine Metropolis

!************************************************************************
!*                                                                      *
!*     MCUpdate                                                         *
!*                                                                      *
!************************************************************************

! ... update energies and coordinates after accepted mc trial move

subroutine MCUpdate

   use MCModule
   use CellListModule, only: UpdateCellIp
   implicit none

   integer(4) :: ip, iploc

                                 u%tot            = u%tot            + du%tot
                                 u%twob(0:nptpt)  = u%twob(0:nptpt)  + du%twob(0:nptpt)
   if (lcharge .and. lewald)     u%rec            = u%rec            + du%rec
   if (lweakcharge .and. lewald) u%rec            = u%rec            + du%rec
   if (ldipole .or. ldipolesph)  u%stat           = u%stat           + du%stat
   if (ldieldis)                 u%oneb           = u%oneb           + du%oneb
   if (lchain)                   u%bond           = u%bond           + du%bond
   if (lchain)                   u%angle          = u%angle          + du%angle
   if (lclink)                   u%crosslink      = u%crosslink      + du%crosslink
   if (luext)                    u%external       = u%external       + du%external

   if (lewald) call EwaldUpdateArray
   if (luext)  call UExternalUpdate(iptmove)

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ro(1:3,ip) = rotm(1:3,iploc)                             ! position
      if (lpolyatom .or. lellipsoid .or. lsuperball .or. lfixedori) ori(1:3,1:3,ip) = oritm(1:3,1:3,iploc)    ! orientation
!  not sure that lfixedori is needed in the line above Jos
      call SetAtomProp(ip, ip, .false.)                        ! atom and dipole
      drostep(1:3,ip) = drostep(1:3,ip)+drotm(1:3,iploc)       ! displacement
      if(lCellList) call UpdateCellIp(ip)
   end do

end subroutine MCUpdate

!************************************************************************
!*                                                                      *
!*     Restorelptm                                                      *
!*                                                                      *
!************************************************************************

! ... restore lptm

subroutine Restorelptm(nptm, ipnptm, lptm)

   implicit none

   integer(4), intent(in)  :: nptm
   integer(4), intent(in)  :: ipnptm(*)
   logical,    intent(out) :: lptm(*)
   integer(4) :: iploc

   do iploc = 1, nptm
      lptm(ipnptm(iploc)) =.false.
   end do

end subroutine Restorelptm

!************************************************************************
!*                                                                      *
!*     MCAver                                                           *
!*                                                                      *
!************************************************************************

! ... calculate means of mc specific variables

subroutine MCAver(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCAver'
   integer(8), allocatable, save :: nys1(:,:,:), nys2(:,:,:)        ! for non-cluster move
   integer(8), allocatable, save :: npcl1s1(:,:,:), npcl1s2(:,:,:)  ! for cluster1 move
   integer(8), allocatable, save :: npcl2s1(:,:), npcl2s2(:,:)      ! for cluster2 move
   integer(4), allocatable, save :: nprimpartclmax(:)               ! for cluster2 move
   integer(8), allocatable, save :: nazs1(:), nazs2(:)              ! for charge move
   integer(4) :: ipt

   if (slave) return

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      if (.not.allocated(nys1)) then
         allocate(nys1(0:nevent,0:npt,nmovetype), nys2(0:nevent,0:npt,nmovetype), &
      npcl1s1(2,npt,nmovetype), npcl1s2(2,npt,nmovetype), npcl2s1(2,npt), npcl2s2(2,npt), nprimpartclmax(npt))
         nys1 = 0
         nys2 = 0
         npcl1s1 = 0
         npcl1s2 = 0
         npcl2s1 = 0
         npcl2s2 = 0
         nprimpartclmax = 0
      end if
      if (.not.allocated(nazs1)) then
         allocate(nazs1(1:nat), nazs2(1:nat))
         nazs1 = 0
         nazs2 = 0
      end if

      nys1           = 0
      npcl1s1        = 0
      npcl2s1        = 0
      nprimpartclmax = 0
      if (txstart == 'continue') read(ucnf) nys1, npcl2s1, nprimpartclmax, npcl1s1, nazs1

   case (iBeforeMacrostep)

      nys2    = 0
      npcl1s2 = 0
      npcl2s2 = 0

   case (iSimulationStep)

      nys2(ievent,iptmove,imovetype) = nys2(ievent,iptmove,imovetype) + 1                ! sample event
      if (ievent == imcaccept) then
         npcl1s2(1,iptmove,imovetype) = npcl1s2(1,iptmove,imovetype) + nptm              ! sample number of displaced particles
         npcl1s2(2,iptmove,imovetype) = npcl1s2(2,iptmove,imovetype) + (npclnew-npclold) ! sample change of net joining members
         if ((imovetype == ispartcl2move) .or. (imovetype == ibrushcl2move)) then        ! additional for cluster2 move
            npcl2s2(1,iptmove) = npcl2s2(1,iptmove) + nptm                               ! sample number of displaced particles
            npcl2s2(2,iptmove) = npcl2s2(2,iptmove) + nprimpartcl                        ! sample number of primary particles
            if (nprimpartcl > nprimpartclmax(iptmove)) nprimpartclmax(iptmove) = nprimpartcl
         end if
      end if

   case (iAfterMacrostep)

      do imovetype = 1, nmovetype
         do ievent = 1, nevent
            nys2(ievent,0,imovetype) = sum(nys2(ievent,1:npt,imovetype))  ! sum over ipt
         end do
         do ipt = 0, npt
            nys2(0,ipt,imovetype) = sum(nys2(1:nevent,ipt,imovetype))     ! sum over ievent
         end do
      end do

      nys1    = nys1    + nys2
      npcl1s1 = npcl1s1 + npcl1s2
      npcl2s1 = npcl2s1 + npcl2s2

      call MCAverWrite(nys2, npcl1s2, npcl2s2)

      write(ucnf) nys1, npcl2s1, nprimpartclmax, npcl1s1, nazs1

   case (iAfterSimulation)

      call MCAverWrite(nys1, npcl1s1, npcl2s1)

      if (nys1(0,0,ispartcl2move) > 0) &
         write(uout,'(a,t50,6(5x,i10,5x))') 'max. num. of clusters in aggr. ', nprimpartclmax

!     call MCAverData
      deallocate(nys1, nys2, npcl1s1, npcl1s2, npcl2s1, npcl2s2,  nprimpartclmax)

   end select

contains

!........................................................................

subroutine MCAverWrite(ny, npcl, nprimpartcl)

   use MollibModule, only: InvInt
   integer(8), intent(in) :: ny(0:nevent,0:npt,1:nmovetype)
   integer(8), intent(in) :: npcl(2,npt,1:nmovetype)
   integer(8), intent(in) :: nprimpartcl(2,npt)

   integer(4) :: ipt
   real(8)    :: any(0:nevent,0:npt), ancl1(2,npt), ancl2(2,npt)

   call WriteHead(2, 'mc statistics', uout)
   do imovetype = 1, nmovetype
      if (ny(0,0,imovetype) == 0) cycle
      if (imovetype > 1) write(uout,'()')
      do ipt = 0, npt
         any(0:nevent,ipt) = ny(0:nevent,ipt,imovetype)*InvInt(ny(0,ipt,imovetype))
      end do
      do ipt = 1, npt
         ancl1(1:2,ipt) = npcl(1:2,ipt,imovetype)*InvInt(ny(imcaccept,ipt,imovetype))
      end do
      write(uout,'(a,t50,a,6(10x,a))') txmovetype(imovetype), 'total', txpt(1:npt)
      write(uout,'(a,t50,a,6(10x,a))') '-------------------', '-----', ('----------',ipt = 1,npt)
      do ievent = 0, nevent
         if (ny(ievent,0,imovetype) > 0) &
         write(uout,'(a,t40,6(i10,f5.2,5x))') txevent(ievent), (ny(ievent,ipt,imovetype),any(ievent,ipt),ipt = 0,npt)
      end do
      write(uout,'()')
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of displaced part.  = ', ancl1(1,1:npt)
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of net joining part.= ', ancl1(2,1:npt)
   end do

   if (ny(0,0,ispartcl2move) > 0) then
      do ipt = 1, npt
         ancl2(1:2,ipt) = nprimpartcl(1:2,ipt)*InvInt(ny(imcaccept,ipt,ispartcl2move))
      end do
      write(uout,'()')
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of diplacad part.   = ', ancl2(1,1:npt)
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of primary part.    = ', ancl2(2,1:npt)
   end if
   if (ny(0,0,ibrushcl2move) > 0) then
      do ipt = 1, npt
         ancl2(1:2,ipt) = nprimpartcl(1:2,ipt)*InvInt(ny(imcaccept,ipt,ibrushcl2move))
      end do
      write(uout,'()')
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of diplacad part.   = ', ancl2(1,1:npt)
      write(uout,'(a,t60,5(5x,f10.2,5x))') 'average no of primary part.    = ', ancl2(2,1:npt)
   end if

end subroutine MCAverWrite

!........................................................................

subroutine MCAverData
   use MollibModule, only: InvInt
   open(90, file = 'mcaver.data', position = 'append')
   write(90,'(15g10.3)') dtran(1:npt), pspart(1:npt), radcl1(1:npt), &
                         nys1(imcaccept,1:npt,ispartmove)*InvInt(nys1(0,1:npt,ispartmove)), &
                         npcl1s1(1,1:npt,ispartmove)*InvInt(nys1(imcaccept,1:npt,ispartmove))
   close(90)
end subroutine MCAverData

!........................................................................

end subroutine MCAver

!************************************************************************
!*                                                                      *
!*     MCAllAver                                                        *
!*                                                                      *
!************************************************************************

! ... calculate means of mc specific variables for mcall

subroutine MCAllAver(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCAllAver'
   integer(4), save :: nys1(0:nevent), nys2(0:nevent)

   if (slave) return

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      nys1 = 0
      if (txstart == 'continue') read(ucnf) nys1

   case (iBeforeMacrostep)

      nys2 = 0

   case (iSimulationStep)

      nys2(ievent) = nys2(ievent) + 1                                     ! sample

   case (iAfterMacrostep)

      nys2(0) = sum(nys2(1:nevent))                                       ! get sum over ievent
      nys1 = nys1 + nys2

      call MCAllAversub(nys2)

      write(ucnf) nys1

   case (iAfterSimulation)

      call MCAllAversub(nys1)

   end select

contains

!........................................................................

subroutine MCAllAversub(ny)

   use MollibModule, only: InvInt
   integer(4), intent(in) :: ny(0:nevent)
   real(8)    :: any(0:nevent)
   integer(4) :: ievent

   any = ny*InvInt(ny(0))

   call WriteHead(2, 'mc statistics', uout)
   write(uout,'(a,t50,a,6(10x,a))') 'all-particle move', 'total'
   write(uout,'(a,t50,a,6(10x,a))') '-----------------', '-----'
   do ievent = 0, nevent
     if (ny(ievent) > 0) write(uout,'(a,t40,i10,f5.2)') txevent(ievent), ny(ievent), any(ievent)
   end do
end subroutine MCAllAversub

!........................................................................

end subroutine MCAllAver

!************************************************************************
!*                                                                      *
!*     MCWeightIO                                                       *
!*                                                                      *
!************************************************************************

! ... perform i/o for an umbrella sampling with a distance dependent weight
!     umbrella sampling: weighting function dependent on separation between two particles

subroutine MCWeightIO(iStage)

   use MCModule
   implicit none

   integer(4), intent(in)  :: iStage

   character(40), parameter :: txroutine ='MCWeightIO'

   namelist /nmlMCWeight/ ipmcw1, ipmcw2, txpotmcw, npolmcw, acoeffmcw

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      npolmcw = 0
      acoeffmcw = Zero

      rewind(uin)
      read(uin,nmlMCWeight)

      call LowerCase(txpotmcw)

      if (min(ipmcw1,ipmcw2) < 1) call Stop(txroutine, 'error in ipmcw1 and/or ipmcw2', uout)
      if (max(ipmcw1,ipmcw2) > np) call Stop(txroutine, 'error in ipmcw1 and/or ipmcw2', uout)
      if (npolmcw > mnpolmcw) call Stop(txroutine, 'npolmcw > mnpolmcw', uout)

   case (iWriteInput)

      call WriteHead(2, 'MCWeight', uout)
      write(uout, *) 'weighting function is applied on the separation between two particles'
      write(uout,'(a)') 'note!!! at this stage averages are not generally reweighted'
      write(uout,'()')
      write(uout,'(a,t35,i10)') 'partice 1                      = ', ipmcw1
      write(uout,'(a,t35,i10)') 'partice 2                      = ', ipmcw2
      if (txpotmcw == 'polynomial') then
         write(uout, '(a,t35,i10)') 'polynomial weighting function (acoeffmcw(0) + acoeffmcw(1)*r + ...)'
         write(uout,'(a,t35,i10)') 'degree of the polynomial       = ', npolmcw
         write(uout,'(a,t35,6g12.5)') 'coefficients of the polynomial = ', acoeffmcw(0:npolmcw)
      else if (txpotmcw == 'exponential') then
         write(uout, '(a,t35,i10)') 'exponential weighting function (acoeffmcw(0)*exp(-acoeffmcw(1)*(r-acoeffmcw(2))))'
         write(uout,'(a,t35,6g12.5)') 'coefficients                   = ', acoeffmcw(0:2)
      else
         call Stop(txroutine, 'error in txpotmcw', uout)
      end if

   end select

end subroutine MCWeightIO

!************************************************************************
!*                                                                      *
!*     MCWeight                                                         *
!*                                                                      *
!************************************************************************

! ... calculate weight based on trial and current position for mcweight

real(8) function MCWeight()

   use MCModule
   implicit none

   integer(4) :: iploc
   real(8)    :: potold, potnew, r2new, r2old, PolVal

   if ((ipmove /= ipmcw1) .and. (ipmove /= ipmcw2)) then    ! particles ipmcw1 and ipmcw2 are not moved
      MCWeight = One
   else                                                     ! either particle ipmcw1 or ipmcw2 is moved
      iploc = 1
      if (ipmove == ipmcw1) then
         r2new = sum((rotm(1:3,iploc)-ro(1:3,ipmcw2))**2)
      else if (ipmove == ipmcw2) then
         r2new = sum((rotm(1:3,iploc)-ro(1:3,ipmcw1))**2)
      end if
      r2old = sum((ro(1:3,ipmcw1)-ro(1:3,ipmcw2))**2)
      if (txpotmcw == 'polynomial') then
         potnew = PolVal(npolmcw,acoeffmcw,sqrt(r2new))
         potold = PolVal(npolmcw,acoeffmcw,sqrt(r2old))
      else if (txpotmcw == 'exponential') then
         potnew = acoeffmcw(0)*exp(-acoeffmcw(1)*(sqrt(r2new)-acoeffmcw(2)))
         potold = acoeffmcw(0)*exp(-acoeffmcw(1)*(sqrt(r2old)-acoeffmcw(2)))
      end if
      MCWeight = exp(-beta*(potnew-potold))
   end if

end function MCWeight


!************************************************************************
!*                                                                      *
!*     MCWeightInverse                                                  *
!*                                                                      *
!************************************************************************

! ... calculate the inverse weight based on current position for mcweight

real(8) function MCWeightInverse()
   use MCModule
   implicit none
   real(8) :: r2, pot, PolVal
   r2 = (ro(1,ipmcw1)-ro(1,ipmcw2))**2 + (ro(2,ipmcw1)-ro(2,ipmcw2))**2 + (ro(3,ipmcw1)-ro(3,ipmcw2))**2
   if (txpotmcw == 'polynomial') then
      pot = PolVal(npolmcw,acoeffmcw,sqrt(r2))
   else if (txpotmcw == 'exponential') then
      pot = acoeffmcw(0)*exp(-acoeffmcw(1)*(sqrt(r2)-acoeffmcw(2)))
   end if
   MCWeightInverse = exp(beta*pot)
end function MCWeightInverse

!************************************************************************
!*                                                                      *
!*     UmbrellaIO                                                       *
!*                                                                      *
!************************************************************************

! ... perform i/o of the umbrella potential sampling

subroutine UmbrellaIO(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='UmbrellaIO'
   integer(4) :: igrid, l19
   integer(4) :: UmbrellaBin
   real(8) :: dum
   character(15) :: fumbin

   namelist /nmlUmbrella/ typeumb, ipumb1, ipumb2, iaumb1, iaumb2, rminumbrella, delumb, numbgrid, cupdate, umbcoord, lreadumb

   if (.not.lmc .and. .not.lmcall) return      ! only for lmc or lmcall

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

     typeumb = '      '
     ipumb1  = 0
     ipumb2  = 0
     iaumb1  = 0
     iaumb2  = 0
     rminumbrella = 3.0d0
     delumb  = 1.0d0
     numbgrid = 0
     cupdate = '    '
     umbcoord = 'r'
     lreadumb = .false.

     rewind(uin)
     read(uin,nmlUmbrella)

     call LowerCase(typeumb)
     call LowerCase(cupdate)
     call LowerCase(umbcoord)

     if ((umbcoord/= 'r') .and. (umbcoord/= 'x') .and. (umbcoord/= 'y') .and. (umbcoord/= 'z')) then
        write(uout,'(2a)') 'the specified value of umbcoord is not supported: ', umbcoord
        stop 1
     end if

! ... set lradumb; if the umbrella potential coordinate is a distance, we should calculate g(r)

     lradumb =.false.
     if ((typeumb == 'ppdist' .or. typeumb == 'aadist') .and. umbcoord == 'r') lradumb =.true.

   case (iWriteInput)

      if (ipumb1 == 0 .or. ipumb2 == 0) call Stop(txroutine, 'ipumb1 == 0 .or. ipumb2 == 0', uout)

! ... write input data

      call WriteHead(2, 'generating umbrella potential', uout)
      if (typeumb == 'ppdist') then
         write(uout,'(a)') 'umbrella potential calculated for particle-particle distance'
         write(uout,'(a,2i5)') 'use particles labeled          = ', ipumb1, ipumb2
      else if (typeumb == 'aadist') then
         write(uout,'(a)') 'umbrella potential calculated for atom-atom distance'
         write(uout,'(a,2i5)') 'use particles labeled          = ', ipumb1, ipumb2
         write(uout,'(a,2i5)') 'use atoms     labeled          = ', iaumb1, iaumb2
      else
         write(uout,'(a)') 'type of umbrella potential not defined'
         stop 1
      end if
      if (umbcoord /= 'r') then
         write(uout,'(2a)') 'the considered particles are allowed to move along axis: ', umbcoord
         write(uout,'(a)')  '               they are not allowed to rotate'
      end if
      if (lradumb) write(uout,'(a)') 'the probability is scaled to give g(r)'
      write(uout,'(a,f8.4)') 'minimum distance for pmf       = ', rminumbrella
      write(uout,'(a,f8.4)') 'size of each grid for pmf      = ', delumb
      write(uout,'(a,i5)')   'number of grid points          = ', numbgrid
      if (cupdate == 'enga') then
         write(uout,'(a)') 'c is calculated as n**3/max(x)'
      else if (cupdate == 'engb') then
         write(uout,'(a)') 'c is calculated as 0.00004/min(x)**0.27'
      end if

! ... initiate values

      if (lreadumb) then
         l19 = 19
         fumbin = 'fumbin'
         call FileOpen(l19, fumbin, 'form/noread')
         read(l19,'(2f18.8)') (dum,xumb(igrid),igrid = 1,numbgrid)
         close(l19)
         xumbmax = maxval(xumb(1:numbgrid))
         xumbmin = minval(xumb(1:numbgrid))
      else
         xumb(1:numbgrid) = One
         if (typeumb == 'ppdist' .or. typeumb == 'aadist') then
            xumb(1) = 1.0d6
            xumb(numbgrid) = 1.0d6
         end if
         xumbmax = One
         xumbmin = One
      end if
      write(uout,'()')
      write(uout,'(a)') 'initial values for the potential of mean force'
      write(uout,'()')
      write(uout,'(2f18.8)') (rminumbrella+delumb*(real(igrid)-0.5),xumb(igrid),igrid = 1,numbgrid)

   case (iBeforeSimulation)

      iumb = UmbrellaBin(0)

   case (iAfterSimulation)

! ... write potential of mean force

     call WriteHead(2, 'umbrella potential', uout)
     write(uout,'(a)') 'potential of mean force'
     write(uout,'()')
     write(uout,'(2f18.8)') (rminumbrella+delumb*(real(igrid)-0.5),xumb(igrid),igrid = 1,numbgrid)

   end select

end subroutine UmbrellaIO

!************************************************************************
!*                                                                      *
!*     UmbrellaWeight                                                   *
!*                                                                      *
!************************************************************************

! ... calculate weight for umbrella sampling

real(8) function UmbrellaWeight(npmove)

   use MCModule
   implicit none

   integer(4), intent(in) :: npmove   ! parameter passing to UmbrellaBin

   integer(4) :: UmbrellaBin
   real(8)    :: rw, rwnew, wrad, wradnew

   iumbnew = UmbrellaBin(npmove)                             ! get umbrella bin for trial move
   UmbrellaWeight = xumb(iumb)/xumb(iumbnew)                 ! weight from xumb

   if (lradumb) then                                         ! include volume part of weight
     rw        = real(iumb   -1)*delumb + rminumbrella
     rwnew     = real(iumbnew-1)*delumb + rminumbrella
     wrad      = (rw    + delumb)**3 - rw**3
     wradnew   = (rwnew + delumb)**3 - rwnew**3
     UmbrellaWeight = UmbrellaWeight*wrad/wradnew
   end if

end function UmbrellaWeight

!************************************************************************
!*                                                                      *
!*     UmbrellaBin                                                      *
!*                                                                      *
!************************************************************************

! ... calculate bin for umbrella potential

integer(4) function UmbrellaBin(npmove)

   use MCModule
   implicit none

   integer(4), intent(in)  :: npmove        ! =0 implies that ro should be used for both particles

   integer(4) :: iploc
   real(8)    :: r2, rp1(3), rp2(3), dx, dy, dz

! ... several types of coordinates for the umbrella potential

   if (typeumb == 'ppdist') then      ! one or all particle are moved

     if (npmove == 0) then            ! all particles are moved
        rp1(1:3) = ro(1:3,ipumb1)
        rp2(1:3) = ro(1:3,ipumb2)
     else                            ! one particle is moved
        iploc = 1
        if (ipmove == ipumb1) then
           rp1(1:3) = rotm(1:3,iploc)
           rp2(1:3) = ro(1:3,ipumb2)
        else if (ipmove == ipumb2) then
           rp1(1:3) = ro(1:3,ipumb1)
           rp2(1:3) = rotm(1:3,iploc)
        else
           rp1(1:3) = ro(1:3,ipumb1)
           rp2(1:3) = ro(1:3,ipumb2)
        end if
     end if
     dx = rp2(1) - rp1(1)
     dy = rp2(2) - rp1(2)
     dz = rp2(3) - rp1(3)

   else if (typeumb == 'aadist') then      ! works only for lmcall (npmove == 0)

     if (npmove == 0) then
        dx = r(1,ianpn(ipumb2)+iaumb2-1) - r(1,ianpn(ipumb1)+iaumb1-1)
        dy = r(2,ianpn(ipumb2)+iaumb2-1) - r(2,ianpn(ipumb1)+iaumb1-1)
        dz = r(3,ianpn(ipumb2)+iaumb2-1) - r(3,ianpn(ipumb1)+iaumb1-1)
     else
        write(uout, '(a)') 'warning: atom-atom ump only implemented for lmcall'
     end if

   end if

   call PBCr2(dx,dy,dz,r2)
   UmbrellaBin = min(int(numbgrid),max(1,int((sqrt(r2)-rminumbrella)/delumb)+1))

end function UmbrellaBin

!************************************************************************
!*                                                                      *
!*     UmbrellaUpdate                                                   *
!*                                                                      *
!************************************************************************

! ... update weight function for updated umbrella potential

subroutine UmbrellaUpdate

   use MCModule
   implicit none

   real(8) :: cumb

   if (ievent == imcaccept) iumb = iumbnew

   if (cupdate == 'enga') then                             ! calculate cumb
     cumb = real(((istep1-1)*nstep2+istep2)**3)/xumbmax
   else if (cupdate == 'engb') then
     cumb = 0.00004d0/xumbmin**0.27d0
   end if

   xumb(iumb) = xumb(iumb)*(One+cumb)                     ! update xumb
   xumbmax = max(xumbmax,xumb(iumb))                      ! update xumbmax
   xumbmin = minval(xumb(1:numbgrid))                     ! update xumbmin

end subroutine UmbrellaUpdate

!************************************************************************
!*                                                                      *
!*     MCPmfIO                                                          *
!*                                                                      *
!************************************************************************

! ... perform i/o of the mc pmf sampling

subroutine MCPmfIO(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MCPmfIO'
   integer(4) :: ibin
   integer(4) :: MCPmfBin

   namelist /nmlMCPmf/ iptmcpmf, nbinmcpmf, rlowmcpmf, ruppmcpmf, termmcpmf

   if (.not.lmc) return      ! only for lmc

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      iptmcpmf  = 1
      nbinmcpmf = 100
      rlowmcpmf = Zero
      ruppmcpmf = 100.0d0

      rewind(uin)
      read(uin,nmlMCPmf)

      binmcpmf = (ruppmcpmf - rlowmcpmf)/nbinmcpmf
      binimcpmf = One/binmcpmf

   case (iWriteInput)

      if (nppt(iptmcpmf) /= 2) call Stop(txroutine, 'nppt(iptmcpmf) /= 2', uout) ! 2 particles of type iptpmcpmf

      mcpmf(1:nbinmcpmf) = One
      mcpmfmax = One
      mcpmfmin = One

      call WriteHead(2, 'mc potential of mean force', uout)
      write(uout,'(a,i5)')   'number of bins                 = ', nbinmcpmf
      write(uout,'(a,f8.4)') 'lower limit of mc pmf          = ', rlowmcpmf
      write(uout,'(a,f8.4)') 'upper limit of mc pmf          = ', ruppmcpmf
      write(uout,'(a,f8.4)') 'bin of mc pmf                  = ', binmcpmf
      write(uout,'(a,f8.4)') 'termmcpmf                      = ', termmcpmf
      write(uout,'()')
      write(uout,'(a)') 'initial values of the mc pmfi (kT)'
      write(uout,'()')
      write(uout,'(3f18.4)') (rlowmcpmf+binmcpmf*(real(ibin)-0.5), mcpmf(ibin), &
         -log(mcpmf(ibin))+log(mcpmf(nbinmcpmf)), ibin = 1,nbinmcpmf)
   case (iBeforeSimulation)

      ibinmcpmf = MCPmfBin(0)
      if (txstart == 'continue') read(ucnf)  mcpmf(1:nbinmcpmf)

   case (iAfterMacrostep)

      write(ucnf) mcpmf(1:nbinmcpmf)

   case (iAfterSimulation)

      call WriteHead(2, 'mc potential of mean force', uout)
      write(uout,'(a)') 'final potential of mean force (kT)'
      write(uout,'()')
      write(uout,'(3f18.4)') (rlowmcpmf+binmcpmf*(real(ibin)-0.5), mcpmf(ibin), &
         -log(mcpmf(ibin))+log(mcpmf(nbinmcpmf)), ibin = 1,nbinmcpmf)
      write(94,'(f12.3,a,f12.3)') (rlowmcpmf+binmcpmf*(real(ibin)-0.5), tab, &
         -log(mcpmf(ibin))+log(mcpmf(nbinmcpmf)), ibin = 1,nbinmcpmf)
      write(95,'(f12.3,a,f12.3)') (rlowmcpmf+binmcpmf*(real(ibin)-0.5), tab, &
         -log((nbinmcpmf+0.5-ibin)/0.5), ibin = 1,nbinmcpmf)

   end select

end subroutine MCPmfIO

!************************************************************************
!*                                                                      *
!*     MCPmfWeight                                                      *
!*                                                                      *
!************************************************************************

! ... calculate weight for mcpmf

real(8) function MCPmfWeight(npmove)

   use MCModule
   implicit none

   integer(4), intent(in) :: npmove   ! parameter passing to MCPmfBin
   integer(4) :: MCPmfBin

   ibinnewmcpmf = MCPmfBin(npmove)                                ! get mcpmf bin for trial move
   MCPmfWeight = mcpmf(ibinmcpmf)/mcpmf(ibinnewmcpmf)             ! weight from mcpmf

   if (itestmc == 9) &
   write(20,'(a,2i5,f15.5)') 'MCPmfWeight: ibinnewmcpmf, ibinmcpmf, MCPmfWeight', &
    ibinnewmcpmf, ibinmcpmf, MCPmfWeight
end function MCPmfWeight

!************************************************************************
!*                                                                      *
!*     MCPmfBin                                                         *
!*                                                                      *
!************************************************************************

! ... calculate bin for umbrella potential

integer(4) function MCPmfBin(npmove)

   use MCModule
   implicit none

   integer(4), intent(in)  :: npmove        ! =0 implies that ro should be used for both particles

   integer(4) :: ip1, ip2, iploc
   real(8)    :: r2, rp1(3), rp2(3), dx, dy, dz

   ip1 = ipnpt(iptmcpmf)
   ip2 = ipnpt(iptmcpmf)+1

   if (npmove == 0) then            ! all particles are moved
      rp1(1:3) = ro(1:3,ip1)
      rp2(1:3) = ro(1:3,ip2)
   else                            ! one particle is moved
      iploc = 1
      if (ipmove == ip1) then
         rp1(1:3) = rotm(1:3,iploc)
         rp2(1:3) = ro(1:3,ip2)
      else if (ipmove == ip2) then
         rp1(1:3) = ro(1:3,ip1)
         rp2(1:3) = rotm(1:3,iploc)
      else
         rp1(1:3) = ro(1:3,ip1)
         rp2(1:3) = ro(1:3,ip2)
      end if
   end if
   dx = rp2(1) - rp1(1)
   dy = rp2(2) - rp1(2)
   dz = rp2(3) - rp1(3)

   call PBCr2(dx,dy,dz,r2)
   MCPmfBin = min(int(nbinmcpmf),max(1,int(binimcpmf*(sqrt(r2)-rlowmcpmf))+1))

   if (itestmc == 9) write(20,'(a,i5,f10.3,i5)') 'MCPmfBin: ipmove, r1, MCPmfBin', ipmove, sqrt(r2), MCPmfBin

end function MCPmfBin

!************************************************************************
!*                                                                      *
!*     MCPmfUpdate                                                      *
!*                                                                      *
!************************************************************************

! ... update weight function for mc pmf

subroutine MCPmfUpdate

   use MCModule
   implicit none

   real(8) :: term

   if (iptmove /= iptmcpmf) return

   if (ievent == imcaccept) ibinmcpmf = ibinnewmcpmf

!  term = real((istep2+(istep1-1)*nstep2)**3)/mcpmfmax
!  term = 0.00004d0/mcpmfmin**0.27d0
   term = log(real(istep2+(istep1-1)*nstep2))/log(real(nstep1*nstep2))
   term = termmcpmf*(1-0.99*term)

   mcpmf(ibinmcpmf) = mcpmf(ibinmcpmf)*(One+term)        ! update mcpmf
   mcpmfmax = max(mcpmfmax, mcpmf(ibinmcpmf))            ! update mcpmfmax
   mcpmfmin = minval(mcpmf(1:nbinmcpmf))                 ! update mcpmfmin

   if (itestmc == 9) &
    write(20,'(a,i10,f15.5)') 'MCPmfUpdate: ibinmcpmf, mcpmf(ibinmcpmf)',ibinmcpmf, mcpmf(ibinmcpmf)
!  write(33,*) ipass+(istep2-1)*np, ibinmcpmf
end subroutine MCPmfUpdate

!************************************************************************
!*                                                                      *
!*     TestMCMove                                                       *
!*                                                                      *
!************************************************************************

! ... test output for MC move

subroutine TestMCMove(unit)

   use MCModule
   implicit none

   character(40), parameter :: txroutine ='TestMCMove'
   integer(4),   intent(in) :: unit
   integer(4) :: iploc, ip, ipt, ialoc, ia

   call WriteHead(3, txroutine, unit)
   write(unit,'(a,i5)') 'ipass        =', ipass
   write(unit,'(a,a)')  'txmovetype   =', txmovetype(imovetype)
   write(unit,'(a,i5)') 'ipmove       =', ipmove
   write(unit,'(a,i5)') 'iptmove      =', iptmove
   write(unit,'(a,i5)') 'nptm         =', nptm

   write(unit,'()')
   write(unit,'(a,t20,a,t40,a,t65,a,t95,a)') 'local particle no', 'particle no', 'particle type', 'old com', 'new com'
   write(unit,'(a,t20,a,t40,a,t65,a,t95,a)') '-----------------', '-----------', '-------------', '-------', '-------'
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      write(unit,'(2(i10,5x),10x,a,2(5x,3F8.3))') iploc, ip, txpt(iptpn(ip)), ro(1:3,ip), rotm(1:3,iploc)
   end do

   ialoc = 0
   write(unit,'()')
   write(unit,'(a,t20,a,t40,a,t65,a,t95,a)') 'local atom no', 'atom no', 'atom type', 'old com', 'new com'
   write(unit,'(a,t20,a,t40,a,t65,a,t95,a)') '-------------', '-------', '---------', '-------', '-------'
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do ia = ianpn(ip), ianpn(ip) + napt(ipt) - 1
         ialoc = ialoc + 1
         write(unit,'(2(i10,5x),10x,a,2(5x,3F8.3))') ialoc, ia, txat(iatan(ia)), r(1:3,ia), rtm(1:3,ialoc)
      end do
   end do

end subroutine TestMCMove
