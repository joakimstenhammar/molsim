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



!> \page nmlPotential
!! The namelist  \ref nmlPotential contains variables that describe the potentials between the particles.  Tabulated potential energy and
!! force evaluation, which often gives a fast code by avoiding square root and divisions. Procedure according to Andrea, Swope, and
!! Andersen, JCP 79, 4576 (1983).
!! * Variables:
!!  * \subpage r2uminin
!!  * \subpage utoltab
!!  * \subpage ftoltab
!!  * \subpage umaxtab
!!  * \subpage fmaxtab
!!  * \subpage rcut
!!  * \subpage txpot
!!  * \subpage npot
!!  * \subpage ipot
!!  * \subpage ucoff
!!  * \subpage relpermitt
!!  * \subpage lscrc
!!  * \subpage scrlen
!!  * \subpage lscrhs
!!  * \subpage lewald
!!  * \subpage txewaldrec
!!  * \subpage iewaldopt
!!  * \subpage uewaldtol
!!  * \subpage ualpha
!!  * \subpage ualphared
!!  * \subpage ncut
!!  * \subpage ncutregion
!!  * \subpage lsurf
!!  * \subpage lewald2dlc
!!  * \subpage order
!!  * \subpage nmesh
!!  * \subpage lrf
!!  * \subpage epsrf
!!  * \subpage epsimg
!!  * \subpage radimg
!!  * \subpage epsi1
!!  * \subpage epsi2
!!  * \subpage boundaryrad
!!  * \subpage lmaxdiel
!!  * \subpage lbg
!!  * \subpage lljcut
!!  * \subpage ljrcut
!!  * \subpage ljushift
!!  * \subpage lambda_sw
!!  * \subpage epsilon_sw
!!  * \subpage alpha_sw
!!  * \subpage lambda_ramp
!!  * \subpage epsilon_ramp
!!  * \subpage alpha_ramp
!!  * \subpage rad_dep
!!  * \subpage rho_dep
!!  * \subpage factor_dep
!!  * \subpage lellipsoid
!!  * \subpage radellipsoid
!!  * \subpage aellipsoid
!!  * \subpage lsuperball
!!  * \subpage radsuperball
!!  * \subpage qsuperball
!!  * \subpage txmethodsuperball
!!  * \subpage nitersuperball
!!  * \subpage tolsuperball
!!  * \subpage meshdepthsuperball
!!  * \subpage dl_damp
!!  * \subpage dl_cut
!!  * \subpage dr_damp
!!  * \subpage dr_cut
!!  * \subpage lstatsuperball
!!  * \subpage luext
!!  * \subpage lmonoatom
!!  * \subpage itestpot

!************************************************************************
!> \page potential potential.F90
!! **PotentialModule**
!! *module for potential*
!************************************************************************
module PotentialModule

   use MolModule

   integer(4), parameter :: mninpot  = 12  ! largest exponent of an 1/r**(-exponent) energy term
!> \page txpot
!! `character(20)`
!! * Text label used for selecting potential of each pair of particle type. The potentials may either be already existing in the program or supplied by the user. The available options are:
!! * '`(1,6,12)`': Charge plus Lennard-Jones interaction \f$ u_{ij}(r) = q_iq_j/4\pi\epsilon_0 r + 4 \epsilon_{ij} ((
!!    \sigma_{ij}/r)^{12}( \sigma_{ij} /r)^6) \f$. The parameters \f$ q_i \f$, \f$ \sigma_{ii} \f$, and \f$ \epsilon_{ii}\f$ are given in
!!     namelist \ref nmlParticle. Lorentz-Berthelot mixing rules are applied for cross terms.
!! * '`setx`': Spherical Ewald potential. As the default potential, but the 1/r term is multiplied with erfc(r*\ref ualphared/\ref rcut).
!! * '`mcy`': The MCY water potential, ref. JCP 64, 1351 (1976).
!! * '`nemo:xxx`': The two-body part of the Nemo potential 'xxx'. The potential form and the coefficients are read from file
!!    molsim.lib. Present installed potentials include nemo:ww (water-water) and nemo:uw(urea-water). Note, the specification of the
!!    many-body polarization interaction is given in namelist \ref nmlParticle.
!! * '`sw`': Square-well potential. Parameters: \ref lambda_sw, \ref epsilon_sw and \ref alpha_sw
!! * '`ramp`': Ramp potential. Parameters
!! * '`asakura-oosawa`': Asukura-Oosawa potential
!! * '`xxx`': Search for user-provided potential labeled 'xxx' called from routine PotentialUser in file moluser.F90.
!! * If there is no match, the default potential form sum[\ref ucoff (1:\ref npot)/r**\ref ipot(1:\ref npot)] and the variables \ref npot, \ref ipot, and \ref ucoff,
!!   which are read in this namelist, are used. If \ref zat /= 0, the Coulomb term need not to be specified.
   character(20), allocatable :: txpot(:)  ! label of two-body potential
!> \page npot
!! `integer`(1:natat)
!! * \ref npot (iatjat) denotes the number of terms of atom type pair iatjat.
   integer(4), allocatable :: npot(:)
   integer(4)    :: npotm                  ! maximum value of npot
!> \page ipot
!! `integer`(1:\ref npot ,1:natat)
!! * \ref ipot (m,iatjat) denotes the exponent of term m of atom type pair iatjat.
   integer(4), allocatable :: ipot(:,:)
   integer(4)    :: ipotm                  ! maximum value of ipot
!> \page ucoff
!! `real`(1:\ref npot ,1:natat)
!! * \ref ucoff (m,iatjat) denotes the coefficient of term m of atom type pair iatjat.
   real(8), allocatable    :: ucoff(:,:)
   real(8), allocatable,save :: ucoffx(:,:)! coefficient of a term
!> \page relpermitt
!! `real`
!! **default:** `1.0`
!! * Relative permittivity.
   real(8)       :: relpermitt
!> \page r2uminin
!! `real`
!! **deault:** `0.1`
!! * Square of the lower end of the tabulated potential.
   real(8)       :: r2uminin
!> \page utoltab
!! `real`
!! **default:** `10d-5`
!! * Energy tolerance of the tabulated potential.
   real(8)       :: utoltab
!> \page ftoltab
!! `real`
!! **default:** `10d-5`
!! * Force tolerance of the tabulated potential.
   real(8)       :: ftoltab
!> \page umaxtab
!! `real`
!! **default:** `2*10d4`
!! * Energy at which the table is cut off at small separation.
   real(8)       :: umaxtab
!> \page fmaxtab
!! `real`
!! **default:** `2*10d4`
!! * Force at which the table is cut off at small separation.
   real(8)       :: fmaxtab                ! maximum value of force
   integer(4),allocatable,save :: nugrid(:)! number of grid points for iatjat
   real(8), allocatable, save :: rumin(:)  ! lower limit of potential table
   real(8), allocatable, save :: rumax(:)  ! upper limit of potential table
   real(8), allocatable, save :: r2umax(:) ! upper limit of potential table squared
   real(8)       :: rumaxfac               ! factor for extending rumax to avoid round-off problems
   logical       :: lmorememory            ! flag for increasing memory
   logical, allocatable :: lsetatat(:)     ! for termorary use
!> \page itestpot
!! `integer`
!! **default:** `0`
!! * Flag for test output. This possibility is for maintenance purposes.
!! * `0`: Nothing. The normal option.
!! * '1': Intermediate potential variables are written
!! * `2`: Examination of accuracy of two-body potentials
!! * `3`: Plot of two-body potentials
!! * `4`: Examination of truncation error of Ewald summation (routine PotTwoBodyTab).
   integer(4)    :: itestpot
!> \page itestpotchain
!! `integer`
!! **default:** `0`
!! * Flag for test output. This possibility is for maintenance purposes.
!! * `0`: Nothing. The normal option.
!! * `1`: Bond length table and bond angle table (routines BondLengthTab and BondAngleTab).
   integer(4)    :: itestpotchain
end module PotentialModule

! relation among potential routines

!     PotentialDriver
!       !
!       !
!       !-------------- IOPotTwoBody
!       !
!       !-------------- PotTwoBodyTab1
!       !                      !
!       !                      ! -----  PotTwoBodyTab2
!       !                      !              !
!       !                                     ! -----  PotTwoBodyTab3
!       !                                     !              !
!       !   lchain                                           !--- potsub
!       !-------------- IOPotChain                                   !
!       !   lchain                                                   !--- SetUBuffer
!       !-------------- ChainTab                                     !
!       !   luext                                                    !--- CheckUBuffer
!       !-------------- IOPotExternal
!       !   lpolarization
!       !-------------- IOPolarizationIter
!       !

!************************************************************************
!> \page potential potential.F90
!! **PotentialDriver**
!! *driver of potential routines*
!************************************************************************


subroutine PotentialDriver(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='PotentialDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

                   call IOPotTwoBody(iStage)
      if (lchain)  call IOPotChain(iStage)
      if (luext)   call IOPotExternal(iStage)
      if (lpolarization)  call IOPolarizationIter(iStage)
      if (lstatsuperball) call SuperballStat(iStage)

   case (iWriteInput)

                   call IOPotTwoBody(iStage)
                   call PotTwoBodyTab1(.true.)
      if (lchain)  call IOPotChain(iStage)
      if (lchain)  call ChainTab
      if (luext)   call IOPotExternal(iStage)
      if (lpolarization)  call IOPolarizationIter(iStage)
      if (lstatsuperball) call SuperballStat(iStage)

   case (iBeforeSimulation)

      if (itestpot == 2) call PlotPotTwoBodyTab
      if (master .and. itestpot == 3) call TestPotTwoBodyTab(iStage,uout)
      if (itestpot == 4) call TestEwald(iStage)
      if (lstatsuperball) call SuperballStat(iStage)

   case (iAfterMacrostep)

      if (itestpot == 4) call TestEwald(iStage)
      if (lstatsuperball) call SuperballStat(iStage)

   case (iAfterSimulation)

      if (itestpot == 4) call TestEwald(iStage)
      if (lstatsuperball) call SuperballStat(iStage)

   end select

contains

!........................................................................

subroutine SuperballStat(iStage)
   integer(4), intent(in) :: iStage
   if (txmethodsuperball == 'nr') call SuperballStatNR(iStage, iaux)
   if (txmethodsuperball == 'mesh') call SuperballStatMesh(iStage, vaux(1:3,1), laux, raux)
end subroutine SuperballStat

!........................................................................

end subroutine PotentialDriver

!************************************************************************
!> \page potential potential.F90
!! **IOPotTwoBody**
!! *perform i/o on potential variables for two body interactions*
!************************************************************************


subroutine IOPotTwoBody(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOPotTwoBody'
   character(80), parameter :: txheading ='two-body potential data'
   integer(4) :: iptjpt, iat, jat, iatjat, m, Getnkvec, Getnkvec2d
   real(8)    :: ucoffaim, fac1, fac2, NetCharge, raddiag2
   logical, allocatable    :: lucoffmod(:)
   real(8)    :: EwaldErrorUReal, EwaldErrorURec
   namelist /nmlPotential/ r2uminin, utoltab, ftoltab, umaxtab, fmaxtab,            &
                           rcut, txpot, npot, ipot, ucoff, relpermitt,              &
                           lscrc, scrlen, lscrhs,                                   &
                           lewald, txewaldrec, iewaldopt, uewaldtol, ualpha, ualphared, ncut, ncutregion, lsurf, lewald2dlc,  &
                           order, nmesh,                                            &
                           lrf, epsrf,                                              &
                           epsimg, radimg,                                          &
                           epsi1, epsi2, boundaryrad, lmaxdiel, lbg,                &
                           lljcut, ljrcut, ljushift,                                &
                           lambda_ramp, epsilon_ramp, alpha_ramp,                   &
                           lambda_sw, epsilon_sw, alpha_sw,                         &
                           rad_dep, rho_dep, factor_dep,                            &
                           lellipsoid, radellipsoid, aellipsoid,                    &
                           lsuperball, radsuperball, qsuperball, txmethodsuperball, nitersuperball, tolsuperball, meshdepthsuperball, &
                           dl_damp, dl_cut, dr_damp, dr_cut, lstatsuperball,        &
                           luext,                                                   &
                           lmonoatom, itestpot

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(txpot)) then
         allocate(txpot(nptpt), npot(natat), ipot(mninpot,natat), ucoff(mninpot,natat))
         txpot = ""
         npot = 0
         ipot = 0
         ucoff = 0.0E+00
      end if

! ... set initial values

      r2uminin       = 0.1d0
      utoltab        = 1e-5
      ftoltab        = 1e-5
      umaxtab        = 20000.0
      fmaxtab        = 20000.0
      rcut           = Zero
      txpot          = 'default'
      ucoff          = Zero
      relpermitt     = One
      lscrc          =.false.
      lscrhs         =.false.
      lewald         =.false.
      txewaldrec     ='std'
      iewaldopt      = -1
      uewaldtol      = Zero
      ualphared      = 3.0d0
      ncutregion     ='sphere'
      lsurf          =.true.
      lewald2dlc     =.false.
      order          = 5
      nmesh          = 48
      lrf            =.false.
      epsrf          = 0.0       ! implies infinite
   !   lmaxdiel = 10
   !   lbg = .false.
      lljcut         =.false.
      ljrcut         = Two**(1.0/6.0)
      ljushift       = One
      lambda_ramp    = 1.1d0     ! end of the ramp potential (divided by hs diamter) (length)
      epsilon_ramp   = 1.0d0     ! depth of the ramp potential (energy)
      alpha_ramp     = 1.0d-3    ! steepness of the rise (lower alpha => steeper rise)
      lambda_sw      = 1.1d0     ! end of the sw potential (divided by hs diamter) (length)
      epsilon_sw     = 1.0d0     ! depth of the sw potential (energy)
      alpha_sw       = 1.0d-3    ! steepness of the rise (lower alpha => steeper rise)
      factor_dep     = 1.0d0     ! depletion-thickness factor
      lellipsoid     = .false.   !
      radellipsoid   = One       ! radius of degenerated axes
      aellipsoid     = One       ! aspect ratio
      lsuperball     = .false.   !
      radsuperball   = One       ! radius of superball
      qsuperball     = One       ! power (=1 => sphere, =inf => cub)
      txmethodsuperball = 'nr'
      nitersuperball = 10
      tolsuperball   = 1.0d-4
      meshdepthsuperball = 4
      dl_cut         = 1d10
      dl_damp        = One
      dr_cut         = 1d10
      dr_damp        = One
      lstatsuperball = .false.
      luext          =.false.

      itestpot       = 0

! ... consistence with earlier versions

      radimg = Zero

! ... read input data

      rewind(uin)
      read(uin,nmlPotential)

      do iptjpt = 1, nptpt
         call LowerCase(txpot(iptjpt))
      end do
      call LowerCase(ncutregion)
      call LowerCase(txewaldrec)

! ... set rumaxfac

      rumaxfac = 1.1d0

! ... set EpsiFourPi

      EpsiFourPi = Epsi0FourPi/relpermitt

! ... set rcut and rcut2

      if (rcut <= Zero) then
         if (lbcbox) then
            if (txbc == 'xyz') then
               rcut = sqrt(boxlen2(1)**2+boxlen2(2)**2+boxlen2(3)**2)
            else if (txbc == 'xy') then
               rcut = sqrt(boxlen2(1)**2+boxlen2(2)**2+boxlen(3)**2)
            else if (txbc == 'z') then
               rcut = sqrt(boxlen(1)**2+boxlen(2)**2+boxlen2(3)**2)
            end if
         else if (lbcrd) then
            rcut = sqrt(Two/Three)*cellside
         else if (lbcto) then
            rcut = sqrt(Three/Two)*cellside
         else if (lbcsph) then
            rcut = Two*sphrad
         else if (lbccyl) then
            rcut = sqrt((Two*cylrad)**2+cyllen**2)
         else if (lbcell) then
            rcut = Two*maxval(ellrad(1:3))
         end if
      end if
      rcut2 = rcut**2

      cyllen2 = half*cyllen

! ... prepare for ewald summation (rcut might be changed)

      if (lewald) call EwaldSetup
!     rcut2 = rcut**2

! ... set rffac and rffac2

      if (lrf) then
         if (epsrf == Zero) then
            rffac = One/(rcut**3)
         else
            rffac = Two*(epsrf-One)/((Two*epsrf+One)*rcut**3)
         end if
         rffac2 = Half*rffac
      end if

! ... prepare for superball

      if (lsuperball) then
         qsuperball2 = Two*qsuperball
         qsuperballi = one/qsuperball
         radsuperball2 = radsuperball**2
         radsuperball2i = one/radsuperball2
         raddiag2 = three**(1-qsuperballi)*radsuperball2  ! radius**2 in (1,1,1)-direction
         rcut2superball(1) = four*radsuperball2           ! dist**2 interval where overlap check is non-trivial
         rcut2superball(2) = four*raddiag2
         call HeapSort(2,rcut2superball)
         if (qsuperball > qsuperball_max_nr) txmethodsuperball = 'mesh'  ! nr algorithm is not stable
         superBallMesh = buildSuperBall(radsuperball,qsuperball2,meshdepthsuperball)

!         volfac = 1.0d0/(TwoPi*(qsuperball)**2)*exp(3*gammaln(half/qsuperball)-gammaln(1.5d0/qsuperball))
!         lenfac = volfac**(Third)
!         boxlen = boxlen*lenfac ! change boxlen to conserve volume fraction of superballs (input: boxlen indep of qsuperball)
!         call SetBoxParam
      end if

   case (iWriteInput)

! ... check conditions

      if (ldipole) then
         if (relpermitt > One) call Stop(txroutine, 'ldipole .and. relpermitt > One not allowed', uout)
         if (lewald2dlc) call Stop(txroutine, 'ldipole .and. lewald2dlc not allowed', uout)
      end if

      if (lpolarization) then
         if (relpermitt > One) call Stop(txroutine, 'lpolarization .and. relpermitt > One not allowed', uout)
      end if

      if (ldipolesph) then
         if (relpermitt > One) call Stop(txroutine, 'ldipolesph and relpermitt > One not allowed', uout)
         if (txbc /= 'sph') call Stop(txroutine, 'ldipolesph and txbc /= sph not allowed', uout)
         if (lscrc .or. lrf .or. lewald) call Stop(txroutine, 'ldipolesph and screened Coulomb / RF / Ewald not supported', uout)
         if (count(napt(1:npt) > 1) > 0) call Stop(txroutine,'ldipolesph .and. lpolyatom',uout)
      end if
      if (laimage .and. (sphrad > radimg)) call Stop(txroutine, 'sphrad > radimg', uout)

      if (lscrc  .and. lewald) call Stop(txroutine, 'both screened coulomb and ewald summation', uout)
      if (lscrc  .and. lrf   ) call Stop(txroutine, 'both screened coulomb and reaction field', uout)
      if (lewald .and. lrf   ) call Stop(txroutine, 'both ewald summation and reaction field', uout)

      if (lewald) then
         if (.not.lPBC)  call Stop(txroutine, 'lewald and non-periodic boundary conditions', uout)  !Check conditions (Jocke)
!123     if (lbd) call Stop(txroutine,'lbd with ewald summation is not implemented',uout)
         do iptjpt = 1, nptpt
            if ((txpot(iptjpt)(1:3) == 'mcy') .or.        &
               (txpot(iptjpt)(1:3) == 'set')) then
               call Stop(txroutine, 'unsupported potential for ewald summation', uout)
            end if
         end do
      end if

      if (lrf) then
         if (lmc) call Stop(txroutine, 'both reaction field and lmc', uout)
         if (lbd) call Stop(txroutine, 'both reaction field and lbd', uout)
         do iptjpt = 1, nptpt
            if ((txpot(iptjpt)(1:3) == 'mcy') .or.        &
               (txpot(iptjpt)(1:3) == 'set')) then
               call Stop(txroutine, 'unsupported potential for reaction field', uout)
            end if
         end do
      end if

      if (lljcut) then
         if (count(txpot=='(1,6,12)') == 0) lljcut =.false.  ! cutoff of lj potential only for (1,6,12) potential
      end if

      if (lellipsoid) radellipsoid2 = radellipsoid**2

! ... check that ucoff(1,:) is consistent to zat(:)

      if(.not.allocated(lucoffmod)) then
         allocate(lucoffmod(natat))
         lucoffmod = .false.
      end if

      lucoffmod =.false.
      if (count(zat(1:nat) /= Zero) > 0) then
         do iat = 1, nat
            do jat = iat, nat
               iatjat = iatat(iat,jat)

               ucoffaim = EpsiFourPi*zat(iat)*zat(jat)     ! calculate ucoffaim

               if (lscrc .and. lscrhs) then                 ! check if modify ucoffaim
                  fac1 = radat(iat)/scrlen
                  fac2 = radat(jat)/scrlen
                  ucoffaim = ucoffaim * exp(fac1)/(1+fac1) * exp(fac2)/(1+fac2)
               end if

               if (ipot(1,iatjat) == 1) then                ! modify coefficient of the 1/r term
                  if (abs(ucoff(1,iatjat)-ucoffaim)/max(abs(ucoff(1,iatjat)),abs(ucoffaim)) > 1d-6) then
                     ucoff(1,iatjat) = ucoffaim
                     lucoffmod(iatjat) =.true.
                  end if
               else                                        ! add a 1/r term
                  do m = npot(iatjat), 1, -1
                     ipot(m+1,iatjat) = ipot(m,iatjat)
                     ucoff(m+1,iatjat) = ucoff(m,iatjat)
                  end do
                  npot(iatjat) = npot(iatjat)+1
                  ipot(1,iatjat) = 1
                  ucoff(1,iatjat) = ucoffaim
                  lucoffmod(iatjat) =.true.
               end if
            end do
         end do
      end if

! ... allocate memory for ucoffx

      npotm = maxval(npot)
      ipotm = 1
      do iatjat = 1, natat
         if (npot(iatjat) > 0) ipotm = max(ipotm,ipot(npot(iatjat),iatjat))
      end do

      if (npotm > mninpot) call Stop(txroutine, 'npotm > mninpot', uout) ! largest number of terms
      if (ipotm > mninpot) call Stop(txroutine, 'ipotm > mninpot', uout) ! largest exponent

      if(.not.allocated(ucoffx)) then
         allocate(ucoffx(ipotm,natat))
         ucoffx = 0.0E+00
      end if

! ... set ucoffx from ucoff

      ucoffx = Zero
      do iat = 1, nat
         do jat = iat, nat
            iatjat = iatat(iat,jat)
            do m = 1, npot(iatjat)
               ucoffx(ipot(m,iatjat),iatjat) = ucoff(m,iatjat)
            end do
         end do
      end do

! ... write input data

      if (master) then

         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f10.3)') 'rcut                           = ', rcut
         write(uout,'(a,t35,g10.3)') 'utoltab                        = ', utoltab
         write(uout,'(a,t35,g10.3)') 'ftoltab                        = ', ftoltab
         write(uout,'(a,t35,g10.3)') 'umaxtab                        = ', umaxtab
         write(uout,'(a,t35,g10.3)') 'fmaxtab                        = ', fmaxtab

         write(uout,'()')
         if (lmonoatom) write(uout,'(a)') 'monoatomic energy and force routines'
         if (lpolyatom) write(uout,'(a)') 'polyatomic energy and force routines'

         if (relpermitt > One) then
            write(uout,'()')
            write(uout,'(a,t35,f10.6)') 'relative permittivity          = ', relpermitt
         end if

         if (maxval(abs(zat(1:nat))) > Zero) then
            write(uout,'()')
            NetCharge = sum(zat(1:nat)*naat(1:nat)*nppt(iptat(1:nat)))
            write(uout,'(a,t35,g10.3)') 'system net charge based on zat =',NetCharge
            if (abs(NetCharge)>1d-10) call Warn('potential' ,'non-zero net charge', uout)
         end if

         if (lewald) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the ewald summation'
            write(uout,'(a)') '------------------------------------------'
            write(uout,'(a,t35,g15.5 )')    'txewaldrec                     = ', txewaldrec
            if (txewaldrec == 'std') then
               write(uout,'(a,t35,3g15.5)') 'iewaldopt                      = ', iewaldopt
               if (iewaldopt > 0) write(uout,'(a,t35,3g15.5)')'uewaldtol                      = ', uewaldtol
               write(uout,'(a,t35,3g15.5)') 'error estimate in real space   = ', EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi)
               write(uout,'(a,t35,3g15.5)') 'error estimate in rec. space   = ', EwaldErrorURec(lq2sum,q2sum,boxlenshort,ualpha,ncut,EpsiFourPi)
               write(uout,'(a,t35,3g15.5)') 'ualpha                         = ', ualpha
               write(uout,'(a,t35,3g15.5)') 'ualphared                      = ', ualphared
               write(uout,'(a,t35,3g15.5)') 'ncut                           = ', ncut
               write(uout,'(a,t35,3g15.5)') 'ncut2                          = ', ncut2
               write(uout,'(a,t35,3g15.5)') 'ncutregion                     = ', ncutregion
               write(uout,'(a,t35,3g15.5)') 'nkvec (one octant)             = ', Getnkvec()
               write(uout,'(a,t35,3g15.5)') 'surface term                   = ', lsurf
               if (lewald2dlc) then
                  write(uout,'()')
                  write(uout,'(a,t35,3g15.5)') '2d layer correction included'
                  write(uout,'(a,t35,3g15.5)') 'ncutd2,                        = ', ncut2d
                  write(uout,'(a,t35,3g15.5)') 'ncut2d2                        = ', ncut2d2
                  write(uout,'(a,t35,3g15.5)') 'nkvec2d (one octant)           = ', Getnkvec2d()
               end if
            else if (txewaldrec == 'spm') then
               write(uout,'(a,t35,3g15.5)') 'ualpha                         = ', ualpha
               write(uout,'(a,t35,3g15.5)') 'ualphared                      = ', ualphared
               write(uout,'(a,t35,3g15.5)') 'error estimate in real space   = ', EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi)
               write(uout,'(a,t35,3g15.5)') 'order                          = ', order
               write(uout,'(a,t35,3g15.5)') 'nmesh                          = ', nmesh
               write(uout,'(a,t35,3g15.5)') 'surface term                   = ', lsurf
            end if
         end if

         if (lrf) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the reaction field method'
            write(uout,'(a)') '------------------------------------------------'
            if (epsrf /= Zero) then
            write(uout,'(a,t35,f15.5)')  'epsrf                          = ', epsrf
            else
            write(uout,'(a,t35,f15.5)')  'epsrf                          =  infinite'
            end if
         end if

         if (laimage) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the image charge approximation'
            write(uout,'(a)') '-----------------------------------------------------'
            write(uout,'(a,t35,f15.5)') 'epsimg                          = ', epsimg
            write(uout,'(a,t35,f15.5)') 'radimg                          = ', radimg
         end if

         if (ldieldis) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the dielectric discontinuity'
            write(uout,'(a)') '---------------------------------------------------'
            write(uout,'(a,t35,f15.5)') 'epsi1                           = ', epsi1
            write(uout,'(a,t35,f15.5)') 'epsi2                           = ', epsi2
            if(txbc == 'sph') then
               write(uout,'(a,t35,f15.5)') 'boundaryrad                     = ', boundaryrad
               write(uout,'(a,t35,3i15)')  'lmaxdiel                        = ', lmaxdiel
               write(uout,'(a,t35,l    )') 'lbg (neutralizing background)   = ', lbg
            end if
         end if

         if (lljcut) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the truncation of lj-potential'
            write(uout,'(a)') '--------------------------------_--------------------'
            write(uout,'(a,t35,g10.3)') 'ljrcut   (in sigma unit)       = ', ljrcut
            write(uout,'(a,t35,g10.3)') 'ljushift (in epsilon unit)     = ', ljushift
         end if

         if (lscrc) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the screened coulomb potential'
            write(uout,'(a)') '-----------------------------------------------------'
            write(uout,'(a,t35,3g15.5)') 'scrlen                         = ', scrlen
            write(uout,'(a,t35,3g15.5)') 'hard-core flag (lscrhs)        = ', lscrhs
         end if

         if (count(txpot=='setx') + count(txpot=='setx') + count(txpot=='setx') > 0) then
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the spherical ewald potential'
            write(uout,'(a)') '----------------------------------------------------'
            write(uout,'(a,t35,3g10.3)') 'ualphared                      = ', ualphared
         end if

         if (count(txpot=='ramp') > 0) then
            if (lambda_ramp < One) call Stop(txroutine//': ramp', 'lambda_ramp < One', uout)
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the ramp potential'
            write(uout,'(a)') '-----------------------------------------'
            write(uout,'(a,t35,3g10.3)') 'lambda_ramp                    = ', lambda_ramp
            write(uout,'(a,t35,3g10.3)') 'epsilon_ramp                   = ', epsilon_ramp
            write(uout,'(a,t35,3g10.3)') 'alpha_ramp                     = ', alpha_ramp
         end if

         if (count(txpot=='sw') > 0) then
            if (lambda_sw < One) call Stop(txroutine//': sw', 'lambda_sw < One', uout)
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the sw potential'
            write(uout,'(a)') '------------------------------------------------'
            write(uout,'(a,t35,3g10.3)') 'lambda_sw                      = ', lambda_sw
            write(uout,'(a,t35,3g10.3)') 'epsilon_sw                     = ', epsilon_sw
            write(uout,'(a,t35,3g10.3)') 'alpha_sw                       = ', alpha_sw
         end if

         if (count(txpot=='asakura-oosawa') > 0) then
            if (npt > 1) call Stop(txroutine//': asakura-oosawa depletion potential', 'npt  > 1', uout)
            if (rho_dep < Zero) call Stop(txroutine//': asakura-oosawa depletion potential', 'rho_dep < Zero', uout)
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the asakura-oosawa potential'
            write(uout,'(a)') '---------------------------------------------------'
            write(uout,'(a,t35,3g10.3)') 'radius of penetrable hs        = ', rad_dep
            write(uout,'(a,t35,3g10.3)') 'number density of penetrable hs= ', rho_dep
         end if

         if (count(txpot=='lekkerkerker-tuinier') > 0) then
            if (npt > 1) call Stop(txroutine//': lekerkerker-tuinier depletion potential', 'npt  > 1', uout)
            if (rho_dep < Zero) call Stop(txroutine//': lekerkerker-tuinier depletion potential', 'rho_dep < Zero', uout)
            write(uout,'()')
            write(uout,'(a)') 'parameters controlling the lekkerkerker-tuinier potential'
            write(uout,'(a)') '---------------------------------------------------------'
            write(uout,'(a,t35,3g10.3)') 'radius of penetrable hs        = ', rad_dep
            write(uout,'(a,t35,3g10.3)') 'number density of penetrable hs= ', rho_dep
            write(uout,'(a,t35,3g10.3)') 'depletion-thickness factor     = ', factor_dep
            if (txuser == 'samia') call WritePotLekkerkerkerTuinier
         end if

         if (lellipsoid) then
            write(uout,'()')
            write(uout,'(a)') 'particles are ellipsoids'
            write(uout,'(a)') '------------------------'
            write(uout,'(a,t35,3g15.3)') 'radius of degerated axes       = ', radellipsoid
            write(uout,'(a,t35,3g15.3)') 'aspect ratio (>1 pro, <1 ob)   = ', aellipsoid
            write(uout,'()')
         end if

         if(lsuperball) then
            write(uout,'()')
            write(uout,'(a)') 'particles are superballs'
            write(uout,'(a)') '------------------------'
            write(uout,'(a,t35,3g15.3)') 'radsuperball                   = ', radsuperball
            write(uout,'(a,t35,3g15.3)') 'qsuperball                     = ', qsuperball
            write(uout,'(a,t35,3g15.3)') 'txmethodsuperball              = ', txmethodsuperball
            if (txmethodsuperball == 'nr') then
            write(uout,'(a,t35,3g15.3)') 'nitersuperball                 = ', nitersuperball
            write(uout,'(a,t35,3g15.3)') 'tolsuperball                   = ', tolsuperball
            else if (txmethodsuperball == 'mesh') then
            write(uout,'(a,t35,3g15.3)') 'meshdepthsuperball             = ', meshdepthsuperball
            else if (txmethodsuperball == 'test') then
            write(uout,'(a,t35,3g15.3)') 'meshdepthsuperball             = ', meshdepthsuperball
            write(uout,'(a,t35,3g15.3)') 'nitersuperball                 = ', nitersuperball
            write(uout,'(a,t35,3g15.3)') 'tolsuperball                   = ', tolsuperball
            end if
            write(uout,'(a,t35,3g15.3)') 'lstatsuperball                 = ', lstatsuperball
            write(uout,'(a,t35,2f10.5)') 'rcutsuperball                  = ', sqrt(rcut2superball(1:2))
            write(uout,'()')
         end if

         write(uout,'()')
         write(uout,'(a)') 'particle pair              potential'
         write(uout,'(a)') '-------------              ---------'
         write(uout,'(a,t30,a)') (txptpt(iptjpt),txpot(iptjpt),iptjpt = 1,nptpt)

         write(uout,'()')
         write(uout,'(a)') 'atom pair         energy coefficients (ipot and ucoff)'
         write(uout,'(a)') '---------         ------------------------------------'
         do iatjat = 1, natat
            write(uout,'(a,t35,t22,6(i3,e14.6))') txatat(iatjat), (ipot(m,iatjat),ucoff(m,iatjat),m = 1,npot(iatjat))
         end do

         if (count(lucoffmod) > 0) then
            write(uout,'()')
            write(uout,'(a,i3)') 'number of 1/r coefficients than has been changed/added based on zat =',count(lucoffmod)
         end if
         deallocate(lucoffmod)

      end if

  end select

contains

!........................................................................

subroutine WritePotLekkerkerkerTuinier
   use PotentialModule
   implicit none

   integer(2), parameter :: nltpot = 2
   character(2), parameter :: txltpot(nltpot) = [ 'ss', 'sw' ]
   integer(2), parameter :: nhdiv = 100
   real(8) :: hlow, hupp, dh, h, u0, u1, u2
   integer(4) :: iat, iltpot, ih

   iat = 1
   hlow = zero
   hupp = three*rad_dep
   dh = (hupp-hlow)/(nhdiv-1)

   call FileOpen(uuser, fuser, 'form/noread')
   write(uuser,'(a)') 'LekkerkerkerTuinier depletion potential'
   write(uuser,'(i5)') nltpot
   do iltpot = 1, nltpot
      write(uuser,'(a)') txltpot(iltpot)
      write(uuser,'(i5)') nhdiv
      do ih = 1, nhdiv
         h = hlow + (ih-1)*dh
         call calc_LT_depletion(txltpot(iltpot), radat(iat), rad_dep, rho_dep, factor_dep, h, u0, u1, u2)
         write(uuser,'(4g15.5)') h/rad_dep, u0/beta
      end do
   end do
   close(uuser)
end subroutine WritePotLekkerkerkerTuinier

!........................................................................

end subroutine IOPotTwoBody

!************************************************************************
!> \page potential potential.F90
!! **PotTwoBodyTab1**
!! *set tables for two-body interactions*
!************************************************************************

! ... set tables for two-body interactions

!       input :     utoltab, ftoltab
!                   umaxtab, fmaxtab
!                   txpot
!                   rcut
!                   ucoffx
!                   r2atat(iatjat)
!       output:     ubuf
!                   nugrid(iatjat)
!                   iubuflow(iatjat)

subroutine PotTwoBodyTab1(lwrite)

   use PotentialModule
   implicit none

   character(40), parameter :: txroutine ='PotTwoBodyTab1'
   logical, intent(in) :: lwrite

   integer(4), save :: nbufinit = 1000, nbuffac = 3, nbufmax = 100000
   external PotDefault
   external Pot_1_6_12
   external MCY
   external Ramp
   external SquareWell
   external SphEwaldTrunc
   external Nemo
   external AsakuraOosawa
   external LekkerkerkerTuinier
   integer(4) :: ipt, jpt, iptjpt, iatjat, ierr, ibuf
   logical    :: lsetconf

   if (ltrace) call WriteTrace(2, txroutine, iWriteInput)

   if(.not.allocated(nugrid)) then
      allocate(iubuflow(natat), nugrid(natat), rumin(natat), rumax(natat), r2umin(natat), r2umax(natat))
      iubuflow = 0
      nugrid = 0
      rumin = 0.0E+00
      rumax = 0.0E+00
      r2umin = 0.0E+00
      r2umax = 0.0E+00
   end if
   if(.not.allocated(lsetatat)) then
      allocate(lsetatat(natat))
      lsetatat = .false.
   end if

   nbuf = nbufinit                              ! initial memory size for potential table

   do while (nbuf <= nbufmax)
      if (allocated(ubuf)) then
         deallocate(ubuf)
         nbuf = nbuffac*nbuf                    ! increase memory for potential table
      end if
      allocate(ubuf(nbuf), stat = ierr)
      ubuf = 0.0E+00
      if(ierr /= 0) call WriteIOStat(txroutine, 'memory allocation of ubuf failed', ierr, 2, 6)

! ... loop over all types of particle pairs

      lmorememory = .false.

      ibuf = 1                                  ! initiate ibuf
      lsetatat =.false.                         ! initiate lsetatat
      do ipt = 1, npt
         do jpt = ipt, npt
            iptjpt = iptpt(ipt,jpt)

! ... check for special potentials

            lsetconf =.true.
            if (txpot(iptjpt) == '(1,6,12)') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, Pot_1_6_12)
            else if (txpot(iptjpt) == 'mcy') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, MCY)
            else if (txpot(iptjpt) == 'setx' .or. txpot(iptjpt) == 'setc' .or. txpot(iptjpt) == 'setd') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, SphEwaldTrunc)
            else if (txpot(iptjpt)(1:5) == 'nemo:') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, Nemo)
            else if (txpot(iptjpt) == 'ramp') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, Ramp)
            else if (txpot(iptjpt) == 'sw') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, SquareWell)
            else if (txpot(iptjpt) == 'asakura-oosawa') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, AsakuraOosawa)
            else if (txpot(iptjpt) == 'lekkerkerker-tuinier') then
               call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, LekkerkerkerTuinier)
            else
               call PotentialUser(ipt, jpt, ibuf, lsetatat, lsetconf)
            end if

! ... default potential

            if (.not.lsetconf) call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat, PotDefault)

            if (lmorememory) goto 10
         end do
      end do
      if(.not.lmorememory) exit
10    continue
   end do
   if(nbuf > nbufmax) call stop(txroutine, 'memory allocation to potenetial table is insufficient', uout)

! ... write data about the two-body potential tables

   if (master .and. lwrite) then
      write(uout,'()')
      write(uout,'(a)') 'specification of two-body potential tables'
      write(uout,'(a)') '------------------------------------------'
      write(uout,'(a,t20,a,t40,a,t55,a,t65,a)') 'atom pair', 'no of grids', 'iubuflow', 'rumin', 'rumax'
      write(uout,'(a,t20,a,t40,a,t55,a,t65,a)') '---------', '-----------', '--------', '-----', '-----'
      write(uout,'(a,t15,i10,t35,i10,t50,2f10.3)')       &
      (txatat(iatjat),nugrid(iatjat),iubuflow(iatjat),rumin(iatjat),rumax(iatjat),iatjat = 1,natat)
      write(uout,'()')
      write(uout,'(a,t35,3g10.3)') 'nbuf (needed, allocated)        = ', iubuflow(natat)+12*nugrid(natat), nbuf
   end if

end subroutine PotTwoBodyTab1

!************************************************************************
!> \page potential potential.F90
!! **PotTwoBodyTab2**
!! *set tables for two-body interactions for particle pair iptjpt*
!************************************************************************

subroutine PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat2, potsub)

   use PotentialModule
   implicit none
   integer(4), intent(in)    :: ipt, jpt
   integer(4), intent(inout) :: ibuf
   logical,    intent(inout) :: lsetatat2(*)

   character(40), parameter :: txroutine ='PotTwoBodyTab2'
   integer(4) :: iptjpt, iat, jat, iatjat, idum1, idum2
   real(8)    :: dum
   external potsub

   if (ltrace) call WriteTrace(2, txroutine, iWriteInput)

   call potsub('init', ipt, jpt, idum1, idum2, dum, dum, dum, dum) ! get initial values of r2umin and r2umax for interaction iptjpt

   iptjpt = iptpt(ipt,jpt)
   do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
      do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
         iatjat = iatat(iat,jat)
         if (lsetatat2(iatjat)) cycle                  ! interaction iatjat is already set
         rumax(iatjat) = sqrt(r2umax(iatjat))          ! maximum value of r for interaction iatjat
         rumin(iatjat) = sqrt(r2umin(iatjat))          ! lower minimum of r for interaction iatjat
         call PotTwoBodyTab3(iat, jat, ibuf, potsub)   ! set iubuflow, nugrid, and ubuf for interaction iatjat
         if (lmorememory) return
         r2umin(iatjat) = rumin(iatjat)**2             ! save r2umin for ineraction iatjat
         lsetatat2(iatjat) =.true.                     ! interaction iatjat is now set
      end do
   end do

   if (itestpot == 1 .and. master) call TestPotTwoBodyTab2(ipt, jpt, potsub, uout)

contains

!........................................................................

subroutine TestPotTwoBodyTab2(ipt, jpt, potsub, unit)
   use PotentialModule
   implicit none
   integer(4), intent(in)    :: ipt, jpt
   integer(4), intent(in)    :: unit
   integer(4) ::  iat, jat, iatjat
   external  potsub
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   if (ilist == 1) then
      write(ulist,'(a,i5,a)') 'twobody potential: ',iptpt(ipt,jpt),txptpt(iptpt(ipt,jpt))
      write(ulist,'(i5)') natpt(ipt)*natpt(jpt)
   end if
   do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
      do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
         iatjat = iatat(iat,jat)
         if (jat < iat) cycle
         call WriteUBuffer(iatjat)
         call TestUBuffer(iat, jat, rumin(iatjat), rumax(iatjat), potsub, unit)
      end do
   end do
end subroutine TestPotTwoBodyTab2

!........................................................................

end subroutine PotTwoBodyTab2

!************************************************************************
!> \page potential potential.F90
!! **PotTwoBodyTab3**
!! *set table for two-body interaction iatjat*
!************************************************************************

!     interpolation procedure is taken from andrea, swope, and andersen JCP 79, 4576 (1983)

subroutine PotTwoBodyTab3(iat, jat, ibuf, potsub)

   use PotentialModule
   implicit none

   integer(4), intent(in)    :: iat, jat
   integer(4), intent(inout) :: ibuf

   character(40), parameter :: txroutine ='PotTwoBodyTab3'
   integer(4),    parameter :: mnugrid = 1200 ! maximum number of potential grids
   integer(4),    parameter :: ndr = 100      ! maximum nuber of iterations in dr
   real(8),       parameter :: drfac = 0.9d0  ! factor by which dr is reduced in succesive iteractions

   logical    :: ok, repul
   integer(4) :: iatjat, igrid, idr, idum
   real(8)    :: rlow, rupp, zlow, zupp, dr
   real(8)    :: u0low, u1low, u2low, u0upp, u1upp, u2upp, ubuft(12)
   external potsub

   if (ltrace) call WriteTrace(3, txroutine, iWriteInput)

   iatjat = iatat(iat,jat)
   iubuflow(iatjat) = ibuf                    ! save location of start of table for iatjpt interaction
   rupp = rumax(iatjat)                       ! initiate rupp
   dr = Half*(rumax(iatjat)-rumin(iatjat))    ! initiate dr
   do igrid = 1, mnugrid                      ! loop over intervals

      zupp = rupp**2                          ! set zupp
      do idr = 1, ndr                         ! loop over decreasing dr
        rlow = rupp-dr                        ! calculate lower end of interval
        rlow = max(rumin(iatjat),rlow)        ! check against rumin
        zlow = rlow**2                        ! set zlow

! ... calculate u0, u1, and u2 for rlow and rupp

        call potsub('calc', idum, idum, iat, jat, rlow, u0low, u1low, u2low)
        call potsub('calc', idum, idum, iat, jat, rupp, u0upp, u1upp, u2upp)

! ... transform to r**2, calculate coefficients, and store coefficients

         call SetUBuffer(zlow, u0low, u1low, u2low, zupp, u0upp, u1upp, u2upp, ubuft)

! ... check accurary of the used dr

         call CheckUBuffer(iat, jat, rlow, rupp, ubuft, potsub, ok, repul)
         if (ok) then
            rupp = rlow                                        ! set rupp for next interval
            if (txpot(1) == 'ramp') dr = Two*dr/drfac          ! inserted 2010-10-31
            if (txpot(1) == 'sw') dr = dr/drfac                ! inserted 2002-12-13
            if (txpot(1) == 'ewald_square_well') dr = dr/drfac ! inserted 2002-12-13
            if (lljcut) dr = Two*dr                            ! for LJ potential: overcome the force break
            exit
         end if

         dr = dr*drfac                                         ! accuracy is not sufficient, reduce dr, and try again
      end do

      if (idr >= ndr+1) write(*,*) 'iatjat, igrid', iatjat, igrid
      if (idr >= ndr+1) call Stop(txroutine, 'idr >= ndr+1', uout)   ! no convergence in dr?
      if (ibuf+11 > nbuf) then                                       ! no space enough to store current set of table data?
         lmorememory = .true.
         return
      end if
      ubuf(ibuf:ibuf+11) = ubuft(1:12)                               ! save table data
      ibuf = ibuf+12                                                 ! update ibuf to next empty slot of ubuf
      if (repul) rumin(iatjat) = rlow                                ! if the strongly repulsive region is reached, update rumin

      if (rlow <= rumin(iatjat)) exit                                ! exit if the total interval is calculated
      if (repul) exit                                                ! exit if the energy or force reaches a max
   end do

   if (igrid >= mnugrid+1) then
      write(uout,'(a,2i5,2e12.5)') 'igrid, idr, rlow, dr', igrid, idr, rlow, dr
      call Stop(txroutine, 'igrid >= mnugrid+1', uout) ! check if converged in fining rlow
   end if
   nugrid(iatjat) = igrid                                            ! store nugrid of interaction iatjat

end subroutine PotTwoBodyTab3

!************************************************************************
!> \page potential potential.F90
!! **SetUBuffer**
!! *transform to \f$ r^2 \f$, calculate coefficients, and store them in a table*
!************************************************************************

subroutine SetUBuffer(zlow, u0low, u1low, u2low, zupp, u0upp, u1upp, u2upp, ubuft)

   use PotentialModule
   implicit none

   real(8), intent(in)  :: zlow, u0low, u1low, u2low
   real(8), intent(in)  :: zupp, u0upp, u1upp, u2upp
   real(8), intent(out) :: ubuft(12)

   real(8) :: dz1, dz2, dz3
   real(8) :: rlow, rupp
   real(8) :: w0low, w1low, w2low
   real(8) :: w0upp, w1upp, w2upp
   real(8) :: a, b, c

!   write(*,'(a,4g15.5)') 'SetUBuffer, zlow, u0, u1, u2', zlow, u0low, u1low, u2low
!   write(*,'(a,4g15.5)') 'SetUBuffer, zupp, u0, u1, u2', zupp, u0upp, u1upp, u2upp

! ... a non-interacting pair has been reached

   if ((u0low == Zero) .and. (u0upp == Zero)) then
      ubuft(1) = zlow
      ubuft(2:12) = Zero
      return
   end if

   rlow = sqrt(zlow)
   rupp = sqrt(zupp)
   dz1 = zupp-zlow
   dz2 = dz1**2
   dz3 = dz2*dz1

!   write(*,'(a,4g15.5)') 'SetUBuffer, rlow, rupp', rlow, rupp
!   write(*,'(a,4g15.5)') 'SetUBuffer, dz1, dz2, dz3', dz1, dz2, dz3

! ... change independent variable from r to z = r**n(power), n = 2

!     w(z) = u(z**(1/n)) = u(r)
!     w'(z) = 1/n*z**(1/n-1)*u'(z**(1/n)) = r/(n*z)*u'(r)
!     w''(z) = 1/(n*z)**2*r*r*u''(r)+(1/n)*(1/n-1)*z**(1/n-2)*u'(r)
!           = (r/(n*z))**2*(u''(r)+(1-n)/r*u'(r))

   w0low = u0low
   w0upp = u0upp
   w1low = 0.50d0*u1low/rlow
   w1upp = 0.50d0*u1upp/rupp
   w2low = 0.25d0*(u2low-u1low/rlow)/zlow
   w2upp = 0.25d0*(u2upp-u1upp/rupp)/zupp

!   write(*,'(a,4g15.5)') 'SetUBuffer, w0low, w1, w2', w0low, w1low, w2low
!   write(*,'(a,4g15.5)') 'SetUBuffer, w0upp, w1, w2', w0upp, w1upp, w2upp

   a = 6.0d0*(w0upp-w0low-w1low*dz1-  0.5d0*w2low*dz2)/dz3
   b = 2.0d0*(w1upp-      w1low            -w2low*dz1)/dz2
   c =      (w2upp-                        w2low     )/dz1

!   write(*,'(a,4g15.5)') 'SetUBuffer, a, b, c', a, b, c

! ... store table in the temporary buffert ubuft

   ubuft(1) = zlow
   ubuft(2) = w0low
   ubuft(3) = w1low
   ubuft(4) = 0.5d0*w2low
   ubuft(5) = ( 10.0d0*a-12.0d0*b+3.0d0*c)/(6.0d0)
   ubuft(6) = (-15.0d0*a+21.0d0*b-6.0d0*c)/(6.0d0*dz1)
   ubuft(7) = (  2.0d0*a- 3.0d0*b+      c)/(2.0d0*dz2)
   ubuft(8) =-2.0d0*ubuft(3)
   ubuft(9) =-4.0d0*ubuft(4)
   ubuft(10) =-6.0d0*ubuft(5)
   ubuft(11) =-8.0d0*ubuft(6)
   ubuft(12) =-10.0d0*ubuft(7)

!   write(*,'(a,12g15.5)') 'buf, ubuft', ubuft(1)
!   write(*,'(a,12g15.5)') 'buf, ubuft', ubuft(2:7)
!   write(*,'(a,12g15.5)') 'buf, ubuft', ubuft(8:12)

end subroutine SetUBuffer

!************************************************************************
!> \page potential potential.F90
!! **CheckUBuffer**
!! *check the accuracy of the tabulation of one interval*
!************************************************************************

subroutine CheckUBuffer(iat, jat, rlow, rupp, ubuft, potsub, ok, repul)

   use PotentialModule
   implicit none

   integer(4), intent(in)  :: iat, jat
   real(8),    intent(in)  :: rlow, rupp
   real(8),    intent(in)  :: ubuft(12)
   logical,    intent(out) :: ok, repul

   integer(4), parameter :: ncheck = 11
   integer(4):: icheck, idum, iatjat
   real(8)   :: dr, dz, r1
   real(8)   :: u0, u1, u2
   real(8)   :: usum, fsum

   iatjat = iatat(iat,jat)
   ok    = .false.
   repul = .false.
   dr = (rupp-rlow)/(ncheck-1)
   do icheck = 1, ncheck
      r1 = rlow+dr*(icheck-1)
      call potsub('calc', idum, idum, iat, jat, r1, u0, u1, u2)
      dz = r1**2-rlow**2
      usum = ubuft(2)+dz*(ubuft(3)+dz*(ubuft(4)+dz*(ubuft(5)+dz*(ubuft(6)+dz*ubuft(7)))))
      fsum = ubuft(8)+dz*(ubuft(9)+dz*(ubuft(10)+dz*(ubuft(11)+dz*ubuft(12))))
      fsum = r1*fsum
      if (abs(usum-u0) > utoltab) return
      if (abs(fsum+u1) > ftoltab) return
      if (abs(usum) > umaxtab) repul = .true.
      if (abs(fsum) > fmaxtab) repul = .true.
   end do
   ok =.true.

end subroutine CheckUBuffer

!************************************************************************
!> \page potential potential.F90
!! **PotDefault**
!! *default potential*
!************************************************************************

subroutine PotDefault(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none
   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='PotDefault'
   integer(4) :: iatjat, m
   real(8) :: c0fac, c1fac, c2fac
   real(8) :: r1i, rni, w0, w1, w2

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2uminin
            if ((zatalpha(iat) > zero) .or. (zatalpha(jat) > zero)) r2umin(iatjat) = 1e-10               ! GCD
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      r1i = One/r1
      rni = One
      w0 = Zero
      w1 = Zero
      w2 = Zero
      do m = 1, 1
         if (lscrc) then
            call CalcScrCUfac(r1, scrlen, c0fac, c1fac, c2fac)
         else if (lcharge .and. lewald) then
            call CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)
            if ((zatalpha(iat) > zero) .or. (zatalpha(jat) > zero)) &                                   ! GCD
               call CalcGCDUfac(r1, zatalpha(iat), zatalpha(jat), c0fac, c1fac, c2fac)
         else if (lweakcharge .and. lewald) then
            call CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)
            if ((zatalpha(iat) > zero) .or. (zatalpha(jat) > zero)) &                                   ! GCD
               call CalcGCDUfac(r1, zatalpha(iat), zatalpha(jat), c0fac, c1fac, c2fac)
         else if (lcharge .and. lrf) then
            call CalcRFUfac(r1, rffac, c0fac, c1fac, c2fac)
         else
            c0fac = One
            c1fac = One
            c2fac = One
            if ((zatalpha(iat) > zero) .or. (zatalpha(jat) > zero)) &                                   ! GCD
               call CalcGCDUfac(r1, zatalpha(iat), zatalpha(jat), c0fac, c1fac, c2fac)
         end if
         rni  = rni * r1i
         w0   = w0 + c0fac*(ucoffx(m,iatjat) * rni)
         w1   = w1 + c1fac*(-m) * (ucoffx(m,iatjat) * rni)
         w2   = w2 + c2fac*(-m) * (-m-1) * (ucoffx(m,iatjat) * rni)
      end do
      do m = 2, ipotm
         rni  = rni * r1i
         w0   = w0 + (ucoffx(m,iatjat) * rni)
         w1   = w1 + (-m) * (ucoffx(m,iatjat) * rni)
         w2   = w2 + (-m) * (-m-1) * (ucoffx(m,iatjat) * rni)
      end do
      u0 = w0
      u1 = w1 * r1i
      u2 = w2 * r1i * r1i

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine PotDefault

!************************************************************************
!> \page potential potential.F90
!! **Pot_1_6_12**
!! *1, 6, 12-potential*
!************************************************************************

subroutine Pot_1_6_12(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='Pot_1_6_12'
   integer(4) :: iatjat
   real(8) :: r1i, qq, sig, eps, fac, c0fac, c1fac, c2fac
   real(8) :: ushort0, ushort1, ushort2

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
        do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
           iatjat = iatat(iat,jat)
!           r2umin(iatjat) = r2atat(iatjat)-1.0d-4
            r2umin(iatjat) = r2uminin
           r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
        end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      qq = zat(iat)*zat(jat)*EpsiFourPi
      if (ldipole .or. ldipolesph) qq = Zero
      if (lpolarization) qq = Zero
      sig = 0.5*(sigat(iat)+sigat(jat))
      eps = sqrt(epsat(iat)*epsat(jat))

      r1i = One/r1
      fac = (sig*r1i)**6

      if (lscrc) then
         call CalcScrCUfac(r1, scrlen, c0fac, c1fac, c2fac)
      else if (lcharge .and. lewald) then
         call CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)
      else if (lweakcharge .and. lewald) then
         call CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)
      else if (lcharge .and. lrf) then
         call CalcRFUfac(r1, rffac, c0fac, c1fac, c2fac)
      else
         c0fac = One
         c1fac = One
         c2fac = One
      end if

      ushort0 = +Four*eps*(fac**2-fac)
      ushort1 = -24.0d0*eps*(Two*fac**2-fac)*r1i
      ushort2 = +24.0d0*eps*(26.0d0*fac**2-7.0d0*fac)*r1i**2
      if (lljcut) then
         if (r1 < ljrcut*sig) then
            ushort0 = ushort0 + ljushift*eps
         else
            ushort0 = Zero
            ushort1 = Zero
            ushort2 = Zero
         end if
      end if

      u0 = qq*c0fac*r1i        + ushort0
      u1 =-qq*c1fac*r1i**2     + ushort1
      u2 = Two*qq*c2fac*r1i**3 + ushort2

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine Pot_1_6_12

!************************************************************************
!> \page potential potential.F90
!! **MCY**
!! *MCY water potential*
!************************************************************************

! ... MCY water potential

!     ref: JCP 64, 1351 (1976). the Four sites are assumed to be ordered o, h1, h2, and m.

subroutine MCY(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='MCY'
   integer(4) :: iatjat
   real(8), save    :: a1, a2, a3, a4, qq, b1, b2, b3, b4
   real(8)    :: r1i, term1, term2, term3, term4, term5

   integer(4), save :: ioo, ioh, ihh, ihm, imm

   if (str(1:4) == 'init') then

      a1 = CalToJ * 1088213.2
      a2 = CalToJ *     666.3373
      a3 = CalToJ *    1455.427
      a4 = CalToJ *     273.5954
      qq = CalToJ *     170.9389
      b1 = - 5.152712
      b2 = - 2.760844
      b3 = - 2.961895
      b4 = - 2.233264

      iat = iatpt(ipt)
      ioo = iatat(iat  ,iat  )
      ioh = iatat(iat  ,iat+1)
      ihh = iatat(iat+1,iat+1)
      ihm = iatat(iat+1,iat+2)
      imm = iatat(iat+2,iat+2)

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      r1i   = One/r1
      if (iatjat == ioo) then
         term1 = a1*exp(max(-expmax,b1*r1))
         u0 = term1
         u1 = b1*term1
         u2 = b1**2*term1
      else if (iatjat == ioh) then
         term3 = a3*exp(max(-expmax,b3*r1))
         term4 = a4*exp(max(-expmax,b4*r1))
         u0 = term3 - term4
         u1 = b3*term3 - b4*term4
         u2 = b3**2*term3   - b4**2*term4
      else if (iatjat == ihh) then
         term2 = a2*exp(max(-expmax,b2*r1))
         term5 = qq*r1i
         u0 = term5 + term2
         u1 =-term5*r1i + b2*term2
         u2 = 2*term5*r1i**2 + b2**2*term2
      else if (iatjat == ihm) then
         term5 = qq*r1i
         u0 =-2*term5
         u1 = 2*term5*r1i
         u2 =-4*term5*r1i**2
      else if (iatjat == imm) then
         term5 = qq*r1i
         u0 = 4*term5
         u1 =-4*term5*r1i
         u2 = 8*term5*r1i**2
      else
         u0 = Zero
         u1 = Zero
         u2 = Zero
      end if

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine MCY

!************************************************************************
!> \page potential potential.F90
!! **SphEwaldTrunc**
!! *spherical ewald truncated potential*
!************************************************************************

!          txpot  =   'setx'  1/r -> (1-erf(r*ualpha))/r
!                     'setc'  1/r ->  1/r
!                     'setd'  1/r ->  erf(r*ualpha)/r

!          where ualpha = ualphared/rcut

subroutine SphEwaldTrunc(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8),      intent(in)    :: r1
   real(8),      intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='SphEwaldTrunc'
   real(8), parameter :: fac = 1.128379167d0
   integer(4) :: iptjpt, iatjat, m
   real(8)    :: er, ex, r1i, rni, w0, w1, w2, term, rb
   real(8), external :: ErfLocal
   real(8), save    :: coul, sett

   if (str(1:4) == 'init') then

      ualpha = ualphared/rcut
      iptjpt = iptpt(ipt,jpt)

      if (txpot(iptjpt) == 'setx') then
         coul = One
         sett = One
      else if (txpot(iptjpt) == 'setc') then
         coul = One
         sett = Zero
      else if (txpot(iptjpt) == 'setd') then
         coul = Zero
         sett = One
      else
         call Stop(txroutine, 'unsupported value of txpot', uout)
      end if

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
!           r2umin(iatjat) = r2atat(iatjat)-1.0d-4
            r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      r1i = One/r1
      rni = r1i
      w0 = Zero
      w1 = Zero
      w2 = Zero
      do m = 2, ipotm
         rni  = rni * r1i
         term = ucoffx(m,iatjat) * rni
         w0   = w0 + term
         w1   = w1 + (-m) * term
         w2   = w2 + (-m) * (-m-1) * term
      end do
      u0 = w0
      u1 = w1 * r1i
      u2 = w2 * r1i * r1i

      rb = ualpha*r1
      er = ErfLocal(rb)
      ex = fac*exp(-rb**2)
      u0 = u0 + ucoffx(1,iatjat)*(coul-sett*(er))*r1i
      u1 = u1 - ucoffx(1,iatjat)*(coul-sett*(er-rb*ex))*r1i**2
      u2 = u2 + ucoffx(1,iatjat)*(coul-sett*(er-rb*(One+rb**2)*ex))*Two*r1i**3

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine SphEwaldTrunc

!************************************************************************
!> \page potential potential.F90
!! **Ramp**
!! *ramp potential (with soft slope change)*
!************************************************************************

!     potential: u(r) =  slope*(r-r_ramp)*(0.5 + (1/pi)*atan(alpha*(r-r_ramp))
!                slope = epsilon_ramp/r1atat*(lambda_ramp-One)
!                r_ramp = r1atat*(lambda_ramp)
!          for
!                r1atat < r < r1atat*(lambda_ramp-One)

subroutine Ramp(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='Ramp'
   real(8)    :: r_ramp, slope, x, f0_ramp, f1_ramp, f2_ramp
   integer(4) :: iatjat

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r_ramp = lambda_ramp*r1atat(iatjat)
            if (r_ramp > rcut) call stop(txroutine,' r_ramp > rcut', uout)  ! check that ramp inside rcut
            r2umin(iatjat) = r2atat(iatjat)-1.0d-4
!           r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      r_ramp = lambda_ramp*r1atat(iatjat)
      slope = epsilon_ramp/(r1atat(iatjat)*(lambda_ramp-One))
      x = r1 - r_ramp
      f0_ramp = (Half - (One/Pi)*atan(x/alpha_ramp))
      f1_ramp = -alpha_ramp / (alpha_ramp**2 + x**2) / Pi
      f2_ramp = -Two*x* f1_ramp / (alpha_ramp**2 + x**2)
      u0 = slope*(x*f0_ramp)
      u1 = slope*(f0_ramp + x*f1_ramp)
      u2 = slope*(2*f1_ramp + x*f2_ramp)

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine Ramp

!************************************************************************
!> \page potential potential.F90
!! **SquareWell**
!! *Square-well potential (with a soft rise)*
!************************************************************************

!     potential: u(r) =  epsilon_sw*(0.5 + (1/pi)*atan(alapha_sw*(r-r_sw))
!                r_sw = r1atat*(lambda_sw)

subroutine SquareWell(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='SquareWell'
   real(8)    :: r_sw, x
   integer(4) :: iatjat

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2atat(iatjat)-1.0d-4
!           r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2
        end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      r_sw = lambda_sw*r1atat(iatjat)
      x = r1 - r_sw
      u0 = -epsilon_sw*(Half - (One/Pi)*atan(x/alpha_sw))
      u1 = epsilon_sw/Pi * alpha_sw / (alpha_sw**2 + x**2)
      u2 = -Two*x* u1 / (alpha_sw**2 + x**2)

      if(txuser == 'jurij_lys') then
         if(iatjat /= 3) then               ! only sw for iatjat == 3
            u0 = zero
            u1 = zero
            u2 = zero
         end if
      end if

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine SquareWell

!************************************************************************
!> \page potential potential.F90
!! **AsakuraOosawa**
!! *Asakura-Oosawa depletion potential*
!************************************************************************


! depletion potential due to nonadsorbing ideal polymers represented as
! penetrable hard spheres with radius of gyration rad_dep and at number density rho_dep
!
! Asakura and Oosawa, J. Polymer Sci. 1958, 33, 183

subroutine AsakuraOosawa(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='AsakuraOosawa'
   integer(4) :: iatjat
   real(8)    :: rd, fac

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2atat(iatjat)
            r2umax(iatjat) = (rcut+racom(ipt)+racom(jpt))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      rd=radat(1)+rad_dep
      fac=r1/rd
      if (r1 <= Two*rd) then
         u0 = -rho_dep*FourPiThird*rd**3*(one-0.75d0*fac+fac**3/16.d0)/beta
         u1 = -rho_dep*Pi*rd**2*(fourth*fac**2-one)/beta
         u2 = -rho_dep*Pi*r1*half/beta
      else
         u0 = Zero
         u1 = Zero
         u2 = Zero
      endif

   else

      call Stop('AsakuraOosawa', 'str value unsupported', uout)

   end if

end subroutine AsakuraOosawa

!************************************************************************
!> \page potential potential.F90
!! **LekkerkerkerTuinier**
!! *Lekkerkerker-Tuinier depletion potential*
!************************************************************************


! depletion potential due to the nonadsorbing ideal polymer modelled as a
! "penetrable hard sphere" with radius of gyration rgyr and density rhophs

! Tuinier, Aarts, Wensink, Lekkerkerker, Phys. Chem. Chem. Phys. 2003, 5, 3707, eq 2, 3, and 16

! and

! Lecture Notes in Physics Volume 833 (2011)
! Colloids and the Depletion Interaction
! Authors: Lekkerkerker, Tuinier
! ISBN: 978-94-007-1222-5 (Print) 978-94-007-1223-2 (Online)
! Chapter 2, p. 147, eq 4.25

subroutine LekkerkerkerTuinier(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='LekkerkerkerTuinier'
   integer(4) :: iatjat
   real(8)    :: h

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2atat(iatjat)
            r2umax(iatjat) = ((rcut+racom(ipt)+racom(jpt)))**2
         end do
      end do

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      h = r1 - two*radat(1)      ! surface-surface separation
      call calc_LT_depletion('ss', radat(iat), rad_dep, rho_dep, factor_dep, h, u0, u1, u2)
      u0 = u0/beta
      u1 = u1/beta
      u2 = u2/beta

   else

      call Stop('LekkerkerkerTuinier', 'str value unsupported', uout)

   end if

end subroutine LekkerkerkerTuinier

!************************************************************************
!> \page potential potential.F90
!! **calc_LT_depletion**
!! *calculate Lekkerkerker-Tuinier depletion potential (in kT units)*
!************************************************************************


!                   phi1
! pot = - int (1/phi * dP * G) dphi
!                    0

subroutine calc_LT_depletion(geometry, R, Rg, phi1, factor_dep, h, pot0, pot1, pot2)
   implicit none
   character(2), intent(in) :: geometry         ! wall-wall, sphere-wall or sphere-sphere
   real(8), intent(in)      :: R                ! colloid radius
   real(8), intent(in)      :: Rg               ! polymer radius of gyration
   real(8), intent(in)      :: phi1             ! polymer concentration divided by the overlap concentration
   real(8), intent(in)      :: factor_dep       ! depletion-layer factor
   real(8), intent(in)      :: h                ! surface-surface separation
   real(8), intent(out)     :: pot0, pot1, pot2 ! depletion potential and its two first derivatives
   real(8), parameter :: phi0 = 1.0d-10         ! lower number density boundary
   integer(4), parameter :: n = 1000            ! number of points
   real(8) :: y0, dy, phi(n), value(n,3), TrapNe
   integer(4) :: i

   y0 = log(phi0)
   if (phi1 > phi0) then
      dy = (log(phi1)-log(phi0))/(n-1)            ! logarithmic division
      do i = 1, n
         phi(i) = exp((i-1)*dy+y0)
         call integrand_LT_depletion(phi(i), value(i,1), value(i,2), value(i,3))
      enddo

      pot0 = TrapNe(n,phi,value(1,1))
      pot1 = TrapNe(n,phi,value(1,2))
      pot2 = TrapNe(n,phi,value(1,3))
   else
      pot0 = 0.0d0
      pot1 = 0.0d0
      pot2 = 0.0d0
   end if
contains

!........................................................................

subroutine integrand_LT_depletion(phi, value0, value1, value2)
   implicit none
   real(8), parameter ::  Pi = 3.1415926535897932d0
   real(8), parameter :: FourPiThird = 4.0d0*Pi/3.0d0
   real(8), intent(in)  :: phi                      ! polymer concentration divided by the overlap concentration
   real(8), intent(out) :: value0, value1, value2   ! integrand and its two first derivatives
   real(8) :: q, dP, Dw, Ds, G0, G1, G2
   real(8) :: fac0, fac1, fac2, fac1d, fac2d, fac1dd, fac

   dP = 1 + 2.63*phi*((1+3.25*phi+4.15*phi**2)/(1+1.48*phi))**0.309         ! osmotic compressbility
   Dw = factor_dep*(1.07*Rg)/((1+3.95*phi**1.54)**0.5)                      ! depletion thickness near a plane, good solvent
   q  = Rg/R                                                                ! polymer-colloid size ratio
   Ds = factor_dep*R*(0.865* (q**(-2)+3.95*(q**(-2))*phi**1.54)**(-0.44))   ! depletion thickness near a sphere,  good solvent

   G0 = 0                                                                   ! adsorbed amount of polymer segments / phi
   G1 = 0
   G2 = 0
   if (geometry == 'ss') then
      if (h < 2*Ds) then
         fac0 = (2*Pi/3)*(Ds**3)
         fac1 = (1 - h/(2*Ds))**2
         fac2 = 2 + 3*R/Ds + h/(2*Ds)
         fac1d = 2*(1 - h/(2*Ds))/(-2*Ds)
         fac2d = 1/(2*Ds)
         fac1dd = 2/(-2*Ds)**2
         G0 = fac0*fac1*fac2
         G1 = fac0*(fac1d*fac2 + fac1*fac2d)
         G2 = fac0*(fac1dd*fac2 + 2*fac1d*fac2d)
      endif
   else if (geometry == 'sw') then
      if (h < Dw+Ds) then
         fac0 = Pi/3
         fac1 = (Ds + Dw - h)**2
         fac2 = 3*R + 2*Ds - Dw + h
         fac1d =-2*(Ds + Dw - h)
         fac2d = 1
         fac1dd = 2
         G0 = fac0*fac1*fac2
         G1 = fac0*(fac1d*fac2 + fac1*fac2d)
         G2 = fac0*(fac1dd*fac2 + 2*fac1d*fac2d)
      endif
   else if (geometry == 'ww') then
      if (h < 2*Dw) then
         G0 = 2*Dw - h
         G1 =-1
         G2 = 0
      endif
   else
      write(*,*) 'wrong option of geometry'
      stop 1
   endif

   fac = 1.0d0/(FourPiThird*Rg**3) * dP
   value0 = -fac*G0
   value1 = -fac*G1
   value2 = -fac*G2

end subroutine integrand_LT_depletion

!........................................................................

end subroutine calc_LT_depletion

!************************************************************************
!> \page potential potential.F90
!! **Nemo**
!! *nemo potentials*
!************************************************************************

!     modified for new potential tables

!     valid potential types:

!  1)  'exponential repulsion'

!       u(r) = qaqb/(4*pi*eps0*r) + aab*exp(-bab*r) + (cab/r)**20
!            + dab*s(r)/r**6
!       s(r) = 1.0 - exp( -(r/asw)**nsw )

!  2)  'r-7 repusion'

!       u(r) = qaqb/(4*pi*eps0*r) + aab/r**7 + dab*s(r)/r**6
!       s(r) = 1.0 - exp( -(r/asw)**nsw )

!  3)  'modified interactions'

!       u(r) = qaqb/(4*pi*eps0*r) - (1-aab*r**6*exp(-bab*r)*cab/r**6
!            + dab/r**nab

!  4) 'damping exponential'

!       u(r) = qaqb/(4*pi*eps0*r) - (1-aab*r**6*exp(-bab*r)*cab/r**6
!            + (dab/r)**nab + eab*exp(-fab*r)

!  5) 'full damping'

!       u(r) = qaqb/(4*pi*eps0*r) - damping*cab/r**6
!            + (dab/r)**nab + eab*exp(-fab*r)

!  6) 'full damping chtr'

!       u(r) = -kcht*exp(-acht*r) - damping*cab/r**6
!            + (dab/r)**nab + eab*exp(-fab*r)

!     it is assumed that each particle is described by the standard
!     number and order of atoms

subroutine Nemo(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='PotTwoBodyTab2'
   character(20) :: strd
   character(40), save :: ptype
   character(40) :: nemop(7) = ['exponential repulsion            ',    &
                                'r-7 repulsion                    ',    &
                                'modified interactions            ',    &
                                'damping exponential              ',    &
                                'full damping                     ',    &
                                'full damping chtr                ',    &
                                'full damping chtr gaussian       ']
   integer(4) :: iptjpt, iatjat, ioerr, idum
   logical    :: foundpot
   integer(4), save :: nsw
   real(8),    save :: asw
   integer(4), allocatable, save ::  nab(:)
   real(8),    allocatable, save :: qa(:), qb(:), aab(:), bab(:), cab(:), dab(:), eab(:), fab(:), acht(:), kcht(:)

   if (str(1:4) == 'init') then

      if(.not.allocated(nab)) then
         allocate(nab(natat), &
        qa(natat), qb(natat), aab(natat), bab(natat), cab(natat), dab(natat), eab(natat), fab(natat), &
        acht(natat), kcht(natat) )
         nab = 0
         qa = 0.0E+00
         qb = 0.0E+00
         aab = 0.0E+00
         bab = 0.0E+00
         cab = 0.0E+00
         dab = 0.0E+00
         eab = 0.0E+00
         fab = 0.0E+00
         acht = 0.0E+00
         kcht = 0.0E+00
      end if

      iptjpt = iptpt(ipt,jpt)

! ... set asw (in m) and nsw which are idential for all nemo potentials except in type 3

      asw = 1.2e-10/scllen
      nsw = 4

      if (master) then

! ... open potential library file

         call FileOpen(ulib, flib, 'form/noread')

         do

! ... read a line of the library file

            read(ulib,'(a)',iostat = ioerr) strd
            if (ioerr/= 0) then
               close (ulib)
               foundpot = .false.
               exit
            end if

! ... look for potential match, if match read potential type

            if (txpot(iptjpt) == strd) then
               read(ulib,'(a)') ptype
               foundpot =.true.
               exit
            end if

         end do

      end if

#if defined (_PAR_)
     call par_bc_characters(ptype, 1)
     call par_bc_logical(foundpot)
#endif

      if (.not.foundpot) call Stop(txroutine, 'nemo potential not found in library', uout)

! ... loop over atom type-atom type pairs and set up the tables

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            if (jat < iat) cycle
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2

            if (master) then

            if (ptype == nemop(1)) then

! ... nemo type 1
               read(ulib,*,iostat = ioerr) idum, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), cab(iatjat), dab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential', uout)
               qa(iatjat) = qa(iatjat)*EpsiFourPi
               aab(iatjat) = aab(iatjat)*1e+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**20
               dab(iatjat) = dab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**6

            else if (ptype == nemop(2)) then

! ... nemo type 2

               read(ulib,*,iostat = ioerr) idum, qa(iatjat), qb(iatjat), aab(iatjat), dab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential', uout)
               qa(iatjat) = qa(iatjat)*EpsiFourPi
               aab(iatjat) = aab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**7
               dab(iatjat) = dab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**6

            else if (ptype == nemop(3)) then

! ... nemo type 3

               read(ulib,*,iostat = ioerr) idum, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), cab(iatjat), dab(iatjat), nab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential', uout)
               qa(iatjat) = qa(iatjat)*EpsiFourPi
               aab(iatjat) = aab(iatjat)*1e+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**6
               dab(iatjat) = dab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**nab(iatjat)

            else if (ptype == nemop(4)) then

! ... nemo type 4

               read(ulib,*,iostat = ioerr) idum, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                                          cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential', uout)
               qa(iatjat) = qa(iatjat)*EpsiFourPi
               aab(iatjat) = aab(iatjat)*1e+3/sclene
               eab(iatjat) = eab(iatjat)*1e+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               fab(iatjat) = fab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**6
               dab(iatjat) = dab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**nab(iatjat)

            else if (ptype == nemop(5)) then

! ... nemo type 5

               read(ulib,*,iostat = ioerr) idum, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                                          cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential 5', uout)
               qa(iatjat) = qa(iatjat)*EpsiFourPi
               aab(iatjat) = aab(iatjat)*1e+3/sclene
               eab(iatjat) = eab(iatjat)*1e+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               fab(iatjat) = fab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**6
               dab(iatjat) = dab(iatjat)*((1e+3/sclene)*(1d-10/scllen))**nab(iatjat)

! ... nemo type 6
            else if (ptype == nemop(6)) then

               read(ulib,*,iostat = ioerr) idum, kcht(iatjat), acht(iatjat), aab(iatjat), bab(iatjat), &
                                          cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential 6', uout)
               acht(iatjat) = acht(iatjat)*scllen/1d-10
               kcht(iatjat) = kcht(iatjat)*1d+3/sclene
               aab(iatjat) = aab(iatjat)*1d+3/sclene
               eab(iatjat) = eab(iatjat)*1d+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               fab(iatjat) = fab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1d+3/sclene)*(1d-10/scllen))**6
               dab(iatjat) = dab(iatjat)*((1d+3/sclene)*(1d-10/scllen))**nab(iatjat)

! ... nemo type 7
            else if (ptype == nemop(7)) then

               read(ulib,*,iostat = ioerr) idum, kcht(iatjat), acht(iatjat), aab(iatjat), bab(iatjat), &
                                          cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat)
               if (ioerr/= 0) call Stop(txroutine, 'error in reading nemo potential 6', uout)
               acht(iatjat) = acht(iatjat)*scllen/1d-10
               kcht(iatjat) = kcht(iatjat)*1d+3/sclene
               aab(iatjat) = aab(iatjat)*1d+3/sclene
               eab(iatjat) = eab(iatjat)*1d+3/sclene
               bab(iatjat) = bab(iatjat)*scllen/1d-10
               fab(iatjat) = fab(iatjat)*scllen/1d-10
               cab(iatjat) = cab(iatjat)*((1d+3/sclene)*(1d-10/scllen))**6
               dab(iatjat) = dab(iatjat)*((1d+3/sclene)*(1d-10/scllen))**nab(iatjat)

            else

               call Stop(txroutine, 'unsupported form of nemo pot.', uout)

            end if
            end if

         end do
      end do
      close(ulib)

#if defined (_PAR_)
       call par_bc_reals(qa , natat)
       call par_bc_reals(qb , natat)
       call par_bc_reals(aab, natat)
       call par_bc_reals(bab, natat)
       call par_bc_reals(cab, natat)
       call par_bc_reals(dab, natat)
       call par_bc_reals(eab, natat)
       call par_bc_reals(fab, natat)
       call par_bc_ints(nab, natat)
#endif

   else if (str(1:4) == 'calc') then

      iatjat = iatat(iat,jat)
      if (ptype == nemop(1)) then
         call NemoType1(r1, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), asw, nsw, u0, u1, u2)
      else if (ptype == nemop(2)) then
         call NemoType2(r1, qa(iatjat), qb(iatjat), aab(iatjat), dab(iatjat), &
                       asw, nsw, u0, u1, u2)
      else if (ptype == nemop(3)) then
         call NemoType3(r1, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), nab(iatjat), u0, u1, u2)
      else if (ptype == nemop(4)) then
         call NemoType4(r1, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat), u0, u1, u2)
      else if (ptype == nemop(5)) then
         call NemoType5(r1, qa(iatjat), qb(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat), u0, u1, u2)
      else if (ptype == nemop(6)) then
         call NemoType6(r1, acht(iatjat), kcht(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), nab(iatjat), u0, u1, u2)
      else if (ptype == nemop(7)) then
         call NemoType7(r1, acht(iatjat), kcht(iatjat), aab(iatjat), bab(iatjat), &
                       cab(iatjat), dab(iatjat), eab(iatjat), fab(iatjat), u0, u1, u2) !nab(iatjat) is not needed
      else
         call Stop(txroutine, 'illegal value of ptype in calc', uout)
      end if

   end if

end subroutine Nemo

!************************************************************************
!> \page potential potential.F90
!! **NemoType1**
!! *nemo potential of type 1*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = qaqb/r + aab*exp(-bab*r) + (cab/r)**20 + dab*s(r)/r**6
!     s(r) = 1.0 - exp( -(r/asw)**nsw )

subroutine NemoType1(r, qa, qb, aab, bab, cab, dab, asw, nsw, ur, urp, urpp)

   implicit none

   integer(4), intent(in)  :: nsw
   real(8),    intent(in)  :: r, qa, qb, aab, bab, cab, dab, asw
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp, ur20, udis, udisp, udispp
   real(8) :: ss, ssp, sspp

   ss (r,asw,nsw) = One - exp( -min(expmax,(r/asw)**nsw) )
   ssp(r,asw,nsw) = real(nsw)*(r**(nsw-1)/asw**nsw)*            &
                    exp( -min(expmax,(r/asw)**nsw) )
!  sspp(r,asw,ns) = real(nsw)*float(nsw-1)*(r**(nsw-2)/asw**nsw)&
!                  *exp( -min(expmax,(r/asw)**nsw) )             &
!                  - ( real(nsw)*( r**(nsw-1)/asw**nsw ) )**2   &
!                  *exp( -min(expmax,(r/asw)**nsw) )
! corrected call of sspp, 97-03-13 /per l
   sspp(r,asw,nsw) = real(nsw)*real(nsw-1)*(r**(nsw-2)/asw**nsw)&
                   *exp( -min(expmax,(r/asw)**nsw) )             &
                   - ( real(nsw)*( r**(nsw-1)/asw**nsw ) )**2   &
                   *exp( -min(expmax,(r/asw)**nsw) )

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uele   = qa*qb*r1i
   uexp   = aab*exp(-min(expmax,bab*r))
   ur20   = (cab*r1i)**20
   udis   = dab*ss  (r,asw,nsw)*r6i
   udisp  = dab*ssp (r,asw,nsw)*r6i
   udispp = dab*sspp(r,asw,nsw)*r6i

   ur   = uele + uexp + ur20 + udis

   urp  =-uele*r1i       - bab*uexp     - 20.0*ur20*r1i    &
        - 6.0*udis*r1i   + udisp

   urpp = 2.0*uele*r2i  + bab*bab*uexp + 420.0*ur20*r2i    &
        + 42.0*udis*r2i - 12.0*udisp*r1i + udispp

end subroutine NemoType1

!************************************************************************
!> \page potential potential.F90
!! **NemoType2**
!! *nemo potential of type 2*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = qaqb/r + aab/r**7 + dab*s(r)/r**6
!     s(r) = 1.0 - exp( -(r/asw)**nsw )

subroutine NemoType2(r, qa, qb, aab, dab, asw, nsw, ur, urp, urpp)

   implicit none

   real(8),    intent(in)  :: r, qa, qb, aab, dab, asw
   integer(4), intent(in)  :: nsw
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp, udis, udisp, udispp
   real(8) :: ss, ssp, sspp

   ss (r,asw,nsw) = One - exp( -min(expmax,(r/asw)**nsw) )
   ssp(r,asw,nsw) = real(nsw)*(r**(nsw-1)/asw**nsw)*             &
                    exp( -min(expmax,(r/asw)**nsw) )
!  sspp(r,asw,ns) = real(nsw)*real(nsw-1)*(r**(nsw-2)/asw**nsw) &
!                  *exp( -min(expmax,(r/asw)**nsw) )              &
!                  - ( real(nsw)*( r**(nsw-1)/asw**nsw ) )**2    &
!                  *exp( -min(expmax,(r/asw)**nsw) )
! corrected call of sspp, 97-03-13 /per l
   sspp(r,asw,nsw) = real(nsw)*real(nsw-1)*(r**(nsw-2)/asw**nsw) &
                   *exp( -min(expmax,(r/asw)**nsw) )               &
                   - ( real(nsw)*( r**(nsw-1)/asw**nsw ) )**2     &
                   *exp( -min(expmax,(r/asw)**nsw) )

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uele   = qa*qb*r1i
   uexp   = aab*r1i*r6i
   udis   = dab*ss  (r,asw,nsw)*r6i
   udisp  = dab*ssp (r,asw,nsw)*r6i
   udispp = dab*sspp(r,asw,nsw)*r6i

   ur   = uele + uexp + udis

   urp  =-uele*r1i       - 7.0*uexp*r1i - 6.0*udis*r1i    + udisp

   urpp = 2.0*uele*r2i   + 56.0*uexp*r2i                &
        + 42.0*udis*r2i  - 12.0*udisp*r1i + udispp

end subroutine NemoType2

!************************************************************************
!> \page potential potential.F90
!! **NemoType3**
!! *nemo potential of type 3*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = qaqb/(4*pi*eps0*r) + aab*exp(-bab*r) - cab/r**6
!          + dab/r**nab

subroutine NemoType3(r, qa, qb, aab, bab, cab, dab, nab, ur, urp, urpp)

   implicit none

   integer(4), intent(in)  :: nab
   real(8),    intent(in)  :: r, qa, qb, aab, bab, cab, dab
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp, udis1, udis2

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uele   = qa*qb*r1i
   uexp   = dab*r1i**nab
   udis1  =-cab*r6i
   udis2  = aab*exp(-bab*r)

   ur   = uele + uexp + udis1 + udis2

   urp  =-uele*r1i       - nab*uexp*r1i  - 6.0*udis1*r1i   - bab*udis2

   urpp = 2.0*uele*r2i   + nab*(nab+One)*uexp*r2i         &
        + 42.0*udis1*r2i  + bab*bab*udis2

end subroutine NemoType3

!************************************************************************
!> \page potential potential.F90
!! **NemoType4**
!! *nemo potential of type 4*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = qaqb/(4*pi*eps0*r) + aab*exp(-bab*r) - cab/r**6
!          + (dab/r)**nab + eab*exp(-fab*r)

subroutine NemoType4(r, qa, qb, aab, bab, cab, dab, eab, fab, nab, ur, urp, urpp)

   implicit none

   integer(4), intent(in)  :: nab
   real(8),    intent(in)  :: r, qa, qb, aab, bab, cab, dab, eab, fab
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp1, uexp2, udis1, udis2

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uele   = qa*qb*r1i
   uexp1  = eab*exp(-min(expmax,fab*r))
   uexp2  = Zero
   if (nab/= 0) uexp2  = dab*r1i**nab
   udis1  =-cab*r6i
   udis2  = aab*exp(-min(expmax,bab*r))

   ur   = uele + uexp1 + uexp2  + udis1 + udis2

   urp  =-uele*r1i       - fab*uexp1    - nab*uexp2*r1i             &
        - 6.0*udis1*r1i  - bab*udis2

   urpp = 2.0*uele*r2i   + fab*fab*uexp1 + nab*(nab+One)*uexp2*r2i  &
        + 42.0*udis1*r2i  + bab*bab*udis2

end subroutine NemoType4

!************************************************************************
!> \page potential potential.F90
!! **NemoType5**
!! *nemo potential of type 5*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = qaqb/(4*pi*eps0*r) + aab*'damping'(r)*exp(-bab*r) - cab/r**6
!          + (dab/r)**nab + eab*exp(-fab*r)
!

subroutine NemoType5(r, qa, qb, aab, bab, cab, dab, eab, fab, nab, ur, urp, urpp)

   implicit none

   integer(4), intent(in)  :: nab
   real(8),    intent(in)  :: r, qa, qb, aab, bab, cab, dab, eab, fab
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp1, uexp2, udis1, udis2
   real(8) :: bri, ud1, ud2, ud3, ud4, ud5, ud6, udd, uddp, uddpp

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uele   = qa*qb*r1i
   uexp1  = eab*exp(-min(expmax,fab*r))
   uexp2  = Zero
   if (nab/= 0) uexp2  = dab*r1i**nab
   udis1  =-cab*r6i
   udis2  = aab*exp(-min(expmax,bab*r))
!
! new for NemoType5
!
   bri = r1i/bab
   ud1 = 6.0*bri
   ud2 = 5.0*bri*ud1
   ud3 = 4.0*bri*ud2
   ud4 = 3.0*bri*ud3
   ud5 = 2.0*bri*ud4
   ud6 =     bri*ud5
   udd = One + ud1 + ud2 + ud3 + ud4 + ud5 + ud6
   uddp  =-r1i*(    ud1 + 2.0*ud2 +  3.0*ud3 +  4.0*ud4 +  5.0*ud5 +  6.0*ud6)
   uddpp = r2i*(2.0*ud1 + 6.0*ud2 + 12.0*ud3 + 20.0*ud4 + 30.0*ud5 + 42.0*ud6)
!
! modified for NemoType5
!
   ur   = uele + uexp1 + uexp2  + udis1 + udd*udis2

   urp  =-uele*r1i       - fab*uexp1      - nab*uexp2*r1i                          &
        - 6.0*udis1*r1i  - bab*udd*udis2  + uddp*udis2

   urpp = 2.0*uele*r2i    + fab*fab*uexp1 + nab*(nab+One)*uexp2*r2i                &
        + 42.0*udis1*r2i  + bab*bab*udd*udis2 - 2.0*bab*udis2*uddp + uddpp*udis2

end subroutine NemoType5

!************************************************************************
!> \page potential potential.F90
!! **NemoType6**
!! *nemo potential of type 6*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = -kcht*exp(-r*acht) + aab*'damping'(r)*exp(-bab*r) - cab/r**6
!          + (dab/r)**nab + eab*exp(-fab*r)
!

subroutine NemoType6(r, acht, kcht, aab, bab, cab, dab, eab, fab, nab, ur, urp, urpp)

   implicit none

   integer(4), intent(in)  :: nab
   real(8),    intent(in)  :: r, acht, kcht, aab, bab, cab, dab, eab, fab
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp1, uexp2, udis1, udis2, uchtexp
   real(8) :: bri, ud1, ud2, ud3, ud4, ud5, ud6, udd, uddp, uddpp

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uchtexp = -kcht*exp(-min(expmax,acht*r))
   uele   = Zero                          !! 2007-06-05 PL
   uexp1  = eab*exp(-min(expmax,fab*r))
   uexp2  = Zero
   if (nab/= 0) uexp2  = dab*r1i**nab
   udis1  =-cab*r6i
   udis2  = aab*exp(-min(expmax,bab*r))
!
! new for NemoType6
!
   bri = r1i/bab
   ud1 = 6.0*bri
   ud2 = 5.0*bri*ud1
   ud3 = 4.0*bri*ud2
   ud4 = 3.0*bri*ud3
   ud5 = 2.0*bri*ud4
   ud6 =     bri*ud5
   udd = One + ud1 + ud2 + ud3 + ud4 + ud5 + ud6
   uddp  =-r1i*(    ud1 + 2.0*ud2 +  3.0*ud3 +  4.0*ud4 +  5.0*ud5 +  6.0*ud6)
   uddpp = r2i*(2.0*ud1 + 6.0*ud2 + 12.0*ud3 + 20.0*ud4 + 30.0*ud5 + 42.0*ud6)
!
! modified for NemoType6
!
   ur   = uele + uexp1 + uexp2  + udis1 + udd*udis2 + uchtexp

   urp  =-uele*r1i       - fab*uexp1      - nab*uexp2*r1i                          &
        - 6.0*udis1*r1i  - bab*udd*udis2  + uddp*udis2 - acht*uchtexp

   urpp = 2.0*uele*r2i    + fab*fab*uexp1 + nab*(nab+One)*uexp2*r2i                &
        + 42.0*udis1*r2i  + bab*bab*udd*udis2 - 2.0*bab*udis2*uddp + uddpp*udis2   &
        + acht*acht*uchtexp
!  Write(*,'(A,6F18.8)') 'From NemoType6',uchtexp, uexp1,ur, r

end subroutine NemoType6

!************************************************************************
!> \page potential potential.F90
!! **NemoType7**
!! *nemo potential of type 7*
!************************************************************************

!     unit matters are handled by the calling routine

!     u(r) = -kcht*exp(-acht*(r-dab)^2) + aab*'damping'(r)*exp(-bab*r) - cab/r**6
!          + eab*exp(-fab*r)
!

!subroutine NemoType7(r, acht, kcht, aab, bab, cab, dab, eab, fab, nab, ur, urp, urpp) !nap is not used
subroutine NemoType7(r, acht, kcht, aab, bab, cab, dab, eab, fab, ur, urp, urpp)

   implicit none

   !integer(4), intent(in)  :: nab
   real(8),    intent(in)  :: r, acht, kcht, aab, bab, cab, dab, eab, fab
   real(8),    intent(out) :: ur, urp, urpp

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0, expmax = 87.0d0
   real(8) :: r1i, r2i, r6i, uele, uexp1, uexp2, udis1, udis2, uchtexp
   real(8) :: bri, ud1, ud2, ud3, ud4, ud5, ud6, udd, uddp, uddpp

   r1i    = One/r
   r2i    = r1i**2
   r6i    = r2i**3

   uchtexp = -kcht*exp(-min(expmax,acht*((r-dab)**2)))
   uele   = Zero                          !! 2007-06-05 PL
   uexp1  = eab*exp(-min(expmax,fab*r))
   uexp2  = Zero
   udis1  =-cab*r6i
   udis2  = aab*exp(-min(expmax,bab*r))
!
! new for NemoType7
!
   bri = r1i/bab
   ud1 = 6.0*bri
   ud2 = 5.0*bri*ud1
   ud3 = 4.0*bri*ud2
   ud4 = 3.0*bri*ud3
   ud5 = 2.0*bri*ud4
   ud6 =     bri*ud5
   udd = One + ud1 + ud2 + ud3 + ud4 + ud5 + ud6
   uddp  =-r1i*(    ud1 + 2.0*ud2 +  3.0*ud3 +  4.0*ud4 +  5.0*ud5 +  6.0*ud6)
   uddpp = r2i*(2.0*ud1 + 6.0*ud2 + 12.0*ud3 + 20.0*ud4 + 30.0*ud5 + 42.0*ud6)
!
! modified for NemoType7
!
   ur   = uele + uexp1 + udis1 + udd*udis2 + uchtexp

   urp  =-uele*r1i       - fab*uexp1                                &
        - 6.0*udis1*r1i  - bab*udd*udis2  + uddp*udis2 - acht*dble(2)*(r-dab)*uchtexp

   urpp = 2.0*uele*r2i    + fab*fab*uexp1                &
        + 42.0*udis1*r2i  + bab*bab*udd*udis2 - 2.0*bab*udis2*uddp + uddpp*udis2   &
        + 2.0d0*acht*(2.0d0*acht*((r-dab)**2)-1.0d0)*uchtexp

end subroutine NemoType7

!************************************************************************
!> \page potential potential.F90
!! **WriteUBuffer**
!! *write the content of ubuf for pair iatjat*
!************************************************************************


subroutine WriteUBuffer(iatjat)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iatjat

   integer(4) :: igrid, ibuf, i

   write(uout,'()')
   write(uout,'(a,i5)')     'iatjat            = ', iatjat
   write(uout,'()')
   write(uout,'(a,i5)')     'nugrid            = ', nugrid(iatjat)
   write(uout,'(a,i5)')     'iubuflow          = ', iubuflow(iatjat)
   write(uout,'(a,2g15.5)') 'rumin, rumax      = ', rumin(iatjat), rumax(iatjat)
   write(uout,'(a,2g15.5)') 'tolerance in u, f = ', utoltab, ftoltab
   write(uout,'()')
   write(uout,'(a)')    'igrid    ibuf        rlow        rupp'
   write(uout,'(a)')    '-----    ----        ----        ----'
   igrid = 1
   ibuf = iubuflow(iatjat)
   write(uout,'(a5,i8,4x,8g12.5)') '1', ibuf, sqrt(ubuf(ibuf)), rumax(iatjat), (ubuf(ibuf+i),i = 1,6)
   do igrid = 2, nugrid(iatjat)
      ibuf = ibuf+12
      write(uout,'(i5,i8,4x,8g12.5)') igrid, ibuf, sqrt(ubuf(ibuf)), sqrt(ubuf(ibuf-12)), (ubuf(ibuf+i),i = 1,6)
   end do

end subroutine WriteUBuffer

!************************************************************************
!> \page potential potential.F90
!! **TestUBuffer**
!! *calculate the accuracy of the table for pair iatjat*
!************************************************************************


subroutine TestUBuffer(iat, jat, rlow, rupp, potsub, unit)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iat, jat
   real(8),    intent(in) :: rlow
   real(8),    intent(in) :: rupp
   integer(4), intent(in) :: unit

   character(40), parameter :: txroutine ='UBuffer'
   integer(4), parameter :: ncheck = 101
   integer(4):: ibuf, icheck, idum, iatjat
   real(8)   :: dr, r1, r2, dz
   real(8)   :: u0, u1, u2
   real(8)   :: usum, fsum
   real(8)   :: uerrmax, ferrmax
   external potsub

   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a)') '       r             u              uerr           f              ferr'
   write(unit,'(a)') '       -             -              ----           -              ----'
   if (ilist > 0) then
      write(ulist,'(2a)') txatat(iatat(iat,jat))
      write(ulist,'(i5)') ncheck
   end if

   dr = (rupp-rlow)/(ncheck-1)
   uerrmax = Zero
   ferrmax = Zero

   iatjat = iatat(iat,jat)
   do icheck = ncheck, 1, -1
      r1 = rlow+dr*(icheck-1)
      if (r1 < rumin(iatjat)) cycle
      if (r1 > rumax(iatjat)) cycle
      ibuf = iubuflow(iatjat)
      r2 = r1**2
      do
         if (r2 >= ubuf(ibuf)) exit
         ibuf = ibuf+12
      end do
      dz = r2-ubuf(ibuf)
      usum = ubuf(ibuf+1)+dz*(ubuf(ibuf+2)+dz*(ubuf(ibuf+3)+dz*(ubuf(ibuf+4)+dz*(ubuf(ibuf+5)+dz*ubuf(ibuf+6)))))
      fsum = ubuf(ibuf+7)+dz*(ubuf(ibuf+8)+dz*(ubuf(ibuf+9)+dz*(ubuf(ibuf+10)+dz*ubuf(ibuf+11))))
      fsum = r1*fsum
      call potsub('calc', idum, idum, iat, jat, r1, u0, u1, u2)
      uerrmax = max(uerrmax,abs(usum-u0))
      ferrmax = max(ferrmax,abs(fsum+u1))
      write(unit,'(5g15.5)') r1, u0, usum-u0, -u1, fsum+u1
      if (ilist > 0) write(ulist,'(5g15.5)') r1, u0, usum-u0
   end do

   write(unit,'()')
   write(unit,'(a,g15.5)') 'max error in u   = ', uerrmax
   write(unit,'(a,g15.5)') 'max error in f   = ', ferrmax

end subroutine TestUBuffer

!************************************************************************
!> \page potential potential.F90
!! **PlotPotTwoBodyTab**
!! *plot two-body potential*
!************************************************************************


subroutine PlotPotTwoBodyTab

   use PotentialModule

   character(40), parameter :: txroutine ='PlotPotTwoBodyTab'
   integer(4), parameter :: ncheck = 101
   integer(4) :: ipt, jpt, iat, jat, iatjat, icheck, ibuf
   real(8)    :: dr, rlow, rupp, r1, r2, d, u0, dum(1)
   real(8)    :: distance(1:ncheck), potential(1:ncheck)

   if (slave) return   ! master only

   call WriteHead(3, txroutine, uout)

   if(.not.allocated(lsetatat)) then
      allocate(lsetatat(1:natat))
      lsetatat = .false.
   end if
   lsetatat = .false.

   write(uout,'()')
   write(uout,'(a)') 'plot of atom-atom two-body potential'
   write(uout,'(a)') '------------------------------------'
   do ipt = 1, npt
      do jpt = ipt, npt
         do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
            do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
               iatjat = iatat(iat,jat)
               if (.not.lsetatat(iatjat)) then
                  if (jat < iat) cycle
                  rupp = rumax(iatjat)
                  rlow = rumin(iatjat)
                  dr = (rupp-rlow)/(ncheck-1)
                  do icheck = ncheck, 1, -1
                     r1 = rlow+dr*(icheck-1)
                     r2 = r1**2
                     ibuf = iubuflow(iatjat)
                     do
                        if (r2 >= ubuf(ibuf)) exit
                        ibuf = ibuf+12
                        if (ibuf > nbuf) call Stop(txroutine, 'ibuf > nbuf', uout)
                     end do
                     d = r2-ubuf(ibuf)
                     u0 = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                          d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
                     distance(icheck) = r1
                     potential(icheck) = u0
                     if (potential(icheck) > Two) exit
                  end do
                  if (icheck < 1) icheck = 1
                  rlow = distance(icheck)
                  call Plot(txatat(iatjat), ncheck-icheck+1, potential(icheck), 'none', rlow, rupp, dum(1), dum(1), uout)
                  lsetatat(iatjat) =.true.
               end if
            end do
         end do
      end do
   end do
   write(uout,'()')
   write(uout,'(a)') 'done'
   stop

end subroutine PlotPotTwoBodyTab

!************************************************************************
!> \page potential potential.F90
!! **TestPotTwoBodyTab**
!! *test of utwobodytab and utwobody*
!************************************************************************


subroutine TestPotTwoBodyTab(iStage, unit)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: unit

   character(40), parameter :: txroutine ='TestPotTwoBodyTab'
   integer(4) :: m
   real(8) :: t0, t1
   real(8) :: SecondsSinceStart, tsetnn

   call WriteHead(3, txroutine, unit)

   lmd   =.true.
   lmc   =.false.
   lmcall =.false.

   t0 = SecondsSinceStart()
   call NList(1)
   tsetnn = SecondsSinceStart()-t0
   write(unit,'(a,f12.3)') 'Neighbour list: cpu time ', tsetnn

   write(unit,'()')
   write(unit,'(a,t25,a,t50,a,t70,a)') 'tolerance', 'energy', 'pressure', 'cpu time (s)'
   write(unit,'(a,t25,a,t50,a,t70,a)') '---------', '------', '--------', '------------'
   do m = 7, 3, -1
      utoltab = 10.0d0**(-m)
      ftoltab = 10.0d0**(-m)
      call PotTwoBodyTab1(.false.)
      t0 = SecondsSinceStart()
      if (lintsite) call SetAtomPos(1, np, lintsite)
      call UTotal(iStage)
      if (lintsite) call SetAtomPos(1, np, .false.)
      t1 = SecondsSinceStart()-t0
      write(uout,'(es8.1,t20,f15.8,t45,f15.8,t70,f8.3)') utoltab, u%tot/np, prsr, t1
   end do
   write(unit,'()')
   write(unit,'(a)') 'done'
   stop

end subroutine TestPotTwoBodyTab

!************************************************************************
!> \page potential potential.F90
!! **CalcScrCUfac**
!! *calculate exp(r1/\ref scrlen) and two related derivatives*
!************************************************************************


subroutine CalcScrCUfac(r1, scrlen, c0fac, c1fac, c2fac)

   implicit none

   real(8), intent(in)  :: r1      ! distance
   real(8), intent(in)  :: scrlen  ! screening length
   real(8), intent(out) :: c0fac   ! exp(r1/srclen)
   real(8), intent(out) :: c1fac   ! (-r**2)d[exp()/r]/dr
   real(8), intent(out) :: c2fac   ! (r**3/2)d2[exp()/r]dr2

   real(8), parameter :: One = 1.0d0, Half = 0.50d0
   real(8) :: ra, fac

   ra = r1/scrlen
   fac = exp(-ra)
   c0fac = fac
   c1fac = (One+ra)*fac
   c2fac = (One+ra+Half*ra**2)*fac

end subroutine CalcScrCUfac

!************************************************************************
!> \page potential potential.F90
!! **CalcGCDUfac**
!! *calculate erfc(1/sqrt(1/a1**2+1/a2**2)*r1) and two related derivatives*
!************************************************************************


subroutine CalcGCDUfac(r1, a1, a2, c0fac, c1fac, c2fac)

   implicit none

   real(8), intent(in)  :: r1      ! distance
   real(8), intent(in)  :: a1      ! 1/width of Gaussian charged distribution a1 = 0 => a1 == inf (point ion)
   real(8), intent(in)  :: a2      ! 1/width of Gaussian charged distribution a2 = 0 => a2 == inf (point ion)
   real(8), intent(inout) :: c0fac ! erfc(gamma*r1)
   real(8), intent(inout) :: c1fac ! (-r**2)d[erf()/r]/dr
   real(8), intent(inout) :: c2fac ! (r**3/2)d2[erf()/r]dr2

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0
   real(8), parameter :: fac  = 1.12837916709d0
   real(8) :: gamma, ra, er, ex
   real(8), external :: ErfLocal

   gamma = Zero
   if (a1 > Zero) gamma = gamma + one/a1**2
   if (a2 > Zero) gamma = gamma + one/a2**2
   gamma = one/sqrt(gamma)

   ra = gamma*r1
   er = ErfLocal(ra)
   ex = fac*exp(-ra**2)
   c0fac = c0fac - One + er
   c1fac = c1fac - One + (er-ra*ex)
   c2fac = c2fac - One + (er-ra*(One+ra**2)*ex)

end subroutine CalcGCDUfac

!************************************************************************
!> \page potential potential.F90
!! **CalcEwaldUfac**
!! *calculate erfc(\ref ualpha*r1) and two related derivatives*
!************************************************************************


subroutine CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)

   implicit none

   real(8), intent(in)  :: r1      ! distance
   real(8), intent(in)  :: ualpha  ! damping factor
   real(8), intent(out) :: c0fac   ! erfc(ualpha*r1)
   real(8), intent(out) :: c1fac   ! (-r**2)d[erf()/r]/dr
   real(8), intent(out) :: c2fac   ! (r**3/2)d2[erf()/r]dr2

   real(8), parameter :: One = 1.0d0
   real(8), parameter :: fac  = 1.12837916709d0
   real(8) :: ra, er, ex
   real(8), external :: ErfLocal

   ra = ualpha*r1
   er = ErfLocal(ra)
   ex = fac*exp(-ra**2)
   c0fac = One-(er)
   c1fac = One-(er-ra*ex)
   c2fac = One-(er-ra*(One+ra**2)*ex)

end subroutine CalcEwaldUfac

!************************************************************************
!> \page potential potential.F90
!! **CalcRFUfac**
!! *calculate modified interactions according to nymand and linse, eqs(55)*
!************************************************************************


subroutine CalcRFUfac(r1, rffac, c0fac, c1fac, c2fac)

   implicit none

   real(8), intent(in)  :: r1      ! distance
   real(8), intent(in)  :: rffac   ! reaction-field factor
   real(8), intent(out) :: c0fac   ! (1 + Half*rffac*r^3)
   real(8), intent(out) :: c1fac   ! (1 -      rffac*r^3)
   real(8), intent(out) :: c2fac   ! (1 + Half*rffac*r^3)
   real(8), parameter :: Half = 0.5d0, One = 1.0d0
   real(8) :: r3

   r3 = r1**3
   c0fac = One+Half*rffac*r3
   c1fac = One-rffac*r3
   c2fac = One+Half*rffac*r3

end subroutine CalcRFUfac

!************************************************************************
!> \page potential potential.F90
!! **EwaldErrorUReal**
!! *estimate error of real space potential energy*
!************************************************************************


real(8) function EwaldErrorUReal(lq2sum, q2sum, l, alpha, rcut, unit)

   implicit none

   integer(4), intent(in) :: lq2sum  ! =0 charge, =1 dipole and no charge
   real(8),    intent(in) :: q2sum   ! sum(q**2) (in units making energy in default unit)
   real(8),    intent(in) :: l       ! box length
   real(8),    intent(in) :: alpha   ! ewald parameter
   real(8),    intent(in) :: rcut    ! cutoff in real space
   real(8),    intent(in) :: unit    ! unit factor

   real(8),    parameter :: Pi = 3.14159265359d0
   real(8) :: fac, facbc

   fac = (alpha*rcut)**2
   if (lq2sum == 0) EwaldErrorUReal = unit*q2sum*sqrt(rcut/(2.0d0*l**3))/fac*exp(-fac) ! kolafa and perram
   if (lq2sum == 1) then
!    b = 2*fac+1                                                        !
!    c = fac*(4*fac+6)+3                                                ! exact
!    facbc = (sqrt(15.0d0)/4.0d0)*sqrt((b*b/4 + c*c/15 - b*c/6))/fac**2 !
     facbc = 1.0d0                              ! leading order - consistent with subroutine AlphaToRcut and RcutTo Alpha
     EwaldErrorUReal = unit*q2sum*4.0d0/sqrt(15.0d0*(rcut*l)**3)*facbc*fac*exp(-fac)       ! wang and holm
  end if
end function EwaldErrorUReal

!************************************************************************
!> \page potential potential.F90
!! **EwaldErrorURec**
!! *estimate error of reciprocal space potential energy*
!************************************************************************


real(8) function EwaldErrorURec(lq2sum, q2sum, l, alpha, ncut, unit)

   implicit none

   integer(4), intent(in) :: lq2sum  ! =0 charge, =1 dipole and no charge
   real(8),    intent(in) :: q2sum   ! sum(q**2) (in units makeing energy in default unit)
   real(8),    intent(in) :: l       ! box length
   real(8),    intent(in) :: alpha   ! ewald parameter
   integer(4), intent(in) :: ncut    ! cutoff in reciprocal space
   real(8),    intent(in) :: unit    ! unit factor

   real(8),    parameter :: Pi = 3.14159265359d0
   real(8) :: fac

   fac = (Pi*ncut/(alpha*l))**2
   if (lq2sum == 0) EwaldErrorURec = unit*q2sum*(ncut/l)/fac*exp(-fac)*(1.0+1.0/(alpha*l*sqrt(real(ncut))))             ! kolafa and perram
   if (lq2sum == 1) EwaldErrorURec = unit*q2sum*(4.0d0/3.0d0)*Pi**2*(ncut/l)**3/fac*exp(-fac)*(1.0+1.0/(alpha*l*sqrt(real(ncut)))) ! wang and holm

end function EwaldErrorURec

!************************************************************************
!> \page potential potential.F90
!! **TestEwald**
!! *examine truncation error of ewald summation (energies in kT)*
!************************************************************************


subroutine TestEwald(iStage)

   use MolModule
   integer(4), intent(in) :: iStage
   if (txewaldrec == 'std') call TestEwaldStd(iStage)
   if (txewaldrec == 'spm') call TestEwaldSPM(iStage)

end subroutine TestEwald

!************************************************************************
!> \page potential potential.F90
!! **TestEwaldStd**
!! *examine truncation error of standard ewald summation*
!************************************************************************


subroutine TestEwaldStd(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='TestEwaldStd'
   real(8) :: x, xlow, xupp, dx, fac
   integer(4) :: i,nnstep
   integer(4) :: iewaldoptsave, ncutsave
   real(8) :: uewaldtolsave, ualphasave, ualpharedsave, rcutsave
   real(8) :: uexact

   call WriteHead(3, txroutine, uout)

! ... determine the "exact" electrostatic energy

   iewaldoptsave = iewaldopt
   uewaldtolsave = uewaldtol
   ualphasave = ualpha
   ualpharedsave = ualphared
   rcutsave = rcut
   ncutsave = ncut

  iewaldopt = 2                      ! for obtaining the 'exact' potential energy (maybe weakly system dependent)
  uewaldtol = 1e-12                  !
  ualpha = 10.5/boxlenshort          !

   call EwaldSetup
   call UTotal(iStage)
   call GetPotEnergy(uexact)         ! save for later use

   call WriteExactSub(ustdout)
   call WriteExactSub(uout)

! ... loop over one variable and calculate the error

   iewaldopt = iewaldoptsave
   uewaldtol = uewaldtolsave
   ualpha = ualphasave
   ualphared = ualpharedsave
   rcut = rcutsave
   ncut = ncutsave

   call WriteHeadSub(ustdout)
   call WriteHeadSub(uout)

   if (iewaldopt == 0) then
      x = ualphared
      fac = 4
      nnstep = 50
   else if (iewaldopt == 1) then
      x = rcut
      fac = 4
      nnstep = 50
   else if (iewaldopt == 2) then
      x = ualpha
      fac = 2
      nnstep = 40
   end if
!  fac = 8           ! older versions
!  nnstep = 23       ! older versions
   xlow = x
   xupp = fac*x
   dx = (xupp-xlow)/(nnstep-1)

   do i = 0, nnstep-1
      x = xlow + i*dx
      if (iewaldopt == 0) ualphared = x
      if (iewaldopt == 1) rcut = x
      if (iewaldopt == 2) ualpha = x
      call EwaldSetup
      if (rcut > boxlenshort) cycle
      call UTotal(iStage)
      call WriteBodySub(ustdout)
      call WriteBodySub(uout)
      call WriteBodySub(ulist)
   end do

   stop

contains

!........................................................................

subroutine GetPotEnergy(uuu)
   real(8), intent(out) :: uuu
   if(txelec == 'charge') uuu = u%tot         ! total energy
   if(txelec == 'dip') uuu = u%stat           ! electrostatic energy
end subroutine GetPotEnergy

subroutine WriteExactSub(unit)
   integer(4), intent(in) :: unit
   real(8) :: EwaldErrorUReal, EwaldErrorURec
   write(unit,'(a)')     '''exact'' calculation'
   write(unit,'(a)')     '-------------------'
   write(unit,'(a,t35,3g15.5)') 'ualpha                         = ', ualpha
   write(unit,'(a,t35,3g15.5)') 'ualphared                      = ', ualphared
   write(unit,'(a,t35,3g15.5)') 'rcut                           = ', rcut
   write(unit,'(a,t35,3g15.5)') 'ncut                           = ', ncut
   write(unit,'(a,t35,3g20.10)') 'u_exact/np                     = ', uexact/np
   write(unit,'(a,t35,3g15.5)') 'error estimate in real space   = ', EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi)
   write(unit,'(a,t35,3g15.5)') 'error estimate in rec. space   = ', EwaldErrorURec(lq2sum,q2sum,boxlenshort,ualpha,ncut,EpsiFourPi)
   write(unit,*)
end subroutine WriteExactSub

!........................................................................

subroutine WriteHeadSub(unit)
   integer(4), intent(in) :: unit
   write(unit,'(10(a,2x))')   &
  'iewaldopt', 'uewaldtol', '   ualphared', '       ualpha  ', '      rcut ', '       ncut', '        u_approx/np', &
  '       ErrorUReal', '    ErrorURec', ' |(u_approx-u_exact)/np|'
end subroutine WriteHeadSub

!........................................................................

subroutine WriteBodySub(unit)
   integer(4), intent(in) :: unit
   real(8) :: EwaldErrorUReal, EwaldErrorURec, uapprox
   call GetPotEnergy(uapprox)
   write(unit,'((i4,a),(g12.5,a),3(f10.5,a),(i4,a),(f20.12,a),3(g12.5,a))') &
   iewaldopt, tab, uewaldtol, tab, ualphared, tab, ualpha, tab,  &
   rcut, tab, ncut, tab, uapprox/np, tab, &
   EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi), tab, &
   EwaldErrorURec(lq2sum,q2sum,boxlenshort,ualpha,ncut,EpsiFourPi), tab, &
   abs(uapprox-uexact)/np
end subroutine WriteBodySub

!........................................................................

end subroutine TestEwaldStd

!************************************************************************
!> \page potential potential.F90
!! **TestEwaldSPM**
!! *examine truncation error of the reciprocal space; smooth particle mesh ewald summation*
!************************************************************************


subroutine TestEwaldSPM(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='TestEwaldSPM'
   integer(4), parameter :: n_order_loop = 6, n_nmesh_loop = 5
   integer(4), save :: iSamp                                     ! sample the number of samplings
   integer(4), save :: order_loop(2) =[4, n_order_loop]          ! indices of order loop
   integer(4), save :: nmesh_loop(2) =[16, n_nmesh_loop]         ! indices of nmesh loop
   real(8), save :: urecsum2(n_order_loop,n_nmesh_loop), tsum(n_order_loop,n_nmesh_loop)
   character(4) :: txewaldrecsave
   integer(4) :: iorder, inmesh
   integer(4) :: iewaldoptsave, ordersave, nmeshsave
   real(8) :: uewaldtolsave, ualphasave, ualpharedsave, rcutsave
   real(8) :: Gettime_ewald
   real(8) :: uexact
   character(2) :: str

   if(master) call WriteHead(3, txroutine, uout)

   select case (iStage)
   case (iReadInput)

! .............. initialize ..............

      iSamp = 0
      urecsum2 = zero

   case (iAfterMacrostep)

! .............. calculate electrostatic energies ..............

      iSamp = iSamp + 1

      txewaldrecsave = txewaldrec
      iewaldoptsave = iewaldopt
      uewaldtolsave = uewaldtol
      ualphasave = ualpha
      ualpharedsave = ualphared
      rcutsave = rcut
      ordersave = order
      nmeshsave = nmesh

! ... determine the "exact" total and reciprocal electrostatic energies

      txewaldrec ='std'
      iewaldopt = 1
      ncut = 25

      txewaldrec ='spm'
      iewaldopt = 0

      ualphared = ualpharedsave
      rcut = rcutsave
      order = 11
      nmesh = 256

      call EwaldSetup
      call UTotal(iStage)
      uexact = u%rec/np       ! save "exact" reciprocal energy

     if(master) call WriteExactSub(uout)

! ... loop over order and nmesh, and write the error of the reciprocal energies

      txewaldrec = txewaldrecsave
      iewaldopt = iewaldoptsave
      uewaldtol = uewaldtolsave
      ualpha = ualphasave
      ualphared = ualpharedsave
      rcut = rcutsave
      order = ordersave
      nmesh = nmeshsave

      if(master) call WriteHeadSub(uout)

      do iorder = 1, order_loop(2)
         order = order_loop(1) + (iorder-1)
         do inmesh = 1, nmesh_loop(2)
            nmesh = nmesh_loop(1)*2**(inmesh-1)
            call EwaldSetup
            if (rcut > minval(boxlen2(1:3))) cycle
            call UTotal(iStage)
            if(master) call WriteBodySub(uout)
            urecsum2(iorder,inmesh) = urecsum2(iorder,inmesh) + (u%rec/np-uexact)**2 ! sample error of reciprocal energy
            tsum(iorder,inmesh) = tsum(iorder,inmesh) + Gettime_ewald()              ! sample cpu time
         end do
      end do

      txewaldrec = txewaldrecsave
      iewaldopt = iewaldoptsave
      uewaldtol = uewaldtolsave
      ualpha = ualphasave
      ualphared = ualpharedsave
      rcut = rcutsave
      order = ordersave
      nmesh = nmeshsave

   case (iAfterSimulation)

! ............... form ratios of rms energies .............

      if(master) then
         do iorder = 1, order_loop(2)
            order = order_loop(1) + (iorder-1)
            write(str,'(i1)') order
            open(unit=30, file = 'order_'//str)
            do inmesh = 1, nmesh_loop(2)
               nmesh = nmesh_loop(1)*2**(inmesh-1)
               urecsum2(iorder,inmesh) = urecsum2(iorder,inmesh)/iSamp
               tsum(iorder,inmesh) = tsum(iorder,inmesh)/nstep
               write(30,'(3(i7,a),3(g12.4,a),2g15.5)') &
                  na, tab, order, tab, nmesh, tab, &
                  nmesh/(float(na)**Third),tab, real(na*order**3), tab, nmesh**3*log(float(nmesh**3)), tab, &
                  sqrt(urecsum2(iorder,inmesh)), tsum(iorder,inmesh)/nstep
            end do
         end do
      end if

   end select

contains

!........................................................................

subroutine WriteExactSub(unit)
   integer(4), intent(in) :: unit
   real(8) :: EwaldErrorUReal
   write(unit,'(a)')     '''exact'' calculation'
   write(unit,'(a)')     '-------------------'
   write(unit,'(a,t35,a)')      'txewaldrec                     = ', txewaldrec
   write(unit,'(a,t35,3g15.5)') 'ualpha                         = ', ualpha
   write(unit,'(a,t35,3g15.5)') 'ualphared                      = ', ualphared
   write(unit,'(a,t35,3g15.5)') 'rcut                           = ', rcut
   write(unit,'(a,t35,3g15.5)') 'ncut                           = ', ncut
   write(unit,'(a,t35,3g15.5)') 'order                          = ', order
   write(unit,'(a,t35,3g15.5)') 'nmesh                          = ', nmesh
   write(unit,'(a,t35,3g20.10)') 'u%rec ''exact''                 = ', uexact
   write(unit,'(a,t35,3g15.5)') 'error estimate in real space   = ', EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi)
   write(unit,*)
end subroutine WriteExactSub

!........................................................................

subroutine WriteHeadSub(unit)
   integer(4), intent(in) :: unit
   write(unit,'(10(a,2x))')   &
  'iewaldopt', 'uewaldtol', '   ualphared', '       ualpha  ', '      rcut  ', '   order', '  nmesh', '        u%rec/np', &
  '       ErrorUReal', '  u%rec/np-uexact'
end subroutine WriteHeadSub

!........................................................................

subroutine WriteBodySub(unit)
   integer(4), intent(in) :: unit
   real(8) :: EwaldErrorUReal
   write(unit,'((i4,a),(g12.5,a),3(f10.5,a),2(i4,a),(f20.12,a),2(g12.5,a))') &
   iewaldopt, tab, uewaldtol, tab, ualphared, tab, ualpha, tab,  &
   rcut, tab, order, tab, nmesh, tab, u%rec/np, tab, &
   EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi), tab, u%rec/np-uexact
end subroutine WriteBodySub

!........................................................................

end subroutine TestEwaldSPM


!************************************************************************
!> \page potential potential.F90
!! **SetImageSph**
!! *Setup of image charges and dipoles according to Friedman*
!************************************************************************

!     Mol. Phys. 1975

subroutine SetImageSph(iplow, ipupp, mode)

   use MolModule
   implicit none

   integer(4), intent(in) :: mode
   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp

   character(40), parameter :: txroutine ='SetImageSph'
   integer(4)             :: ip, iploc
   logical, save          :: first = .true.
   real(8), save          :: rad2img, zfac

   if (first) then
      if (.not.allocated(zimg)) then
         allocate(zimg(na_alloc))
         zimg = 0.0E+00
      end if
      if (.not.allocated(rimg)) then
         allocate(rimg(3,na_alloc))
         rimg = 0.0E+00
      end if
      if (.not.allocated(dipimg)) then
         allocate(dipimg(3,na_alloc))
         dipimg = 0.0E+00
      end if
      rad2img = radimg**2
      zfac = ((epsimg-1)/(epsimg+1))*radimg
      first = .false.
   end if

   if (mode == 1) then        ! set image particles from iplow to ipupp
      do ip = iplow, ipupp
         call SetImageSphSub(rad2img, zfac, ro(1:3,ip), az(ip), dip(1:3,ip), rimg(1:3,ip), zimg(ip), dipimg(1:3,ip))
      end do
   else if (mode == 2) then     ! update image particles for moved particles
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         call SetImageSphSub(rad2img, zfac, rotm(1:3,iploc), az(ip), diptm(1:3,iploc), rimg(1:3,ip), zimg(ip), dipimg(1:3,ip))
      end do
   else if (mode == 3) then     ! update moved image particles to old configuration
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         call SetImageSphSub(rad2img, zfac, ro(1:3,ip), az(ip), dip(1:3,ip), rimg(1:3,ip), zimg(ip), dipimg(1:3,ip))
      end do
   else
      call Stop(txroutine,'Error in mode',uout)
   end if

contains

!........................................................................

subroutine SetImageSphSub(rad2img, zfac, ro, az, dip, rimg, zimg, dipimg)
   real(8), intent(in)  :: rad2img
   real(8), intent(in)  :: zfac
   real(8), intent(in)  :: ro(3)
   real(8), intent(in)  :: az
   real(8), intent(in)  :: dip(3)
   real(8), intent(out) :: rimg(3)
   real(8), intent(out) :: zimg
   real(8), intent(out) :: dipimg(3)
   real(8) :: rfac, dipfac, r2, r3i, r2i, r1i, dipm2, ac, angcos, rot(3,3)
   r2 = ro(1)**2 + ro(2)**2 + ro(3)**2
   r2i = One/r2
   r1i = sqrt(r2i)
   r3i = r1i*r2i
   rfac = rad2img*r2i
   rimg(1:3) = ro(1:3)*rfac
   dipm2 = dip(1)**2 + dip(2)**2 + dip(3)**2
   if (dipm2 < 1d-12) then
      zimg = -zfac*az*r1i
      dipimg(1:3) = Zero
   else
      ac = angcos(dip(1),dip(2),dip(3),ro(1),ro(2),ro(3))
      zimg = zfac*(-az*r1i + r2i*ac*sqrt(dipm2))
      dipfac = zfac*rad2img*r3i
      call AxisAngToOri(ro(1:3), Pi, rot)
      dipimg(1) = dipfac*(dip(1)*rot(1,1) + dip(2)*rot(1,2) + dip(3)*rot(1,3))
      dipimg(2) = dipfac*(dip(1)*rot(2,1) + dip(2)*rot(2,2) + dip(3)*rot(2,3))
      dipimg(3) = dipfac*(dip(1)*rot(3,1) + dip(2)*rot(3,2) + dip(3)*rot(3,3))
   endif
end subroutine SetImageSphSub

!........................................................................

end subroutine SetImageSph

!************************************************************************
!> \page potential potential.F90
!! **IOPotChain**
!! *perform i/o on bond and angle interaction variables*
!************************************************************************


!> \page nmlPotentialChain
!! The namelist  \ref nmlPotentialChain contains variables that describe the potentials involving particles belonging to the same chain.
!! There is a bond potential between two consecutive particles in a chain and an angular potential between three consecutive particles
!! in a chain. The potentials are of the type (k/p)(x-x0)**p.
!! * Variables:
!!  * \subpage bond
!!  * \subpage angle
!!  * \subpage clink
!!  * \subpage itestpotchain

subroutine IOPotChain(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOPotChain'
   character(80), parameter :: txheading ='ipotchain data'

   namelist /nmlPotentialChain/ bond, angle, clink, itestpotchain

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(bond)) then
         allocate(bond(nct))
      end if
      if (.not.allocated(angle)) then
         allocate(angle(nct))
      end if

! ... set initial values

      bond  = bond_var(Zero, 2, Zero)
      angle = bond_var(Zero, 2, 180.0d0)
      clink = bond_var(Zero, 2, Zero)
      itestpotchain = 0

! ... read input data

      rewind(uin)
      read(uin,nmlPotentialChain)

! ... check some conditions

      if (count(bond(1:nct)%k < Zero) > 0)  call Stop(txroutine, 'bond%k < 0', uout)
      if (count(bond(1:nct)%p < 1) > 0)     call Stop(txroutine, 'bond%p < 1', uout)
      if (count(bond(1:nct)%eq < Zero) > 0) call Stop(txroutine, 'bond%eq < 0', uout)
      if (count(angle(1:nct)%k < Zero) > 0) call Stop(txroutine, 'angle%k < 0', uout)
      if (count(angle(1:nct)%p < 1) > 0)    call Stop(txroutine, 'angle%p < 1', uout)
      if (clink%k < Zero)  call Stop(txroutine, 'clink%k < 0', uout)
      if (clink%p < 1)     call Stop(txroutine, 'clink%p < 1', uout)
      if (clink%eq < Zero) call Stop(txroutine, 'clink%eq < 0', uout)

   case (iWriteInput)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a,6f15.5)') 'bond: force constant           = ', bond(1:nct)%k
         write(uout,'(a,6g15.5)') 'bond: power                    = ', bond(1:nct)%p
         write(uout,'(a,6f15.5)') 'bond: eqilibrium separation    = ', bond(1:nct)%eq
         write(uout,'(a,6f15.5)') 'angle: force constant          = ', angle(1:nct)%k
         write(uout,'(a,6g15.5)') 'angle: power                   = ', angle(1:nct)%p
         write(uout,'(a,6f15.5)') 'angle: eqilibrium separation   = ', angle(1:nct)%eq
         if (lclink) write(uout,'(a,2f15.5)') 'crosslink: force constant      = ', clink%k
         if (lclink) write(uout,'(a,2g15.5)') 'crosslink: power               = ', clink%p
         if (lclink) write(uout,'(a,2f15.5)') 'crosslink: eqilibrium separ    = ', clink%eq
         write(uout,'(a,2g15.5)') 'itestpotchain                  = ', itestpotchain
      end if

      angle(1:nct)%k = angle(1:nct)%k/sclang**2     ! angle force constant in rad
      angle(1:nct)%eq = angle(1:nct)%eq*sclang      ! equilibrium angle in rad

      bond(1:nct)%k = bond(1:nct)%k/bond(1:nct)%p   ! divide force konstant by varible for fewer operations in energy evaluations
      angle(1:nct)%k = angle(1:nct)%k/angle(1:nct)%p
      clink%k = clink%k/clink%p

   end select

end subroutine IOPotChain

!************************************************************************
!> \page potential potential.F90
!! **ChainTab**
!! *call routines setting up potential tables for chains*
!************************************************************************


subroutine ChainTab

   use PotentialModule
   implicit none

   integer(4) :: idum, iseedsave
   real(8)    :: dum

   call BondLengthTab('setup', idum, dum)
   call BondAngleTab('setup', idum, dum)

   if (itestpotchain == 1) then
      iseedsave = iseed
      call TestChainTab(nct,uout)
      iseed = iseedsave
   end if

contains

!........................................................................

! ... test BondLengthTab and BondAngleTab

subroutine TestChainTab(nct, uout)

   implicit none

   integer(4), intent(in) :: nct   ! number of chain types
   integer(4), intent(in) :: uout  ! output unit

   integer(4) :: n, ict, i
   real(8)    :: x, xsum, xsum2

   write(uout,'()')
   write(uout,'(a)') 'quantity    chain type      <x>      <x**2>**0.5  sqrt(<x**2>-<x>**2)'

   n = 100000

   do ict = 1, nct
      xsum = 0.0d0
      xsum2 = 0.0d0
      do i = 1, n
         call BondLengthTab('value', ict, x)
         xsum = xsum + x
         xsum2 = xsum2 + x**2
      end do
      write(uout,'(a,i6,3f15.4)') 'bond length', ict, xsum/n, sqrt(xsum2/n), sqrt(xsum2/n-(xsum/n)**2)
   end do

   do ict = 1, nct
      xsum = 0.0d0
      xsum2 = 0.0d0
      do i = 1, n
         call BondAngleTab('value', ict, x)
         xsum = xsum + x
         xsum2 = xsum2 + x**2
      end do
      write(uout,'(a,i6,3f15.4)') 'bond angle ', ict, xsum/n/sclang, sqrt(xsum2/n)/sclang, sqrt(xsum2/n-(xsum/n)**2)/sclang
   end do

end subroutine TestChainTab

!........................................................................

end subroutine ChainTab

!************************************************************************
!> \page potential potential.F90
!! **BondLengthTab**
!! *make and use of lookup table to get random bond lengths*
!************************************************************************


!     p(b) = const * exp(-bond_k*(b-bond_eq)**2)*b**bond_p
!     ptab(b) = integ(0_to_b) p(b') db'

subroutine BondLengthTab(action, ictx, b)

   use PotentialModule
   implicit none

   character(5), intent(in)    :: action           ! 'setup' or 'value'
   integer(4),   intent(in) :: ictx             ! chain type
   real(8),      intent(out)   :: b                ! random bond length

   character(40), parameter :: txroutine ='BondLengthTab'
   integer(4)   , parameter :: nbin = 100          ! number of bins of the potential tables btab and ptab
   real(8), allocatable, save :: btab(:,:)         ! table of bond lengths
   real(8), allocatable, save :: ptab(:,:)         ! integrated probability of corresponding bond length
   real(8)       :: bond_k                         ! force constant/kT
   integer(4)    :: bond_p                         ! power of bond length
   real(8)       :: bond_eq                        ! equilibrium bond length
   real(8)       :: blow                           ! lower end of table
   real(8)       :: bupp                           ! upper end of table

   integer(4) :: ict, i
   real(8)    :: db, p, Random, VInter

   if (action == 'setup') then

      if(.not.allocated(btab)) then
         allocate(btab(0:nbin,nct))
         btab = 0.0E+00
      end if
      if(.not.allocated(ptab)) then
         allocate(ptab(0:nbin,nct))
         ptab = 0.0E+00
      end if

      do ict = 1, nct

         bond_k  = bond(ict)%k*beta
         bond_p  = bond(ict)%p
         bond_eq = bond(ict)%eq

         if (bond_k <= 0.0d0) call Stop(txroutine, 'bond_k <= 0', 6)

         blow = max(0.0d0,bond_eq-sqrt(30.0d0/bond_k))
         bupp = bond_eq+sqrt(30.0d0/bond_k)
         db = (bupp-blow)/nbin

         btab(0,ict) = blow
         ptab(0,ict) = 0.0d0
         do i = 1, nbin
            b = blow+i*db                                       ! calculate b
            btab(i,ict) = b                                     ! store b(i)
            b = b-0.5d0*db                                      ! shift b to get right p increment
            p = exp(-bond_k*(b-bond_eq)**2)*b**bond_p           ! calculate p
            ptab(i,ict) = ptab(i-1,ict)+p                       ! store ptab(i)
         end do
         ptab(0:nbin,ict) = ptab(0:nbin,ict)/ptab(nbin,ict)

      end do

      if (itestpotchain == 1) call TestBondLengthTab(uout)

   else if (action == 'value') then

      b = VInter(nbin+1,ptab(0,ictx),btab(0,ictx),Random(iseed))

   else

      call Stop(txroutine, 'action out of range', 6)

   end if

contains

!........................................................................

subroutine TestBondLengthTab(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(10a)') (' bond distance    probability   ', ict = 1, nct)
   do i = 0, nbin
      write(unit,'(5(2g15.5,2x))') (btab(i,ict),ptab(i,ict),ict = 1,nct)
   end do
end subroutine TestBondLengthTab

!........................................................................

end subroutine BondLengthTab

!************************************************************************
!> \page potential potential.F90
!! **BondAngleTab**
!! *make and use of lookup table to get random bond angles*
!************************************************************************


!     p(t) = const * exp(-angle_k*(t-angle_eq)**2)*sin(t)
!     ptab(t) = integ(0_to_t) p(t') dt'

subroutine BondAngleTab(action, ictx, a)

   use PotentialModule
   implicit none

   character(5), intent(in)    :: action           ! 'setup' or 'value'
   integer(4),   intent(in) :: ictx             ! chain type
   real(8),      intent(out)   :: a                ! random bond angle

   character(40), parameter :: txroutine ='BondAngleTab'
   integer(4),    parameter :: nbin = 1800         ! number of bins of the potential tables btab and ptab
   real(8), allocatable, save :: atab(:,:)         ! table of bond angles
   real(8), allocatable, save :: ptab(:,:)         ! integrated probability of corresponding bond angle
   real(8)       :: angle_k                        ! force constant/kT
   integer(4)    :: angle_p                        ! power of bond angle
   real(8)       :: angle_eq                       ! equilibrium bond angle
   real(8)       :: alow                           ! lower end of table
   real(8)       :: aupp                           ! upper end of table

   integer(4) :: ict, i
   real(8)    :: da, p, Random, VInter

   if (action == 'setup') then

      if(.not.allocated(atab)) then
         allocate(atab(0:nbin,nct))
         atab = 0.0E+00
      end if
      if(.not.allocated(ptab)) then
         allocate(ptab(0:nbin,nct))
         ptab = 0.0E+00
      end if

      do ict = 1, nct

         angle_k  = angle(ict)%k*beta
         angle_p  = angle(ict)%p
         angle_eq = angle(ict)%eq
         alow = Zero
         aupp = Pi
         if (angle_k > Zero) alow = max(alow,angle_eq-sqrt(30.0d0/angle_k))
         if (angle_k > Zero) aupp = min(aupp,angle_eq+sqrt(30.0d0/angle_k))
         da = (aupp-alow)/nbin

         atab(0,ict) = alow
         ptab(0,ict) = Zero
         do i = 1, nbin
            a = alow+i*da                                        ! calculate a
            atab(i,ict) = a                                      ! store a(i)

            a = a-Half*da                                        ! shift a to get right p increment
            p = exp(-angle_k*(a-angle_eq)**angle_p)*sin(a)       ! calculate p
            ptab(i,ict) = ptab(i-1,ict)+p                        ! store ptab(i)
         end do
         ptab(0:nbin,ict) = ptab(0:nbin,ict)/ptab(nbin,ict)

      end do

      if (itestpotchain == 1) call TestBondAngleTab(uout)

   else if (action == 'value') then

      a = VInter(nbin+1,ptab(0,ictx),atab(0,ictx),Random(iseed))

   else

      call Stop(txroutine, 'action out of range', 6)

   end if

contains

!........................................................................

subroutine TestBondAngleTab(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(uout,'(20a)') (' bond angle       probability   ', ict = 1, nct)
   do i = 0, nbin
      write(unit,'(5(2g15.5,2x))') (atab(i,ict)/sclang,ptab(i,ict),ict = 1,nct)
   end do
end subroutine TestBondAngleTab

!........................................................................

end subroutine BondAngleTab

!************************************************************************
!> \page potential potential.F90
!! **IOPotExternal**
!! *perform i/o on external potential variables*
!************************************************************************

!> \page nmlPotentialExternal
!! The namelist  \ref nmlPotentialExternal contains variables that describe the potentials involving atoms and an external potential.
!! * Variables:
!!  * \subpage txuext
subroutine IOPotExternal(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOPotExternal'
   character(80), parameter :: txheading ='external potential data'
   integer(4) :: ipt, iat

   namelist /nmlPotentialExternal/ txuext,                                         &
                                   wall_z_ext,                                     &
                                   gravitation_force,                              &
                                   lambda_ramp_ext, epsilon_ramp_ext,              &
                                   lambda_sw_ext, epsilon_sw_ext, alpha_sw_ext,    &
                                   sigma_ext, epsilon_ext,                         &
                                   z3coeff_ext, z9coeff_ext, u1_drift,             &
                                   c_ext, lx_ext, ly_ext,                          &
                                   efield_ext,                                     &
                                   surfchargeden, llongrangecontr, zdist, chden,   &
                                   rInsSphere, zInsSphere, rChargeIn, rChargeOut, rInSphere, rOutSphere, &
                                   ruext, auext, nuext,                            &
                                   rcap, dcap,                                     &
                                   boundaryrad, epsi1, epsi2, lmaxdiel,            &
                                   ofstrength, ofaxis,                             &
                                   rCylinder, zCylinder

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(sigma_ext)) then
         allocate(sigma_ext(nat), epsilon_ext(nat), z3coeff_ext(nat), z9coeff_ext(nat))
         sigma_ext = 0.0E+00
         epsilon_ext = 0.0E+00
         z3coeff_ext = 0.0E+00
         z9coeff_ext = 0.0E+00
      end if
      if (.not.allocated(txuext)) then
         allocate(txuext(npt), ruext(3,npt), ruexti(3,npt))
         txuext = ""
         ruext = 0.0E+00
         ruexti = 0.0E+00
      end if

! ... set initial values

      txuext      = ''
      wall_z_ext  = Zero
      ruext       = Zero
      sigma_ext   = Zero
      epsilon_ext = Zero
      z3coeff_ext = Zero
      z9coeff_ext = Zero
      u1_drift    = Zero
      c_ext       = 0.0d0
      lx_ext      = 5.0d0
      ly_ext      = 5.0d0
      rInsSphere  = 0.0d0
      zInsSphere  = 0.0d0
      rInSphere   = 0.0d0
      rOutSphere  = 0.0d0
      efield_ext  = [Zero, Zero, Zero]
      epsilon_ramp_ext = Two*epsilon_ramp   ! thus assuming equal matter
      lambda_ramp_ext  = lambda_ramp        ! thus assuming equal matter
      alpha_sw_ext= 1.0d-3                  ! steepness of the rise (lower alpha => steeper rise)
      epsi1       = 1
      epsi2       = 2
      lmaxdiel    = 10
      lbg         = .false.
      ofstrength  = Zero
      ofaxis(1:3) = [Zero, Zero, One]
      rCylinder   = Zero
      zCylinder   = Zero                    ! number of charges per unit length

! ... read input data

      rewind(uin)
      read(uin,nmlPotentialExternal)

      do ipt = 1, npt
         call LowerCase(txuext(ipt))
      end do

! ... check some conditions

      if (wall_z_ext > boxlen2(3)) call stop(txroutine, 'wall_z_ext > boxlen2(3)', uout)

      do ipt = 1, npt
         if (txuext(ipt) == 'i_soft_sphere') then
            if (nuext < 2) call Stop(txroutine, 'nuext < 2', uout)
            if (1.2*ruext(1,ipt) > sphrad) call Warn(txroutine, '1.2*ruext(1,ipt) > sphrad', uout)
         end if
         if (txuext(ipt) == 'estat_field_z') then
            if (txbc /= 'sph') call Warn(txroutine, 'txbc /= ''sph''', uout)
         end if
         if (txuext(ipt) =='core_shell') then
            if (txbc /= 'sph') call Warn(txroutine, 'txbc /= ''sph''', uout)
            if (rChargeOut > sphrad) call Stop(txroutine, 'rChargeOut > simulation cell', uout)
            if (rChargeIn > rChargeOut) call Stop(txroutine, 'rChargeIn > rChargeOut', uout)
         end if
         if (txuext(ipt) =='insulating_sphere') then
            if (txbc /= 'sph') call Warn(txroutine, 'txbc /= ''sph''', uout)
            if (rInsSphere > sphrad) call Stop(txroutine, 'rInsSphere > simulation cell', uout)
         end if
         if (txuext(ipt) =='hollow_sphere') then
            if (txbc /= 'sph') call Warn(txroutine, 'txbc /= ''sph''', uout)
            if (rInSphere > sphrad) call Stop(txroutine, 'rInSphere > simulation cell', uout)
            else if (txuext(ipt) == 'orienting_field') then
               write(uout,'(a)') 'orienting field in z-direction'
               write(uout,'(a,t35,3g10.3)') 'strength                       = ', ofstrength
               write(uout,'(a,t35,3g10.3)') 'axis                           = ', ofaxis(1:3)
         end if
         if (txuext(ipt)=='hard_cylinder') then
            if (txbc == 'cyl') then
               if (rCylinder > cylrad) call Stop(txroutine, 'rCylinder > cell radius', uout)
            endif
         endif
      end do

! ... handle simplificaion

      where (ruext(2,1:npt) == Zero) ruext(2,1:npt) = ruext(1,1:npt)
      where (ruext(3,1:npt) == Zero) ruext(3,1:npt) = ruext(1,1:npt)

! ... set ruexti

       where (ruext(1:3,1:npt) > Zero) ruexti(1:3,1:npt) = One/ruext(1:3,1:npt)

   case (iWriteInput)

      do ipt = 1, npt
         iat = ipt
         if (z3coeff_ext(iat) == Zero .and. z9coeff_ext(iat) == Zero) then
            z3coeff_ext(iat) = -(Two*Pi/Three*epsilon_ext(iat))*sigma_ext(iat)**3
            z9coeff_ext(iat) = +(Two*Pi/Three*epsilon_ext(iat))*(Two/15.0d0)*sigma_ext(iat)**9
         else if (sigma_ext(iat) == Zero .and. epsilon_ext(iat) == Zero) then
            sigma_ext(iat) = (7.5d0*z9coeff_ext(iat)/(-z3coeff_ext(iat)))**(1.0d0/6.0d0)
            epsilon_ext(iat)  = Three/(Two*Pi)*(-z3coeff_ext(iat))/sigma_ext(iat)**3
         endif
      end do

      if(.not.allocated(zmin_ext)) then
         allocate(zmin_ext(nat), delta_ext(nat))
         zmin_ext = 0.0E+00
         delta_ext = 0.0E+00
      end if

      do ipt = 1, npt
         if (txuext(ipt) == 'lj_wall_z_mod') then
            iat = ipt
            zmin_ext(iat) = 0.4d0**(1.0d0/6.0d0)*sigma_ext(iat)
            delta_ext(iat) = Two*Pi*sqrt(10.0d0)/9.0d0*epsilon_ext(iat)
            TwoPilxi_ext = TwoPi/lx_ext
            TwoPilyi_ext = TwoPi/ly_ext
         else if (txuext(ipt) == 'lj_wall_z_ts') then
            iat = ipt
            zmin_ext(iat) = 0.4d0**(1.0d0/6.0d0)*sigma_ext(iat)
            delta_ext(iat) = Two*Pi*sqrt(10.0d0)/9.0d0*epsilon_ext(iat)
         end if
      end do

! ... write input data

      if (master) then
         call WriteHead(2, txheading, uout)
         do ipt = 1, npt
            write(uout,'(a,t35,a)')      'particle type                  = ', txpt(ipt)
            if (txuext(ipt) == 'wall_z') then
               write(uout,'(a)') 'hard walls at z = -wall_z_ext and +wall_z_ext'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
            else if (txuext(ipt) == 'gravitation_wall_z') then
               write(uout,'(a)') 'hard wall + gravitation potential'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,3g10.3)') 'gravitation force              = ', gravitation_force
            else if (txuext(ipt) == 'superball_wall_z') then
               write(uout,'(a)') 'hard wall + superball nonpenetrating potential'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
            else if (txuext(ipt) == 'laura') then
               write(uout,'(a)') 'hard wall + superball nonpenetrating potential, gravitation, and estat field'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,3g10.3)') 'gravitation force              = ', gravitation_force
               write(uout,'(a,t35,3g10.3)') 'efield_ext                     = ', efield_ext(1:3)
            else if (txuext(ipt) == 'ramp_wall_z') then
               write(uout,'(a)') 'hard wall + ramp potential extending to lambda_ramp_ext*r1atat'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,3g10.3)') 'lambda_ramp_ext                = ', lambda_ramp_ext
               write(uout,'(a,t35,3g10.3)') 'epsilon_ramp_ext               = ', epsilon_ramp_ext
            else if (txuext(ipt) == 'sw_wall_zlow') then
               write(uout,'(a)') 'hard wall + sw potential extending to lambda_sw_ext*r1atat'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,5g10.3)') 'lambda_sw_ext                  = ', lambda_sw_ext
               write(uout,'(a,t35,5g10.3)') 'epsilon_sw_ext                 = ', epsilon_sw_ext
               write(uout,'(a,t35,5g10.3)') 'alpha_sw_ext                   = ', alpha_sw_ext
            else if (txuext(ipt) == 'lj_wall_z') then
               write(uout,'(a)') 'lj wall at z = +-wall_z_ext'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,5g10.3)') 'sigma_ext       (atom type)    = ', sigma_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'epsilon_ext     (atom type)    = ', epsilon_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z3 coefficients (atom type)    = ', z3coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z9 coefficients (atom type)    = ', z9coeff_ext(1:nat)
            else if (txuext(ipt) == 'lj_wall_z_ts') then
               write(uout,'(a)') 'lj wall at z = +-wall_z_ext, truncated and shifted to zero at minimum.'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,5g10.3)') 'sigma_ext       (atom type)    = ', sigma_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'epsilon_ext     (atom type)    = ', epsilon_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z3 coefficients (atom type)    = ', z3coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z9 coefficients (atom type)    = ', z9coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'cutoff distance (atom type)    = ', zmin_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'energy shift    (atom type)    = ', delta_ext(1:nat)
            else if (txuext(ipt) == 'lj_wall_z_mod') then
               write(uout,'(a)') 'modulated lj wall at z = +-wall_z_ext'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,5g10.3)') 'sigma_ext       (atom type)    = ', sigma_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'epsilon_ext     (atom type)    = ', epsilon_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z3 coefficients (atom type)    = ', z3coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z9 coefficients (atom type)    = ', z9coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'c_ext                          = ', c_ext
               write(uout,'(a,t35,5g10.3)') 'lx_ext                         = ', lx_ext
               write(uout,'(a,t35,5g10.3)') 'ly_ext                         = ', ly_ext
            else if (txuext(ipt) == 'lj_wall_zlow') then
               u0_drift = -u1_drift*(-boxlen2(3))                  ! Zero contribution at -boxlen2(3)
               write(uout,'(a)') 'lj wall at z = -wall_z_ext and hard wall at +wall_z_ext'
               write(uout,'(a)') 'with the additional potential u_drift = u0_drift + u1_drift * z'
               write(uout,'(a,t35,3g10.3)') 'wall_z_ext                     = ', wall_z_ext
               write(uout,'(a,t35,5g10.3)') 'sigma_ext       (atom type)    = ', sigma_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'epsilon_ext     (atom type)    = ', epsilon_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z3 coefficients (atom type)    = ', z3coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'z9 coefficients (atom type)    = ', z9coeff_ext(1:nat)
               write(uout,'(a,t35,5g10.3)') 'u0_drift                       = ', u0_drift
               write(uout,'(a,t35,5g10.3)') 'u1_drift                       = ', u1_drift
            else if (txuext(ipt) == 'estat_field') then
               write(uout,'(a)') 'electrical field in z-direction'
               write(uout,'(a,t35,3g10.3)') 'efield_ext                     = ', efield_ext(1:3)
            else if (txuext(ipt) == 'hom_charged_walls') then
               write(uout,'(a)') 'homgeneously charged hard walls at z = +-boxlen(3)/2'
               write(uout,'(a,t35,3g10.3)') 'surface charge density (C/m**2) = ', surfchargeden
               write(uout,'(a,t35,3g10.3)') 'long-range correction           = ', llongrangecontr
            else if (txuext(ipt) == 'i_soft_sphere') then
               write(uout,'(a)') 'external potential u(r) = 0                            , r < ruext(1,ipt)'
               write(uout,'(a)') 'external potential u(r) = auext*(r-ruext(1,ipt))**nuext, r > ruext(1,ipt)'
               write(uout,'(a,t35,3g10.3)') 'ruext(1,ipt)                   = ', ruext(1,ipt)
               write(uout,'(a,t35,3g10.3)') 'auext                          = ', auext
               write(uout,'(a,t35,3g10.3)') 'nuext                          = ', nuext
            else if (txuext(ipt) == 'Gunnar_soft_sphere') then
               write(uout,'(a)') 'external potential u(r) = auext*(r/ruext(1,ipt))**nuext'
               write(uout,'(a,t35,3g10.3)') 'ruext(1,ipt)                   = ', ruext(1,ipt)
               write(uout,'(a,t35,3g10.3)') 'auext                          = ', auext
               write(uout,'(a,t35,3g10.3)') 'nuext                          = ', nuext
            else if (txuext(ipt) == 'out_hard_ellipsoid') then
               write(uout,'(a)') 'external potential u(r)         = 0              , outside'
               write(uout,'(a)') 'external potential u(r)         = infinify       , inside'
               write(uout,'(a,t35,3g10.3)') 'length of principal axes       = ', ruext(1:3,ipt)
            else if (txuext(ipt) == 'capsid_shell') then
               write(uout,'(a)') 'hard wall between rcap and rcap+dcap'
               write(uout,'(a,t35,3g10.3)') 'capsid radius                  = ', rcap
               write(uout,'(a,t35,3g10.3)') 'capsid radius + thickness      = ', rcap+dcap
            else if (txuext(ipt) == 'uniform_shell') then
               write(uout,'(a)') 'hard wall between rcap and rcap+dcap and uniform charge'
               write(uout,'(a,t35,3g10.3)') 'capsid radius                  = ', rcap
               write(uout,'(a,t35,3g10.3)') 'capsid radius + thickness      = ', rcap+dcap
            else if (txuext(ipt) == 'sphdielboundary_q') then
               write(uout,'(a,t35,a)')      'spherical dielectric boundary'
               write(uout,'(a)')            'multipole expansion'
               write(uout,'(a,t35,3g10.3)') 'boundaryrad                    = ', boundaryrad
               write(uout,'(a,t35,3g10.3)') 'epsi1                          = ', epsi1
               write(uout,'(a,t35,3g10.3)') 'epsi2                          = ', epsi2
               write(uout,'(a,t35,3i5)')    'lmaxdiel                       = ', lmaxdiel
            else if (txuext(ipt)(1:17) == 'sphdielboundary_p') then
               write(uout,'(a,t35,a)')      'spherical dielectric boundary'
               if (txuext(ipt)(18:18) == '1')  write(uout,'(a)') 'pairwise summation: shorter eta expansion'
               if (txuext(ipt)(18:18) == '2')  write(uout,'(a)') 'pairwise summation: longer eta expansion'
               write(uout,'(a,t35,3g10.3)') 'boundaryrad                    = ', boundaryrad
               write(uout,'(a,t35,3g10.3)') 'epsi1                          = ', epsi1
               write(uout,'(a,t35,3g10.3)') 'epsi2                          = ', epsi2
               write(uout,'(a,t35,3i5)')    'lmaxdiel                       = ', lmaxdiel
            else if (txuext(ipt) == 'sphdielboundary') then
               write(uout,'(a,t35,a)')      'spherical dielectric boundary'
               write(uout,'(a,t35,3g10.3)') 'boundaryrad                    = ', boundaryrad
               write(uout,'(a,t35,3g10.3)') 'epsi1 (inside)                 = ', epsi1
               write(uout,'(a,t35,3g10.3)') 'epsi2 (outside)                = ', epsi2
               write(uout,'(a,t35,3i5)')    'lmaxdiel                       = ', lmaxdiel
               if(lbg) write(uout,'(a)')            'neutralizing background is invoked, be sure of what you are doing'
            else if (txuext(ipt) == 'insulating_sphere') then
               write(uout,'(a,t35,a)')      'particle type                  = ', txpt(ipt)
               write(uout,'(a)') 'penetrable homogeneously charged sphere at the origin'
               write(uout,'(a,t35,3g10.3)') 'Radius                         = ', rInsSphere
               write(uout,'(a,t35,3g10.3)') 'Charge                         = ', zInsSphere
            else if (txuext(ipt) == 'core_shell') then
               write(uout,'(a,t35,a)')      'particle type                  = ', txpt(ipt)
               write(uout,'(a)') 'hard walls at R = rChargeIn and rChargeIn'
               write(uout,'(a,t35,3g10.3)') 'rChargeIn                      = ', rChargeIn
               write(uout,'(a,t35,3g10.3)') 'rChargeOut                     = ', rChargeOut
            else if (txuext(ipt) == 'hollow_sphere') then
               write(uout,'(a,t35,a)')      'particle type                  = ', txpt(ipt)
               write(uout,'(a)') 'penetrable homogeneously charged Sphere at the origin'
               write(uout,'(a,t35,3g10.3)') 'Radius der inneren Schale      = ', rInSphere
               write(uout,'(a,t35,3g10.3)') 'Radius der ausseren Schale     = ', rOutSphere
               write(uout,'(a,t35,3g10.3)') 'Charge                         = ', zInsSphere
            else if (txuext(ipt) == 'hard_cylinder') then
               write(uout,'(a)') 'hard cylinder wih a line charge along the z-direction '
               write(uout,'(a)') 'external potential: psi = EpsiFourPi*zcylinder*log((sqrt(r^2+L/2)+L/2)/(sqrt(r^2+L/2)-L/2))'
               write(uout,'(a,t35,3g10.3)') 'radius of the hard cylinder    = ', rCylinder
               write(uout,'(a,t35,3g10.3)') 'line charge density (e/length) = ', zCylinder
            else if (txuext(ipt) == 'lekkerkerker-tuinier') then
               write(uout,'(a)') 'external potential: lekkerkerker-tuinier depletion interaction'
               write(uout,'(a,t35,3g10.3)') 'radius of penetrable hs        = ', rad_dep
               write(uout,'(a,t35,3g10.3)') 'number density of penetrable hs= ', rho_dep
               write(uout,'(a,t35,3g10.3)') 'depletion-thickness factor     = ', factor_dep
            end if
         end do
      end if

   end select

end subroutine IOPotExternal

!************************************************************************
!> \page potential potential.F90
!! **IOPolarizationIter**
!! *perform i/o on polarization iteration variables*
!************************************************************************


!> \page nmlPolarizationIter
!! The namelist  \ref nmlPolarizationIter contains variables that control the calculation of the many-body polarization contribution to the potential energy and forces, ref. JPC 94, 1649 (1990).
!! * Variables:
!!  * \subpage tpolit
!!  * \subpage mpolit
!!  * \subpage npolit
!!  * \subpage ldamping

subroutine IOPolarizationIter(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOPolarizationIter'
   character(80), parameter :: txheading ='polarization calculation data'

   namelist /nmlPolarizationIter/ tpolit, mpolit, npolit, ldamping

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      tpolit  = 0.0001
      mpolit  = 15
      npolit  = 5
      ldamping = .false.

      rewind(uin)
      read(uin,nmlPolarizationIter)

! .. setup for analysis, force iteration of induced dipole moments

      if (lana) npolit = 1

      if (tpolit <= Zero) call Stop(txroutine, 'tpolit <= 0.0', uout)
      if (mpolit < 1    ) call Stop(txroutine, 'mpolit < 1', uout)
      if (npolit < 1    ) call Stop(txroutine, 'npolit < 1', uout)

   case (iWriteInput)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f10.6)') 'abs. tol. of polarization iter.= ', tpolit
         write(uout,'(a,t35,i10)')   'maximum number of iterations   = ', mpolit
         write(uout,'(a,t35,i10)')   'interval of iteration          = ', npolit
         if (ldamping) write(uout,'(a)') 'electrostatics are dampted'
      end if

   end select

end subroutine IOPolarizationIter

