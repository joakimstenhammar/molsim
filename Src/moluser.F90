!************************************************************************
!> \page moluser moluser.F90
!! **PotentialUser**
!! *driver of user-provided routines for potentials*
!************************************************************************


subroutine PotentialUser(ipt, jpt, ibuf, lsetatat2, lsetconf)

   use PotentialModule
   implicit none

   integer(4), intent(in)    :: ipt, jpt
   integer(4), intent(inout) :: ibuf
   logical,    intent(inout) :: lsetatat2(*)
   logical,    intent(out)   :: lsetconf   ! =.true. if match and generation of tables
                                           ! =.false. otherwise
   integer(4) :: iptjpt

   external SodiumChloride
!   external SquareWell
   external EwaldSquareWell
   external ElstatScreen

! ... see if user provided routine exist for txpot(iptjpt)

   iptjpt = iptpt(ipt,jpt)
   lsetconf =.true.
   if (txpot(iptjpt) == 'nana' .or. txpot(iptjpt) == 'nacl' .or. txpot(iptjpt) == 'clcl') then
      call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat2, SodiumChloride)
   else if (txpot(iptjpt) == 'ewald_square_well') then
      call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat2, EwaldSquareWell)
   else if (txpot(iptjpt) == 'elstat_screen') then
      call PotTwoBodyTab2(ipt, jpt, ibuf, lsetatat2, ElstatScreen)
   else
      lsetconf =.false.
   end if

end subroutine PotentialUser

!************************************************************************
!> \page moluser moluser.F90
!! **SodiumChloride**
!! *na-na, na-cl, and cl-cl potentials*
!************************************************************************


!     ref: JCP 84, 867 (1986)

subroutine SodiumChloride(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8),      intent(in)    :: r1
   real(8),      intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='SodiumChloride'
   integer(4) :: iptjpt, iatjat
   real(8), save :: a1, a2, a3, a4, a5, a6, qq
   real(8)    :: r1i, term1, term2, term3, term4, term5


   if (str(1:4) == 'init') then

      iptjpt = iptpt(ipt,jpt)
      if (txpot(iptjpt) == 'nana') then
         a2 =                        (-  4.61112)
         a1 = AuTokJ / Bohr**(+a2) * ( 140.38399)
         a3 = AuTokJ / Bohr        * (-  0.25701)
         a4 = One / Bohr           * (   1.18166)
         a6 =                        (   5.57162)
         a5 = AuTokJ / Bohr**(-a6) * (  13.19009)
         qq = AuTokJ / Bohr**(-1)  * (   1.0    )
      else if (txpot(iptjpt) == 'nacl') then
         a2 =                        (-  3.33344)
         a1 = AuTokJ / Bohr**(+a2) * (  76.82255)
         a3 = AuTokJ / Bohr        * (-  0.25408)
         a4 = One / Bohr           * (   0.66391)
         a6 =                        (+  3.99429)
         a5 = AuTokJ / Bohr**(-a6) * (  16.75816)
         qq = AuTokJ / Bohr**(-1)  * (-  1.0    )
      else if (txpot(iptjpt) == 'clcl') then
         a2 =                        (-  1.03475)
         a1 = AuTokJ / Bohr**(+a2) * (   2.14167)
         a3 = AuTokJ / Bohr        * (-  0.06119)
         a4 =   One / Bohr         * (   0.58659)
         a6 =                        (  11.34935)
         a5 = AuTokJ / Bohr**(-a6) * ( 627519.57840)
         qq = AuTokJ / Bohr**(-1)  * (   1.0    )
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

      r1i   = One/r1
      term1 = a1*r1**( a2)
      term2 = a3*r1
      term3 = exp(max(-expmax,-a4*r1))
      term4 = a5*r1**(-a6)
      term5 = qq*r1i
      u0 = (term1+term2)*term3         +  term4 + term5
      u1 = (term1*( a2)*r1i+a3)*term3  + (term1+term2)*(-a4)*term3 &
         +  term4*(-a6)*r1i            +  term5*(-One)*r1i
      u2 = (term1*( a2)*( a2-1.)*r1i**2)*term3  &
         + Two*(term1*( a2)*r1i+a3)*(-a4)*term3 &
         + (term1+term2)*(-a4)*(-a4)*term3      &
         +  term4*(-a6)*(-a6-1.)*r1i**2         &
         +  term5*(-1.)*(-2.)*r1i**2

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine SodiumChloride

!************************************************************************
!> \page moluser moluser.F90
!! **EwaldSquareWell**
!! *real space ewald plus square-well potential (with a soft rise)*
!************************************************************************


subroutine EwaldSquareWell(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='EwaldSquareWell'
   real(8) :: u0x, u1x, u2x

   if (str(1:4) == 'init') then

     call EwaldRealSpacePot(str, ipt, jpt, iat, jat, r1, u0x, u1x, u2x)
     call SquareWell(str, ipt, jpt, iat, jat, r1, u0x, u1x, u2x)

   else if (str(1:4) == 'calc') then

      u0 = Zero
      u1 = Zero
      u2 = Zero
      call EwaldRealSpacePot(str, ipt, jpt, iat, jat, r1, u0x, u1x, u2x)
      u0 = u0 + u0x
      u1 = u1 + u1x
      u2 = u2 + u2x
      call SquareWell(str, ipt, jpt, iat, jat, r1, u0x, u1x, u2x)
      u0 = u0 + u0x
      u1 = u1 + u1x
      u2 = u2 + u2x

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

contains

!........................................................................

!     real ewald term: 1/r -> (1-erf(r*ualpha))/r      with ualpha = ualphared/rcut

subroutine EwaldRealSpacePot(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='EwaldRealSpacePot'
   real(8), parameter :: fac = 1.128379167d0
   integer(4) :: iatjat, m
   real(8) :: c0fac, c1fac, c2fac
   real(8) :: r1i, rni, w0, w1, w2, term

   if (str(1:4) == 'init') then

      do iat = iatpt(ipt), iatpt(ipt)+natpt(ipt)-1
         do jat = iatpt(jpt), iatpt(jpt)+natpt(jpt)-1
            iatjat = iatat(iat,jat)
            r2umin(iatjat) = r2atat(iatjat)-1.0d-4
!           r2umin(iatjat) = r2uminin
            r2umax(iatjat) = (rumaxfac*(rcut+racom(ipt)+racom(jpt)))**2

! ... check existence of 1/r-term

            if (ipot(1,iatjat) /= 1) call Stop(txroutine, 'no 1/r-term', uout)

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

      c0fac = One
      c1fac = One
      c2fac = One
      if (lewald .and. lcharge) then
         call CalcEwaldUfac(r1, ualpha, c0fac, c1fac, c2fac)
      else
         c0fac = One
         c1fac = One
         c2fac = One
      end if
      u0 = u0 + ucoffx(1,iatjat)*c0fac*r1i
      u1 = u1 - ucoffx(1,iatjat)*c1fac*r1i**2
      u2 = u2 + ucoffx(1,iatjat)*c2fac*Two*r1i**3

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine EwaldRealSpacePot

!........................................................................

!     square-well term: u(r) =  epsilon_well*(0.5 + (1/pi)*atan(alapha_well*(r-r_well))

subroutine SquareWell(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='SquareWell'
   real(8), save :: epsilon_well
   real(8), save :: r_well
   real(8), save :: alpha_well
   integer(4) :: iatjat
   real(8)    :: x

   if (str(1:4) == 'init') then

      epsilon_well = -One     ! depth of the square well (energy)
      r_well       = 3.0d0      ! end of the square well   (length)
      alpha_well   = 1.0d-3     ! steepness of the rise (lower alpha => steeper rise)

       write(uout,*) 'square-well potential'
       write(uout,*) '---------------------'
       write(uout,'(a,t35,3g15.5)') 'epsilon_well                   = ', epsilon_well
       write(uout,'(a,t35,3g15.5)') 'r_well                         = ', r_well
       write(uout,'(a,t35,3g15.5)') 'alpha_well                     = ', alpha_well
       write(uout,*)

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

      x = r1 - r_well
      u0 = epsilon_well*(Half - (One/Pi)*atan(x/alpha_well))
      u1 = -epsilon_well/Pi * alpha_well / (alpha_well**2 + x**2)
      u2 = -Two*x* u1 / (alpha_well**2 + x**2)

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine SquareWell

!........................................................................

end subroutine EwaldSquareWell

!************************************************************************
!> \page moluser moluser.F90
!! **ElstatScreen**
!! *potential where the first term is screened if it is a 1/r term*
!************************************************************************

!     other terms may also be present
!     notice that epsilonr has to be set in this routine
!     ucoff for the 1/r term is not used instead zat is used

subroutine ElstatScreen(str, ipt, jpt, iat, jat, r1, u0, u1, u2)

   use PotentialModule
   implicit none

   character(*), intent(in)    :: str
   integer(4),   intent(in)    :: ipt, jpt
   integer(4),   intent(inout) :: iat, jat
   real(8)   ,   intent(in)    :: r1
   real(8)   ,   intent(out)   :: u0, u1, u2

   character(40), parameter :: txroutine ='ElstatScreen'
   real(8), save :: facscr
   integer(4) :: iatjat, m
   real(8)    :: r1i, rni, w0, w1, w2, term, qq, expfac

   if (str(1:4) == 'init') then

      if (scrlen <= Zero) call Stop(txroutine, 'scrlen <= 0', uout)
      facscr = One/scrlen

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
      if (ldipole) qq = Zero
      if (lpolarization) qq = Zero
      r1i = One/r1
      expfac = exp(-facscr*r1)
      rni = One
      w0 = Zero
      w1 = Zero
      w2 = Zero
      do m = 1, ipotm
         rni  = rni * r1i
         if (ipot(m,iatjat) == 1) cycle
         term = ucoffx(m,iatjat) * rni
         w0   = w0 + term
         w1   = w1 + (-m) * term
         w2   = w2 + (-m) * (-m-1) * term
      end do

      u0 = qq*r1i*expfac                                        +w0
      u1 =-qq*r1i*(facscr+r1i)*expfac                           +w1*r1i
      u2 = qq*r1i*(facscr**2+Two*facscr*r1i+Two*r1i**2)*expfac  +w2*r1i*r1i

   else

      call Stop(txroutine, 'str value unsupported', uout)

   end if

end subroutine ElstatScreen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **SetConfigurationUser**
!! *driver of user-provided routines for generating a start configuration*
!************************************************************************


subroutine SetConfigurationUser(ipt, lsetconf)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type
   logical,    intent(out) :: lsetconf     ! =.t, if match and call of user supplied lattice routine
                                           ! =.false. otherwise
   external setpa3

   lsetconf =.true.

   if (txsetconf(ipt) == 'setcnf') then
      call SetCNF(ipt)
   else if (txsetconf(ipt) == 'onehelix') then
      call SetOnehelix(ipt)
   else if (txsetconf(ipt) == 'cubic2d1surf') then
      call SetCubic2D1Surf(ipt)
   else if (txsetconf(ipt) == 'cubic2d2surf') then
      call SetCubic2D2Surf(ipt)
   else if (txsetconf(ipt) == 'random_trans') then
      call SetChainRandomTrans() !ipt is not needed
   else if (txsetconf(ipt) == 'capsid') then
      call SetCapsid(ipt)
   else if (txsetconf(ipt) == 'capsid_250') then
      call SetCapsid_250(ipt)
   else if (txsetconf(ipt) =='random_capsid_inside') then
      call SetInsideCapsid(ipt)
   else if (txsetconf(ipt) =='random_chain_inside') then
      call SetChainInsideCapsid() !ipt is not needed
   else if (txsetconf(ipt) =='randomcapsid') then
      call SetCapsidRandom(ipt)
   else if (txsetconf(ipt) =='spool') then
      call SetSpool(ipt)
   else if (txsetconf(ipt) =='spoolm') then
      call SetSpoolm(ipt)
   else if (txsetconf(ipt) =='spoolm1') then
      call SetSpoolm1(ipt)
   else if (txsetconf(ipt) =='sheet') then
      call setsheet(ipt,1.0d0,1.0d0,.true.)
   else if (txsetconf(ipt) =='randomcylindershell') then
      call SetRandomCylinderShell(ipt)
   else if (txsetconf(ipt) =='chainlinepma') then
      call SetChainLinePMA() !ipt is not needed
   else
      lsetconf =.false.
   end if

end subroutine SetConfigurationUser

!************************************************************************
!> \page moluser moluser.F90
!! **SetCNF**
!! *generate a configuration by reading ro and ori from the begining of fcnf*
!************************************************************************


subroutine SetCNF(ipt)

   use CoordinateModule
   implicit none

   logical      :: ldum
   integer(4)   :: ip, ipt, iplow, ipupp
   integer(4)   :: idum
   real(8)      :: dum

   if (master) then
      iplow = ipnpt(ipt)
      ipupp = ipnpt(ipt)+nppt(ipt)-1
      rewind ucnf
      read(ucnf) idum, idum, dum, dum, dum
      read(ucnf) idum, idum, dum, dum, ldum, ldum
      read(ucnf) (ro(1:3,ip),qua(0:3,ip),ip = iplow,ipupp)
   end if

#if defined (_PAR_)
   call par_bc_reals(ro(1,iplow)       , 3*(nppt(ipt)))
   call par_bc_reals(qua(1,iplow)      , 4*(nppt(ipt)))
#endif

   call QuaToOri(np, iplow, ipupp, qua, ori)

   call SetAtomPos(iplow, ipupp, .false.)

end subroutine SetCNF

!************************************************************************
!> \page moluser moluser.F90
!! **SetOneHelix**
!! *generate a helical conformation*
!************************************************************************


subroutine SetOneHelix(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt

   character(40), parameter :: txroutine ='SetOneHelix'
   integer(4) :: ict, ic, iseg, ip
   real(8)    :: pitch, radhelix, phi, nbeadperturn

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   ict = ictpt(ipt)          ! chain type
   if (ict == 0) return      ! see if chain
   ic = icnct(ict)           ! first chain of chain type ict
   if (ncct(ict) /= 1) call Stop(txroutine, 'ncct(ict) /= 1', uout)

   pitch = 2.5*radatset(ipt) ! set value of pitch
   radhelix = 20.0           ! set value of radhelix

   phi = Two*asin((bond(ict)%eq/Two)/radhelix)
   nbeadperturn = TwoPi/phi

!  write(*,*) 'ipt', ipt
!  write(*,'(a,g15.5)') 'picth       ', pitch
!  write(*,'(a,g15.5)') 'radhelix    ', radhelix
!  write(*,'(a,g15.5)') 'phi       ', phi
!  write(*,'(a,g15.5)') 'nbeadperturn', nbeadperturn

   do iseg = 1, npct(ict)
      ip = ipnsegcn(iseg,ic)
      ro(1,ip) = radhelix*(cos((iseg-1)*phi))
      ro(2,ip) = radhelix*(sin((iseg-1)*phi))
      ro(3,ip) = pitch*(iseg-1)/nbeadperturn
!     write(*,'(4g15.5)') ip, ro(1:3,ip)
   end do

end subroutine SetOneHelix

!************************************************************************
!> \page moluser moluser.F90
!! **SetCubic2D1Surf**
!! *generate a cubic configuration with particles at z = -(boxlenz/2 + delta)*
!************************************************************************


subroutine SetCubic2D1Surf(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetCubic2D1Surf'
   integer(4) :: nset, nsurf, m, ix, iy, isurf, ip
   real(8) :: xoffset, yoffset, dx, dy

   nsurf = 1                               ! number of surfaces (1 or 2)
   m = sqrt(real(nppt(ipt)/nsurf))
   if (nsurf*m**2 /= nppt(ipt)) call Stop(txroutine, 'nsurf*m**2 /= nppt(ipt)', uout)

   dx = boxlen(1)/m
   dy = boxlen(2)/m

   xoffset = -boxlen2(1) + Half*dx
   yoffset = -boxlen2(2) + Half*dx

   nset = 0
   do ix = 0, m-1
      do iy = 0, m-1
         do isurf = 1, nsurf
            nset = nset + 1
            ip = nset - 1 + ipnpt(ipt)
            ro(1,ip) = xoffset + ix*dx
            ro(2,ip) = yoffset + iy*dy
            ro(3,ip) = (2*isurf-3)*(boxlen2(3) + Two)
            call SetAtomPos(ip, ip, .false.)
         end do
      end do
   end do

end subroutine SetCubic2D1Surf

!************************************************************************
!> \page moluser moluser.F90
!! **SetCubic2D2Surf**
!! *generate a cubic configuration with particles at z = +-(boxlenz/2 + delta)*
!************************************************************************


subroutine SetCubic2D2Surf(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetCubic2D2Surf'
   integer(4) :: nset, nsurf, m, ix, iy, isurf, ip
   real(8) :: xoffset, yoffset, dx, dy

   nsurf = 2                               ! number of surfaces (1 or 2)
   m = sqrt(real(nppt(ipt)/nsurf))
   if (nsurf*m**2 /= nppt(ipt)) call Stop(txroutine, 'nsurf*m**2 /= nppt(ipt)', uout)

   dx = boxlen(1)/m
   dy = boxlen(2)/m

   xoffset = -boxlen2(1) + Half*dx
   yoffset = -boxlen2(2) + Half*dx

   nset = 0
   do ix = 0, m-1
      do iy = 0, m-1
         do isurf = 1, nsurf
            nset = nset + 1
            ip = nset - 1 + ipnpt(ipt)
            ro(1,ip) = xoffset + ix*dx
            ro(2,ip) = yoffset + iy*dy
            ro(3,ip) = (2*isurf-3)*(boxlen2(3) + Two)
            call SetAtomPos(ip, ip, .false.)
         end do
      end do
   end do

end subroutine SetCubic2D2Surf

!************************************************************************
!> \page moluser moluser.F90
!! **SetChainRandomTrans**
!! *generate a random start configurations for chain particles + chain translation*
!************************************************************************

!     for Niklas

subroutine SetChainRandomTrans!(iptset) iptset is not needed

   use CoordinateModule
   implicit none

   !integer(4), intent(in) :: iptset

   character(40), parameter :: txroutine ='SetChainRandomTrans'
   integer(4) :: ntry, nset, itry, ic, ict, ip, ipt, jp
   real(8)    :: zaim, zlow
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   if (.not.first) return                                ! should be called only once
   first =.false.

   do ic = 1, nc                                        ! loop over chains
      ict = ictcn(ic)

      ntry = 100*npct(ict)
      nset = 1
      do itry = 1, ntry                                 ! loop over attempts to set a particle
         ip = ipnsegcn(nset,ic)
         ipt = iptpn(ip)

         if (nset == 1) then                             ! a first segment (at the origin)
            ro(1:3,ip) = Zero
         else                                           ! remaining segments
            jp = ipnsegcn(nset-1,ic)
            call SetPartPosRandomN(ip, jp, bond(ict)%eq)
         end if
         if (CheckPartOutsideBox(ip)) cycle
         if (CheckTooFoldedChain(ip,bond(ict)%eq)) cycle! check if a too folded chain
         call SetPartOriRandom(iseed,ori(1,1,ip))       ! set trial random particle orientation
         call SetAtomPos(ip,ip,.false.)                 ! set trial atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
         lpset(ip) =.true.                              ! trial configuration accepted
         if (nset == npct(ict)) exit                     ! check if finnished
         nset = nset+1                                  ! update nset

      end do

      if (itry > ntry) then                         ! number of trial attempts exceeds the maximal one ?
         if (master) write(uout,'(i5,a,i5,a,i5,a,i5,a)') nset, 'of', npct(ict), 'particles in chain', ic, 'set after', ntry, 'attemps'
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
     end if

! .. translate the chain so that the particle with smallest z-coordinate is positioned at zaim

     zaim = -boxlen2(3) + 12.0                       ! Niklas 2005-05-05
     zlow = minval(ro(3,1:np))                    ! Niklas 2005-05-05
     ro(3,1:np) = (zaim-zlow) + ro(3,1:np)        ! Niklas 2005-05-05
! write(*,*) sum(ro(3,1:np))/np                   ! Niklas 2004-04-01
   end do

end subroutine SetChainRandomTrans

!************************************************************************
!> \page moluser moluser.F90
!! **SetCapsid**
!! *generate a uniform distribution on the surface of a sphere of rcap+2*
!************************************************************************

!     radius for particles of type ipt

subroutine SetCapsid(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetCapsid'
   integer(4) :: ntry, nset, itry, ip
   real(8)    :: x ,y, z
   logical    :: lWarnHCOverlap

   if (nppt(ipt)==0) return

   ntry =  10*nppt(ipt)
   nset = 1
   do itry=1,ntry

      ip=nset-1+ipnpt(ipt)

      call SphRandom(iseed,x,y,z)
      ro(1,ip)=(rcap+2)*x
      ro(2,ip)=(rcap+2)*y
      ro(3,ip)=(rcap+2)*z

      call SetPartOriRandom(iseed,ori(1,1,ip))  ! set trial random particle orientation
      call SetAtomPos(ip,ip,.false.)            ! set trial atom positions
      if (lWarnHCOverlap(ip, radatset, .true.)) cycle  ! check if atom-atom hard-core overlap
      lpset(ip)=.true.                          ! trial configuration accepted
      if (nset==nppt(ipt)) exit                 ! check if finnished
      nset=nset+1                               ! update nset

   end do

   if (itry>10*nppt(ipt)) then            ! number of trial attempts exceeds the maximal one ?
      if (master) write(uout,'(i5,a,i5,a,i5,a,i5,a)') nset, 'of', nppt(ipt), 'particles of type', ipt,' set after', ntry, 'attemps'
      call Stop(txroutine, 'random configuration for capsid failed, itry > ntry', uout)
   end if
 end subroutine SetCapsid

!************************************************************************
!> \page moluser moluser.F90
!! **SetCapsid_250**
!! *generate a random configuration with 250 particles on a spherical surface*
!************************************************************************


subroutine SetCapsid_250(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   integer(4) :: nset, ip
   integer(4) :: itry                      ! number of circle
   integer(4) :: jtry
   integer(4) :: n(19)                     ! number of particles on each circle

   if (nppt(ipt)/=250) return

   nset=1

   do itry=1,19
      n(itry)=int(20*sin(pi*itry/20))    ! provided that int(x) returns [x]
      do jtry=1,n(itry)
         ip=nset-1+ipnpt(ipt)
         ro(1,ip)=(rcap+2)*cos(twopi*jtry/n(itry))*sin(pi*itry/20)
         ro(2,ip)=(rcap+2)*sin(twopi*jtry/n(itry))*sin(pi*itry/20)
         ro(3,ip)=(rcap+2)*cos(pi*itry/20)
         call SetPartOriRandom(iseed,ori(1,1,ip))      ! set trial random particle orientation
         call SetAtomPos(ip,ip,.false.)                ! set trial atom positions
         lpset(ip)=.true.
         nset=nset+1                                   ! update nset
      end  do
   end  do

   ro(1,ip+1)=0
   ro(2,ip+1)=0
   ro(3,ip+1)=rcap+2
   call SetPartOriRandom(iseed,ori(1,1,ip+1))         ! set trial random particle orientation
   call SetAtomPos(ip+1,ip+1,.false.)                 ! set trial atom positions
   lpset(ip+1)=.true.

   ro(1,ip+2)=0
   ro(2,ip+2)=0
   ro(3,ip+2)=-(rcap+2)
   call SetPartOriRandom(iseed,ori(1,1,ip+2))         ! set trial random particle orientation
   call SetAtomPos(ip+2,ip+2,.false.)                 ! set trial atom positions
   lpset(nppt(ipt))=.true.

  end subroutine  SetCapsid_250

!************************************************************************
!> \page moluser moluser.F90
!! **SetCapsidRandom**
!! *generate a random start configurations for paticles of type ipt*
!************************************************************************

!     in a spherical cell containing the capsid

subroutine SetCapsidRandom(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetCapsidRandom'
   integer(4) :: ntry, nset, itry, ip
   real(8)    :: Random
   logical    :: lWarnHCOverlap

   if (nppt(ipt)==0) return

   ntry = 10*nppt(ipt)
   nset=1
   do itry=1,ntry

      ip=nset-1+ipnpt(ipt)

! ... trial particle position

      do
        ro(1,ip)=Two*sphrad*(Random(iseed)-0.5)
        ro(2,ip)=Two*sphrad*(Random(iseed)-0.5)
        ro(3,ip)=Two*sphrad*(Random(iseed)-0.5)
        if ((ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2>=(rcap-radat(ipt))**2).and.   &
        (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2<=(rcap+dcap+radat(ipt))**2))  cycle
        if (ro(1,ip)<-Two*sphrad.or.ro(1,ip)>Two*sphrad) cycle
        if (ro(2,ip)<-Two*sphrad.or.ro(2,ip)>Two*sphrad) cycle
        if (ro(3,ip)<-Two*sphrad.or.ro(3,ip)>Two*sphrad) cycle
        if (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2<=(rcap-radat(ipt))**2-0.1d0)  then
           exit
        else if ((ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2>=(rcap+dcap+radat(ipt))**2-0.1d0).and. &
            (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2<=sphrad**2-0.1d0))  then
           exit
        endif
      end do

      call SetPartOriRandom(iseed,ori(1,1,ip))  ! set trial random particle orientation
      call SetAtomPos(ip,ip,.false.)            ! set trial atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
      lpset(ip)=.true.                          ! trial configuration accepted
      if (nset==nppt(ipt)) exit                  ! check if finnished
      nset=nset+1                               ! update nset

   end do

if (itry > ntry) then                            ! number of trial attempts exceeds the maximal one ?
      if (master) write(uout,'(i5,a,i5,a,i5,a,i5,a)') &
         nset, 'of', nppt(ipt), 'particles of type', ipt,' set after', ntry, 'attemps'
      call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
   end if

end subroutine SetCapsidRandom

!************************************************************************
!> \page moluser moluser.F90
!! **SetInsideCapsid**
!! *generate a random configuration inside a capsid*
!************************************************************************


subroutine SetInsideCapsid(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetInsideCapsid'
   integer(4) :: ntry, nset, itry, ip
   real(8)    :: Random
   logical    :: lWarnHCOverlap

   if (nppt(ipt)==0) return

   ntry = 10*nppt(ipt)
   nset = 1

   do itry=1,ntry

      ip=nset-1+ipnpt(ipt)

! ... particle position

      do
         ro(1,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
         ro(2,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
         ro(3,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
         if (ro(1,ip)<-Two*(rcap-radat(ipt)).or.ro(1,ip)>Two*(rcap-radat(ipt))) cycle
         if (ro(2,ip)<-Two*(rcap-radat(ipt)).or.ro(2,ip)>Two*(rcap-radat(ipt))) cycle
         if (ro(3,ip)<-Two*(rcap-radat(ipt)).or.ro(3,ip)>Two*(rcap-radat(ipt))) cycle
         if (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2<=(rcap-radat(ipt))**2-0.1d0) exit
      end do

      call SetPartOriRandom(iseed,ori(1,1,ip))  ! set trial random particle orientation
      call SetAtomPos(ip,ip,.false.)            ! set trial atom positions
      if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
      lpset(ip)=.true.                          ! trial configuration accepted
      if (nset==nppt(ipt)) exit                  ! check if finnished
      nset=nset+1                               ! update nset

   end do

if (itry > ntry) then                     ! number of trial attempts exceeds the maximal one ?
      if (master) write(uout,'(i5,a,i5,a,i5,a,i5,a)') nset, 'of', nppt(ipt), 'particles of type', ipt,' set after', ntry, 'attemps'
      call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
   end if

end subroutine SetInsideCapsid

!************************************************************************
!> \page moluser moluser.F90
!! **SetChainInsideCapsid**
!! *generate a random chain configuration inside a capsid*
!************************************************************************


subroutine SetChainInsideCapsid!(iptset) iptset is not needed

   use CoordinateModule
   implicit none

   !integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainInsideCapsid'
   integer(4) :: ntry, nset, itry, ic, ict, ip, ipt, jp
   real(8)    :: Random
   logical    :: first =.true.
   logical    :: CheckTooFoldedChain, lWarnHCOverlap

   if (.not.first) return                           ! should be called only once
   first =.false.

   do ic = 1, nc                                   ! loop over chains
      ict = ictcn(ic)                              ! chain type

      ntry = 1000*npct(ict)
      nset = 1

      do itry = 1, ntry                            ! loop over attempts to set a particle
         ip = ipnsegcn(nset,ic)
         ipt = iptpn(ip)

         if (nset == 1) then                        ! a first segment
            do
               ro(1,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
               ro(2,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
               ro(3,ip)=Two*(rcap-radat(ipt))*(Random(iseed)-0.5)
               if (ro(1,ip)<-Two*(rcap-radat(ipt)).or.ro(1,ip)>Two*(rcap-radat(ipt))) cycle
               if (ro(2,ip)<-Two*(rcap-radat(ipt)).or.ro(2,ip)>Two*(rcap-radat(ipt))) cycle
               if (ro(3,ip)<-Two*(rcap-radat(ipt)).or.ro(3,ip)>Two*(rcap-radat(ipt))) cycle
               if (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2<=(rcap-radat(ipt))**2-0.1d0) exit
            end do
         else                                            ! remaining segments
            jp = ipnsegcn(nset-1,ic)
            call SetPartPosRandomN(ip, jp, bond(ict)%eq)
            if (ro(1,ip)<-Two*(rcap-radat(ipt)).or.ro(1,ip)>Two*(rcap-radat(ipt))) cycle
            if (ro(2,ip)<-Two*(rcap-radat(ipt)).or.ro(2,ip)>Two*(rcap-radat(ipt))) cycle
            if (ro(3,ip)<-Two*(rcap-radat(ipt)).or.ro(3,ip)>Two*(rcap-radat(ipt))) cycle
            if (ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2>=(rcap-radat(ipt))**2-0.1d0) cycle
         end if
         if (CheckTooFoldedChain(ip,bond(ict)%eq)) cycle ! check if a too folded chain
         call SetPartOriRandom(iseed,ori(1,1,ip))        ! set trial random particle orientation
         call SetAtomPos(ip,ip,.false.)                  ! set trial atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
         lpset(ip) =.true.                               ! trial configuration accepted
         if (nset == npct(ict)) exit                     ! check if finnished
         nset = nset+1                                   ! update nset

      end do

      if (itry > ntry) then                         ! number of trial attempts exceeds the maximal one ?
         if (master) write(uout,'(i5,a,i5,a,i5,a,i5,a)') &
            nset, 'of', npct(ict), 'particles in chain', ic, 'set after', ntry, 'attemps'
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if

  end do

end subroutine SetChainInsideCapsid

!************************************************************************
!> \page moluser moluser.F90
!! **SetSpool**
!! *generate a spool-like structure inside a capsid*
!************************************************************************


!     The distance between the projections of two adjacent circles on Oz axisis constant and is 6.4

subroutine SetSpool(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetSpool'

   integer(4) :: nset, ip
   integer(4) :: ncircle, n, nbeads        ! number of circles in one semi-sphere
                                           ! number of beads on circle icircle
   integer(4) :: icircle, jtry, ntry
   real(8)    :: bb, delta, alpha          ! rms bead-to-bead distance, circle-to-circle distance
                                           ! angle at centre
   logical    :: lWarnHCOverlap

   alpha = twopi*8/360
   bb = 6.0
   delta = (rcap-2*radat(ipt))*sin(alpha)
   ncircle = 0
   nbeads  = 0
   if (nppt(ipt) < int(twopi*(rcap-2*radat(ipt))/bb)) then
       ncircle = 0
   else
      do  n = 1,10
         nbeads =  nbeads + int(twopi*(rcap-2*radat(ipt))*cos(n*alpha)/bb)
         if (nbeads >  (nppt(ipt)-int(twopi*(rcap-2*radat(ipt))/bb))/2)  then
           ncircle = n
           exit
        end if
     end do
   end if
   if (nbeads < (nppt(ipt)-int(twopi*(rcap-radat(ipt))/bb))/2) then
     call Stop(txroutine, 'spool configuration failed, ncircle > nmax', uout)
   end if

   nset = 1

   do icircle  = -ncircle, ncircle
      ntry = int(twopi*(rcap-2*radat(ipt))*cos(icircle*alpha)/bb)
      do jtry = 1, ntry
         ip=nset-1+ipnpt(ipt)
         ro(1,ip)= (rcap-2*radat(ipt))*cos(icircle*alpha)*cos(twopi*jtry/ntry)
         ro(2,ip)= (rcap-2*radat(ipt))*cos(icircle*alpha)*sin(twopi*jtry/ntry)
         ro(3,ip)= (rcap-2*radat(ipt))*sin(icircle*alpha)
         call SetPartOriRandom(iseed,ori(1,1,ip))    ! set trial random particle orientation
         call SetAtomPos(ip,ip,.false.)              ! set trial atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
         lpset(ip)=.true.
         if (nset==nppt(ipt)) exit
         nset=nset+1                   ! update nset
      end do
      if (nset==nppt(ipt)) exit
   end do

end subroutine SetSpool

!************************************************************************
!> \page moluser moluser.F90
!! **SetSpoolm**
!! *generate a multilayered spool-like structure inside a capsid*
!************************************************************************


subroutine SetSpoolm(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetSpoolm'
   integer(4) :: nset, ip
   integer(4) :: ncircle, n, nbeads, m     ! number of circles in one semi-sphere
                                           ! number of beads on circle icircle
                                           ! number of layer
   integer(4) :: icircle, jtry, ntry, k, n0
   real(8)    :: bb, delta, alpha          ! rms bead-to-bead distance, circle-to-circle distance
                                        ! angle at centre
   logical    :: lWarnHCOverlap

   alpha = twopi*8/360
   bb = 6.0
   delta = (rcap-2*radat(ipt))*sin(alpha)
   m = 2
   ncircle = 0
   n0  = 0
   nbeads = 0

   do k = 1, m
      n0 = n0 + int(twopi*(rcap-2*k*radat(ipt))/bb)
   end do

   if (nppt(ipt) < n0)  then
       ncircle = 0
   else
      do  n = 1,7
         do k = 1, m
         nbeads =  nbeads + int(twopi*(rcap-2*k*radat(ipt))*cos(n*alpha)/bb)
         write(*,*) n, k, nbeads
         if (nbeads >  (nppt(ipt)- n0)/2)  then
            ncircle = n
            exit
         end if
        end do
         if (nbeads >  (nppt(ipt)- n0)/2) then
            exit
         end if
      end do
   end if

   if (nbeads < (nppt(ipt)-n0)/2) then
     call Stop(txroutine, 'spool configuration failed, ncircle > nmax', uout)
   end if

   nset = 1

   do icircle  = -ncircle, ncircle
      do  k = 1, m
         ntry = int(twopi*(rcap-3*k*radat(ipt))*cos(icircle*alpha)/bb)
         do jtry = 1, ntry
             ip=nset-1+ipnpt(ipt)
             ro(1,ip)= (rcap-2*k*radat(ipt))*cos(icircle*alpha)*cos(twopi*jtry/ntry)
             ro(2,ip)= (rcap-2*k*radat(ipt))*cos(icircle*alpha)*sin(twopi*jtry/ntry)
             ro(3,ip)= (rcap-2*k*radat(ipt))*sin(icircle*alpha)
             call SetPartOriRandom(iseed,ori(1,1,ip))    ! set trial random particle orientation
             call SetAtomPos(ip,ip,.false.)              ! set trial atom positions
             if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
             lpset(ip)=.true.
             if (nset==nppt(ipt)) exit
             nset=nset+1                   ! update nset
         end do
         if (nset ==nppt(ipt)) exit
      end do
    if (nset ==nppt(ipt)) exit
   end do

  end subroutine SetSpoolm

!************************************************************************
!> \page moluser moluser.F90
!! **SetSpoolm1**
!! *generate a multilayered spool-like structure inside a capsid*
!************************************************************************


subroutine SetSpoolm1(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt           ! particle type

   character(40), parameter :: txroutine ='SetSpoolm1'
   integer(4) :: nset, ip
   integer(4) :: ncircle, n, nbeads, m     ! number of circles in one semi-sphere
                                           ! number of beads on circle icircle
                                           ! number of layer
   integer(4) :: icircle, jtry, ntry, k, n0
   real(8)    :: bb, delta, alpha               ! rms bead-to-bead distance, circle-to-circle distance
                                                                                        ! angle at centre
   logical    :: lWarnHCOverlap

   alpha = twopi*8/360
   bb = 6.0
   delta = (rcap-2*radat(ipt))*sin(alpha)
   m = 2
   ncircle = 0
   n0  = 0
   nbeads = 0

   do k = 1, m
     n0 = n0 + int(twopi*(rcap-2*k*radat(ipt))/bb)
   end do

   if (nppt(ipt) < n0)  then
       ncircle = 0
   else
      do n = 1,7
         do k = 1, m
            nbeads =  nbeads + int(twopi*(rcap-2*k*radat(ipt))*cos(n*alpha)/bb)
            write(*,*) n, k, nbeads
            if (nbeads >  (nppt(ipt)- n0)/2)  then
               ncircle = n
               exit
           end if
        end do
        if (nbeads >  (nppt(ipt)- n0)/2) then
           exit
        end if
     end do
   end if

   if (nbeads < (nppt(ipt)-n0)/2) then
     call Stop(txroutine, 'spool configuration failed, ncircle > nmax', uout)
   end if

   nset = 1

   do  k = 1, m
      do icircle  = -ncircle, ncircle
         ntry = int(twopi*(rcap-3*k*radat(ipt))*cos(icircle*alpha)/bb)
         do jtry = 1, ntry
            ip=nset-1+ipnpt(ipt)
            ro(1,ip)= (rcap-2*k*radat(ipt))*cos(icircle*alpha)*cos(twopi*jtry/ntry)
            ro(2,ip)= (rcap-2*k*radat(ipt))*cos(icircle*alpha)*sin(twopi*jtry/ntry)
            if (k > 1) then
               ro(3,ip)= -(rcap-2*k*radat(ipt))*sin(icircle*alpha)
            else
               ro(3,ip)= (rcap-2*k*radat(ipt))*sin(icircle*alpha)
            end if
            call SetPartOriRandom(iseed,ori(1,1,ip)) ! set trial random particle orientation
            call SetAtomPos(ip,ip,.false.)           ! set trial atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle       ! check if atom-atom hard-core overlap
            lpset(ip)=.true.
            if (nset==nppt(ipt)) exit
            nset=nset+1                   ! update nset
            end do
         if (nset ==nppt(ipt)) exit
      end do
      if (nset ==nppt(ipt)) exit
   end do

end subroutine SetSpoolm1

!************************************************************************
!> \page moluser moluser.F90
!! **setsheet**
!! *set bcc 'sheet' or 'tube'*
!************************************************************************


subroutine setsheet(ipt,scfac1,scfac2,first)

   use MolModule
   use CoordinateModule
   implicit none

   logical, save :: tube, twoD
   logical      :: lWarnHCOverlap,CheckPartOutsideBox,first
   integer(4), save :: nx, nz
   integer(4)   :: ix, iz, ip, ipt, iplow, ipupp
   real(8)      :: com(3), dalpha, lx, lz, tmprad, scfac1, scfac2

   namelist /nmlSheet/ nx, nz, tube, twoD

   if (first) then
      tube = .true.
      twoD = .false.
      nx = 6
      nz = 12
      rewind(uin)
      read(uin,nmlSheet)
   end if

   iplow = ipnpt(ipt)
   ipupp = ipnpt(ipt)+nppt(ipt)-1


!   tmpang=atan(sqrt(1.0d0-0.6d0*aellipsoid**2/(aellipsoid**2-1.0d0))/(0.6d0*aellipsoid**2/sqrt(aellipsoid**2-1.0d0)))

!   tmpang=atan(sqrt(1.0d0-0.6d0*(aellipsoid**2-1.0d0)/aellipsoid)/sqrt(0.6d0*(aellipsoid**2-1.0d0)))

!   lx = 9.0d0

!   lz = 4.0d0*radellipsoid*sqrt(0.6d0*(aellipsoid**2-1.0d0))

!   lz = 1.02d0*4.0d0*radellipsoid*sqrt(0.6d0*(aellipsoid**2-1.0d0))

   lz = scfac1*4.0d0*radellipsoid*sqrt(0.6d0)*aellipsoid**2/sqrt(aellipsoid**2-1.0d0)

   lx = sqrt(16.0d0*radellipsoid**2-lz**2/aellipsoid**2)

   lz = lz*scfac2

   lx = lx*scfac2

!   write(6,*) lx, lz

   if (mod(nppt(ipt),nx*nz) == 0.and.nppt(ipt)/(nx*nz) == 1) then

      ip=iplow-1
      com(1:3)=0.0d0

      if (tube) then
         dalpha=4.0d0*asin(1.0d0)/nx
         tmprad=lx/(4.0d0*sin(dalpha/4.0d0))
         do iz=1, nz
            do ix=1, nx
               ip=ip+1
               ro(1,ip)= tmprad*sin(dalpha*ix+0.5d0*dalpha*mod(iz,2))
               ro(2,ip)= tmprad*cos(dalpha*ix+0.5d0*dalpha*mod(iz,2))
               ro(3,ip)= 0.5d0*iz*lz
               call SetPartOriLab(ori(1,1,ip))
               com(1)=com(1) + ro(1,ip)
               com(2)=com(2) + ro(2,ip)
               com(3)=com(3) + ro(3,ip)
            end do
         end do

      else if (.not.twoD) then

         do iz=1, nz
            do ix=1, nx
               ip=ip+1
               ro(1,ip)= ix*lx+0.5d0*lx*mod(iz,2)
               ro(2,ip)= 0.0d0
               ro(3,ip)= 0.5d0*iz*lz
               call SetPartOriLab(ori(1,1,ip))
               com(1)=com(1) + ro(1,ip)
               com(2)=com(2) + ro(2,ip)
               com(3)=com(3) + ro(3,ip)
            end do
         end do

      else

         do iz=1, nz
            do ix=1, nx
               ip=ip+1
               ro(1,ip)= ix*lx+0.5d0*lx*mod(iz,2)
               ro(2,ip)= 0.5d0*iz*lz
               ro(3,ip)= 0.0d0

               ori(1,1,ip) = 1.0d0
               ori(2,1,ip) = 0.0d0
               ori(3,1,ip) = 0.0d0

               ori(1,2,ip) = 0.0d0
               ori(2,2,ip) = 0.0d0
               ori(3,2,ip) = 1.0d0

               ori(1,3,ip) = 0.0d0
               ori(2,3,ip) = 1.0d0
               ori(3,3,ip) = 0.0d0

               com(1)=com(1) + ro(1,ip)
               com(2)=com(2) + ro(2,ip)
               com(3)=com(3) + ro(3,ip)
            end do
         end do


      end if
   else

      call Stop('SetSheet', 'Inconsistency btween nx, nz and nppt.', uout)

   end if

!   If (twoD) then
!      Boxlen(1)=nx*lx
!      Boxlen(2)=0.5d0*nz*lz
!   Else
!      Boxlen(3)=0.5d0*nz*lz
!   End If

!   call SetBoxParam
!   if (lewald) call EwaldSetup

   do ip=iplow, ipupp
      ro(1,ip)=ro(1,ip)-com(1)/nppt(ipt)
      ro(2,ip)=ro(2,ip)-com(2)/nppt(ipt)
      ro(3,ip)=ro(3,ip)-com(3)/nppt(ipt)
      call SetAtomPos(ip, ip, .false.)
      if (lWarnHCOverlap(ip, radatset, .true.)) call Stop('SetSheet', 'Hard spere overlap.', uout)
!      if (CheckHCOverlap(ip)) call Stop('SetSheet', 'Hard spere overlap.', uout)
      if (CheckPartOutsideBox(ip)) call Stop('SetSheet', 'Particle outside box.', uout)
   end do

end subroutine setsheet

!************************************************************************
!> \page moluser moluser.F90
!! **SetRandomCylinderShell**
!! *generate a random configuration inside a cylinder with a hard cylinder within*
!************************************************************************


subroutine SetRandomCylinderShell(ipt)

   use CoordinateModule
   implicit none
   integer(4), intent(in) :: ipt                   ! particle type

   character(40), parameter :: txroutine ='SetRandomCylinderShell'
   integer(4) :: ntry, itry, ip, iploc
   real(8) :: rmin2
   logical :: lWarnHCOverlap

   if (nppt(ipt) == 0) return
   ntry = ntrydef
   rmin2 = Radlimit(1)**2

   do iploc = 1, nppt(ipt)                                      ! loop over particles of type ipt
      ip = iploc-1+ipnpt(ipt)                                   ! particle to be set

      do itry = 1, ntry                                         ! loop over attempts to set the particle
         call SetPartPosRandom(ip)                              ! set particle position
         call SetPartOriRandom(iseed,ori(1,1,ip))               ! set particle orientation
         call SetAtomPos(ip,ip,.false.)                         ! set atom positions
         if (ro(1,ip)**2+ro(2,ip)**2 < rmin2) cycle             ! particle appeared in the central core
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle        ! check if atom-atom hard-core overlap
         lpset(ip) = .true.                                     ! configuration accepted
         exit
      end do

      if (itry > ntry) then                         ! number of  attempts exceeds the maximal one ?
         if (master) write(uout,'(4(a,i5,2x))') 'particle type =', ipt,'iploc =', iploc, 'nppt(ipt) =', nppt(ipt)
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if

   end do

end subroutine SetRandomCylinderShell

!************************************************************************
!> \page moluser moluser.F90
!! **SetChainLinePMA**
!! *generate a line configuration for PMA chain particles*
!************************************************************************


subroutine SetChainLinePMA!(iptset) iptset is not used

   use CoordinateModule
   implicit none

   !integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainLinePMA'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, jp
   integer(4) :: n,i,j,k
   real(8)    :: dtranx(20),dtrany(20),dtranz(20)
   real(8)    :: bondloc
   real(8)    :: alpha
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   if (lclink) call Stop(txroutine, 'lclink is true', uout)
   ntry = ntrydef

   if (.not.first) return                                   ! should be called only once
   first =.false.

! setup of translational vectors for up to 20 chains

   if (nc>20) call Stop(txroutine, 'too many PMA chains', uout)

   n=20   ! empirically estimated
   k=0
      do i=1,5
         do j=1,4
            k=k+1
            dtranx(k)=(float(i)-0.5)*boxlen(1)/n*(-1.)**k
            dtrany(k)=(float(j)-0.5)*boxlen(2)/n*(-1.)**k
            dtranz(k)=(float(k)-0.5)*boxlen(3)/n
         enddo
      enddo


   alpha=acos(-1.d0/3.d0)/2.d0                              ! half of the tethraedral angle (109.47 deg)

   do ic = 1, nc                                            ! loop over chains
      ict = ictcn(ic)                                       ! chain type
      bondloc = bondscl(ict)*bond(ict)%eq                   ! bond length to be used
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)                             ! particle to be set

         do itry = 1, ntry                                  ! loop over attempts to set the particle
            if (iseg == 1) then                             ! a first segment
               ro(1,ip)=dtranx(ic)
               ro(2,ip)=dtrany(ic)
               ro(3,ip)=dtranz(ic)
            else                                            ! a remaining segment
               jp = ipnsegcn(iseg-1,ic)
               if(mod(iseg,2) == 0)then
                 ro(3,ip)=zero                              + dtranz(ic)  !start building chains at different locations in the MC cell
                 ro(2,ip)=bondloc*cos(alpha)                + dtrany(ic)
                 ro(1,ip)=ro(1,jp)+sqrt(2.d0/3.d0)*bondloc
               else
                 ro(3,ip)=zero                              + dtranz(ic)
                 ro(2,ip)=zero                              + dtrany(ic)
                 ro(1,ip)=ro(1,jp)+sqrt(2.d0/3.d0)*bondloc
               endif
            end if
            if (CheckPartOutsideBox(ip)) cycle              ! check if particle is outside the box
            if (CheckTooFoldedChain(ip, bondloc)) cycle     ! check if a too folded chain

            call SetPartOriLab(ori(1,1,ip))                 ! set particle orientation

            call SetAtomPos(ip,ip,.false.)                  ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap

            lpset(ip) =.true.                               ! configuration accepted
            exit
         end do

         if (itry > ntry) then                              ! number of attempts exceeds the maximal one ?
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            call Stop(txroutine, 'linear PMA chain configuration failed, itry > ntry', uout)
         end if

      end do

   end do

end subroutine SetChainLinePMA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **DumpUser**
!! *driver of user-provided routines for dumping*
!************************************************************************


subroutine DumpUser(iStage)

   use DumpModule
   implicit none
   integer(4), intent(in) :: iStage
   if (master) then
!     call BBDump(iStage)
!     call ChainCOMDump(iStage)
!     call ChainReeDump(iStage)
!     call DipMomTotDump(iStage)
     if (txuser == 'jasper_xy_coord') call xy_coordinates_Jasper(iStage)
   end if

end subroutine DumpUser

!************************************************************************
!> \page moluser moluser.F90
!! **BBDump**
!! *dumps three angles and two energies*
!************************************************************************


!     adapted for two benzene moleucules in water

!     phi          angle formed by x1' and x2' axis
!     theta1       angles formed by r(b-b) and x1'-axis
!     theta2       angle formed by r(b1-b2) and x2'-axis
!     u%twob(2)   bw energy
!     u%twob(3)   bb energy

subroutine BBDump(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   integer(4) :: m, ip1, ip2
   real(8) :: dum, dx, dy, dz, phi, theta1, theta2
   real(8), external :: angle_deg
   save

   select case (iStage)
   case (iWriteInput)

      call FileOpen(uuser, fuser, 'form/noread')

! ... get particle numbers of the benzene molecules

      ip1 = ipnpt(2)
      ip2 = ip1+1

   case (iBeforeSimulation)

      if (txstart == 'continue') then               ! advance fuser
         do m = 1, (nstep2/idump)*(nstep1beg-1)+1
            read(uuser,'(f10.4,3f12.6)') dum, dum, dum, dum, dum
         end do
      else if (txstart/= 'continue') then           ! dump initial value
         dx = ro(1,ip1)-ro(1,ip2)
         dy = ro(2,ip1)-ro(2,ip2)
         dz = ro(3,ip1)-ro(3,ip2)
         phi = angle_deg(ori(1,3,ip1),ori(2,3,ip1),ori(3,3,ip1),ori(1,3,ip2),ori(2,3,ip2),ori(3,3,ip2))
         theta1 = angle_deg(ori(1,3,ip1),ori(2,3,ip1),ori(3,3,ip1),-dx,-dy,-dz)
         theta2 = angle_deg(-dx,-dy,-dz,ori(1,3,ip2),ori(2,3,ip2),ori(3,3,ip2))
         write(uuser,'(f10.4,3f12.6)') phi, theta1, theta2, u%twob(2), u%twob(3)
      end if

   case (iSimulationStep)

      if (mod(istep2,idump) == 0) then            ! dump
         dx = ro(1,ip1)-ro(1,ip2)
         dy = ro(2,ip1)-ro(2,ip2)
         dz = ro(3,ip1)-ro(3,ip2)
         phi = angle_deg(ori(1,3,ip1),ori(2,3,ip1),ori(3,3,ip1),ori(1,3,ip2),ori(2,3,ip2),ori(3,3,ip2))
         theta1 = angle_deg(ori(1,3,ip1),ori(2,3,ip1),ori(3,3,ip1),-dx,-dy,-dz)
         theta2 = angle_deg(-dx,-dy,-dz,ori(1,3,ip2),ori(2,3,ip2),ori(3,3,ip2))
         write(uuser,'(f10.4,3f12.6)') phi, theta1, theta2, u%twob(2), u%twob(3)
      end if

   case (iAfterSimulation)

      close(uuser)

   end select

end subroutine BBDump

!************************************************************************
!> \page moluser moluser.F90
!! **ChainCOMDump**
!! *dump center of mass of chains*
!************************************************************************


subroutine ChainCOMDump(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   real(8), allocatable, save :: rcom(:,:)
   real(8)    :: dum
   integer(4) :: m, idum, ic
   integer(4), save :: idumplocal
   type(chainprop_var) :: ChainProperty

   select case (iStage)
   case (iReadInput)

      call FileOpen(upos, fpos, 'unform/noread')    ! borrow upos and fpos
      idumplocal = 1

   case (iBeforeSimulation)

      if (.not.allocated(rcom)) then
         allocate(rcom(3,nc))
         rcom = 0.0E+00
      end if

      if (txstart /= 'continue') then                  ! dump initial value
         call GetChainData
         write(upos) idum, boxlen(1:3)
         write(upos) rcom(1:3,1:nc)
         call FileFlush(upos)
      else if (txstart == 'continue') then             ! advance fpos
         do m = 1, (nstep2/idumplocal)*(nstep1beg-1)+1
         read(upos) idum, dum, dum, dum
         read(upos) (dum,dum,dum,ic = 1,nc)
         end do
      end if

   case (iSimulationStep)

      if ((lsim .and. mod(istep2,idumplocal) == 0) .or. lana) then  ! dump
         call GetChainData
         write(upos) idum, boxlen(1:3)
         write(upos) rcom(1:3,1:nc)
         call FileFlush(upos)
      end if

   case (iAfterMacrostep)

      deallocate(rcom)

      call FileFlush(upos)

   case (iAfterSimulation)

      close(upos)

   end select

contains

!........................................................................

subroutine GetChainData
   do ic = 1, nc
      call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
      call CalcChainProperty(ic, vaux, ChainProperty)
      rcom(1:3,ic) = ChainProperty%ro(1:3)
   end do
end subroutine GetChainData

!........................................................................

end subroutine ChainCOMDump

!************************************************************************
!> \page moluser moluser.F90
!! **ChainReeDump**
!! *dump end-to-end separations of chains*
!************************************************************************


subroutine ChainReeDump(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   real(8), allocatable :: ree(:,:)
   real(8)    :: dum
   integer(4) :: m, ic
   integer(4), save :: idumplocal
   type(chainprop_var) :: ChainProperty

   select case (iStage)
   case (iReadInput)

      call FileOpen(uuser, fuser, 'unform/noread')
      idumplocal = 1

   case (iBeforeSimulation)

      if (.not.allocated(ree)) then
         allocate(ree(3,nc))
         ree = 0.0E+00
      end if

      if (txstart /= 'continue') then                  ! dump initial value
         call GetChainData
         write(uuser) ree(1:3,1:nc)
         call FileFlush(uuser)
      else if (txstart == 'continue') then             ! advance fuser
         do m = 1, (nstep2/idumplocal)*(nstep1beg-1)+1
         read(uuser) (dum,dum,dum,ic = 1,nc)
         end do
      end if

   case (iSimulationStep)

      if ((lsim .and. mod(istep2,idumplocal) == 0) .or. lana) then  ! ... dump
         call GetChainData
         write(uuser) ree(1:3,1:nc)
         call FileFlush(uuser)
      end if

   case (iAfterMacrostep)

      deallocate(ree)

      call FileFlush(uuser)

   case (iAfterSimulation)

      close(uuser)

   end select

contains

!........................................................................

subroutine GetChainData
   do ic = 1, nc
      call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
      call CalcChainProperty(ic, vaux, ChainProperty)
      ree(1:3,ic) = ChainProperty%ree(1:3)
   end do
end subroutine GetChainData

!........................................................................

end subroutine ChainReeDump

!************************************************************************
!> \page moluser moluser.F90
!! **DipMomTotDump**
!! *dump the total dipole moment*
!************************************************************************


subroutine DipMomTotDump(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   real(8)    :: dum, dipmomtot(1:3)
   integer(4) :: m, ic
   integer(4), save :: idumplocal

   select case (iStage)
   case (iReadInput)

      call FileOpen(uuser, fuser, 'unform/noread')
!     idumplocal = 1                                ! Jocke
      idumplocal = 100

   case (iBeforeSimulation)

      if (txstart /= 'continue') then                  ! dump initial value
         call DipMomTotDumpSub
      else if (txstart == 'continue') then             ! advance fuser
         do m = 1, (nstep2/idumplocal)*(nstep1beg-1)+1
            read(uuser) (dum,dum,dum,ic = 1,nc)
         end do
      end if

   case (iSimulationStep)

      if ((lsim .and. mod(istep2,idumplocal) == 0) .or. lana) then  ! dump
         call DipMomTotDumpSub
      end if

   case (iAfterMacrostep)

   case (iAfterSimulation)

      close(uuser)

   end select

contains

!........................................................................

subroutine DipMomTotDumpSub
   call CalcSysDipMom(dipmomtot)
   write(uuser) dipmomtot(1:3)
   call FileFlush(uuser)
end subroutine DipMomTotDumpSub

!........................................................................

end subroutine DipMomTotDump

!************************************************************************
!> \page moluser moluser.F90
!! **xy_coordinates_Jasper**
!! *read x,y coordinates*
!************************************************************************


subroutine xy_coordinates_Jasper(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='xy_coordinates_Jasper'
   integer(4), save :: npixel(1:2), ip, ipt, iploc, ierr
   real(8), save, allocatable :: xy(:,:)

   select case (iStage)
   case (iReadInput)

      call FileOpen(uuser, fuser, 'form/noread')

   case (iBeforeSimulation)

     allocate(xy(1:2,1:np))
     xy = 0.0E+00

   case (iSimulationStep)

     read(uuser,*) npixel(1)                        ! number of pixels in x-dir
     npixel(2)=npixel(1)                            ! number of pixels in y-dir
     ip = 0
     do ipt = 1, npt                                ! loop over particle types
        read(uuser,*)
        read(uuser,*, iostat = ierr) nppt(ipt)
        if (ierr /= 0) call WriteIOStat(txroutine, 'read of number of particles failed', ierr, 2, uout)
        do iploc = 1, nppt(ipt)                     ! loop over particles of particle type ipt
           ip = ip + 1
           read(uuser,*,iostat = ierr) xy(1,ip), xy(2,ip)
           write(*,*) 'ip ok', ip
           if (ierr /= 0) call WriteIOStat(txroutine, 'read of particles failed', ierr, 2, uout)
        end do
     end do

! ... copy to global particle position

     call SetObjectParam1
     call SetObjectParam2
     ro(1,1:np) = xy(1,1:np) - npixel(1)/2
     ro(2,1:np) = xy(2,1:np) - npixel(2)/2
     ro(3,1:np) = Zero
     r(1:3,1:np) = ro(1:3,1:np)

   case (iAfterMacrostep)

   case (iAfterSimulation)

      close(uuser)
      write(*,*) 'iStage', iStage
      write(*,'(i5,2f10.2)') (iptpn(ip), xy(1:2,ip), ip = 1, np)

   end select

end subroutine xy_coordinates_Jasper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **GroupUser**
!! *driver of user-provided routines for particle group division*
!************************************************************************


subroutine GroupUser(iStage, m, txtype, lsetconf)

   use MolModule

   integer(4),    intent(in)  :: iStage     != 1 call for setting ngr
                                           ! ngr(m)       number of groups
                                           != 2 call for setting iptgr, and grvar()%label
                                           ! iptgr(igr,m)       partcle type associated with group igr
                                           ! grvar(m,igr)%label label of members in gropup igr
                                           ! = 4 call for calculating igrpn and grvar()%value
                                           ! igrpn(i,m)         group number of particle i
                                           ! grvar(m,igr)%value number of particles belonging to group igr
   integer(4),    intent(in)  :: m         ! = 1 classification of reference particles
                                           ! = 2 classification of field particles
   character(20), intent(in)  :: txtype(2) ! lable of type of group division
   logical,       intent(out) :: lsetconf  ! =.true. if match and use of user supplied group section
                                           ! =.false. otherwise

   lsetconf =.true.
   if (txtype(m) == 'benzene_in_water_1') then
      call GroupBW1(iStage, m)
   else if (txtype(m) == 'benzene_in_water_2') then
      call GroupBW2(iStage, m)
   else if (txtype(m) == 'benzene_in_water_3') then
      call GroupBW3(iStage, m)
   else if (txtype(m) == 'bb(aq)_w1') then
      call GroupBBW1(iStage, m)
   else if (txtype(m) == 'bb(aq)_w2') then
      call GroupBBW2(iStage, m)
   else if (txtype(m) == 'bb(aq)_w3') then
      call GroupBBW3(iStage, m)
   else if (txtype(m) == 'bb(aq)_b') then
      call GroupBBW4(iStage, m)
   else if (txtype(m) == 'surf') then
      call GroupSurface(iStage, m)
   else if (txtype(m) == 'surf2') then
      call GroupSurface2(iStage, m)
   else if (txtype(m) == 'ads1') then   ! previously 'niklas_ads1'
      call GroupAds1(iStage, m)
   else if (txtype(m) == 'ads_layer1_ramp') then
      call GroupAds_layer1_ramp(iStage, m)
   else if (txtype(m) == 'networkgenerations') then
      call GroupNetworkGenerations(iStage, m)
   else if (txtype(m) == 'weakcharge') then
      call GroupWeakCharge(iStage, m)
   else
      lsetconf =.false.
   end if

end subroutine GroupUser

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBW1**
!! *group division for one benzene molecule in water*
!************************************************************************


subroutine GroupBW1(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr
   real(8)    :: dx, dy, dz, r2, rmax1, rmax2

   select case (iStage)
   case (iReadInput)
      ngr(m) = 4
   case (iWriteInput)
      iptgr(1,m) = 2
      iptgr(2,m) = 1
      iptgr(3,m) = 1
      iptgr(4,m) = 1
      grvar(igrpnt(m,1))%label = 'benzene'
      grvar(igrpnt(m,2))%label = '1:st hydration shell'
      grvar(igrpnt(m,3))%label = '2:nd hydration shell'
      grvar(igrpnt(m,4))%label = 'other water molecules'
   case (iSimulationStep)
      rmax1 = 6.0
      rmax2 = 9.0
      do ip = 1, np
         if (ip == np) then
            igr = 1
         else
            igr = 2
            dx = ro(1,ip)-ro(1,np)
            dy = ro(2,ip)-ro(2,np)
            dz = ro(3,ip)-ro(3,np)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rmax1**2) igr = 3
            if (r2 > rmax2**2) igr = 4
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBW1

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBW2**
!! *group division for one benzene molecule in water*
!************************************************************************


subroutine GroupBW2(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr
   real(8)    :: dx, dy, dz, r2, acp, angcos

   select case (iStage)
   case (iReadInput)
      ngr(m) = 3
   case (iWriteInput)
      iptgr(1,m) = 2
      iptgr(2,m) = 1
      iptgr(3,m) = 1
      grvar(igrpnt(m,1))%label = 'benzene'
      grvar(igrpnt(m,2))%label = 'water top'
      grvar(igrpnt(m,3))%label = 'water side'
   case (iSimulationStep)
      do ip = 1, np
         if (ip == np) then
            igr = 1
         else
            igr = 0
            dx = ro(1,ip)-ro(1,np)
            dy = ro(2,ip)-ro(2,np)
            dz = ro(3,ip)-ro(3,np)
            call PBCr2(dx,dy,dz,r2)
            acp = angcos(ori(1,3,np),ori(2,3,np),ori(3,3,np),-dx,-dy,-dz)
            if (abs(acp) > 0.8660254 .and. r2 < 27.65) igr = 2
            if (abs(acp) < 0.5000000 .and. r2 < 40.32) igr = 3
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBW2

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBW3**
!! *group division for one benzene molecule in water*
!************************************************************************


subroutine GroupBW3(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr
   real(8)    :: dx, dy, dz, r2, acp, absacp, angcos
   real(8)    :: t0, t1, t2, t3, t4

   select case (iStage)
   case (iReadInput)
      ngr(m) = 5
   case (iWriteInput)
      iptgr(1,m) = 2
      iptgr(2,m) = 1
      iptgr(3,m) = 1
      iptgr(4,m) = 1
      iptgr(5,m) = 1
      grvar(igrpnt(m,1))%label = 'benzene'
      grvar(igrpnt(m,2))%label = 'water reg 1'
      grvar(igrpnt(m,3))%label = 'water reg 2'
      grvar(igrpnt(m,4))%label = 'water reg 3'
      grvar(igrpnt(m,5))%label = 'water reg 4'
   case (iSimulationStep)
      t0 = 1.0
      t1 = 0.8660254
      t2 = 0.7071068
      t3 = 0.5000000
      t4 = 0.0000000
      do ip = 1, np
         if (ip == np) then
            igr = 1
         else
            igr = 0
            dx = ro(1,ip)-ro(1,np)
            dy = ro(2,ip)-ro(2,np)
            dz = ro(3,ip)-ro(3,np)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > 100) cycle
            acp = angcos(ori(1,3,np),ori(2,3,np),ori(3,3,np),-dx,-dy,-dz)
            absacp = abs(acp)
            if (absacp < t0 .and. absacp > t1) igr = 2
            if (absacp < t1 .and. absacp > t2) igr = 3
            if (absacp < t2 .and. absacp > t3) igr = 4
            if (absacp < t3 .and. absacp > t4) igr = 5
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBW3

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBBW1**
!! *group division for two benzene molecules in water*
!************************************************************************


subroutine GroupBBW1(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, jp, igr
   real(8)    :: dx, dy, dz, r2, rmax1, rmax2

   select case (iStage)
   case (iReadInput)
      ngr(m) = 3
   case (iWriteInput)
      iptgr(1,m) = 1
      iptgr(2,m) = 1
      iptgr(3,m) = 1
      grvar(igrpnt(m,1))%label = '1:st hydration shell'
      grvar(igrpnt(m,2))%label = '2:nd hydration shell'
      grvar(igrpnt(m,3))%label = 'other water molecules'
   case (iSimulationStep)
      rmax1 = 6.0
      rmax2 = 9.0
      do ip = 1, np
         if (ip == np-1 .or. ip == np) then
            igr = 0
         else
            igr = 3
            do jp = np-1, np
               dx = ro(1,ip)-ro(1,jp)
               dy = ro(2,ip)-ro(2,jp)
               dz = ro(3,ip)-ro(3,jp)
               call PBCr2(dx,dy,dz,r2)
               if (r2 < rmax2**2) igr = min(int(igr),2)
               if (r2 < rmax1**2) igr = 1
            end do
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBBW1

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBBW2**
!! *group division for two benzene molecules in water*
!************************************************************************


subroutine GroupBBW2(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, jp, igr
   real(8)    :: dx, dy, dz, r2, acp, angcos, absz, rho2

   select case (iStage)
   case (iReadInput)
      ngr(m) = 3
   case (iWriteInput)
      iptgr(1,m) = 1
      iptgr(2,m) = 1
      iptgr(3,m) = 1
      grvar(igrpnt(m,1))%label = 'water between'
      grvar(igrpnt(m,2))%label = 'water top'
      grvar(igrpnt(m,3))%label = 'water side'
   case (iSimulationStep)
      do ip = 1, np
         if (ip == np-1 .or. ip == np) then
            igr = 0
         else
            igr = 0
            do jp = np-1, np
               dx = ro(1,ip)-ro(1,jp)
               dy = ro(2,ip)-ro(2,jp)
               dz = ro(3,ip)-ro(3,jp)
               call PBCr2(dx,dy,dz,r2)
               acp = angcos(ori(1,3,jp),ori(2,3,jp),ori(3,3,jp),-dx,-dy,-dz)
               if (abs(acp) > 0.8660254 .and. r2 < 27.65) igr = 2
               if (abs(acp) < 0.5000000 .and. r2 < 40.32) igr = 3
            end do
            absz = abs(ro(3,ip))
            rho2 = ro(1,ip)**2+ro(2,ip)**2
            if (absz < 3.0 .and. rho2 < 6.25) igr = 1
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBBW2

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBBW3**
!! *group division for two benzene molecules in water*
!************************************************************************


subroutine GroupBBW3(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr
   real(8)    :: absz, rho2

   select case (iStage)
   case (iReadInput)
      ngr(m) = 1
   case (iWriteInput)
      iptgr(1,m) = 1
      grvar(igrpnt(m,1))%label = 'water between'
   case (iSimulationStep)
      do ip = 1, np
         if (ip == np-1 .or. ip == np) then
            igr = 0
         else
            igr = 0
            absz = abs(ro(3,ip))
            rho2 = ro(1,ip)**2+ro(2,ip)**2
            if (absz < 3.0 .and. rho2 < 6.25) igr = 1
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBBW3

!************************************************************************
!> \page moluser moluser.F90
!! **GroupBBW4**
!! *group division for two benzene molecules in water*
!************************************************************************


subroutine GroupBBW4(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr

   select case (iStage)
   case (iReadInput)
      ngr(m) = 2
   case (iWriteInput)
      iptgr(1,m) = 2
      iptgr(2,m) = 2
      grvar(igrpnt(m,1))%label = 'benzene 1'
      grvar(igrpnt(m,2))%label = 'benzene 2'
   case (iSimulationStep)
      do ip = 1, np
         if (ip == np-1) then
            igr = 1
         else if (ip == np) then
            igr = 2
         else
            igr = 0
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupBBW4

!************************************************************************
!> \page moluser moluser.F90
!! **GroupSurface**
!! *group division for particles of type ipt = 2  at surface and in bulk*
!************************************************************************


subroutine GroupSurface(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr

   select case (iStage)
   case (iReadInput)
      ngr(m) = 2
   case (iWriteInput)
      iptgr(1,m) = 2
      iptgr(2,m) = 2
      grvar(igrpnt(m,1))%label = 'prot_surf'
      grvar(igrpnt(m,2))%label = 'prot_bulk'
   case (iSimulationStep)
      do ip = 1, np
         if (iptpn(ip) /= 2) cycle
         if (boxlen2(3)-abs(ro(3,ip)) < 0) then
            igr = 0
         else if (boxlen2(3)-abs(ro(3,ip)) < 23.54) then
            igr = 1
         else
            igr = 2
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupSurface

!************************************************************************
!> \page moluser moluser.F90
!! **groupsurface2**
!! *group division for particles only at surface*
!************************************************************************


subroutine GroupSurface2(iStage, m)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: m

   integer(4) :: ip, igr

   select case (iStage)
   case (iReadInput)
      ngr(m) = 1
   case (iWriteInput)
      iptgr(1,m) = 2
      grvar(igrpnt(m,1))%label = 'prot_surf'
   case (iSimulationStep)
      do ip = 1, np
         if (iptpn(ip) /= 2) cycle
         if (boxlen2(3)-abs(ro(3,ip)) < 23.54) then
            igr = 1
         else
            igr = 0
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupSurface2

!************************************************************************
!> \page moluser moluser.F90
!! **GroupAds1**
!! *group division for particles belonging to adsorbed chains*
!************************************************************************


subroutine GroupAds1(iStage, m)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: m

   type(adscond_var), parameter :: adscond = adscond_var('plane_xy','-|+',6.0d0)
   logical,   allocatable, save :: ladschain(:)    ! .true. if chain is adsorbed
   logical,   allocatable, save :: ladsseg(:)      ! .true. if segment is adsorbed
   integer(4) :: ip, igr, ic

   if (.not.allocated(ladschain)) then
      allocate(ladschain(nc), ladsseg(maxval(npct(1:nct))))
      ladschain = .false.
      ladsseg = .false.
   end if

   select case (iStage)
   case (iReadInput)
      ngr(m) = 1
   case (iWriteInput)
      iptgr(1,m) = 1
      grvar(igrpnt(m,1))%label = 'part in ads chains'
   case (iSimulationStep)
      do ic = 1, nc
         call CheckAdsChainSeg(ic, adscond, ladschain(ic), ladsseg)
      end do
      do ip = 1, np
         ic = icnpn(ip)
         if (ladschain(ic)) then
            igr = 1
         else
            igr = 0
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupAds1

!************************************************************************
!> \page moluser moluser.F90
!! **GroupAds_layer1_ramp**
!! *group division for adsorbed (layer 1) particles*
!************************************************************************


subroutine GroupAds_layer1_ramp(iStage, m)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: m

   integer(4) :: ip, igr

   real(8)         , save :: adsoffset = One

   select case (iStage)
   case (iReadInput)
      ngr(m) = 1
   case (iWriteInput)
      iptgr(1,m) = 1
      grvar(igrpnt(m,1))%label = 'ads particles (layer 1)'
   case (iSimulationStep)
      do ip = 1, np
         if (boxlen2(3)-abs(ro(3,ip)) < adsoffset) then
            igr = 1
         else
            igr = 0
         end if
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupAds_layer1_ramp

!************************************************************************
!> \page moluser moluser.F90
!! **GroupNetworkGenerations**
!! *group division according to generations of non-periodic network*
!************************************************************************


subroutine GroupNetworkGenerations(iStage,m)

   use MolModule
   implicit none

   integer(4),    intent(in)     :: iStage
   integer(4),    intent(in)     :: m

   character(len=*), parameter   :: txroutine = 'GroupNetworkGenerations'

   integer(4), allocatable, save :: igencn(:)    ! chain                   -> its generation
   integer(4), allocatable, save :: icnclpn(:,:) ! crosslink and particle  -> crosslinked chain
   integer(4), allocatable, save :: ipnclcn(:,:) ! crosslink and chain     -> crosslinked particle
   integer(4), allocatable, save :: nclcn(:)     ! chain number            -> its number of cross-links

   integer(4)                    :: ip, igr, ic, ict

   integer(4), parameter         :: inwt = 1 ! Currently works for systems with one network type only

   character(10)                 :: str

   select case (iStage)

   case (iReadInput)

      ! ... for the allocations in subroutine Group: Set ngr(m) high enough
      ! ... the actual number of generations cannot exceed nc/4 for networks as set by SetNetwork
      ! ... add 1 in order to account for the cross-links of the network
      ngr(m) = nc/4 + 1

      if (.not.lnetwork) then
         call Warn(txroutine,'ref/field == ''networkgenerations'' .and. .not.lnetwork',uout)
      else if (nnw > 1) then
         call Warn(txroutine,'ref/field == ''networkgenerations'' .and. nnw > 1',uout)
      end if

      ! ... alloctate cross-link related pointers
      if (.not.allocated(icnclpn)) then
         allocate(icnclpn(maxvalnbondcl,np_alloc))
         icnclpn = 0
      end if
      if (.not.allocated(nclcn)) then
         allocate(nclcn(nc))
         nclcn = 0
      end if
      if (.not.allocated(ipnclcn)) then
         allocate(ipnclcn(2,nc))
         ipnclcn = 0
      end if

   case (iWriteInput)

      ! ... Evaluate required information
      call SetNetworkGenerationPointer

      ! ... Determine ngr(m) and assign generations to chains
      if(.not.allocated(igencn)) allocate(igencn(0:nc))
      igencn(0:nc) = 0
      do ic = 1, nc
         if (nclcn(ic) == 1) then ! dangling chain ic
            call AssignChainGeneration(ic,0,igencn)
         end if
      end do

      ! ... add 1 in order to account for the cross-links of the network as an own group
      ngr(m) = maxval(igencn(1:nc)) + 1

      ! ... Determine iptgr(igr,m)
      do ict = 1, nct
         if (inwtct(ict) == inwt) exit
      end do
      iptgr(1:ngr(m)-1,m) = iptpn(ipnsegcn(1,icnct(ict))) ! particle type constituting chains of network
      iptgr(ngr(m),m) = iptclnwt(inwt) ! particle type constituting cross-links of network

      ! ... Determine grvar(igrpnt(m,igr))%label
      str = '      '
      do igr = 1, ngr(m)
         write(str,'(i10)') igr
         grvar(igrpnt(m,igr))%label = "nw-gen:"//trim(adjustl(str))
      end do

      ! Set igrpn(ip,m) in order to know it for image generation - Group assignment of particles is fixed throughout the simulation
      do ip = 1, np
         igrpn(ip,m) = merge(ngr(m), igencn(icnpn(ip)), iptpn(ip) == iptclnwt(inwt))
         grvar(igrpnt(m,igrpn(ip,m)))%value = grvar(igrpnt(m,igrpn(ip,m)))%value + 1
      end do

   case (iSimulationStep)

      do ip = 1, np
         ! particle ip is either:
         ! -> part of a chain of network (igrpn is then the ID of its respective chain generation igrpn = igencn)
         ! -> part of the group of cross-links of the network (igrpn is then the highest group number: igrpn = ngr(m))
         ! -> no part of the network (igrpn is then the ID of its respective chain generations: igrpn = igencn = 0)
         igrpn(ip,m) = merge(ngr(m), igencn(icnpn(ip)), iptpn(ip) == iptclnwt(inwt))
         grvar(igrpnt(m,igrpn(ip,m)))%value = grvar(igrpnt(m,igrpn(ip,m)))%value + 1
      end do

   end select

contains

!........................................................................

recursive subroutine AssignChainGeneration(ic,jc,igencn)

   implicit none

   integer(4), intent(in)     :: ic, jc
   integer(4), intent(inout)  :: igencn(0:nc)

   integer(4)  :: jccl
   integer(4)  :: icl, ibondcl
   integer(4)  :: ipnode

   if ((igencn(ic) > (igencn(jc)+1)) .or. (igencn(ic) == 0)) then
      igencn(ic) = igencn(jc) + 1
      do icl = 1, nclcn(ic)
         ipnode = ipnclcn(icl,ic)
         do ibondcl = 1, nbondcl(ipnode)
            jccl = icnclpn(ibondcl,ipnode)
            if (ic == jccl) cycle ! ... don't go back to where you came from
            call AssignChainGeneration(jccl,ic,igencn)
         end do
      end do
   else
      continue
   end if

end subroutine AssignChainGeneration

!........................................................................

subroutine SetNetworkGenerationPointer

   use MolModule
   implicit none

   integer(4)  :: ip
   integer(4)  :: ibondcl
   integer(4)  :: ic

   do ip = 1, np
      if ((nbondcl(ip) > 0) .and. (icnpn(ip) == 0)) then ! ip is a node: It got crosslinks but it's no part of a chain
         do ibondcl = 1, nbondcl(ip)
            ic = icnpn(bondcl(ibondcl,ip))
            nclcn(ic) = nclcn(ic) + 1
            icnclpn(ibondcl,ip) = ic
            ipnclcn(nclcn(ic),ic) = ip
         end do
      end if
   end do

end subroutine SetNetworkGenerationPointer

end subroutine GroupNetworkGenerations

!************************************************************************
!> \page moluser moluser.F90
!! **GroupWeakCharge**
!! *group division according to titratable species - respective division in*
!************************************************************************

! ... charged and uncharged state
! ...
! ... currently works for monoatomic systems only

subroutine GroupWeakCharge(iStage,m)

   use MolModule
   implicit none

   integer(4),    intent(in)     :: iStage
   integer(4),    intent(in)     :: m

   character(len=*), parameter   :: txroutine = 'GroupWeakCharge'

   character(len=9), parameter   :: txchargestate(2) = [ 'charged  ', 'uncharged' ]

   integer(4), allocatable, save :: igrref(:,:) ! group reference for fast assignment of particles to groups

   integer(4)                    :: ichargestate, ip, ipt, igr

   select case (iStage)

   case (iReadInput)

      if (.not.lweakcharge) then
         call Stop(txroutine,'ref/field == ''weakcharge'' .and. .not.lweakcharge',uout)
      else if (.not.lmonoatom) then
         call Stop(txroutine,'ref/field == ''weakcharge'' .and. .not.lmonoatom',uout)
      end if

      ngr(m) = 0
      do ipt = 1, npt   !                   weak charge?            counterion?
         ngr(m) = merge(ngr(m)+2, ngr(m)+1, latweakcharge(ipt) .or. any(jatweakcharge(1:npt) == ipt))
      end do

      if (.not.allocated(igrref)) then
         allocate(igrref(2,npt))
         igrref = 0
      end if

   case (iWriteInput)

      ! ... Determine iptgr(igr,m), grvar(igrpnt(m,igr))%label
      igr = 0
      do ipt = 1, npt
         if (latweakcharge(ipt) .or. any(jatweakcharge(1:npt) == ipt)) then ! weak charge or counterion
            do ichargestate = 1, 2
               igr = igr + 1
               iptgr(igr,m) = ipt
               grvar(igrpnt(m,igr))%label = trim(adjustl(txchargestate(ichargestate)))//' '//&
                                           &trim(adjustl(txpt(ipt)))
               igrref(ichargestate,ipt) = igr
            end do
         else ! no charge or fixed charge
            igr = igr + 1
            iptgr(igr,m) = ipt
            grvar(igrpnt(m,igr))%label = trim(adjustl(txpt(ipt)))
            igrref(1:2,ipt) = igr
         end if
      end do

   case (iSimulationStep)

      igrpn(1:np,m) = 0
      do ipt = 1, npt
         do ip = ipnpt(ipt), ipnpt(ipt)+nppt(ipt)-1
            igrpn(ip,m) = merge(igrref(1,ipt), igrref(2,ipt), laz(ip))
            grvar(igrpnt(m,igrpn(ip,m)))%value = grvar(igrpnt(m,igrpn(ip,m)))%value + 1
         end do
      end do

   end select

end subroutine GroupWeakCharge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **StaticUser**
!! *driver of user-provided routines for static analysis*
!************************************************************************


subroutine StaticUser(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

! ... place routines for serial computation here

   if (master) then
      if (txuser == 'scalardemo') call ScalarDemo1(iStage, u%tot, uout)
      if (txuser == 'scalardemo') call ScalarDemo2(iStage, uout)
      if (txuser(1:6) == 'jasper') call BondOrder(iStage)
      if (txuser == 'anna') call ChainBeadCylContact(iStage)

      if (txuser(1:4) == 'comp') call DoComplexation(iStage)

      if(txuser == 'ellipsoid') then     ! ellipsoid Erik W
!        call Min(iStage)
          call S1(iStage)
!        call S2(iStage)
!        call Q(iStage)
!        call C2(iStage)
          call CylDistFunc(iStage)
          call SFPBC2D(iStage)
!        call SFPBC2DNoAv(iStage)
      end if

      if(txuser == 'twocomponentmodel_corr') then
         if (lbcsph) then                 ! spherical geomentry, one macroion
            call PosMeanSph(iStage)        ! means of positions
            call MacroionOneSph(iStage)    ! projection of position onto the surface of one particle
         else if (lbccyl) then            ! cylindrical geomentry, two macroion
            call MacroionTwoCyl(iStage)    ! various df for system of two macroions + ions, cyl geometry
            call MeanElFieldZCyl(iStage)   ! mean electrostastic field and force in the z-direction
         end if
      end if

      if (txuser == 'chainads') then      ! adsorption of chains
         if (lmc) then                    ! Monte Carlo
            call SCDF(iStage)              ! chain distribution functions
            call AdsRadGyr(iStage)         ! radius of gyration df of adsorbed and adsorbing chains
            call AdsBondOrder(iStage)      ! bond order df of adsorbed chains
         else if (lbd) then               ! Brownian dynamics
            call AdsPropDyn(iStage)        ! time dependence of adsorbed chain properties
            call AdsEventDyn(iStage)       ! write time and occurrance of adsorbtion events
            if (iStage == iAfterSimulation) call AdsExam ! analysis of written time and occurrance of adsorbtion event
         end if
      end if

      if (txuser == 'daniel' .and. lclink) then        ! Daniel Angelescu
         call COMBAver(iStage)         ! averages of comb chain quantities l
         call SPDF_COMB(iStage)        ! radial distrib. vs. COM
         call COMB_DF(iStage)          ! type distribution functions of comb polymer (coarse model)
         call OCF(iStage)              !  orientational correlation function(OCF)  for more chains
         call OCF_DF(iStage)           ! distribution functions of OCF at severasl separations along the chain contour
         call ChainBeadBeadContact(iStage)  ! distance between the backbone  ict =1 and linear chain ict =3
         if (txbc == 'xyz' .or. txbc == 'xy') then
            call SFPBC_COMB(iStage)    ! partial structure factor; loop over hierarchy, stored in ipt =2, works only for one comb
         else
            call Stop('SFDriver','Unsupported boundary condition',uout)
         end if
      end if

      if (txuser(1:11) == 'xyprojectdf') call XYProjectDF(iStage) ! projection onto z' = 0 plane

      if (txuser == 'jurij') call Z_DF_Slit(iStage)     ! Jurij Rescic  distribution in a slit
      if (txuser == 'rudi') call ElMom(iStage)     ! Rudi Podgorni  electrostatic moments

   end if

! ... place routines for parallel computation here

   if (txuser(1:9) == 'md_dipole') then   ! Dipole project with Gunnar Karlström
       call DomainDriver(iStage)           ! domain analysis based on Kirkwoods gk-factor
   endif

end subroutine StaticUser

!************************************************************************
!> \page moluser moluser.F90
!! **ScalarDemo1**
!! *calculate an average of a scalar quantity*
!************************************************************************


subroutine ScalarDemo1(iStage, value, unit)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), intent(in) :: iStage
   real(8),    intent(in) :: value     ! value of the scalar quantity
   integer(4), intent(in) :: unit

   character(80), parameter :: txheading ='scalar calc: average potential energy'
   integer(4)   , parameter :: nvar = 1
   type(scalar_var), save :: var(nvar)

   if (slave) return

   select case (iStage)
   case (iReadInput)
      var(1)%label = 'potential energy               = '
      var(1)%norm = One/np
   case (iBeforeSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
   case (iBeforeMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
   case (iSimulationStep)
      var(1)%value = value
      call ScalarSample(iStage, 1, nvar, var)
   case (iAfterMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      if (lsim .and. master) write(ucnf) var
   case (iAfterSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, unit)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', unit)
   end select

end subroutine ScalarDemo1

!************************************************************************
!> \page moluser moluser.F90
!! **ScalarDemo2**
!! *calculate average of scalar quantities*
!************************************************************************


subroutine ScalarDemo2(iStage,unit)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: unit

   character(40), parameter :: txroutine ='ScalarDemo2'
   character(80), parameter :: txheading ='scalar calc: average potential energy'
   integer(4)   , parameter :: nvar = 2
   type(scalar_var), save :: var(nvar)
   real(8)                :: value

!  namelist /nmlScalarDemo/ xxxx

   if (slave) return

!   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)
!     rewind(uin)
!     read(uin,nmlScalarDemo)
      var(1)%label = 'scalar 1:  1*upot               = '
      var(2)%label = 'scalar 2:  2*upot               = '
      var%norm = One/np
   case (iBeforeSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
   case (iBeforeMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
   case (iSimulationStep)
      value = u%tot                            ! evaluate value of scalar quantities
      var(1)%value = 1*value
      var(2)%value = 2*value
      call ScalarSample(iStage, 1, nvar, var)
   case (iAfterMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      if (lsim .and. master) write(ucnf) var
!     call WriteHead(2, txheading, unit)
!     call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', unit)
   case (iAfterSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, unit)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', unit)
   end select

end subroutine ScalarDemo2

!************************************************************************
!> \page moluser moluser.F90
!! **BondOrder**
!! *calculate an average of a scalar quantity*
!************************************************************************


subroutine BondOrder(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='BondOrder'
   character(80), parameter :: txheading ='bond order'
   integer(4)   , parameter :: nvar = 7
   integer(4)   , parameter :: lmax = nvar-1
   type(scalar_var), save :: var(nvar)
   real(8), save :: npartsum, rcutoff
   real(8) :: dx, dy, dz, r2, theta, phi, ql(0:nvar), thetann(50), phinn(50), qlsum(0:nvar)
   integer(4), save :: npart, nframe
   integer(4) :: m, l, ip, jp, nneigh, iploc
   complex(8) :: cclm, qlm


   if (slave) return

   select case (iStage)
   case (iReadInput)
      var(1)%label = 'l=0'
      var(2)%label = 'l=1'
      var(3)%label = 'l=2'
      var(4)%label = 'l=3'
      var(5)%label = 'l=4'
      var(6)%label = 'l=5'
      var(7)%label = 'l=6'
      var(1:nvar)%norm = One
      rcutoff = 2.8*radat(1)                        ! define cut-off length
   case (iBeforeSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
   case (iBeforeMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
   case (iSimulationStep)

      qlsum = 0
      npart = 0
      nframe = nframe + 1
      do ip=1,np
          nneigh = 0
!         if(abs(ro(1,ip))>boxlen2(1)-3*radat(1)) cycle
!         if(abs(ro(2,ip))>boxlen2(2)-3*radat(1)) cycle
!         if(abs(ro(3,ip))>boxlen2(3)-3*radat(1)) cycle

         do jp=1,np
            if (jp==ip) cycle                          ! exclude center particle
            dx = ro(1,jp)-ro(1,ip)                     ! define x-distance
            dy = ro(2,jp)-ro(2,ip)                     ! define y-distance
            dz = ro(3,jp)-ro(3,ip)                     ! define z-distance
            call PBC(dx,dy,dz)                         ! applying periodic boundary condition
            r2 = dx**2+dy**2+dz**2                     ! calculate cathesian distance
            if (r2 > rcutoff**2) cycle                 ! exclude large distance particles
            nneigh = nneigh + 1
            if (nneigh > 50) call stop(txroutine, 'nneigh > 50', uout)
            call CarToSph('rad',dx,dy,dz,r2,theta,phi) ! calculate spherical coordinates
            thetann(nneigh)=theta
            phinn(nneigh)=phi
         end do
         write(*,*) 'ip, nneigh',ip, nneigh

         if (nneigh == 0) then                         ! check if no neighbours
            call warn(txroutine, 'no neighbours detected', uout)
            cycle
         end if
         npart = npart + 1

         do l=0,lmax
            ql(l) = zero
            do m=-l,l                                        ! sum m over all values
               qlm = zero
               do iploc=1,nneigh
                  theta=thetann(iploc)
                  phi=phinn(iploc)
                  qlm =  qlm + CCLM(l,m,theta,phi,0)         ! use CCLM subroutine for calculation Qlm
               end do     ! iploc
               qlm = qlm/nneigh
               ql(l) = ql(l) + real(qlm*conjg(qlm))                ! complex conjugate
            end do ! m-values
            ql(l)=sqrt(ql(l))                                      ! end equation
            qlsum(l) = qlsum(l) + ql(l)
         end do                                              ! l-values
         call TabulationQl(lmax,ql(0),drotm(1,ip))          ! store order parameter in drotm
         write(*,*) 'ip, op',npart,force(1,ip)
      end do  ! first, iploop
      var(1:nvar)%value = qlsum(0:nvar-1)/npart
      call ScalarSample(iStage, 1, nvar, var)

      npartsum = npartsum + npart

   case (iAfterMacrostep)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      if (lsim .and. master) write(ucnf) var
   case (iAfterSimulation)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)
      write(uout,*) 'npartsum', npartsum/nframe
      write(uout,*) 'frame', nframe
      write(uout,*) 'rcutoff', rcutoff
   end select

end subroutine BondOrder

!************************************************************************
!> \page moluser moluser.F90
!! **TabulationQl**
!! *calculate an average of a scalar quantity*
!************************************************************************


subroutine TabulationQl(lmax,qlsum,op)

   implicit none

   integer(4), intent(in) :: lmax
   real(8), intent(in)    :: qlsum(0:lmax)
   real(8), intent(out)   :: op

   integer(4) :: l
   real(8) :: qlhex2d(0:lmax)

   qlhex2d(0:6) = 0.0d0
   qlhex2d(0) = 1.0d0
   qlhex2d(2) = 0.5d0
   qlhex2d(4) = 0.375d0
   qlhex2d(6) = 0.74083d0

   op = 0.0d0
   do l = 0,lmax
      op = op + abs(qlsum(l)-qlhex2d(l))
       write(*,*) 'l, qlsum(l), qlhex2d(l)',l, qlsum(l),qlhex2d(l)
   end do
   op = op/lmax
    write(*,'(a,7f12.6)') 'Tab..., op', qlsum(1:6), op
end subroutine TabulationQl

!************************************************************************
!> \page moluser moluser.F90
!! **ChainBeadCylContact**
!! *calculate probability of chain bead-cylinder contact*
!************************************************************************


   subroutine ChainBeadCylContact(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainBeadCylContact'
   character(80), parameter :: txheading ='probability of chain bead-cylinder contacts'
   integer(4)   , parameter :: ntype = 10
   integer(4),                 save :: nbin
   integer(4),                 save :: nvar
   real(8),                    save :: rcontact
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:)
   integer(4) :: iseg, ic, ict, ip, ivar, ibin
   real(8) ::  r2
   character(6) :: str1

   namelist /nmlCBCC/ rcontact

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      rewind(uin)
      read(uin,nmlCBCC)

      if (nct > 1) call Stop(txroutine, 'nct > 1', uout)

   case (iWriteInput)

      nbin = npct(1)
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

! ... set nvar as well as allocate memory

      nvar = nc
      allocate(var(nvar), ipnt(nc, ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do ic = 1, nc                         ! loop over all chains
         ict = ictcn(ic)
         ivar = ivar+1
         ipnt(ic,1) = ivar
         write(str1,'(i6)') ic
         var(ivar)%label = 'c.'//trim(adjustl(str1))
         var(ivar)%min = 1-0.5
         var(ivar)%max = npct(ict)+0.5
         var(ivar)%nbin = nbin
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc
         ict = ictcn(ic)
         do iseg = 1, npct(ict)
            ip = ipnsegcn(iseg,ic)
               r2 = ro(1,ip)**2 + ro(2,ip)**2
               if (r2 < rcontact**2) then
                  ivar = ipnt(ic,1)
                  ibin = iseg-1
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if
         end do
      end do

   case (iAfterMacrostep)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

    write(*,*) '1'
      call DistFuncSample(iStage, nvar, var)
    write(*,*) '2'
      call WriteHead(2, txheading, uout)
    write(*,*) '3'
      write(uout,'(a,t35,f8.2)') 'upper contact distance         = ', rcontact
      call DistFuncHead(nvar, var, uout)
    write(*,*) '4'
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
    write(*,*) '5'

     deallocate(var, ipnt)
    write(*,*) '6'

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainBeadCylContact

!************************************************************************
!> \page moluser moluser.F90
!! **Min**
!! *calculate minimum energy configuration of a bcc "sheet" or "tube"*
!************************************************************************


subroutine Min(iStage)

   use MolModule
   implicit none

   integer(4) :: iStage, sgn
   real(8) :: delta, scfac, U1, U2

   select case (iStage)
   case (iReadInput)

   case (iWriteInput)

   case (iBeforeSimulation)

      U1 = u%tot
      scfac=1.0d0
      delta =0.01d0
      sgn = 1

      do
         scfac= scfac + sgn*delta
         Call setsheet(1,scfac,1.0d0,.false.)
         Call UTotal(iStage)
         U2 = u%tot
         if (U2 > U1) then
            sgn = -sgn
            scfac= scfac + sgn*delta
            delta = 0.5d0*delta
         else
            U1=U2
         end if
         If (delta.lt.1.0d-7) exit
      end do

   case (iBeforeMacrostep)

   case (iSimulationStep)

   case (iAfterMacrostep)

   case (iAfterSimulation)

   end select

end subroutine Min

!************************************************************************
!> \page moluser moluser.F90
!! **S1**
!! *calculate orientational order parameter S1 in all dimensions*
!************************************************************************

! ... with respect to particle x-axis.

subroutine S1(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), parameter :: mnvar = 3
   integer(4), intent(in) :: iStage

   character(60),    save :: head = 'S1 order parameter'
   integer(4),       save :: nvar
   type(scalar_var), save :: var(mnvar)
   real(8)                :: tmp(mnvar)
   integer                :: ip, j

   if (ltrace) call WriteTrace(2,'S1', iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      if (nvar > mnvar) call Stop('ScalarDemo', 'nvar > mnvar', uout)
      var(1)%label = 'S1 order parameter (x)         = '
      var(2)%label = 'S1 order parameter (y)         = '
      var(3)%label = 'S1 order parameter (z)         = '
      nvar=mnvar
      var(1:3)%norm = One/np

   case (iWriteInput)

      call WriteHead(2, head, uout)

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      tmp(1:mnvar) = Zero
      do ip=1, np
         do j=1, mnvar
            tmp(j) = tmp(j) + ori(j,1,ip)
         end do
      end do
       do j=1, mnvar
          var(j)%value = tmp(j)
       end do

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         if (lsim .and. master) write(ucnf) var(1:nvar) ! possition ?
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   end select

end subroutine S1

!************************************************************************
!> \page moluser moluser.F90
!! **S2**
!! *calculate orientational order parameter S2 in all dimensions*
!************************************************************************

! ... with respect to particle z-axis.

subroutine S2(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), parameter :: mnvar = 3
   integer(4), intent(in) :: iStage

   character(60),    save :: head = 'S2 order parameter'
   integer(4),       save :: nvar
   type(scalar_var), save :: var(mnvar)
   real(8)                :: tmp(mnvar)
   integer                :: ip, j

   if (ltrace) call WriteTrace(2,'S2', iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      if (nvar > mnvar) call Stop('ScalarDemo', 'nvar > mnvar', uout)
      var(1)%label = 'S2 order parameter (x)         = '
      var(2)%label = 'S2 order parameter (y)         = '
      var(3)%label = 'S2 order parameter (z)         = '
      var(1:3)%norm = One/np
      nvar=mnvar

   case (iWriteInput)

      call WriteHead(2, head, uout)

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      tmp(1:nvar) = Zero
      do ip=1, np
         do j=1, nvar
            tmp(j) = tmp(j) + 1.5d0*ori(j,3,ip)**2-0.5d0
         end do
      end do
       do j=1, nvar
          var(j)%value = tmp(j)
       end do

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         if (lsim .and. master) write(ucnf) var(1:nvar) ! possition ?
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   end select

end subroutine S2

!************************************************************************
!> \page moluser moluser.F90
!! **C2**
!! *calculate C2 nematic order parameter of projection on xy-plane*
!************************************************************************


subroutine C2(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), parameter :: mnvar = 1
   integer(4), intent(in) :: iStage

   character(60),    save :: head = 'nematic order parameter'
   integer(4),       save :: nvar
   type(scalar_var), save :: var(mnvar)
   real(8)                :: tmp
   integer                :: ip, jp

   if (ltrace) call WriteTrace(2,'Q', iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      if (nvar > mnvar) call Stop('ScalarDemo', 'nvar > mnvar', uout)
      var(1)%label = 'C2                     = '
      var(1)%norm = One
      nvar=mnvar

   case (iWriteInput)

      call WriteHead(2, head, uout)

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var(1)%nsamp2 = var(1)%nsamp2+Half*np*(np-1)
      var(1)%value = Zero
      do ip=1, np
         do jp=ip+1, np
            tmp=sqrt(ori(1,3,ip)**2+ori(2,3,ip)**2)*sqrt(ori(1,3,jp)**2+ori(2,3,jp)**2)
            tmp = (ori(1,3,ip)*ori(1,3,jp)+ori(2,3,ip)*ori(2,3,jp))/tmp
            tmp = 2.0e0*tmp**2-1.0e0
            var(1)%value = var(1)%value + tmp
         end do
      end do
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         if (lsim .and. master) write(ucnf) var(1:nvar) ! possition ?
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   end select

end subroutine C2

!************************************************************************
!> \page moluser moluser.F90
!! **Q**
!! *calculate nematic order parameter of projection on xy-plane with x*
!************************************************************************

! as reference direction.

subroutine Q(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), parameter :: mnvar = 2
   integer(4), intent(in) :: iStage

   character(60),    save :: head = 'nematic order parameter'
   integer(4),       save :: nvar
   type(scalar_var), save :: var(mnvar)
   real(8)                :: tmp(mnvar), x, y
   integer                :: ip, j

   if (ltrace) call WriteTrace(2,'Q', iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      if (nvar > mnvar) call Stop('ScalarDemo', 'nvar > mnvar', uout)
      var(1)%label = '<sin(2a)>              = '
      var(2)%label = '<cos(2a)>              = '
      var(1:2)%norm = One/np
      nvar=mnvar

   case (iWriteInput)

      call WriteHead(2, head, uout)

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      tmp(1:nvar) = Zero
      do ip=1, np
         x = ori(1,3,ip)/sqrt(ori(1,3,ip)**2+ori(2,3,ip)**2)
         y = ori(2,3,ip)/sqrt(ori(1,3,ip)**2+ori(2,3,ip)**2)
         tmp(1) = tmp(1) + 2.0d0*x*y
         tmp(2) = tmp(2) + x**2 - y**2
       end do
       do j=1, nvar
          var(j)%value = tmp(j)
       end do

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         if (lsim .and. master) write(ucnf) var(1:nvar) ! possition ?
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call ScalarNorm(iStage, 1, nvar, var, 1)
         call WriteHead(2, head, uout)
         call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5)', uout)
      end if

   end select

end subroutine Q

subroutine CylDistFunc(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), intent(in) :: iStage

   character(60), parameter :: txheading = 'cylindrical 2d distribution functions'
   integer(4),    parameter :: ntype = 3
   type(static2D_var), save :: vtype(ntype)
   integer(4),    save :: nvar
   integer(4), allocatable,     save :: ipnt(:,:)
   type(df2d_var), allocatable, save :: var(:)

   character(60), parameter :: txheadingC = '2D nematic order parameter'
   integer(4),    parameter :: nvarC = 1
   type(scalar_var), save :: varC(nvarC)

   character(60), parameter :: txheadingS = '3D nematic order parameter'
   integer(4),    parameter :: nvarS = 1
   type(scalar_var), save :: varS(nvarS)

   integer(4)    :: ip, jp, ipt, itype, ibin1, ibin2, ivar
   real(8)       :: dz,dx,dy,drho,S2,C2,tmp

   namelist /nmlCylDistFunc/ vtype

   if(ltrace) call WriteTrace(1,'CylDistFunc', iStage)

   if (ltime) call CpuAdd('start', 'CylDistFunc', 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype(1:ntype)%label = [ 'cyldf', 'C2df ', 'S2df ' ]
      do itype = 1, ntype
         vtype(itype)%l =.true.
         vtype(itype)%min(1:2) = [ Zero, Zero ]
         vtype(itype)%max(1:2) = [ 50.0d0, 50.0d0 ]
         vtype(itype)%nbin(1:2) = [ 10, 10 ]
      end do

      rewind(uin)
      read(uin,nmlCylDistFunc)

! ... check condition

      if(vtype(1)%nbin(1) > mnbin_df2d) call Stop('CylDistFunc', 'vtype(1)%nbin(1) > mnbin_df2d', uout)
      if(vtype(1)%nbin(2) > mnbin_df2d) call Stop('CylDistFunc', 'vtype(1)%nbin(2) > mnbin_df2d', uout)

   case (iWriteInput)

! ... set nvar and allocate memory

      nvar = ntype*npt
      if (.not.allocated(var)) then
         allocate(var(nvar), ipnt(ntype,nvar))
         ipnt = 0
      end if

! ... set ipnt, label, vlow, vupp, and nbin

      ivar = 0
      do itype = 1, ntype
         do ipt = 1, npt
            ivar = ivar + 1
            ipnt(ipt,itype) = ivar
            var(ivar)%label = vtype(itype)%label
            var(ivar)%min(1:2) = vtype(1)%min(1:2)  ! note itype = 1
            var(ivar)%max(1:2) = vtype(1)%max(1:2)  ! note itype = 1
            var(ivar)%nbin(1:2)= vtype(1)%nbin(1:2) ! note itype = 1
         end do
      end do

      call DistFunc2DSample(iStage, nvar, var)

      varC(1)%label = 'C2                     = '
      varC(1)%norm = One
      varS(1)%label = 'S2                     = '
      varS(1)%norm = One

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)
      call ScalarSample(iStage, 1, nvarC, varC)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) varC(1:nvar)
      call ScalarSample(iStage, 1, nvarS, varS)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) varS(1:nvar)

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvarC, varC)
      call ScalarSample(iStage, 1, nvarS, varS)

   case (iSimulationStep)

      varC(1)%nsamp2 = varC(1)%nsamp2+Half*np*(np-1)
      varS(1)%nsamp2 = varS(1)%nsamp2+Half*np*(np-1)
      var(1:nvar)%nsamp2 = var(1:nvar)%nsamp2+np
      varC(1)%value = Zero
      varS(1)%value = Zero

      itype = 1
      do ip = 1, np
         do jp=ip+1, np
            if (.not.vtype(itype)%l) cycle
            if (ip==jp) cycle
            ipt = iptpn(jp)
            ivar = ipnt(itype,ipt)
            dx = ro(1,jp)-ro(1,ip)
            dy = ro(2,jp)-ro(2,ip)
            dz = ro(3,jp)-ro(3,ip)
            call PBC(dx,dy,dz)
            drho = sqrt(dx**2+dy**2)
            dz = abs(dz)
            tmp=sqrt(ori(1,3,ip)**2+ori(2,3,ip)**2)*sqrt(ori(1,3,jp)**2+ori(2,3,jp)**2)
            C2 = (ori(1,3,ip)*ori(1,3,jp)+ori(2,3,ip)*ori(2,3,jp))/tmp
            C2 = 2.0e0*C2**2-1.0e0
            S2 = ori(1,3,ip)*ori(1,3,jp)+ori(2,3,ip)*ori(2,3,jp)+ori(3,3,ip)*ori(3,3,jp)
            S2 = 1.5e0*S2**2-0.5e0
            ibin1 = max(-1,min(floor(var(1)%bini(1)*(dz-var(1)%min(1))),var(1)%nbin(1)))
            ibin2 = max(-1,min(floor(var(1)%bini(2)*(drho-var(1)%min(2))),var(1)%nbin(2)))
            var(1)%avs2(ibin1,ibin2) = var(1)%avs2(ibin1,ibin2)+One
            ibin1 = max(-1,min(floor(var(2)%bini(1)*(dz-var(2)%min(1))),var(2)%nbin(1)))
            ibin2 = max(-1,min(floor(var(2)%bini(2)*(drho-var(2)%min(2))),var(2)%nbin(2)))
            var(2)%avs2(ibin1,ibin2) = var(2)%avs2(ibin1,ibin2)+C2
            varC(1)%value = varC(1)%value + C2
            var(3)%avs2(ibin1,ibin2) = var(3)%avs2(ibin1,ibin2)+S2
            varS(1)%value = varS(1)%value + S2
         end do
      end do

      call ScalarSample(iStage, 1, nvarC, varC)

      call ScalarSample(iStage, 1, nvarS, varS)

   case (iAfterMacrostep)

! ... normalization

      if(master) call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var(1:nvar)

      if (master) then
         call ScalarSample(iStage, 1, nvarC, varC)
         call ScalarNorm(iStage, 1, nvarC, varC, 1)
         if (lsim .and. master) write(ucnf) varC(1:nvar) ! possition ?
         call WriteHead(2, txheadingC, uout)
         call ScalarWrite(iStage, 1, nvarC, varC, 0, '(a,t35,2f15.5)', uout)
      end if

      if (master) then
         call ScalarSample(iStage, 1, nvarS, varS)
         call ScalarNorm(iStage, 1, nvarS, varS, 1)
         if (lsim .and. master) write(ucnf) varS(1:nvar) ! possition ?
         call WriteHead(2, txheadingS, uout)
         call ScalarWrite(iStage, 1, nvarS, varS, 0, '(a,t35,2f15.5)', uout)
      end if

   case (iAfterSimulation)

      if(master) then
         call DistFunc2DSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,2i5,5x,a,i5,2x,a)') 'number of grid points          = ', var(1)%nbin, '(',mnbin_df2d,')'
         call DistFunc2DHead(nvar, var, uout)

!         inquire (uuser, opened=connected)
!         if (.not.connected) call FileOpen(uuser, fuser, 'form/noread')
!         write(uuser,'(a)') '2d CylDistFunc'
!         write(uuser,'(a)')  '1'
!         write(uuser,'(a)') txtype(1)
!         write(uuser,'(a)')  '1'
!         do ibin1=0, var(1)%nbin(1)-1
!            write(uuser,'')
!            do ibin2=0, var(1)%nbin(2)-1
!              tmp = Boxlen(1)*Boxlen(2)*Boxlen(3)/(TwoPi*var(1)%bin(1)*(Two*(var(1)%min(2)+ibin2*var(1)%bin(2))*var(1)%bin(2) + var(1)%bin(2)**2))/(Half*(np-1))
!              norm = InvFlt(var(1)%avs1(ibin1,ibin2))
!              ratio = var(2)%avs1(ibin1,ibin2)*norm
!              dratio = var(2)%avsd(ibin1,ibin2)*norm
!              ratio2 = var(3)%avs1(ibin1,ibin2)*norm
!              dratio2 = var(3)%avsd(ibin1,ibin2)*norm
!              write(uuser,'(8e13.4)') var(1)%min(1)+(ibin1+0.5)*var(1)%bin(1), var(1)%min(2)+(ibin2+0.5)*var(1)%bin(2), var(1)%avs1(ibin1,ibin2)*tmp, var(1)%avsd(ibin1,ibin2)*tmp, ratio, dratio, rati
!            end do
!            write(uuser,'(a1)') ' '
!         end do

! ... normalize

         do ibin1=0, var(1)%nbin(1)-1
            do ibin2=0, var(1)%nbin(2)-1
              tmp = Boxlen(1)*Boxlen(2)*Boxlen(3)/(TwoPi*var(1)%bin(1)*(Two*(var(1)%min(2)+ibin2*var(1)%bin(2))*var(1)%bin(2) + var(1)%bin(2)**2))/(Half*(np-1))
              var(1)%avs1(ibin1,ibin2)=var(1)%avs1(ibin1,ibin2)*tmp
              var(1)%avsd(ibin1,ibin2)=var(1)%avsd(ibin1,ibin2)*tmp
            end do
         end do

         call DistFunc2DShow(1, txheading, nvar, var, uout)
         call DistFunc2DList(1, txheading, nvar, var, ulist)

         call ScalarSample(iStage, 1, nvarC, varC)
         call ScalarNorm(iStage, 1, nvarC, varC, 1)
         call WriteHead(2, txheadingC, uout)
         call ScalarWrite(iStage, 1, nvarC, varC, 0, '(a,t35,2f15.5)', uout)

         call ScalarSample(iStage, 1, nvarS, varS)
         call ScalarNorm(iStage, 1, nvarS, varS, 1)
         call WriteHead(2, txheadingS, uout)
         call ScalarWrite(iStage, 1, nvarS, varS, 0, '(a,t35,2f15.5)', uout)

      end if

   end select

    if (ltime) call CpuAdd('stop', 'CylDistFunc', 1, uout)

end subroutine CylDistFunc

!************************************************************************
!> \page moluser moluser.F90
!! **SFPBC2D**
!! *documentation_missing*
!************************************************************************

subroutine SFPBC2D(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), intent(in) :: iStage

   character(60), parameter :: txheading = '2D structure factor (cylindrical symmetry)'
   integer(4),    parameter :: ntype = 6
   character(14), parameter :: txtype(ntype) = [ "SI 10,1       ", "SI 10,1 points", "SI 10,1 shape ", "SI 11,1       ", "SI 11,1 points", "SI 11,1 shape " ]
   integer(4)   , save :: dirlow(ntype), dirupp(ntype)
   integer(4)   , save :: nbin(2), nvar
   type(df2d_var), allocatable, save :: var(:), Withff(:), ff(:)
   integer(4),     allocatable, save :: ipnt(:,:)
   real(8),        allocatable :: sfpar(:,:,:), sfparsd(:,:)
   real(8)    :: norm, kx, ky, kz, kk
   integer(4) :: ip, ipt, jpt, ivar, iptjpt, ibin, jbin, itype, m, dupp, dlow
   real(8)    :: sffac(ntype),  uu, fff, ffstore(-1:mnbin_df2d,-1:mnbin_df2d,1:13,1:mnpt)
   complex(8) :: eikr(-1:mnbin_df2d,-1:mnbin_df2d,1:4,1:mnpt), eikrwff(-1:mnbin_df2d,-1:mnbin_df2d,1:4,1:mnpt), eikr1(1:4), eikrtemp, eikrz
   real(8), save :: scrad, scasp

   namelist /nmlSF2D/ nbin, scrad, scasp

   if(ltrace) call WriteTrace(1,'SFPBC2D', iStage)

   if (ltime) call CpuAdd('start', 'SFPBC2D', 1, uout)

   select case (iStage)
   case (iReadInput)

      nbin(1:2) = 40
      scrad = one
      scasp = one
      rewind(uin)
      read(uin,nmlSF2D)

      dirlow(1:ntype) = [ 1, 1, 1, 3, 3, 3 ]
      dirupp(1:ntype) = [ 2, 2, 2, 4, 4, 4 ]

      if(boxlen(1) /= boxlen(2)) call Stop('SFPBC2D', 'boxlen(1) /= boxlen(2)', uout)
      if(npt /= 1) call Stop('SFPBC2D', 'npt /= 1', uout)
      if(nbin(1) > mnbin_df2d) call Stop('SFPBC2D', 'nbin(1) > mnbin_df2d', uout)
      if(nbin(2) > mnbin_df2d) call Stop('SFPBC2D', 'nbin(2) > mnbin_df2d', uout)

      sffac(1:ntype) = [ One, One, One, sqrt(Two), sqrt(Two), sqrt(Two) ]

   case (iWriteInput)

! ... set nvar and allocate memory

      nvar = ntype*npt**2
      if (.not.allocated(var)) then
         allocate(var(nvar), Withff(nvar), ff(nvar))
      end if
      if (.not.allocated(ipnt)) then
         allocate(ipnt(ntype,nptpt), sfpar(2*mnbin_df2d,2*mnbin_df,nptpt), sfparsd(2*mnbin_df2d,nptpt))
         ipnt = 0
         sfpar = 0.0E+00
         sfparsd = 0.0E+00
      end if

! ... set ipnt, label, vlow, vupp, and nbin

      nvar = 0
      do itype = 1, ntype
         do ipt  = 1, npt
            do jpt = 1, npt
               iptjpt = iptpt(ipt,jpt)
               nvar = nvar+1
               ipnt(itype,iptjpt) = nvar
               var(nvar)%label = trim(txtype(itype))//' '//txpt(ipt)
               var(nvar)%min(1) = Half*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%min(2) = Half*TwoPiBoxi(3)
               var(nvar)%max(1) = (nbin(1)+Half)*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%max(2) = (nbin(2)+Half)*TwoPiBoxi(3)
               var(nvar)%nbin(1:2)= nbin(1:2)
            end do
         end do
      end do

      call DistFunc2DSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)

   case (iSimulationStep)

      var(1:nvar)%nsamp2 = var(1:nvar)%nsamp2+1

      eikr(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)
      eikrwff(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)
      ffstore(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)

      do ip = 1, np
         ipt = iptpn(ip)
         kx = TwoPiBoxi(1)*ro(1,ip)
         ky = TwoPiBoxi(2)*ro(2,ip)
         kz = TwoPiBoxi(3)*ro(3,ip)
         eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
         eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
         eikr1(3) = cmplx(cos(kx+ky),sin(kx+ky))        ! ( 1 1 0) direction
         eikr1(4) = cmplx(cos(kx-ky),sin(kx-ky))        ! ( 1-1 0) direction
         eikrz = cmplx(cos(kz),sin(kz))                 ! ( 0 0 1) direction
         do m=1,  dirupp(ntype)
            do ibin = -1, nbin(2)-1
               eikrtemp = eikrz**(ibin+1)
               do jbin = -1, nbin(1)-1
                  eikr(jbin,ibin,m,ipt) = eikr(jbin,ibin,m,ipt)+eikrtemp
                  if (m == 1) then
                     kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                     uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                  else if (m == 2) then
                     kk = sqrt((TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                     uu = (TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                  else if (m == 3) then
                     kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                     uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) + TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                  else if (m == 4) then
                     kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                     uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) - TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                  else
                     call Stop('SFPBC2D', 'invalid direction (m)', uout)
                  end if
                  uu = radellipsoid*scrad*kk*sqrt((1.0-uu**2)+(aellipsoid*scasp)**2*uu**2)
                  fff = 3.0*(sin(uu)-uu*cos(uu))/uu**3
                  eikrwff(jbin,ibin,m,ipt) = eikrwff(jbin,ibin,m,ipt)+eikrtemp*fff
                  ffstore(jbin,ibin,m,ipt) = ffstore(jbin,ibin,m,ipt)+ fff**2
                  eikrtemp = eikrtemp*eikr1(m)
               end do
            end do
         end do
      end do

      do itype = 1, ntype
         dlow=dirlow(itype)
         dupp=dirupp(itype)
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(itype,iptpt(ipt,jpt))
               do ibin = -1, nbin(2)-1
                  do jbin = -1, nbin(1)-1
                     if (mod(itype,3) == 1) then
                        var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ &
                             sum(real(eikrwff(jbin,ibin,dlow:dupp,ipt))*real(eikrwff(jbin,ibin,dlow:dupp,jpt)) &
                          +aimag(eikrwff(jbin,ibin,dlow:dupp,ipt))*aimag(eikrwff(jbin,ibin,dlow:dupp,jpt)))
                     else if (mod(itype,3) == 2) then
                        var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ &
                             sum(real(eikr(jbin,ibin,dlow:dupp,ipt))*real(eikr(jbin,ibin,dlow:dupp,jpt)) &
                             +aimag(eikr(jbin,ibin,dlow:dupp,ipt))*aimag(eikr(jbin,ibin,dlow:dupp,jpt)))
                     else if (mod(itype,3) == 0) then
                        var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ sum(ffstore(jbin,ibin,dlow:dupp,ipt))
                     end if
                  end do
               end do
            end do
         end do
      end do

   case (iAfterMacrostep)

! ... normalization

      do itype = 1, ntype
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(itype,iptpt(ipt,jpt))
               norm = One/((dirupp(itype)-dirlow(itype)+1)*sqrt(real(nppt(ipt))*real(nppt(jpt))))
               var(ivar)%avs2(-1:nbin(1),-1:nbin(2)) = var(ivar)%avs2(-1:nbin(1),-1:nbin(2))*norm
            end do
         end do
      end do

      call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var(1:nvar)

   case (iAfterSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,2i5,5x,a,i5,2x,a)') 'number of grid points          = ', nbin, '(',mnbin_df2d,')'
      call DistFunc2DHead(nvar, var, uout)
      call DistFunc2DShow(1, txheading, nvar, var, uout)
      call DistFunc2DList(1, txheading, nvar, var, ulist)

!      inquire (uuser, opened=connected)
!      if (.not.connected) call FileOpen(uuser, fuser, 'form/noread')
!      write(uuser,*) 'from SFPBC2D: SI_1 data'
!      write(uuser,'(a1)') ' '
!      do ibin=-nbin(2), nbin(2)
!         do jbin=-nbin(1), nbin(1)
!            if ((ibin == 0) .and. (jbin == 0)) then
!               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),Zero,Zero,Zero,Zero,Zero,Zero
!            else
!               write(uuser,'(8e13.4)') jbin*var(1)%bin(1), ibin*var(1)%bin(2),&
!                var(1)%avs1(abs(jbin)-1,abs(ibin)-1), var(1)%avsd(abs(jbin)-1,abs(ibin)-1), &
!                var(2)%avs1(abs(jbin)-1,abs(ibin)-1), var(2)%avsd(abs(jbin)-1,abs(ibin)-1), &
!                var(3)%avs1(abs(jbin)-1,abs(ibin)-1), var(3)%avsd(abs(jbin)-1,abs(ibin)-1)
!            End if
!         end do
!         write(uuser,'(a1)') ' '
!      end do
!
!      write(uuser,*) 'from SFPBC2D: SI_2 data'
!      write(uuser,'(a1)') ' '
!      do ibin=-nbin(2), nbin(2)
!         do jbin=-nbin(1), nbin(1)
!            if ((ibin == 0) .and. (jbin == 0)) then
!               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),Zero,Zero,Zero,Zero,Zero,Zero
!            else
!               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),&
!                var(4)%avs1(abs(jbin)-1,abs(ibin)-1), var(4)%avsd(abs(jbin)-1,abs(ibin)-1), &
!                var(5)%avs1(abs(jbin)-1,abs(ibin)-1), var(5)%avsd(abs(jbin)-1,abs(ibin)-1), &
!                var(6)%avs1(abs(jbin)-1,abs(ibin)-1), var(6)%avsd(abs(jbin)-1,abs(ibin)-1)
!            End if
!         end do
!         write(uuser,'(a1)') ' '
!      end do

   end select

   if (ltime) call CpuAdd('stop', 'SFPBC2D', 1, uout)

end subroutine SFPBC2D

!************************************************************************
!> \page moluser moluser.F90
!! **SFPBC2DNoAv**
!! *Does not average over x and y, compatible with boxlength in x \= y*
!************************************************************************


subroutine SFPBC2DNoAv(iStage)

   use MolModule
   use StatisticsModule
   implicit none

   integer(4), parameter :: mntype = 6
   integer(4), parameter :: mnvar = mntype*mnpt

   integer(4), intent(in) :: iStage
   !character(60), save :: head = '2D structure factor'
   character(18),  save :: txtype(mntype) = [ "SI '10,1'         ", "SI (points) '10,1'", "SI (shape) '10,1' ", "SI '11,1'         ", "SI (points) '11,1'", "SI (shape) '11,1' " ]
   integer(4)   , save :: dirlow(mntype), dirupp(mntype)
   integer(4)   , save :: nbin(2), ntype, nvar
   integer(4), allocatable, save :: ipnt(:,:)
   real(8)    :: norm, kx, ky, kz, kk
   real(8), allocatable :: sfpar(:,:,:), sfparsd(:,:)
   type(df2d_var),  save :: var(mnvar)
   integer(4) :: ip, ipt, jpt, ivar, iptjpt, ibin, jbin, itype, m, dupp, dlow
   real(8)       :: sffac(mntype),  uu, fff, ffstore(0:mnbin_df2d,0:mnbin_df2d,1:13,1:mnpt)
   complex(8) :: eikr(0:mnbin_df2d,0:mnbin_df2d,1:4,1:mnpt), eikrwff(0:mnbin_df2d,0:mnbin_df2d,1:4,1:mnpt), eikr1(1:4), eikrtemp, eikrz
   integer(4), save :: isf
   logical :: connected

   namelist /nmlSF2D/ nbin, isf

   if(ltrace) call WriteTrace(1,'SFPBC2DNoAv', iStage)

   if (ltime) call CpuAdd('start', 'SFPBC2DNoAv', 1, uout)

   select case (iStage)
   case (iReadInput)

      nbin(1:2) = 40
      ntype=mntype
      isf = 1

      rewind(uin)
      read(uin,nmlSF2D)

      dirlow(1:ntype) = [ 1, 1, 1, 2, 2, 2 ]
      dirupp(1:ntype) = [ 1, 1, 1, 2, 2, 2 ]

      if(npt /= 1) call Stop('SFPBC2DNoAv', 'npt /= 1', uout)
      if(nbin(1) > mnbin_df2d) call Stop('SFPBC2DNoAv', 'nbin(1) > mnbin_df2d', uout)
      if(nbin(2) > mnbin_df2d) call Stop('SFPBC2DNoAv', 'nbin(2) > mnbin_df2d', uout)
      sffac(1:ntype) = [ One, One, One, sqrt(Two), sqrt(Two), sqrt(Two) ]

   case (iWriteInput)

      allocate(ipnt(mntype,nptpt), sfpar(2*mnbin_df2d,2*mnbin_df,nptpt), sfparsd(2*mnbin_df2d,nptpt))
      ipnt = 0
      sfpar = 0.0E+00
      sfparsd = 0.0E+00

! ... set nvar, ipnt, label, vlow, vupp, and nbin

      nvar = 0
      do itype = 1, ntype/2
         do ipt  = 1, npt
            do jpt = 1, npt
               iptjpt = iptpt(ipt,jpt)
               nvar = nvar+1
               ipnt(itype,iptjpt) = nvar
               var(nvar)%label = trim(txtype(itype))//' '//txpt(ipt)
               var(nvar)%min(1) = Half*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%min(2) = Half*TwoPiBoxi(3)
               var(nvar)%max(1) = (nbin(1)+Half)*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%max(2) = (nbin(2)+Half)*TwoPiBoxi(3)
               var(nvar)%nbin(1:2)= nbin(1:2)
            end do
         end do
      end do

      do itype = ntype/2+1, ntype
         do ipt  = 1, npt
            do jpt = 1, npt
               iptjpt = iptpt(ipt,jpt)
               nvar = nvar+1
               ipnt(itype,ipt) = nvar
               var(nvar)%label = trim(txtype(itype))//' '//txpt(ipt)
               var(nvar)%min(1) = Half*TwoPiBoxi(2)*sffac(itype)
               var(nvar)%min(2) = Half*TwoPiBoxi(3)
               var(nvar)%max(1) = (nbin(1)+Half)*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%max(2) = (nbin(2)+Half)*TwoPiBoxi(3)
               var(nvar)%nbin(1:2)= nbin(1:2)
            end do
         end do
      end do
      if(nvar > mnvar) call Stop('SFPBC2D', 'nvar > mnvar', uout)

      call DistFunc2DSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master .and. txstart == 'continue') read(ucnf) var(1:nvar)

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)

   case (iSimulationStep)

      If (mod(istep,isf)==0) then

         var(1:nvar)%nsamp2 = var(1:nvar)%nsamp2+1
         eikr(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)
         eikrwff(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)
         ffstore(-1:nbin(1),-1:nbin(2),1:dirupp(ntype),1:npt) = cmplx(Zero,Zero)

         do ip = 1, np
            ipt = iptpn(ip)
            kx = TwoPiBoxi(1)*ro(1,ip)
            ky = TwoPiBoxi(2)*ro(2,ip)
            kz = TwoPiBoxi(3)*ro(3,ip)
            eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
            eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
            eikrz = cmplx(cos(kz),sin(kz)) ! ( 0 0 1) direction
            do m=1,  dirupp(ntype)
               do ibin = -1, nbin(2)-1
                  eikrtemp = eikrz**(ibin+1)
                  do jbin = -1, nbin(1)-1
                     eikr(jbin,ibin,m,ipt) = eikr(jbin,ibin,m,ipt)+eikrtemp
                     if (m == 1) then
                        kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                        uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                     else if (m == 2) then
                        kk = sqrt((TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                        uu = (TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                     else if (m == 3) then
                        kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                        uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) + TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                     else if (m == 4) then
                        kk = sqrt((TwoPiBoxi(1)*(jbin+1))**2 + (TwoPiBoxi(2)*(jbin+1))**2 + (TwoPiBoxi(3)*(ibin+1))**2)
                        uu = (TwoPiBoxi(1)*(jbin+1)*ori(1,3,ip) - TwoPiBoxi(2)*(jbin+1)*ori(2,3,ip) + TwoPiBoxi(3)*(ibin+1)*ori(3,3,ip))/kk
                     else
                        call Stop('SFPBC2D', 'invalid direction (m)', uout)
                     end if
                     uu = radellipsoid*kk*sqrt((1.0-uu**2)+aellipsoid**2*uu**2)
                     fff = 3.0*(sin(uu)-uu*cos(uu))/uu**3
                     eikrwff(jbin,ibin,m,ipt) = eikrwff(jbin,ibin,m,ipt)+eikrtemp*fff
                     ffstore(jbin,ibin,m,ipt) = ffstore(jbin,ibin,m,ipt)+ fff**2
                     eikrtemp = eikrtemp*eikr1(m)
                  end do
               end do
            end do
         end do

         do itype = 1, ntype
            dlow=dirlow(itype)
            dupp=dirupp(itype)
            do ipt = 1, npt
               do jpt = ipt, npt
                  ivar = ipnt(itype,iptpt(ipt,jpt))
                  do ibin = -1, nbin(2)-1
                     do jbin = -1, nbin(1)-1
                        if (mod(itype,3) == 1) then
                           var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ &
                                sum(real(eikrwff(jbin,ibin,dlow:dupp,ipt))*real(eikrwff(jbin,ibin,dlow:dupp,jpt)) &
                             +aimag(eikrwff(jbin,ibin,dlow:dupp,ipt))*aimag(eikrwff(jbin,ibin,dlow:dupp,jpt)))
                        else if (mod(itype,3) == 2) then
                           var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ &
                                sum(real(eikr(jbin,ibin,dlow:dupp,ipt))*real(eikr(jbin,ibin,dlow:dupp,jpt)) &
                                +aimag(eikr(jbin,ibin,dlow:dupp,ipt))*aimag(eikr(jbin,ibin,dlow:dupp,jpt)))
                        else if (mod(itype,3) == 0) then
                           var(ivar)%avs2(jbin,ibin) = var(ivar)%avs2(jbin,ibin)+ sum(ffstore(jbin,ibin,dlow:dupp,ipt))
                        end if
                     end do
                  end do
               end do
            end do
         end do

      endif

   case (iAfterMacrostep)

! ... normalization

      do itype = 1, ntype
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(itype,iptpt(ipt,jpt))
               norm = One/((dirupp(itype)-dirlow(itype)+1)*sqrt(real(nppt(ipt))*real(nppt(jpt))))
               var(ivar)%avs2(-1:nbin(1),-1:nbin(2)) = var(ivar)%avs2(-1:nbin(1),-1:nbin(2))*norm
            end do
         end do
      end do

      call DistFunc2DSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var(1:nvar)

   case (iAfterSimulation)

      call DistFunc2DSample(iStage, nvar, var)

      inquire (uuser, opened=connected)
      if (.not.connected) call FileOpen(uuser, fuser, 'form/noread')
!     Open(unit=666,file='SI_1',status='unknown')
      write(uuser,*) 'from SFPBC2DNoAv: SI_1 data'
      write(uuser,'(a1)') ' '
      do ibin=-nbin(2), nbin(2)
         do jbin=-nbin(1), nbin(1)
            if ((ibin == 0) .and. (jbin == 0)) then
               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),Zero,Zero,Zero,Zero,Zero,Zero
            else
               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),&
                  var(1)%avs1(abs(jbin)-1,abs(ibin)-1), var(1)%avsd(abs(jbin)-1,abs(ibin)-1), &
                  var(2)%avs1(abs(jbin)-1,abs(ibin)-1), var(2)%avsd(abs(jbin)-1,abs(ibin)-1), &
                  var(3)%avs1(abs(jbin)-1,abs(ibin)-1), var(3)%avsd(abs(jbin)-1,abs(ibin)-1)
              End if
         end do
         write(uuser,'(a1)') ' '
      end do
!     Close(unit=666,status='keep')

!     Open(unit=666,file='SI_2',status='unknown')
      write(uuser,*) 'from SFPBC2DNoAv: SI_2 data'
      write(uuser,'(a1)') ' '
      do ibin=-nbin(2), nbin(2)
         do jbin=-nbin(1), nbin(1)
            if ((ibin == 0) .and. (jbin == 0)) then
               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),Zero,Zero,Zero,Zero,Zero,Zero
            else
               write(uuser,'(8e13.4)') jbin*var(4)%bin(1), ibin*var(4)%bin(2),&
                  var(4)%avs1(abs(jbin)-1,abs(ibin)-1), var(4)%avsd(abs(jbin)-1,abs(ibin)-1), &
                  var(5)%avs1(abs(jbin)-1,abs(ibin)-1), var(5)%avsd(abs(jbin)-1,abs(ibin)-1), &
                  var(6)%avs1(abs(jbin)-1,abs(ibin)-1), var(6)%avsd(abs(jbin)-1,abs(ibin)-1)
              End if
         end do
         write(uuser,'(a1)') ' '
      end do
!     Close(unit=666,status='keep')

   end select

   if (ltime) call CpuAdd('stop', 'SFPBC2DNoAv', 1, uout)

end subroutine SFPBC2DNoAv

!************************************************************************
!> \page moluser moluser.F90
!! **PosMeanSph**
!! *calculate means of positions*
!************************************************************************


subroutine PosMeanSph(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='PosMeanSph'
   character(80), parameter :: txheading ='center-of-mass of particles of different types'
   integer(4)   , parameter :: ntype = 6
   character(6) , parameter :: txtype(ntype) = ['<xo>  ', '<yo>  ', '<zo>  ', '<xo^2>', '<yo^2>', '<zo^2>']
   type(scalar_var), allocatable, save :: var(:)
   integer(4),       allocatable, save :: ipnt(:,:)
   integer(4),                    save :: nvar

   integer(4) :: ip, ipt, ivar, itype
   real(8), allocatable, save :: rosum(:,:)

   select case (iStage)
   case (iWriteInput)

! ... set nvar and allocate memory

     nvar = npt*ntype
     allocate(var(nvar), ipnt(npt,ntype), rosum(3,npt))
     ipnt = 0
     rosum = 0.0E+00

! ... set ipnt, label, and norm

     ivar = 0
     do ipt = 1, npt
        do itype = 1, ntype
           ivar = ivar + 1
           ipnt(ipt,itype) = ivar
           var(ivar)%label = txpt(ipt)//' '//txtype(itype)
           var(ivar)%norm = One/nppt(ipt)
        end do
     end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      rosum(1:3,1:npt) = Zero
      do ip = 1, np
         ipt = iptpn(ip)
         rosum(1:3,ipt) = rosum(1:3,ipt) + ro(1:3,ip)
      end do
      do ipt = 1, npt
         do itype = 1, 6
             ivar = ipnt(ipt,itype)
             var(ivar)%value = rosum(1+mod(itype-1,3),ipt)**(1+itype/4)
         end do
      end do
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 0)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 0)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

     deallocate(var, ipnt, rosum)

   end select

end subroutine PosMeanSph

!************************************************************************
!> \page moluser moluser.F90
!! **MacroionOneSph**
!! *calculate various df for system of one macroion + ions, sph geometry*
!************************************************************************


!     type  label  quantity
!     ----  -----  --------
!     1     proj   projection of the pos of part of one type onto the surface of one particle
!     2     |idm|  magnitude of induced dipole moment
!     3     (idm)z projection of the induced dipole moment onto the z-axis

subroutine MacroionOneSph(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MacroionOneSph'
   character(80), parameter :: txheading ='various df for system of one macroion + ions, sph geometry'
   integer(4)   , parameter :: ntype = 3
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: iptmacroion
   integer(4),                 save :: iptproject
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),                    save :: radsph
   integer(4) :: itype, ivar, ibin, ip, ipt
   real(8) ::  dx, dy, dz, normi, idmm(1:3), varrms

   namelist /nmlMacroionOneSph/ vtype, iptmacroion, iptproject

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.false.
      vtype%min = [-sphrad, Zero,  -5.0d0]
      vtype%max = [+sphrad, 50.0d0, 5.0d0]
      vtype%nbin   = 100

      rewind(uin)
      read(uin,nmlMacroionOneSph)

      if (txbc /= 'sph') call Stop (txroutine,'txbc /= ''sph''',uout)                  ! spherical geometry
      if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)
      if (nppt(iptmacroion) /= 1) call Stop(txroutine, 'nnpt(iptmacroion) /= 1', uout)
      if (iptproject == iptmacroion) call Stop(txroutine, 'iptproject == iptmacroion', uout)
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

      radsph = radat(iptpn(iptmacroion))

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['projection', '|idm|     ', 'idm_z     ']
      vtype%nvar = 1

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(1,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ivar + 1
            ipnt(1,itype) = ivar
            var(ivar)%label = vtype(itype)%label
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      idmm(1:3) = Zero
      do ip = 1, np
         ipt = iptpn(ip)
         idmm(1:3) = idmm(1:3)+zat(ipt)*r(1:3,ip)
      end do

! ... sample type 1

      itype = 1
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         do ip = ipnpt(iptproject), ipnpt(iptproject) + nppt(iptproject) - 1
            dx = ro(1,ip)-ro(1,iptmacroion)
            dy = ro(2,ip)-ro(2,iptmacroion)
            dz = ro(3,ip)-ro(3,iptmacroion)
            call PBC(dx,dy,dz)
            normi = One/sqrt(dx**2+dy**2+dz**2)
            dx = dx*normi
            dy = dy*normi
            dz = dz*normi
            dz = dz*radsph
            var%nsamp2 = var%nsamp2 + 1
            ibin = max(-1,min(floor(var(ivar)%bini*(dz-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
         end do
      end if

! ... sample type 2

      itype = 2
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         varrms = sqrt(idmm(1)**2+idmm(2)**2+idmm(3)**2)
         ibin = max(-1,min(floor(var(ivar)%bini*(varrms-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 3

      itype = 3
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(idmm(3)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,i8  )') 'iptmacroion                      = ', iptmacroion
      write(uout,'(a,i8  )') 'iptproject                       = ', iptproject
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate (var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine MacroionOneSph

!************************************************************************
!> \page moluser moluser.F90
!! **MacroionTwoCyl**
!! *calculate various df for system of two macroions + ions, cyl geometry*
!************************************************************************


!     type  label        quantity
!     ----  -----        --------
!     1     z-sep        z-separation between the two macroions
!     2     z-mean       mean z-coordinate of the two macroions
!     3     net charge 1 net change of z < 0 volume
!     4     idm 1        magnitude of induced dipole moment of z < 0 volume
!     5     idm 2        magnitude of induced dipole moment of z > 0 volume
!     6     idm z-axis 1 projection of induced dipole moment on the z-axis of z < 0 volume
!     7     idm z-axis 2 projection on induced dipole moment on the z-axis of z > 0 volume
!     8     idm-dim cosa cos of angle of induced dipole moment 1 - induced dipole moment 2

subroutine MacroionTwoCyl(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MacroionTwoCyl'
   character(80), parameter :: txheading ='various df for system of two macroions + ions, cyl geometry'
   integer(4)   , parameter :: ntype = 8
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   integer(4)   ,              save :: iptmacroion
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4) :: itype, ivar, ibin, m
   real(8)    :: value, netcharge, idmv(1:3,2), idms(2), cosang

   namelist /nmlMacroionTwoCyl/ vtype, iptmacroion

   select case (iStage)
   case (iReadInput)

      vtype%l   =.false.
      vtype%min = [Zero, -cyllen/2, -Two, Zero, Zero, -One, -One, -One]
      vtype%max = [cyllen, cyllen/2, Two,  One,  One,  One,  One,  One]
      vtype%nbin = 100
!     vtype(1:2)%nbin = 20*cyllen+0.5                             ! z10 etc
!     vtype(1:2)%nbin = 0.5*cyllen+0.5                            ! 60-1 etc
      vtype%nbin = 100
      iptmacroion = 1                                             ! select particle type of macroions

      rewind(uin)
      read(uin,nmlMacroionTwoCyl)

      if (txbc /= 'cyl') call Stop (txroutine,'txbc /= ''cyl''',uout)                  ! cylindrical geometry
      if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)
      if (nppt(iptmacroion) /= 2) call Stop (txroutine, 'nppt(iptmacroion) /=2',uout)  ! 2 particles of type iptmacroion
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['z-sep   ', 'z-mean  ', 'net ch 1', 'idm 1   ', 'idm 2   ', 'idm_z 1 ', 'idm_z 2 ', 'cosang  ']
      vtype%nvar = 1

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(1,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ivar + 1
            ipnt(1,itype) = ivar
            var(ivar)%label = vtype(itype)%label
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      netcharge = sum(zat(iptpn(1:np)), 1, (r(3,1:np) < Zero))               ! net charge of z < 0 volume
      do m = 1, 3
         idmv(m,1) = sum(zat(iptpn(1:np))*r(m,1:np), 1, (r(3,1:np) < Zero))  ! idm of z < 0 volume
         idmv(m,2) = sum(zat(iptpn(1:np))*r(m,1:np), 1, (r(3,1:np) > Zero))  ! idm of z > 0 volume
      end do
      idms(1) = sqrt(idmv(1,1)**2+idmv(2,1)**2+idmv(3,1)**2)                 ! magnetude of idm ov z < volume
      idms(2) = sqrt(idmv(1,2)**2+idmv(2,2)**2+idmv(3,2)**2)                 ! magnetude of idm of z > volume
      cosang = sum(idmv(1:3,1)*idmv(1:3,2))/(idms(1)*idms(2))                ! cos(angle) of idm 1 and idm 2

! ... sample type 1

      itype = 1
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         value = abs(ro(3,ipnpt(iptmacroion))-ro(3,ipnpt(iptmacroion)+1))
         ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 2

      itype = 2
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         value = sum(ro(3,ipnpt(iptmacroion):ipnpt(iptmacroion)+nppt(iptmacroion)-1))/nppt(iptmacroion)
         ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 3

      itype = 3
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(netcharge-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 4

      itype = 4
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(idms(1)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 5

      itype = 5
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(idms(2)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 6

      itype = 6
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(idmv(3,1)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 7

      itype = 7
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(idmv(3,2)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

! ... sample type 8

      itype = 8
      if(vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(cosang-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate (var, ipnt)

   end select

end subroutine MacroionTwoCyl

!************************************************************************
!> \page moluser moluser.F90
!! **MeanElFieldZCyl**
!! *calculate mean electrostatic field in the z-direction*
!************************************************************************


!     system: cylinder with particles 1 and 2 fixed, z(1) < z(2) is needed

subroutine MeanElFieldZCyl(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MeanElFieldZCyl'
   character(80), parameter :: txheading ='mean electrostatic field in the z-direction and mean z-force'
   integer(4)   , parameter :: nlayer = 3
   integer(4)   , parameter :: npoint = 21
   integer(4)   , parameter :: nvar = 3*nlayer*npoint
   integer(4)   , parameter :: nvar_2 = 2
   type(scalar_var),  allocatable, save :: var(:)             ! sampling of electrostatic field of the np particles
   type(scalar_var),  allocatable, save :: var_2(:)           ! sampling of the force between the two particles
   integer(4)      ,  allocatable, save :: ipnt(:,:,:)
   integer(4), save :: isamp2
   real(8), save :: xsamp(npoint), ysamp(npoint), wsamp(npoint), elfield2s2(3,nlayer,npoint), zsamp(3), dzcut
   character(4) :: chila, chipn, chidir
   integer(4) :: iptmacroion, ipoint, ip, ilayer, ivar, idir
   real(8) :: elfield(3), elfieldsum(3), integrand(nlayer,npoint), mffelec(nlayer), dmffelec(2), fac
   real(8) :: dx, dy, dz, dr2, dr, norm

   select case (iStage)
   case (iReadInput)

      iptmacroion = 1                                                                 ! select particle type
      if (txbc /= 'cyl') call Stop (txroutine,'txbc /= ''cyl''',uout)                 ! cylindrical geometry
      if (nppt(iptmacroion) /= 2) call Stop (txroutine, 'nppt(iptmacroion) /=2',uout) ! 2 particles of type iptmacroion

      allocate(var(nvar), var_2(nvar_2), ipnt(3, npoint, nlayer))
      ipnt = 0

      ivar = 0                                                               ! set var()%label, var()%norm
      do ilayer = 1, nlayer
         write(chila,'(i4)') ilayer
         do ipoint = 1, npoint
            write(chipn,'(i4)') ipoint
            do idir = 1, 3
               write(chidir,'(i4)') idir
               ivar = ivar + 1
               var(ivar)%label = 'layer'//' '//chila//'   '//'point'//' '//chipn//' '//'dir'//' '//chidir
               var(ivar)%norm = One
               ipnt(idir,ipoint,ilayer) = ivar
            end do
         end do
      end do

      var_2%label = ['force(z) (from mean el field) at -z', 'force(z) (from mean el field) at +z']
                                                      ! set var_2()%label
      dzcut = 0.1d0                                   ! cutoff between integration ponts and part
      zsamp(1:3) = [Zero, -cyllen/2, cyllen/2]        ! z-values of the integration ponts

   case (iWriteInput)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarSample(iStage, 1, nvar_2, var_2)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var, var_2

      call IntegCirPoints(npoint, cylrad, xsamp, ysamp, wsamp)  ! get integration poins and their weights

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarSample(iStage, 1, nvar_2, var_2)

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarSample(iStage, 1, nvar_2, var_2)

      isamp2 = 0                                                ! initialize variable for sumation of
      elfield2s2(1:3,1:nlayer,1:npoint) = Zero                  ! square of the electrical field from the part

   case (iSimulationStep)

! ... sample the electrical field at layer (1) z = 0, (2) z = -cyllen/2, (3) z = +cyllen/2

      isamp2 = isamp2 + 1
      do ipoint = 1, npoint
         do ilayer = 1, nlayer
            elfieldsum(1:3) = Zero
            do ip = 1, np
               dx = xsamp(ipoint) - ro(1,ip)
               dy = ysamp(ipoint) - ro(2,ip)
               dz = zsamp(ilayer) - ro(3,ip)
               if (abs(dz) < dzcut) then                        ! to avoid singularities or near singularities
     !            write(*,*) 'ipoint, ilayer, ip, dz', ipoint, ilayer, ip, dz
               else
                  dr2 = dx**2 + dy**2 + dz**2
                  dr = sqrt(dr2)
                  fac = az(ip)/(dr2*dr)
                  elfield(1) = dx*fac
                  elfield(2) = dy*fac
                  elfield(3) = dz*fac
                  elfieldsum(1:3) = elfieldsum(1:3) + elfield(1:3)
               end if
            end do
            var(ipnt(1:3,ipoint,ilayer))%value = elfieldsum(1:3)! sample electrical field of the particles
            elfield2s2(1:3,ilayer,ipoint) = elfield2s2(1:3,ilayer,ipoint) + elfieldsum(1:3)**2
         end do

      end do

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      norm = InvInt(isamp2)                                  ! square of the electrical field from the particles
      elfield2s2(1:3,1:nlayer,1:npoint) = elfield2s2(1:3,1:nlayer,1:npoint)*norm ! ... averaged over macrosteps
      integrand(1:nlayer,1:npoint) = &
           Half*(elfield2s2(3,1:nlayer,1:npoint) - elfield2s2(1,1:nlayer,1:npoint) - elfield2s2(2,1:nlayer,1:npoint))
      do ilayer = 1, nlayer
         mffelec(ilayer) = Pi*cylrad**2 * sum(integrand(ilayer,1:npoint)*wsamp(1:npoint))/FourPi
      end do
      dmffelec(1) = mffelec(1)-mffelec(2)                          ! sum over surfaces that contribute
      dmffelec(2) =-mffelec(1)+mffelec(3)                          ! sum over surfaces that contribute

      call WriteHead(2, 'integration of electrostatic field squared', uout)
      write(uout,'(a)') 'plane, point   x          y           w          elfield2s2(1:3)      integrand'
      write(uout,'(a)') '-----  -----   -          -           -          ---------------      ---------'
      write(uout,'(2i5,7f10.4)') ((ilayer,ipoint,xsamp(ipoint),ysamp(ipoint),wsamp(ipoint), &
                 elfield2s2(1:3,ilayer,ipoint),integrand(ilayer,ipoint),ipoint = 1,npoint), ilayer = 1, nlayer)
      write(uout,*)
      write(uout,'(a,3g15.5)') 'mean-field el stress (zz) (of the three layers) =', mffelec(1:nlayer)
      write(uout,'(a,g15.5)') (var_2(idir)%label//' =', dmffelec(idir), idir = 1, 2)

      var_2(1:2)%nsamp2 = var_2(1:2)%nsamp2 + 1
      var_2(1:2)%avs2 = dmffelec(1:2)                              ! store force
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarSample(iStage, 1, nvar_2, var_2)
      if (lsim .and. master) write(ucnf) var, var_2

   case (iAfterSimulation)

      call WriteHead(2, txheading, uout)
      call ScalarSample(iStage, 1, nvar, var)
      call ScalarSample(iStage, 1, nvar_2, var_2)
      call ScalarWrite(iStage, 1, nvar_2, var_2, 1, '(a,t35,4f15.5,f15.0)', uout)    ! write force
      var_2(1:2)%label = var_2(1:2)%label//' (kT/Å)'
      var_2(1:2)%norm = EpsiFourPi*beta
      call ScalarNorm(iStage, 1, nvar_2, var_2, 0)                                   ! change normalization
      write(uout,*)
      call ScalarWrite(iStage, 1, nvar_2, var_2, 1, '(a,t35,4f15.5,f15.0)', uout)    ! write force again
      write(uout,*)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)        ! write electrical field

      deallocate(var, var_2, ipnt)

   end select

end subroutine MeanElFieldZCyl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **DomainModule**
!! *module for domain analysis*
!************************************************************************


module DomainModule
   use MolModule
   integer(4),    parameter :: mndommax = 50
   integer(4),    parameter :: ntype = 1
   type(static1D_var), save :: vtype(ntype)
   character(10)            :: txdomaindef            ! selection of definition of domains
   integer(4)               :: npdomainsearch         ! number of particles used for establishing domian properties
   real(8)                  :: gkfac                  ! factor for determining lower limit of gk of domains
   integer(4)               :: ndommax                ! maximum number of domains to be considered
   integer(4)               :: ndom                   ! number of domains
   integer(4)               :: ipdom(mndommax)        ! id of particle defining a domain
   real(8)                  :: gkdom(mndommax)        ! Kirkwood Gk factor of the domain
   real(8)                  :: raddom(mndommax)       ! radius of the domain
   real(8)                  :: dipdom(3,mndommax)     ! dipole vector of the domain
end module DomainModule

!************************************************************************
!*                                                                      *
!*     DomainDriver                                                     *
!*                                                                      *
!************************************************************************

! ... driver of domain analyse based on Kirkwood gk-factor

subroutine DomainDriver(iStage)

   use DomainModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DomainDriver'

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)
      call DomainIO(iStage)
      call DomainDistFunc(iStage)
      call TCFDipDomain(iStage)
   case (iWriteInput)
      call CalcDomain
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   case (iBeforeSimulation)
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   case (iBeforeMacrostep)
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   case (iSimulationStep)
      call CalcDomain
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   case (iAfterMacrostep)
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   case (iAfterSimulation)
      call DomainIO(iStage)
      call DomainDominating(iStage)
      call DomainDistFunc(iStage)
      call UDomain(iStage)
      call TCFDipDomain(iStage)
   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine DomainDriver

!************************************************************************
!> \page moluser moluser.F90
!! **DomainIO**
!! *domain analyse based on Kirkwood gk-factor*
!************************************************************************


subroutine DomainIO(iStage)

   use DomainModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DomainIO'
   character(80), parameter :: txheading ='IO domain analysis'
   integer(4), save :: ndomsum, nsamp2, nsum
   integer(4) :: i, idom

   namelist /nmlDomain/ vtype, txdomaindef, npdomainsearch, gkfac, ndommax

   select case (iStage)
   case (iReadInput)

      vtype(1)%min = Zero
      vtype(1)%max =+10.0d0
      vtype(1)%nbin = 100
      txdomaindef = 'first_max'
      npdomainsearch = 0
      gkfac = Half
      ndommax = 20

      rewind(uin)
      read(uin,nmlDomain)

      call LowerCase(txdomaindef)
      if (npdomainsearch == 0) npdomainsearch = np

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)
      if (npdomainsearch > np) call Stop(txroutine, 'npdomainsearch > np', uout)
      if (ndommax > mndommax) call Stop(txroutine, 'ndommax > mndommax', uout)

   case (iWriteInput)

   case (iBeforeSimulation)

      ndomsum = 0
      nsamp2 = 0
      nsum = 0
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) ndomsum, nsamp2, nsum

      if (master) then
         open(unit=70,file='gkr_extensive')
         open(unit=71,file='gkr_short')
         open(unit=72,file='gkr_short_dom')
         open(unit=73,file='gk_sorted')
      end if

      if (lsim .and. master .and. txstart == 'continue') then
         do i = 1, ndomsum
            read(70,*)
            read(71,*)
         end do
         do i = 1, nsamp2
            read(72,*)
         end do
         do i = 1, nsum
            read(73,*)
         end do
      end if

   case (iBeforeMacrostep)

   case (iSimulationStep)

      ndomsum = ndomsum + ndom
      nsamp2 = nsamp2 + 1
      if (ndom > 3) nsum = nsum + 1

      if (master) then
         write(70,'(2i8,8f12.3)') &
         (idom, ipdom(idom), ro(1:3,ipdom(idom)), raddom(idom), gkdom(idom), dipdom(1:3,idom), idom = 1, ndom)
         write(71,'(f12.3,a,f12.3)') (raddom(idom), char(9), gkdom(idom), idom = 1, ndom)
         write(72,'(f12.3,a,f12.3)') raddom(1), char(9), gkdom(1)
         if (ndom > 3) write(73,'(i8,4(a,f12.3))') ndom, (char(9), gkdom(idom), idom = 1,  4)
      end if

   case (iAfterMacrostep)

      if (lsim .and. master) write(ucnf) ndomsum, nsamp2, nsum

   case (iAfterSimulation)

      if (master) then
         close(70)
         close(71)
         close(72)
         close(73)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,a   )') 'txdomaindef                      = ', txdomaindef
         write(uout,'(a,i8  )') 'npdomainsearch                   = ', npdomainsearch
         write(uout,'(a,f8.2)') 'gkfac                            = ', gkfac
         write(uout,'(a,i8  )') 'ndommax                          = ', ndommax
         write(uout,*)
         write(uout,'(a,f8.2)') 'average number of domains        = ', float(ndomsum)/float(nsamp2)
         write(uout,*)
      end if

   end select

end subroutine DomainIO

!************************************************************************
!> \page moluser moluser.F90
!! **CalcDomain**
!! *calculate dipole domains based on Kirkwoods gk-factor (G Karlstrom)*
!************************************************************************


subroutine CalcDomain

   use DomainModule
   implicit none

   character(40), parameter :: txroutine ='CalcDomain'
   type(df_var),  allocatable, save :: var(:)
   integer(4) :: nvar, ivar, ip, jp, ibin, i, idom, iptemp(1), maxbin(1)
   real(8), allocatable, save :: gmax(:), rmax(:)
   real(8)    :: dx, dy, dz, r2, ac, angcos, gmaxval(1)
   !integer(4) :: ipstep

   nvar = 1
   allocate(var(nvar))
   var(1)%label = trim('Krik g2 domain')
   var(1)%min = vtype(1)%min
   var(1)%max = vtype(1)%max
   var(1)%nbin = vtype(1)%nbin

   call DistFuncSample(iWriteInput, nvar, var)
   if (.not.allocated(gmax)) then
      allocate(gmax(np), rmax(np))
      gmax = 0.0E+00
      rmax = 0.0E+00
   end if

#if defined (_PAR_)
    gmax = Zero
    rmax = Zero
#endif

   ivar = 1

   do ip = ipmyid(1), ipmyid(2), max(1, np/npdomainsearch)
!  do ip = ipmyid(1), ipmyid(2)

! ... calculate gk(r) for particle ip

      call DistFuncSample(iBeforeMacrostep, nvar, var)
      do jp = 1, np
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
         ibin = max(-1,min(floor(var(ivar)%bini*(sqrt(r2)-var(ivar)%min)),var(ivar)%nbin))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ac
      end do
      do ibin = 0, var(ivar)%nbin
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1) + var(ivar)%avs2(ibin)
      end do

! ... determine the maximum value of gk(r) and its location

      if (txdomaindef == 'first_max') then   ! first maximum
         do ibin = 1, var(ivar)%nbin - 3              ! continue until three following values all are smaller
            if (maxval(var(ivar)%avs2(ibin+1:ibin+3)) < var(ivar)%avs2(ibin)) then
               gmax(ip) = var(ivar)%avs2(ibin)
               rmax(ip) = var(ivar)%min +(ibin+0.5)*var(ivar)%bin
               exit
            end if
         end do
         if (ibin == var(ivar)%nbin - 2) then      ! should not come here
            gmax(ip) = var(ivar)%avs2(ibin)
            rmax(ip) = var(ivar)%min + (ibin+0.5)* var(ivar)%bin
         end if
      else if (txdomaindef == 'global_max') then  ! global maximum
         if (maxval(var(ivar)%avs2(1:var(ivar)%nbin)) >= Zero) then
            maxbin = maxloc(var(ivar)%avs2(1:var(ivar)%nbin))
            gmax(ip) = var(ivar)%avs2(maxbin(1))
            rmax(ip) = var(ivar)%min +(maxbin(1)+0.5)*var(ivar)%bin
         end if
      else
         call stop('CalcDomain', 'error in txdomaindef', uout)
      end if

      rmax(ip) = min(rmax(ip), half*sqrt(boxlen(1)**2+boxlen(2)**2+boxlen(3)**2)) ! limit rmax to half box diagonal
      if (txbc /= 'xyz') call Stop(txroutine, 'txbc /= ''xyz''', uout)

   end do

#if defined (_PAR_)
   call par_allreduce_reals(gmax, vaux, np)
   call par_allreduce_reals(rmax, vaux, np)
#endif

! ... eliminate those particles for which gmax < gkfac*gmax

   gmaxval = maxval(gmax(1:np))
   where (gmax(1:np) < gkfac*gmaxval(1)) gmax(1:np) = Zero

! ... save domain properties

  call CalcPartDipMom(vaux)

   idom = 0
   do i = 1, ndommax
      iptemp = maxloc(gmax(1:np))       ! particle id
      ip = iptemp(1)
      if (gmax(ip) <= Zero) exit         ! no candiages left
      idom = idom + 1                   ! update number of domains
      ipdom(idom) = ip                  ! save particle id
      gkdom(idom) = gmax(ip)            ! save gmax
      raddom(idom) = rmax(ip)           ! save rmax
      dipdom(1:3,idom) = Zero
      do jp = 1, np
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rmax(ip)**2) then                             ! particle jp inside domain
            gmax(jp) = Zero                                    ! eliminate particle jp from further domain considerations
            dipdom(1:3,idom) = dipdom(1:3,idom) + vaux(1:3,jp) ! add contribution to domain dipole moment
         end if
      end do
   end do
   ndom = idom

   deallocate(var)

 end subroutine CalcDomain

!************************************************************************
!> \page moluser moluser.F90
!! **DomainDominating**
!! *sample gk for the dominating domain of each configuration*
!************************************************************************


subroutine DomainDominating(iStage)

   use DomainModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine =''
   character(80), parameter :: txheading ='gk for dominating domain'
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4) :: ivar, ibin, ip, jp
   real(8)    :: dx, dy, dz, r2, ac, angcos

   if (slave) return   ! only master

   select case (iStage)
   case (iReadInput)

   case (iWriteInput)

      nvar = 1
      allocate(var(nvar))
      var(1)%label = trim('gk(r) dominating dom')
      var(1)%min = vtype(1)%min
      var(1)%max = vtype(1)%max
      var(1)%nbin = vtype(1)%nbin
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

! ... sample gk(r) for the particle identifying the most dominating domain

      var%nsamp2 = var%nsamp2 + 1
      ivar = 1
      ip = ipdom(1)
      do jp = 1, np
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
         ibin = max(-1,min(floor(var(ivar)%bini*(sqrt(r2)-var(ivar)%min)),var(ivar)%nbin))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ac
      end do

   case (iAfterMacrostep)

      ivar = 1
      do ibin = 0, var(ivar)%nbin
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1) + var(ivar)%avs2(ibin)
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

end subroutine DomainDominating

!************************************************************************
!> \page moluser moluser.F90
!! **DomainDistFunc**
!! *Domain distribution functions*
!************************************************************************


subroutine DomainDistFunc(iStage)

   use DomainModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine =''
   character(80), parameter :: txheading ='domain distribution functions'
   integer(4),                save :: nvar
   type(df_var), allocatable, save :: var(:)
   real(8)       :: R_hig, G_hig
   integer(4)    :: ivar, idom, ibin, nbinloc

   if (slave) return   ! only master

   select case (iStage)
   case (iReadInput)

   case (iWriteInput)

! ... determine R_hig and G_hig for the histograms

      R_hig = vtype(1)%max
      G_hig = GetG_hig()

! ... initiate

      nvar = 5
      allocate(var(nvar))
      nbinloc = 200
      var%label = ['number of domains    ', 'P(R_max) domains     ', 'P(G_max) domains     ', 'P(R_max) imp. domains', 'P(G_max) imp. domains']
      var%min   = [Half,         Zero,  Zero,  Zero,  Zero]
      var%max   = [nbinloc+Half, R_hig, G_hig, R_hig, G_hig]
      var%nbin  = 5*(nbinloc)

      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1
      ivar = 1
      ibin = max(-1,min(floor(var(ivar)%bini*(ndom-0.5-var(ivar)%min)),int(var(ivar)%nbin)))
      var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      do idom = 1, ndom
         ivar = 2
         ibin = max(-1,min(floor(var(ivar)%bini*(raddom(idom)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
         ivar = 3
         ibin = max(-1,min(floor(var(ivar)%bini*(gkdom(idom)-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      end do
      ivar = 4
      ibin = max(-1,min(floor(var(ivar)%bini*(raddom(1)-var(ivar)%min)),int(var(ivar)%nbin)))
      var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
      ivar = 5
      ibin = max(-1,min(floor(var(ivar)%bini*(gkdom(1)-var(ivar)%min)),int(var(ivar)%nbin)))
      var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call WriteHead(2, txheading, uout)
      call DistFuncSample(iStage, nvar, var)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

contains

!........................................................................

! ... return the value of GetG_hig (the values given here are based on previous simulations)

real(8) function GetG_hig()
   real(8) :: dipmom

!  dipmom = 0.5d0
   dipmom = sqrt(dipa(1,1,1)**2+dipa(2,1,1)**2+dipa(3,1,1)**2)

   if (lewald) then
      if (dipmom < 0.11) then
         GetG_hig = 100.0
      else if ((0.11 < dipmom) .and. (dipmom < 0.24)) then
         GetG_hig = 200.0
      else
         GetG_hig = 1000.0
      endif
   else ! (MI in mind)
      if (dipmom < 0.11) then
         GetG_hig = 100.0
      else if ((0.11 < dipmom) .and. (dipmom < 0.24)) then
         if (np == 1000) GetG_hig = 100.0
         if (np == 3000) GetG_hig = 200.0
         if (np == 10000) GetG_hig = 300.0
         if (np == 30000) GetG_hig = 500.0
         if (np == 100000) GetG_hig = 800.0
         if (np == 300000) GetG_hig = 1000.0
      else
         if (np == 1000) GetG_hig = 100.0
         if (np == 3000) GetG_hig = 300.0
         if (np == 10000) GetG_hig = 1000.0
         if (np == 30000) GetG_hig = 3000.0
         if (np == 100000) GetG_hig = 5000.0
         if (np == 300000) GetG_hig = 10000.0
      end if
   end if
end function GetG_hig

!........................................................................

end subroutine DomainDistFunc

!************************************************************************
!> \page moluser moluser.F90
!! **UDomain**
!! *domain-domain dipole-dipole interaction energy*
!************************************************************************


subroutine UDomain(iStage)

   use DomainModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine =''
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: svar(:)
   integer(4) :: ivar, i, j, ip, jp
   real(8)    :: dx, dy, dz, r2, r1i, dotip, dotjp, dot, udom

   if (slave) return   ! only master

   select case (iStage)
   case (iReadInput)

   case (iWriteInput)

      nvar = 1
      allocate(svar(nvar))
      svar%label = 'dom-dom dipole energy'
      svar%norm = One/np

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, svar)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) svar

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, svar)

   case (iSimulationStep)

      udom = Zero
      do i = 1, ndom
         ip = ipdom(i)
         do j = i+1, ndom
            jp = ipdom(j)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            r1i = One/sqrt(r2)
            dotip = dipdom(1,i)*dx + dipdom(2,i)*dy + dipdom(3,i)*dz
            dotjp = dipdom(1,j)*dx + dipdom(2,j)*dy + dipdom(3,j)*dz
            dot = dipdom(1,i)*dipdom(1,j) + dipdom(2,i)*dipdom(2,j) + dipdom(3,i)*dipdom(3,j)
            udom = udom - (Three*dotip*dotjp*r1i**2 - dot)*r1i**3
         end do
      end do
      udom = Epsi0FourPi*udom

      ivar = 1
      svar(ivar)%value = udom
      call ScalarSample(iStage, 1, nvar, svar)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, svar)
      if (lsim .and. master) write(ucnf) svar
      call ScalarNorm(iStage, 1, nvar, svar, 0)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, svar)
      call ScalarNorm(iStage, 1, nvar, svar, 0)
      call ScalarWrite(iStage, 1, nvar, svar, 1, '(a,t35,2f15.4,f15.0)', uout)

      deallocate(svar)

   end select

end subroutine UDomain

!************************************************************************
!> \page moluser moluser.F90
!! **TCFDipDomain**
!! *sample gk for the dominating domain of each configuration*
!************************************************************************


subroutine TCFDipDomain(iStage)

   use DomainModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine =''
   character(80), parameter :: txheading ='domain dipole moment tcf'
   real(8),    save :: raddom0(mndommax), dipdom0(1:3,mndommax)
   integer(4), save :: ndom0, ipdom0(mndommax)
   integer(4), save :: nvar
   integer(4), save :: its                 ! time step
   integer(4), save :: nti                 ! number of time intervals of a tcf
   integer(4), save :: iti                 ! time interval of tcf
   real(8),    allocatable, save :: c0(:), c(:,:), ctemp(:,:)
   integer(4), allocatable, save :: norm(:)
   real(8),    save :: rlow, rupp, dr, dri, dx, dy, dz, r2
   real(8) :: dipdomt(1:3,mndommax), anorm
   integer(4) :: idom, ip, jp, ivar

   if (slave) return   ! only master

   select case (iStage)
   case (iReadInput)

      nvar = 5
      nti = 20
      rlow = 5.0
      rupp = 30.0
      dr = (rupp-rlow)/nvar
      dri = One/dr

   case (iWriteInput)

      allocate(c0(nvar), c(nvar,0:nti), norm(nvar),ctemp(mndommax,0:nti))
      c0 = 0.0E+00
      c = 0.0E+00
      norm = 0
      ctemp = 0.0E+00

      c = Zero
      norm = 0

      its = 0
      iti = 0

      call Initialize
      call Accumulate

   case (iBeforeSimulation)

   case (iBeforeMacrostep)

   case (iSimulationStep)

      its = its + 1
      iti = iti + 1

! ... calculate the dipole moments of the domains at iti = 0

      call CalcPartDipMom(vaux)

      do idom = 1, ndom0
         ip = ipdom0(idom)
         dipdomt(1:3,idom) = Zero
         do jp = 1, np
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < raddom(idom)**2) then                          ! particle jp inside domain
               dipdomt(1:3,idom) = dipdomt(1:3,idom) + vaux(1:3,jp) ! add contribution to domain dipole moment
            end if
         end do
      end do
      call Accumulate

      iti = mod(its,nti)                      ! update iti
      if (iti == 0) then
         call Accumulate2
         call Initialize                      ! start of a new tcf
         call Accumulate
      end if

   case (iAfterMacrostep)

   case (iAfterSimulation)

! ... divide by norm

      do ivar = 1, nvar
         anorm = InvInt(norm(ivar))
         c(ivar,0:nti) = anorm*c(ivar,0:nti)
      end do

! ... normalize

      do ivar = 1, nvar
         c0(ivar) = c(ivar,0)
         if (c0(ivar) > Zero) then
            anorm = One/c0(ivar)
            c(ivar,0:nti) = anorm*c(ivar,0:nti)
         end if
      end do

      call WriteHead(2, txheading, uout)
      write(uout,'(a,i8  )') 'nvar (number of tcfs)            = ', nvar
      write(uout,'(a,f8.2)') 'rlow                             = ', rlow
      write(uout,'(a,f8.2)') 'rupp                             = ', rupp

      write(uout,*)
      write(uout,'(12x,5i10)') norm
      write(uout,'(12x,5f10.3)') c0
      do iti = 0, nti
         write(uout,'(i10,2x,5f10.3)') iti, c(1:nvar,iti)
      end do

      if (ilist > 0) then
         write(ulist,'(a)') 'volume dipole tcf'
         write(ulist,'(i3)') nvar
         do ivar = 1, nvar
            write(ulist,'(i3)') ivar
            write(ulist,'(i3)') nti+1
            do iti = 0, nti
                                  ! assumes that tstep*istatic = 0.1
               write(ulist,'(f12.4,a,f12.4,a,f12.4)') float(iti)*0.1, char(9), c(ivar,iti), char(9), Zero
 !             write(ulist,'(f12.4,a,f12.4,a,f12.4)') float(iti)*One, char(9), c(ivar,iti), char(9), Zero   ! Jocke
            end do
         end do
      end if

      deallocate(c0, c, norm, ctemp)

   end select

contains

!........................................................................

subroutine Initialize
   ndom0 = ndom
   ipdom0(1:ndom) = ipdom(1:ndom)
   raddom0(1:ndom) = raddom(1:ndom)
   dipdom0(1:3,1:ndom) = dipdom(1:3,1:ndom)
   ctemp = Zero
   dipdomt(1:3,1:ndom0) = dipdom0(1:3,1:ndom0)
end subroutine Initialize

subroutine Accumulate
   integer(4) :: idom
   character(40), parameter :: txroutine ='Accumuate'
   if ((iti < 0) .or. (iti > nti)) call Stop(txroutine, 'error in iti', uout)
   do idom = 1, ndom0
      ctemp(idom,iti) = ctemp(idom,iti) + sum(dipdom0(1:3,idom)*dipdomt(1:3,idom))
   end do
end subroutine Accumulate

subroutine Accumulate2
   integer(4) :: ivar, idom
   do idom = 1, ndom0
      ivar = One + (raddom0(idom)-rlow)*dri
      if (ivar >= 1 .and. ivar <= nvar) then
         norm(ivar) = norm(ivar) + 1
         c(ivar,0:nti) = c(ivar,0:nti) + ctemp(idom,0:nti)
      end if
   end do
end subroutine Accumulate2

!........................................................................

end subroutine TCFDipDomain


!************************************************************************
!> \page moluser moluser.F90
!! **SCDF**
!! *calculate single chain distribution functions as a function of z coordinate*
!************************************************************************


!     type  label  quantity
!     ----  -----  --------
!     1     zden   number density in the z-direction
!     2     zden2  number density in the z-direction
!     3     rg     radius of gyration as a function of zcom
!     4     rg_xy  parallel radius of gyration as a function of zcom
!     5     rg_z   perpendicular radius of gyration as a function of zcom

!     min, max, and nbin of types 3-5 are coupled to type 2

subroutine SCDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SCDF'
   character(80), parameter :: txheading ='single chain distribution functions'
   integer(4)   , parameter :: ntype = 5
   type(static1D_var),         save :: vtype(ntype)
   logical,                    save :: lsym
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),                    save :: facsym
   type(chainprop_var),  save :: ChainProperty
   integer(4)  :: itype, ivar, ibin, ic, ict
   real(8)     :: zcom, norm, InvFlt

   namelist /nmlSCDF/ vtype, lsym

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      vtype%l =.true.
      vtype%min =-boxlen2(3)
      vtype%max =+boxlen2(3)
      vtype%nbin = 100
      lsym = .false.

      rewind(uin)
      read(uin,nmlSCDF)

      do itype = 3, 5
         if (vtype(itype)%l .and. .not.vtype(2)%l) call Stop(txroutine, 'zden (itype = 2) is needed for normalization', uout)
      end do
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

      facsym = One
      if (lsym) facsym = Two    ! chains with negative z-coordinate are reflected at z = 0

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['zden ','zden2','rg   ','rg_xy','rg_z ']
      vtype%nvar = nct

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nctct,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ivar + 1
               ipnt(ict,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

      do ic = 1, nc
         ict = ictcn(ic)
         if (ict <= 0) cycle
         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)

   ! ... summation of type 1

         itype = 1
         if(vtype(itype)%l) then
            zcom = ChainProperty%ro(3)
            if (lsym) zcom = abs(zcom)
            ivar = ipnt(ict,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(zcom-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

   ! ... summation of type 2

         itype = 2
         if(vtype(itype)%l) then
            zcom = ChainProperty%ro(3)
            if (lsym) zcom = abs(zcom)
            ivar = ipnt(ict,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(zcom-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

   ! ... summation of type 3

         itype = 3
         if(vtype(itype)%l) then
            ivar = ipnt(ict,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(zcom-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + sqrt(ChainProperty%rg2)
         end if

   ! ... summation of type 4

         itype = 4
         if(vtype(itype)%l) then
            ivar = ipnt(ict,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(zcom-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + sqrt(ChainProperty%rg2xy)
         end if

   ! ... summation of type 5

         itype = 5
         if(vtype(itype)%l) then
            ivar = ipnt(ict,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(zcom-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + sqrt(ChainProperty%rg2z)
         end if
      end do

   case (iAfterMacrostep)

      do ict = 1, nct
         do itype = 1, 2
            if (vtype(itype)%l) then
               ivar = ipnt(ict,itype)
               norm = One/(facsym*boxlen(1)*boxlen(2)*var(ivar)%bin)
               do ibin = 0, var(ivar)%nbin
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm
               end do
            end if
         end do
         do itype = 3, 5
            if (vtype(itype)%l) then
               ivar = ipnt(ict,itype)
               norm = var(ivar)%nsamp2/(facsym*boxlen(1)*boxlen(2)*var(ivar)%bin)
               do ibin = 0, var(ivar)%nbin
                  var(ivar)%avs2(ibin) = (var(ivar)%avs2(ibin)*norm)*InvFlt(var(ipnt(ict,2))%avs2(ibin))
               end do
            end if
         end do
      end do

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

end subroutine SCDF

!************************************************************************
!> \page moluser moluser.F90
!! **AdsRadGyr**
!! *calculate radius of gyration distribution functions of adsorbed and adsorbing chains*
!************************************************************************

!
!     type  label  quantity
!     ----  -----  --------
!     1     rg A    radius of gyration df for adsorbED chains
!     2     rg_z A  perpendicular radius of gyration df for adsorbED chains
!     3     rg_xy A parallel radius of gyration df for adsorbED chains
!     4     rg F    radius of gyration df for adsorbING chains
!     5     rg_z F  perpendicular radius of gyration df for adsorbING chains
!     6     rg_xy F parallel radius of gyration df for adsorbING chains

subroutine AdsRadGyr(iStage)

   use Molmodule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AdsRadGyr'
   character(80), parameter :: txheading ='radius of gyration df of adsorbed and adsorbing chains'
   integer(4)   , parameter :: ntype = 6
   type(static1D_var),         save :: vtype(ntype)
   type(adscond_var),          save :: adscond
   integer(4)    ,             save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   logical,       allocatable, save :: ladschainold(:)
   logical,       allocatable, save :: ladsseg(:)
   logical       , save       :: ladschain
   type(chainprop_var)        :: ChainProperty
   integer(4)                 :: ic, ict, itype, ivar, ibin
   real(8)                    :: value

   namelist /nmlAdsRadGyr/ vtype, adscond

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = Zero
      vtype%max = 100.0d0
      vtype%nbin = 200
      adscond = adscond_var('plane_xy', '-|+', 6.0d0)

      rewind(uin)
      read(uin,nmlAdsRadGyr)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

      if (.not.allocated(ladschainold)) then
         allocate(ladschainold(nc), ladsseg(maxval(npct(1:nct))))
         ladschainold = .false.
         ladsseg = .false.
      end if

! ... set remaining elements of vtype

      vtype%label = ['rg A   ','rg_z A ','rg_xy A','rg F   ','rg_z F ','rg_xy F']
      vtype%nvar = nct

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nct,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ivar+1
               ipnt(ict,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      ladschainold = .false.
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) ladschainold, var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc
         ict = ictcn(ic)
         call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
         if (ladschain) then
            call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
            call CalcChainProperty(ic, vaux, ChainProperty)
         end if

         if (ladschain) then                                         ! adsorbed chain
            do itype = 1, 3
               if(vtype(itype)%l) then
                  ivar = ipnt(ict,itype)
                  if (itype == 1) then
                     value = sqrt(ChainProperty%rg2)
                  else if (itype == 2) then
                     value = sqrt(ChainProperty%rg2z)
                  else if (itype == 3) then
                     value = sqrt(ChainProperty%rg2xy)
                  end if
                  ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
               end if
            end do
         end if

         if (ladschain .and. .not.ladschainold(ic)) then                  ! adsorbING chain
            ladschainold(ic) = .true.
            do itype = 4, 6
               if(vtype(itype)%l) then
                  ivar = ipnt(ict,itype)
                  if (itype == 4) then
                     value = sqrt(ChainProperty%rg2)
                  else if (itype == 5) then
                     value = sqrt(ChainProperty%rg2z)
                  else if (itype == 6) then
                     value = sqrt(ChainProperty%rg2xy)
                  end if
                  ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
               end if
            end do
         end if
         if (.not.ladschain .and. ladschainold(ic)) then                  ! desorbING chain
            ladschainold(ic) = .false.
         end if

      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) ladschainold, var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t40,a    )') 'adsorbing object               = ', adscond%txobject
      write(uout,'(a,t40,a    )') 'adsorbing surface              = ', adscond%txsurface
      write(uout,'(a,t35,f10.3)') 'adsorbing distance             = ', adscond%dist
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

      deallocate(var, ipnt)
      deallocate(ladschainold, ladsseg)

   end select

end subroutine AdsRadGyr

!************************************************************************
!> \page moluser moluser.F90
!! **AdsBondOrder**
!! *calculate bond order distribution function for adsorbed chains*
!************************************************************************


subroutine AdsBondOrder(iStage)

   use Molmodule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AdsBondOrder'
   character(80), parameter :: txheading ='bond order distribution function for adsorbed chains'
   integer(4)   , parameter :: ntype = 3
   integer(4)   , parameter :: mnbond = 1000
   type(static1D_var),         save :: vtype(ntype)
   real(8)           ,         save :: radius(ntype)
   type(adscond_var) ,         save :: adscond
   integer(4)        ,         save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   logical,       allocatable, save :: ladschain(:), ladsseg(:)
   integer(4),    allocatable, save :: bondlist(:,:), nq(:)
   real(8),       allocatable, save :: bondorder(:), rdir(:,:), rc(:,:), q(:,:)
   integer(4)                 :: ic, ict, itype, ivar, ibin, ib, nbond, isum
   type(scalar_var), allocatable, save :: var_s(:)

   namelist /nmlAdsBondOrder/ vtype, adscond, radius

   if (slave) return   ! only master

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = Zero
      vtype%max = One
      vtype%nbin = 100
      radius(1:ntype) = [10.0, 20.0, 50.0]
      adscond = adscond_var('plane_xy', '-|+', 6.0d0)

      rewind(uin)
      read(uin,nmlAdsBondOrder)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

   if (.not.allocated(ladschain)) then
      allocate(ladschain(nc), ladsseg(maxval(npct(1:nct))))
      ladschain = .false.
      ladsseg = .false.
   end if
   if (.not.allocated(bondorder)) then
      allocate(bondorder(mnbond), bondlist(2,mnbond), nq(mnbond), rdir(3,mnbond), rc(3,mnbond), q(6,mnbond))
      bondorder = 0.0E+00
      bondlist = 0
      nq = 0
      rdir = 0.0E+00
      rc = 0.0E+00
      q = 0.0E+00
   end if

! ... set remaining elements of vtype

      vtype%nvar = nct

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nct,ntype), var_s(nvar))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            write(vtype(itype)%label,'(f5.1)') radius(itype)
            vtype(itype)%label = 'rad ='//trim(vtype(itype)%label)
            do ict = 1, nct
               ivar = ivar+1
               ipnt(ict,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//', '//trim(txct(ict))
               var(ivar)%min = vtype(ict)%min
               var(ivar)%max = vtype(ict)%max
               var(ivar)%nbin = vtype(ict)%nbin
               var_s(ivar)%label = var(ivar)%label
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, var_s)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var_s

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, var_s)

   case (iSimulationStep)

      do ic = 1, nc
         call CheckAdsChainSeg(ic, adscond, ladschain(ic), ladsseg)
      end do

      var%nsamp2 = var%nsamp2 + 1
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ipnt(ict,itype)
               call CalcBondOrder(ict, ladschain, radius(itype), nbond, mnbond, bondorder, bondlist, nq, rdir, rc, q)
               do ib = 1, nbond
                  if (bondorder(ib) > Zero) then
                     ibin = max(-1,min(floor(var(ivar)%bini*(bondorder(ib)-var(ivar)%min)),int(var(ivar)%nbin)))
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
                  end if
               end do
               isum = 0
               var_s(ivar)%value = Zero
               do ib = 1, nbond
                  if (bondorder(ib) > Zero) then
                     isum = isum + 1
                     var_s(ivar)%value = var_s(ivar)%value + bondorder(ib)
                  end if
               end do
               var_s(ivar)%value = var_s(ivar)%value/isum
            end do
         end if
      end do
      call ScalarSample(iStage, 1, nvar, var_s)

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, var_s)
      if (lsim .and. master) write(ucnf) var
      if (lsim .and. master) write(ucnf) var_s

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, var_s)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t40,a    )') 'adsorbing object               = ', adscond%txobject
      write(uout,'(a,t40,a    )') 'adsorbing surface              = ', adscond%txsurface
      write(uout,'(a,t35,f10.3)') 'adsorbing distance             = ', adscond%dist
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      write(uout,*)
      call ScalarWrite(iStage, 1, nvar, var_s, 1, '(a,t35,2f15.4,f15.0)', uout)

      deallocate(var, ipnt, var_s)
      deallocate(ladschain, ladsseg)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AdsBondOrder

!************************************************************************
!> \page moluser moluser.F90
!! **AdsPropDyn**
!! *write chain adsorption dynamic data on external unit*
!************************************************************************


subroutine AdsPropDyn(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AdsPropDyn'
   character(80), parameter :: txheading ='write dynamic properties of adsorbed chains on external unit'
   integer(4)   , parameter :: ntype = 18
   integer(4)   , parameter :: mnobj = 100
   integer(4)   , parameter :: unit = 30
   real(8)      , parameter :: ScaleTime = 1.0d-3
   character(8) , parameter :: txtype(ntype) = ['time    ', 'Nads_p  ', 'Nads_seg', &
                                                'r_g     ', 'r_gxy   ', 'r_gz    ', &
                                                'n_loop  ', 'n_tail  ', 'n_train ', &
                                                'l_loop  ', 'l_tail  ', 'l_train ', &
                                                'ns_loop ', 'ns_tail ', 'ns_train', &
                                                'b.o._10 ', 'b.o._20 ', 'b.o._50 ']
   type(adscond_var),    save :: adscond
   real(8), allocatable, save :: var(:,:)             ! summation array of the variables
   integer(4)          , save :: nvar                 ! number of variables
   real(8)             , save :: blocklen             ! number of samplings in a block
   real(8)             , save :: facblocklen          ! factor for increasing blocklen
   integer(4)          , save :: isum                 ! current sampling number
   integer(4)          , save :: nblock               ! number of blocks sampled
   logical, allocatable, save :: ladschain(:)         ! .true. if chain is adsorbed
   logical, allocatable, save :: ladsseg(:)           ! .true. if segment is adsorbed
   type(chainprop_var) :: ChainProperty
   integer(4)  :: nobj(3)                             ! number of objects
   integer(4)  :: lenobj(3,mnobj)                     ! length of objects
   integer(4)  :: nsegobj(3)                          ! number of segments involved
   real(8), save :: radius(3) = [ 10.0d0, 20.0d0, 50.0d0 ]
   real(8)     :: avbondorder
   integer(4)  :: ic, ict, i, m
   real(8)     :: GetTime, MomTime

   namelist /nmlAdsPropDyn/ adscond, blocklen, facblocklen

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      adscond = adscond_var('plane_xy', '-|+', 6.0d0)
      blocklen = 10
      facblocklen = 1.1d0

      rewind(uin)
      read(uin,nmlAdsPropDyn)
      if (lana .and. lbd) call IOBD(iStage)                                     ! to read tstep

      if (.not.allocated(ladschain)) then
         allocate(ladschain(nc), ladsseg(maxval(npct(1:nct))))
         ladschain = .false.
         ladsseg = .false.
      end if

      nvar = ntype
      allocate(var(nvar,nct))
      var = 0.0E+00

   case (iBeforeSimulation)

      isum = 0
      nblock = 0
      var = Zero
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) isum, blocklen, nblock, var
      if (master) then
         do ict = 1, nct
            open(unit+ict-1, file = 'ads_'//trim(txct(ict)))
            if (lsim .and. txstart == 'continue') then                    ! advance unit 1 + nblock rows
               do m = 1, 1 + nblock
                  read(unit+ict-1,*)
               end do
            else                                                          ! write heading
               write(unit+ict-1,'(18(a15,a))') (txtype(i),tab,i = 1,nvar)
            end if
         end do
      end if

   case (iSimulationStep)

      isum = isum + 1
      MomTime = ScaleTime*GetTime()
      var(1,1:nct) = var(1,1:nct) + MomTime                                ! sample time
      do ic = 1, nc                                                        ! loop over chains
         ict = ictcn(ic)
         call CheckAdsChainSeg(ic, adscond, ladschain(ic), ladsseg)
         if (ladschain(ic)) then
            call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)           ! undo the periodic boundary condtion
            call CalcChainProperty(ic, vaux, ChainProperty)                ! calculate chain property
            call CalcLTT(npct(ict),ladsseg,mnobj,nobj,lenobj,nsegobj)      ! calculate nobj, lenobj, and nsegobj
            var(2,ict) = var(2,ict) + 1                                    ! sample number of adsorbed chains
            var(3,ict) = var(3,ict) + count(ladsseg(1:npct(ict)))          ! sample number of adsorbed segments
            var(4,ict) = var(4,ict) + ChainProperty%rg2
            var(5,ict) = var(5,ict) + ChainProperty%rg2xy
            var(6,ict) = var(6,ict) + ChainProperty%rg2z
            do i = 1, 3                                                    ! loop over the different characteristics
               var(i+6,ict) = var(i+6,ict) + nobj(i)
               var(i+9,ict) = var(i+9,ict) + sum(lenobj(i,1:nobj(i)))
               var(i+12,ict) = var(i+12,ict) + nsegobj(i)
            end do
         end if
      end do
      do ict = 1, nct
         do i = 1, 3
            call CalcAvBondOrder(ict, ladschain, radius(i), avbondorder)
            var(i+15,ict) = var(i+15,ict) + avbondorder
         end do
      end do

      if (isum >= blocklen-1.0d-20) then                                          ! sample if current isum > current block length
         blocklen = blocklen*facblocklen
         nblock = nblock + 1
         do ict = 1, nct
            var(1,ict) = var(1,ict)/isum
            if (var(2,ict) > 0) var(4:15,ict) = var(4:15,ict) / var(2,ict)        ! make average for adsorbed chain (rg and ltt)
            var(2:3,ict) = var(2:3,ict)/isum                                      ! normalize over number of samplings
            var(4:6,ict) = sqrt(var(4:6,ict))                                     ! make rms
            where(var(7:9,ict) > 0) var(10:12,ict) = var(10:12,ict)/var(7:9,ict)  ! make average length of objects
            var(16:18,ict) = var(16:18,ict)/isum                                  ! make average of bond order
            write(unit+ict-1,'(f15.4,a,17(f15.3,a))') (var(i,ict),tab,i = 1,nvar)
            var(1:nvar,ict) = Zero
         end do
         isum = 0
      end if

   case (iAfterMacrostep)

      if (lsim .and. master) write(ucnf) isum, blocklen, nblock, var

   case (iAfterSimulation)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t40,a    )') 'adsorbing object               = ', adscond%txobject
      write(uout,'(a,t40,a    )') 'adsorbing surface              = ', adscond%txsurface
      write(uout,'(a,t35,f10.3)') 'adsorbing distance             = ', adscond%dist
      write(uout,*)
      write(uout,'(a,4a25)') 'generated file(s): ', ('ads_'//trim(txct(ict)),ict = 1, nct)
      write(uout,'(a,t40,f10.3)') 'final block length             = ', blocklen
      write(uout,'(a,t40,f10.3)') 'block length factor            = ', facblocklen
      write(uout,'(a,t40,i10  )') 'number of averages             = ', nblock

      deallocate(var)
      deallocate(ladschain, ladsseg)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AdsPropDyn

!************************************************************************
!> \page moluser moluser.F90
!! **CalcAvBondOrder**
!! *calculation of the average bond order among provided chains*
!************************************************************************


subroutine CalcAvBondOrder(ictbond, ladschain, radius, avbondorder)

   use Molmodule
   implicit none

   integer(4)   , parameter :: mnbond = 1000
   character(40), parameter :: txroutine =''

   integer(4), intent(in)  :: ictbond          ! chain type
   logical,    intent(in)  :: ladschain(*)     ! chains to consider
   real(8),    intent(in)  :: radius           ! radius of region to consider
   real(8),    intent(out) :: avbondorder      ! average bond order of the bonds

   integer(4), allocatable, save :: bondlist(:,:), nq(:)
   real(8), allocatable, save :: bondorder(:), rdir(:,:), rc(:,:), q(:,:)
   integer(4) :: nbond, ib, ibsum

   if (.not.allocated(bondorder)) then
      allocate(bondorder(mnbond), bondlist(2,mnbond), nq(mnbond), rdir(3,mnbond), rc(3,mnbond), q(6,mnbond))
      bondorder = 0.0E+00
      bondlist = 0
      nq = 0
      rdir = 0.0E+00
      rc = 0.0E+00
      q = 0.0E+00
   end if

   call CalcBondOrder(ictbond, ladschain, radius, nbond, mnbond, bondorder, bondlist, nq, rdir, rc, q)
   ibsum = 0
   avbondorder = Zero
   do ib = 1, nbond
      if (bondorder(ib)  > Zero) then
        ibsum = ibsum + 1
        avbondorder = avbondorder + bondorder(ib)
      end if
   end do
   if (ibsum > 0) avbondorder = avbondorder/ibsum

end subroutine CalcAvBondOrder

!************************************************************************
!> \page moluser moluser.F90
!! **CalcBondOrder**
!! *calculation of the bond order among bonds in chains of type ict*
!************************************************************************


subroutine CalcBondOrder(ictbond, ladschain, radius, nbond, mnbond, bondorder, bondlist, nq, rdir, rc, q)

   use Molmodule
   implicit none
   character(40), parameter :: txroutine ='CalcBondOrder'

   integer(4), intent(in)  :: ictbond          ! chain type
   logical,    intent(in)  :: ladschain(*)     ! chains to consider
   real(8),    intent(in)  :: radius           ! radius of region to consider
   integer(4), intent(out) :: nbond            ! number of bonds considered
   integer(4), intent(in)  :: mnbond           ! maximum number of bonds
   real(8),    intent(out) :: bondorder(*)     ! bond order of the bonds
   integer(4) :: bondlist(2,*), nq(*)          ! work
   real(8)    :: rdir(3,*), rc(3,*), q(6,*)    ! work

   integer(4), save :: itestloc = 0
   integer(4) :: ic, ict, iseg, ip, jp, ib, jb, nrot
   real(8)    :: norm, dx, dy, dz, r2, qq(3,3), diagonal(3), eivr(3,3)

! ... make bondlist of adsorbed chains

      nbond = 0
      do ic = 1, nc
         if (ladschain(ic)) then
            ict = ictcn(ic)
            do iseg = 1, npct(ict)-1
               ip = ipnsegcn(iseg,ic)
               jp = ipnsegcn(iseg+1,ic)
               nbond = nbond+1
               if (nbond > mnbond) call Stop(txroutine, 'nbond > mnbond', 6)
               bondlist(1,nbond) = ip
               bondlist(2,nbond) = jp
            end do
         end if
      end do
      if (itestloc == 1) then
         write(*,*) 'bondlist'
         write(*,'((2i10))') (bondlist(1:2,ib),ib=1,nbond)
      end if

! ... calculate bond direction and location of bonds of particles that belongs to adsorbed chains

   do ib = 1, nbond
      ip = bondlist(1,ib)
      jp = bondlist(2,ib)
      rdir(1:3,ib) = r(1:3,jp)-r(1:3,ip)
      call PBC(rdir(1,ib),rdir(2,ib),rdir(3,ib))
      rc(1:3,ib) = r(1:3,ip) + Half*rdir(1:3,ib)
      call PBC(rc(1,ib),rc(2,ib),rc(3,ib))
      norm = 1/dsqrt(rdir(1,ib)**2+rdir(2,ib)**2+rdir(3,ib)**2)
      rdir(1:3,ib) = rdir(1:3,ib)*norm
   end do
   if (itestloc == 1) then
      write(*,*) ' ip        r                   '
      write(*,'((i5,3f10.3))') (ip, r(1:3,ip), ip = 1, na)
      write(*,*) ' ip        rc                   rdir'
      write(*,'((i5, 3f10.3,4x,3f10.3))') (ib, rc(1:3,ib), rdir(1:3,ib), ib = 1, nbond)
   end if

! ... calculate bond alingment tensor

   q(1:6,1:nbond) = zero
   nq(1:nbond) = 0
   do ib = 1, nbond
      ip = bondlist(1,ib)
      do jb = ib, nbond
         jp = bondlist(1,jb)
!        if (icnpn(ip) == icnpn(jp)) cycle                  ! exclude bonds in same chain
         dx = rc(1,ib) - rc(1,jb)
         dy = rc(2,ib) - rc(2,jb)
         dz = rc(3,ib) - rc(3,jb)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < radius**2) then
            if (itestloc == 1) write(*,'(a,2i5,2f10.3)') 'ib, jb, r2, radius**2' , ib, jb, r2, radius**2
            if (ictpn(ip) == ictbond) then
               nq(ib) = nq(ib) + 1
               q(1,ib) = q(1,ib) + Three*rdir(1,jb)**2-One
               q(2,ib) = q(2,ib) + Three*rdir(2,jb)**2-One
               q(3,ib) = q(3,ib) + Three*rdir(3,jb)**2-One
               q(4,ib) = q(4,ib) + Three*rdir(1,jb)*rdir(2,jb)
               q(5,ib) = q(5,ib) + Three*rdir(1,jb)*rdir(3,jb)
               q(6,ib) = q(6,ib) + Three*rdir(2,jb)*rdir(3,jb)
            end if
            if (ib == jb) cycle
            if (ictpn(jp) == ictbond) then
               nq(jb) = nq(jb) + 1
               q(1,jb) = q(1,jb) + Three*rdir(1,ib)**2-One
               q(2,jb) = q(2,jb) + Three*rdir(2,ib)**2-One
               q(3,jb) = q(3,jb) + Three*rdir(3,ib)**2-One
               q(4,jb) = q(4,jb) + Three*rdir(1,ib)*rdir(2,ib)
               q(5,jb) = q(5,jb) + Three*rdir(1,ib)*rdir(3,ib)
               q(6,jb) = q(6,jb) + Three*rdir(2,ib)*rdir(3,ib)
            end if
         end if
      end do
   end do
   if (itestloc == 1) write(*,'((a,i4,6f8.3))') ('nq, q',nq(ib),q(1:6,ib),ib=1,nbond)

! ... diagolalize and collect the largest eigenvalue

   do ib = 1, nbond
      bondorder(ib) = Zero
      if (nq(ib) == 0) cycle
      norm = Half/nq(ib)
      qq(1,1) = norm*q(1,ib)
      qq(2,2) = norm*q(2,ib)
      qq(3,3) = norm*q(3,ib)
      qq(1,2) = norm*q(4,ib)
      qq(2,1) = norm*q(4,ib)
      qq(1,3) = norm*q(5,ib)
      qq(3,1) = norm*q(5,ib)
      qq(2,3) = norm*q(6,ib)
      qq(3,2) = norm*q(6,ib)
      call Diag(3, qq, diagonal, eivr, nrot)
      bondorder(ib) = min(One-1.d-10, max(diagonal(1),diagonal(2),diagonal(3)))
      if (itestloc == 1) then
         write(*,'(a,2i5,3f10.3,2x,f10.3)') 'ib, nq(ib), eigenvalues, bond order', ib, nq(ib), qq(1,1), qq(2,2), qq(3,3), bondorder(ib)
      end if
   end do

end subroutine CalcBondOrder

!************************************************************************
!> \page moluser moluser.F90
!! **AdsEventDyn**
!! *write time and occurrence of adsorption events on FUSER*
!************************************************************************


subroutine AdsEventDyn(iStage)

   use Molmodule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AdsEventDyn'
   character(80), parameter :: txheading ='write time and occurrence of adsorption events'
   integer(4)   , parameter :: mnadsevent = 200
   integer(4)   , parameter :: unit = 25
   real(8)      , parameter :: ScaleTime = 1.0d-3
   type(adscond_var) , save :: adscond
   real(8)           , save :: iinterval
   real(8)           , save :: isum
   logical, allocatable, save :: ladschainold(:)              ! .true. if chain is adsorbed
   logical, allocatable, save :: ladsseg(:)                   ! .true. if segment is adsorbed
   logical             , save :: ladschain
   integer(4), allocatable , save :: nAdsEvent(:)
   real(8), allocatable, save :: AdsStart(:,:), AdsLength(:,:)
   integer(4)               :: i, ic, ict, icnlow, icnupp
   real(8), allocatable, save :: ResTime(:)
   real(8)                  :: GetTime, MomTime

   namelist /nmlAdsEventDyn/ adscond, iinterval

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      adscond = adscond_var('oplane_xy', '-|+', 6.0d0)
      iinterval = 1

      rewind(uin)
      read(uin,nmlAdsEventDyn)
      if (lana .and. lbd) call IOBD(iStage)                                     ! to read tstep

   case (iWriteInput)

      if (.not.allocated(ladschainold)) then
         allocate(ladschainold(nc),ladsseg(maxval(npct(1:nct))))
         ladschainold = .false.
         ladsseg = .false.
      end if
      if (.not.allocated(nAdsEvent)) then
         allocate(nAdsEvent(nc))
         nAdsEvent = 0
      end if
      if (.not.allocated(AdsStart)) then
         allocate(AdsStart(mnadsevent,nc), AdsLength(mnadsevent,nc))
         AdsStart = 0.0E+00
         AdsLength = 0.0E+00
      end if
      if (.not.allocated(ResTime)) then
         allocate(ResTime(nc))
         ResTime = 0.0E+00
      end if

   case (iBeforeSimulation)

      ladschainold = .false.
      nAdsEvent = 0
      AdsStart(1,1:nc) = Zero
      AdsLength(1,1:nc) = Zero
      if (txstart /= 'continue') then
         MomTime = Zero
         do ic = 1, nc
            ict = ictcn(ic)
            call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
            if (ladschain) then                         ! adsorbING chain
               ladschainold(ic) = .true.
               if (nAdsEvent(ic) <= mnadsevent) then                      ! prepare data to be written on uuser
                  nAdsEvent(ic) = nAdsEvent(ic) + 1
                  AdsStart(nAdsEvent(ic),ic) = MomTime
               end if
            end if
         end do
      end if
      if (lsim .and. master .and. txstart == 'continue') then
         read(ucnf) (ladschainold(ic), nAdsEvent(ic), AdsStart(1:nadsEvent(ic),ic), AdsLength(1:nAdsEvent(ic),ic), ic = 1, nc)
         do ic = 1, nc
         call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
         if (ladschain) ladschainold(ic) = .true.
         end do
      end if

   case (iBeforeMacrostep)

   case (iSimulationStep)

      isum = isum + 1
      if (isum >= iinterval) then
         MomTime = ScaleTime*GetTime()
         do ic = 1, nc
            ict = ictcn(ic)
            call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
            if (ladschain .and. .not.ladschainold(ic)) then                    ! adsorbING chain
               ladschainold(ic) = .true.
               if (nAdsEvent(ic) <= mnadsevent) then                      ! prepare data to be written on uuser
                  nAdsEvent(ic) = nAdsEvent(ic) + 1
                  AdsStart(nAdsEvent(ic),ic) = MomTime
               end if
            end if
            if (.not.ladschain .and. ladschainold(ic)) then                    ! desorbING chain
               ladschainold(ic) = .false.
               if (nAdsEvent(ic) <= mnadsevent) &                         ! prepare data to be written on uuser
               AdsLength(nAdsEvent(ic),ic) = MomTime - AdsStart(nAdsEvent(ic),ic)
            end if
         end do
         isum = 0
      end if

   case (iAfterMacrostep)

      if (master) write(ucnf) (ladschainold(ic), nAdsEvent(ic), AdsStart(1:nadsEvent(ic),ic), AdsLength(1:nAdsEvent(ic),ic), ic = 1, nc)

   case (iAfterSimulation)

! ... add contribution from chains adsorbed at the end of the simulation

      MomTime = ScaleTime*GetTime()
      do ic = 1, nc
         call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
         if (ladschain) then
             if(nAdsEvent(ic) > 0) then
                AdsLength(nAdsEvent(ic),ic) = MomTime - AdsStart(nAdsEvent(ic),ic)   ! prepare data to be written on uuser
             else
                call warn('txroutine', 'soft exception appeard; nAdsEvent(ic) < 1', uout)
             end if
         end if
      end do

! ... write data on uuser

      call FileOpen(uuser, fuser, 'form/noread')
      write(uuser,'(a,i8)') 'nct =',nct
      write(uuser,'(a,i8)') 'nc  =',nc
      write(uuser,*)
      write(uuser,'(7(a,2x))') 'chain type', 'chain nr', 'ads events', 'total ads time', 'ads start', 'ads length', '...'
      write(uuser,'(7(a,2x))') '----------', '--------', '----------', '--------------', '---------', '----------', '---'
      do ict = 1, nct
         icnlow = icnct(ict)
         icnupp = icnlow + ncct(ict)-1
         do ic = icnlow, icnupp
            write(uuser,'(3i10,5x,1002f9.3)') &
            ict, ic, nAdsEvent(ic), sum(AdsLength(1:nAdsEvent(ic),ic)), (AdsStart(i,ic), AdsLength(i,ic), i = 1, nAdsEvent(ic))
         end do
      end do
      close(uuser)

! ... calculate residence time

      do ict = 1, nct
         icnlow = icnct(ict)
         icnupp = icnct(ict) + ncct(ict) - 1
         ResTime(ict) = Zero
         do ic = icnlow, icnupp
            ResTime(ict) = ResTime(ict) + sum(AdsLength(1:nAdsEvent(ic),ic))
         end do
         ResTime(ict) = ResTime(ict)*InvInt(sum(nAdsEvent(icnlow:icnupp)))
      end do

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t40,a    )') 'adsorbing object               = ', adscond%txobject
      write(uout,'(a,t40,a    )') 'adsorbing surface              = ', adscond%txsurface
      write(uout,'(a,t35,f10.3)') 'adsorbing distance             = ', adscond%dist
      write(uout,*)
      write(uout,*) 'average residence time of chain type at surface'
      write(uout,'()')
      write(uout,'(a,3x,3(a10,5x))') 'chain type                = ', txct(1:nct)
      write(uout,'(a,3(f10.4,4x))')  'residence time            = ', ResTime(1:nct)

      deallocate(ladschainold, ladsseg)
      deallocate(nAdsEvent)
      deallocate(AdsStart, AdsLength)
      deallocate(ResTime)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AdsEventDyn

!************************************************************************
!> \page moluser moluser.F90
!! **AdsModule**
!! *module for adsorption module*
!************************************************************************


module AdsModule

   use MolModule

   integer(4), parameter :: mnAdsEvent = 500

   integer(4), allocatable :: Idt(:)             ! chain type
   integer(4), allocatable :: Id(:)              ! chain number
   integer(4), allocatable :: nAdsEvent(:)       ! number of adsorption events
   real(8)   , allocatable :: AdsStart(:,:)      ! start of adsorption event
   real(8)   , allocatable :: AdsLength(:,:)     ! end of adsorption event
   integer(4), allocatable :: IdtSave(:)         ! chain type
   integer(4), allocatable :: IdSave(:)          ! chain number
   integer(4), allocatable :: nAdsEventSave(:)   ! number of adsorption events
   real(8)   , allocatable :: AdsStartSave(:,:)  ! start of adsorption event
   real(8)   , allocatable :: AdsLengthSave(:,:) ! end of adsorption event

end module AdsModule

!************************************************************************
!> \page moluser moluser.F90
!! **AdsExam**
!! *call of AdsExam routines*
!************************************************************************


subroutine AdsExam

   use AdsModule
   implicit none

   character(80), parameter :: txheading ='analysis of adsorption data'

   call FileOpen(uin, fin , 'form/noread')
   call FileOpen(uuser, fuser, 'form/noread')
   call WriteHead(1, txheading, uout)
   call AdsExam1
   call AdsExam2
   call AdsExam3
   call AdsExam4
   close(uin)
   close(uuser)

end subroutine AdsExam

!************************************************************************
!> \page moluser moluser.F90
!! **AdsExam1**
!! *generate curvz file using the primary adsorption data*
!************************************************************************


subroutine AdsExam1

   use AdsModule
   implicit none

   character(80), parameter :: txheading ='generate curvz file after sorting primary ads. data'
   integer(4), save :: itestloc = 1

   call WriteHead(2, txheading, uout)
   call ReadPrimAdsData
   if (itestloc == 1) call WritePrimAdsData('primary adsorption data before time sorting:')
   call sort
   if (itestloc == 1) call WritePrimAdsData('primary adsoroption data after time sorting:')
   call PrepCurvz
   write(uout,*)
   write(uout,'(a)') 'Curvz file have been generated'

contains

!........................................................................

subroutine sort                                ! sort in ascending AdsStart(1,i)
   integer(4), allocatable :: index(:)
   integer(4) :: i, ict, ic, ilow, ioff, idx

   allocate(IdtSave(nc), IdSave(nc), nAdsEventSave(nc), AdsStartSave(mnAdsEvent,nc), AdsLengthSave(mnAdsEvent,nc))
   IdtSave = 0
   IdSave = 0
   nAdsEventSave = 0
   AdsStartSave = 0.0E+00
   AdsLengthSave = 0.0E+00
   allocate(index(nc))
   index = 0

   do i = 1, nc
      IdtSave(i) = Idt(i)
      IdSave(i) = Id(i)
      nAdsEventSave(i) = nAdsEvent(i)
      AdsStartSave(i,1) = Zero                 ! zero is assigned to chains never adsorbed
      AdsStartSave(i,1:nAdsEvent(i)) = AdsStart(1:nAdsEvent(i),i)
      AdsLengthSave(i,1:nAdsEvent(i)) = AdsLength(1:nAdsEvent(i),i)
   end do
   ilow = 1
   do ict = 1, nct
      call HeapSortIndex(ncct(ict), AdsStartSave(ilow,1), index(ilow))
      ilow = ilow + ncct(ict)
   end do
   i = 0
   ioff = 0
   do ict = 1, nct
      do ic = 1, ncct(ict)
         i = i + 1
         idx = ioff + index(ioff+ic)
         Idt(i) = IdtSave(idx)
         Id(i) = IdSave(idx)
         nAdsEvent(i) = nAdsEventSave(idx)
         AdsStart(1:nAdsEvent(i),i) = AdsStartSave(idx,1:nAdsEventSave(idx))
         AdsLength(1:nAdsEvent(i),i) = AdsLengthSave(idx,1:nAdsEventSave(idx))
      end do
      ioff = ioff + ncct(ict)
   end do
   deallocate(IdtSave, IdSave, nAdsEventSave, AdsStartSave, AdsLengthSave)
   deallocate(index)
end subroutine sort

!........................................................................

subroutine PrepCurvz                     ! generate input file to Curvz
   real(8),      parameter :: Zero = 0.0d0
   character(1), parameter :: tab = char(9)
   integer(4) ::  i, j, ic
   open(20,file = 'chainadsevent.curvz.txt')
   ic = 0
   do i = 1, nc
      if (nAdsEvent(i) > 0) then         ! write only chains that have been adsorbed
         ic = ic + 1                     ! order adsorbed chains from 1
         write(20,'(a)') '%%%%%%%%%%%%'
         write(20,'(f8.3,a,f8.1)') Zero, tab, Zero
         write(20,'(f8.3,a,f8.1)') AdsStart(1,i), tab, Zero
         do j = 1, nAdsEvent(i)-1
            write(20,'(f8.3,a,i8)') AdsStart(j,i), tab, ic
            write(20,'(f8.3,a,i8)') AdsStart(j,i)+AdsLength(j,i), tab, ic
            write(20,'(f8.3,a,f8.1)') AdsStart(j,i)+AdsLength(j,i), tab, Zero
            write(20,'(f8.3,a,f8.1)') AdsStart(j+1,i), tab, Zero
         end do
         write(20,'(f8.3,a,i8)') AdsStart(j,i), tab, ic
         write(20,'(f8.3,a,i8)') AdsStart(j,i)+AdsLength(j,i), tab, ic
         write(20,'(a)') '%%%%%%%%%%%%'
      end if
   end do
   close(20)
end subroutine PrepCurvz

!........................................................................

end subroutine AdsExam1

!************************************************************************
!> \page moluser moluser.F90
!! **AdsExam2**
!! *analysis of primary adsorption data; adsorption events*
!************************************************************************


subroutine AdsExam2

   use AdsModule
   implicit none

   character(40), parameter :: txroutine ='AdsExam2'
   character(80), parameter :: txheading ='adsorption related distribution function'
   integer(4), parameter :: ntype = 3
   character(10), parameter :: txtype(ntype) = ['adsorbing', 'desorbing', 'adsorbed ']
   real(8)        :: vmin, vmax
   integer(4)     :: nbin, nvar
   integer(4),   allocatable :: ipnt(:,:)
   type(df_var), allocatable :: var(:)
   integer(4)     :: itype, ivar, ibin, ivara, ivard, j, ict, ic
   character(2)   :: str

   namelist /nmlAdsExam2/ vmin, vmax, nbin

   call WriteHead(2, txheading, uout)
   call ReadPrimAdsData

   vmin = Zero
   vmax = 10.0
   nbin = 100

   rewind(uin)
   read(uin,nmlAdsExam2)

   if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

   nvar = ntype*nct
   allocate(var(nvar),ipnt(nct,ntype))
   ipnt = 0

   nvar = 0
   do itype = 1, ntype
      do ict = 1, nct
         nvar = nvar + 1
         ipnt(ict,itype) = nvar
         write(str,'(i2)') ict
         var(nvar)%label = trim(txtype(itype))//' '//'ict='//trim(adjustl(str))
         var(nvar)%min = vmin
         var(nvar)%max = vmax
         var(nvar)%nbin = nbin
      end do
   end do

   call DistFuncSample(2, nvar, var)
   call DistFuncSample(3, nvar, var)
   call DistFuncSample(4, nvar, var)

   var%nsamp2 = var%nsamp2+1

! ... sampling of type 1

   itype = 1
   do ic = 1, nc
      ict = Idt(ic)
      ivar = ipnt(ict,itype)
      do j = 1, nAdsEvent(ic)
         ibin = max(-1,min(floor(var(ivar)%bini*(AdsStart(j,ic)-var(ivar)%min)),var(ivar)%nbin))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
      end do
   end do

! ... sampling of type 2

   itype = 2
   do ic = 1, nc
      ict = Idt(ic)
      ivar = ipnt(ict,itype)
      do j = 1, nAdsEvent(ic)
         ibin = max(-1,min(floor(var(ivar)%bini*(AdsStart(j,ic)+AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
      end do
   end do

! ... sampling of type 3

   itype = 3
   do ict = 1, nct
      ivara = ipnt(ict,1)
      ivard = ipnt(ict,2)
      ivar = ipnt(ict,itype)
      var(ivar)%avs2(-1) = (var(ivara)%avs2(-1) - var(ivard)%avs2(-1))
      do ibin = 0, nbin
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1) + (var(ivara)%avs2(ibin) - var(ivard)%avs2(ibin))
      end do
   end do

! ... normalize to the number of adsorbing and desorpbing items per unit time

   do itype = 1, 2
      do ict = 1, nct
         ivar = ipnt(ict,itype)
         var(ivar)%avs2(-1:nbin) = var(ivar)%avs2(-1:nbin)*var(ivar)%bini
      end do
   end do

   call DistFuncSample(6, nvar, var)
   call DistFuncSample(7, nvar, var)
   call DistFuncHead(nvar, var, uout)
   call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

   deallocate(var, ipnt)

end subroutine AdsExam2

!************************************************************************
!> \page moluser moluser.F90
!! **AdsExam3**
!! *analysis of primary adsorption data; adsorption length*
!************************************************************************


subroutine AdsExam3

   use AdsModule
   implicit none

   character(40), parameter :: txroutine ='AdsExam2'
   character(80), parameter :: txheading ='ads. length dist. function; first, early, late, all'
   integer(4), parameter :: ntype = 4
   character(15), parameter :: txtype(ntype) = ['first ads event', 'early ads event', 'late ads event ', 'all ads event  ']

   logical        :: ltype(ntype)
   real(8)        :: tdiv, tmax, vmin, vmax
   integer(4)     :: nbin, nvar
   integer(4),   allocatable :: ipnt(:,:)
   type(df_var), allocatable :: var(:)
   integer(4)     :: itype, ivar, ibin, ict, ic, j
   character(2)   :: str

   namelist /nmlAdsExam3/ ltype, tdiv, tmax, vmin, vmax, nbin

   call WriteHead(2, txheading, uout)
   call ReadPrimAdsData


   ltype =.false.
   tdiv = 10.0
   vmin = Zero
   vmax = 100.0
   nbin = 100

   rewind(uin)
   read(uin,nmlAdsExam3)

   nvar = ntype*nct
   allocate(var(nvar),ipnt(nct,ntype))
   ipnt = 0

   if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

   nvar = 0
   do itype = 1, ntype
      if (ltype(itype)) then
         do ict = 1, nct
            nvar = nvar + 1
            ipnt(ict,itype) = nvar
            write(str,'(i2)') ict
            var(nvar)%label = trim(txtype(itype))//' '//'ict='//trim(adjustl(str))
            var(nvar)%min = vmin
            var(nvar)%max = vmax
            var(nvar)%nbin = nbin
         end do
      end if
   end do

   call DistFuncSample(2, nvar, var)
   call DistFuncSample(3, nvar, var)
   call DistFuncSample(4, nvar, var)

   var%nsamp2 = var%nsamp2+1

! ... sampling of itype = 1

   itype = 1
   if (ltype(itype)) then
      do ic = 1, nc
         ict = Idt(ic)
         ivar = ipnt(itype,ict)
         do j = 1, min(nAdsEvent(ic),1)
            if (AdsStart(j,ic) > tmax) cycle
            ibin = max(-1,min(floor(var(ivar)%bini*(AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end do
   end if

! ... sampling of itype = 2

   itype = 2
   if (ltype(itype)) then
      do ic = 1, nc
         ict = Idt(ic)
         ivar = ipnt(ict,itype)
         do j = 1, nAdsEvent(ic)
            if (AdsStart(j,ic) <= tdiv) then
               if (AdsStart(j,ic) > tmax) cycle
               ibin = max(-1,min(floor(var(ivar)%bini*(AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do
      end do
   end if

! ... sampling of itype = 3

   itype = 3
   if (ltype(itype)) then
      do ic = 1, nc
         ict = Idt(ic)
         ivar = ipnt(ict,itype)
         do j = 1, nAdsEvent(ic)
            if (AdsStart(j,ic) >= tdiv) then
               if (AdsStart(j,ic) > tmax) cycle
               ibin = max(-1,min(floor(var(ivar)%bini*(AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do
      end do
   end if

! ... sampling of itype = 4

   itype = 4
   if (ltype(itype)) then
      do ic = 1, nc
         ict = Idt(ic)
         ivar = ipnt(ict,itype)
         do j = 1, nAdsEvent(ic)
            if (AdsStart(j,ic) > tmax) cycle
            ibin = max(-1,min(floor(var(ivar)%bini*(AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end do
   end if

   call DistFuncNorm(1, nvar, var)

   call DistFuncSample(6, nvar, var)
   call DistFuncSample(7, nvar, var)
   write(uout,'(a,t35,f8.2)') 'tdiv (division early - late)   = ', tdiv
   write(uout,'(a,t35,f8.2)') 'tmax (largest AdsStart used)   = ', tmax
   call DistFuncHead(nvar, var, uout)
   call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

   deallocate(var, ipnt)

end subroutine AdsExam3

!************************************************************************
!> \page moluser moluser.F90
!! **AdsExam4**
!! *analysis of primary adsorption data; adsorption length*
!************************************************************************


subroutine AdsExam4

   use AdsModule
   implicit none

   integer(4), parameter :: mntype = 100
   character(40), parameter :: txroutine ='AdsExam4'
   character(80), parameter :: txheading ='adsorption length dist. function; ads. start separated'
   integer(4)     :: ntype, nbin, nvar
   real(8)        :: tlow, tupp, dt, tmin(mntype), tmax(mntype), vmin, vmax
   integer(4),   allocatable :: ipnt(:,:)
   type(df_var), allocatable :: var(:)
   integer(4)     :: itype, ivar, ibin, ict, ic, j
   character(8)   :: strlow, strupp
   character(2)   :: str

   namelist /nmlAdsExam4/ ntype, tlow, tupp, vmin, vmax, nbin

   call WriteHead(2, txheading, uout)
   call ReadPrimAdsData

   ntype = 10
   tlow = Zero
   tupp = 10.0
   vmin = Zero
   vmax = 100.0
   nbin = 100

   rewind(uin)
   read(uin,nmlAdsExam4)

   nvar = ntype*nct
   allocate(var(nvar),ipnt(nct,ntype))
   ipnt = 0

   if (ntype > mntype) call Stop(txroutine, 'ntype > mntype', uout)
   if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

   nvar = 0
   dt = (tupp-tlow)/ntype
   do itype = 1, ntype
      tmin(itype) = tlow + (itype-1)*dt
      tmax(itype) = tmin(itype) + dt
      write(strlow,'(f8.2)') tmin(itype)
      write(strupp,'(f8.2)') tmax(itype)
      do ict = 1, nct
         nvar = nvar + 1
         ipnt(ict,itype) = nvar
         write(str,'(i2)') ict
         var(nvar)%label = trim(adjustl(strlow))//'-'//trim(adjustl(strupp))//' '//'ict='//trim(adjustl(str))
         var(nvar)%min = vmin
         var(nvar)%max = vmax
         var(nvar)%nbin = nbin
      end do
   end do

   call DistFuncSample(2, nvar, var)
   call DistFuncSample(3, nvar, var)
   call DistFuncSample(4, nvar, var)

   var%nsamp2 = var%nsamp2+1

   do ic = 1, nc
      ict = Idt(ic)
      do j = 1, nAdsEvent(ic)
         do itype = 1, ntype
            if ((AdsStart(j,ic) > tmin(itype)) .and. (AdsStart(j,ic) < tmax(itype))) then
               ivar = ipnt(ict,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(AdsLength(j,ic)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do
      end do
   end do

   call DistFuncNorm(1, nvar, var)

   call DistFuncSample(6, nvar, var)
   call DistFuncSample(7, nvar, var)
   call DistFuncHead(nvar, var, uout)
   call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

   deallocate(var, ipnt)

end subroutine AdsExam4

!************************************************************************
!> \page moluser moluser.F90
!! **ReadPrimAdsData**
!! *read user-provided primary adsorption data*
!************************************************************************


subroutine ReadPrimAdsData

   use AdsModule
   implicit none

   character(40), parameter :: txroutine ='ReadPrimAdsData'
   integer(4) :: ic, j
   real(8)    :: dum

   rewind(uuser)
   read(uuser,'(5x,i8)') nct
   read(uuser,'(5x,i8)') nc

   if (.not.allocated(Idt)) then
      allocate(Idt(nc), Id(nc), nAdsEvent(nc), AdsStart(mnAdsEvent,nc), AdsLength(mnAdsEvent,nc))
      Idt = 0
      Id = 0
      nAdsEvent = 0
      AdsStart = 0.0E+00
      AdsLength = 0.0E+00
   end if
   ncct = 0

   read(uuser,*)
   read(uuser,*)
   read(uuser,*)
   do ic = 1, nc
      read(uuser,*,end=999) Idt(ic), Id(ic), nAdsEvent(ic), dum, (AdsStart(j,ic),AdsLength(j,ic), j = 1, nAdsEvent(ic))
      ncct(idt(ic)) = ncct(idt(ic)) + 1
   end do
 999 continue
   if (maxval(nadsEvent(1:nc)) > mnAdsEvent) call Stop(txroutine, 'maxval(nAdsEvent(1:nc) > mnAdsEvent', uout)

   write(uout,'(a,t40,a    )') 'primary adsorption data read from FUSER'
   write(uout,'(a,t40,i8   )') 'number of chain types          = ', nct
   write(uout,'(a,t40,i8   )') 'number of chains               = ', nc
   write(uout,'(a,t40,4i8  )') 'number of chains of diff types = ', ncct(1:nct)

end subroutine ReadPrimAdsData

!************************************************************************
!> \page moluser moluser.F90
!! **WritePrimAdsData**
!! *write user-provided primary adsorption data*
!************************************************************************


subroutine WritePrimAdsData(txheading)

   use AdsModule
   implicit none

   character(*), intent(in) :: txheading
   integer(4) :: i, j

   write(uout,*)
   write(uout,'(a)') txheading
   write(uout,*)
   write(uout,'(8(a,2x))') 'order nr', 'chain type', 'chain nr', 'ads events', 'total ads time', 'ads start', 'ads length', '...'
   write(uout,'(8(a,2x))') '--------', '----------', '--------', '----------', '--------------', '---------', '----------', '---'
   do i = 1, nc
      write(uout,'(4(i8,2x),5x,400f9.3)') &
      i, Idt(i), Id(i), nAdsEvent(i), sum(AdsLength(1:nAdsEvent(i),i)), (AdsStart(j,i), AdsLength(j,i), j = 1, max(0,nAdsEvent(i)))
   end do

end subroutine WritePrimAdsData

!************************************************************************
!> \page moluser moluser.F90
!! **Z_DF_Slit**
!! *distribution function based on z-coordinates, planar geometry*
!************************************************************************


subroutine Z_DF_Slit(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Z_DF_Slit'
   character(80), parameter :: txheading ='1d z-distribution function, one particle type'
   integer(4)   , parameter :: ntype = 3
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:)
   integer(4),    save :: iptz1,iptz2

   integer(4) :: itype, ivar, ibin, i, j

   select case (iStage)
   case (iReadInput)

      vtype%l =.true.
      vtype%min =-boxlen2(3)
      vtype%max =+boxlen2(3)
      vtype%nbin = 200
      iptz1 = 1
      iptz2 = 2

      if (txbc /= 'xy') call Stop (txroutine,'txbc /= ''xy''',uout)
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['atom of p1', 'atom of p2','a1 of p1  ']
      vtype%nvar = 1

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ivar + 1
            ipnt(itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

   ! ... summation of type 1

      itype = 1
      if(vtype(itype)%l) then
         ivar = ipnt(itype)
         do i = ipnpt(iptz1),ipnpt(iptz1)+nppt(iptz1)-1
            do j = ianpn(i),ianpn(i)+1
               ibin = max(-1,min(floor(var(ivar)%bini*(r(3,j)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end do
         end do
      end if

   ! ... summation of type 2

      itype = 2
      if(vtype(itype)%l) then
         ivar = ipnt(itype)
         do i=ipnpt(iptz2),ipnpt(iptz2)+nppt(iptz2)-1
            do j=ianpn(i),ianpn(i)+1
               ibin = max(-1,min(floor(var(ivar)%bini*(r(3,j)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end do
         end do
      end if

   ! ... summation of type 3

      itype = 3
      if(vtype(itype)%l) then
         ivar = ipnt(itype)
         do i=1,nppt(iptz1)
            ibin = max(-1,min(floor(var(ivar)%bini*(r(3,ianpn(i))-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
         end do
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,3x,i5)') 'particle type 1           = ', iptz1
      write(uout,'(a,3x,i5)') 'particle type 2           = ', iptz2
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var,ipnt)

   end select

end subroutine Z_DF_Slit

!************************************************************************
!> \page moluser moluser.F90
!! **ElMom**
!! *electrostatic moment of a particle*
!************************************************************************


subroutine ElMom(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ElMom'
   character(80), parameter :: txheading ='Electrostatic moments of a particle'
   integer(4), parameter :: ntype = 1       ! maximum number of properties
   integer(4)      , save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4) :: nrot
   real(8) :: qr(3), qrr(3,3), diagonal(3), eivr(3,3)

!  namelist /nmlElMom/ xxxx

   if (slave) return

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)
!     rewind(uin)
!     read(uin,nmlScalarDemo)

   case (iWriteInput)

! ... set nvar and allocate memory

      nvar = 5
      allocate(var(nvar))

! ... set label

      var(1)%label = 'net charge                      = '
      var(2)%label = 'dipole moment (com)             = '
      var(3)%label = 'quadrupole moment (com, min)    = '
      var(4)%label = 'quadrupole moment (com, mid)    = '
      var(5)%label = 'quadrupole moment (com, max)    = '
      var%norm = One

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
      call FileOpen(uuser, fuser , 'form/noread')
      write(uuser,'(a)') 'net charge,  dip mom, theta,   phi,  qq(min), qq(mid),qq(max)'

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

!      do ia = 1, na
!         write(*,'(l5,f8.3,2x,3f8.3)') laz(ia), az(ia), (r(m,ia)-ro(m,1), m = 1,3)
!      end do
      var(1)%value = sum(az,1,laz)
      qr(1) = sum(az(:)*(r(1,:)-ro(1,1)), 1, laz(:))
      qr(2) = sum(az(:)*(r(2,:)-ro(2,1)), 1, laz(:))
      qr(3) = sum(az(:)*(r(3,:)-ro(3,1)), 1, laz(:))
!      call writevec(3,qr,'qr',1,3,1,6)
!      call CarToSph('deg',qr(1),qr(2),qr(3),rr,theta,phi)
      var(2)%value = sqrt(sum(qr(:)**2))
      qrr(1,1) = sum(az(:)*(r(1,:)-ro(1,1))**2, 1, laz(:))
      qrr(2,2) = sum(az(:)*(r(2,:)-ro(2,1))**2, 1, laz(:))
      qrr(3,3) = sum(az(:)*(r(3,:)-ro(3,1))**2, 1, laz(:))
      qrr(1,2) = sum(az(:)*(r(1,:)-ro(1,1))*(r(2,:)-ro(2,1)), 1, laz(:))
      qrr(1,3) = sum(az(:)*(r(1,:)-ro(1,1))*(r(3,:)-ro(3,1)), 1, laz(:))
      qrr(2,3) = sum(az(:)*(r(2,:)-ro(2,1))*(r(3,:)-ro(3,1)), 1, laz(:))
      qrr(2,1) = qrr(1,2)
      qrr(3,1) = qrr(1,3)
      qrr(3,2) = qrr(1,2)
!      call writemat(3,3,qrr,'qrr',1,3,1,1,3,1,6)
      call Diag(3, qrr, diagonal, eivr, nrot)
!      call writevec(3,diagonal,'diag',1,3,1,6)
      diagonal=diagonal-third*sum(diagonal)
!      call writevec(3,diagonal,'diag',1,3,1,6)
!      call writemat(3,3,eivr,'eivr',1,3,1,1,3,1,6)
      var(3)%value = minval(diagonal(:))
      var(5)%value = maxval(diagonal(:))
      var(4)%value = sum(diagonal(:))-(var(3)%value+var(5)%value)
!      write(*,'(a,10f10.3)') 'var(1:nvar)', var(1:nvar)%value
!      if(mod(istep1,10) == 0)      write(uuser,'(f8.0,2x,3f8.2,2x,3f8.2)') var(1)%value,rr,theta,phi, var(3:5)%value
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

   end select

end subroutine ElMom


!************************************************************************
!> \page moluser moluser.F90
!! **XYProjectDF**
!! *calculate on 'z' = 0 projected normalized density df*
!************************************************************************

!     here adapted for pure benezene fluid (cf skipper JACS 132, 5735, 2010, fig7)

subroutine XYProjectDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='XYProjectDF'
   character(80), parameter :: txheading ='on ''z'' = 0 projected normalized density df'
   integer(4)   , parameter :: ntype = 1
   type(static2D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df2d_var),allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),                    save :: rcutdist2
   real(8),                    save :: zlim = 2*Half
   integer(4) :: itype, ivar, ibin1, ibin2, ip, jp
   real(8)    :: dx, dy, dz, r2, norm, ddvol, xp, yp, zp

   namelist /nmlXYProjectDF/ vtype

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.false.
      vtype(1)%min(1:2) = [ -5.0, -5.0 ]
      vtype(1)%max(1:2) = [ +5.0, +5.0 ]
      vtype(1)%nbin(1:2) = [ 20, 20 ]

      rewind(uin)
      read(uin,nmlXYProjectDF)

      rcutdist2 = 3*vtype(1)%min(1)**2

! ... check condition

      if (maxval(vtype%nbin(1)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(1) > mnbin_df2d', uout)
      if (maxval(vtype%nbin(2)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(2) > mnbin_df2d', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['in-plane']
      vtype%nvar = 1

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(1,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ivar+1
            ipnt(1,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end if
      end do

      call DistFunc2DSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ipnt(1,itype)
            do ip = 1, np
               do jp = 1, np
                  if (ip == jp) cycle
                  dx = ro(1,ip)-ro(1,jp)
                  dy = ro(2,ip)-ro(2,jp)
                  dz = ro(3,ip)-ro(3,jp)
                  call PBCr2(dx,dy,dz,r2)
                  if (r2 < rcutdist2) then
                     xp = dx*ori(1,1,jp)+dy*ori(2,1,jp)+dz*ori(3,1,jp)  ! dr's projection on x'-axis
                     yp = dx*ori(1,2,jp)+dy*ori(2,2,jp)+dz*ori(3,2,jp)  ! dr's projection on y'-axis
                     zp = dx*ori(1,3,jp)+dy*ori(2,3,jp)+dz*ori(3,3,jp)  ! dr's projection on y'-axis
                     if (abs(zp) < zlim) then                           ! select neighbours in the z' = 0 plane
                        ibin1 = max(-1,min(floor(var(ivar)%bini(1)*(xp-var(ivar)%min(1))),int(var(ivar)%nbin(1))))
                        ibin2 = max(-1,min(floor(var(ivar)%bini(2)*(yp-var(ivar)%min(2))),int(var(ivar)%nbin(2))))
                        var(ivar)%avs2(ibin1,ibin2) = var(ivar)%avs2(ibin1,ibin2)+One
                     end if
                  end if
               end do
            end do
         end if
      end do

   case (iAfterMacrostep)

! ... normalization

      do itype = 1, ntype
         ivar = ipnt(1,itype)
         ddvol = (var(ivar)%max(1)-var(ivar)%min(1))*(var(ivar)%max(2)-var(ivar)%min(2))*(Two*zlim)/(var(ivar)%nbin(1)*var(ivar)%nbin(2))
         norm = (vol/ddvol)/np**2
         var(ivar)%avs2(-1:var(ivar)%nbin(1),-1:var(ivar)%nbin(2)) = var(ivar)%avs2(-1:var(ivar)%nbin(1),-1:var(ivar)%nbin(2))*norm
      end do
      call DistFunc2DSample(iStage, nvar, var)

      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,f10.3)') 'rcutdist2                      = ', rcutdist2
      write(uout,'(a,t35,f10.3)') 'zlim                           = ', zlim
      call DistFunc2DHead(nvar, var, uout)
      call DistFunc2DShow(1, txheading, nvar, var, uout)
      call DistFunc2DList(1, txheading, nvar, var, ulist)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine XYProjectDF

!************************************************************************
!> \page moluser moluser.F90
!! **SPDF_COMB**
!! *calculate single particle distribution functions for comb polymers (coarse_grained)*
!************************************************************************

!     parallel version (master only) By Daniel A.
!     works for PBC

!     ltype(1)  COM of each comb copolymer
!     ltype(2)  COM of PECS (comb  + linear)

subroutine SPDF_COMB(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   integer(4)   , parameter :: ntype = 3
   character(40), parameter :: txroutine ='SPDF_COMB'
   character(80), parameter :: txheading ='single particle distribution functions for comb polymer COM whole comb polymer'
   character(5),  parameter :: txtype(ntype) = ['rrONE','rrTWO','rrTHR']

   logical,       save :: ltype(ntype)
   real(8),       save :: vmin(ntype), vmax(ntype)
   integer(4),    save :: nbin, nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),   allocatable, save :: ipnt(:,:)
   real(8),      allocatable, save :: rotemp(:,:), ro_loc(:,:), ro_temp(:,:)

   integer(4) :: ivar, ibin, itype
   integer(4) :: ip, igr, ic, ict, n_ic(nh), iseg, ih, igen, n_chain
   real(8)    :: r1, r12, r2, r22, norm, norm1, norm2,vsum, vsum1, vsum2, dvol
   real(8)    :: xcom1, ycom1, zcom1, xcom(nh), ycom(nh), zcom(nh)
   real(8)    ::  dinf, dx, dy, dz, r_dist
   integer(4) ::  jct, jc, jseg, jp

   namelist /nmlSPDF_COMB/ ltype, vmin, vmax, nbin

   if(slave) return

   if(ltrace) call WriteTrace(2,txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ltype(1:ntype) =.false.
      vmin(1) = Zero
      vmax(1) = 100.0
      vmin(2) = Zero
      vmax(2) = 100.0
      vmin(3) = Zero
      vmax(3) = 100.0
      nbin    = 100

      rewind(uin)
      read(uin,nmlSPDF_COMB)

      if(nbin > mnbin_df) call Stop('txroutine', 'nbin > mnbin_df', uout)

   case (iWriteInput)

     if(.not.lclink) call stop(txroutine, '.not.lclink', uout)

     if (.not.allocated(rotemp)) then
        allocate(rotemp(3,np), ro_loc(3,np), ro_temp(3,np))
        rotemp = 0.0E+00
        ro_loc = 0.0E+00
        ro_temp = 0.0E+00
     end if

! ... set nvar as well as allocate memory

      nvar  = ngr(1)*count(ltype(1:ntype))
      allocate(var(nvar), ipnt(ngrgr,ntype))
      ipnt = 0

      nvar = 0
    do itype = 1, ntype
       if(ltype(itype)) then
          if(itype == 1) then
             do igr = 1, ngr(1)
                nvar = nvar+1
                ipnt(igr,itype) = nvar
                var(nvar)%label = trim(txtype(itype))//' '//txgr(igr)
                var(nvar)%min = vmin(itype)
                var(nvar)%max = vmax(itype)
                var(nvar)%nbin = nbin
             end do
          else if(itype == 2) then
             do igr = 1, ngr(1)
                nvar = nvar+1
                ipnt(igr,itype) = nvar
                var(nvar)%label = trim(txtype(itype))//' '//txgr(igr)
                var(nvar)%min = vmin(itype)
                var(nvar)%max = vmax(itype)
                var(nvar)%nbin = nbin
             end do
          else if(itype == 3) then
             do igr = 1, ngr(1)
                nvar = nvar+1
                ipnt(igr,itype) = nvar
                var(nvar)%label = trim(txtype(itype))//' '//txgr(igr)
                var(nvar)%min = vmin(itype)
                var(nvar)%max = vmax(itype)
                var(nvar)%nbin = nbin
             end do
          end if
       end if
    end do

    call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if(lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

 !! ....    calculate COM
!    ... COM_all   it includes backbone and strands

     xcom(1:nh) = Zero
     ycom(1:nh) = Zero
     zcom(1:nh) = Zero
     n_ic(1:nh) = 0
     vaux(1:3,1:np) = 0.0
     ro_loc(1:3,1:np) = 0.0

! ... undo hierarchical structures

   do ih = 1, nh                              ! loop over number of hierarchic structures
      do igen = 0, ngen                                                    ! loop over generations
         ict = ictgen(igen)                                                ! chain type
         do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1       ! loop over chains of the structure
            do iseg = 1, npct(ict)                                         ! loop over segments
               ip = ipnsegcn(iseg,ic)
               if (ict ==1) then                                           ! UNDO backbone (first chain type)
                  call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                  ro_temp(1:3,1:ip) = vaux(1:3,1:ip)                       ! store undo backbone
              end if
              if (nbondcl(ip) ==1.AND.ict > 1 ) call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)  ! UNDO grafted chains, taking  UNDO positions
            !linked beads belonging to the backbone as reference
            end do
         end do
      end do
   end do
! ...  calculate COM of each comb copolymer
   do ih = 1, nh
      do igen = 0, ngen                                                      ! loop over generations
        ict = ictgen(igen)                                                   ! chain type
        do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
           do iseg = 1, npct(ict)                                            ! loop over segments
              ip = ipnsegcn(iseg,ic)
              ro_loc(1:3,ip) = vaux(1:3,ip)
              xcom(ih) = xcom(ih) + ro_loc(1,ip)
              ycom(ih) = ycom(ih) + ro_loc(2,ip)
              zcom(ih) = zcom(ih) + ro_loc(3,ip)
              n_ic(ih) = n_ic(ih) + 1
           end do
        end do
      end do
      xcom(ih) = xcom(ih)/n_ic(ih)
      ycom(ih) = ycom(ih)/n_ic(ih)
      zcom(ih) = zcom(ih)/n_ic(ih)
   end do

!  undo non-hierarchical (linear) chains taking as reference the backbone ic = 1
    do ic = 1, nc
       ict = ictcn(ic)
       ip = ipnsegcn(1,ic)
       if(ihnpn(ip) /= 0) cycle           ! exclude hierarchical structures
       call UndoPBCChain(ro_temp(1,ipnsegcn(1,1)), ic, -1, vaux)
       do iseg = 1, npct(ict)             ! loop over segments
          ip = ipnsegcn(iseg,ic)
          ro_loc(1:3,ip) = vaux(1:3,ip)
       end do
    end do
    dinf = boxlen(1)
    ict =1
    ic = icnct(ict)
    do iseg = 1, npct(ict)
       ip = ipnsegcn(iseg,ic)
       jct =3
       jc = icnct(jct)
       do jseg = 1, npct(jct)
          jp = ipnsegcn(jseg,jc)
          dx = ro_loc(1,ip) - ro_loc(1,jp)
          dy = ro_loc(2,ip) - ro_loc(2,jp)
          dz = ro_loc(3,ip) - ro_loc(3,jp)
          r_dist = dx**2+dy**2+dz**2
          if (sqrt(r_dist) < dinf) dinf = sqrt(r_dist)
       end do
    end do
    if (dinf > boxlen(1)/3.0)  then
       do ic = 1, nc
          ict = ictcn(ic)
          ip = ipnsegcn(1,ic)
          if(ihnpn(ip) /= 0) cycle       ! exclude hierarchical structures
          call UndoPBCChain(ro_temp(1,ipnsegcn(1,1)), ic, 1, vaux)
          do iseg = 1, npct(ict)         ! loop over segments
             ip = ipnsegcn(iseg,ic)
             ro_loc(1:3,ip) = vaux(1:3,ip)
          end do
       end do
   end if

! ...  calculate COM of  combed PECS (ALL comb copolymers plus linear chain)

    xcom1 = 0.0
    ycom1 = 0.0
    zcom1 = 0.0
    n_chain= 0
    do ic = 1, nc                      ! loop over ALL CHAINS (they all made one PECS)
       ict = ictcn(ic)
       ip = ipnsegcn(1,ic)
       do iseg = 1, npct(ict)          ! loop over segments
          ip = ipnsegcn(iseg,ic)
          ro_loc(1:3,ip) = vaux(1:3,ip)
          xcom1 = xcom1 + ro_loc(1,ip)
          ycom1 = ycom1 + ro_loc(2,ip)
          zcom1 = zcom1 + ro_loc(3,ip)
          n_chain = n_chain+1
       end do
    end do

    xcom1 = xcom1/n_chain
    ycom1 = ycom1/n_chain
    zcom1 = zcom1/n_chain

     do ip = 1,np
        igr = igrpn(ip,1)
        if(igr <= 0 ) cycle
        if(ltype(1)) then       !! ... spdf  (COM_comb) of each branched polymer
                ivar = ipnt(igr,1)
                                ih = ihnpn(ip)
                                if (ih /= 0) then
                                        rotemp(1,ip)=ro_loc(1,ip)-xcom(ih)  !   spdf for ih hierarchical structures
                                        rotemp(2,ip)=ro_loc(2,ip)-ycom(ih)
                                        rotemp(3,ip)=ro_loc(3,ip)-zcom(ih)
                                else
                                        rotemp(1,ip)=ro_loc(1,ip)-xcom(1)  !   spdf for the linear chain complexed with one hierarchical structure
                                        rotemp(2,ip)=ro_loc(2,ip)-ycom(1)
                                        rotemp(3,ip)=ro_loc(3,ip)-zcom(1)
                                end if
                                r12 = rotemp(1,ip)**2+rotemp(2,ip)**2+rotemp(3,ip)**2
                                r1 = sqrt(r12)
                        ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                        var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
        if(ltype(2)) then        !! ... spdf  for whole combed PECS
                        ivar = ipnt(igr,2)
                        rotemp(1,ip)=ro_loc(1,ip)-xcom1
                        rotemp(2,ip)=ro_loc(2,ip)-ycom1
                        rotemp(3,ip)=ro_loc(3,ip)-zcom1
                        r22 = rotemp(1,ip)**2+rotemp(2,ip)**2+rotemp(3,ip)**2
                        r2 = sqrt(r22)
                        ibin = max(-1,min(floor(var(ivar)%bini*(r2-var(ivar)%min)),int(var(ivar)%nbin)))
                        var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
        if(ltype(3)) then
            end if
     end do

   case (iAfterMacrostep)

   do igr = 1,ngr(1)
         if(ltype(1)) then
            ivar = ipnt(igr,1)
            vsum = sum(var(ivar)%avs2(0:var(ivar)%nbin-1))
            norm = Zero
            if(vsum/= Zero) norm = var(ivar)%nsamp2*(var(ivar)%max**3-var(ivar)%min**3)/vsum
            do ibin = 0, nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
         if(ltype(2)) then
            ivar = ipnt(igr,2)
            vsum1 = sum(var(ivar)%avs2(0:var(ivar)%nbin-1))
            norm1 = Zero
            if(vsum1/= Zero) norm1 = var(ivar)%nsamp2*(var(ivar)%max**3-var(ivar)%min**3)/vsum1
            do ibin = 0, nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm1/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
         if(ltype(3)) then
            ivar = ipnt(igr,3)
            vsum2 = sum(var(ivar)%avs2(0:var(ivar)%nbin-1))
            norm2 = Zero
            if(vsum2/= Zero) norm2 = var(ivar)%nsamp2*(var(ivar)%max**3-var(ivar)%min**3)/vsum2
            do ibin = 0, nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm2/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
        end do

      call DistFuncSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, ipnt,rotemp, ro_loc, ro_temp)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine SPDF_COMB


!************************************************************************
!> \page moluser moluser.F90
!! **COMB_DF**
!! *calculate  type distribution functions of comb polymer type (coarse model)*
!************************************************************************

!     works for PBC
!     By Daniel A

!     type  label  quantity
!     ----  -----  --------
!     1     rg_comb     radius of gyration of the comb polymer(include backbone and strands), COM of comb particle; includes all hierarchical structures
!     2     asp_comb    asphericty of the comb polymer(include backbone and strands), COM of comb particle; includes all hierarchical structures
!    3-5   rg_comb(ih)  radius of gyration of up to three comb polymers(ih >1)
!    6-8   asp_comb(ih) asphericity of up to three comb polymers(ih >1)

subroutine COMB_DF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='COMB_DF'
   character(60), parameter :: txheading = 'distribution functions of the radius of gyration and asphericity of comb copolymer'
   integer(4), parameter :: ntype = 8
   character(5),  parameter :: txtype(ntype) = [ 'rg_co','as_co','rg_c1','rg_c2','rg_c3','as_c1','as_c2','as_c3']

   logical,       save :: ltype(ntype)
   real(8),       save :: vmin(ntype), vmax(ntype)
   integer(4),    save :: nbin, nvar

   type(df_var),  allocatable, save :: var(:)
   integer(4),   allocatable, save :: ipnt(:,:)
   real(8),      allocatable, save :: ro_loc(:,:), ro_temp(:,:)

   integer(4)     :: ic, ip, ict, itype, ivar, ibin, iseg, nrot, ih, n_ic(nh), igen
   real(8)        :: properties(13), value, xcom(nh), ycom(nh), zcom(nh), dx(nh), dy(nh), dz(nh), r2, vsumr(nh)
   real(8)        :: l2_small2, l2_mid2, l2_large2, mimat(3,3,nh), mimat_loc(3,3), dummat_loc(3,3), diagonal(3)

   namelist /nmlCOMB_DF/ ltype, vmin, vmax, nbin

   if(ltrace) call WriteTrace(2,txroutine, iStage)

   if(slave) return

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

    if(.not.lclink) call stop(txroutine, '.not.lclink', uout)

      ltype(1:ntype) =.false.
      vmin(1) =   0.0
      vmax(1) = 100.0
      vmin(2) =   0.0
      vmax(2) =   1.0
      vmin(3) =   0.0
      vmax(3) = 100.0
      vmin(4) =   0.0
      vmax(4) = 100.0
      vmin(5) =   0.0
      vmax(5) = 100.0
      vmin(6) =   0.0
      vmax(6) =   1.0
      vmin(7) =   0.0
      vmax(7) =   1.0
      vmin(8) =   0.0
      vmax(8) =   1.0

      nbin  =  100
      rewind(uin)
      read(uin,nmlCOMB_DF)

      if(nbin > mnbin_df) call Stop('txroutine', 'nbin > mnbin_df', uout)

   case (iWriteInput)

          if(.not.lclink) call stop(txroutine, '.not.lclink', uout)

      if (.not.allocated(ro_loc)) then
         allocate(ro_loc(3,np), ro_temp(3,np))
         ro_loc = 0.0E+00
         ro_temp = 0.0E+00
      end if

      nvar = nct*count(ltype(1:ntype))
      allocate(var(nvar), ipnt(nct,ntype))
      ipnt = 0

      nvar = 0
      do itype = 1, ntype
         if(ltype(itype)) then
            do ict = 1, nct
               nvar = nvar+1
               ipnt(ict,itype) = nvar
               var(nvar)%label = trim(txtype(itype))//' '//txct(ict)
               var(nvar)%min = vmin(itype)
               var(nvar)%max = vmax(itype)
               var(nvar)%nbin = nbin
            end do
         end if
      end do

    call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

     call DistFuncSample(iStage, nvar, var)
     if(lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

         xcom(1:nh) = Zero
     ycom(1:nh) = Zero
     zcom(1:nh) = Zero
     n_ic(1:nh) = 0
         vaux(1:3,1:np) = 0.0
         ro_loc(1:3,np) = 0.0

!...  undo conditions for hierarchical structure

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
   do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
             if (ict ==1) then                                                 ! UNDO backbone (first chain type)
                                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                                ro_temp(1:3,1:na) = vaux(1:3,1:na)                         ! store undo backbone
                         end if
                         if (nbondcl(ip) ==1.AND.ict > 1 )     call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)      ! UNDO grafted chains, taking  UNDO positions
                                                                                                                                                          !linked beads belonging to the backbone as reference
         end do
      end do
   end do
   end do

 !    ... COM_comb polymer, xcom(ih), ycom(ih), zcom(ih)

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
                        xcom(ih) = xcom(ih) + ro_loc(1,ip)
                        ycom(ih) = ycom(ih) + ro_loc(2,ip)
                        zcom(ih) = zcom(ih) + ro_loc(3,ip)
                        n_ic(ih) =    n_ic(ih) + 1
         end do
      end do
    end do
        xcom(ih) = xcom(ih)/n_ic(ih)
    ycom(ih) = ycom(ih)/n_ic(ih)
    zcom(ih) = zcom(ih)/n_ic(ih)
   end do


 ! ... properties(1)    radius of gyration

   vsumr(1:nh) = Zero
   mimat(1:3,1:3,1:nh) = Zero

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
   do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
                        dx(ih) = ro_loc(1,ip) -xcom(ih)
                        dy(ih) = ro_loc(2,ip) -ycom(ih)
                        dz(ih) = ro_loc(3,ip) -zcom(ih)
                        r2 = dx(ih)**2+dy(ih)**2+dz(ih)**2
                        vsumr(ih) = vsumr(ih)+r2
                         mimat(1,1,ih) = mimat(1,1,ih)+dx(ih)**2        !        +dy1**2+dz1**2
             mimat(1,2,ih) = mimat(1,2,ih)+dx(ih)*dy(ih)    !        -dx1*dy1
             mimat(1,3,ih) = mimat(1,3,ih)+dx(ih)*dz(ih)    !        -dx1*dz1
             mimat(2,2,ih) = mimat(2,2,ih)+dy(ih)**2        !        +dx1**2+dz1**2
             mimat(2,3,ih) = mimat(2,3,ih)+dy(ih)*dz(ih)    !        -dy1*dz1
             mimat(3,3,ih) = mimat(3,3,ih)+dz(ih)**2        !        +dx1**2+dy1**2
         end do
      end do
          end do
         vsumr(ih) = sqrt(vsumr(ih)/n_ic(ih))
   end do

        properties(1) = sum(vsumr(1:nh))/nh

!... properties(3-5)    radius of gyration of each comb polymer
    if (nh >1 .AND. nh <=3) properties(3:5) = vsumr(1:3)

! ... properties(2) asphericity, comb polymer, COM comb polymer

        properties(2) = Zero

    do ih = 1, nh                                                            ! loop over number of hierarchic structures
                mimat(2,1,ih) = mimat(1,2,ih)
                mimat(3,1,ih) = mimat(1,3,ih)
                mimat(3,2,ih) = mimat(2,3,ih)
                mimat_loc(1:3,1:3) = mimat(1:3,1:3,ih)
                call Diag(3, mimat_loc, diagonal, dummat_loc, nrot)
                diagonal(1) = diagonal(1)/n_ic(ih)
                diagonal(2) = diagonal(2)/n_ic(ih)
                diagonal(3) = diagonal(3)/n_ic(ih)
                l2_small2 = min(diagonal(1),diagonal(2),diagonal(3))
                l2_large2 = max(diagonal(1),diagonal(2),diagonal(3))
                l2_mid2  = diagonal(1)+diagonal(2)+diagonal(3)-(l2_small2+l2_large2)
                properties(2) = properties(2) +((l2_small2-l2_mid2)**2+(l2_small2-l2_large2)**2+(l2_large2-l2_mid2)**2)/(2.0d0*(l2_small2+l2_mid2+l2_large2)**2  )
    if (nh >1 .AND. nh <=3)  properties(5+ih) = ((l2_small2-l2_mid2)**2+(l2_small2-l2_large2)**2+(l2_large2-l2_mid2)**2)/(2.0d0*(l2_small2+l2_mid2+l2_large2)**2  )
        end do

        properties(2) = properties(2)/nh


! ... summation of type 1 to type 8

         do itype = 1, 8
            if(ltype(itype)) then
               ivar = ipnt(1,itype)
               if(itype == 1) then
                  value = properties(1)
               else if(itype == 2) then
                  value = properties(2)
               else if(itype == 3) then
                  value = properties(3)
               else if(itype == 4) then
                  value = properties(4)
               else if(itype == 5) then
                  value = properties(5)
               else if(itype == 6) then
                  value = properties(6)
               else if(itype == 7) then
                  value = properties(7)
               else if(itype == 8) then
                  value = properties(8)
               end if
               ibin = max(-1,min(int(var(ivar)%bini*(value-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
         end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

          deallocate(var, ipnt, ro_loc, ro_temp)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine COMB_DF

!***********************************************************************
!> \page moluser moluser.F90
!! **COMBAver**
!! *calculate averages of comb chain quantities*
!************************************************************************

!     By Daniel A.

!     type  label          quantity
!     ----  -----          --------
!     1     <r(g)**2>**0.5
!     2     <A>
!     3                    hydrodynamic radius
!     4-6   <r(g)**2>**0.5 for up to three separate comb chains
!     7-9   <A>                  for up to three separate comb chains

subroutine COMBAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine =''
   character(80), parameter :: txheading ='comb structure quantities'
   integer(4), parameter :: ntype = 9               ! maximum number of properties
   integer(4)      , save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   real(8)         , allocatable, save :: mimat(:,:,:), ro_temp(:,:), ro_loc(:,:)
   integer(4)     :: ic, ict, nrot, ih, n_ic(nh), jct, igen, jgen, iseg, jseg, jc, ip, im
   real(8)        :: xcom(nh), ycom(nh), zcom(nh), dx(nh), dy(nh), dz(nh), r2, vsumr(nh), vsum_rh(nh)
   real(8)        :: mimat_loc(3,3), dummat_loc(3,3), l2_small2, l2_mid2, l2_large2, diagonal(3)


   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2,txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = ntype
      allocate(var(nvar),mimat(3,3,nh),ro_temp(3,na),ro_loc(3,np))
      mimat = 0.0E+00
      ro_temp = 0.0E+00
      ro_loc = 0.0E+00

      var(1)%label = '<r(g)**2>**0.5                 = ' ! rms radius of gyration
      var(1)%norm = One
          var(2)%label = '<Asphericity>                  = ' ! asphericity
      var(2)%norm = One
          var(3)%label = '<Hydrodynamic radius>          = ' ! hydrodynamic radius (Arun Yethiraj, J. Chem. Phys. 125, 204901, 2006)
      var(3)%norm = One
          var(4)%label = '<r(g)**2>**0.5  1st comb       = ' ! rms radius of gyration first comb
      var(4)%norm = One
          var(5)%label = '<r(g)**2>**0.5  2nd comb       = ' ! rms radius of gyration first comb
      var(5)%norm = One
          var(6)%label = '<r(g)**2>**0.5  3rd comb       = ' ! rms radius of gyration first comb
      var(6)%norm = One
          var(7)%label = '<Asphericity>   1st comb       = ' ! asphericity
      var(7)%norm = One
          var(8)%label = '<Asphericity>   2nd comb       = ' ! asphericity
      var(8)%norm = One
          var(9)%label = '<Asphericity>   3rd comb       = ' ! asphericity
      var(9)%norm = One

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)
!..............................................

         xcom(1:nh) = Zero
     ycom(1:nh) = Zero
     zcom(1:nh) = Zero
     n_ic(1:nh) = 0
         vaux(1:3,1:np) = 0.0
         ro_loc(1:3,1:np) = 0.0

!...  undo conditions for hierarchical structure

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
            if (ict ==1) then                                                 ! UNDO backbone (first chain type)
                                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                                ro_temp(1:3,1:na) = vaux(1:3,1:na)                         ! store undo backbone
                        end if
                        if (nbondcl(ip) ==1.AND.ict > 1 ) call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)         ! UNDO grafted chains, taking  UNDO positions
                                                                                                                                                                                                !linked beads belonging to the backbone as reference
         end do
      end do
    end do
   end do

 !    ... COM_comb polymer, xcom(ih), ycom(ih), zcom(ih)

    do ih = 1, nh                                                            ! loop over number of hierarchic structures
                do igen = 0, ngen                                                       ! loop over generations
                        ict = ictgen(igen)                                                   ! chain type
                        do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
                                do iseg = 1, npct(ict)                                            ! loop over segments
                                        ip = ipnsegcn(iseg,ic)
                                        ro_loc(1:3,ip) = vaux(1:3,ip)
                                        xcom(ih) = xcom(ih) + ro_loc(1,ip)
                                        ycom(ih) = ycom(ih) + ro_loc(2,ip)
                                        zcom(ih) = zcom(ih) + ro_loc(3,ip)
                                        n_ic(ih) =    n_ic(ih) + 1
                                end do
                        end do
                end do
                xcom(ih) = xcom(ih)/n_ic(ih)
        ycom(ih) = ycom(ih)/n_ic(ih)
        zcom(ih) = zcom(ih)/n_ic(ih)
        end do

 ! ... var(1)%value =    radius of gyration

   vsumr(1:nh) = Zero
   mimat(1:3,1:3,1:nh) = Zero

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
                ict = ictgen(igen)                                                   ! chain type
                do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
                        do iseg = 1, npct(ict)                                            ! loop over segments
                                ip = ipnsegcn(iseg,ic)
                                ro_loc(1:3,ip) = vaux(1:3,ip)
                                dx(ih) = ro_loc(1,ip) -xcom(ih)
                                dy(ih) = ro_loc(2,ip) -ycom(ih)
                                dz(ih) = ro_loc(3,ip) -zcom(ih)
                                r2 = dx(ih)**2+dy(ih)**2+dz(ih)**2
                                vsumr(ih) = vsumr(ih)+r2
                                mimat(1,1,ih) = mimat(1,1,ih)+dx(ih)**2        !        +dy1**2+dz1**2
                                mimat(1,2,ih) = mimat(1,2,ih)+dx(ih)*dy(ih)    !        -dx1*dy1
                                mimat(1,3,ih) = mimat(1,3,ih)+dx(ih)*dz(ih)    !        -dx1*dz1
                                mimat(2,2,ih) = mimat(2,2,ih)+dy(ih)**2        !        +dx1**2+dz1**2
                                mimat(2,3,ih) = mimat(2,3,ih)+dy(ih)*dz(ih)    !        -dy1*dz1
                                mimat(3,3,ih) = mimat(3,3,ih)+dz(ih)**2        !        +dx1**2+dy1**2
                        end do
                end do
        end do
        vsumr(ih) = sqrt(vsumr(ih)/n_ic(ih))
   end do

        var(1)%value = sum(vsumr(1:nh))/nh

 ! ... var(4:6)%value =    radius of gyration, 1st, 2nd and 3rd comb
        if (nh >1 .AND. nh <=3)  var(4:6)%value = vsumr(1:3)

! ... var(2)%value = asphericity, comb polymer, COM comb polymer

        var(2)%value = Zero

    do ih = 1, nh                                                            ! loop over number of hierarchic structures
                mimat(2,1,ih) = mimat(1,2,ih)
                mimat(3,1,ih) = mimat(1,3,ih)
                mimat(3,2,ih) = mimat(2,3,ih)
                mimat_loc(1:3,1:3) = mimat(1:3,1:3,ih)
                call Diag(3, mimat_loc, diagonal, dummat_loc, nrot)
                diagonal(1) = diagonal(1)/n_ic(ih)
                diagonal(2) = diagonal(2)/n_ic(ih)
                diagonal(3) = diagonal(3)/n_ic(ih)
                l2_small2 = min(diagonal(1),diagonal(2),diagonal(3))
                l2_large2 = max(diagonal(1),diagonal(2),diagonal(3))
                l2_mid2  = diagonal(1)+diagonal(2)+diagonal(3)-(l2_small2+l2_large2)
                var(2)%value = var(2)%value +((l2_small2-l2_mid2)**2+(l2_small2-l2_large2)**2+(l2_large2-l2_mid2)**2)/(2.0d0*(l2_small2+l2_mid2+l2_large2)**2  )
! ... var(7:9)%value =    radius of gyration, 1st, 2nd and 3rd comb
                if (nh >1 .AND. nh <=3) var(6+ih)%value = ((l2_small2-l2_mid2)**2+(l2_small2-l2_large2)**2+(l2_large2-l2_mid2)**2)/(2.0d0*(l2_small2+l2_mid2+l2_large2)**2)
        end do

        var(2)%value = var(2)%value/nh

! ... var(3)%value = hydrodynamic radius, comb polymer, COM comb polymer

   vsum_rh(1:nh) = Zero
   do ih = 1, nh                                                            ! loop over number of hierarchic structures
        do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
                        do jgen = 0, ngen
                                jct = ictgen(jgen)
                                do jc = icihigen(ih,jgen), icihigen(ih,jgen) + nch(jgen) -1          ! loop over chains of the structure
                                        do jseg = 1, npct(ict)
                                                im = ipnsegcn(jseg,jc)
                                                if(ip >= im) cycle
                                                dx(ih) = ro_loc(1,ip)-ro_loc(1,im)
                                                dy(ih) = ro_loc(2,ip)-ro_loc(2,im)
                                                dz(ih) = ro_loc(3,ip)-ro_loc(3,im)
                                                r2 = dx(ih)**2+dy(ih)**2+dz(ih)**2
                                                vsum_rh(ih) = vsum_rh(ih) + 1/sqrt(r2)
                                        end do
                                end do
                        end do
         end do
      end do
   end do
   vsum_rh(ih) = vsum_rh(ih)/(n_ic(ih)**2)
  end do                                                                  !   end loop over number of hierarchic structures

    vsum_rh(1:nh) = 1/vsum_rh(1:nh)                                       ! get hydrodynamic radius for each hierarchic structures
        var(3)%value = sum(vsum_rh(1:nh))/nh                                    ! mediate over hierarchic structures
!..............................................

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)

      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5,f15.0)', uout)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)

      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5,f15.0)', uout)

      deallocate(var,mimat,ro_temp,ro_loc)

   end select

end subroutine COMBAver

!************************************************************************
!> \page moluser moluser.F90
!! **SFPBC_COMB**
!! *calculate partial structure factors*
!************************************************************************

!     calculate actually the form factor as it works only for nh =1
!     form factor for backbone stored in ipt = 1
!     evaluated by summing over hierarchy and stored in ipt = 2
!     form factor for all side chains stored in ipt = 3
!     bugs: does not calcuate scaterring intensities
!     does not calculate scattering intensities
!     By Daniel A.

subroutine SFPBC_COMB(iStage)

   use MolModule
   implicit none

   integer(4), intent(in):: iStage

   character(40), parameter :: txroutine ='SFPBC_COMB'
   character(80), parameter :: txheading ='structure factor'
   integer(4)   , parameter :: ntype = 3
   character(14), parameter :: txtype(ntype) = ['100 directions', '110 directions', '111 directions']

   integer(4)   , save :: dirlow(ntype), dirupp(ntype)
   integer(4)   , save :: ndim, nbin, nvar

   type(df_var) , allocatable, save :: var(:)
   integer(4)   , allocatable, save :: ipnt(:,:)
   complex(8)   , allocatable, save :: eikr(:,:,:)
   real(8)      , allocatable, save :: q(:)
   real(8)      , allocatable, save :: sfpar(:,:), sfparsd(:,:)

   logical      , save :: lqsorted, lsi
   real(8)      , save :: sffac(3)
   integer(4)   , save :: ndiv = 1  ! factor by which the number of bins is reduced (in FLIST)
   integer(4)   , save :: n_ic(100)
   logical,       save :: ltype
   integer(4) :: ip, ipt, jpt, ivar, iptjpt, ibin, jbin(ntype), itype, nq, m, dupp, dlow
   integer(4) :: ih, ict, igen, iseg, ic, ilocal
   real(8)    :: norm, kx, ky, kz
   complex(8) :: eikr1(1:13), eikrtemp

   namelist /nmlSFPBC_COMB/ ndim, nbin, lqsorted, lsi, ltype

   if ((txbc /= 'xyz') .and. (txbc /= 'xy')) return

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2,txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ndim = 3
      nbin = 100
      lqsorted =.true.
      lsi =.false.
      ltype =.false.

      rewind(uin)
      read(uin,nmlSFPBC_COMB)

      if (ndim == 3) then
         dirlow(1:ndim) = [ 1, 4, 10 ]
         dirupp(1:ndim) = [ 3, 9, 13 ]
         if (boxlen(1) /= boxlen(2) .or. boxlen(1) /= boxlen(3)) call Stop(txroutine, 'boxlen(1) /= boxlen(2) .or. boxlen(1) /= boxlen(3)', uout)
      else if (ndim == 2) then
         dirlow(1:ndim) = [ 1, 3 ]
         dirupp(1:ndim) = [ 2, 4 ]
         if (boxlen(1) /= boxlen(2)) call Stop(txroutine, 'boxlen(1) /= boxlen(2)', uout)
      else
         call Stop(txroutine, 'ndim out of range', uout)
      end if
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

      sffac(1:3) = [ One, sqrt(Two), sqrt(Three) ]

   case (iWriteInput)

      nvar = nptpt*ndim
      allocate(var(nvar), ipnt(nptpt,ntype), eikr(0:nbin,1:13,1:npt), &
              q(3*nbin), sfpar(3*nbin,nptpt), sfparsd(3*nbin,nptpt))
      ipnt = 0
      eikr = cmplx(Zero,Zero)
      q = 0.0E+00
      sfpar = 0.0E+00
      sfparsd = 0.0E+00

      nvar = 0
      do itype = 1, ndim
         if(ltype) then
         do ipt = 1, npt
            do jpt = ipt, npt
               iptjpt = iptpt(ipt,jpt)
               nvar = nvar + 1
               ipnt(iptjpt,itype) = nvar
               var(nvar)%label = trim(txptpt(iptjpt))//' '//txtype(itype)
               var(nvar)%min = Half*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%max = (nbin+Half)*TwoPiBoxi(1)*sffac(itype)
               var(nvar)%nbin = nbin
            end do
         end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

      eikr(0:nbin,1:dirupp(ndim),1:npt) = cmplx(Zero,Zero)

! ... calculate sum( exp(i (k * r ))
!                i               i

    n_ic(1:nh) = 0

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        kx = TwoPiBoxi(1)*ro(1,ip)
                        ky = TwoPiBoxi(2)*ro(2,ip)
                        kz = TwoPiBoxi(3)*ro(3,ip)
                        if (ndim == 3) then
                                eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
                                eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
                                eikr1(3) = cmplx(cos(kz),sin(kz))              ! ( 0 0 1) direction
                                eikr1(4) = cmplx(cos(kx+ky),sin(kx+ky))        ! ( 1 1 0) direction
                                eikr1(5) = cmplx(cos(ky+kz),sin(ky+kz))        ! ( 0 1 1) direction
                                eikr1(6) = cmplx(cos(kz+kx),sin(kz+kx))        ! ( 1 0 1) direction
                                eikr1(7) = cmplx(cos(kx-ky),sin(kx-ky))        ! ( 1-1 0) direction
                                eikr1(8) = cmplx(cos(ky-kz),sin(ky-kz))        ! ( 0 1-1) direction
                                eikr1(9) = cmplx(cos(kz-kx),sin(kz-kx))        ! (-1 0 1) direction
                                eikr1(10) = cmplx(cos(kx+ky+kz),sin(kx+ky+kz)) ! ( 1 1 1) direction
                                eikr1(11) = cmplx(cos(kx+ky-kz),sin(kx+ky-kz)) ! ( 1 1-1) direction
                                eikr1(12) = cmplx(cos(kx-ky+kz),sin(kx-ky+kz)) ! ( 1-1 1) direction
                                eikr1(13) = cmplx(cos(kz+ky-kx),sin(kz+ky-kx)) ! (-1 1 1) direction
                        else if (ndim == 2) then
                                eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
                                eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
                                eikr1(3) = cmplx(cos(kx+ky),sin(kx+ky))        ! ( 1 1 0) direction
                                eikr1(4) = cmplx(cos(kx-ky),sin(kx-ky))        ! ( 1-1 0) direction
                        end if

                        do m = 1, dirupp(ndim)
                                eikrtemp = cmplx(One,Zero)
                                do ibin = 0, nbin-1
                                        eikrtemp = eikrtemp*eikr1(m)
                                        if(iptpn(ip)==1) then                                 ! backbone, stored in ipt = 1
                                                eikr(ibin,m,1) = eikr(ibin,m,1)+eikrtemp
                                        end if
                                        if(iptpn(ip)>1)  then                    ! store all side chains in ipt 3
                                                eikr(ibin,m,2) = eikr(ibin,m,2)+eikrtemp
                                        end if
                                end do
                        end do
                        n_ic(ih) =    n_ic(ih) + 1
         end do
      end do
    end do
  end do

           do m = 1, dirupp(ndim)
        do ibin = 0, nbin-1
                        eikr(ibin,m,2) = eikr(ibin,m,2) +eikr(ibin,m,1)
                end do
      end do

! ... sum over different equivalent k-vectors

      do itype = 1, ndim
           if(ltype) then
         dlow=dirlow(itype)
         dupp=dirupp(itype)
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(iptpt(ipt,jpt),itype)
               do ibin = 0, nbin-1
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ &
                  sum(real(eikr(ibin,dlow:dupp,ipt))*real(eikr(ibin,dlow:dupp,jpt))  &
                      +aimag(eikr(ibin,dlow:dupp,ipt))*aimag(eikr(ibin,dlow:dupp,jpt)))
               end do
            end do
         end do
                 end if
      end do

   case (iAfterMacrostep)

! ... normalizing

      do itype = 1, ndim
           if(ltype) then
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(iptpt(ipt,jpt),itype)
                           if (ipt ==1) then
               norm = One/((dirupp(itype)-dirlow(itype)+1)*sqrt(real(nppt(ipt))*real(nppt(jpt))))
                           else if (ipt ==2.AND.jpt ==2) then
                           norm = One/((dirupp(itype)-dirlow(itype)+1)*sqrt(real(n_ic(1)-nppt(1))*real(n_ic(1)-nppt(1))))
                           end if
               var(ivar)%avs2(-1:nbin) = var(ivar)%avs2(-1:nbin)*norm
            end do
         end do
                 end if
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      if (ndim == 3) write(uout,'(a,t35,i5,5x,a)')         'number of dimensions           = ', ndim, &
       'averaged over 3x(100), 6x(110), and 4x(111) directions'
      if (ndim == 2) write(uout,'(a,t35,i5,5x,a)')         'number of dimensions           = ', ndim, &
       'averaged over 2x(100) and 2x(110) directions'
      write(uout,'(a,t35,i5,5x,a,i5,2x,a)') 'number of grid points          = ', nbin
      write(uout,'(a,t35,g10.3)')           'q-sorted                       = ', lqsorted
      write(uout,'(a,t35,f8.2)')            'lower wavevector limit         = ', var(1)%min+0.5*TwoPiBoxi(1)
      write(uout,'(a,t35,f8.2)')            'upper wavevector limit         = ', var(1)%max-0.5*TwoPiBoxi(1)

      if (.not.lqsorted) then

         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      else

! ... store sf data in q-sorted order
         if(ltype) then
         do ipt = 1, npt
            do jpt = ipt, npt
               iptjpt = iptpt(ipt,jpt)
               jbin(1:ndim) = 1
               nq = 0
               if (ndim == 3) then
                  do ibin = 1, ndim*nbin
                     if (jbin(1)*sffac(1) < jbin(2)*sffac(2) .and. jbin(1)*sffac(1) < jbin(3)*sffac(3) .and. jbin(1) <= nbin) then
                        ilocal=1
                        call sortsub(ilocal)
                     else if (jbin(2)*sffac(2) < jbin(3)*sffac(3) .and. jbin(2) <= nbin) then
                        ilocal=2
                        call sortsub(ilocal)
                     else if (jbin(3) <= nbin) then
                        ilocal=3
                        call sortsub(ilocal)
                     end if
                  end do
               else if (ndim == 2) then
                  do ibin = 1, ndim*nbin
                     if (jbin(1)*sffac(1) < jbin(2)*sffac(2) .and. jbin(1) <= nbin) then
                        ilocal=1
                        call sortsub(ilocal)
                     else if (jbin(2) <= nbin) then
                        ilocal=2
                        call sortsub(ilocal)
                     end if
                  end do
               end if
            end do
         end do

! ... write sf data in q-sorted order

         if (ishow > 0) then
            write(uout,'()')
            write(uout,'(a)') 'q-sorted structure factor'
            do ipt = 1, npt
               do jpt = ipt, npt
                  iptjpt = iptpt(ipt,jpt)
                  write(uout,'(a)') txptpt(iptjpt)
                  write(uout,'(i5,3g15.5)') (ibin,q(ibin),sfpar(ibin,iptjpt),sfparsd(ibin,iptjpt),ibin = 1,nq)
               end do
            end do
         end if
         if (ilist > 0) then
            write(ulist,*) 'structure factor'
            write(ulist,*) nptpt
            do ipt = 1, npt
               do jpt = ipt, npt
                  iptjpt = iptpt(ipt,jpt)
                  write(ulist,*) 'sf '//txptpt(iptjpt)
                  write(ulist,*) nq/ndiv
                  write(ulist,'(g15.5,a,g15.5,a,g15.5)') (q(ibin),tab,sfpar(ibin,iptjpt),tab,sfparsd(ibin,iptjpt),ibin = 1,nq/ndiv)
               end do
            end do
         end if
                 end if

         if (lsi) call ScatIntens(nq, q, sfpar)

      end if

          deallocate(var, ipnt, eikr, q, sfpar, sfparsd)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine sortsub(itype)
   integer(4), intent(in) :: itype                           ! type of direction
   ivar = ipnt(iptjpt,itype)                                 ! get variable id
   nq = nq+1                                                 ! number of items in q-sorted list
   q(nq) = jbin(itype)*TwoPiBoxi(itype)*sffac(itype)         ! get wave vector
   sfpar(ibin,iptjpt) = var(ivar)%avs1(jbin(itype)-1)        ! store value of sf
   sfparsd(ibin,iptjpt) = var(ivar)%avsd(jbin(itype)-1)      ! store uncertainty of sf
   jbin(itype) = jbin(itype)+1                               ! update jbin
end subroutine sortsub

!..................................................

end subroutine SFPBC_COMB

!************************************************************************
!> \page moluser moluser.F90
!! **OCF**
!! *calculate bond-bond orientational correlation function and where particle i belongs to group igr*
!************************************************************************

!     By Daniel A.

subroutine OCF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine = 'OCF'
   character(60), parameter :: txheading = 'cos(tetha_i,i+n) between i and i+n bonds'
   integer(4), parameter :: ntype = 3
   character(5),  parameter :: txtype(ntype) = ['cos_t','cos_t','cos_t']

   logical,       save :: ltype(ntype)
   integer(4),    save :: nvar

   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:), ip_loc(:,:), ip_low(:,:), ip_high(:,:)
   real(8),       allocatable, save :: ro_loc(:,:), ro_temp(:,:)
   real(8),        allocatable, save :: n(:,:,:,:), corr(:,:,:), corrh(:,:,:), corr2(:,:,:), corrp(:,:,:,:), corrmed(:,:,:)

   integer(4),    save  ::  n_excl
   integer(4) :: itype
   integer(4) :: ip, s,  k, ih, igen, ict, ic,iseg, p_macro
   integer(4) :: unit = 99

   namelist /nmlOCF/ ltype, n_excl

   if(slave) return

   if(ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ltype(1:ntype) =.false.
          n_excl = 0

      rewind(uin)
      read(uin,nmlOCF)

   case (iWriteInput)

        allocate(var(nvar), ipnt(nptpt, ntype), ro_loc(3,np), ro_temp(3,na), n(3,np,nc,nc), corr(nct,nc,np), corrh(nct,nc,np), corr2(nc,nc,np), corrp(nc,nc,20,np), corrmed(nc,nc,np) )
    allocate (ip_loc(nc,nc), ip_low(nc,nc), ip_high(nc,nc))

        do ict = 1, nct                        ! loop over chain types,   ntype == nct!
                if (ltype(ict) .eqv. .true.) then
                        do ic = icnct(ict), icnct(ict)+ncct(ict)-1        ! loop over chains
                                do s = 1, npct(ict)-1-2*n_excl                         ! loop over segments
                                        corr(ict,ic,s) = 0
                                end do
                        end do
                end if
        end do

   case (iBeforeSimulation)

    if(lsim .and. master .and. txstart == 'continue')  read(ucnf)  var

   case (iBeforeMacrostep)

        do ict = 1, nct                        ! loop over chain types
                if (ltype(ict) .eqv. .true.) then
                        do ic = icnct(ict), icnct(ict)+ncct(ict)-1        ! loop over chains
                                do s = 1, npct(ict)-1-2*n_excl                         ! loop over segments
                                        corrh(ict,ic,s) = 0
                                end do
                        end do
                end if
        end do


   case (iSimulationStep)

         vaux(1:3,1:np) = 0.0
         ro_loc(1:3,1:np) = 0.0

!...  undo conditions for hierarchical structure

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
        do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
            if (ict ==1) then                                                 ! UNDO backbone (first chain type)
                                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                                ro_temp(1:3,1:na) = vaux(1:3,1:na)                         ! store undo backbone
                        end if
                          if (nbondcl(ip) ==1.AND.ict > 1 ) call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)      ! UNDO grafted chains, taking  UNDO positions
                                                                                                                                                          !linked beads belonging to the backbone as reference
         end do
      end do
        end do
   end do

 ! ... store position of the undo structure in ro_loc(1:3,ip)
   do ih = 1, nh                                                            ! loop over number of hierarchic structures
        do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
         end do
      end do
        end do
   end do

 ! ... undo conditions for linear chains

    do ic = 1, nc                            ! undo the linear chains
                ict = ictcn(ic)
                ip = ipnsegcn(1,ic)
                if(ihnpn(ip) /= 0) cycle                                            ! exclude hierarchical structures
                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, -1, vaux)
                do iseg = 1, npct(ict)
                        ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
                end do
    end do


 !.... calculate orientational correlations

        do itype = 1, ntype                    !   ntype == nct
          if(ltype(itype)) then
                do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                        ip_low(itype,ic)  = ipnsegcn(1 + n_excl,ic)
                        ip_high(itype,ic) = ipnsegcn(1 + npct(itype) -2-n_excl,ic)
                        ip_loc(itype,ic) = 0
                end do

        do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                 do ip = ip_low(itype,ic), ip_high(itype,ic)           ! calculate bond orientation
                        ip_loc(itype,ic) = ip_loc(itype,ic)+1
                        n(1,ip_loc(itype,ic),itype,ic) = (ro_loc(1,ip+1)-ro_loc(1,ip))/sqrt((ro_loc(1,ip+1)-ro_loc(1,ip))**2+(ro_loc(2,ip+1)-ro_loc(2,ip))**2+(ro_loc(3,ip+1)-ro_loc(3,ip))**2)
                        n(2,ip_loc(itype,ic),itype,ic) = (ro_loc(2,ip+1)-ro_loc(2,ip))/sqrt((ro_loc(1,ip+1)-ro_loc(1,ip))**2+(ro_loc(2,ip+1)-ro_loc(2,ip))**2+(ro_loc(3,ip+1)-ro_loc(3,ip))**2)
                        n(3,ip_loc(itype,ic),itype,ic) = (ro_loc(3,ip+1)-ro_loc(3,ip))/sqrt((ro_loc(1,ip+1)-ro_loc(1,ip))**2+(ro_loc(2,ip+1)-ro_loc(2,ip))**2+(ro_loc(3,ip+1)-ro_loc(3,ip))**2)
                        end do
        end do

                do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                        do s = 1, npct(itype)-2-2*n_excl + 1  !   s = 1, npct(ict)-1-2*n_excl
                                do k = 1, npct(itype)-s-2*n_excl
                                        corr(itype,ic,s)  = corr(itype,ic,s)  + n(1,k,itype,ic)*n(1,s+k-1,itype,ic)+n(2,k,itype,ic)*n(2,s+k-1,itype,ic)+n(3,k,itype,ic)*n(3,s+k-1,itype,ic)
                                        corrh(itype,ic,s) = corrh(itype,ic,s) + n(1,k,itype,ic)*n(1,s+k-1,itype,ic)+n(2,k,itype,ic)*n(2,s+k-1,itype,ic)+n(3,k,itype,ic)*n(3,s+k-1,itype,ic)
                                end do
                        end do
        end do
           end if
        end do

   case (iAfterMacrostep)

        do itype = 1, ntype
          if(ltype(itype)) then
                        do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                                do s = 1, npct(itype)-2-2*n_excl + 1
                                        corr(itype,ic,s)  = corr(itype,ic,s)/(npct(itype)-s-2*n_excl)
                                        corrh(itype,ic,s) = corrh(itype,ic,s)/(npct(itype)-s-2*n_excl)
                                end do
                                do s = 1, npct(itype)-2-2*n_excl +1
                                        corrp(itype,ic,istep1,s) = corrh(itype,ic,s)/(nstep2/istatic)
                                end do
                        end do
                open(unit, status = 'unknown')
                write(unit,*) 'macrostep ' , istep1
                        do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                                do s = 1, npct(itype)-2-2*n_excl +1
                                        write(unit,*) itype,ic, s, corrp(itype,ic,istep1,s)
                                end do
                        end do
          end if
        end do

   case (iAfterSimulation)

        do itype = 1, ntype
                if(ltype(itype)) then
                        do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                                do s = 1, npct(itype)-2-2*n_excl +1     !!!!! replace  ncct(itype)  with  nptct(itype)  !!!!!
                                        corr2(itype,ic,s) = 0
                                        corrmed(itype,ic,s) = 0
                                        do p_macro= 1, nstep1
                                                corrmed(itype,ic,s) = corrmed(itype,ic,s) + corrp(itype,ic,p_macro,s)
                                        end do
                                        corrmed(itype,ic,s) = corrmed(itype,ic,s)/nstep1
                                        do p_macro= 1, nstep1
                                                corr2(itype,ic,s) = corr2(itype,ic,s)+ (corrp(itype,ic,p_macro,s)- corrmed(itype,ic,s))**2
                                        end do
                                end do
                        end do
                        write(ulist,*)   'final orientational correlation function'
                        do ic = icnct(itype), icnct(itype)+ncct(itype)-1
                                do s = 1, npct(itype)-2-2*n_excl +1
                                        write(ulist,'(3i3,2f12.4)') itype,ic, s, corrmed(itype,ic,s), sqrt(corr2(itype,ic,s)/nstep1/(nstep1-1))
                                end do
                        end do
                        close(unit)
                end if
        end do

        deallocate(var,ipnt, ro_loc, ro_temp, n, corr, corrh, corr2,  corrp,  corrmed, ip_loc, ip_low, ip_high)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OCF

!************************************************************************
!> \page moluser moluser.F90
!! **OCF_DF**
!! *calculate  type distribution functions of comb polymer type (coarse model)*
!************************************************************************

!     calculate distribution functions for bond-bond orientational correlations cos(theta_i,i+n)
!     for nine given values of n: 5, 10, 15, 20, 25, 30, 40, 50, 60
!     By Daniel A.

subroutine OCF_DF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine = 'OCF_DF'
   character(60), parameter :: txheading = 'distribution functions of the bond-bond orientational correlations '
   integer(4), parameter :: ntype = 9
   character(5),  parameter :: txtype(ntype) = ['oc(1)','oc(2)','oc(3)','oc(4)','oc(5)','oc(6)','oc(7)','oc(8)','oc(9)']

   logical,       save :: ltype(ntype)
   real(8),       save :: vmin(ntype), vmax(ntype)
   integer(4),    save  ::    pipt, n_excl
   integer(4),    save :: nbin, nvar

   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),       allocatable, save :: ro_temp(:,:), ro_loc(:,:), n(:,:), corr(:), corrh(:)

   integer(4)     :: ic, ip, ict, itype, ivar, ibin
   real(8)        :: properties(13), value
   integer(4)     :: iseg, ih, igen, s,  k, ip_loc, ip_low, ip_high

   namelist /nmlOCF_DF/ ltype, vmin, vmax, nbin, pipt, n_excl

   if(ltrace) call WriteTrace(2,txroutine, iStage)

   if(slave) return

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ltype(1:ntype) =.false.
      vmin(1:9) =   -1.0
      vmax(1:9) =   1.0

      nbin  =  100

      pipt = 1
          n_excl = 0
      rewind(uin)
      read(uin,nmlOCF_DF)

      if(nbin > mnbin_df) call Stop('OCF_DF', 'nbin > mnbin_df', uout)

   case (iWriteInput)

          nvar =  nct*count(ltype(1:ntype) )
      allocate(var(nvar), ipnt(nptpt, ntype), ro_loc(3,np), ro_temp(3,na), n(3,np), corr(np), corrh(np))
      ipnt = 0
      ro_loc = 0.0E+00
      ro_temp = 0.0E+00
      n = 0.0E+00
      corr = 0.0E+00
      corrh = 0.0E+00

          do s = 1, np
                corr(s) = 0
      end do

      nvar = 0
    do itype = 1, ntype
         if(ltype(itype)) then
            do ict = 1, nct
               nvar = nvar+1
               ipnt(ict,itype) = nvar
               var(nvar)%label = trim(txtype(itype))//' '//txct(ict)
               var(nvar)%min = vmin(itype)
               var(nvar)%max = vmax(itype)
               var(nvar)%nbin = nbin
            end do
         end if
    end do

    call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

     call DistFuncSample(iStage, nvar, var)
     if(lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

        do s = 1, nppt(pipt)-1-2*n_excl
                corrh(s) = 0
          end do
    call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

         vaux(1:3,1:np) = 0.0
         ro_loc(1:3,np) = 0.0

!...  undo conditions for hierarchical structure

   do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
             if (ict ==1) then                                                 ! UNDO backbone (first chain type)
                                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                                ro_temp(1:3,1:na) = vaux(1:3,1:na)                         ! store undo backbone
                         end if
                          if (nbondcl(ip) ==1.AND.ict > 1 )  call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)      ! UNDO grafted chains, taking  UNDO positions
                                                                                                                                                          !linked beads belonging to the backbone as reference
         end do
      end do
    end do
   end do

 ! ... store position of the undo structure in ro_loc(1:3,ip)
  do ih = 1, nh                                                            ! loop over number of hierarchic structures
    do igen = 0, ngen                                                       ! loop over generations
      ict = ictgen(igen)                                                   ! chain type
      do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
         do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                        ro_loc(1:3,ip) = vaux(1:3,ip)
         end do
      end do
   end do
 end do

 !.... calculate orientational correleations

        do itype = 1, ntype
          if(ltype(itype)) then
                 ip_low  = ipnpt(pipt) + n_excl
                 ip_high = ipnpt(pipt)+nppt(pipt)-2-n_excl
                 ip_loc = 0
                do ip = ip_low, ip_high              ! calculate bond orientation
                    ip_loc = ip_loc+1
                        n(1,ip_loc) = (ro(1,ip+1)-ro(1,ip))/sqrt((ro(1,ip+1)-ro(1,ip))**2+(ro(2,ip+1)-ro(2,ip))**2+(ro(3,ip+1)-ro(3,ip))**2)
                        n(2,ip_loc) = (ro(2,ip+1)-ro(2,ip))/sqrt((ro(1,ip+1)-ro(1,ip))**2+(ro(2,ip+1)-ro(2,ip))**2+(ro(3,ip+1)-ro(3,ip))**2)
                        n(3,ip_loc) = (ro(3,ip+1)-ro(3,ip))/sqrt((ro(1,ip+1)-ro(1,ip))**2+(ro(2,ip+1)-ro(2,ip))**2+(ro(3,ip+1)-ro(3,ip))**2)
                end do
        do s = 1, np
                 corr(s) = 0
                end do

                do s = 1, nppt(pipt)-2-2*n_excl + 1
                        do k = 1, nppt(pipt)-s-2*n_excl
                        corr(s)  = corr(s)  + n(1,k)*n(1,s+k-1)+n(2,k)*n(2,s+k-1)+n(3,k)*n(3,s+k-1)
                        corrh(s) = corrh(s) + n(1,k)*n(1,s+k-1)+n(2,k)*n(2,s+k-1)+n(3,k)*n(3,s+k-1)
                        end do
                end do
      end if
        end do

        do s = 1, nppt(pipt)-2-2*n_excl + 1
                corr(s) = corr(s)/(nppt(pipt)-s-2*n_excl)
        end do

        properties(1) = corr(5)
    properties(2) = corr(10)
        properties(3) = corr(15)
        properties(4) = corr(20)
        properties(5) = corr(25)
        properties(6) = corr(30)
        properties(7) = corr(40)
        properties(8) = corr(50)
        properties(9) = corr(60)

! ... summation of type 1 to type 9

         do itype = 1, 9
            if(ltype(itype)) then
               ivar = ipnt(1,itype)
               if(itype == 1) then
                  value = properties(1)
               else if(itype == 2) then
                  value = properties(2)
               else if(itype == 3) then
                  value = properties(3)
               else if(itype == 4) then
                  value = properties(4)
               else if(itype == 5) then
                  value = properties(5)
               else if(itype == 6) then
                  value = properties(6)
               else if(itype == 7) then
                  value = properties(7)
               else if(itype == 8) then
                  value = properties(8)
                        else if(itype == 9) then
                  value = properties(9)
               end if
               ibin = max(-1,min(int(var(ivar)%bini*(value-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
         end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if(lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

         deallocate(var, ipnt, ro_loc, ro_temp, n, corr, corrh)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OCF_DF


!***********************************************************************
!> \page moluser moluser.F90
!! **ChainBeadBeadContact**
!! *documentation_missing*
!************************************************************************

!     By Daniel A.

subroutine ChainBeadBeadContact(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40)   , parameter :: txroutine = 'ChainBeadBeadContact'
   character(60)   , parameter :: txheading = 'cbead-cbead complexation'

   integer(4)      , save :: nvar

   type(scalar_var), allocatable, save :: var(:)
   real(8),       allocatable, save :: ro_temp(:,:), ro_loc(:,:), r_dist(:)
   integer(4),    save :: iptpart_bead
   real(8),       save :: rcontact_bead
   logical,       save :: ltype

   integer(4)     :: ih, ic, ict, ip, jp, jploc, iseg, igen
   real(8)        :: dx, dy, dz, r2

   namelist /nmlCBCB/  ltype, iptpart_bead, rcontact_bead

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2,txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)
      ltype = .false.
      if(.not.lclink) call stop(txroutine, '.not.lclink',uout)

          rewind(uin)
      read(uin,nmlCBCB)

        case (iWriteInput)

          nvar = npct(1)

      if(.not.lclink) call stop(txroutine, '.not.lclink',uout)
          allocate(var(nvar), ro_temp(3,np), ro_loc(3,np), r_dist(np))
          ro_temp = 0.0E+00
          ro_loc = 0.0E+00
          r_dist = 0.0E+00

      var(1:nvar)%label = 'complexation                = '
      var(1:nvar)%norm = One

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

        r_dist(1:nvar) = Zero
        if(ltype) then
        ! ... undo hierarchical structure

        do ih = 1, nh                                                            ! loop over number of hierarchic structures
                do igen = 0, ngen                                                       ! loop over generations
                ict = ictgen(igen)                                                   ! chain type
                        do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1          ! loop over chains of the structure
                                do iseg = 1, npct(ict)                                            ! loop over segments
                                        ip = ipnsegcn(iseg,ic)
                                        if (ict ==1) then                                                 ! UNDO backbone (first chain type)
                                                call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
                                                ro_temp(1:3,1:ip) = vaux(1:3,1:ip)                         ! store undo backbone
                                        end if
                                        if (nbondcl(ip) ==1.AND.ict > 1 ) call UndoPBCChain(ro_temp(1,bondcl(1,ip)), ic, 1, vaux)      ! UNDO grafted chains, taking  UNDO positions
                                                                                                                                                          !linked beads belonging to the backbone as reference
                                end do
                        end do
                end do
        end do

!  undo non-hierarchical (linear) chains taking as reference the backbone ic = 1
    do ic = 1, nc
                ict = ictcn(ic)
                ip = ipnsegcn(1,ic)
                if(ihnpn(ip) /= 0) cycle                                            ! exclude hierarchical structures
                call UndoPBCChain(ro_temp(1,ipnsegcn(1,1)), ic, -1, vaux)        ! it assumes ic =1 is the backbone of the hierarchy
                do iseg = 1, npct(ict)                                            ! loop over segments
            ip = ipnsegcn(iseg,ic)
                    ro_loc(1:3,ip) = vaux(1:3,ip)
                end do
    end do


! ... calculate complexation between each beads of ic =1 (backbone) and beads belonging to ict=3 (linear polyion)

    do ic = 1, nc
        ict = ictcn(ic)
                if(ict/=1) cycle
        do iseg = 1, npct(ict)
            ip = ipnsegcn(iseg,ic)
            do jploc = 1, nppt(iptpart_bead)
               jp = ipnpt(iptpart_bead)+(jploc-1)
               if (ip == jp) call Stop('ChainBeadBeadContact', 'ip == jp', uout)
               dx = ro_temp(1,ip)-ro_loc(1,jp)
               dy = ro_temp(2,ip)-ro_loc(2,jp)
               dz = ro_temp(3,ip)-ro_loc(3,jp)
               r2 = dx**2+dy**2+dz**2
               if (r2 < rcontact_bead**2) then
                              r_dist(ip) = r_dist(ip) +1
               end if
            end do
        end do
    end do

          end if    !    end if(ltype)

        var(1:nvar)%value = r_dist(1:nvar)

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)

      if(ltype) then
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 0, '(a,t35,2f15.5,f15.0)', uout)

      end if
      deallocate(var, ro_temp, ro_loc, r_dist)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainBeadBeadContact

!************************************************************************
!> \page moluser moluser.F90
!! **ElPot**
!! *calculate electrostatic potential as a function of distance from the center of lab frame*
!************************************************************************

!     by Stefanie Schneider

subroutine ElPot(iStage)

   use MolModule
   implicit none

   integer, intent(in)   :: iStage

   character(40), parameter :: txroutine ='ElPot'
   character(80), parameter :: txheading ='potential radial distr function'
   integer(4)   , parameter :: ntype = 1
   real(8), save         :: vmin,vmax
   integer(4), save      :: nbin, nvar
   integer(4), save      :: nsamp1(-1:1000)
   type(df_var), allocatable, save :: var(:)
   integer(4)            :: ibin, ip, jp, ipt
   real(8)               :: ui, dx,dy, dz, r2, dr, uuu,pot, charge, norm, norm1, fsum(3)

   namelist /nmlElPot/ vmin, vmax, nbin

   if (slave) return

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

! ... set default values

      vmin = Zero
      vmax = 10.0
      nbin = 100

      rewind(uin)
      read(uin,nmlElPot)

      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

   case (iWriteInput)

      nvar = ntype
      allocate(var(nvar))

      var(1)%label = 'radial potential distribution'
      var(1)%min   = vmin
      var(1)%max   = vmax
      var(1)%nbin  = nbin

      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)
      nsamp1(-1:1000)=0
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      var(1)%nsampbin(-1:var(1)%nbin)=zero

    case (iSimulationStep)

      var(1)%nsamp2 = var(1)%nsamp2 + 1

      do ip = 1, np
         ui=0.0d0
         do jp = 1, np
            if (ip == jp) cycle
            call UTwoBodyPair(ip, jp, uuu, fsum)
            ui = ui+uuu
         end do
         ipt=iptpn(ip)
         charge=zat(ipt)
         if (abs(charge) > 1d-10) then
            pot=ui/charge/ech
         else
            pot=0.0d0
         end if

         dx = ro(1,ip)
         dy = ro(2,ip)
         dz = ro(3,ip)
         r2 = dx**2+dy**2+dz**2
         dr=sqrt(r2)

         if (luext) then
            if (txuext(ipt) == 'insulating_sphere') then
              if (dr >= rInsSphere) then
                 pot=pot+zInsSphere*EpsiFourPi/dr/ech
               else
                 pot=pot+zInssphere*EpsiFourPi/ech/(2.0d0*rInsSphere)*(3.0d0-r2/rInsSphere**2)
               end if
            end if
         end if
         pot=pot/avNo*sclene

         ibin = max(-1,min(int(var(1)%bini*(dr-var(1)%min)),var(1)%nbin))
         var(1)%nsampbin(ibin)=var(1)%nsampbin(ibin) + 1
         var(1)%avs2(ibin) = var(1)%avs2(ibin)+pot
      end do

   case (iAfterMacrostep)

      do ibin=-1,var(1)%nbin
         norm =zero
         if (var(1)%nsampbin(ibin) > 0) then
            nsamp1(ibin)=nsamp1(ibin) + 1
            norm=one/var(1)%nsampbin(ibin)
         else
         end if
         var(1)%avs2(ibin)=var(1)%avs2(ibin)*norm
         var(1)%avs1(ibin) = var(1)%avs1(ibin)+var(1)%avs2(ibin)
         var(1)%avsd(ibin) = var(1)%avsd(ibin)+var(1)%avs2(ibin)**2
      end do
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)
      norm = Zero
      norm1 = Zero

      do ibin = -1, var(1)%nbin
         if (nsamp1(ibin) > 0) norm = One/nsamp1(ibin)
         if (nsamp1(ibin) > 1) norm1 = One/(nsamp1(ibin)-1)

         var(1)%avs1(ibin) = var(1)%avs1(ibin)*norm
         var(1)%avsd(ibin) = sqrt(max(Zero,var(1)%avsd(ibin)*norm-var(1)%avs1(ibin)**2)*norm1)
      end do
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ElPot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **ImageUser**
!! *driver of user-provided routines for generating image input files*
!************************************************************************


subroutine ImageUser(iStage)

   use DumpModule
   use MolModule, only : ltrace
   implicit none
   character(40), parameter :: txroutine ='ImageUser'

   integer(4), intent(in) :: iStage

   if (ltrace) call WriteTrace(2,txroutine, iStage)
!  call ximage(iStage)
   call Stop(txroutine, 'no user-provided image files is available', uout)

end subroutine ImageUser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
!> \page moluser moluser.F90
!! **JosUser**
!! *auxillary routine for syncronizing input data for Jos' project*
!************************************************************************


subroutine JosUser(iMode)

   use MolModule
   implicit none

   integer(4), intent(in) :: iMode

   character(40), parameter :: txroutine ='JosUser'
   real(8), parameter :: conc_salt = 0.001D0                 ! salt concentration in M
   real(8), parameter :: rho_salt = conc_salt*AvNo*1d-27     ! salt concentration in Å**-3

   write(*,'(2a,i4)') trim(txroutine), '  iMode = ', iMode

   if (iMode == 1) then ! adjust nppt(1) and nppt(2) according to input variables (call from particle.F90)

      nppt(1) = nppt(1) - nppt(2)
      nppt(2) = rho_salt*cyllen*pi*(cylrad**2 - rcylinder**2)
      nppt(1) = nppt(1) + nppt(2)
      write(*,'(a,f10.5)') 'cylrad (in)    =', cylrad
      write(*,'(a,i10  )') 'nppt(2) (new)  =', nppt(2)

   else if (iMode == 2) then ! adjust cylrad to nppt(2) (call from molsim.F90)

      cylrad = sqrt(nppt(2)/(rho_salt*pi*cyllen)+rcylinder**2)
      write(*,'(a,i10  )') 'nppt(2) (in)   =', nppt(2)
      write(*,'(a,f10.5)') 'cylrad (new)   =', cylrad
      write(*,*)
      call SetBoxParam

   end if

end subroutine JosUser

!************************************************************************
!> \page moluser moluser.F90
!! **ComplexationModule**
!! *Module for analysing the Complexation*
!************************************************************************


!> \page nmlComplexation
!! The namelist \ref nmlComplexation contains variables that can be used for the analysis of interparticle complexation.
!! * Variables:
!!  * \subpage rcut_complexation
!!  * \subpage lClusterDF
!!  * \subpage lComplexFraction
!!  * \subpage lComplexDist
!!  * \subpage lSegmentComplex

!> \page lClusterDF
!! `logical`
!! **default:** `.false.`
!! * `.true.`: The cluster size distribution  will be calculated.
!! * `.false.`: Nothing.

!> \page lComplexFraction
!! `logical`
!! **default:** `.false.`
!! * `.true.`: The fraction of complexation will be calculated.
!! * `.false.`: Nothing.

!> \page lComplexDist
!! `logical`
!! **default:** `.false.`
!! * _documentation missing_

!> \page lSegmentComplex
!! `logical`
!! **default:** `.false.`
!! * _documentation missing_

module ComplexationModule
   implicit none
   private
   public  ComplexationDriver

!> \page rcut_complexation
!! `real`
!! **default:** `0.0`
!! * Cut-off distance for complexation.
   real(8)  :: rcut_complexation
   real(8)  :: r2cut_cmplx

   logical, allocatable :: lcmplx_ipjp(:,:)

   contains

      !************************************************************************
      !*                                                                      *
      !*     ComplexationDriver                                               *
      !*                                                                      *
      !************************************************************************

      ! ... Driver for the Complexation Analysis

      subroutine ComplexationDriver(iStage)
         use MolModule, only: ltrace, ltime, uout, master, uin
         use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
         use MolModule, only: np
         implicit none
         integer(4), intent(in)  :: iStage
         character(40), parameter :: txroutine ='ComplexationDriver'
         character(80), parameter :: txheading ='complexation Analysis'
         logical,       save :: lClusterDF, lComplexFraction, lSegmentComplex, lComplexDist

         namelist /nmlComplexation/ rcut_complexation, lClusterDF, lComplexFraction, lComplexDist, lSegmentComplex

         if (ltrace) call WriteTrace(2, txroutine, iStage)
         if (ltime) call CpuAdd('start', txroutine, 0, uout)

         select case (iStage)
         case (iReadInput)
            rcut_complexation = 0.0
            lClusterDF = .false.
            lComplexFraction = .false.
            lComplexDist    = .false.
            lSegmentComplex = .false.
            rewind(uin)
            read(uin,nmlComplexation)

            call ComplexationDriverSub

         case (iWriteInput)
            if(.not. allocated(lcmplx_ipjp)) then
               allocate(lcmplx_ipjp(np,np))
            end if
            r2cut_cmplx = rcut_complexation**2
            call ComplexationDriverSub

         case (iBeforeSimulation)
            call ComplexationDriverSub

         case (iBeforeMacrostep)
            call ComplexationDriverSub

         case (iSimulationStep)
            call GetComplex(iStage)
            call ComplexationDriverSub

         case (iAfterMacrostep)
            call ComplexationDriverSub

         case (iAfterSimulation)
            if (master) then
               call WriteHead(2, txheading, uout)
               write(uout,'(a,t35,e13.6)')     'cutoff-distance                = ', rcut_complexation
               write(uout,'(a)') 'static analysis routines used'
               write(uout,'(a)') '-----------------------------'
               if (lClusterDF)        write(uout,'(a)') '   ClusterDF      '
               if (lComplexFraction)  write(uout,'(a)') '   ComplexFraction'
               if (lComplexDist)      write(uout,'(a)') '   ComplexDist    '
               if (lSegmentComplex)   write(uout,'(a)') '   SegmentComplex '
            end if

            call ComplexationDriverSub

         end select

         if (ltime) call CpuAdd('stop', txroutine, 1, uout)

         contains
            subroutine ComplexationDriverSub
               if (lComplexFraction)  call ComplexFraction(iStage)
               if (lComplexDist)      call ComplexDistribution(iStage)
               if (lSegmentComplex)   call SegmentComplex(iStage)
               if (lClusterDF)        call ClusterDF(iStage)
               continue
            end subroutine ComplexationDriverSub
      end subroutine

      !************************************************************************
      !*                                                                      *
      !*     GetComplex                                                       *
      !*                                                                      *
      !************************************************************************

      ! ... Detect which particles form a Complex

      subroutine GetComplex(iStage)
         use MolModule, only: ltrace, ltime, uout, ro, np
         implicit none
         integer(4), intent(in)  :: iStage
         character(40), parameter :: txroutine ='GetComplex'

         integer(4)  :: ip, jp
         real(8)  :: d(3), r2

         if (ltrace) call WriteTrace(3, txroutine, iStage)
         if (ltime) call CpuAdd('start', txroutine, 0, uout)

         if(.not. allocated(lcmplx_ipjp)) then
            call Stop(txroutine, 'lcmplx_ipjp is not allocated!', uout)
         end if

         lcmplx_ipjp = .false.

         do ip = 1, np-1
            do jp = ip + 1, np
               d(1:3) = ro(1:3,ip) - ro(1:3,jp)
               call PBCr2(d(1), d(2), d(3), r2)
               if (r2 .le. r2cut_cmplx) then
                  lcmplx_ipjp(ip,jp) = .true.
                  lcmplx_ipjp(jp,ip) = .true.
               end if
            end do
         end do

      end subroutine GetComplex

      !************************************************************************
      !*                                                                      *
      !*     ComplexFraction                                                  *
      !*                                                                      *
      !************************************************************************

      ! ... calculate how many particles are complexed

      subroutine ComplexFraction(iStage)

         use MolModule, only: ltrace, ltime, uout, master, lsim, txstart, ucnf
         use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
         use MolModule, only: np, npt, txpt, nppt, ipnpt, iptpn
         use StatisticsModule, only: scalar_var
         implicit none

         integer(4), intent(in) :: iStage
! These are documented in the manual in Chapter 7 (file datastructures.md)
         type cluster_var
            integer(4), allocatable :: ip
            integer(4), allocatable :: np
         end type cluster_var

         character(40), parameter :: txroutine ='ComplexFraction'
         character(80), parameter :: txheading ='Fraction of Complexation'
         integer(4), save         :: nvar
         type(scalar_var), allocatable, save :: var(:)

         integer(4)  :: ipt, jpt, ivar, ip

         if (ltrace) call WriteTrace(3, txroutine, iStage)
         if (ltime) call CpuAdd('start', txroutine, 1, uout)

         select case (iStage)
         case (iReadInput)
            nvar = npt**2
            if(.not. allocated(var)) then
               allocate(var(nvar))
            end if

         case (iWriteInput)

            do ipt = 1, npt
               do jpt = 1, npt
                  ivar = ivar_ptpt(ipt,jpt)
                  var(ivar)%label = 'w(cmplx): '//trim(txpt(ipt))//' - '//trim(txpt(jpt))
                  var(ivar)%norm = 1.0d0/nppt(ipt)
               end do
            end do

         case (iBeforeSimulation)

            call ScalarSample(iStage, 1, nvar, var)
            if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

         case (iBeforeMacrostep)

            call ScalarSample(iStage, 1, nvar, var)

         case (iSimulationStep)

            var%value = 0.0d0
            do ip = 1, np
               ipt = iptpn(ip)
               do jpt = 1, npt
                  if( any(lcmplx_ipjp( ipnpt(jpt):(ipnpt(jpt) + nppt(jpt) - 1), ip ))) then !if any particle of type jpt is complexed withparticle ip
                     ivar = ivar_ptpt(ipt,jpt)
                     var(ivar)%value = var(ivar)%value + 1.0d0
                  end if
               end do
            end do

            call ScalarSample(iStage, 1, nvar, var)


         case (iAfterMacrostep)

            call ScalarSample(iStage, 1, nvar, var)
            if (lsim .and. master) write(ucnf) var
            call ScalarNorm(iStage, 1, nvar, var, 1)

         case (iAfterSimulation)

            call ScalarSample(iStage, 1, nvar, var)
            call ScalarNorm(iStage, 1, nvar, var, 1)
            if(master) then
               call WriteHead(2, txheading, uout)
               call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)
            endif
            deallocate(var)

         end select

         if (ltime) call CpuAdd('stop', txroutine, 1, uout)

      contains

         pure function ivar_ptpt(ipt, jpt) result(ivar)
            use MolModule, only: npt
            implicit none
            integer(4), intent(in)  :: ipt
            integer(4), intent(in)  :: jpt
            integer(4)  :: ivar
            ivar = (ipt-1)*npt + jpt
         end function ivar_ptpt

      end subroutine ComplexFraction

      !************************************************************************
      !*                                                                      *
      !*     ComplexDistribution                                              *
      !*                                                                      *
      !************************************************************************

      ! ... calculate the distribution functions of the complexation

!> \page nmlComplexDist
!! The namelist  \ref nmlComplexDist contains variables that control the calculation of complexation distribution functions.
!! Any combination of the types of distribution functions listed below may be selected
!! through vtype\%l.
!!
!!    | type | label | quantity                 |
!!    | ---- | ----- | -----------------------  |
!!    | 1    | w     | fraction of complexation |
!!
!! * Variables:
!!  * \subpage nmlComplexDist_vtype

!> \page nmlComplexDist_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`
!! * Flag for engagement, lower end, upper end, number of bins. Other flags are not used.
!! * Min: /-0.005,-0.005,-0.005/ Max: /1,005,1.005,1.005/

      subroutine ComplexDistribution(iStage)

         use MolModule, only: ltrace, ltime, uin, uout, master, lsim, txstart, ucnf, ulist, ishow, iplot, ilist
         use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
         use MolModule, only: npt, ipnpt, nppt, txpt
         use StatisticsModule, only: df_var, mnbin_df
         use MolModule, only: static1D_var
         use MollibModule, only: InvInt
         implicit none

         integer(4), intent(in) :: iStage

         character(40), parameter :: txroutine ='ComplexDist'
         character(80), parameter :: txheading ='Distribution of the rate of complexation'
         integer(4), parameter    :: ntype = 1
         type(static1D_var), save :: vtype(ntype)
         integer(4), save         :: nvar
         type(df_var), allocatable, save :: var(:)
         integer(4), allocatable, save :: ivariptjptitype(:,:,:)

         integer(4)  :: ipt, jpt
         integer(4)  :: itype, ivar, ibin
         real(8)  :: w !fraction of complexation

         namelist /nmlComplexDist/ vtype

         if (ltrace) call WriteTrace(3, txroutine, iStage)
         if (ltime) call CpuAdd('start', txroutine, 1, uout)

         select case (iStage)
         case (iReadInput)

            vtype(1)%l      = .false.
            vtype(1)%min    = -0.005d0
            vtype(1)%max    = 1.005d0+epsilon(vtype(1)%max)
            vtype(1)%nbin   = 101

            rewind(uin)
            read(uin,nmlComplexDist)

            if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

         case (iWriteInput)

            vtype%label = ['w']
            vtype(1)%nvar = npt**2
            nvar = sum(vtype(:)%nvar, 1, vtype(:)%l)
            if(.not. allocated(var)) then
               allocate(var(nvar))
            end if
            if(.not. allocated(ivariptjptitype)) then
               allocate(ivariptjptitype(npt,npt,ntype))
            end if

            ivar = 0
            do itype = 1, ntype
               if (vtype(itype)%l) then
                  do ipt = 1, npt
                     do jpt = 1, npt
                        ivar = ivar + 1
                        ivariptjptitype(ipt,jpt,itype) = ivar
                        var(ivar)%label = trim(vtype(itype)%label)//': pt: '//trim(txpt(ipt))//'; pt:'//trim(txpt(jpt))
                        var(ivar)%min   = vtype(itype)%min
                        var(ivar)%max   = vtype(itype)%max
                        var(ivar)%nbin  = vtype(itype)%nbin
                        var(ivar)%norm  = InvInt(nppt(ipt))
                     end do
                  end do
               end if
            end do

            call DistFuncSample(iStage, nvar, var)

         case (iBeforeSimulation)

            call DistFuncSample(iStage, nvar, var)
            if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

         case (iBeforeMacrostep)

            call DistFuncSample(iStage, nvar, var)

         case (iSimulationStep)

            var%nsamp2 = var%nsamp2 + 1

            do ipt = 1, npt
               do jpt = 1, npt

                  !sample of type 1
                  itype = 1
                  if (vtype(itype)%l) then
                     ivar = ivariptjptitype(ipt,jpt,itype)
                     w = count(any(lcmplx_ipjp( ipnpt(jpt):(ipnpt(jpt) + nppt(jpt) - 1), ipnpt(ipt):(ipnpt(ipt)+nppt(ipt)-1)),DIM=1 ))
                     w=w*var(ivar)%norm
                     ibin = max(-1,min(floor(var(ivar)%bini*(w-var(ivar)%min)),var(ivar)%nbin))
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+1.0d0
                  end if

               end do
            end do

         case (iAfterMacrostep)

            call DistFuncNorm(1, nvar, var)
            call DistFuncSample(iStage, nvar, var)
            if (lsim .and. master) write(ucnf) var

         case (iAfterSimulation)

            call DistFuncSample(iStage, nvar, var)
            call WriteHead(2, txheading, uout)
            call DistFuncHead(nvar, var, uout)
            call DistFuncAverValue(nvar, var, uout)
            call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

         deallocate(var)

      end select

      if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   end subroutine ComplexDistribution

      !************************************************************************
      !*                                                                      *
      !*     SegmentComplex                                                   *
      !*                                                                      *
      !************************************************************************

      ! ... calculate how many particles are complexed

      subroutine SegmentComplex(iStage)

         use MolModule, only: ltrace, ltime, uout, master, lsim, txstart, ucnf, ulist, ishow, iplot, ilist
         use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
         use MolModule, only: nct, npct, ncct, nc, ictcn, ipnsegcn, txct
         use MolModule, only: npt, ipnpt, nppt, txpt
         use StatisticsModule, only: df_var
         use MollibModule, only: InvInt
         implicit none

         integer(4), intent(in) :: iStage

         character(40), parameter :: txroutine ='SegmentComplex'
         character(80), parameter :: txheading ='Fraction of Complexation of each segment of the chains'
         integer(4), save         :: nvar
         type(df_var), allocatable, save :: var(:)
         integer(4), allocatable, save :: ivarictipt(:,:)

         integer(4)  :: ict, ic, iseg
         integer(4)  :: ip, ipt
         integer(4)  :: ivar

         if (ltrace) call WriteTrace(3, txroutine, iStage)
         if (ltime) call CpuAdd('start', txroutine, 1, uout)

         select case (iStage)
         case (iReadInput)

            nvar = nct*npt
            if(.not. allocated(var)) then
               allocate(var(nvar))
            end if
            if(.not. allocated(ivarictipt)) then
               allocate(ivarictipt(nct,npt))
            end if

         case (iWriteInput)

            ivar = 0
            do ict = 1, nct
               do ipt = 1, npt
                  ivar = ivar + 1
                  ivarictipt(ict,ipt) = ivar
                  var(ivar)%label = 'ct: '//txct(ict)//'; pt:'//txpt(ipt)
                  var(ivar)%min   = 0.5d0
                  var(ivar)%max   = npct(ict) + 0.5d0
                  var(ivar)%nbin  = npct(ict)
                  var(ivar)%norm  = InvInt(ncct(ict))
               end do
            end do
            call DistFuncSample(iStage, nvar, var)

         case (iBeforeSimulation)

            call DistFuncSample(iStage, nvar, var)
            if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

         case (iBeforeMacrostep)

            call DistFuncSample(iStage, nvar, var)

         case (iSimulationStep)

            var%nsamp2 = var%nsamp2 + 1

            do ic = 1, nc
               ict = ictcn(ic)
               do iseg = 1, npct(ict)
                  ip = ipnsegcn(iseg, ic)
                  do ipt = 1, npt
                     if( any(lcmplx_ipjp( ipnpt(ipt):(ipnpt(ipt) + nppt(ipt) - 1), ip ))) then !if any particle of type ipt is complexed with particle ip
                        ivar = ivarictipt(ict,ipt)
                        var(ivar)%avs2(iseg-1) = var(ivar)%avs2(iseg-1) + 1.0d0
                     end if
                  end do
               end do
            end do

         case (iAfterMacrostep)

            do ivar = 1, nvar
               var(ivar)%avs2(:) = var(ivar)%avs2(:)*var(ivar)%norm
            end do

            call DistFuncSample(iStage, nvar, var)
            if (lsim .and. master) write(ucnf) var

         case (iAfterSimulation)

         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

         deallocate(var)

      end select

      if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   end subroutine SegmentComplex

   !************************************************************************
   !*                                                                      *
   !*     ClusterDistribution
   !*                                                                      *
   !************************************************************************

   ! ... Calculate clusters and their distribution (first and second moment)

   subroutine ClusterDF(iStage)

      use MolModule, only: ltrace, ltime, uout, lsim, master, txstart, ucnf, ulist, ishow, iplot, ilist
      use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
      use MolModule, only: npt, txpt, nppt, ipnpt, np
      use StatisticsModule, only: df_var, mnbin_df
      implicit none

      integer(4), intent(in) :: iStage

      character(40), parameter :: txroutine ='ClusterDF'
      character(80), parameter :: txheading ='Size distribution of the binary clusters'
      integer(4),    save :: nvar
      integer(4)   , parameter :: nmoment = 2
      type(df_var),  allocatable, save              :: var(:,:)
      integer, allocatable, save              :: ivar_ptpt(:,:)

      logical(4), allocatable, save :: linclstr(:)        !logical if a particle in a cluster
      integer(4)  :: nbead_clstr                    !number of clusters

      integer(4)  :: ip, ipt, jpt, ivar, ibin, imoment

      if (ltrace) call WriteTrace(3, txroutine, iStage)
      if (ltime) call CpuAdd('start', txroutine, 1, uout)

      select case (iStage)
      case (iReadInput)
         if(.not. allocated(ivar_ptpt)) then
            allocate(ivar_ptpt(npt,npt))
         end if
         ivar = 0
         do ipt = 1, npt
            do jpt = ipt + 1, npt
               ivar = ivar + 1
               ivar_ptpt(ipt,jpt) = ivar
            end do
         end do

         nvar = ivar_ptpt(npt-1, npt)
         if(.not. allocated(var)) then
            allocate(var(nvar,nmoment))
         end if

      case (iWriteInput)

         do imoment = 1, nmoment
            do ipt = 1, npt
               ivar = ivar_ptpt(ipt,ipt)
               do jpt = ipt + 1, npt
                  ivar = ivar_ptpt(ipt,jpt)
                  write(var(ivar,imoment)%label,'(4a,i1)') trim(txpt(ipt)), ' - ', trim(txpt(jpt)), '; Moment = ', imoment
                  var(ivar,imoment)%min = -0.5d0
                  var(ivar,imoment)%max = nppt(ipt) + nppt(jpt) + 0.5
                  var(ivar,imoment)%nbin = min(mnbin_df,nppt(ipt) + nppt(jpt) + 1)
                  var(ivar,imoment)%norm = 1.0d0
               end do
            end do
         end do

         do imoment = 1, nmoment
            call DistFuncSample(iStage, nvar, var(:,imoment))
         end do

      case (iBeforeSimulation)

         do imoment = 1, nmoment
            call DistFuncSample(iStage, nvar, var(:,imoment))
         end do
         if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var
         if(.not. allocated(linclstr)) then
            allocate(linclstr(np))
            linclstr = .false.
         endif

      case (iBeforeMacrostep)

         call DistFuncSample(iStage, nvar, var)
         do imoment = 1, nmoment
            call distfuncsample(iStage, nvar, var(:,imoment))
         end do

      case (iSimulationStep)

         var%nsamp2 = var%nsamp2 + 1

         do ipt = 1, npt
            do jpt = ipt + 1, npt
               ivar = ivar_ptpt(ipt, jpt)

               !initialize linclstr
               linclstr = .true.
               linclstr((ipnpt(ipt)):(ipnpt(ipt)+nppt(ipt)-1)) = .false.
               linclstr(ipnpt(jpt):(ipnpt(jpt)+nppt(jpt)-1)) = .false.

               !count complexes (it is sufficient to count the ones where ipt is inside, the rest un in a cluster of size 1
               do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
                  if(.not. linclstr(ip)) then
                     nbead_clstr = 0
                     call get_nbead_clstr(ip,ipt, jpt,nbead_clstr)
                     if(nbead_clstr < 1) then
                        call Stop(txroutine, 'found particle which has wrong cluster size', uout)
                     endif
                     do imoment = 1, nmoment
                        ibin = max(-1,min(floor(var(ivar,imoment)%bini*(nbead_clstr-var(ivar,imoment)%min)),int(var(ivar,imoment)%nbin)))
                        var(ivar,imoment)%avs2(ibin) = var(ivar,imoment)%avs2(ibin) + nbead_clstr**(imoment - 1)
                     end do
                  end if
               end do

               !all particles of jpt which are not in a cluster with particles of ipt are individial particles (ibin 1)
               var(ivar,:)%avs2(1) = var(ivar,:)%avs2(1) + count(.not. linclstr(ipnpt(jpt):(ipnpt(jpt)+nppt(jpt)-1)))

            end do
         end do

      case (iAfterMacrostep)

         do imoment = 1, nmoment
            call DistFuncSample(iStage, nvar, var(:,imoment))
         end do
         if (lsim .and. master) write(ucnf) var

      case (iAfterSimulation)

         do imoment = 1, nmoment
            call DistFuncSample(iStage, nvar, var(:,imoment))
         end do
         call WriteHead(2, txheading, uout)
         do imoment = 1, nmoment
            call DistFuncHead(nvar, var(:,imoment), uout)
            call DistFuncWrite(txheading, nvar, var(:,imoment), uout, ulist, ishow, iplot, ilist)
         end do

         deallocate(var, linclstr)

      end select

      if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   contains

      recursive subroutine get_nbead_clstr(ip, iptc, jptc, nbead)

         implicit none

         integer(4), intent(in)  :: ip
         integer(4), intent(in)  :: iptc
         integer(4), intent(in)  :: jptc
         integer(4), intent(inout)  :: nbead

         integer(4)  :: jp

         if(linclstr(ip)) then
            return
         else
            nbead = nbead + 1
            linclstr(ip) = .true.
            do jp = ipnpt(jptc), ipnpt(jptc) + nppt(jptc) - 1
               if( (.not. linclstr(jp)) .and. (lcmplx_ipjp(ip,jp))) then
                     call get_nbead_clstr(jp, jptc, iptc, nbead)
               endif
            end do
         end if


      end subroutine get_nbead_clstr

   end subroutine ClusterDF

end module ComplexationModule

!*******************************************************************
!> \page moluser moluser.F90
!! **DoComplexation**
!! *documentation_missing*
!*******************************************************************

subroutine DoComplexation(iStage)
   use ComplexationModule, only: ComplexationDriver
   implicit none
   integer(4), intent(in)  :: iStage      ! event of SSO-Move
   call ComplexationDriver(iStage)
end subroutine
