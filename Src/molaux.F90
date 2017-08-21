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

module MolauxModule

   implicit none
   private
   public CalcCOM

   contains

      !************************************************************************
      !*                                                                      *
      !*     CalcCOM                                                              *
      !*                                                                      *
      !************************************************************************

      ! ... Calculate the Center of Mass of given coordinates.
      !     If MASK is given, use only the coordinates for which MASK is true

      function CalcCOM(r,MASK, MASS) result(com)
         implicit none
         real(8), intent(in), dimension (:,:)   :: r
         logical, intent(in), optional, dimension (:)   :: MASK
         real(8), intent(in), optional, dimension (:)   :: MASS
         real(8)   :: com(3)

         logical, dimension (size(r,DIM=2))   :: l
         real(8), dimension (size(r,DIM=2))   :: m

         integer(4)  :: i, n, first
         real(8)  :: d(3)

         if(present(MASK)) then
            l = MASK
         else
            l = .true.
         end if
         if(present(MASS)) then
            m = MASS
         else
            m = 1.0d0
         end if
         n = size(l)
         com = 0.0d0

         do i = 1, n   ! locate first particle which is to be used
            if(l(i)) then
               first = i
               exit
            end if
         end do

         do i = 1, n
            if (.not. l(i)) cycle
            d(1:3) = r(1:3,i) - r(1:3,first)
            call PBC(d(1),d(2),d(3))
            com(1:3) = com(1:3) + m(i)*d(1:3)
         end do
         com(1:3) = com(1:3)/(sum(m,MASK=l)) + r(1:3,first)

         call PBC(com(1),com(2),com(3))

      end function CalcCOM

end module MolauxModule

!************************************************************************
!*                                                                      *
!*     PBC                                                              *
!*                                                                      *
!************************************************************************

! ... apply periodic boundary conditions

subroutine PBC(dx,dy,dz)

   use MolModule
   implicit none

   real(8), intent(inout) :: dx, dy, dz

   if (lPBC) then                                                              ! periodic boundary condition
      if (lbcbox) then                                                         ! box-like cell
         ! if (abs(dx) > boxlen2(1)) dx = dx - sign(dpbc(1),dx)
         ! if (abs(dy) > boxlen2(2)) dy = dy - sign(dpbc(2),dy)
         ! if (abs(dz) > boxlen2(3)) dz = dz - sign(dpbc(3),dz)
         if (abs(dx) > boxlen2(1)) dx = dx - dpbc(1) * ANINT(dx*boxleni(1))
         if (abs(dy) > boxlen2(2)) dy = dy - dpbc(2) * ANINT(dy*boxleni(2))
         if (abs(dz) > boxlen2(3)) dz = dz - dpbc(3) * ANINT(dz*boxleni(3))
      else if (lbcrd) then                                                     ! rhombic dodecahedral cell
         if (abs(dx) > boxlen2(1)) dx = dx - boxlen(1) * ANINT(dx*boxleni(1))
         if (abs(dy) > boxlen2(2)) dy = dy - boxlen(2) * ANINT(dy*boxleni(2))
         if (abs(dz) > boxlen2(3)) dz = dz - boxlen(3) * ANINT(dz*boxleni(3))
         if (abs(dx) + abs(dy) + SqTwo*abs(dz) > boxlen(1)) then
            dx = dx - sign(boxlen2(1),dx)
            dy = dy - sign(boxlen2(2),dy)
            dz = dz - sign(boxlen2(3),dz)
         end if
      else if (lbcto) then                                                     ! truncated octahedral cell
         if (abs(dx) > boxlen2(1)) dx = dx - boxlen(1) * ANINT(dx*boxleni(1))
         if (abs(dy) > boxlen2(2)) dy = dy - boxlen(2) * ANINT(dy*boxleni(2))
         if (abs(dz) > boxlen2(3)) dz = dz - boxlen(3) * ANINT(dz*boxleni(3))
         if (abs(dx) + abs(dy) + abs(dz) > ThreeHalf*boxlen2(1)) then
            dx = dx - sign(boxlen2(1),dx)
            dy = dy - sign(boxlen2(2),dy)
            dz = dz - sign(boxlen2(3),dz)
         end if
      end if
   end if

end subroutine PBC

!************************************************************************
!*                                                                      *
!*     PBCr2                                                            *
!*                                                                      *
!************************************************************************

! ... apply periodic boundary conditions and calculate r**2

subroutine PBCr2(dx,dy,dz,r2)

   use MolModule
   implicit none

   real(8), intent(inout) :: dx, dy, dz
   real(8), intent(out) :: r2

   if (lPBC) then                                                              ! periodic boundary condition
      if (lbcbox) then                                                         ! box-like cell
         if (abs(dx) > boxlen2(1)) dx = dx - sign(dpbc(1),dx)
         if (abs(dy) > boxlen2(2)) dy = dy - sign(dpbc(2),dy)
         if (abs(dz) > boxlen2(3)) dz = dz - sign(dpbc(3),dz)
      else if (lbcrd) then                                                     ! rhombic dodecahedral cell
         if (abs(dx) > boxlen2(1)) dx = dx - sign(boxlen(1),dx)
         if (abs(dy) > boxlen2(2)) dy = dy - sign(boxlen(2),dy)
         if (abs(dz) > boxlen2(3)) dz = dz - sign(boxlen(3),dz)
         if (abs(dx) + abs(dy) + SqTwo*abs(dz) > boxlen(1)) then
            dx = dx - sign(boxlen2(1),dx)
            dy = dy - sign(boxlen2(2),dy)
            dz = dz - sign(boxlen2(3),dz)
         end if
      else if (lbcto) then                                                     ! truncated octahedral cell
         if (abs(dx) > boxlen2(1)) dx = dx - sign(boxlen(1),dx)
         if (abs(dy) > boxlen2(2)) dy = dy - sign(boxlen(2),dy)
         if (abs(dz) > boxlen2(3)) dz = dz - sign(boxlen(3),dz)
         if (abs(dx) + abs(dy) + abs(dz) > ThreeHalf*boxlen2(1)) then
            dx = dx - sign(boxlen2(1),dx)
            dy = dy - sign(boxlen2(2),dy)
            dz = dz - sign(boxlen2(3),dz)
         end if
      end if
   end if
   r2 = dx**2+dy**2+dz**2

end subroutine PBCr2

!************************************************************************
!*                                                                      *
!*     PBC2                                                             *
!*                                                                      *
!************************************************************************

! ... calculate periodic boundary condition displacement

subroutine PBC2(dx,dy,dz,dxpbc,dypbc,dzpbc)

   use MolModule
   implicit none

   real(8), intent(in)  :: dx, dy, dz
   real(8), intent(out) :: dxpbc, dypbc, dzpbc

   real(8)              :: dxaux, dyaux, dzaux

   if (lPBC) then                                                          ! periodic boundary condition
      if (lbcbox) then                                                     ! box-like cell
         dxpbc = Zero
         dypbc = Zero
         dzpbc = Zero
         if (abs(dx) > boxlen2(1)) dxpbc = sign(dpbc(1),dx)
         if (abs(dy) > boxlen2(2)) dypbc = sign(dpbc(2),dy)
         if (abs(dz) > boxlen2(3)) dzpbc = sign(dpbc(3),dz)
      else if (lbcrd) then                                                 ! rhombic dodecahedral cell
         dxaux = dx
         dyaux = dy
         dzaux = dz
         if (abs(dxaux) > boxlen2(1)) dxaux = dxaux - sign(boxlen(1),dxaux)
         if (abs(dyaux) > boxlen2(2)) dyaux = dyaux - sign(boxlen(2),dyaux)
         if (abs(dzaux) > boxlen2(3)) dzaux = dzaux - sign(boxlen(3),dzaux)
         if (abs(dxaux) + abs(dyaux) + SqTwo*abs(dzaux) > boxlen(1)) then
            dxaux = dxaux - sign(boxlen2(1),dxaux)
            dyaux = dyaux - sign(boxlen2(2),dyaux)
            dzaux = dzaux - sign(boxlen2(3),dzaux)
         end if
         dxpbc = dx - dxaux
         dypbc = dy - dyaux
         dzpbc = dz - dzaux
      else if (lbcto) then                                                 ! truncated octahedral cell
         dxaux = dx
         dyaux = dy
         dzaux = dz
         if (abs(dxaux) > boxlen2(1)) dxaux = dxaux - sign(boxlen(1),dxaux)
         if (abs(dyaux) > boxlen2(2)) dyaux = dyaux - sign(boxlen(2),dyaux)
         if (abs(dzaux) > boxlen2(3)) dzaux = dzaux - sign(boxlen(3),dzaux)
         if (abs(dxaux) + abs(dyaux) + abs(dzaux) > ThreeHalf*boxlen2(1)) then
            dxaux = dxaux - sign(boxlen2(1),dxaux)
            dyaux = dyaux - sign(boxlen2(2),dyaux)
            dzaux = dzaux - sign(boxlen2(3),dzaux)
         end if
         dxpbc = dx - dxaux
         dypbc = dy - dyaux
         dzpbc = dz - dzaux
      end if
   else                                                                 ! non-periodic boundary condition
      dxpbc = Zero
      dypbc = Zero
      dzpbc = Zero
   end if

end subroutine PBC2

!************************************************************************
!*                                                                      *
!*     FileOpen                                                         *
!*                                                                      *
!************************************************************************

! ... open external files

subroutine FileOpen(unit, fname, txopt)

   implicit none

   integer(4),   intent(in) :: unit           ! unit number
   character(*), intent(in) :: fname          ! file name
   character(*), intent(in) :: txopt          ! mode of action

   integer(4) :: ios

   if (txopt == 'unform/noread') then
      open(unit, file = fname, form = 'unformatted', status = 'unknown', err = 999, iostat = ios)
   else if (txopt == 'unform/read') then
      open(unit, file = fname, form = 'unformatted', status = 'unknown', err = 999, iostat = ios)
10    read(unit,end = 900)
      goto 10
   else if (txopt == 'form/noread') then
      open(unit, file = fname, form = 'formatted', status = 'unknown', err = 999, iostat = ios)
   else if (txopt == 'form/read') then
      open(unit, file = fname, form = 'formatted', status = 'unknown', err = 999, iostat = ios)
20    read(unit,*,end = 900)
      goto 20
   else
      call Warn('FileOpen', 'unsupported value of txopt', 6)
   end if
999 if (ios/= 0) then
      call Warn('FileOpen', ' ', 6)
      write(*,*) 'unit = ', unit, '  fname = ', fname, '  ios = ', ios
   end if
900 continue

end subroutine FileOpen

!************************************************************************
!*                                                                      *
!*     FileFlush                                                        *
!*                                                                      *
!************************************************************************

! ... flush external files   (contains a system call)

subroutine FileFlush(unit)
   implicit none
   integer(4),   intent(in) :: unit
   call flush(unit)
end subroutine FileFlush

!************************************************************************
!*                                                                      *
!*     SetBoxParam                                                      *
!*                                                                      *
!************************************************************************

! ... set box related parameters

!     periodical boundary conditions should be applied as
!     if (abs(dx) > boxlen2(1)) dx = dx - sign(dpbc(1),dx)
!     if (abs(dy) > boxlen2(2)) dy = dy - sign(dpbc(2),dy)
!     if (abs(dx) > boxlen2(3)) dz = dz - sign(dpbc(3),dz)

subroutine SetBoxParam

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='SetBoxParam'

   lbcbox = .false.
   lbcrd  = .false.
   lbcto  = .false.
   lbcsph = .false.
   lbccyl = .false.
   lbcell = .false.
   lPBC   = .false.

   if ((txbc == 'xyz') .or. (txbc == 'xy') .or. (txbc == 'z')) then    ! box
      lbcbox = .true.
      boxlen2 = boxlen/Two
      boxleni = One/boxlen
      boxlen2i = One/boxlen2
      TwoPiBoxi = TwoPi*boxleni
      vol = boxlen(1)*boxlen(2)*boxlen(3)
      if (txbc == 'xyz') then                                          ! pbc in x-, y-, and z-dir
         dpbc = boxlen
      else if (txbc == 'xy') then                                      ! pbc in x- and y-dir
         dpbc(1:2) = boxlen(1:2)
         dpbc(3) = Zero
      else if (txbc == 'z') then                                       ! pbc in z-dir
         dpbc(1:2) = Zero
         dpbc(3) = boxlen(3)
      end if
   else if (txbc == 'rd') then                                         ! rhombic dodecahedral cell
      lbcrd  = .true.
      vol = (16.0/9.0)*sqrt(Three)*cellside**3
      boxlen(1:2) = Two*sqrt(Two/Three)*cellside                       ! sides of box inscribing simulation cell
      boxlen(3)   = SqTwo*boxlen(1)
      boxlen2 = boxlen/Two
      boxleni = One/boxlen
      boxlen2i = One/boxlen2
      TwoPiBoxi = TwoPi*boxleni
   else if (txbc == 'to') then                                         ! truncated octahedral cell
      lbcto  = .true.
      vol = 8.0*SqTwo*cellside**3
      boxlen = Two*SqTwo*cellside                                      ! sides of box inscribing simulation cell
      boxlen2 = boxlen/Two
      boxleni = One/boxlen
      boxlen2i = One/boxlen2
      TwoPiBoxi = TwoPi*boxleni
   else if (txbc == 'sph') then                                        ! spherical cell
      lbcsph = .true.
      sphrad2 = sphrad**2
      vol = FourPiThird*sphrad**3
      dpbc = Zero
   else if (txbc == 'cyl') then                                        ! cylindrical cell
      lbccyl = .true.
      cylrad2 = cylrad**2
      vol = Pi*cylrad2*cyllen
      dpbc = Zero
   else if (txbc == 'ell') then                                        ! ellipsoidal cell
      lbcell = .true.
      ellradi = One/ellrad
      vol = FourPiThird*ellrad(1)*ellrad(2)*ellrad(3)
   else
      call Stop(txroutine, 'wrong option of txbc', uout)
   end if

   if (lbcbox .or. lbcrd .or. lbcto) lPBC = .true.        ! periodic boundary condition(s)
   if (lPBC) boxlenshort = minval(boxlen(1:3))            ! shortest boxlen(1:3)

end subroutine SetBoxParam

!************************************************************************
!*                                                                      *
!*     SetPartOriRandom                                                 *
!*                                                                      *
!************************************************************************

! ... set particle orientation matrix corresponding to a random orientation

subroutine SetPartOriRandom(iseed, ori)

   implicit none

   integer(4), intent(inout)  :: iseed         ! seed of random number generator
   real(8),    intent(out)    :: ori(1:3,1:3)  ! oritentation matrix

   real(8), parameter :: half = 0.5d0
   real(8) :: Random

   ori(1,1) = Random(iseed)-half                                       ! x'-axis
   ori(2,1) = Random(iseed)-half
   ori(3,1) = Random(iseed)-half

   ori(1,2) = Random(iseed)-half                                       ! y'-axis
   ori(2,2) = Random(iseed)-half
   ori(3,2) =-(ori(1,1)*ori(1,2)+ori(2,1)*ori(2,2))/ori(3,1)

   ori(1,3) = ori(2,1)*ori(3,2) - ori(3,1)*ori(2,2)                    ! z'-axis
   ori(2,3) = ori(3,1)*ori(1,2) - ori(1,1)*ori(3,2)
   ori(3,3) = ori(1,1)*ori(2,2) - ori(2,1)*ori(1,2)

   ori(1:3,1) = ori(1:3,1)/sqrt(ori(1,1)**2+ori(2,1)**2+ori(3,1)**2)   ! normalize
   ori(1:3,2) = ori(1:3,2)/sqrt(ori(1,2)**2+ori(2,2)**2+ori(3,2)**2)
   ori(1:3,3) = ori(1:3,3)/sqrt(ori(1,3)**2+ori(2,3)**2+ori(3,3)**2)
end subroutine SetPartOriRandom

!************************************************************************
!*                                                                      *
!*     SetPartOriLab                                                    *
!*                                                                      *
!************************************************************************

! ... set particle orientation matrix to match lab frame orientation

subroutine SetPartOriLab(ori)
   implicit none
   real(8),    intent(out)    :: ori(3,3)  ! oritentation matrix

   ori(1,1) = 1.0d0
   ori(2,1) = 0.0d0
   ori(3,1) = 0.0d0

   ori(1,2) = 0.0d0
   ori(2,2) = 1.0d0
   ori(3,2) = 0.0d0

   ori(1,3) = 0.0d0
   ori(2,3) = 0.0d0
   ori(3,3) = 1.0d0

end subroutine SetPartOriLab

!************************************************************************
!*                                                                      *
!*     SetPartOriBond                                                   *
!*                                                                      *
!************************************************************************

! ... set particle orientation matrix to match bond directions of the particle

!     zero angle problem remains

subroutine SetPartOriBond(txloc, r21, r23, ori)

   implicit none

   character(40), parameter :: txroutine ='SetPartOriBond'
   character(3), intent(in)    :: txloc     ! 'mid' => noend particle, 'end' => end particle
   real(8),      intent(in)    :: r21(3)    ! bond vector
   real(8),      intent(in)    :: r23(3)    ! bond vector
   real(8),      intent(out)   :: ori(3,3)  ! oritentation matrix
   real(8) :: r21n(3), r23n(3)

   if (txloc == 'end') then           ! end particle

! ... set z' axis in the direction of the bond vector

      ori(1,3) = r21(1)
      ori(2,3) = r21(2)
      ori(3,3) = r21(3)
      ori(1:3,3) = ori(1:3,3)/sqrt(ori(1,3)**2+ori(2,3)**2+ori(3,3)**2)             ! normalize

! ... y' = rz x (0,1,0)

      ori(1,2) = -ori(3,3)
      ori(2,2) = 0
      ori(3,2) = ori(1,3)
      ori(1:3,2) = ori(1:3,2)/sqrt(ori(1,2)**2+ori(2,2)**2+ori(3,2)**2)             ! normalize

! ... x' = rz x ry

      ori(1,1) = ori(2,3)*ori(3,2) - ori(3,3)*ori(2,2)
      ori(2,1) = ori(3,3)*ori(1,2) - ori(1,3)*ori(3,2)
      ori(3,1) = ori(1,3)*ori(2,2) - ori(2,3)*ori(1,2)
      ori(1:3,1) = ori(1:3,1)/sqrt(ori(1,1)**2+ori(2,1)**2+ori(3,1)**2)             ! normalize

   else if (txloc == 'mid') then     ! non-end particle

      r21n = r21/sqrt(sum(r21**2))   ! normalize bond vector
      r23n = r23/sqrt(sum(r23**2))   ! normalize bond vector

      ori(1,1) = r21n(1)+r23n(1)
      ori(2,1) = r21n(2)+r23n(2)
      ori(3,1) = r21n(3)+r23n(3)

      ori(1,2) = r21n(2)*r23n(3)-r21n(3)*r23n(2)
      ori(2,2) = r21n(3)*r23n(1)-r21n(1)*r23n(3)
      ori(3,2) = r21n(1)*r23n(2)-r21n(2)*r23n(1)

      ori(1,3) = ori(2,1)*ori(3,2) - ori(3,1)*ori(2,2)                    ! z'-axis
      ori(2,3) = ori(3,1)*ori(1,2) - ori(1,1)*ori(3,2)
      ori(3,3) = ori(1,1)*ori(2,2) - ori(2,1)*ori(1,2)

      ori(1:3,1) = ori(1:3,1)/sqrt(ori(1,1)**2+ori(2,1)**2+ori(3,1)**2)   ! normalize
      ori(1:3,2) = ori(1:3,2)/sqrt(ori(1,2)**2+ori(2,2)**2+ori(3,2)**2)
      ori(1:3,3) = ori(1:3,3)/sqrt(ori(1,3)**2+ori(2,3)**2+ori(3,3)**2)

   else

      call Stop(txroutine, 'error in txloc', 6)

   end if

end subroutine SetPartOriBond

!************************************************************************
!*                                                                      *
!*     SetAtomPos                                                       *
!*                                                                      *
!************************************************************************

! ... set atom positions in the laboratory frame

!  r = ro + ori . ra

!  r(1,ia)   ro(1,ip)   ori(1,1,ip) ori(1,2,ip) ori(1,3,ip)   ra(1,ialoc,ipt)
!  r(2,ia) = ro(2,ip) + ori(2,1,ip) ori(2,2,ip) ori(2,3,ip) x ra(2,ialoc,ipt)
!  r(3,ia)   ro(3,ip)   ori(3,1,ip) ori(3,2,ip) ori(3,3,ip)   ra(3,ialoc,ipt)

subroutine SetAtomPos(iplow, ipupp, lint)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp
   logical,    intent(in) :: lint

   integer(4) :: ip, ipt, ia, ialoc, ialow, iaupp

   ialow = ianpn(iplow)
   iaupp = ianpn(ipupp)+napt(iptpn(ipupp))-1

   if (.not.lint) then   ! use ra

      do ia = ialow, iaupp
         ip = ipnan(ia)
         ipt = iptan(ia)
         ialoc = ia-(ianpn(ip)-1)
         r(1,ia) = ro(1,ip)+ori(1,1,ip)*ra(1,ialoc,ipt)+ori(1,2,ip)*ra(2,ialoc,ipt)+ori(1,3,ip)*ra(3,ialoc,ipt)
         r(2,ia) = ro(2,ip)+ori(2,1,ip)*ra(1,ialoc,ipt)+ori(2,2,ip)*ra(2,ialoc,ipt)+ori(2,3,ip)*ra(3,ialoc,ipt)
         r(3,ia) = ro(3,ip)+ori(3,1,ip)*ra(1,ialoc,ipt)+ori(3,2,ip)*ra(2,ialoc,ipt)+ori(3,3,ip)*ra(3,ialoc,ipt)
      end do

   else                 ! use rasite

      do ia = ialow, iaupp
         ip = ipnan(ia)
         ipt = iptan(ia)
         ialoc = ia-(ianpn(ip)-1)
         r(1,ia) = ro(1,ip)+ori(1,1,ip)*rasite(1,ialoc,ipt)+ori(1,2,ip)*rasite(2,ialoc,ipt)+ori(1,3,ip)*rasite(3,ialoc,ipt)
         r(2,ia) = ro(2,ip)+ori(2,1,ip)*rasite(1,ialoc,ipt)+ori(2,2,ip)*rasite(2,ialoc,ipt)+ori(2,3,ip)*rasite(3,ialoc,ipt)
         r(3,ia) = ro(3,ip)+ori(3,1,ip)*rasite(1,ialoc,ipt)+ori(3,2,ip)*rasite(2,ialoc,ipt)+ori(3,3,ip)*rasite(3,ialoc,ipt)
      end do

   end if

end subroutine SetAtomPos

!************************************************************************
!*                                                                      *
!*     SetAtomDipMom                                                    *
!*                                                                      *
!************************************************************************

! ... set atom dipole moments in the laboratory frame

!  dip =  ori . dipa

!  dip(1,ia)   ori(1,1,ip) ori(1,2,ip) ori(1,3,ip)   dipa(1,ialoc,ipt)
!  dip(2,ia) = ori(2,1,ip) ori(2,2,ip) ori(2,3,ip) x dipa(2,ialoc,ipt)
!  dip(3,ia)   ori(3,1,ip) ori(3,2,ip) ori(3,3,ip)   dipa(3,ialoc,ipt)

subroutine SetAtomDipMom(iplow, ipupp)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp

   integer(4) :: ip, ipt, ia, ialoc, ialow, iaupp

   ialow = ianpn(iplow)
   iaupp = ianpn(ipupp)+napt(iptpn(ipupp))-1
   do ia = ialow, iaupp
      ip = ipnan(ia)
      ipt = iptan(ia)
      ialoc = ia-(ianpn(ip)-1)
      dip(1,ia) = ori(1,1,ip)*dipa(1,ialoc,ipt)+ori(1,2,ip)*dipa(2,ialoc,ipt)+ori(1,3,ip)*dipa(3,ialoc,ipt)
      dip(2,ia) = ori(2,1,ip)*dipa(1,ialoc,ipt)+ori(2,2,ip)*dipa(2,ialoc,ipt)+ori(2,3,ip)*dipa(3,ialoc,ipt)
      dip(3,ia) = ori(3,1,ip)*dipa(1,ialoc,ipt)+ori(3,2,ip)*dipa(2,ialoc,ipt)+ori(3,3,ip)*dipa(3,ialoc,ipt)
   end do

end subroutine SetAtomDipMom

!************************************************************************
!*                                                                      *
!*     SetAtomPolTens                                                   *
!*                                                                      *
!************************************************************************

! ... set atom polarization tensors in the laboratory frame

!                                T
!  poltens =  ori . poltens . ori

subroutine SetAtomPolTens(iplow, ipupp)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp

   integer(4) :: ip, ipt, ia, ialoc, ialow, iaupp, m, n
   real(8)    :: atemp(3,3), btemp(3,3)

   ialow = ianpn(iplow)
   iaupp = ianpn(ipupp)+napt(iptpn(ipupp))-1

   do ia = ialow, iaupp
      ip = ipnan(ia)
      ipt = iptan(ia)
      ialoc = ia-(ianpn(ip)-1)

      atemp(1,1) = poltensa(1,ialoc,ipt)
      atemp(1,2) = poltensa(4,ialoc,ipt)
      atemp(1,3) = poltensa(5,ialoc,ipt)
      atemp(2,1) = poltensa(4,ialoc,ipt)
      atemp(2,2) = poltensa(2,ialoc,ipt)
      atemp(2,3) = poltensa(6,ialoc,ipt)
      atemp(3,1) = poltensa(5,ialoc,ipt)
      atemp(3,2) = poltensa(6,ialoc,ipt)
      atemp(3,3) = poltensa(3,ialoc,ipt)

      do m = 1, 3
         btemp(m,1) = ori(m,1,ip)*atemp(1,1)+ori(m,2,ip)*atemp(2,1)+ori(m,3,ip)*atemp(3,1)
         btemp(m,2) = ori(m,1,ip)*atemp(1,2)+ori(m,2,ip)*atemp(2,2)+ori(m,3,ip)*atemp(3,2)
         btemp(m,3) = ori(m,1,ip)*atemp(1,3)+ori(m,2,ip)*atemp(2,3)+ori(m,3,ip)*atemp(3,3)
      end do

      do n = 1, 3
         atemp(1,n) = btemp(1,1)*ori(n,1,ip)+btemp(1,2)*ori(n,2,ip)+btemp(1,3)*ori(n,3,ip)
         atemp(2,n) = btemp(2,1)*ori(n,1,ip)+btemp(2,2)*ori(n,2,ip)+btemp(2,3)*ori(n,3,ip)
         atemp(3,n) = btemp(3,1)*ori(n,1,ip)+btemp(3,2)*ori(n,2,ip)+btemp(3,3)*ori(n,3,ip)
      end do

      poltens(1,ia) = atemp(1,1)
      poltens(2,ia) = atemp(2,2)
      poltens(3,ia) = atemp(3,3)
      poltens(4,ia) = atemp(1,2)
      poltens(5,ia) = atemp(1,3)
      poltens(6,ia) = atemp(2,3)

   end do

end subroutine SetAtomPolTens

!************************************************************************
!*                                                                      *
!*     CalcIndDipMom                                                    *
!*                                                                      *
!************************************************************************

! ... calculate induced dipole moment: total and particle

subroutine CalcIndDipMom

   use MolModule
   implicit none

   integer(4) :: ip, ipt, ia, ialoc, ialow

   do ip = 1, np
      ipt = iptpn(ip)
      ialow = ianpn(ip)-1
      idmo(1:3,ip) = Zero
      do ialoc = 1, napt(ipt)
         ia = ialow+ialoc
         idmo(1,ip) = idmo(1,ip) + idm(1,ia)
         idmo(2,ip) = idmo(2,ip) + idm(2,ia)
         idmo(3,ip) = idmo(3,ip) + idm(3,ia)
      end do
   end do
   idmsys(1) = sum(idmo(1,1:np))
   idmsys(2) = sum(idmo(2,1:np))
   idmsys(3) = sum(idmo(3,1:np))

end subroutine CalcIndDipMom

!************************************************************************
!*                                                                      *
!*     CalcSysDipMom                                                    *
!*                                                                      *
!************************************************************************

! ... calculate system dipole moment from charges, static, and induced dipoles

subroutine CalcSysDipMom(dipsys)

   use MolModule
   implicit none

   real(8), intent(out) :: dipsys(3)
   integer(4) :: ia

   dipsys(1:3) = Zero
   if (lcharge) then
      do ia = 1, na
         dipsys(1) = dipsys(1) + az(ia)*r(1,ia)
         dipsys(2) = dipsys(2) + az(ia)*r(2,ia)
         dipsys(3) = dipsys(3) + az(ia)*r(3,ia)
      end do
   else if (ldipole .or. ldipolesph) then
      do ia = 1, na
         dipsys(1) = dipsys(1) + az(ia)*r(1,ia) + dip(1,ia)
         dipsys(2) = dipsys(2) + az(ia)*r(2,ia) + dip(2,ia)
         dipsys(3) = dipsys(3) + az(ia)*r(3,ia) + dip(3,ia)
      end do
   else if (lpolarization) then
      do ia = 1, na
         dipsys(1) = dipsys(1) + az(ia)*r(1,ia) + diptot(1,ia)
         dipsys(2) = dipsys(2) + az(ia)*r(2,ia) + diptot(2,ia)
         dipsys(3) = dipsys(3) + az(ia)*r(3,ia) + diptot(3,ia)
      end do
   end if

end subroutine CalcSysDipMom

!************************************************************************
!*                                                                      *
!*     CalcPartDipMom                                                   *
!*                                                                      *
!************************************************************************

! ... calculate particle dipole moment from charges, static, and induced dipoles

subroutine CalcPartDipMom(dipo)

   use MolModule
   implicit none

   integer(4) :: ip, ipt, ia, ialoc, ialow
   real(8)    :: dipo(3,np)

   do ip = 1, np
      ipt = iptpn(ip)
      ialow = ianpn(ip)-1
      dipo(1:3,ip) = Zero
      if (lcharge) then
         do ialoc = 1, napt(ipt)
            ia = ialow+ialoc
            dipo(1,ip) = dipo(1,ip) + az(ia)*r(1,ia)
            dipo(2,ip) = dipo(2,ip) + az(ia)*r(2,ia)
            dipo(3,ip) = dipo(3,ip) + az(ia)*r(3,ia)
         end do
      else if (ldipole .or. ldipolesph) then
         do ialoc = 1, napt(ipt)
            ia = ialow+ialoc
            dipo(1,ip) = dipo(1,ip) + az(ia)*r(1,ia) + dip(1,ia)
            dipo(2,ip) = dipo(2,ip) + az(ia)*r(2,ia) + dip(2,ia)
            dipo(3,ip) = dipo(3,ip) + az(ia)*r(3,ia) + dip(3,ia)
         end do
      else if (lpolarization) then
         do ialoc = 1, napt(ipt)
            ia = ialow+ialoc
            dipo(1,ip) = dipo(1,ip) + az(ia)*r(1,ia) + diptot(1,ia)
            dipo(2,ip) = dipo(2,ip) + az(ia)*r(2,ia) + diptot(2,ia)
            dipo(3,ip) = dipo(3,ip) + az(ia)*r(3,ia) + diptot(3,ia)
         end do
      end if
   end do

end subroutine CalcPartDipMom

!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

!************************************************************************
!*                                                                      *
!*     dvol                                                             *
!*                                                                      *
!************************************************************************

! ... return the volume between two concentric shells, 4*pi/3 is not included

real(8) function dvol(i, rl, dr)
   implicit none
   integer(4), intent(in) :: i      ! lower radius: rl + i*dr,  upper radius: rl + (i+1)*dr
   real(8), intent(in)    :: rl
   real(8), intent(in)    :: dr
   dvol = ( ( (3*i*(i+1)+1)*dr + 3.0d0*(2*i+1)*rl )*dr + 3.0d0*rl**2 )*dr
   if (dvol <= 0.0d0) call Warn('dvol', 'dvol <= Zero', 6)
end function dvol

!************************************************************************
!*                                                                      *
!*     darea                                                            *
!*                                                                      *
!************************************************************************

! ... return the area between two concentric spheres, pi is not included

real(8) function darea(i, rl, dr)
   implicit none
   integer(4), intent(in) :: i      ! lower radius: rl + i*dr,  upper radius: rl + (i+1)*dr
   real(8), intent(in)    :: rl
   real(8), intent(in)    :: dr
   darea = ( (2*i+1)*dr + 2.0d0*rl )*dr
   if (darea <= 0.0d0) call Warn('darea', 'darea <= Zero', 6)
end function darea

!************************************************************************
!*                                                                      *
!*     angcos                                                           *
!*                                                                      *
!************************************************************************

! ... return the angle cosine between two 3d vectors

real(8) function angcos(x1, y1, z1, x2, y2, z2)
   implicit none
   real(8), intent(in) :: x1, y1, z1, x2, y2, z2
   real(8) :: term1, term2, term3
   term1 = x1*x2+y1*y2+z1*z2
   term2 = x1**2+y1**2+z1**2
   term3 = x2**2+y2**2+z2**2
   angcos = term1/sqrt(term2*term3)
end function angcos

!************************************************************************
!*                                                                      *
!*     angle_rad                                                        *
!*                                                                      *
!************************************************************************

! ... return the angle between two 3d vectors

real(8) function angle_rad(x1, y1, z1, x2, y2, z2)
   implicit none
   real(8), parameter :: One = 1.0d0
   real(8), intent(in) :: x1, y1, z1, x2, y2, z2
   real(8) :: term1, term2, term3, term
   term1 = x1*x2+y1*y2+z1*z2
   term2 = x1**2+y1**2+z1**2
   term3 = x2**2+y2**2+z2**2
   term = term1/sqrt(term2*term3)
   angle_rad = acos(max(-One,min(term,One)))
end function angle_rad

!************************************************************************
!*                                                                      *
!*     angle_deg                                                        *
!*                                                                      *
!************************************************************************

! ... return the angle between two 3d vectors

real(8) function angle_deg(x1, y1, z1, x2, y2, z2)
   implicit none
   real(8), parameter :: One = 1.0d0, radgrd = 57.29577951d0
   real(8), intent(in) :: x1, y1, z1, x2, y2, z2
   real(8) :: term1, term2, term3, term
   term1 = x1*x2+y1*y2+z1*z2
   term2 = x1**2+y1**2+z1**2
   term3 = x2**2+y2**2+z2**2
   term = term1/sqrt(term2*term3)
   angle_deg = radgrd*acos(max(-One,min(term,One)))
end function angle_deg

!************************************************************************
!*                                                                      *
!*     WritePartAtomProp                                                *
!*                                                                      *
!************************************************************************

! ... write particle and atom properties

subroutine WritePartAtomProp

   use MolModule
   implicit none

   integer(4) :: ip, ia

   if (slave) return

   if (ipart /= 0 .or. iatom /= 0) call WriteHead(2, 'coordinates', uout)

   if (ipart /= 0) then
      write(uout,'()')
      write(uout,'(a)') 'particle coordinates'
      write(uout,'(a)') '--------------------'
      write(uout,'(a,t15,a,t35,a,t60,a,t83,a,t106,a)')                  &
       'particle no', 'particle type', 'com (x,y,z)', 'x''-axis (x,y,z)', 'y''-axis (x,y,z)', 'z''-axis (x,y,z)'
      write(uout,'(a,t15,a,t35,a,t60,a,t83,a,t106,a)')                  &
       '-----------', '-------------', '-----------', '---------------', '---------------', '---------------'
      write(uout,'(i6,t15,a,t27,3f9.3,t55,3f7.3,t78,3f7.3,t101,3f7.3)') &
        (ip,txpt(iptpn(ip)),ro(1:3,ip),ori(1:3,1:3,ip),ip = 1,np,ipart)
   end if

   if (lmd .and. ipart /= 0) then
      write(uout,'()')
      write(uout,'(a)') 'particle velocities'
      write(uout,'(a)') '-------------------'
      write(uout,'(a,t15,a,t33,a,t60,a)')                                                 &
        'particle no', 'particle type', 'linear vel. (x,y,z)', 'angular vel. (x'',y'',z'')', &
        '-----------', '-------------', '-------------------', '-----------------------'
      write(uout,'(i6,t15,a,t28,3f9.3,t56,3f9.3)')                                        &
       (ip,txpt(iptpn(ip)),rod(1:3,ip),angvelo(1:3,ip),ip = 1,np,ipart)
   end if

   if (iatom /= 0) then
      write(uout,'()')
      write(uout,'(a)') 'atom coordinates'
      write(uout,'(a)') '----------------'
      write(uout,'(a,t15,a,t30,a,t45,a,t70,a,t85,a,t95,a)')                        &
        'atom no', 'atom type', 'particle no', 'particle type', 'coordinate (x,y,z)'
      write(uout,'(a,t15,a,t30,a,t45,a,t70,a,t85,a,t95,a)')                        &
        '-------', '---------', '-----------', '-------------', '------------------'
      write(uout,'(i6,t15,a,t30,i6,t45,a,t60,3f12.5)')                             &
       (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),r(1:3,ia),ia = 1,na,iatom)
   end if

    if (iatom /= 0 .and. lweakcharge) then
       write(uout,'()')
       write(uout,'(a)') 'atom charge'
       write(uout,'(a)') '-----------'
       write(uout,'(a,t15,a,t30,a)') 'atom no', 'weak charge', 'charged'
       write(uout,'(a,t15,a,t30,a)') '-------', '-----------', '-------'
       write(uout,'(i6,t10,2l12  )') (ia,latweakcharge(iatan(ia)),laz(ia),ia = 1,na,iatom)
    end if

    if (iatom /= 0 .and. (ldipole .or. ldipolesph)) then
       write(uout,'()')
       write(uout,'(a)') 'atom dipole moments'
       write(uout,'(a)') '-------------------'
       write(uout,'(a,t15,a,t55,a)') 'atom no', 'dipole moment (x,y,z)'
       write(uout,'(a,t15,a,t55,a)') '-------', '---------------------'
       write(uout,'(i6,t10,3f10.5)') (ia,dip(1:3,ia),ia = 1,na,iatom)
    end if

    if (iatom /= 0 .and. lpolarization) then
       write(uout,'()')
       write(uout,'(a)') 'atom dipole moments and polarizabilities'
       write(uout,'(a)') '----------------------------------------'
       write(uout,'(a,t15,a,t55,a)') 'atom no', 'dipole moment (x,y,z)', 'polarizability (xx,yy,zz,xy,xz,yz)'
       write(uout,'(a,t15,a,t55,a)') '-------', '---------------------', '----------------------------------'
       write(uout,'(i6,t10,3f10.5,5x,6f10.5)') (ia,dip(1:3,ia),poltens(1:6,ia),ia = 1,na,iatom)
    end if

   if ((iatom /= 0) .and. lpolarization .and. (istep1 > nstep1)) then   ! ind dip mom are not calculated before simulation
      write(uout,'()')
      write(uout,'(a)') 'induced atom dipole moments'
      write(uout,'(a)') '---------------------------'
      write(uout,'(a,t15,a,t30,a,t45,a,t70,a,t85,a,t95,a)')                             &
        'atom no', 'atom type', 'particle no', 'particle type', 'dipole moment (x,y,z)'
      write(uout,'(a,t15,a,t30,a,t45,a,t70,a,t85,a,t95,a)')                             &
        '-------', '---------', '-----------', '-------------', '---------------------'
      write(uout,'(i4,t15,a,t30,i6,t45,a,t60,3f12.5)')                                  &
       (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),idm(1:3,ia),ia = 1,na,iatom)
   end if

end subroutine WritePartAtomProp

!************************************************************************
!*                                                                      *
!*     lWarnHCOverlap                                                   *
!*                                                                      *
!************************************************************************

! ... warn if hard-core overlap between particle ip and other particles

logical function lWarnHCOverlap(ip,radaths,lflag)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip            ! the particle to be considered
   real(8),    intent(in) :: radaths(*)    ! hard-core radii of atom types
   logical,    intent(in) :: lflag         !=.true. consider only other particles with lpset ==.true.
                                           !=.false. consider all other particles except ip
   logical :: EllipsoidOverlap, SuperballOverlap

   integer(4) :: ipt, jp, jpt, ia, iat, ja, jat
   real(8)    :: dx, dy, dz, r2, dro(3), dropbc(3)

   lWarnHCOverlap =.true.
   ipt = iptpn(ip)
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      iat = iatan(ia)
      do jp = 1, np
         if(lflag) then
            if (.not.lpset(jp)) cycle  ! exclude particles not set
         else
            if (jp == ip) cycle        ! exclude jp == ip
         end if
         jpt = iptpn(jp)
         dro(1) = ro(1,ip)-ro(1,jp)
         dro(2) = ro(2,ip)-ro(2,jp)
         dro(3) = ro(3,ip)-ro(3,jp)
         call PBC2(dro(1),dro(2),dro(3),dropbc(1),dropbc(2),dropbc(3))
         do ja = ianpn(jp), ianpn(jp)+napt(jpt)-1
            jat = iatan(ja)
            dx = r(1,ia)-r(1,ja)-dropbc(1)
            dy = r(2,ia)-r(2,ja)-dropbc(2)
            dz = r(3,ia)-r(3,ja)-dropbc(3)
            r2 = dx**2+dy**2+dz**2
            if (r2 < (radaths(iat)+radaths(jat))**2) return
            if (lellipsoid) then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ipnan(ia)),ori(1,1,ipnan(ja)),radellipsoid2,aellipsoid)) return
            end if
            if (lsuperball) then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ipnan(ia)),ori(1,1,ipnan(ja)))) return
            end if
         end do
      end do
   end do
   lWarnHCOverlap =.false.

end function lWarnHCOverlap

!************************************************************************
!*                                                                      *
!*     WarnHCOverlap                                                    *
!*                                                                      *
!************************************************************************

! ... warn if hard-core overlap among atoms

subroutine WarnHCOverlap(iplow, ipupp)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow   ! lower limit of particles to be considered
   integer(4), intent(in) :: ipupp   ! upper limit of particles to be considered

   logical    :: loverlap, ltext, SuperballOverlap, ltemp
   integer(4) :: ip, jp, ia, ialow, iaupp, iat, ja, jat, iatjat
   real(8)    :: dxx, dyy, dzz, dx, dy, dz, dxpbc, dypbc, dzpbc, r2

   if (sum(r2atat(1:natat)) > Zero) then

     ltemp = lstatsuperball             ! to exlcude call of SuperballStat
     lstatsuperball = .false.

      ialow = ianpn(iplow)
      iaupp = ianpn(ipupp)+napt(iptpn(ipupp))-1

      ltext =.true.
      do ia = ialow, iaupp-1
         iat = iatan(ia)
         ip = ipnan(ia)
         do ja = ia+1, iaupp
            jat = iatan(ja)
            iatjat = iatat(iat,jat)
            jp = ipnan(ja)
            if (ip == jp) cycle
            dxx = ro(1,ip)-ro(1,jp)
            dyy = ro(2,ip)-ro(2,jp)
            dzz = ro(3,ip)-ro(3,jp)
            call PBC2(dxx,dyy,dzz,dxpbc,dypbc,dzpbc)
            dx = r(1,ia)-r(1,ja)-dxpbc
            if (dx > r1atat(iatjat)) cycle
            dy = r(2,ia)-r(2,ja)-dypbc
            if (dy > r1atat(iatjat)) cycle
            dz = r(3,ia)-r(3,ja)-dzpbc
            if (dz > r1atat(iatjat)) cycle
            r2 = dx**2+dy**2+dz**2
            loverlap = .false.
            if (lsuperball) then
               loverlap = SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ia),ori(1,1,ja))
            else
               if (r2 < r2atat(iatjat)) loverlap = .true.
            end if
            if (loverlap) then
               if (ltext) then
                  call Warn('WarnHCOverlap', 'hard-core overlap between atoms', uout)
                  write(uout,'(a,5x,a,5x,a,35x,a,25x,a)') 'particle no', 'atom no', 'atom type', 'coordinates', 'separation'
               end if
               ltext =.false.
               write(uout,'(2i5,t15,2i5,3x,a,1x,a,t50,3(3f8.3,5x))')      &
               ip, jp, ia, ja, txat(iat), txat(jat), r(1:3,ia), r(1:3,ja), sqrt(r2)
            end if
         end do
      end do

      lstatsuperball = ltemp             ! restore value of lusperballstat

   end if

end subroutine WarnHCOverlap

!************************************************************************
!*                                                                      *
!*     WarnAtomOutsideBox                                               *
!*                                                                      *
!************************************************************************

! ... warn if atoms are outside box

subroutine WarnAtomOutsideBox(ialow, iaupp)

   use MolModule
   implicit none

   integer(4), intent(in) :: ialow, iaupp

   logical    :: ltext
   integer(4) :: ia, iat
   real(8)    :: r2

   ltext =.true.
   if (lbcbox) then
      if (txbc == 'xyz') then
         continue
      else if (txbc == 'xy') then
         do ia = ialow, iaupp
            iat = iatan(ia)
            if (abs(r(3,ia)) > boxlen2(3)-radat(iat)) then
               if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the box', uout)
               ltext =.false.
               write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'z = ', r(3,ia)
            end if
         end do
      else if (txbc == 'z') then
         do ia = ialow, iaupp
            iat = iatan(ia)
            if (abs(r(1,ia)) > boxlen2(1)-radat(iat)) then
               if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the box', uout)
               ltext =.false.
               write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'x = ', r(1,ia)
            end if
            if (abs(r(2,ia)) > boxlen2(2)-radat(iat)) then
               if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the box', uout)
               ltext =.false.
               write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'x = ', r(2,ia)
            end if
         end do
      end if
   else if (lbcrd .or. lbcto) then
      continue
   else if (lbcsph) then
      do ia = ialow, iaupp
         iat = iatan(ia)
         r2 = r(1,ia)**2+r(2,ia)**2+r(3,ia)**2
         if (r2 > (sphrad+radat(iat))**2) then
            if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the sphere', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'r = ', sqrt(r2)
         end if
      end do
   else if (lbccyl) then
      do ia = ialow, iaupp
         iat = iatan(ia)
         r2 = r(1,ia)**2+r(2,ia)**2
         if (r2 > (cylrad+radat(iat))**2) then
            if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the cylinder', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'r = ', sqrt(r2)
         end if
         if (abs(r(3,ia)) > Half*cyllen-radat(ia)) then
            if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the cylinder', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'z = ', r(3,ia)
         end if
      end do
   else if (lbcell) then
      do ia = ialow, iaupp
         iat = iatan(ia)
         r2 = sum((r(1:3,ia)*ellradi(1:3))**2)                ! radat not yet involved
         if (r2 > One ) then
            if (ltext) call Warn('WarnAtomOutsideBox', 'atom outside the ellipsoid', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'atom', ia, 'r = ', sqrt(r2)
         end if
      end do
   end if

end subroutine WarnAtomOutsideBox

!************************************************************************
!*                                                                      *
!*     WarnPartOutsideBox                                               *
!*                                                                      *
!************************************************************************

! ... warn if particles are outside box

subroutine WarnPartOutsideBox(iplow, ipupp)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow, ipupp

   logical    :: ltext
   integer(4) :: ip
   real(8)    :: r2, dxpbc, dypbc, dzpbc

   ltext =.true.
   if (lbcbox) then
      do ip = iplow, ipupp
         if (abs(ro(1,ip)) > boxlen2(1)) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the box', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'x = ', ro(1,ip)
            write(uout,'(a,f10.5)') 'boxlen2', boxlen2(1)
         end if
         if (abs(ro(2,ip)) > boxlen2(2)) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the box', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'y = ', ro(2,ip)
         end if
         if (abs(ro(3,ip)) > boxlen2(3)) then
!           if (abs(ro(3,ip)) == (boxlen2(3)+2)) exit   ! alow mica charges 2 length units outside the box (FREDRIK C)
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the box', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'z = ', ro(3,ip)
         end if
      end do
   else if (lbcrd .or. lbcto) then
      do ip = iplow, ipupp
         call PBC2(ro(1,ip),ro(2,ip),ro(3,ip),dxpbc,dypbc,dzpbc)
         if (abs(dxpbc) > Zero .or. abs(dypbc) > Zero .or. abs(dzpbc) > Zero) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the cell', uout)
            ltext = .false.
            write(uout,'(a,i4,5x,a,3f10.5)') 'particle', ip, 'r(1:3) = ', ro(1:3,ip)
         end if
      end do
   else if (lbcsph) then
      do ip = iplow, ipupp
         r2 = ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2
         if (r2 > sphrad2) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the sphere', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'r = ', sqrt(r2)
         end if
      end do
   else if (lbccyl) then
      do ip = iplow, ipupp
         r2 = ro(1,ip)**2+ro(2,ip)**2
         if (r2 > cylrad2) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the cylinder', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'r = ', sqrt(r2)
         end if
         if (abs(ro(3,ip)) > Half*cyllen) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the cylinder', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'z = ', ro(3,ip)
         end if
      end do
   else if (lbcell) then
      do ip = iplow, ipupp
         r2 = sum((ro(1:3,ip)*ellradi(1:3))**2)
         if (r2 > One ) then
            if (ltext) call Warn('WarnPartOutsideBox', 'particle outside the ellipsoid', uout)
            ltext =.false.
            write(uout,'(a,i4,5x,a,f10.5)') 'particle', ip, 'r = ', sqrt(r2)
         end if
      end do
   end if

end subroutine WarnPartOutsideBox

!************************************************************************
!*                                                                      *
!*    UndoPBCChain                                                      *
!*                                                                      *
!************************************************************************

! ... undo the periodical boundary conditions for a chain molecule

subroutine UndoPBCChain(rref, ic, idir, rotemp)

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='UndoPBCChain'

   real(8),    intent(in)  :: rref(3)     ! reference point for the undo
   integer(4), intent(in)  :: ic          ! chain number
   integer(4), intent(in)  :: idir        ! direction: -1 or +1
   real(8),    intent(out) :: rotemp(3,*) ! coordinate of the undone chain

   integer(4) :: ip, jp, iseg, ipfirst, iplast
   real(8) :: dx, dy, dz

! ... determine first and last particle from idir

   if (idir == -1) then
      ipfirst = npct(ictcn(ic))
      iplast = 1
   else if (idir == 1) then
      ipfirst = 1
      iplast = npct(ictcn(ic))
   else
     call Stop(txroutine, 'invalid value of idir', uout)
   end if

! ... undo pbc of the first segment with respect to the reference point

   iseg = ipfirst
   ip = ipnsegcn(iseg,ic)
   dx = ro(1,ip) - rref(1)
   dy = ro(2,ip) - rref(2)
   dz = ro(3,ip) - rref(3)
   call PBC(dx,dy,dz)
   rotemp(1,ip) = rref(1) + dx
   rotemp(2,ip) = rref(2) + dy
   rotemp(3,ip) = rref(3) + dz


!   WRITE(*,*) 'UNDO, ip = ',ip
!   WRITE(*,'(a,3f10.3)') 'rref(1:3)= ',rref(1:3)
!   WRITE(*,'(a,3f10.3)') 'ro(1:3,ip)=',ro(1:3,ip)
!   WRITE(*,'(a,3f10.3)') 'dx, dy, dz=',dx, dy, dz

! ... undo pbc of following segments with respect to the previous one

   do iseg = ipfirst+idir, iplast, idir
      jp = ip
      ip = ipnsegcn(iseg,ic)
      dx = ro(1,ip) - ro(1,jp)
      dy = ro(2,ip) - ro(2,jp)
      dz = ro(3,ip) - ro(3,jp)
      call PBC(dx,dy,dz)
      rotemp(1,ip) = rotemp(1,jp) + dx
      rotemp(2,ip) = rotemp(2,jp) + dy
      rotemp(3,ip) = rotemp(3,jp) + dz
   end do

end subroutine UndoPBCChain

!************************************************************************
!*                                                                      *
!*    CalcChainProperty                                                 *
!*                                                                      *
!************************************************************************

! ... calculate properties of a single chain

subroutine CalcChainProperty(ic, rotemp, ChainProperty)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4),          intent(in)  :: ic            ! chain number
   real(8),             intent(in)  :: rotemp(3,*)   ! coordinate of the undone chain
   type(chainprop_var), intent(out) :: ChainProperty ! chain properties

   real(8) :: vsum, vsum2, vsumr, vsumz, vsumxy, dx, dy, dz, r2, xcom, ycom, zcom, norm, angle_deg, rbb
   real(8) :: dxm, dym, dzm, dxp, dyp, dzp, theta, mimat(3,3), diagonal(3), eivr(3,3)
   real(8) :: l2_small, l2_mid, l2_large, toroidparam
   !real(8) :: i_small, i_mid, i_large ! used in older code snippet below - currently not needed
   integer(4) :: iseg, ip, jp, jp_m, jp_p, ict, nrot
   real(8) :: hx, hy, hz, hxsum, hysum, hzsum
   real(8) :: InvFlt, PerLengthRg, Asphericity

   ict = ictcn(ic)

! .. center of mass

   xcom = sum(rotemp(1,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   ycom = sum(rotemp(2,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   zcom = sum(rotemp(3,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   call PBC(xcom,ycom,zcom)
   ChainProperty%ro(1) = xcom
   ChainProperty%ro(2) = ycom
   ChainProperty%ro(3) = zcom

! ... bead-to-bead distance squared

   vsum = Zero
   do iseg = 1, npct(ict)-1
      ip = ipnsegcn(iseg,ic)
      jp = ipnsegcn(iseg+1,ic)
      dx = rotemp(1,ip)-rotemp(1,jp)
      dy = rotemp(2,ip)-rotemp(2,jp)
      dz = rotemp(3,ip)-rotemp(3,jp)
      r2 = dx**2+dy**2+dz**2
      vsum = vsum + r2
   end do
   ChainProperty%rbb2 = vsum*InvInt(npct(ict)-1)

! ... angle between consecutive beads

   vsum = Zero
   vsum2 = Zero
   do iseg = 2, npct(ict)-1
      jp_m = ipnsegcn(iseg-1,ic)
      ip = ipnsegcn(iseg,ic)
      jp_p = ipnsegcn(iseg+1,ic)
      dxm = ro(1,ip)-ro(1,jp_m)
      dym = ro(2,ip)-ro(2,jp_m)
      dzm = ro(3,ip)-ro(3,jp_m)
      call PBC(dxm,dym,dzm)
      dxp = ro(1,ip)-ro(1,jp_p)
      dyp = ro(2,ip)-ro(2,jp_p)
      dzp = ro(3,ip)-ro(3,jp_p)
      call PBC(dxp,dyp,dzp)
      theta = angle_deg(dxm,dym,dzm,dxp,dyp,dzp)
      vsum = vsum+theta
      vsum2 = vsum2+cos(Pi*(One-theta/180.0d0))
   end do
   ChainProperty%angle = vsum*InvInt(npct(ict)-2)
   ChainProperty%cos = vsum2*InvInt(npct(ict)-2)

! ... end-to-end vector and end-to-end distance squared

   ip = ipnsegcn(1,ic)
   jp = ipnsegcn(npct(ict),ic)
   dx = rotemp(1,ip)-rotemp(1,jp)
   dy = rotemp(2,ip)-rotemp(2,jp)
   dz = rotemp(3,ip)-rotemp(3,jp)
   r2 = dx**2+dy**2+dz**2
   ChainProperty%ree(1) = dx
   ChainProperty%ree(2) = dy
   ChainProperty%ree(3) = dz
   ChainProperty%ree2 = r2

! ... radius of gyration squared and projections on the z-axis, the xy-plane and the prinicpal axes

!     rg**2 = l2_small**2 + l2_mid**2 + l2_large**2

   xcom = sum(rotemp(1,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   ycom = sum(rotemp(2,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   zcom = sum(rotemp(3,ipnsegcn(1:npct(ict),ic)))/npct(ict)
   vsumr = Zero
   vsumz = Zero
   vsumxy = Zero
   mimat(1:3,1:3) = Zero
   do iseg = 1, npct(ict)
      ip = ipnsegcn(iseg,ic)
      dx = rotemp(1,ip)-xcom
      dy = rotemp(2,ip)-ycom
      dz = rotemp(3,ip)-zcom
      r2 = dx**2+dy**2+dz**2
      vsumr = vsumr + r2
      vsumz = vsumz + dz**2
      vsumxy = vsumxy + dx**2+dy**2
      mimat(1,1) = mimat(1,1)+dx**2    !        +dy**2+dz**2
      mimat(1,2) = mimat(1,2)+dx*dy    !        -dx*dy
      mimat(1,3) = mimat(1,3)+dx*dz    !        -dx*dz
      mimat(2,2) = mimat(2,2)+dy**2    !        +dx**2+dz**2
      mimat(2,3) = mimat(2,3)+dy*dz    !        -dy*dz
      mimat(3,3) = mimat(3,3)+dz**2    !        +dx**2+dy**2
   end do
   ChainProperty%rg2 = vsumr/npct(ict)
   ChainProperty%rg2z = vsumz/npct(ict)
   ChainProperty%rg2xy = vsumxy/npct(ict)

   mimat(2,1) = mimat(1,2)
   mimat(3,1) = mimat(1,3)
   mimat(3,2) = mimat(2,3)
   call Diag(3, mimat, diagonal, eivr, nrot)
   l2_small = min(diagonal(1),diagonal(2),diagonal(3))/npct(ict)
   l2_large = max(diagonal(1),diagonal(2),diagonal(3))/npct(ict)
   l2_mid  = (diagonal(1)+diagonal(2)+diagonal(3))/npct(ict)-(l2_small+l2_large)
   ChainProperty%rg2s = max(Zero,l2_small)
   ChainProperty%rg2m = max(Zero,l2_mid  )
   ChainProperty%rg2l = max(Zero,l2_large)

! ... principal moments of inertia  ( i_i = mass*(rg**2-l2_i) )

!  i_small = npct(ict)*(ChainProperty%rg2-ChainProperty%rg2l)
!  i_mid   = npct(ict)*(ChainProperty%rg2-ChainProperty%rg2m)
!  i_large = npct(ict)*(ChainProperty%rg2-ChainProperty%rg2s)

! ... persistence lengths (ref JCP 107, 1279 (1997))

   rbb = sqrt(ChainProperty%rbb2)
   ChainProperty%lpree = ChainProperty%ree2 *InvFlt(Two*(npct(ict)-1)*rbb) + Half*rbb   ! based on end-to-end separation
   ChainProperty%lprg = PerLengthRg(ChainProperty%rg2, (npct(ict)-1)*rbb)          ! based on radius of gyration

! ... shape ratio  ( = 6 for gaussian chain and 12 for a rod)

   ChainProperty%shape = ChainProperty%ree2 * InvFlt(ChainProperty%rg2)

! ... asphericity  ( = 0 for sphere, 0.526 for gaussian chain, and 1 for a rod)

   ChainProperty%asph = Asphericity(l2_small,l2_mid,l2_large)

! ... toroidicity (sum of vector products between successive bond directions)

   hxsum = Zero
   hysum = Zero
   hzsum = Zero
   do iseg = 2, npct(ict)-1
      jp_m = ipnsegcn(iseg-1,ic)
      ip = ipnsegcn(iseg,ic)
      jp_p = ipnsegcn(iseg+1,ic)
      dxm = rotemp(1,ip)-rotemp(1,jp_m)
      dym = rotemp(2,ip)-rotemp(2,jp_m)
      dzm = rotemp(3,ip)-rotemp(3,jp_m)
      call PBC(dxm,dym,dzm)
      dxp = rotemp(1,ip)-rotemp(1,jp_p)
      dyp = rotemp(2,ip)-rotemp(2,jp_p)
      dzp = rotemp(3,ip)-rotemp(3,jp_p)
      call PBC(dxp,dyp,dzp)
      hx = dym*dzp - dzm*dyp
      hy = dzm*dxp - dxm*dzp
      hz = dxm*dyp - dym*dxp
      norm = InvFlt(sqrt(hx**2+hy**2+hz**2))
      hx = hx*norm
      hy = hy*norm
      hz = hz*norm
      hxsum = hxsum + hx
      hysum = hysum + hy
      hzsum = hzsum + hz
   end do
   toroidparam = sqrt(hxsum**2+hysum**2+hzsum**2)*InvInt((npct(ict)-2))
   ChainProperty%torp = toroidparam
end subroutine CalcChainProperty

real(8) function PerLengthRg(rg2, l)
   implicit none
   real(8), intent(in)  :: rg2         ! radius of gyration squared
   real(8), intent(in)  :: l           ! contour length
   integer(4) ::i
   real(8) :: PerLengthRgOld, fac

   PerLengthRg = 0.0d0
   if (l == 0.0d0) return
   do i = 1, 500
      PerLengthRgOld = PerLengthRg
      fac = PerLengthRg/l
      if(fac == 0.0d0) then ! prevent division by zero
         PerLengthRg = (3.0/l)*(rg2 + PerLengthRg**2.0)
      else
         PerLengthRg = (3.0/l)*(rg2 + PerLengthRg**2.0*(1-2*fac*(1.0-fac*(1.0-exp(-1.0/fac)))))
      end if
      if (abs(PerLengthRg-PerLengthRgOld)/PerLengthRg < 1e-3) return
   end do
end function PerLengthRg

!************************************************************************
!*                                                                      *
!*    CalcNetworkProperty                                               *
!*                                                                      *
!************************************************************************

! ... calculate properties of network number inw
! ... no prior undoing of the periodic boundary conditions necessary

subroutine CalcNetworkProperty(inw, NetworkProperty)

   use MolModule
   use MolauxModule, only: CalcCOM
   implicit none

   integer(4),            intent(in)  :: inw             ! network number
   type(networkprop_var), intent(out) :: NetworkProperty ! network properties

! ... properties
   real(8)     :: rcom(3)                    ! center of mass
   real(8)     :: rg2x, rg2y, rg2z           ! radius of gyration squared projected on the x, y and z-axes
   real(8)     :: l2_small, l2_mid, l2_large ! small, middle and large extension along principal axes
   real(8)     :: eivr(3,3)                  ! eigenvectors of the principal frame
   real(8)     :: Asphericity                ! asphericity
   real(8)     :: theta(3)                   ! angles of axes of largest extension and x-, y-, and z-axes of main frame
   real(8)     :: alpha                      ! degree of ionization

! ... processing variables
   real(8)     :: dr(3)
   real(8)     :: r2
   real(8)     :: vsumr
   real(8)     :: mimat(3,3)
   real(8)     :: diagonal(3)
   integer(4)  :: nrot
   integer(4)  :: npcharged

! ... counter
   integer(4)  :: inwt, ip, iploc
   integer(4)  :: irow

! ... determine network type

   inwt = inwtnwn(inw)

! ... center of mass

   rcom(1:3) = CalcCOM(ro(1:3,1:np),MASK=lpnnwn(1:np,inw),MASS=massp(1:np))
   NetworkProperty%ro(1:3) = rcom(1:3)

! ... radius of gyration squared and projections on the prinicpal axes

   ! rg**2 = l2_small**2 + l2_mid**2 + l2_large**2

   vsumr          = Zero
   mimat(1:3,1:3) = Zero
   do iploc = 1, npnwt(inwt)
      ip = ipnplocnwn(iploc,inw)
      dr(1:3) = ro(1:3,ip)-rcom(1:3)
      call PBCr2(dr(1),dr(2),dr(3),r2)
      vsumr = vsumr + r2
      do irow = 1, 3
         mimat(irow,1:3) = mimat(irow,1:3) + dr(irow)*dr(1:3)
      end do
   end do
   NetworkProperty%rg2 = vsumr*massinwt(inwt)

! ... radius of gyration squared projected on the x-, y- and z-axes

   rg2x = mimat(1,1)*massinwt(inwt)
   rg2y = mimat(2,2)*massinwt(inwt)
   rg2z = mimat(3,3)*massinwt(inwt)

   NetworkProperty%rg2x = rg2x
   NetworkProperty%rg2y = rg2y
   NetworkProperty%rg2z = rg2z

! ... normalized eigenvectors of the principal frame in descending order of the eigenvalues (due to Eigensort)

   call Diag(3,mimat,diagonal,eivr,nrot)
   call Eigensort(diagonal,eivr,3)

   NetworkProperty%eivr(1:3,1:3) = eivr(1:3,1:3)

! ... small, middle and large square extension along principal axes and eigenvectors of the principal frame

   l2_large = diagonal(1)*massinwt(inwt)
   l2_mid   = diagonal(2)*massinwt(inwt)
   l2_small = diagonal(3)*massinwt(inwt)

   NetworkProperty%rg2s = max(Zero,l2_small)
   NetworkProperty%rg2m = max(Zero,l2_mid  )
   NetworkProperty%rg2l = max(Zero,l2_large)

! ... angle of axes of largest extension (principal frame) and x,y,z-axes of main frame

   theta(1) = RadToDeg*acos(eivr(1,1)) ! eivr(1,1) is the dot product of the eigenvector of the largest extension and the normalized eigenvector of the x-axis
   theta(2) = RadToDeg*acos(eivr(2,1)) ! eivr(2,1) is the dot product of the eigenvector of the largest extension and the normalized eigenvector of the y-axis
   theta(3) = RadToDeg*acos(eivr(3,1)) ! eivr(3,1) is the dot product of the eigenvector of the largest extension and the normalized eigenvector of the z-axis
   NetworkProperty%theta = theta

! ... asphericity  ( = 0 for sphere, 0.526 for gaussian chain, and 1 for a rod)

   NetworkProperty%asph = Asphericity(l2_small,l2_mid,l2_large)

! ... degree of ionization

   if (lweakcharge) then
      npcharged = count(laz(1:np).and.lpnnwn(1:np,inw))
      alpha     = real(npcharged)/npweakchargenwt(inwt)
      NetworkProperty%alpha = alpha
   else
      NetworkProperty%alpha = Zero
   end if

end subroutine CalcNetworkProperty

!************************************************************************
!*                                                                      *
!*     Asphericity                                                      *
!*                                                                      *
!************************************************************************

! ... calculate the asphericity based on the 2:nd moments in the principal frame
! ... 0 for sphere, 0.526 for gaussian chain, and 1 for a rod

real(8) function Asphericity(l1, l2, l3)
   implicit none
   real(8), intent(in)  :: l1, l2, l3    ! 2:nd moments
   real(8) :: InvFlt
   Asphericity = ((l1-l2)**2 + (l2-l3)**2 + (l3-l1)**2) * InvFlt(2.0d0*(l1 + l2 + l3)**2)
end function Asphericity

!************************************************************************
!*                                                                      *
!*     CalcLTT                                                          *
!*                                                                      *
!************************************************************************

! ... calculate loop, tail, and tail characteristics for a single chain

subroutine CalcLTT(nseg, ltag, mnobj, nobj, lenobj, nsegobj)

   implicit none

   integer(4), parameter :: iloop = 1, itail = 2, itrain = 3
   integer(4), intent(in) :: nseg              ! number of segments in a chain
   logical,    intent(in) :: ltag(*)           ! label marking if the segment is adsorbed
   integer(4), intent(in) :: mnobj             ! max number of objects allowed
   integer(4), intent(out) :: nobj(3)          ! number of objects
   integer(4), intent(out) :: lenobj(3,mnobj)  ! length of objects
   integer(4), intent(out) :: nsegobj(3)       ! number of segments involved

   character(40), parameter :: txroutine ='CalcLTT'
   integer(4) :: n
   integer(4) :: iobj, iseg

! ... initiate

   nobj = 0
   lenobj = 0
   nsegobj = 0

! ... first segment

   iseg = 1
   if (ltag(iseg)) then
      iobj = itrain
   else
      iobj = itail
   end if
   n = 1

! ... remaining segments

   do iseg = 2, nseg
      if ((ltag(iseg) .and. ltag(iseg-1)) .or. (.not.ltag(iseg) .and. .not.ltag(iseg-1))) then  ! same object, just advance
         n = n + 1
      else                                           ! different object: store and initiate segment iseg
         nobj(iobj) = nobj(iobj)  + 1
         if (nobj(iobj) > mnobj) call Stop(txroutine, 'nobj(iobj) > mnobj', 6)
         lenobj(iobj,nobj(iobj)) = lenobj(iobj,nobj(iobj)) + n
         if (ltag(iseg)) then
            iobj = itrain
         else
            iobj = iloop
         end if
         n = 1
      end if

   end do

! ... store the end

   if (iobj == iloop) iobj = itail                ! if labelled loop, make last sequence a tail
   if (n == nseg .and. iobj == itail) return       ! chain not adsorbed

   nobj(iobj) = nobj(iobj) + 1
   if (nobj(iobj) > mnobj) call Stop(txroutine, 'nobj(iobj) > mnobj', 6)
   lenobj(iobj,nobj(iobj)) = lenobj(iobj,nobj(iobj)) + n

! ... calculate nsegobj

   nsegobj(iloop) = sum(lenobj(iloop,1:nobj(iloop)))
   nsegobj(itail) = sum(lenobj(itail,1:nobj(itail)))
   nsegobj(itrain) = sum(lenobj(itrain,1:nobj(itrain)))

! ... check

   if (sum(nsegobj(1:3)) /= nseg) call Stop(txroutine, 'sum(nsegobj(1:3)) /= nseg', 6)

!  write(*,*) nobj(1), ' ', nsegobj(1), ' ', lenobj(1,1:nobj(1))
!  write(*,*) nobj(2), ' ', nsegobj(2), ' ', lenobj(2,1:nobj(2))
!  write(*,*) nobj(3), ' ', nsegobj(3), ' ', lenobj(3,1:nobj(3))

end subroutine CalcLTT

!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

! ... routines for cluster analyses

!************************************************************************
!*                                                                      *
!*     CalcChainPairListGeneral                                         *
!*                                                                      *
!************************************************************************

! ... calculate a list of chain pairs based on particle-particle separation

subroutine CalcChainPairListGeneral(nobj, ichain, rcluster2, mnpair, npair, n1, n2)

   use MolModule
   implicit none

   integer(4), intent(in)  :: nobj           ! number of chais to be considered
   integer(4), intent(in)  :: ichain(*)      ! chains to be considered (global number)
   real(8),    intent(in)  :: rcluster2      ! maximal distance for consider two objects as a pair
   integer(4), intent(in)  :: mnpair         ! maximal number of pairs
   integer(4), intent(out) :: npair          ! number of pairs found
   integer(4), intent(out) :: n1(mnpair)     ! particle number of pairs
   integer(4), intent(out) :: n2(mnpair)     ! particle number of pairs

   character(40), parameter :: txroutine ='CalcChainPairListGeneral'
   real(8)    :: dx, dy, dz, r2
   integer(4) :: ict, ic, iobj, jct, jc, jobj, iseg, ip, jseg, jp

   npair = 0
   do iobj = 1, nobj
      ic = ichain(iobj)
      ict = ictcn(ic)
      do jobj = iobj+1, nobj
         jc = ichain(jobj)
         jct = ictcn(jc)
         do iseg = 1, npct(ict)
            ip = ipnsegcn(iseg,ic)
            do jseg = 1, npct(jct)
               jp = ipnsegcn(jseg,jc)
               dx = ro(1,ip)-ro(1,jp)
               dy = ro(2,ip)-ro(2,jp)
               dz = ro(3,ip)-ro(3,jp)
               call PBCr2(dx,dy,dz,r2)
               if (r2 < rcluster2) then                      ! particle-particle separation below rcluster
                  npair = npair+1
                  if (npair > mnpair) call Stop(txroutine, 'npair > mnpair', 6)
                  n1(npair) = iobj
                  n2(npair) = jobj
                  goto 10
               end if
            end do
         end do
10       continue
      end do
   end do

end subroutine CalcChainPairListGeneral

!************************************************************************
!*                                                                      *
!*     CalcChainComPairListGeneral                                      *
!*                                                                      *
!************************************************************************

! ... calculate a list of chain pairs based on com-com separation

subroutine CalcChainComPairListGeneral(nobj, ichain, rcluster2, mnpair, npair, n1, n2)

   use MolModule
   implicit none

   integer(4), intent(in)  :: nobj           ! number of chais to be considered
   integer(4), intent(in)  :: ichain(*)      ! chains to be considered (global number)
   real(8),    intent(in)  :: rcluster2      ! maximal distance for consider two objects as a pair
   integer(4), intent(in)  :: mnpair         ! maximal number of pairs
   integer(4), intent(out) :: npair          ! number of pairs found
   integer(4), intent(out) :: n1(mnpair)     ! particle number of pairs
   integer(4), intent(out) :: n2(mnpair)     ! particle number of pairs

   character(40), parameter :: txroutine ='CalcChainComPairListGeneral'
   real(8)    :: dx, dy, dz, r2, xcom, ycom, zcom, rcom(1:3,nc)
   integer(4) :: ict, ic, iobj, jc, jobj

   do ic = 1, nc                              ! get com for all chains
      ict = ictcn(ic)
      call UndoPBCChain(ro(1:3,ipnsegcn(1,ic)), ic, 1, vaux)
      xcom = sum(vaux(1,ipnsegcn(1:npct(ict),ic)))/npct(ict)
      ycom = sum(vaux(2,ipnsegcn(1:npct(ict),ic)))/npct(ict)
      zcom = sum(vaux(3,ipnsegcn(1:npct(ict),ic)))/npct(ict)
      call PBC(xcom,ycom,zcom)
      rcom(1,ic) = xcom
      rcom(2,ic) = ycom
      rcom(3,ic) = zcom
   end do

   npair = 0
   do iobj = 1, nobj
      ic = ichain(iobj)
      do jobj = iobj+1, nobj
         jc = ichain(jobj)
         dx = rcom(1,ic)-rcom(1,jc)
         dy = rcom(2,ic)-rcom(2,jc)
         dz = rcom(3,ic)-rcom(3,jc)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rcluster2) then
            npair = npair+1
            if (npair > mnpair) call Stop(txroutine, 'npair > mnpair', 6)
            n1(npair) = iobj
            n2(npair) = jobj
         end if
      end do
   end do

end subroutine CalcChainComPairListGeneral

!************************************************************************
!*                                                                      *
!*     CalcPartPairListAll                                              *
!*                                                                      *
!************************************************************************

! ... calculate a list of particle pairs based on their separations

subroutine CalcPartPairListAll(nobj, rcluster2, mnpair, npair, n1, n2)

   use MolModule
   implicit none

   integer(4), intent(in)  :: nobj           ! number of particles to be considered
   real(8),    intent(in)  :: rcluster2      ! maximal distance for consider two objects as a pair
   integer(4), intent(in)  :: mnpair         ! maximal number of pairs
   integer(4), intent(out) :: npair          ! number of pairs found
   integer(4), intent(out) :: n1(mnpair)     ! particle number of pairs
   integer(4), intent(out) :: n2(mnpair)     ! particle number of pairs

   character(40), parameter :: txroutine ='CalcPartPairListAll'
   real(8)    :: dx, dy, dz, r2
   integer(4) :: ip, iobj, jp, jobj

   npair = 0
   do iobj = 1, nobj
      ip = iobj
      do jobj = iobj+1, nobj
         jp = jobj
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rcluster2) then
            npair = npair+1
            if (npair > mnpair) call Stop(txroutine, 'npair > mnpair', 6)
            n1(npair) = iobj
            n2(npair) = jobj
         end if
      end do
   end do

end subroutine CalcPartPairListAll

!************************************************************************
!*                                                                      *
!*     CalcPartPairListGeneral                                          *
!*                                                                      *
!************************************************************************

! ... calculate a list of particle pairs

subroutine CalcPartPairListGeneral(nobj, iparticle, rcluster2, mnpair, npair, n1, n2)

   use MolModule
   implicit none

   integer(4), intent(in)  :: nobj           ! number of particles to be considered
   integer(4), intent(in)  :: iparticle(*)   ! particles to be considered (global number)
   real(8),    intent(in)  :: rcluster2      ! maximal distance for consider two objects as a pair
   integer(4), intent(in)  :: mnpair         ! maximal number of pairs
   integer(4), intent(out) :: npair          ! number of pairs found
   integer(4), intent(out) :: n1(mnpair)     ! particle number of pairs
   integer(4), intent(out) :: n2(mnpair)     ! particle number of pairs

   character(40), parameter :: txroutine ='CalcPartPairListGeneral'
   real(8)    :: dx, dy, dz, r2
   integer(4) :: ip, iobj, jp, jobj

   npair = 0
   do iobj = 1, nobj
      ip = iparticle(iobj)
      do jobj = iobj+1, nobj
         jp = iparticle(jobj)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rcluster2) then
            npair = npair+1
            if (npair > mnpair) call Stop(txroutine, 'npair > mnpair', 6)
            n1(npair) = iobj
            n2(npair) = jobj
         end if
      end do
   end do

end subroutine CalcPartPairListGeneral

!************************************************************************
!*                                                                      *
!*     CalcPartPairListIpt                                              *
!*                                                                      *
!************************************************************************

! ... calculate pair list for particles

subroutine CalcPartPairListIpt(iptcluster, rcluster2, mnpair, npair, n1, n2)

   use MolModule
   implicit none

   integer(4), intent(in)  :: iptcluster  ! particle type for pair search
   real(8),    intent(in)  :: rcluster2   ! max distance for consider two particles as a pair
   integer(4), intent(in)  :: mnpair      ! maximal number of pairs
   integer(4), intent(out) :: npair       ! number of pairs found
   integer(4), intent(out) :: n1(mnpair)  ! particle number of pairs
   integer(4), intent(out) :: n2(mnpair)  ! particle number of pairs

   character(40), parameter :: txroutine ='CalcPartPairListIpt'
   real(8)    :: dx, dy, dz, r2
   integer(4) :: ip, iplow, ipupp, jp

   npair = 0
   iplow = ipnpt(iptcluster)
   ipupp = ipnpt(iptcluster)+nppt(iptcluster)-1
   do ip = iplow, ipupp
      do jp = ip+1, ipupp
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rcluster2) then
            npair = npair+1
            if (npair > mnpair) call Stop(txroutine, 'npair > mnpair', 6)
            n1(npair) = ip-iplow+1
            n2(npair) = jp-iplow+1
         end if
      end do
   end do

end subroutine CalcPartPairListIpt

!************************************************************************
!*                                                                      *
!*     Perculation                                                      *
!*                                                                      *
!************************************************************************

! ... determine if a system of bonded objects is perculated

logical function Perculation(nobjtot, iclusteriobj, npair, n1, n2,  r, boxlen, ltest, unit)

   implicit none

   integer(4), intent(in)  :: nobjtot         ! number of objects to consider
   integer(4), intent(in)  :: iclusteriobj(*) ! pointer: object -> its cluster (cluster id = id of member with lowest id)
   integer(4), intent(in)  :: npair           ! number of object pairs
   integer(4), intent(in)  :: n1(*)           ! object number of a pair of objects
   integer(4), intent(in)  :: n2(*)           ! object number of a pair of objects
   real(8),    intent(in)  :: r(1:3,*)        ! coordinate of object
   real(8),    intent(in)  :: boxlen(3)       ! box lengths
   logical,    intent(in)  :: ltest           ! logical flag for test output
   integer(4), intent(in)  :: unit            ! output unit

   character(40), parameter :: txroutine ='Perculation'
   integer(4), parameter   :: mnobjtot = 100  ! maximum number of objects to consider
   integer(4), parameter   :: mpathlength = 100
   integer(4), parameter   :: mnobj = 100     ! maximum number of objects in a cluster
   integer(4), parameter   :: mnn = 10
   integer(4)              :: nsize(mnobjtot) ! sizes of the clusters
   integer(4)              :: nobj            ! number of objects in a cluseter
   integer(4)              :: obj(mnobj)      ! id of object
   integer(4)              :: nn(mnobj)       ! number of neighbours of object
   integer(4)              :: idnn(mnn,mnobj) ! id of neighbour
   integer(4)              :: objloc(mnobjtot)! local id of object

   integer(4) :: pathlength, pathobj(0:mpathlength), i, icluster, iobj, jobj, iobjloc, jobjloc, istep, ipair
   real(8) :: dr(1:3), boxlen2(3), dx, dy, dz
   logical :: loop
   integer(4) :: itemp(1:1)

   if (ltest) then
      write(unit,*)
      write(unit,*) 'In Perculation'
      write(unit,*) 'nobjtot = ',nobjtot
      write(unit,*) 'iclusteriobj(1:nobjtot) =', iclusteriobj(1:nobjtot)
      write(unit,*) 'npair =', npair
      write(unit,*) 'n1(1:npair) =', n1(1:npair)
      write(unit,*) 'n2(1:npair) =', n2(1:npair)
   end if

   if (nobjtot > mnobjtot) call stop(txroutine, 'nobjtot > mnobjtot', unit)

! ... get identity of the largest cluster

   nsize(1:nobjtot) = 0
   do iobj = 1, nobjtot
      nsize(iclusteriobj(iobj)) = nsize(iclusteriobj(iobj)) + 1
   end do
   itemp(1:1) = maxloc(nsize(1:nobjtot))
   icluster = itemp(1)
   if (ltest) write(unit,*) 'nsize(1:nobjtot) = ', nsize(1:nobjtot)
   if (ltest) write(unit,*) 'icluster =', icluster

! ... calculate nobj, obj, nn, and idnn for the cluster icluster

   iobjloc = 0
   do iobj = 1, nobjtot                                ! loop over all objects
      if (iclusteriobj(iobj) == icluster) then          ! check if object belongs to the cluster
         iobjloc = iobjloc + 1                         ! increase the number of objects in the cluster
         if (iobjloc > mnobj) call stop(txroutine, 'iobjloc > mnobj', unit)
         obj(iobjloc) = iobj                           ! assign pointer iobjloc -> iobj
         objloc(iobj) = iobjloc                        ! assign pointer iobj -> iobjloc
         nn(iobjloc) = 0
         idnn(1:mnn, iobjloc) = 0
         do ipair = 1, npair                           ! loop over all pairs to set nn and idnn
            if (n1(ipair) == iobj) then
               nn(iobjloc) = nn(iobjloc) + 1
               if (nn(iobjloc) > mnn) call stop(txroutine, 'nn(iobjloc) > mnn', unit)
               idnn(nn(iobjloc),iobjloc) = n2(ipair)
            else if (n2(ipair) == iobj) then
               nn(iobjloc) = nn(iobjloc) + 1
               if (nn(iobjloc) > mnn) call stop(txroutine, 'nn(iobjloc) > mnn', unit)
               idnn(nn(iobjloc),iobjloc) = n1(ipair)
            end if
         end do
      end if
   end do
   nobj = iobjloc                                      ! save the number of objects in the cluster

   if (ltest) then
      write(unit,*)
      do iobjloc = 1, nobj
         write(unit,*) 'iobjloc, iobj, nn(), idnn() =', iobjloc, obj(iobjloc), nn(iobjloc), idnn(1:nn(iobjloc),iobjloc)
      end do
   end if

! ... initiate

   pathlength = 0
   iobjloc = 1
   iobj = obj(iobjloc)
   pathobj(pathlength) = iobj
   dr(1:3) = 0.0
   boxlen2(1:3) = 0.5d0*boxlen(1:3)
   loop = .false.
   perculation = .false.

! ... loop over events

   do istep = 1, 1000

      if (ltest) then
         write(unit,*)
         write(unit,*) 'Top of event loop'
         write(unit,*) 'pathlength, pathobj(pathlength), iobj =', pathlength, pathobj(pathlength), iobj
         write(unit,*) 'iobjloc, iobj, nn(), idnn() =', iobjloc, obj(iobjloc), nn(iobjloc), idnn(1:nn(iobjloc),iobjloc)
      end if

! ... check if reverse

      if (loop .or. nn(iobjloc) == 0) then
         if (ltest) write(unit,*) 'reverse'
         do i = pathlength-1, 0, -1
            iobj = pathobj(i)
            iobjloc = objloc(iobj)
            dx = r(1,iobj) - r(1,pathobj(i+1))
            dy = r(2,iobj) - r(2,pathobj(i+1))
            dz = r(3,iobj) - r(3,pathobj(i+1))
            if (abs(dx) > boxlen2(1)) dx = dx - boxlen(1)*sign(1.0d0,dx)
            if (abs(dy) > boxlen2(2)) dy = dy - boxlen(2)*sign(1.0d0,dy)
            if (abs(dz) > boxlen2(3)) dz = dz - boxlen(3)*sign(1.0d0,dz)
            dr(1) = dr(1) + dx
            dr(2) = dr(2) + dy
            dr(3) = dr(3) + dz
            if (ltest) then
               write(unit,*) 'test pathlength, i, iobjloc, iobj, nn() =', i, iobjloc, iobj, nn(iobjloc)
               write(unit,'(a,3f10.2)') 'dr',dr
            end if
            if (nn(iobjloc) > 0) exit      ! object has nonvisited neighbours
         end do
         if (i == -1) exit     ! the full graph has been examined
         pathlength = i
         if (ltest) write(unit,*) 'new pathlength and new obj', pathlength, iobj
         loop = .false.
      end if

! ... advance in the graph

      jobj = idnn(1, iobjloc)                           ! identify new object
      jobjloc = objloc(jobj)                            ! identify new object local
      if (ltest) then
         write(unit,*) 'before bonds removed'
         write(unit,*) 'iobjloc, iobj, nn(), idnn() =', iobjloc, obj(iobjloc), nn(iobjloc), idnn(1:nn(iobjloc),iobjloc)
         write(unit,*) 'jobjloc, jobj, nn(), idnn() =', jobjloc, obj(jobjloc), nn(jobjloc), idnn(1:nn(jobjloc),jobjloc)
      end if
      call updatenn(jobj, nn(iobjloc), idnn(1,iobjloc)) ! remove forward bond
      call updatenn(iobj, nn(jobjloc), idnn(1,jobjloc)) ! remove reverse bond
      dx = r(1,jobj) - r(1,iobj)
      dy = r(2,jobj) - r(2,iobj)
      dz = r(3,jobj) - r(3,iobj)
      call PBC(dx,dy,dz)
      dr(1) = dr(1) + dx
      dr(2) = dr(2) + dy
      dr(3) = dr(3) + dz

      if (ltest) then
         write(unit,*) 'after bonds removed'
         write(unit,*) 'iobjloc, iobj, nn(), idnn() =', iobjloc, obj(iobjloc), nn(iobjloc), idnn(1:nn(iobjloc),iobjloc)
         write(unit,*) 'jobjloc, jobj, nn(), idnn() =', jobjloc, obj(jobjloc), nn(jobjloc), idnn(1:nn(jobjloc),jobjloc)
         write(unit,'(a,3f10.2)') 'dr =',dr
      end if

      pathlength = pathlength + 1          ! update pathlength
      if (pathlength > mpathlength) call stop(txroutine, 'pathlength > mpathlength',6)
      iobj = jobj                          ! update current object
      iobjloc = jobjloc                    ! update current object local
      pathobj(pathlength) = iobj           ! update pathobj with current object

      if (ltest) write(unit,*) 'pathlength, pathobj(pathlength), iobj =', pathlength, pathobj(pathlength), iobj

! ... chech for loop

      do i = pathlength-1, 0, -1
         if (iobj == pathobj(i)) loop = .true.
      end do
      if (ltest) write(unit,*) 'loop', loop

! ... check for perculation

      if ((iobj == pathobj(0)) .and. (count(abs(dr(1:3)) > 0.25*minval(boxlen(1:3))) > 0)) then
         perculation = .true.
         exit
      end if
   end do
   if (istep > 1000) call warn('perculation', 'istep > 1000', unit)

   if (ltest) write(unit,*) 'perculation', perculation

contains

!........................................................................

subroutine updatenn(iobj, nn, idnn)  ! update nn and idnn
   implicit none
   integer(4), intent(in) :: iobj
   integer(4), intent(inout) :: nn
   integer(4), intent(inout) :: idnn(*)
   integer(4) :: i
   logical :: hit
   hit = .false.
   do i = 1, nn
      if (idnn(i) == iobj) then
         idnn(i:nn-1) = idnn(i+1:nn)
         nn = nn - 1
         hit = .true.
         exit
      end if
   end do
   if (.not.hit) then
      write(unit,*)
      write(unit,*) 'no hit in updatenn'
      write(unit,*) 'iobj, nn, idnn(1:nn)', iobj, nn, idnn(1:nn)
      stop 1
   end if
end subroutine updatenn

!........................................................................

end function Perculation

!¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤

!************************************************************************
!*                                                                      *
!*     CorrAnalysis                                                     *
!*                                                                      *
!************************************************************************

! ... make correlation analysis of a single variable

! use: insert the following alls in the routine calculating the property of interest


! at case (iBeforeSimulation):
!       call CorrAnalysis('initialize',xxx,blmin)
!
! at case (iSimulationStep):
!       call CorrAnalysis('sample',xxx,blmin)
!
! at case (iAfterSimulation):
!       call CorrAnalysis('finalize',xxx,blmin)
!
! where xxx is the value of the variable to be examined after a time step/pass (real(8))
! and blmin is smallest block length involved (integer(4))

subroutine CorrAnalysis(whattodo, variable, blmin)

   implicit none

   character(10), intent(in) :: whattodo    ! controls task to be performed
   real(8),       intent(in) :: variable    ! value of variable to be averaged
   integer(4),    intent(in) :: blmin       ! smallest block length to consider

   integer(4), save :: ndata                ! number of data points
   integer(4), save :: isum                 ! counter for primary average
   real(8),    save :: aver                 ! primary average


   if (whattodo(1:10) == 'initialize') then

      open(21, status = 'scratch', form = 'formatted')
      ndata = 0
      isum = 0
      aver = 0.0d0

   else if (whattodo(1:6) == 'sample') then

      ndata = ndata + 1                    ! increment ndata
      isum = isum + 1                      ! increment isum
      aver = aver + variable               ! add variable to aver
      if (isum == blmin) then              ! form average containing blmin values
         aver = aver/isum                  ! make the average
         write(21,*) aver                  ! write the average on unit 21
         isum = 0
         aver = 0.0d0
      end if

   else if (whattodo(1:8) == 'finalize') then

      open(22, file = 'blockaveranalysis_result', status = 'unknown')
      call BlockAverAnalysis(ndata,21,22)  ! averaging using different block lengths
      close(21)
      close(22)

   end if

end subroutine CorrAnalysis

!************************************************************************
!*                                                                      *
!*     BlockAverAnalysis                                                *
!*                                                                      *
!************************************************************************

! ... read data and calculate average and precision for different block lengths

subroutine BlockAverAnalysis(ndata, unitin, unitout)

   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: ndata                                  ! number of data points
   integer(4), intent(in) :: unitin                                 ! unit containing the data
   integer(4), intent(in) :: unitout                                ! unit to write output data

   integer(4), parameter :: mnblocklen = 20
   real(8), parameter :: Zero = 0.0d0
   integer(4) :: nblocklen
   integer(4) :: blocklen(mnblocklen)
   integer(4) :: nblock(mnblocklen)
   real(8)    :: av_s1(mnblocklen)
   real(8)    :: av_sd(mnblocklen)
   real(8)    :: av_s2(mnblocklen)
   integer(4) :: ibl, idata
   real(8)    :: data, norm, norm1
   real(8)    :: xfit(50), yfit(50), wfit(50), afit(0:2), PolVal, dum1, dum2, av_sd_extrap, av_stateff

! ... initiation

   nblock = 0
   av_s1  = Zero
   av_sd  = Zero
   av_s2  = Zero

   blocklen = [ (2**(ibl-1), ibl = 1, mnblocklen) ]

! ... calculation

   nblocklen = 20

   rewind (unitin)
   do idata = 1, ndata                                              ! loop over data points
      read(unitin,*) data

      do ibl = 1, nblocklen                                         ! loop over block lengths
         av_s2(ibl) = av_s2(ibl) + data                             ! add data to av_s2
         if (mod(idata,blocklen(ibl)) == 0) then                    ! check for end of block length
            nblock(ibl) = nblock(ibl) + 1                           ! update number of blocks
            av_s2(ibl) = av_s2(ibl) / blocklen(ibl)                 ! get average of the block
            av_s1(ibl) = av_s1(ibl) + av_s2(ibl)                    ! sum up av_sd for run average
            av_sd(ibl) = av_sd(ibl) + av_s2(ibl)**2                 ! sum up av_s2 precision of run average
            av_s2(ibl) = Zero                                       ! initialize av_s2 for next block
         end if
      end do

   end do

   do ibl = 1, nblocklen                                            ! loop over block lengths
      norm = InvInt(nblock(ibl))
      norm1 = InvInt(nblock(ibl)-1)
      av_s1(ibl) = av_s1(ibl)*norm                                  ! predicted run average for block length blocklen(i)
      av_sd(ibl) = sqrt(max(Zero,av_sd(ibl)*norm-av_s1(ibl)**2)*norm1) ! predicted precision of run average for block length blocklen(i)
      xfit(ibl) = log10(real(nblock(ibl)))                          ! save for linear fit
      yfit(ibl) = av_sd(ibl)                                        ! save for linear fit
      wfit(ibl) = sqrt(real(nblock(ibl)))                           ! save for linear fit
      if (nblock(ibl) == 1) exit                                    ! exit when number of blocks becomes 1
   end do
   nblocklen = ibl                                                  ! number of block lengths used

! ... extrapolation by linear fit of av_sd vs log(nblock) to nblock = 1

   call PolFit(1, nblocklen, xfit, yfit, wfit, 0, 6, afit, dum1, dum2)
   av_sd_extrap = max(Zero,PolVal(1,afit,0.0d0))
   av_stateff = (nblock(1)/(nblock(1)-1))*(av_sd_extrap/av_sd(1))**2

! ... output

   write(unitout,*) 'results from block average routine'
   write(unitout,*) '----------------------------------'
   write(unitout,*)
   write(unitout,'(a)') '   nblock       blocklen         product         average          sd             fit '
   write(unitout,'(a)') '   ------       --------         -------         -------          --             --- '
   write(unitout,'(3(i10,a),g14.7,a,g12.5,a,g12.5)') &
     (nblock(ibl), char(9), blocklen(ibl), char(9), nblock(ibl)*blocklen(ibl), char(9), &
      av_s1(ibl), char(9), av_sd(ibl), char(9), PolVal(1,afit,xfit(ibl)), ibl = 1,nblocklen)
   write(unitout,*)
   write(unitout,'(a,t35,g12.5)') 'average value            = ', av_s1(1)
   write(unitout,'(a,t35,g12.5)') 'fluctuation              = ', av_sd(1)*sqrt(real(nblock(1)-1))
   write(unitout,'(a,t35,g12.5)') 'extrapolated uncertainty = ', av_sd_extrap
   write(unitout,'(a,t35,g12.5)') 'statistical inefficiency = ', av_stateff

end subroutine BlockAverAnalysis

!************************************************************************
!*                                                                      *
!*     WriteTrace                                                       *
!*                                                                      *
!************************************************************************

! ... write trace information

subroutine WriteTrace(ilevel,name , iStage)

   use MolModule
   implicit none

   integer(4),   intent(in) :: ilevel
   character(*), intent(in) :: name
   integer(4),   intent(in) :: iStage

   logical, save :: first = .true.
   integer(4) :: unit

   if (first) then
      open(40, file = 'trace.master.data', status = 'unknown')
      if (nproc > 1) open(41, file = 'trace.slave.data', status = 'unknown')
      first = .false.
   end if

   if (master) unit = 40
   if (slave) unit = 41
   if (ilevel == 1) write(unit,'(i3,a,i2,t15,a,i5)') myid, ' iStage =', iStage, name
   if (ilevel == 2) write(unit,'(i3,a,i2,t25,a,i5)') myid, ' iStage =', iStage, name
   if (ilevel == 3) write(unit,'(i3,a,i2,t35,a,i5)') myid, ' iStage =', iStage, name, ipass
   if (ilevel == 4) write(unit,'(i3,a,i2,t45,a,i5)') myid, ' iStage =', iStage, name
   call fileflush(unit)

end subroutine WriteTrace

!************************************************************************
!*                                                                      *
!*     WriteVecAppend                                                   *
!*                                                                      *
!************************************************************************

! ... write real array at the end of an external unit

subroutine WriteVecAppend(ndata, data, funit)

   implicit none

   integer(4),   intent(in) :: ndata
   real(8),      intent(in) :: data(ndata)
   character(*), intent(in) :: funit

   character(40), parameter :: txroutine ='WriteVecApped'
   integer(4), save :: unit = 90
   integer(4) :: i

   if (ndata > 12) call Stop(txroutine, 'ndata too large', 6)
   open(unit, file = trim(adjustl(funit)), position = 'append', status = 'unknown')
   write(unit,'(12(g13.5,a))') (data(i), char(9), i = 1,ndata)
   close(unit)

end subroutine WriteVecAppend

module XModule
   implicit none
   real(8) :: rad2(3), rad2i, r12(1:3), ori1(1:3), ori2(1:3), e2, a(1:3,1:3,1:2)
end module XModule

!************************************************************************
!*                                                                      *
!*     EllipsoidOverlap                                                 *
!*                                                                      *
!************************************************************************

! ... check overlap between two ellipsoids with equal semi-axes
!     Ref. Perram & Wertheim J. Comp. Phys. 58, 409 (1984)

logical function EllipsoidOverlap(r2, r12_, ori1_, ori2_, rad2_, e_)

   use XModule
   implicit none

   real(8), intent(in) :: r2        ! distance squared between ellipsoid 1 and 2
   real(8), intent(in) :: r12_(3)   ! vector between ellipsoid 1 and 2
   real(8), intent(in) :: ori1_(3,3)! orientation matrix of ellipsoid 1
   real(8), intent(in) :: ori2_(3,3)! orientation matrix of ellipsoid 2
   real(8), intent(in) :: rad2_     ! square of degenerated semi-axes
   real(8), intent(in) :: e_        ! aspect ratio (>1 prolate, <1 oblate)

   real(8),    parameter :: zero = 0.0d0, half = 0.5d0, one = 1.0d0, four = 4.0d0
   character(40), parameter :: txroutine ='EllipsoidOverlap'
   real(8),    parameter :: functol = 1.0d-2              ! function tolerance
   real(8),    parameter :: lamthr = one                  ! threshold for abortion of maximum searching
   real(8) :: lam1, lam2, lam3, lammax, funcmax, BrentMod
   real(8), external :: Ellipsoidfunc, Ellipsoidfunc2
   real(8) :: func
   integer(4) :: imethod = 1
   logical :: ltest = .false.

     if (e_ > one) then
        if (r2 > four*rad2_*e_**2) then
           EllipsoidOverlap = .false.
           return
        end if
        if (r2 < four*rad2_) then
           EllipsoidOverlap = .true.
           return
        end if
     else
        if (r2 > four*rad2_) then
           EllipsoidOverlap = .false.
           return
        end if
        if (r2 < four*rad2_/e_**2) then
           EllipsoidOverlap = .true.
           return
        end if
     end if

   if(ltest) then
      call WriteHead(3, 'Test'//trim(txroutine), 6)
      write(*,'(a,4f10.5)') 'r12',r12_, sqrt(sum(r12_(1:3)**2))
      write(*,'(a,9f10.5)') 'ori1z''',ori1_
      write(*,'(a,9f10.5)') 'ori2z''',ori2_
      write(*,'(a,3f10.5)') 'rad2',rad2_
      write(*,'(a,3f10.5)') 'e',e_
   end if

   if (imethod == 0 .or. imethod == 1) then  ! PW

! ... make connection to Ellipsoidfunc

   r12 = r12_
   ori1 = ori1_(1:3,3)
   ori2 = ori2_(1:3,3)
   rad2i = one/rad2_
   e2 = e_**2

   func = Ellipsoidfunc(0.5d0)

! ... ininitate lambda

   lam1 = zero
   lam2 = half
   lam3 = one

! ... determine maximal value of Ellipsoidfunc, abort if larger than lamthr

   funcmax = -BrentMod(lam1, lam2, lam3, Ellipsoidfunc, functol, -lamthr, lammax)

   if (ltest)  write(*,'(2a,f25.18)') trim(txroutine), ': overlap function at 0.50 ',func
   if (ltest) write(*,'(2a,f25.18)') trim(txroutine), ': contact function ',funcmax

   endif

   if (imethod == 0 .or. imethod == 2) then ! MatInv

! ... make connection to Ellipsoidfunc2

   r12 = r12_
   rad2 = rad2_
   rad2(3) = rad2(1)*e_**2

   call MatDiagMatT(ori1_(1,1),rad2(1),a(1,1,1))
   call MatDiagMatT(ori2_(1,1),rad2(1),a(1,1,2))
   func = Ellipsoidfunc2(0.5d0)

! ... ininitate lambda

   lam1 = zero
   lam2 = half
   lam3 = one

! ... determine maximal value of Ellipsoidfunc, abort if larger than lamthr

   funcmax = -BrentMod(lam1, lam2, lam3, Ellipsoidfunc2, functol, -lamthr, lammax)

   if (ltest)  write(*,'(2a,f25.18)') trim(txroutine), ': overlap function at 0.50 ',func
   if (ltest) write(*,'(2a,f25.18)') trim(txroutine), ': contact function ',funcmax

   end if

   if (ltest) stop 1

   if (funcmax < one) then
      EllipsoidOverlap = .true.
   else
      EllipsoidOverlap = .false.
   end if
end function EllipsoidOverlap

!........................................................................

real(8) function Ellipsoidfunc(lam)

   use XModule
   implicit none
   real(8), intent(in) :: lam

   real(8),    parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0
   integer :: i, j
   real(8) :: oriprod, x, y, di, alpha, beta, gamma

   oriprod = ori1(1)*ori2(1) + ori1(2)*ori2(2) + ori1(3)*ori2(3)
   x = (one - lam)*(one - e2)
   y = lam*(one - e2)
   di = one/((one - x)*(one - y) - x*y*oriprod**2)
   alpha = x*(one - y)*di
   beta = y*(one - x)*di
   gamma = x*y*oriprod*di

   ellipsoidfunc = r12(1)**2 + r12(2)**2 + r12(3)**2
    do i = 1, 3
       do j = 1, 3
          ellipsoidfunc = ellipsoidfunc &
          + r12(i)*r12(j)*(alpha*ori1(i)*ori1(j) + beta*ori2(i)*ori2(j) + gamma*(ori2(i)*ori1(j) + ori1(i)*ori2(j)))
       end do
    end do
!    ellipsoidfunc = ellipsoidfunc &                       ! with O3: not faster
!      + r12(1)*r12(1)*(alpha*ori1(1)*ori1(1) + beta*ori2(1)*ori2(1) + two*gamma*ori2(1)*ori1(1)) &
!      + r12(2)*r12(2)*(alpha*ori1(2)*ori1(2) + beta*ori2(2)*ori2(2) + two*gamma*ori1(2)*ori2(2)) &
!      + r12(3)*r12(3)*(alpha*ori1(3)*ori1(3) + beta*ori2(3)*ori2(3) + two*gamma*ori1(3)*ori2(3)) &
!      + two*(r12(1)*r12(2)*(alpha*ori1(1)*ori1(2) + beta*ori2(1)*ori2(2) + gamma*(ori2(1)*ori1(2) + ori1(1)*ori2(2))) &
!      + r12(1)*r12(3)*(alpha*ori1(1)*ori1(3) + beta*ori2(1)*ori2(3) + gamma*(ori2(1)*ori1(3) + ori1(1)*ori2(3))) &
!      + r12(2)*r12(3)*(alpha*ori1(2)*ori1(3) + beta*ori2(2)*ori2(3) + gamma*(ori2(2)*ori1(3) + ori1(2)*ori2(3))))
   ellipsoidfunc  = -lam*(one-lam)*ellipsoidfunc*rad2i

end function Ellipsoidfunc

real(8) function Ellipsoidfunc2(lam)

   use XModule
   implicit none
   real(8), intent(in) :: lam

   real(8),    parameter :: one = 1.0d0, two = 2.0d0
   real(8) :: mat(3,3), matinv(3,3)

   mat(1:3,1:3) = (one-lam)*a(1:3,1:3,1) + lam*a(1:3,1:3,2)
   call MatInv3(mat,matinv)
   ellipsoidfunc2 = -lam*(one-lam)*(r12(1)*r12(1)*matinv(1,1) + r12(2)*r12(2)*matinv(2,2) + r12(3)*r12(3)*matinv(3,3) &
                         + Two*(r12(1)*r12(2)*matinv(1,2) + r12(1)*r12(3)*matinv(1,3) + r12(2)*r12(3)*matinv(2,3)))

end function Ellipsoidfunc2

subroutine MatDiagMatT(ori,diag,mat)
   implicit none
   real(8),    intent(in)  :: ori(3,3), diag(3)
   real(8),    intent(out) :: mat(3,3)

   mat(1,1) = ori(1,1)*diag(1)*ori(1,1) + ori(1,2)*diag(2)*ori(1,2) + ori(1,3)*diag(3)*ori(1,3)
   mat(1,2) = ori(1,1)*diag(1)*ori(2,1) + ori(1,2)*diag(2)*ori(2,2) + ori(1,3)*diag(3)*ori(2,3)
   mat(1,3) = ori(1,1)*diag(1)*ori(3,1) + ori(1,2)*diag(2)*ori(3,2) + ori(1,3)*diag(3)*ori(3,3)
   mat(2,1) = mat(1,2)
   mat(2,2) = ori(2,1)*diag(1)*ori(2,1) + ori(2,2)*diag(2)*ori(2,2) + ori(2,3)*diag(3)*ori(2,3)
   mat(2,3) = ori(2,1)*diag(1)*ori(3,1) + ori(2,2)*diag(2)*ori(3,2) + ori(2,3)*diag(3)*ori(3,3)
   mat(3,1) = mat(1,3)
   mat(3,2) = mat(2,3)
   mat(3,3) = ori(3,1)*diag(1)*ori(3,1) + ori(3,2)*diag(2)*ori(3,2) + ori(3,3)*diag(3)*ori(3,3)
end subroutine MatDiagMatT

!************************************************************************
!*                                                                      *
!*     SuperballOverlap                                                 *
!*                                                                      *
!************************************************************************

! ... superball overlap check

logical function SuperballOverlap(r2, r21, ori1, ori2) result(loverlap)
   use MolModule
   implicit none

   real(8), intent(in) :: r2            ! distance squared between superball 1 and 2
   real(8), intent(in) :: r21(3)        ! vector from superball 2 to superball 1
   real(8), intent(in) :: ori1(3,3)     ! orientation matrix of superball 1
   real(8), intent(in) :: ori2(3,3)     ! orientation matrix of superball 2

   real(8) :: ori2_1(3,3), r12_1(3)     ! coordinates in ball 1's frame
   real(8) :: t0, t1, of
   logical :: loverlap2, SuperballOverlap_NR

   if (lstatsuperball) then
      if (txmethodsuperball == 'mesh') then
         call cpu_time(t0)
         dopc = 0
         tric = 0
      end if
   end if

!  write(*,*) 'SuperballOverlap: txmethodsuperball', txmethodsuperball
!  write(*,'(a,3f20.10)') 'rcut2superball,r2',rcut2superball, r2

   if (abs(r21(1)**2 + r21(2)**2 + r21(3)**2 - r2) > 1e-10) call stop('SuperballOverlap','error',6)   ! consistency check

   if (r2 > rcut2superball(2)) then         ! separation > curcumscribing sphere
      loverlap = .false.
   else if(r2 < rcut2superball(1)) then     ! separation < inscribing sphere
      loverlap = .true.
   else                                     ! more detailed overlap check is needed
      if (txmethodsuperball == 'nr') then                 ! NR
         loverlap = SuperballOverlap_NR(r21, ori1, ori2, of)
!         write(*,*) 'HÄR 2', loverlap
      else if (txmethodsuperball == 'mesh') then          ! mesh
         r12_1 = matmul(transpose(ori1),-r21)
         ori2_1 = matmul(transpose(ori1),ori2)
         loverlap = OverlapMesh(superBallMesh,r12_1,ori2_1)
      else if (txmethodsuperball == 'test') then          ! both NR and mesh (test purpose)
         loverlap = SuperballOverlap_NR(r21, ori1, ori2, of)
         r12_1 = matmul(transpose(ori1),-r21)
         ori2_1 = matmul(transpose(ori1),ori2)
         loverlap2 = OverlapMesh(superBallMesh,r12_1,ori2_1)
         call WriteHead(3, 'SuperballOverlap', uout)
         write(uout,'(a,2l3,a)') 'overlap (NR, Mesh): (',loverlap,loverlap2,')'
         stop 1
      end if
   end if

   if (lstatsuperball) then
      if (txmethodsuperball == 'mesh') then
         call cpu_time(t1)
         call SuperballStatMesh(iSimulationStep, r21, loverlap, t1-t0)
      end if
   end if

end function SuperballOverlap

!************************************************************************
!*                                                                      *
!*     SuperballOverlapOF                                               *
!*                                                                      *
!************************************************************************

! ... superball overlap check, return overlap function

!real(8) function SuperballOverlapOF(r2, r21, ori1, ori2) result(of) !r2 is not used
real(8) function SuperballOverlapOF(r21, ori1, ori2) result(of)
   use MolModule
   implicit none
   !real(8), intent(in) :: r2            ! distance squared between superball 1 and 2
   real(8), intent(in) :: r21(3)        ! vector from superball 2 to superball 1
   real(8), intent(in) :: ori1(3,3)     ! orientation matrix of superball 1
   real(8), intent(in) :: ori2(3,3)     ! orientation matrix of superball 2
   logical :: loverlap, SuperballOverlap_NR
   loverlap = SuperballOverlap_NR(r21, ori1, ori2, of)
end function SuperballOverlapOF

!************************************************************************
!*                                                                      *
!*     SuperballOverlap_NR                                              *
!*                                                                      *
!************************************************************************

! ... check for overlap between two superballs
!     ref: Ni et al., Soft Matter 2012, 8, 8826

logical function SuperballOverlap_NR(r21, ori1, ori2, of) result(loverlap)
   use MolModule
   implicit none

   real(8), intent(in) :: r21(3)        ! vector between superball 1 and 2
   real(8), intent(in) :: ori1(3,3)     ! orientation matrix of superball 1
   real(8), intent(in) :: ori2(3,3)     ! orientation matrix of superball 2
   real(8), intent(out) :: of           ! overlap function
   integer(4) :: k, iter
   real(8) :: lambda, lambda0, ori1i(3,3), ori2i(3,3), orihelp(3,3,2), rc(3), rc0(3), rsb(3,2), dr(3), dll, dg(3), dl, vtemp(3)
   real(8) :: sf(2), sfd(3,2), sfdd(3,3,2), ofd(3), ofdd(3,3), minv(3,3)
   real(8) :: r1help, r2help, r3help, gp, gpp, fac1, fac1help, gradsftilde(3), grad2sftilde(3)
   real(8) :: qsuperball2m2, qsuperball2qsuperball2m1, qsuperballim1, qsuperballim2

   call MatInv3(ori1,ori1i)
   call MatInv3(ori2,ori2i)
   orihelp(1:3,1:3,1) = ori1
   orihelp(1:3,1:3,2) = ori2

   lambda = 0.5d0                 ! initial value of lambda
   lambda = 0.8d0*lambda
   rc = half*(-r21)               ! initial value of rc
   rc = 0.8d0*rc
!   write(*,*) 'HÄR 1', itest
   if (itest == 5) call TestSuperballOverlap_NR(1)
   lambda0 = lambda
   rc0 = rc

   do iter = 1, nitersuperball

      if (itest == 5) call TestSuperballOverlap_NR(2)
      rsb(1:3,1) = matmul(ori1i,rc)
      rsb(1:3,2) = matmul(ori2i,rc+r21)
      qsuperball2m2 = qsuperball2-two
      qsuperball2qsuperball2m1 = qsuperball2*(qsuperball2-one)
      qsuperballim1 = qsuperballi - one
      qsuperballim2 = qsuperballi - two

! ... calculation of shape functions and their two first partial derivatives

      do k = 1, 2
         r1help = abs(rsb(1,k))**qsuperball2m2
         r2help = abs(rsb(2,k))**qsuperball2m2
         r3help = abs(rsb(3,k))**qsuperball2m2
         grad2sftilde(1) = qsuperball2qsuperball2m1*r1help
         grad2sftilde(2) = qsuperball2qsuperball2m1*r2help
         grad2sftilde(3) = qsuperball2qsuperball2m1*r3help
         gradsftilde(1) = qsuperball2*abs(rsb(1,k))*r1help*sign(one,rsb(1,k))
         gradsftilde(2) = qsuperball2*abs(rsb(2,k))*r2help*sign(one,rsb(2,k))
         gradsftilde(3) = qsuperball2*abs(rsb(3,k))*r3help*sign(one,rsb(3,k))
         fac1 = abs(rsb(1,k))**2*r1help + abs(rsb(2,k))**2*r2help + abs(rsb(3,k))**2*r3help
         fac1help = fac1**qsuperballim2
         gp = qsuperballi*radsuperball2i*fac1*fac1help                          ! g'(x)
         gpp = qsuperballi*qsuperballim1*radsuperball2i*fac1help                ! g''(x)
         sf(k) = radsuperball2i*fac1**2*fac1help                                ! shape function
         sfd(1:3,k) = gp*matmul(orihelp(1:3,1:3,k),gradsftilde)                 ! its gradient
         sfdd(1,1,k) = gp*grad2sftilde(1) + gpp*gradsftilde(1)*gradsftilde(1)   ! its second gradient
         sfdd(2,2,k) = gp*grad2sftilde(2) + gpp*gradsftilde(2)*gradsftilde(2)
         sfdd(3,3,k) = gp*grad2sftilde(3) + gpp*gradsftilde(3)*gradsftilde(3)
         sfdd(1,2,k) = gpp*gradsftilde(1)*gradsftilde(2)
         sfdd(1,3,k) = gpp*gradsftilde(1)*gradsftilde(3)
         sfdd(2,3,k) = gpp*gradsftilde(2)*gradsftilde(3)
         sfdd(2,1,k) = sfdd(1,2,k)
         sfdd(3,1,k) = sfdd(1,3,k)
         sfdd(3,2,k) = sfdd(2,3,k)
         sfdd(1:3,1:3,k) = matmul(matmul(orihelp(1:3,1:3,k),sfdd(1:3,1:3,k)),transpose(orihelp(1:3,1:3,k)))
      end do

! ... overlap functions

      of = lambda*sf(1) + (one-lambda)*sf(2)
      ofd(1:3) = lambda*sfd(1:3,1) + (one-lambda)*sfd(1:3,2)
      ofdd(1:3,1:3) = lambda*sfdd(1:3,1:3,1) + (one-lambda)*sfdd(1:3,1:3,2)
      if (itest == 5) call TestSuperballOverlap_NR(3)

! ... NR matter

      ofdd(1,1) = ofdd(1,1) + 1d-10 ! to avoid singular matrix
      ofdd(2,2) = ofdd(2,2) + 1d-10
      ofdd(3,3) = ofdd(3,3) + 1d-10
      call MatInv3(ofdd,minv)                                                           ! minv
      dg(1:3) = sfd(1:3,1) - sfd(1:3,2)                                                 ! dg
      vtemp = matmul(dg,minv)                                                           ! dg*minv
      dll = vtemp(1)*dg(1) + vtemp(2)*dg(2) + vtemp(3)*dg(3)                            ! dll
      dl = ((sf(1)-sf(2)) - (vtemp(1)*ofd(1) + vtemp(2)*ofd(2) + vtemp(3)*ofd(3)))/dll  ! dl
      dr = matmul(minv,(-dg(1:3)*dl-ofd(1:3)))
      if (abs(dl) > dl_cut) dl = sign(dl_cut,dl)+(dl-sign(dl_cut,dl))*dl_damp
      where (abs(dr) > dr_cut) dr = sign(dr_cut,dr)+(dr-sign(dr_cut,dr))*dr_damp
      lambda = lambda + dl
 !    if (iter > 1) then
         rc = rc + dr
 !    end if
      if (itest == 5) call TestSuperballOverlap_NR(4)
      if (max(abs(dl),maxval(abs(dr))) < tolsuperball) then
         loverlap = .true.
         if (of > one) loverlap = .false.         ! of is here the overlap potential
         if (lstatsuperball) call SuperballStatNR(iSimulationStep, iter)
         if (itest == 5) call TestSuperballOverlap_NR(5)
         goto 999
      end if
   end do
   if (iter == nitersuperball + 1) then
      call warn('SuperballOverlap_NR', 'no convergence in NR iterattion', 6)
      call warn('SuperballOverlap_NR', 'no convergence in NR iterattion', uout)
      call TestSuperballOverlap_NR(4)
      stop 1
   end if
 999  continue

contains

!........................................................................

subroutine TestSuperballOverlap_NR(iMode)
   integer(4), intent(in) :: iMode
   integer(4) :: unit
   unit = uout
   if (iMode == 1) then
      call writehead(3,'SuperballOverlap_NR', unit)
      call WriteVec(3,r21,'r21',1,3,1,unit)
      call WriteMat(3,3,ori1,'ori1',1,3,1,1,3,1,unit)
      call WriteMat(3,3,ori2,'ori2',1,3,1,1,3,1,unit)
      call WriteMat(3,3,ori1i,'ori1i',1,3,1,1,3,1,unit)
      call WriteMat(3,3,ori2i,'ori2i',1,3,1,1,3,1,unit)
   else if(iMode == 2) then
      write(unit,*)
      write(unit,'(a,i10)')             'iteration                    =', iter
      write(unit,'(a,f10.5,4x,3f10.5)') 'lambda, rc                   =', lambda, rc(1:3)
   else if(iMode == 3) then
      write(unit,'(a,2f8.3,5x,f8.3)')   'shape & overlap func, value  =', sf(1:2), of
      write(unit,'(a,3f8.3,4x,3f8.3,5x,3f10.3)') 'shape & overlap func, grad   =', sfd(1:3,1:2), ofd(1:3)
      write(unit,'(a,3(3f8.3,4x))')     'shape func, hessian 1        =', sfdd(1:3,1:3,1)
      write(unit,'(a,3(3f8.3,4x))')     'shape func, hessian 2        =', sfdd(1:3,1:3,2)
      write(unit,'(a,3(3f8.3,4x))')     'overlap func, hessian        =', ofdd(1:3,1:3)
   else if(iMode == 4) then
      write(unit,'(a, (3i8  ,4x))')     'istep                        =', istep
      write(unit,'(a,2f8.3,5x,f8.3)')   'overlap func                 =', of
      write(unit,'(a,3(3f8.3,4x))')     'minv                         =', minv
      write(unit,'(a,3f8.3)')           'dg                           =', dg
      write(unit,'(a,3(3f8.3,4x))')     'vtemp                        =', vtemp
      write(unit,'(a,3(3f8.3,4x))')     'dll                          =', dll
      write(unit,'(a,f10.5,4x,3f10.5)') 'dl, dr                       =', dl, dr
      write(unit,'(a,f10.5,4x,3f10.5)') 'lambda, rc                   =', lambda, rc
      write(unit,'(a,3e13.5)')          'max(abs(dl),abs(dr))         =', max(abs(dl),maxval(abs(dr)))
!      write(99,'(a,2i5,2(f10.5,4x,3f10.5,5x))') 'istep, iter, dl, dr, lambda, rc   =', istep, iter, dl, dr, lambda, rc
   else if(iMode == 5) then
      write(unit,*)
      write(unit,'(a,f10.5,4x,3f10.5)') 'dl, dr (total)               =', lambda - lambda0, rc - rc0
      write(unit,'(a,l5,f10.5)') '-> convergence in nr-iteration: loverlap, overlap function',loverlap, of
   end if
end subroutine TestSuperballOverlap_NR

!........................................................................

end function SuperballOverlap_NR

!************************************************************************
!*                                                                      *
!*     SuperballStatNR                                                  *
!*                                                                      *
!************************************************************************

! ... superball check statistics, NR overlap method

!subroutine SuperballStatNR(iStage, rr, iter) !rr is not used
subroutine SuperballStatNR(iStage, iter)
   use MolModule
   implicit none
   integer(4), intent(in) :: iStage
   !real(8), intent(in)     :: rr(3)
   integer(4), intent(in)  :: iter

   character(40), parameter :: txroutine ='SuperballStatNR'
   character(80), parameter :: txheading ='superball: NR overlap method'
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4) :: ivar

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      allocate(var(nvar))

      var(1)%label = 'number of NR-iterations        = '
      var(1)%norm = One

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      ivar = 1
      var(ivar)%value = iter
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

      deallocate(var)

   end select

end subroutine SuperballStatNR

!************************************************************************
!*                                                                      *
!*     SuperballStatMesh                                                *
!*                                                                      *
!************************************************************************

! ... superball check statistics, mesh overlap check

subroutine SuperballStatMesh(iStage, rr, loverlap, time)
   use MolModule
   integer(4), intent(in) :: iStage
   real(8), intent(in)     :: rr(3)
   logical, intent(in)     :: loverlap
   real(8), intent(in)     :: time

   call SuperballAver(iStage, rr, loverlap, time)
   call SuperballDF(iStage, rr, time)

end subroutine SuperballStatMesh

!************************************************************************
!*                                                                      *
!*     SuperballAver                                                    *
!*                                                                      *
!************************************************************************

! ... superball overlap-check average

subroutine SuperballAver(iStage, rr, loverlap, time)
   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   real(8), intent(in)     :: rr(3)
   logical, intent(in)     :: loverlap
   real(8), intent(in)     :: time

   character(40), parameter :: txroutine ='SuperballAver'
   character(80), parameter :: txheading ='superball averages'
   integer(4),           save :: nvar, nbin
   real(8), allocatable, save :: samp(:), sum1(:,:), sum2(:,:)
   real(8),             save :: rmin, rmax, bin, bini
   real(8)    :: r1, vals(4), rlim
   integer(4) :: ibin, i

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      rmin = sqrt(rcut2superball(1))
      rmax = sqrt(rcut2superball(2))
      nbin = 40
      bin = (rmax - rmin)/nbin
      bini = One/bin
      nvar = 4
      allocate(samp(-1:nbin), sum1(nvar,-1:nbin), sum2(nvar,-1:nbin))
      samp = 0.0E+00
      sum1 = 0.0E+00
      sum2 = 0.0E+00
      samp = 0
      sum1 = 0
      sum2 = 0

   case (iSimulationStep)

      r1 = sqrt(dot_product(rr,rr))
      ibin = int(bini*(r1-rmin))
      ibin = max(-1,min(ibin,nbin))
      vals = [1d0*dopc,1d0*tric, merge(1d0,0d0,loverlap),1000*time]
      samp(ibin) = samp(ibin) + One
      sum1(:,ibin) = sum1(:,ibin) + vals
      sum2(:,ibin) = sum2(:,ibin) + vals**2

   case (iAfterSimulation)

      if (master) then
          call WriteHead(2, txheading, uout)
          write(uout,'(a,t35,i10  )') 'nbin                           = ', nbin
          write(uout,'(a,t35,f10.4)') 'rmin                           = ', rmin
          write(uout,'(a,t35,f10.4)') 'rmax                           = ', rmax
          write(uout,'(a,t35,f10.4)') 'bin                            = ', bin
          write(uout,*)
          write(uout,'(a,4x,a,3x,a,3x,a,4x,a,4x,a,5x,a)') &
                       'bin','rlim','sampl','polytrope tests','triangle tests','fraction overlap','time (ms)'
          write(uout,'(a,4x,a,3x,a,3x,a,4x,a,4x,a,5x,a)') &
                       '---','----','-----','---------------','--------------','----------------','---------'
          do i=1,nvar
              where (samp /= 0) sum1(i,:) = sum1(i,:)/samp
              where (samp > 1) sum2(i,:) = sqrt(max(Zero,sum2(i,:)/samp - sum1(i,:)**2))
          end do
          do ibin = -1, nbin
              rlim = rmin+(ibin+Half)*bin
              write(uout,'(i4,f8.4,i8,2(2f8.1,2x),(2f8.3,2x),2f8.3)') &
                           ibin, rlim, int(samp(ibin)), (sum1(i,ibin),sum2(i,ibin),i = 1,nvar)
          end do
       end if

    end select

end subroutine SuperballAver

!************************************************************************
!*                                                                      *
!*     SuperballDF                                                      *
!*                                                                      *
!************************************************************************

! ... superball overlap-check distribution functions

!     type  label    quantity
!     ----  -----    --------
!     1     polytrope  dop (polytrope)
!     2     triangle   tri (triangle)
!     3     time       cpu time

!subroutine SuperballDF(iStage, rr, loverlap, time) !loverlap is not used
subroutine SuperballDF(iStage, rr, time)
   use MolModule
   implicit none
   integer(4), intent(in) :: iStage
   real(8), intent(in)     :: rr(3)
   !logical, intent(in)     :: loverlap
   real(8), intent(in)     :: time

   character(40), parameter :: txroutine ='SuperballDF'
   character(80), parameter :: txheading ='superball distribution function'
   integer(4)   , parameter :: ntype = 3
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8), save :: radmin, radmax, radbin, radbini
   integer(4), save :: nrad
   integer(4) :: itype, ivar, ibin, irad
   real(8)    :: r1
   character(1) :: str

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.true.
      vtype%max = [2000d0, 200d0, 1.0d0]
      vtype%min = [Zero,  Zero, Zero]
      vtype%nbin= 100
      nrad = 2

      radmin = sqrt(rcut2superball(1))
      radmax = sqrt(rcut2superball(2))
      radbin = (radmax - radmin)/nrad
      radbini = one/radbin

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['polytrope ','triangle  ','time (ms) ']
      vtype%nvar = nrad

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nvar,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do irad = 1, vtype(itype)%nvar
               ivar = ivar+1
               ipnt(irad,itype) = ivar
               write(str,'(i1)') irad
               var(ivar)%label = trim(vtype(itype)%label)//' '//str
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
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

      r1 = sqrt(dot_product(rr,rr))
      irad = 1+int(floor(radbini*(r1-radmin)))
      if(irad > 0 .and. irad < nrad+1) then
         var%nsamp2 = var%nsamp2 + 1

! ... sample type 1

         itype = 1
         if (vtype(itype)%l) then
            ivar = ipnt(irad,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(dopc-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + one
         end if

! ... sample type 2

          itype = 2
          if (vtype(itype)%l) then
             ivar = ipnt(irad,itype)
             ibin = max(-1,min(floor(var(ivar)%bini*(tric-var(ivar)%min)),int(var(ivar)%nbin)))
             var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + one
          end if

! ... sample type 3

         itype = 3
         if (vtype(itype)%l) then
            ivar = ipnt(irad,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(1000*time-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + one
         end if
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f8.2)') 'lower separation distance      = ', sqrt(rcut2superball(1))
         write(uout,'(a,t35,f8.2)') 'upper separation distance      = ', sqrt(rcut2superball(2))
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var, ipnt)

   end select

end subroutine SuperballDF
