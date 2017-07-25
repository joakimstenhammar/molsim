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

! ... lib calls: func  : Erf        , Erf2           , GammaLn     ,
!                        PL         , PLM            , CCLM        , SetWW3j     , WW3j
!                integ : Trap       , TrapNE         , IntegCirPoints
!                lsfit : PolFit     , PolVal
!                matrix: MatVecMul  , VecMatMul      , MatInv2     , MatInv3     , MatInv      , Diag
!                smooth: Smooth     , CsEval         , Vinter      , LinInter    , CardinalBSpline
!                sort  : HeapSort   , HeapSortIndex  , MakeCluster , CalcClusterMember
!                transf: CarToSph   , SphToCar       , CarToStd1   , OriToEuler  , EulerToQua
!                        OriToQua   , QuaToOri       , AxisAngToOri,
!                        LabToPri   , PriToLab       , EulerRot    ,
!                        CheckOriOrtho, QuaNorm,     , AngVelToQuaVel, QuaVelToAngVel,
!               random : Random     , GauRandom      , CirRandom   , SphRandom
!                mixed : InvInt     , InvFlt         , RelDiff     , CpuAdd      , CpuLeft     , CpuTot     , SecondsSinceStart
!            readwrite : WriteVec   , WriteMat       , WriteStd    , WriteFront  , WriteHead   , WriteDateTime, WriteIOStat, Warn       , Stop
!               string : Center     , SpaceOut       , LowerCase   , UpperCase   , SubStr      , Advance    , ReplaceText
!                 plot : SignMagn   , Plot
!               powell : BrentMod

module MollibModule !Starting to migrate to Module

   implicit none
   private
   public InvInt, Center, SpaceOut, ReplaceText

   interface InvInt
      module procedure InvInt_single, InvInt_double
   end interface InvInt

   contains
!************************************************************************
!*                                                                      *
!*     InvInt                                                           *
!*                                                                      *
!************************************************************************

! ... return the inverse of an integer

   pure elemental real(8) function InvInt_single(i)
      implicit none
      integer(4), intent(in) :: i
      InvInt_single = 0.0d0
      if (i .ne. 0) InvInt_single = 1.0d0/real(i)
   end function InvInt_single

   pure elemental real(8) function InvInt_double(i)
      implicit none
      integer(8), intent(in) :: i
      InvInt_double = 0.0d0
      if (i .ne. 0) InvInt_double = 1.0d0/real(i)
   end function InvInt_double

!************************************************************************
!*                                                                      *
!*     Center                                                           *
!*                                                                      *
!************************************************************************

! ... center a string

   pure function Center(nw,string)
      implicit none
      integer(4),   intent(in) :: nw              ! width of line
      character(*), intent(in) :: string          ! string
      character(len=nw) :: Center
      integer(4) :: nchar, noff
      Center = repeat(' ', nw)
      nchar = len_trim(string)
      noff = (nw-nchar)/2
      Center(noff+1:noff+nchar) = string
   end function Center

!************************************************************************
!*                                                                      *
!*     SpaceOut                                                         *
!*                                                                      *
!************************************************************************

! ... space out a string

   pure function SpaceOut(string)
      implicit none
      character(*), intent(in) :: string   ! string to be spaced out
      character(:), allocatable :: SpaceOut
      integer(4) :: i
      allocate( character(len=2*len(string)) :: SpaceOut)
      SpaceOut=repeat(' ', len(SpaceOut))
      do i=1,len_trim(adjustl(string))
         SpaceOut(2*i-1:2*i)=string(i:i)//' '
      end do
   end function SpaceOut

!************************************************************************
!*                                                                      *
!*     ReplaceText                                                      *
!*                                                                      *
!************************************************************************

! ... Replace text in a string
   pure function ReplaceText (string,toReplace,replaceWith)  result(outString)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: toReplace
      character(len=*), intent(in) :: replaceWith
      character(len=:), allocatable :: outString

      integer(4) :: lenToReplace, lenReplaceWith, lenOutString
      integer(4) :: posString, posOutString
      integer(4) :: pos, nrep

      lenToReplace = len(toReplace)
      lenReplaceWith = len(replaceWith)

      posString = 1
      posOutString = 1

      nrep = 0
      do
         pos = index(string(posString:),toReplace) - 1
         if (pos == -1) then
            exit
         endif
         nrep = nrep + 1
         posString = posString + pos + lenToReplace
      end do

      lenOutString = len(string) + nrep * (lenReplaceWith - lenToReplace)
      allocate(character(len=lenOutString) :: outString)
      outString = repeat(" ", lenOutString)

      posString = 1
      posOutString = 1
      do
         pos = index(string(posString:),toReplace) - 1
         if (pos == -1) then
            exit
         endif
         outString(posOutString:(posOutString + pos + lenReplaceWith - 1)) = string(posString:(posString + pos - 1))//replaceWith
         posString = posString + pos + lenToReplace
         posOutString = posOutString + pos + lenReplaceWith
      end do
      outString(posOutString:) = string(posString:)
   end function ReplaceText

end module MollibModule

!************************************************************************
!*                                                                      *
!*     ErfLocal                                                         *
!*                                                                      *
!************************************************************************

! ... return the error function, erf(x)

real(8) function ErfLocal(x)
   implicit none
   real(8), intent(in) :: x
   real(8) :: Erf2

!   ErfLocal = erf(x)       ! intrinsic ≈ 1e-10
!   ErfLocal = Erf1(x)      ! abs(dev) < 3e-7
   ErfLocal = Erf2(x)      ! abs(dev) ≈< 1e-16

end function ErfLocal

!************************************************************************
!*                                                                      *
!*     Erf1                                                             *
!*                                                                      *
!************************************************************************

! ... return the error function, erf(x)
!     abs(dev) < 3 e-7
!     abramowitz and stegun, handbook of mathematical functions

real(8) function Erf1(x)
   implicit none
   real(8), intent(in) :: x
   Erf1 = 1.0d0 - ( 1.0d0 + x * (0.0705230784d0 +                &
                            x * (0.0422820123d0 +                &
                            x * (0.0092705272d0 +                &
                            x * (0.0001520143d0 +                &
                            x * (0.0002765672d0 +                &
                            x *  0.0000430638d0 ))))))**(-16)
end function Erf1

!************************************************************************
!*                                                                      *
!*     Erf2                                                             *
!*                                                                      *
!************************************************************************

! ... return the error function, erf(x)
!     abs(dev) < (ca) 1e-16
!     originating from mark 4 release nag 1974, ...

real(8) function Erf2(x)

   real(8), intent(in) :: x

   real(8), parameter  :: sqrtpi = 1.7724538509055160d0
   real(8), parameter  :: zero   = 0.0d0, half   = 0.5d0, one    = 1.0d0, &
                          two    = 2.0d0, three  = 3.0d0, twenty =20.0d0

   real(8), parameter  :: xup = 6.25d0
   integer(4)          :: ncfc  = 18
   integer(4)          :: ncfd  = 17
   real(8)             :: c(18) =                                             &
         [ 1.9449071068178803d0,  4.20186582324414d-2,  -1.86866103976769d-2, &
           5.1281061839107d-3,   -1.0683107461726d-3,    1.744737872522d-4,   &
          -2.15642065714d-5,      1.7282657974d-6,      -2.00479241d-8,       &
          -1.64782105d-8,         2.0008475d-9,          2.57716d-11,         &
          -3.06343d-11,           1.9158d-12,            3.703d-13,           &
           -5.43d-14,             -4.0d-15,               1.2d-15       ]

   real(8)             :: d(17) =                                             &
         [ 1.4831105640848036d0, -3.010710733865950d-1,  6.89948306898316d-2, &
          -1.39162712647222d-2,   2.4207995224335d-3,   -3.658639685849d-4,   &
           4.86209844323d-5,     -5.7492565580d-6,       6.113243578d-7,      &
          -5.89910153d-8,         5.2070091d-9,         -4.232976d-10,        &
           3.18811d-11,          -2.2361d-12,            1.467d-13,           &
          -9.0d-15,               5.0d-16                              ]

   integer(4)          :: j
   real(8)             :: bj, bjp1, bjp2, x2, xv

   xv = abs(x)

   if (xv <= two) then
      x2 = x**2 - two
      bjp2 = zero
      bjp1 = d(ncfd)
      do j = ncfd - 1, 2, -1
         bj = x2*bjp1 - bjp2 + d(j)
         bjp2 = bjp1
         bjp1 = bj
      end do
      bj = x2*bjp1 - bjp2 + d(1)
      Erf2 = half*(bj-bjp2)*x
   elseif (two < xv .and. xv < xup) then
      x2 = two - twenty/(xv+three)
      bjp2 = zero
      bjp1 = c(ncfc)
      do j = ncfd - 1, 2, -1
         bj = x2*bjp1 - bjp2 + c(j)
         bjp2 = bjp1
         bjp1 = bj
      end do
      bj = x2*bjp1 - bjp2 + c(j)
      x2 = half*(bj-bjp2)/xv*exp(-x*x)/sqrtpi
      Erf2 = (one-x2)*sign(one,x)
   elseif (xv > xup) then
      Erf2 = sign(one,x)
   endif

end function Erf2

!************************************************************************
!*                                                                      *
!*     GammaLn                                                          *
!*                                                                      *
!************************************************************************

! ... return the logarithm of the gamma function
!     relative error < 2 e-10 for xx > 1
!     lanczos, j. s.i.a.m. numerical analysis, ser b, vol 1, pg 86, 1964

real(8) function GammaLn(xx)
   implicit none
   real(8), parameter :: half = 0.5d0, zero = 0.0d0, one = 1.0d0, fpf = 5.5d0
   real(8), parameter :: cof(6) = [ 76.18009173d0, -86.50532033d0,  24.01409822d0, &
                                   -1.231739516d0, 0.120858003d-2, -0.536382d-5 ]
   real(8), parameter :: stp = 2.50662827465d0
   real(8), intent(in) :: xx
   integer(4) :: i
   real(8) :: x, temp, sum
   if (xx <= zero) call stop('GammaLn','xx <= 0.0',6)
   if (xx == one .or. xx == 2.0d0) then
      GammaLn = zero
   else
      x = xx-one
      temp = x+fpf
      temp = (x+half)*log(temp)-temp
      sum = one
      do i = 1,6
         x = x+one
         sum = sum+cof(i)/x
      end do
      GammaLn = temp+log(stp*sum)
   end if
end function GammaLn

!************************************************************************
!*                                                                      *
!*     PL                                                               *
!*                                                                      *
!************************************************************************

! ... return the legendre polynomial p(l) at x.
!     recurrence rel: (l+1)*p(l+1)(x) = (2l+1)*x*p(l)(x) - l*p(l-1)(x)

real(8) function PL(l,x)
   implicit none
   real(8), parameter :: one = 1.0d0 , two = 2.0d0
   integer(4), intent(in) :: l
   real(8),    intent(in) :: x
   integer(4) :: i
   real(8)    :: p0, p1, p2

   if (l < 0)        call stop('PL','l < 0',6)
   if (abs(x) > 1.0) call stop('PL','abs(x) > 1.0',6)
   if (l == 0) then
      PL = one
   else if (l == 1) then
      PL = x
   else if (l >= 2) then
      p0 = one
      p1 = x
      do i = 2,l
         p2 = ((two*i-one)*x*p1-(i-one)*p0)/i
         p0 = p1
         p1 = p2
      end do
      PL = p2
   end if
end function PL

!************************************************************************
!*                                                                      *
!*     PLM                                                              *
!*                                                                      *
!************************************************************************

! ... return the associate legendre polynomial p(l,m) at x.
!     p(m,m)(x) = (2m-1)!!(1-x*x)*(m/2) and the recurrence relation
!     (l-m+1)*p(l+1,m)(x)=(2*l+1)*x*p(l,m)(x)-(l+m)*p(l-1,m)(x)

real(8) function PLM(l,m,x)
   implicit none
   real(8), parameter :: one = 1.0d0 , two = 2.0d0
   integer(4), intent(in) :: l, m
   real(8),    intent(in) :: x
   integer(4) :: i, j
   real(8)    :: pmm, factor, fac, pmp1m, pjm

   if (l < 0) call stop('PLM','l < 0',6)
   if (m > l) call stop('PLM','m > l',6)
   if (m < 0) call stop('PLM','m < 0',6)
   if (abs(x) > one) call stop('PLM','abs(x) > one',6)
   pmm = one
   if (m > 0) then
      factor = one
      fac = sqrt((one-x)*(one+x))
      do i = 1,m
         pmm = pmm*factor*fac
         factor = factor+two
      end do
   end if
   if (l == m) then
      PLM = pmm
   else
      pmp1m = (2*m+1)*x*pmm
      if (l == m+1) then
         PLM = pmp1m
      else
         do j = m+2,l
            pjm = (x*(2*j-1)*pmp1m-(j+m-1)*pmm)/(j-m)
            pmm = pmp1m
            pmp1m = pjm
         end do
         PLM = pjm
      end if
   end if
end function PLM

!************************************************************************
!*                                                                      *
!*     CCLM                                                             *
!*                                                                      *
!************************************************************************

! ... return the spherical harmonics
!             norm = 0,  modified spherical harmonics c(l,m)(theta,phi)
!             nrom = 1,  spherical harmonics y(l,m)(theta,phi)

complex(8) function CCLM(l,m,theta,phi,norm)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   real(8), parameter :: pi = 3.1415926535897932d0 , facpi = 4.0d0*pi
   integer(4), parameter :: lmax = 150, l2max = 2*lmax      !if lmax is larger there is a floating overflow
   integer(4), intent(in) :: l, m, norm
   real(8),    intent(in) :: theta, phi
   logical, save :: first=.true.
   real(8), save :: s(0:l2max), si(0:l2max), vnorm(0:lmax)
   integer(4) :: i, mabs
   real(8) :: x, PLM

   if (l > lmax) call stop('CCLM','l > lmax',6)
   if (l < 0)    call stop('CCLM','l < 0',6)
   if (m > l)    call stop('CCLM','m > l',6)
   if (first) then
      s(0) = one
      si(0) = one
      do i = 1, l2max
         s(i) = s(i-1)*sqrt(dble(i))     ! sqrt(n!)
         si(i) = one/s(i)                ! 1/sqrt(n!)
      end do
      do i = 0, lmax
         vnorm(i) = sqrt((2*i+one)/facpi)
      end do
      first = .false.
   end if
   mabs = abs(m)
   x = cos(theta)
   CCLM = (-1)**mabs*s(l-mabs)*si(l+mabs)*PLM(l,mabs,x)*cdexp(dcmplx(zero,mabs*phi))
   if (m < 0) CCLM = (-1)**mabs*conjg(CCLM)
   if (norm == 1) CCLM = vnorm(l)*CCLM
end function CCLM

!************************************************************************
!*                                                                      *
!*     DDJKM                                                            *
!*                                                                      *
!************************************************************************

! ... return one wigner's rotation matrix element
!      norm = 0, wigner's rotation matrix element d(jkm)(alfa,beta,gamma)
!      norm = 1, normalized wigner's rotation matrix element d(jkm)(alfa,beta,gamma)
!     definition according to Edmonds

complex(8) function DDJKM(j,k,m,alfa,beta,gamma,norm)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0, two = 2.0d0
   real(8), parameter :: pi = 3.1415926535897932d0 , facpi = 8.0d0*pi**2
 !!!  integer(4), parameter :: jmax = 5
   integer(4), parameter :: jmax = 1
   integer(4), intent(in) :: j, k, m, norm
   real(8),    intent(in) :: alfa, beta, gamma
   logical, save :: first = .true.
   real(8), save :: vnorm(0:jmax), del(0:jmax,-jmax:jmax,-jmax:jmax)
   real(8)    :: temp
   integer(4) :: imod, n

   if (j > jmax) call stop('DDJKM','j > jmax', 6)
   if (k > j   ) call stop('DDJKM','k > j'   , 6)
   if (m > j   ) call stop('DDJKM','m > j'   , 6)

   if (first) then
      do n = 0,jmax
         vnorm(n) = sqrt((2*n+one)/facpi)
      end do
      call CalcDel(jmax,del)
      first = .false.
   end if

   imod = mod(4+k-m,4)
   if (imod == 0) then
      temp = +del(j,0,k)*del(j,0,m)
      do n = 1,j
         temp = temp+two*del(j,n,k)*del(j,n,m)*(+cos(n*beta))
      end do
   else if (imod == 1) then
      temp = zero
      do n = 1,j
         temp = temp+two*del(j,n,k)*del(j,n,m)*(+sin(n*beta))
      end do
   else if (imod == 2) then
      temp = -del(j,0,k)*del(j,0,m)
      do n = 1,j
         temp = temp+two*del(j,n,k)*del(j,n,m)*(-cos(n*beta))
      end do
   else if (imod == 3) then
      temp = zero
      do n = 1,j
         temp = temp+two*del(j,n,k)*del(j,n,m)*(-sin(n*beta))
      end do
   end if
   ddjkm = cdexp(dcmplx(zero,k*gamma+m*alfa))*dcmplx(temp,zero)
   if (norm == 1) ddjkm = vnorm(j)*ddjkm

contains

!........................................................................

! ... calculate the del matrix

subroutine CalcDel(jmax,del)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   integer(4), intent(in) :: jmax
   real(8),   intent(out) :: del(0:jmax,-jmax:jmax,-jmax:jmax)
   real(8)    :: dum(0:2*jmax,-2*jmax:2*jmax,-2*jmax:2*jmax)
   real(8)    :: fac(0:20), s(0:20), si(0:20), ffac, sum
   integer(4) :: i, j, jj, k, kk, m, mm

   if (jmax > 10) call stop('CalcDel','jmax > 10', 6)

! ... set fac, s, and si

   fac(0) = one
   s(0) = zero
   ffac = one
   do i = 1,20
      ffac = ffac*dble(i)
      fac(i) = ffac
      s(i) = sqrt(dble(i))
      si(i) = one/s(i)
   end do

! ... generate dum

    dum(0,0,0) = one
    do jj = 1,2*jmax
       mm = jj
       do  kk = jj,-jj,-2
          dum(jj,kk,mm) = (-1)**((jj-kk)/2)*sqrt(fac(jj)/(fac((jj+kk)/2)*fac((jj-kk)/2)*2**jj))
       end do
       do mm = jj-2,-jj,-2
          do kk = jj,-jj,-2
             sum = zero
             if (kk+1 <= jj-1) sum = sum+s((jj-kk)/2)*si((jj-mm)/2)*dum(jj-1,kk+1,mm+1)*si(2)
             if (abs(kk-1) <= jj-1) sum = sum+s((jj+kk)/2)*si((jj-mm)/2)*dum(jj-1,kk-1,mm+1)*si(2)
             dum(jj,kk,mm) = sum
          end do
       end do
    end do

! ... skip half integer indices

    do j = 0,jmax
       jj = 2*j
       do k = j,-j,-1
          kk = 2*k
          do m = j,-j,-1
             mm = 2*m
             del(j,k,m) = dum(jj,kk,mm)
         end do
       end do
    end do

end subroutine CalcDel

!........................................................................

end function DDJKM

!************************************************************************
!*                                                                      *
!*     SetW3J                                                           *
!*                                                                      *
!************************************************************************

! ... precalculate wigner's 3j-symbols

subroutine SetW3J(j1max,j2max,j3max,w3j)
   implicit none
   integer(4), intent(in)  :: j1max, j2max, j3max
   real(8),    intent(out) :: w3j(0:j1max,0:j2max,0:j3max,-j1max:j1max,-j2max:j2max,-j3max:j3max)
   integer(4) :: j1, j2, j3, m1, m2, m3
   real(8)    :: WW3J

   do j1 = 0,j1max
      do j2 = 0,j2max
         do j3 = 0,j3max
            do m1 = -j1,j1
               do m2 = -j2,j2
                  m3 = -(m1+m2)
                  if (abs(m3) > j3) cycle
                  w3j(j1,j2,j3,m1,m2,m3) = WW3J(j1,j2,j3,m1,m2,m3)
               end do
            end do
         end do
      end do
   end do
end subroutine SetW3J

!************************************************************************
!*                                                                      *
!*     WW3J                                                             *
!*                                                                      *
!************************************************************************

! ... return one wigner 3j-symbol

!     OK with (j1,j2,j3) = (45,45,90)
!     overflow with (j1,j2,j3) = (46,46,92)

real(8) function WW3J(j1,j2,j3,m1,m2,m3)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   integer(4), parameter :: jjjmax = 181
   integer(4), intent(in) :: j1, j2, j3, m1, m2, m3
   real(8), save :: s(0:jjjmax), si(0:jjjmax)
   logical, save :: first=.true.
   integer :: i, n, n1, n2, n3, n4, n5, nupp, nlow, isign
   real(8) :: term1, term2, term3

   if (j1+j2+j3+1 > jjjmax) call stop('WW3J','j1+j2+j3+1 > jjjmax',6)

   if (first) then
      s(0) = one
      si(0) = one
      do i = 1,jjjmax
         s(i) = s(i-1)*sqrt(dble(i))     ! sqrt(n!)
         si(i) = one/s(i)                ! 1/sqrt(n!)
      end do
      first = .false.
   end if

   WW3J = zero
   if (j1 < 0) return
   if (j2 < 0) return
   if (j3 < 0) return
   if ((abs(j1-j2) > j3) .or. (j1+j2 < j3)) return
   if ((abs(j2-j3) > j1) .or. (j2+j3 < j1)) return
   if ((abs(j3-j1) > j2) .or. (j3+j1 < j2)) return
   if (m1+m2+m3 /= 0) return
   if (abs(m1) > j1) return
   if (abs(m2) > j2) return
   if (abs(m3) > j3) return

   term1 = s(j1+j2-j3) * s(j2+j3-j1) * s(j3+j1-j2) * si(j1+j2+j3+1)

   term2 = s(j1+m1) * s(j1-m1) * s(j2+m2) * s(j2-m2) * s(j3+m3) * s(j3-m3)

   n1 = j1+j2-j3
   n2 = j1-m1
   n3 = j2+m2
   n4 = j3-j2+m1
   n5 = j3-j1-m2
   nupp = min(n1,n2,n3)
   nlow = max(0,-n4,-n5)
   term3 = zero
   do n = nlow,nupp
      term3 = term3 + (-1)**n * (si(n) * si(n1-n) * si(n2-n) * si(n3-n) * si(n4+n) * si(n5+n))**2
   end do
   isign = -1
   if (abs(mod(j1-j2-m3,2)) == 0) isign = 1
   term3 = isign * term3

   WW3J = term1 * term2 * term3

!  write(55,'(a,6i4,4g15.5)') 'j1,j2,j3,m1,m2,m3,WW3J',j1,j2,j3,m1,m2,m3,WW3J, term1, term2, term3
end function WW3J

!************************************************************************
!*                                                                      *
!*     Trap                                                             *
!*                                                                      *
!************************************************************************

! ... integration by the trapezoidal rule

real(8) function Trap(n,y,h)
   implicit none
   integer(4), intent(in) :: n      ! extent of y
   real(8),    intent(in) :: y(n)   ! function to integrate
   real(8),    intent(in) :: h      ! step size
   Trap = h*(sum(y(2:n-1))+0.5d0*(y(1)+y(n)))
end function Trap

!************************************************************************
!*                                                                      *
!*     TrapNe                                                           *
!*                                                                      *
!************************************************************************

! ... integration by the trapezoidal rule with varying steplength

real(8) function TrapNe(n,x,y)
   implicit none
   integer(4), intent(in) :: n      ! extent of x and y
   real(8), intent(in)    :: x(n)   ! independent variable
   real(8), intent(in)    :: y(n)   ! function to integrate
   TrapNe = 0.5d0*sum(abs(x(2:n)-x(1:n-1))*(y(2:n)+y(1:n-1)))
end function TrapNe

!************************************************************************
!*                                                                      *
!*     IntegCirPoints                                                   *
!*                                                                      *
!************************************************************************

! ... get points and weight for numerical integration over a circle
!     ref: abramowitz and stegun

subroutine IntegCirPoints(npoint,radius,x,y,w)

   implicit none

   real(8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0, six = 6.0d0, nine = 9.0d0, ten = 10.0d0
   real(8), parameter :: pi = 3.1415926535897932d0

   integer(4), intent(in)  :: npoint     ! number of points
   real(8),    intent(in)  :: radius     ! radius of circle
   real(8),    intent(out) :: x(npoint)  ! x-coordinate
   real(8),    intent(out) :: y(npoint)  ! y-coordinate
   real(8),    intent(out) :: w(npoint)  ! weight

   integer(4) :: i
   real(8)    :: k, fac1, fac2, fac3

   if (npoint == 21) then

      x(1) = zero
      y(1) = zero
      w(1) = one/nine

      fac1 = two*pi/ten

      fac2 = sqrt((six-sqrt(six))/ten)
      fac3 = (16.0d0+sqrt(six))/360
      do i = 2, 11
          k = i - 1
          x(i) = radius*fac2*cos(fac1*k)
          y(i) = radius*fac2*sin(fac1*k)
          w(i) = fac3
      end do

      fac2 = sqrt((six+sqrt(six))/ten)
      fac3 = (16.0d0-sqrt(six))/360
      do i = 12, 21
         k = i - 11
         x(i) = radius*fac2*cos(fac1*k)
         y(i) = radius*fac2*sin(fac1*k)
         w(i) = fac3
      end do

   else

      call stop('IntegCirPoints','unsupported value of npoint',6)

   endif

end subroutine IntegCirPoints

!************************************************************************
!*                                                                      *
!*     PolFit                                                           *
!*                                                                      *
!************************************************************************

! ... calculate least square fit to a polynomial

!                        ndata
!       minimization of  sum (y(i)-(a(0)+a(1)*x(i)+...+a(npol)*x(i)**npol
!                        i=1

!       whith respect of a(0),a(1),...,a(npol) gives the system of linea

!       npol     ndata             npol
!       sum a(l) sum x(i)**(k+l) = sum x(i)**(k)*y(i)  k=0,1,...,npol
!       l=0      i=1               i=1

!                <--- t(k,l) -->   <----- c(k) ----->

!       or                t(k,l)*a(l) = c(k)

subroutine PolFit(npol, ndata, x, y, w, iprint, unit, a, dmax, sd)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   integer(4), intent(in)    :: npol      ! degree of the polynomial
   integer(4), intent(in)    :: ndata     ! number of data points
   real(8),    intent(inout) :: x(ndata)  ! x-values
   real(8),    intent(inout) :: y(ndata)  ! y-values
   real(8),    intent(in)    :: w(ndata)  ! weight
   integer(4), intent(in)    :: iprint    ! 0 no output
                                          ! print dmax and sd
                                          ! full output
   integer(4), intent(in)    :: unit      ! output unit
   real(8),    intent(inout) :: a(0:npol) ! coefficients of the polynomial
   real(8),    intent(out)   :: dmax      ! maximal deviation
   real(8),    intent(out)   :: sd        ! standard deviation

   integer(4) :: i, j, k, l, index(2,npol+1)
   real(8) :: t(0:npol,0:npol), c(0:npol), det
   real(8) :: xshift, yshift, facx, facxx, facy, wsumi, PolVal, diff, value

   if (npol > ndata-1)  call stop('PolFit','npol > ndata-1',unit)

! ... shift x and y

   xshift = sum(x(1:ndata))/ndata
   yshift = sum(y(1:ndata))/ndata
   x(1:ndata) = x(1:ndata)-xshift
   y(1:ndata) = y(1:ndata)-yshift

! ... calculate the t matrix and the c vector

   t(0:npol,0     ) = zero
   t(npol,  0:npol) = zero
   c(0:npol) = zero
   do i = 1,ndata
       facx = one*w(i)
       facy = y(i)*w(i)
       do k = 0,npol
          t(k,0) = t(k,0)+facx
          c(k)   = c(k)  +facy
          facx   = facx  *x(i)
          facy   = facy  *x(i)
      end do
      do k = 1,npol
         t(npol,k) = t(npol,k)+facx
         facx = facx  *x(i)
      end do
   end do
   do k = 1,npol
      do l = 1,k
         t(k-l,l) = t(k,0)
         if (k /= npol) t(npol-l,npol+l-k) = t(npol,npol-k)
      end do
   end do
   wsumi = one/sum(w(1:ndata))
   t(0:npol,0:npol) = t(0:npol,0:npol)*wsumi
   c(0:npol) = c(0:npol)*wsumi

! ... solve the system of linear equations

   call MatInv(npol+1, t(0,0), 1, c(0), det)

! ... get the coefficient for the nonshifted data

   index(1,1) = 1
   facx = one
   a(0) = yshift+c(0)
   do i = 1,npol
      index(2,1  ) = 1
      index(2,i+1) = 1
      facx = facx*(-xshift)
      a(0) = a(0)+c(i)*facx
      a(i) = c(i)
      facxx = 1.0
      do j = i-1,1,-1
         index(2,j+1) = index(1,j)+index(1,j+1)
         facxx = facxx*(-xshift)
         a(j) = a(j)+index(2,j+1)*c(i)*facxx
      end do
      do j = 0,i
         index(1,j+1) = index(2,j+1)
      end do
   end do

! ... reshift x and y

   x(1:ndata) = x(1:ndata)+xshift
   y(1:ndata) = y(1:ndata)+yshift

! ... calculate dmax, sd, and print

   if (iprint == 2) write(unit,'(/,t10,a,t25,a,t35,a,t50,a,t70,a/)') 'x','y','weight','y(polyn)','delta'

   dmax = zero
   sd = zero
   do i = 1,ndata
      value = PolVal(npol,a,x(i))
      diff = value-y(i)
      dmax = max(dmax,abs(diff))
      sd = sd+diff**2
      if (iprint == 2)write(unit,'(5g15.5)')x(i),y(i),w(i),value,diff
   end do
   if (ndata == npol+1) sd = zero
   if (ndata > npol+1) sd = sqrt(sd/(ndata-npol-1))
   if (iprint >= 1) then
      write(unit,'(a,t25,5g12.5)') 'maximum deviation   =',dmax
      write(unit,'(a,t25,5g12.5)') 'standard deviation  =',sd
      write(unit,'(a,t25,5g12.5)') 'polyn. coeff. =',a(0:npol)
   end if

end subroutine PolFit

!************************************************************************
!*                                                                      *
!*     PolVal                                                           *
!*                                                                      *
!************************************************************************

! ... calculate the value of the polynomial for a given x-value

real(8) function PolVal(npol,a,x)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   integer(4), intent(in) :: npol      ! degree of polynomial
   real(8),    intent(in) :: a(0:npol) ! coefficients of the polynomial
   real(8),    intent(in) :: x         ! x-value
   integer(4) :: k
   real(8) :: facx

   PolVal = zero
   facx = one
   do k = 0,npol
      PolVal = PolVal+facx*a(k)
      facx = facx*x
   end do
end function PolVal

!************************************************************************
!*                                                                      *
!*     MatVecMul                                                        *
!*                                                                      *
!************************************************************************

! ... multiply a matrix with a vector, c = a*b

subroutine MatVecMul(n1,n2,a,b,c)
   implicit none
   integer(4), intent(in)  :: n1, n2
   real(8),    intent(in)  :: a(n1,n2), b(n2)
   real(8),    intent(out) :: c(n1)
   integer(4) :: i
   do i = 1, n1
      c(i) = dot_product(a(i,1:n2),b(1:n2))
   end do
end subroutine MatVecMul

!************************************************************************
!*                                                                      *
!*     VecMatMul                                                        *
!*                                                                      *
!************************************************************************

! ... multiply a vector with a matrix, c = b*a

subroutine VecMatMul(n1,n2,b,a,c)
   implicit none
   integer(4), intent(in)  :: n1, n2
   real(8),    intent(in)  :: b(n1), a(n1,n2)
   real(8),    intent(out) :: c(n2)
   integer(4) :: i
   do i = 1, n1
      c(i) = dot_product(b(1:n1),a(1:n1,i))
   end do
end subroutine VecMatMul

!************************************************************************
!*                                                                      *
!*     MatInv2                                                          *
!*                                                                      *
!************************************************************************

! ... matrix inversion of a 2x2 matrix

subroutine MatInv2(a,b)
   implicit none
   real(8), intent(in)  :: a(2,2)
   real(8), intent(out) :: b(2,2)
   real(8) :: detinv
   detinv = 1.0d0/(a(1,1)*a(2,2)-a(1,2)*a(2,1))
   if (detinv == 0.d0) call Stop('MatInv2', 'singular determinant', 6)
   b(1,1) = a(2,2)*detinv
   b(1,2) = -a(1,2)*detinv
   b(2,1) = -a(2,1)*detinv
   b(2,2) = a(1,1)*detinv
end subroutine MatInv2

!************************************************************************
!*                                                                      *
!*     MatInv3                                                          *
!*                                                                      *
!************************************************************************

! ... matrix inversion of a 3x3 matrix

subroutine MatInv3(a,b)
   implicit none
   real(8), intent(in)  :: a(3,3)
   real(8), intent(out) :: b(3,3)
   real(8) :: detinv
   detinv = a(1,1)*a(2,2)*a(3,3) + a(2,1)*a(3,2)*a(1,3) + a(3,1)*a(1,2)*a(2,3)  &
                  -a(3,1)*a(2,2)*a(1,3) - a(2,1)*a(1,2)*a(3,3) - a(1,1)*a(3,2)*a(2,3)
   if (detinv == 0.d0) call Stop('MatInv3', 'singular determinant', 6)
   detinv = 1.0d0/detinv
   b(1,1) = (a(2,2)*a(3,3) - a(2,3)*a(3,2))*detinv
   b(1,2) = (a(1,3)*a(3,2) - a(1,2)*a(3,3))*detinv
   b(1,3) = (a(1,2)*a(2,3) - a(1,3)*a(2,2))*detinv
   b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))*detinv
   b(2,2) = (a(1,1)*a(3,3) - a(1,3)*a(3,1))*detinv
   b(2,3) = (a(1,3)*a(2,1) - a(1,1)*a(2,3))*detinv
   b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))*detinv
   b(3,2) = (a(1,2)*a(3,1) - a(1,1)*a(3,2))*detinv
   b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))*detinv
end subroutine MatInv3

!************************************************************************
!*                                                                      *
!*     MatInv                                                           *
!*                                                                      *
!************************************************************************

! ... matrix inversion with solution of linear equations

!     routine solves matrix equation a*x = b where
!         in: a(n,n)
!             b(n,m)
!        out: x(n,m) is returned in the b array.

!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.

subroutine MatInv(n,a,m,b,determ)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   integer(4), intent(in) :: n          ! extent of the a matrix
   integer(4), intent(in) :: m          ! number of equations to be solved
   real(8), intent(inout) :: a(n,n)     ! in: matrix to be inverted
                                        ! out: inverted matrix
   real(8), intent(inout) :: b(n,m)     ! in: known vector(s)
                                        ! out: solution
   real(8), intent(out) :: determ       ! out: determinant
   integer(4), parameter :: nmax = 500
   integer(4) :: ipiv(nmax),indxr(nmax),indxc(nmax)
   integer(4) :: i, k, j, l, ll, irow, icol
   real(8) :: big, dum, pivinv

   determ = one
   ipiv(1:n) = 0

   do i=1,n

! ... search for pivot element

      big = zero
      do j = 1,n
         if (ipiv(j) /= 1) then
            do k = 1,n
               if (ipiv(k) == 0) then
                  if (abs(a(j,k)) >= big) then
                     big = abs(a(j,k))
                     irow = j
                     icol = k
                  endif
               else if (ipiv(k) > 1) then
                  call SubMatInv
               endif
            end do
         endif
      end do
      ipiv(icol) = ipiv(icol)+1

! ... interchange rows to put pivot element on diagonal

      if (irow /= icol) then
         determ = -determ
         do l = 1,n
            dum = a(irow,l)
            a(irow,l) = a(icol,l)
            a(icol,l) = dum
         end do
         do l = 1,m
            dum = b(irow,l)
            b(irow,l) = b(icol,l)
            b(icol,l) = dum
         end do
      endif
      indxr(i) = irow
      indxc(i) = icol
      if (a(icol,icol) == zero) call SubMatInv
      pivinv = one/a(icol,icol)
      determ = determ*a(icol,icol)

! ... divide pivot row by pivot element

      a(icol,icol) = one
      a(icol,1:n) = a(icol,1:n)*pivinv
      b(icol,1:m) = b(icol,1:m)*pivinv

! ... reduce non-pivot rows

      do ll = 1,n
         if (ll /= icol) then
           dum = a(ll,icol)
           a(ll,icol) = zero
           a(ll,1:n) = a(ll,1:n)-a(icol,1:n)*dum
           b(ll,1:m) = b(ll,1:m)-b(icol,1:m)*dum
         endif
      end do
   end do

! ... interchange columns

   do l = n,1,-1
      if (indxr(l) /= indxc(l)) then
         do k = 1,n
            dum = a(k,indxr(l))
            a(k,indxr(l)) = a(k,indxc(l))
            a(k,indxc(l)) = dum
         end do
      endif
   end do

contains

!........................................................................

subroutine SubMatInv
   call warn('SubMatInv', 'singular matrix', 6)
end subroutine SubMatInv

!........................................................................

end subroutine MatInv

!************************************************************************
!*                                                                      *
!*     Diag                                                             *
!*                                                                      *
!************************************************************************

! ... diagonalize a real matrix and calculate eigenvectors

subroutine Diag_old(n,a,eivr)
   implicit none
   real(8), parameter :: zero = 0.0d0, fourth = 0.25d0, one = 1.0d0, two = 2.0d0
   integer(4), intent(in)    :: n          ! extent of the matrices
   real(8),    intent(inout) :: a(n,n)     ! in: matrix to diagonalize
                                           ! out: diagonal contains the eigenv.
   real(8),    intent(inout) :: eivr(n,n)  ! out: eigenvectors stored in columns
   integer(4) :: i, i2, j, ii, jj, irow,jcol, jcol1, iflag
   real(8)    :: atop, aii, ajj, aij, avgf, c, d, dstop, s, t, u, thrsh

   do 10 j = 1,n
   do 20 i = 1,n
20 eivr(i,j) = zero
10 eivr(j,j) = one

   atop = 0.
   do 30 i = 1,n
   do 30 j = i,n
   if (atop-abs(a(i,j))) 31,30,30
31 atop = abs(a(i,j))
30 continue
   if (atop)2,2,3
2  return
3  avgf = (n*(n-1))*.55
   d = 0.0
   do 40 jj = 2,n
   do 40 ii = 2,jj
   s = a(ii-1,jj)/atop
40 d = s*s+d
   dstop = (1.e-06)*d
   thrsh = sqrt(d/avgf)*atop
4  iflag = 0
   do 50 jcol = 2,n
   jcol1 = jcol-1
   do 50 irow = 1,jcol1
   aij = a(irow,jcol)
   if (abs(aij)-thrsh) 50,50,51
51 aii = a(irow,irow)
   ajj = a(jcol,jcol)
   s = ajj-aii
   if (abs(aij)-1.e-09*abs(s)) 50,50,52
52 iflag = 1
   if (1.e-10*abs(aij)-abs(s)) 53,54,54
54 s = .707106781d0
   c = s
   go to 55
53 t = aij/s
   s = fourth/sqrt(fourth+t*t)
   c = sqrt(0.5d0+s)
   s = two*t*s/c
55 do 60 i = 1,irow
   t = a(i,irow)
   u = a(i,jcol)
   a(i,irow) = c*t-s*u
60 a(i,jcol) = s*t+c*u
   i2 = irow+2
   if (i2-jcol)5,5,6
5  continue
   do 70 i = i2,jcol
   t = a(i-1,jcol)
   u = a(irow,i-1)
   a(i-1,jcol) = s*u+c*t
70 a(irow,i-1) = c*u-s*t
6  a(jcol,jcol) = s*aij+c*ajj
   a(irow,irow) = c*a(irow,irow)-s*(c*aij-s*ajj)
   do 80 j = jcol,n
   t = a(irow,j)
   u = a(jcol,j)
   a(irow,j) = c*t-s*u
80 a(jcol,j) = s*t+c*u
   do 90 i = 1,n
   t = eivr(i,irow)
   eivr(i,irow) = c*t-eivr(i,jcol)*s
90 eivr(i,jcol) = s*t+eivr(i,jcol)*c
   s = aij/atop
   d = d-s*s
   if (d-dstop)7,8,8
7  d = 0.
   do 100 jj = 2,n
   do 100 ii = 2,jj
   s = a(ii-1,jj)/atop
100 d = s*s+d
   dstop = (1.e-06)*d
8  thrsh = sqrt(d/avgf)*atop
50 continue
   if (iflag)4,9,4
9  return
end subroutine Diag_old

!************************************************************************
!*                                                                      *
!*     Diag                                                             *
!*                                                                      *
!************************************************************************

! ... diagonalize a real matrix and calculate eigenvectors

!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.
!     -> subroutine "jacobi"

subroutine Diag(n,a,d,v,nrot)
   implicit none
   real(8), parameter :: zero = 0.0d0, fourth = 0.25d0, one = 1.0d0, two = 2.0d0
   integer(4), parameter :: mnrot = 50     ! maximum number sweeps
   integer(4), intent(in)    :: n          ! extent of the matrices
   real(8),    intent(inout) :: a(n,n)     ! in: matrix to diagonalize
                                           ! out: upper right triangle destroyed
   real(8),    intent(out) :: d(n)         ! out: eigenvalues
   real(8),    intent(out) :: v(n,n)       ! out: eigenvectors stored in columns
   integer(4), intent(out) :: nrot         ! number of sweeps (n(n-1)/2 Jacobi rotations)
   real(8) :: b(n),z(n)
   real(8) :: sm, t, tresh, h, g, theta, c, s, tau
   integer(4) :: ip, iq, i, j

   v = zero
   do ip = 1,n
      v(ip,ip) = one
   end do
   do ip = 1,n
      b(ip) = a(ip,ip)
      d(ip) = b(ip)
      z(ip) = zero
   end do

   nrot = 0
   do i = 1,mnrot

! ... check for convergence

      sm = zero
      do  ip = 1,n-1
         do iq = ip+1,n
            sm = sm+abs(a(ip,iq))
         end do
      end do
      if (sm == zero) exit

! ... make n(n-1)/2 Jacobi rotations

      tresh = zero
      if (i < 4) tresh = 0.2d0*sm/n**2
      do ip = 1,n-1
         do iq = ip+1,n
            g = 100.0d0*abs(a(ip,iq))
            if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) .and. (abs(d(iq))+g == abs(d(iq)))) then
               a(ip,iq) = zero
            else if (abs(a(ip,iq)) > tresh) then
               h = d(iq)-d(ip)
               if (abs(h)+g == abs(h)) then
                  t = a(ip,iq)/h
               else
                  theta = 0.5d0*h/a(ip,iq)
                  t = one/(abs(theta)+sqrt(one+theta**2))
                  if (theta < zero) t=-t
               endif
               c = one/sqrt(1+t**2)
               s = t*c
               tau = s/(one+c)
               h = t*a(ip,iq)
               z(ip) = z(ip)-h
               z(iq) = z(iq)+h
               d(ip) = d(ip)-h
               d(iq) = d(iq)+h
               a(ip,iq) = zero
               do j = 1,ip-1
                  g = a(j,ip)
                  h = a(j,iq)
                  a(j,ip) = g-s*(h+g*tau)
                  a(j,iq) = h+s*(g-h*tau)
               end do
               do j = ip+1,iq-1
                  g = a(ip,j)
                  h = a(j,iq)
                  a(ip,j) = g-s*(h+g*tau)
                  a(j,iq) = h+s*(g-h*tau)
               end do
               do j = iq+1,n
                  g = a(ip,j)
                  h = a(iq,j)
                  a(ip,j) = g-s*(h+g*tau)
                  a(iq,j) = h+s*(g-h*tau)
               end do
               do j = 1,n
                  g = v(j,ip)
                  h = v(j,iq)
                  v(j,ip) = g-s*(h+g*tau)
                  v(j,iq) = h+s*(g-h*tau)
               end do
               nrot = nrot+1
            endif
         end do
      end do

      do ip = 1,n
         b(ip) = b(ip)+z(ip)
         d(ip) = b(ip)
         z(ip) = zero
      end do

   end do
   if (nrot > mnrot) call stop('Diag', 'no convergence', 6)

end subroutine Diag

!************************************************************************
!*                                                                      *
!*     Eigensort                                                        *
!*                                                                      *
!************************************************************************

! ... Given the eigenvalues d and eigenvectors v as output from Diag this
! ... routine sorts the eigenvalues into descending order, and rearranges
! ... the columns of v correspondingly. The method is straight insertion.

!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.
!     --> subroutine "eigsrt"

subroutine Eigensort(d,v,n)

   implicit none

   integer(4)  :: n
   real(8)     :: d(n), v(n,n)
   integer(4)  :: i, j, k
   real(8)     :: p

   do i = 1, n-1
      k=i
      p=d(i)
      do j = i+1, n
         if (d(j).ge.p) then
            k=j
            p=d(j)
         end if
      end do
      if (k.ne.i) then
         d(k)=d(i)
         d(i)=p
         do j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
         end do
      end if
   end do

   return

end subroutine

!************************************************************************
!*                                                                      *
!*     Smooth                                                           *
!*                                                                      *
!************************************************************************

! ... smooth by using spline functions. argonne routine e350s
!     see C. H. Reinsch, Num Math 10, 177 (1967), 16, 451 (1971)
!     x(i) are strictly increasing or decreasing (see comment below)

subroutine Smooth(n,x,y,dy,s,a,b,c,d)
   implicit none
   real(8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0, three = 3.0d0
   integer(4), intent(in)  :: n                      ! extent of arrays
   real(8),    intent(in)  :: x(n), y(n), dy(n)      ! function to smooth
   real(8),    intent(in)  :: s                      ! smooth parameter
   real(8),    intent(out) :: a(n), b(n), c(n), d(n) ! spline coeff
   integer(4) :: i, i1, nm1
   real(8)    :: r(n), r1(n), r2(n), t(n), t1(n), u(n), v(n)
   real(8)    :: e, f, f2, g, h, p

   nm1 = n-1

   h = x(2)-x(1)
   f = (y(2)-y(1))/h
   do i = 2,nm1
      g = h
      h = x(i+1)-x(i)
      e = f
      f = (y(i+1)-y(i))/h
      a(i) = f-e
      t(i) = two*(g+h)/three
      t1(i) = h/three
      r2(i) = dy(i-1)/g
      r(i) = dy(i+1)/h
      r1(i) = -dy(i)/g-dy(i)/h
   end do
   r(1) = zero
   r1(n) = zero
   r2(n) = zero

   do i = 2,nm1
      b(i) = r(i)**2+r1(i)**2+r2(i)**2
      c(i) = r(i)*r1(i+1)+r1(i)*r2(i+1)
      if (i /= nm1) d(i) = r(i)*r2(i+2)
   end do
   d(nm1) = zero

   p = zero
   f2 = -s

! ... next iteration

300 e = zero
   u(1) = zero
   u(n) = zero
   do i = 2,nm1
      if (i == 2) go to 350
      r2(i-2) = g*r(i-2)
      e = g*r2(i-2)
350   r1(i-1) = f*r(i-1)
      r(i) = one/(p*b(i)+t(i)-f*r1(i-1)-e)
      u(i) = a(i)-r1(i-1)*u(i-1)
      if (i /= 2) u(i) = u(i)-r2(i-2)*u(i-2)
      f = p*c(i)+t1(i)-h*r1(i-1)
      g = h
      h = d(i)*p
   end do

   do i1 = 2,nm1
      i = nm1+2-i1
      u(i) = r(i)*u(i)-r1(i)*u(i+1)
      if (i /= nm1) u(i) = u(i)-r2(i)*u(i+2)
   end do

   e = zero
   h = zero
   do i = 1,nm1
      g = h
      h = (u(i+1)-u(i))/(x(i+1)-x(i))
      v(i) = (h-g)*dy(i)*dy(i)
      e = e+v(i)*(h-g)
   end do
   v(n) = -h*dy(n)*dy(n)
   e = e-v(n)*h
   g = f2
   f2 = e*p*p
   if (f2 >= s.or.f2 <= g) go to 800

   f = zero
   h = (v(2)-v(1))/(x(2)-x(1))
   do i = 2,nm1
      g = h
      h = (v(i+1)-v(i))/(x(i+1)-x(i))
      g = h-g-r1(i-1)*r(i-1)
      if (i /= 2) g = g-r2(i-2)*r(i-2)
      f = f+g*r(i)*g
      r(i) = g
   end do
   h = e-p*f
   if (h <= zero) go to 800

! ... use negative branch of square root if x(i) are strictly decreasing

   p = p+(s-f2)/((p+sqrt(s/e))*h)
   go to 300

! ... final step after convergence

800 continue
    do i = 1,n
      a(i) = y(i)-p*v(i)
      c(i) = u(i)
    end do
    do i = 1,nm1
       h = x(i+1)-x(i)
       d(i) = (c(i+1)-c(i))/(three*h)
       b(i) = (a(i+1)-a(i))/h-(h*d(i)+c(i))*h
    end do

end subroutine Smooth

!************************************************************************
!*                                                                      *
!*     CsEval                                                           *
!*                                                                      *
!************************************************************************

! ... return the value and the first two derivatives of a spline func.

!       ideriv = 0:  a(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!       ideriv = 1:  b(i)+2*c(i)*(u-x(i))+3*d(i)*(u-x(i))**2
!       ideriv = 2:  2*c(i)+6*d(i)*(u-x(i))
!       ideriv =-1:  integral of the cubic spline funtion

!       where  x(i) < u < x(i+1), using horner's rule

!     if u < x(1) then  i = 1 is used.
!     if u >= x(n) then  i = n is used.

real(8) function CsEval(n,u,ideriv,x,a,b,c,d)
   implicit none
   real(8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0, three = 3.0d0, four = 4.0d0
   integer(4), intent(in) :: n                      ! extent of arrays
   real(8), intent(in) :: u                         ! absicssa
   integer(4), intent(in) :: ideriv                 ! = 0, calculate f(u)
                                                    ! = 1, calculate f'(u)
                                                    ! = 2, calculate f''(u)
                                                    ! =-1, integrate from x(1) to u
   real(8),    intent(in) :: x(n)                   ! array of x-values
   real(8),    intent(in) :: a(n), b(n), c(n), d(n) ! spline coeff
   integer(4), save :: i = 1
   integer(4) :: j, k
   real(8)    :: dx

   if (ideriv < -1.or.ideriv > 2) call stop('CsEval','check ideriv',6)

   if (i >= n) i = 1
   if (u < x(i).or.u > x(i+1)) then

! ... binary search if outside previous interval

      i = 1
      j = n+1
20    k = (i+j)/2
      if (u < x(k)) j = k
      if (u >= x(k)) i = k
      if (j > i+1) go to 20

   end if

! ... evaluate the spline

   if (ideriv == -1) then
      CseVal = zero
      if (u <= x(1)) return
      dx = u-x(i)
      CseVal = dx*(a(i)+dx*(b(i)/two+dx*(c(i)/three+dx*d(i)/four)))
      if (i > 1) then
         do j = 1,i-1
            dx = x(j+1)-x(j)
            CseVal = CseVal+dx*(a(j)+dx*(b(j)/two+dx*(c(j)/three+dx*d(j)/four)))
         end do
      end if
   else
      dx = u-x(i)
      if (ideriv == 0) then
         CseVal = a(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      else if (ideriv == 1) then
         CseVal = b(i)+dx*(two*c(i)+three*dx*d(i))
      else if (ideriv == 2) then
         CseVal = two*c(i)+6.0d0*dx*d(i)
      end if
   end if
end function CsEval

!************************************************************************
!*                                                                      *
!*     VInter                                                           *
!*                                                                      *
!************************************************************************

! ... interpolate in a 1d array

real(8) function VInter(n,x,y,xin)
   implicit none
   integer(4), intent(in) :: n           ! extent of arrays
   real(8),    intent(in) :: x(n), y(n)  ! function to use
   real(8),    intent(in) :: xin         ! x value for which the y value is desired
   integer(4) :: i, j, k
   real(8) :: LinInter

   if (x(1) < x(n)) then
      if (xin < x(1)) then
         VInter = y(1)
      else if (xin > x(n)) then
         VInter = y(n)
      else
         i = 1
         j = n+1
   20    k = (i+j)/2
         if (xin < x(k)) j = k
         if (xin >= x(k)) i = k
         if (j > i+1) goto 20
         VInter = LinInter(x(i),x(i+1),xin,y(i),y(i+1))
      end if
   else
      if (xin < x(n)) then
         VInter = y(n)
       else if (xin > x(1)) then
         VInter = y(1)
      else
         i = 1
         j = n+1
    30   k = (i+j)/2
         if (xin < x(k)) i = k
         if (xin >= x(k)) j = k
         if (j > i+1) goto 30
         VInter = LinInter(x(i),x(i+1),xin,y(i),y(i+1))
     end if
   end if

end function VInter

!************************************************************************
!*                                                                      *
!*     LinInter                                                         *
!*                                                                      *
!************************************************************************

! ... linear interpolation between (x0,y0) and (x1,y1)

real(8) function LinInter(x0,x1,xk,y0,y1)
   implicit none
   real(8), intent(in) :: x0 ! lower x-value
   real(8), intent(in) :: x1 ! upper x-value
   real(8), intent(in) :: xk ! x-value
   real(8), intent(in) :: y0 ! lower y-value
   real(8), intent(in) :: y1 ! upper y-value
   if ((x1-xk)*(xk-x0) < 0.0d0) call Warn('LinInter','(x1-xk)*(xk-x0) < 0.0d0',6)
   LinInter = y0 + (xk-x0) / (x1-x0) * (y1-y0)
end function LinInter

!************************************************************************
!*                                                                      *
!*     CardinalBSpline                                                  *
!*                                                                      *
!************************************************************************

! ... calculate values and derivatives of the n:th order cardinal B-spline
!     ref. Schoenberg, Cardinal Spline Interpolation
!     (Society for Industrial and Applied Mathematics, Philadelphia, PA, 1973)

module CardinalBSplineModule

contains

subroutine CardinalBSpline(order, pos, kmax, values, derivs, derivs2)
    implicit none
    integer(4), intent(in)  :: order                     ! order of the cardinal B-spline
    real(8),    intent(in)  :: pos                       ! position
    integer(4), intent(out) :: kmax                      ! floor(pos)
    real(8),  intent(out)   :: values(0:order-1)         ! values at M_n(pos-kmin)...M_n(pos-(kmin+n-1))
    real(8), intent(out), optional :: derivs(0:order-1)  ! values of first derivative
    real(8), intent(out), optional :: derivs2(0:order-1) ! values of second derivative

    integer(4) :: v, n, i, kmin
    real(8) :: x, a
    real(8) :: lastvals(0:order-1), lastvals2(0:order-1)

    v = floor(pos)
    x = pos-v                   !  pos = v + x: int(v), 0<=x<1
    kmax = v
    kmin = kmax-order+1
    lastvals = 0.0d0
    lastvals2 = 0.0d0

    values(0) = x
    if (order > 1)  then        ! from version 1.6.0
       values(1) = (1.0d0-x)
       do n = 3, order
           lastvals = values
           if ((n == order-1) .and. present(derivs2)) lastvals2 = values
           values(0) = x*lastvals(0)/(n-1)
           do i = 1, n-2
               a = i + x
               values(i) = (a*lastvals(i)+(n-a)*lastvals(i-1))/(n-1)
           end do
           values(n-1) = (1-x)*lastvals(n-2)/(n-1)
       end do

       n = order
       if (present(derivs)) then
           derivs(0) = lastvals(0)
           do i = 1, n-2
               derivs(i) = lastvals(i)-lastvals(i-1)
           end do
           derivs(n-1) = -lastvals(n-2)
       end if

       if (present(derivs2)) then
           derivs2(0) = lastvals2(0)
           derivs2(1) = lastvals2(1)-2*lastvals2(0)
           do i = 2, n-3
               derivs2(i) = lastvals2(i)-2*lastvals2(i-1)+lastvals2(i-2)
           end do
           if (n > 3) derivs2(n-2) = -2*lastvals2(n-3)+lastvals2(n-4)
           if (n > 2) derivs2(n-1) = lastvals2(n-3)
       end if
    end if

end subroutine CardinalBSpline

!........................................................................

end module CardinalBSplineModule

!************************************************************************
!*                                                                      *
!*     HeapSort                                                         *
!*                                                                      *
!************************************************************************

! ... sort an array into ascending order, Heapsort algorithm
!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.

subroutine HeapSort(n,vec)
   implicit none
   integer(4), intent(in) :: n         ! extent of the array
   real(8),    intent(inout) :: vec(n) ! in: unsorted array
                                       ! out: sorted array
   integer(4) :: l, i, j, ir
   real(8) :: vect

   if (n < 0) call stop('HeapSort','n < 0',6)
   if (n == 1) return
   l = n/2+1
   ir = n
10 continue
      if (l>1) then
         l = l-1
         vect = vec(l)
      else
         vect = vec(ir)
         vec(ir) = vec(1)
         ir = ir-1
         if (ir == 1) then
            vec(1) = vect
            return
         end if
      end if
      i = l
      j = l+l
20    if (j <=  ir) then
         if (j < ir) then
            if (vec(j) < vec(j+1)) j = j+1
         end if
         if (vect < vec(j)) then
            vec(i) = vec(j)
            i = j
            j = j+j
         else
            j = ir+1
         end if
      goto 20
      end if
      vec(i) = vect
   goto 10
end subroutine HeapSort

!************************************************************************
!*                                                                      *
!*     HeapSortIndex                                                    *
!*                                                                      *
!************************************************************************

! ... generate an index array such that vec(index(j)) is ascending, Heapsort algorithm
!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.

subroutine HeapSortIndex(n,vec,index)
   implicit none
   integer(4), intent(in)  :: n        ! extent of the array
   real(8),    intent(in)  :: vec(n)   ! unsorted array
   integer(4), intent(out) :: index(n) ! index array

   integer(4) :: l, i, j, ir,indext
   real(8) :: vect

   if (n < 0) call stop('HeapSortIndex','n < 0',6)

   index(1:n) = [ (i,i = 1,n) ]
   if (n == 1) return
   l = n/2+1
   ir = n
10 continue
      if (l>1) then
         l = l-1
         indext = index(l)
         vect = vec(indext)
      else
         indext = index(ir)
         vect = vec(indext)
         index(ir) = index(1)
         ir = ir-1
         if (ir == 1) then
         index(1) = indext
            return
         end if
      end if
      i = l
      j = l+l
   20 if (j <= ir) then
         if (j < ir) then
            if (vec(index(j)) < vec(index(j+1))) j = j+1
         end if
         if (vect < vec(index(j))) then
            index(i) = index(j)
            i = j
            j = j+j
         else
            j = ir+1
         end if
      goto 20
      end if
      index(i) = indext
   goto 10
end subroutine HeapSortIndex

!************************************************************************
!*                                                                      *
!*     MakeCluster                                                      *
!*                                                                      *
!************************************************************************

! ... make clusters

subroutine MakeCluster(nobj, npair, n1, n2, icliobj)

   implicit none

   integer(4), intent(in)  :: nobj                 ! number of objects
   integer(4), intent(in)  :: npair                ! number of object pairs
   integer(4), intent(in)  :: n1(*), n2(*)         ! object number of a pair of objects
   integer(4), intent(out) :: icliobj(*)           ! pointer: object -> its cluster

   if (nobj > npair) then
      call MakeCluster1(nobj, npair, n1, n2, icliobj)
   else
      call MakeCluster2(nobj, npair, n1, n2, icliobj)
   end if

end subroutine MakeCluster

!************************************************************************
!*                                                                      *
!*     MakeCluster1                                                     *
!*                                                                      *
!************************************************************************

! ... make clusters
!     nobj >> npair: order(nobj*npair) operations
!     nobj << npair: order(npair) operations

subroutine MakeCluster1(nobj, npair, n1, n2, icliobj)

   implicit none

   integer(4), intent(in)  :: nobj                 ! number of objects
   integer(4), intent(in)  :: npair                ! number of object pairs
   integer(4), intent(in)  :: n1(*), n2(*)         ! object number of a pair of objects
   integer(4), intent(out) :: icliobj(*)           ! pointer: object -> its cluster

   integer(4) :: ipair, icl, icl1, icl2

   icliobj(1:nobj) = [ (icl,icl = 1,nobj) ]        ! initizate iclobj

   do ipair = 1, npair                             ! start make clusters
      icl1 = icliobj(n1(ipair))
      icl2 = icliobj(n2(ipair))
      if (icl1/= icl2) where(icliobj(1:nobj) == icl2) icliobj(1:nobj) = icl1
!     write(*,'(i5,5x,20i5)') ipair, icliobj(1:nobj)
   end do

end subroutine MakeCluster1

!************************************************************************
!*                                                                      *
!*     MakeCluster2                                                     *
!*                                                                      *
!************************************************************************

! ... make clusters
!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.
!     npair >> nobj: number of operations increases fast with increasing npair

subroutine MakeCluster2(nobj, npair, n1, n2, icliobj)

   implicit none

   integer(4), intent(in)  :: nobj                 ! number of objects
   integer(4), intent(in)  :: npair                ! number of object pairs
   integer(4), intent(in)  :: n1(*), n2(*)         ! object number of a pair of objects
   integer(4), intent(out) :: icliobj(*)           ! pointer: object -> its cluster

   integer(4) :: icl, ipair, iobj, iobj1, iobj2

   icliobj(1:nobj) = [ (icl,icl = 1,nobj) ]        ! initizate iclobj

   do ipair = 1, npair                             ! start make clusters
      iobj1 = n1(ipair)
      iobj2 = n2(ipair)
      call MakeClusterSub(iobj1)                   ! get top ancestor of iobj1
      call MakeClusterSub(iobj2)                   ! get top ancestor of iobj2
      if (iobj1 /= iobj2) icliobj(iobj2) = iobj1   ! if different top ancestor, merge the two families
!     write(*,'(i5,5x,10i5)') ipair, icliobj(1:nobj)
   end do
   do iobj = 1, nobj                               ! fill up
      do
         if (icliobj(iobj) == icliobj(icliobj(iobj))) exit
         icliobj(iobj) = icliobj(icliobj(iobj))
      end do
   end do

contains

! ........................................................................

subroutine MakeClusterSub(iobj)
   integer(4), intent(inout) :: iobj
   do
      if (iobj == icliobj(iobj)) exit
      iobj = icliobj(iobj)
   end do
end subroutine MakeClusterSub

! ........................................................................

end subroutine MakeCluster2

!************************************************************************
!*                                                                      *
!*     CalcClusterMember                                                *
!*                                                                      *
!************************************************************************

! ... get members of the cluster to which object iobj belongs to

subroutine CalcClusterMember(nobj, icliobj, iobj, nmem, imemlis)

   implicit none

   integer(4), intent(in)  :: nobj                ! number of objects
   integer(4), intent(in)  :: icliobj(*)          ! pointer: object -> its cluster
   integer(4), intent(in)  :: iobj                ! object
   integer(4), intent(out) :: nmem                ! number of members
   integer(4), intent(out) :: imemlis(*)          ! members of the cluster to which object iobj belongs to

   integer(4) :: jobj

   if((iobj < 1) .or. (iobj > nobj)) call Stop('CalcClusterMember', 'objt out of range', 6)

   nmem = 0
   do jobj = 1, nobj                              ! loop over objects
      if (icliobj(jobj) == icliobj(iobj)) then    ! check if object jobj belongs to the cluster of object iobj
         nmem = nmem+1
         imemlis(nmem) = jobj
      end if
   end do

end subroutine CalcClusterMember

! ... quaternions:
!     ref: Goldstein 'Classical mechanics' p. 147-155 and Evans & s. Murad, Mol. Phys., 32, 327 (1977).
!     notation:
!               here              Evans & Murad      Goldstein
!     matrix    ori(-1)           a  equ (4)         a (equ 4-46,4-67)
!     quat.     qua(0)            chi                e(0)
!               qua(1)            eta                e(1)
!               qua(2)            xsi               -e(2)
!               qua(3)            zet                e(3)

!************************************************************************
!*                                                                      *
!*     CarToSph                                                         *
!*                                                                      *
!************************************************************************

! ... transform a vector from a cartensian to a spherical polar coordinate system

subroutine CarToSph(txoptin,x,y,z,r,theta,phi)
   implicit none
   real(8), parameter :: radtodeg = 57.29577951d0
   real(8), parameter :: pihalf = 1.570796327d0
   real(8), parameter :: zero = 0.0d0, one = 1.0d0, small = 1.0d-6
   character(*), intent(in)  :: txoptin              ! select 'rad' or 'deg'
   real(8),      intent(in)     :: x, y, z            ! cartesian coordinate
   real(8),      intent(out)    :: r, theta, phi      ! spherical polar coordinate
   real(8) :: norm, xn, yn, zn
   character(LEN=len(txoptin))  :: txopt              ! select 'rad' or 'deg'

   txopt = txoptin
   call LowerCase(txopt)

   r = sqrt(x**2+y**2+z**2)
   if (r == zero) call Stop('CarToSph','r == zero', 6)
   norm = one/r
   xn = max(-one,min(x*norm,one))
   yn = max(-one,min(y*norm,one))
   zn = max(-one,min(z*norm,one))
   theta = acos(zn)
   if (abs(xn) < small) then
      phi =        + sign(pihalf,yn)
   else if (abs(yn) < small) then
      phi = pihalf + sign(pihalf,-xn)
   else
      phi = atan(yn/xn) + pihalf + sign(pihalf,-xn)
      if (phi < 0) phi = phi + 4.0d0*pihalf
   end if
   if (txopt == 'deg') then
      theta = radtodeg*theta
      phi   = radtodeg*phi
   end if

end subroutine CarToSph

!************************************************************************
!*                                                                      *
!*     SphToCar                                                         *
!*                                                                      *
!************************************************************************

! ... transform a vector from a spherical polar to a cartesian coordinate system

subroutine SphToCar(txopt,r,theta,phi,x,y,z)
   implicit none
   real(8), parameter :: degtorad = 0.01745329252d0
   character(*), intent(inout)  :: txopt              ! select 'rad' or 'deg'
   real(8),      intent(in)     :: r, theta, phi      ! spherical polar coordinate
   real(8),      intent(out)    :: x, y, z            ! cartesian coordinate

   call LowerCase(txopt)

   if (txopt == 'deg') then
      x = r*sin((degtorad*theta))*cos((degtorad*phi))
      y = r*sin((degtorad*theta))*sin((degtorad*phi))
      z = r*cos((degtorad*theta))
   else
      x = r*sin(theta)*cos(phi)
      y = r*sin(theta)*sin(phi)
      z = r*cos(theta)
   end if

end subroutine SphToCar

!************************************************************************
!*                                                                      *
!*     CarToStd1                                                        *
!*                                                                      *
!************************************************************************

! ... transform a vector from a cartensian coordinate system to standard form

subroutine CarToStd1(x,y,z,q)
   implicit none
   real(8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0
   real(8), parameter :: sqrt2i = one/dsqrt(two)
   real(8),     intent(in)  :: x, y, z            ! vector in cartesian coordinate system
   complex(8),  intent(out) :: q(-1:1)            ! vector in standard form

   q( 1) =  -cmplx(x, +y)*sqrt2i
   q( 0) =  cmplx(z, zero)
   q(-1) =  cmplx(x, -y)*sqrt2i

end subroutine CarToStd1

!************************************************************************
!*                                                                      *
!*     Std1ToCar                                                        *
!*                                                                      *
!************************************************************************

! ... transform a vector from standard form to a cartensian coordinate system

subroutine Std1ToCar(q,x,y,z)
   implicit none
   real(8), parameter :: zero = 0.0d0, one = 1.0d0, two = 2.0d0
   real(8), parameter :: sqrt2i = one/dsqrt(two)
   complex(8),  intent(in)  :: q(-1:1)            ! vector in standard form
   real(8),     intent(out) :: x, y, z            ! vector in cartesian coordinate system

   x = -real(q(1)-q(-1))*sqrt2i
   y = -imag(q(1)+q(-1))*sqrt2i
   z =  real(q(0))

end subroutine Std1ToCar

!************************************************************************
!*                                                                      *
!*     OriToEuler                                                       *
!*                                                                      *
!************************************************************************

! ... transforming orientation matrix ori to Euler angles a, b, c
!     (zxz or zyz convention)

subroutine OriToEuler(txopt, ori, a, b, c)

   implicit none
   real(8), parameter :: Zero = 0.0d0, one = 1.0d0
   real(8), parameter :: Pi2 = 6.28318530718d0
   character(1), intent(in)  :: txopt
   real(8),      intent(in)  :: ori(3,3)            ! orientation matrix
   real(8),      intent(out) :: a, b, c             ! Euler angles
   real(8) :: ca, sa, sb, cc, sc

   b = acos(max(-one,min(ori(3,3),one)))
   sb = sin(b)

   if (txopt == 'x') then
      if (abs(sb) > 1e-8) then
         ca =-ori(2,3)/sb
         sa = ori(1,3)/sb
         a = acos(max(-one,min(ca,one)))
         if (sa < Zero) a = Pi2-a
         cc = ori(3,2)/sb
         sc = ori(3,1)/sb
         c = acos(max(-one,min(cc,one)))
         if (sc < Zero) c = Pi2-c
      else
         a = acos(max(-one,min(ori(1,1),one)))
         if (ori(2,1) < Zero) a =-a
         c = Zero
      end if
   else if (txopt == 'y') then
      if (abs(sb) > 1e-8) then
         sa = ori(2,3)/sb
         ca = ori(1,3)/sb
         a = acos(max(-one,min(ca,one)))
         if (sa < Zero) a = Pi2-a
         sc = ori(3,2)/sb
         cc =-ori(3,1)/sb
         c = acos(max(-one,min(cc,one)))
         if (sc < Zero) c = Pi2-c
      else
         a = acos(max(-one,min(ori(2,2),one)))
         if (ori(2,1) < Zero) a =-a
         c = Zero
      end if
   else
      write(*,*) txopt
      call Stop('OriToEuler', 'check txopt', 6)
   end if

end subroutine OriToEuler

!************************************************************************
!*                                                                      *
!*     EulerToQua                                                       *
!*                                                                      *
!************************************************************************

! ... transform Euler angles a, b, c to quaternions

subroutine EulerToQua(txopt, a, b, c, qua)

   implicit none
   real(8), parameter :: half = 0.5d0, one = 1.0d0
   character(1), intent(in)  :: txopt
   real(8),      intent(in)  :: a, b, c               ! orientation matrix
   real(8),      intent(out) :: qua(0:3)              ! Euler angles
   real(8) :: cb2, sb2, norm

   cb2 = cos(half*b)
   sb2 = sin(half*b)

   if (txopt == 'x') then
      qua(0) = cb2 * cos(half*(c+a))
      qua(1) = sb2 * cos(half*(c-a))
      qua(2) = sb2 * sin(half*(c-a))
      qua(3) = cb2 * sin(half*(c+a))
   else if (txopt == 'y') then
      qua(0) = cb2 * cos(half*(c+a))
      qua(1) = sb2 * sin(half*(c-a))
      qua(2) =-sb2 * cos(half*(c-a))
      qua(3) = cb2 * sin(half*(c+a))
   else
      call Stop('EulerToQua', 'check txopt', 6)
   end if

   norm = one/sqrt(qua(0)**2+qua(1)**2+qua(2)**2+qua(3)**2)
   qua(0) = qua(0)*norm
   qua(1) = qua(1)*norm
   qua(2) = qua(2)*norm
   qua(3) = qua(3)*norm

end subroutine EulerToQua

!************************************************************************
!*                                                                      *
!*     OriToQua                                                         *
!*                                                                      *
!************************************************************************

! ... transform orientation matrix to quaternions

!     note, the routine selects the branch where the largest quaternion is positive

subroutine OriToQua(np,iplow, ipupp, ori, qua)
   implicit none
   real(8), parameter :: fourth = 0.25d0, one = 1.0d0, two = 2.0d0
   integer(4), intent(in)  :: np
   integer(4), intent(in)  :: iplow
   integer(4), intent(in)  :: ipupp
   real(8),    intent(in)  :: ori(3,3,np)
   real(8),    intent(out) :: qua(0:3,np)
   integer(4) :: ip
   real(8)    :: trace, qua00, qua11, qua22, qua33, quafac

   if (iplow < lbound(qua,2)) call stop('OriToQua', 'iplow < lbound(qua,2)', 6)
   if (ipupp > ubound(qua,2)) call stop('QriToQua', 'ipupp > ubound(qua,2)', 6)
   do ip = iplow, ipupp
      trace = ori(1,1,ip)+ori(2,2,ip)+ori(3,3,ip)
      qua00 = fourth*(one+trace)
      qua11 = fourth*(one+two*ori(1,1,ip)-trace)
      qua22 = fourth*(one+two*ori(2,2,ip)-trace)
      qua33 = fourth*(one+two*ori(3,3,ip)-trace)
      if (qua00 >= qua11 .and. qua00 >= qua22 .and. qua00 >= qua33) then
         qua(0,ip) = sqrt(qua00)
         quafac = fourth/qua(0,ip)
         qua(1,ip) = quafac*(+ori(3,2,ip)-ori(2,3,ip))
         qua(2,ip) = quafac*(+ori(3,1,ip)-ori(1,3,ip))
         qua(3,ip) = quafac*(+ori(2,1,ip)-ori(1,2,ip))
      else if (qua11 >= qua22 .and. qua11 >= qua33 .and. qua11 >= qua00) then
         qua(1,ip) = sqrt(qua11)
         quafac = fourth/qua(1,ip)
         qua(2,ip) = quafac*(-ori(2,1,ip)-ori(1,2,ip))
         qua(3,ip) = quafac*(+ori(3,1,ip)+ori(1,3,ip))
         qua(0,ip) = quafac*(+ori(3,2,ip)-ori(2,3,ip))
      else if (qua22 >= qua33 .and. qua22 >= qua00 .and. qua22 >= qua11) then
         qua(2,ip) = sqrt(qua22)
         quafac = fourth/qua(2,ip)
         qua(3,ip) = quafac*(-ori(3,2,ip)-ori(2,3,ip))
         qua(0,ip) = quafac*(+ori(3,1,ip)-ori(1,3,ip))
         qua(1,ip) = quafac*(-ori(2,1,ip)-ori(1,2,ip))
      else
         qua(3,ip) = sqrt(qua33)
         quafac = fourth/qua(3,ip)
         qua(0,ip) = quafac*(+ori(2,1,ip)-ori(1,2,ip))
         qua(1,ip) = quafac*(+ori(3,1,ip)+ori(1,3,ip))
         qua(2,ip) = quafac*(-ori(3,2,ip)-ori(2,3,ip))
      end if
   end do
   call QuaNorm(np, iplow, ipupp, qua)
end subroutine OriToQua

!************************************************************************
!*                                                                      *
!*     QuaToOri                                                         *
!*                                                                      *
!************************************************************************

! ... transform quaternions qua to orientation matrix

subroutine QuaToOri(np, iplow, ipupp, qua, ori)
   implicit none
   real(8), parameter :: two = 2.0d0
   integer(4), intent(in)  :: np
   integer(4), intent(in)  :: iplow
   integer(4), intent(in)  :: ipupp
   real(8),    intent(in)  :: qua(0:3,np)
   real(8),    intent(out) :: ori(3,3,np)
   integer(4) :: ip
   real(8) :: qua00, qua11, qua22, qua33, qua01, qua02, qua03, qua12, qua13, qua23

   if (iplow < lbound(qua,2)) call stop('QuaToOri', 'iplow < lbound(qua,2)', 6)
   if (ipupp > ubound(qua,2)) call stop('QuaToQri', 'ipupp > ubound(qua,2)', 6)
   do ip = iplow, ipupp
      qua00 = qua(0,ip)*qua(0,ip)
      qua11 = qua(1,ip)*qua(1,ip)
      qua22 = qua(2,ip)*qua(2,ip)
      qua33 = qua(3,ip)*qua(3,ip)
      qua01 = qua(0,ip)*qua(1,ip)
      qua02 = qua(0,ip)*qua(2,ip)
      qua03 = qua(0,ip)*qua(3,ip)
      qua12 = qua(1,ip)*qua(2,ip)
      qua13 = qua(1,ip)*qua(3,ip)
      qua23 = qua(2,ip)*qua(3,ip)
      ori(1,1,ip) = -qua22 + qua11 - qua33 + qua00
      ori(2,1,ip) =  two*( qua03 - qua12 )
      ori(3,1,ip) =  two*( qua13 + qua02 )
      ori(1,2,ip) = -two*( qua12 + qua03 )
      ori(2,2,ip) =  qua22 - qua11 - qua33 + qua00
      ori(3,2,ip) =  two*( qua01 - qua23 )
      ori(1,3,ip) =  two*( qua13 - qua02 )
      ori(2,3,ip) = -two*( qua23 + qua01 )
      ori(3,3,ip) = -qua22 - qua11 + qua33 + qua00
   end do
end subroutine QuaToOri

!************************************************************************
!*                                                                      *
!*     AxisAngToOri                                                     *
!*                                                                      *
!************************************************************************

! ... transform rotation axis and rotation angle to a rotation matrix

! see Goldstein "classical mechanics" p 165-166.

! actually the rotation is made around -rotaxis

subroutine AxisAngToOri(rotaxis, alpha, ori)

   implicit none
   real(8), parameter :: half = 0.5d0, one = 1.0d0, two = 2.0d0
   real(8), intent(in)  :: rotaxis(3)                   ! direction of rotation axis
   real(8), intent(in)  :: alpha                        ! rotation angle
   real(8), intent(out) :: ori(3,3)                     ! rotation matrix
   real(8) :: norm, ca, sa, qua0, qua1, qua2, qua3

   norm = one/sqrt(rotaxis(1)**2+rotaxis(2)**2+rotaxis(3)**2)
   ca = cos((half*alpha))
   sa = sin((half*alpha))

   qua0 = ca
   qua1 = rotaxis(1)*norm*sa
   qua2 = rotaxis(2)*norm*sa
   qua3 = rotaxis(3)*norm*sa

   if (((qua0**2)+(qua1**2)+(qua2**2)+(qua3**2))-one > 1.0d-12) call Warn('AxisAngToOri', 'qua is not normalized', 6)

   ori(1,1) = (qua0**2)+(qua1**2)-(qua2**2)-(qua3**2)
   ori(1,2) = two*((qua1*qua2)+(qua0*qua3))
   ori(1,3) = two*((qua1*qua3)-(qua0*qua2))
   ori(2,1) = two*((qua1*qua2)-(qua0*qua3))
   ori(2,2) = (qua0**2)-(qua1**2)+(qua2**2)-(qua3**2)
   ori(2,3) = two*((qua2*qua3)+(qua0*qua1))
   ori(3,1) = two*((qua1*qua3)+(qua0*qua2))
   ori(3,2) = two*((qua2*qua3)-(qua0*qua1))
   ori(3,3) = (qua0**2)-(qua1**2)-(qua2**2)+(qua3**2)

   call CheckOriOrtho(ori, 6)

end subroutine AxisAngToOri

!************************************************************************
!*                                                                      *
!*     LabToPri                                                         *
!*                                                                      *
!************************************************************************

! ... transform a vector from lab to principal frame

!     x(principal) =  x(lab) * ori

subroutine LabToPri(x, ori)
   implicit none
   real(8), intent(inout) :: x(3)
   real(8), intent(in)    :: ori(3,3)
   real(8) :: xt1, xt2, xt3
   xt1 = x(1)
   xt2 = x(2)
   xt3 = x(3)
   x(1) = xt1*ori(1,1) + xt2*ori(2,1) + xt3*ori(3,1)
   x(2) = xt1*ori(1,2) + xt2*ori(2,2) + xt3*ori(3,2)
   x(3) = xt1*ori(1,3) + xt2*ori(2,3) + xt3*ori(3,3)
end subroutine LabToPri

!************************************************************************
!*                                                                      *
!*     PriToLab                                                         *
!*                                                                      *
!************************************************************************

! ... transform a vector from principal to laboratory frame

!     x(lab) =  ori * x(principal)

subroutine PriToLab(x, ori)
   implicit none
   real(8), intent(inout) :: x(3)
   real(8), intent(in)    :: ori(3,3)
   real(8) :: xt1, xt2, xt3
   xt1 = x(1)
   xt2 = x(2)
   xt3 = x(3)
   x(1) = ori(1,1)*xt1 + ori(1,2)*xt2 + ori(1,3)*xt3
   x(2) = ori(2,1)*xt1 + ori(2,2)*xt2 + ori(2,3)*xt3
   x(3) = ori(3,1)*xt1 + ori(3,2)*xt2 + ori(3,3)*xt3
end subroutine PriToLab

!************************************************************************
!*                                                                      *
!*     EulerRot                                                         *
!*                                                                      *
!************************************************************************

! ... rotation of a column vector (x,y,z) using the Euler angles
!     positive rotation, right-handed coordinate system, zyz convention

!     rotated object          rotation axes          sequence
!     --------------          -------------          --------
!         axes                  rotated           p(+c)p(+b)p(+a)
!         axes                  fixed             p(+a)p(+b)p(+c)
!         function              rotated           p(-a)p(-b)p(-c)
!         function              fixed             p(-c)p(-b)p(-a)

subroutine EulerRot(mode1in,mode2in,alpha,beta,gamma,x,y,z)
   implicit none
   real(8), parameter :: degrad = 0.01745329252d0
   character(*), intent(in) :: mode1in            ! select convention
   character(*), intent(in) :: mode2in            ! select rad/deg
   real(8), intent(in) :: alpha, beta, gamma       ! Euler angles
   real(8), intent(inout) :: x, y, z               ! cartesian coordinate
   real(8) :: a, b, c, ca, cb, cc, sa, sb, sc
   real(8) :: d11, d12, d13, d21, d22, d23, d31, d32, d33
   real(8) :: x1, y1, z1
   character(len=len(mode1in)) :: mode1            ! select convention
   character(len=len(mode2in)) :: mode2            ! select rad/deg

   mode1 = mode1in
   mode2 = mode2in

   call LowerCase(mode1)
   call LowerCase(mode2)

   if (mode1 == 'axes/rot') then
      a = alpha
      b = beta
      c = gamma
   else if (mode1 == 'axes/fix') then
      a = gamma
      b = beta
      c = alpha
   else if (mode1 == 'func/rot') then
      a = -gamma
      b = -beta
      c = -alpha
   else if (mode1 == 'func/fix') then
      a = -alpha
      b = -beta
      c = -gamma
   else
      call Warn('EulerRot','wrong value of mode1',6)
   end if
   if (mode2 == 'deg') then
      a = degrad*a
      b = degrad*b
      c = degrad*c
   else if (mode2 == 'rad') then
      a = a
      b = b
      c = c
   else
      call Warn('EulerRot','wrong value of mode2',6)
   end if
   ca = cos(a)
   sa = sin(a)
   cb = cos(b)
   sb = sin(b)
   cc = cos(c)
   sc = sin(c)
   d11 = cc*cb*ca - sc*sa
   d12 = cc*cb*sa + sc*ca
   d13 = -cc*sb
   d21 = -sc*cb*ca - cc*sa
   d22 = -sc*cb*sa + cc*ca
   d23 = sc*sb
   d31 = sb*ca
   d32 = sb*sa
   d33 = cb
   x1 = x
   y1 = y
   z1 = z
   x = x1*d11 + y1*d12 + z1*d13
   y = x1*d21 + y1*d22 + z1*d23
   z = x1*d31 + y1*d32 + z1*d33
end subroutine EulerRot

!************************************************************************
!*                                                                      *
!*     EulerRotStd                                                      *
!*                                                                      *
!************************************************************************

! ... rotate tensor of standard form using Euler angles (NOTE CONVENTION NOT CLEAR)

subroutine EulerRotStd(mlmax, l, ang, Std)

   implicit none
   integer(4),  parameter :: Zero = 0.0d0
   integer(4),  intent(in)  :: mlmax                          ! upper limit of l in declaration
   integer(4),  intent(in)  :: l                              ! order of tensor
   real(8),     intent(in)  :: ang(3)                         ! Euler angles (in rad)
   complex(8),  intent(inout) :: Std(0:mlmax,-mlmax:mlmax)    ! Std before and after rotation
   integer(4) :: k, m
   complex(8) :: StdHelp(0:mlmax,-mlmax:mlmax), term, DDJKM


!  write(*,*) 'in EulerRotStd:'
   do k = -l, l
      StdHelp(l,k) = cmplx(Zero, Zero)
      do m = -l, l
         term = DDJKM(l, m, k, ang(3), ang(2), ang(1), 0) * Std(l,m)
         StdHelp(l,k) = StdHelp(l,k) + term
!        write(*,'(a,2i5,2f10.5)') 'DDJKM', k,m, DDJKM(l2, k, m, ang(1), ang(2), ang(3), 0)
!        write(*,'(a,2i5,2f10.5)') '                     term', k,m, term
      end do
!      write(*,'(a,3(i5,2f10.5,2x))') 'k, Std(l,k)   ',k, Std(l,k)
!      write(*,'(a,3(i5,2f10.5,2x))') 'k, StdHelp(l,k)',k, StdHelp(l,k)
   end do
   Std(l,-l:l)=StdHelp(l,-l:l)

end subroutine EulerRotStd

!************************************************************************
!*                                                                      *
!*     OrthoOri                                                         *
!*                                                                      *
!************************************************************************

! ... orthogonalize particle frame

subroutine OrthoOri(np, iplow, ipupp, ori, tol, unit)

   implicit none
   real(8), parameter :: One = 1.0d0
   integer(4), intent(in)    :: np
   integer(4), intent(in)    :: iplow
   integer(4), intent(in)    :: ipupp
   real(8),    intent(inout) :: ori(3,3,np)
   real(8),    intent(in)    :: tol
   integer(4), intent(in)   :: unit

   logical    :: ltext
   integer(4) :: ip
   real(8)    :: xn, yn, zn, xy, xz, yz, yt(3)

! ... check

   ltext =.true.
   do ip = iplow, ipupp
      xn = One/sqrt(ori(1,1,ip)**2+ori(2,1,ip)**2+ori(3,1,ip)**2)
      yn = One/sqrt(ori(1,2,ip)**2+ori(2,2,ip)**2+ori(3,2,ip)**2)
      zn = One/sqrt(ori(1,3,ip)**2+ori(2,3,ip)**2+ori(3,3,ip)**2)
      xy = ori(1,1,ip)*ori(1,2,ip)+ori(2,1,ip)*ori(2,2,ip)+ori(3,1,ip)*ori(3,2,ip)
      xz = ori(1,1,ip)*ori(1,3,ip)+ori(2,1,ip)*ori(2,3,ip)+ori(3,1,ip)*ori(3,3,ip)
      yz = ori(1,2,ip)*ori(1,3,ip)+ori(2,2,ip)*ori(2,3,ip)+ori(3,2,ip)*ori(3,3,ip)
      if (abs(xn-One) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad normalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'xn'' = ', xn
      else if (abs(yn-One) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad normalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'yn'' = ', yn
      else if (abs(zn-One) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad normalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'zn'' = ', zn
      else if (abs(xy) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad orthogonalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'x''y'' = ', xy
      else if (abs(xz) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad orthogonalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'x''z'' = ', xz
      else if (abs(yz) > tol) then
         if (ltext) call Warn('OrthoOri', 'bad orthogonalization', unit)
         ltext =.false.
         write(unit,'(a,i4,5x,a,e10.3)') 'particle', ip, 'y''z'' = ', yz
      end if
   end do

! ... make orthogonal

   do ip = iplow, ipupp
! ... x'-axis
      ori(1:3,1,ip) = ori(1:3,1,ip)/sqrt(ori(1,1,ip)**2+ori(2,1,ip)**2+ori(3,1,ip)**2)
! ... y'-axis
      xy = ori(1,1,ip)*ori(1,2,ip) + ori(2,1,ip)*ori(2,2,ip) + ori(3,1,ip)*ori(3,2,ip)
      yt(1:3) = ori(1:3,2,ip) - xy*ori(1:3,1,ip)
      ori(1:3,2,ip) = yt(1:3)/sqrt(yt(1)**2+yt(2)**2+yt(3)**2)
! ... z'-axis
      ori(1,3,ip) = ori(2,1,ip)*ori(3,2,ip) - ori(3,1,ip)*ori(2,2,ip)
      ori(2,3,ip) = ori(3,1,ip)*ori(1,2,ip) - ori(1,1,ip)*ori(3,2,ip)
      ori(3,3,ip) = ori(1,1,ip)*ori(2,2,ip) - ori(2,1,ip)*ori(1,2,ip)
   end do

end subroutine Orthoori

!************************************************************************
!*                                                                      *
!*     CheckOriOrtho                                                    *
!*                                                                      *
!************************************************************************

! ... check that a matrix is orthogonal

subroutine CheckOriOrtho(ori, unit)
   implicit none
   real(8),  intent(in)   :: ori(3,3)
   integer(4), intent(in) :: unit

   real(8) :: det, o11, o22, o33, o12, o13, o23
   real(8) :: acc = 1.0d-10

   det = ori(1,1)*ori(2,2)*ori(3,3)+ori(2,1)*ori(3,2)*ori(1,3) &
        +ori(3,1)*ori(1,2)*ori(2,3)-ori(3,1)*ori(2,2)*ori(1,3) &
        -ori(1,1)*ori(3,2)*ori(2,3)-ori(2,1)*ori(1,2)*ori(3,3)

   o11 = ori(1,1)**2+ori(2,1)**2+ori(3,1)**2
   o22 = ori(1,2)**2+ori(2,2)**2+ori(3,2)**2
   o33 = ori(1,3)**2+ori(2,3)**2+ori(3,3)**2
   o12 = ori(1,1)*ori(1,2)+ori(2,1)*ori(2,2)+ori(3,1)*ori(3,2)
   o13 = ori(1,1)*ori(1,3)+ori(2,1)*ori(2,3)+ori(3,1)*ori(3,3)
   o23 = ori(1,2)*ori(1,3)+ori(2,2)*ori(2,3)+ori(3,2)*ori(3,3)

   if (abs(det-1.0d0) > acc .or. abs(o11-1.0d0) > acc .or. abs(o11-1.0d0) > acc .or. abs(o11-1.0d0) > acc &
       .or. abs(o12) > acc .or. abs(o13) > acc .or. abs(o23) > acc) then
      write(unit,'(a,3g15.5)') 'deteminant - 1 = ', det-1.0d0
      write(unit,'(a,3g15.5)') 'o11 - 1 = ', o11-1.0d0
      write(unit,'(a,3g15.5)') 'o22 - 1 = ', o22-1.0d0
      write(unit,'(a,3g15.5)') 'o33 - 1 = ', o33-1.0d0
      write(unit,'(a,3g15.5)') 'o12 = ', o12
      write(unit,'(a,3g15.5)') 'o13 = ', o13
      write(unit,'(a,3g15.5)') 'o23 = ', o23
      write(unit,'(a,3g15.5)') 'ori(1,1:3) = ', ori(1,1:3)
      write(unit,'(a,3g15.5)') 'ori(2,1:3) = ', ori(2,1:3)
      write(unit,'(a,3g15.5)') 'ori(3,1:3) = ', ori(3,1:3)
      call Warn('CheckOriOrtho', 'ori not orthogonal', unit)
   end if

end subroutine CheckOriOrtho

!************************************************************************
!*                                                                      *
!*     QuaNorm                                                          *
!*                                                                      *
!************************************************************************

! ... normalize quaternions

subroutine QuaNorm(np, iplow, ipupp, qua)
   implicit none
   real(8), parameter :: one = 1.0d0
   integer(4), intent(in)    :: np
   integer(4), intent(in)    :: iplow
   integer(4), intent(in)    :: ipupp
   real(8),    intent(inout) :: qua(0:3,np)
   integer(4) :: ip
   real(8)    :: norm
   if (iplow < lbound(qua,2)) call stop('QuaNorm', 'iplow < lbound(qua,2)', 6)
   if (ipupp > ubound(qua,2)) call stop('QuaNorm', 'ipupp > ubound(qua,2)', 6)
   do ip = iplow, ipupp
      norm = one/sqrt(qua(0,ip)**2+qua(1,ip)**2+qua(2,ip)**2+qua(3,ip)**2)
      qua(0,ip) = qua(0,ip)*norm
      qua(1,ip) = qua(1,ip)*norm
      qua(2,ip) = qua(2,ip)*norm
      qua(3,ip) = qua(3,ip)*norm
   end do
end subroutine QuaNorm

!************************************************************************
!*                                                                      *
!*     AngVelToQuaVel                                                   *
!*                                                                      *
!************************************************************************

! ... calculate quaternion velocities from quaternions and angular velocities (principal frame)

subroutine AngVelToQuaVel(np, iplow, ipupp, qua, angvelo, quad)

   implicit none
   real(8), parameter :: half = 0.5d0
   integer(4), intent(in)  :: np
   integer(4), intent(in)  :: iplow
   integer(4), intent(in)  :: ipupp
   real(8),    intent(in)  :: qua(0:3,np)
   real(8),    intent(in)  :: angvelo(1:3,np)
   real(8),    intent(out) :: quad(0:3,np)
   integer :: ip

   if (iplow < lbound(qua,2)) call stop('AngVelToQuaVel', 'iplow < lbound(qua,2)', 6)
   if (ipupp > ubound(qua,2)) call stop('AngVelToQuaVel', 'ipupp > ubound(qua,2)', 6)
   do ip = iplow, ipupp
      quad(0,ip) = half*(-qua(1,ip)*angvelo(1,ip)+qua(2,ip)*angvelo(2,ip)-qua(3,ip)*angvelo(3,ip))
      quad(1,ip) = half*(+qua(0,ip)*angvelo(1,ip)-qua(3,ip)*angvelo(2,ip)-qua(2,ip)*angvelo(3,ip))
      quad(2,ip) = half*(-qua(3,ip)*angvelo(1,ip)-qua(0,ip)*angvelo(2,ip)+qua(1,ip)*angvelo(3,ip))
      quad(3,ip) = half*(+qua(2,ip)*angvelo(1,ip)+qua(1,ip)*angvelo(2,ip)+qua(0,ip)*angvelo(3,ip))
   end do

end subroutine AngVelToQuaVel

!************************************************************************
!*                                                                      *
!*     QuaVelToAngVel                                                   *
!*                                                                      *
!************************************************************************

! ... calculate angular velocities (principal frame) from quaternions and quaternion velocities

subroutine QuaVelToAngVel(np, iplow, ipupp, qua, quad, angvelo)

   implicit none
   real(8), parameter :: two = 2.0d0
   integer(4), intent(in)  :: np
   integer(4), intent(in)  :: iplow
   integer(4), intent(in)  :: ipupp
   real(8),    intent(in)  :: qua(0:3,np)
   real(8),    intent(in)  :: quad(0:3,np)
   real(8),    intent(out) :: angvelo(1:3,np)
   integer :: ip

   do ip = iplow, ipupp
      angvelo(1,ip) = two*(-qua(3,ip)*quad(2,ip)+qua(0,ip)*quad(1,ip)+qua(2,ip)*quad(3,ip)-qua(1,ip)*quad(0,ip))
      angvelo(2,ip) = two*(-qua(0,ip)*quad(2,ip)-qua(3,ip)*quad(1,ip)+qua(1,ip)*quad(3,ip)+qua(2,ip)*quad(0,ip))
      angvelo(3,ip) = two*(+qua(1,ip)*quad(2,ip)-qua(2,ip)*quad(1,ip)+qua(0,ip)*quad(3,ip)-qua(3,ip)*quad(0,ip))
   end do

end subroutine QuaVelToAngVel

!************************************************************************
!*                                                                      *
!*     Random                                                           *
!*                                                                      *
!************************************************************************

! ... return a random number in the range of 0 < ran < 1
!     iseed should not be equal to zero, and the first seed should be negative
!     "Numerical recipes in Fortran 90" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1992.
!     modified from function ran

!“Minimal” random number generator of Park and Miller combined with a Marsaglia shift sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). This fully portable, scalar generator has the “traditional” (not Fortran 90) calling sequence with a random deviate as the returned function value: call with idum a negative integer to initialize; thereafter, do not alter idum except to reinitialize. The period of this generator is about 3.1 × 10^18 .

module Random_Module
   integer, parameter :: k4b=selected_int_kind(9) ! = 4 on intel fortran and gfortran
   real(8) :: am
   integer(k4b) :: ix=-1,iy=-1
end module Random_Module

function Random(idum)
   use Random_Module
   implicit none
   integer(k4b), intent(inout) :: idum
   real(8) :: Random
   integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
   integer(k4b)   :: k
   if (idum <= 0 .or. iy < 0) then           !initialize.
      am=nearest(1.0,-1.0)/im
      iy=ior(ieor(888889999,abs(idum)),1)
      ix=ieor(777755555,abs(idum))
      idum=abs(idum)+1                          !set idum positive.
   end if
   ix=ieor(ix,ishft(ix,13))                  !marsaglia shift sequence with period 2^32 − 1.
   ix=ieor(ix,ishft(ix,-17))
   ix=ieor(ix,ishft(ix,5))
   k=iy/iq                                   !park-miller sequence by schrage’s method, period 2^31 − 2.
   iy=ia*(iy-k*iq)-ir*k
   if (iy < 0) iy=iy+im
   Random=am*ior(iand(im,ieor(ix,iy)),1)     !combine the two generators with masking to ensure nonzero value.
end function Random

!************************************************************************
!*     Random                                                           *
!*                                                                      *
!************************************************************************

! ... return a random number in the range of 0 to 1
!     iseed should not be equal to zero
!     by Learmonth and Lewis 1973 (32 bits integer)

!real(8) function Random(iseed)
   !implicit none
   !integer(4), parameter :: k = 16087, l = 0, nb = 31
   !real(8), parameter :: f = 2.0d0**(-nb)
   !integer(4), intent(inout) :: iseed
   !iseed = iseed*k+l
   !iseed = ibclr(iseed,nb)
   !Random = f*iseed
!end function Random

!************************************************************************
!*                                                                      *
!*     GauRandom                                                        *
!*                                                                      *
!************************************************************************

! ... return a normally distributed stochastic variable
!     zero mean and unit variance. box-muller mehtod
!     random is explicitly inlined

real(8) function GauRandom(iseed)
   implicit none
   real(8), parameter :: half = 0.5d0 , one = 1.0d0 , two = 2.0d0
   integer(4), parameter :: k = 16087, l = 0, nb = 31
   real(8),    parameter :: f = 2.0d0**(-nb)
   integer(4), intent(inout) :: iseed  ! seed to random number generator
   real(8) :: rx, ry, r, fac
   logical, save :: flag =.true.
   real(8), save :: save

   if (flag) then
      do
         iseed = ibclr(iseed*k+l,nb)
         rx = two*(f*iseed-half)
         iseed = ibclr(iseed*k+l,nb)
         ry = two*(f*iseed-half)
         r = rx**2+ry**2
         if (r < one) exit
      end do
      fac = sqrt(-two*log(r)/r)
      GauRandom = rx*fac
      save = ry*fac
      flag = .false.
   else
      GauRandom = save
      flag = .true.
   end if
end function GauRandom

!************************************************************************
!*                                                                      *
!*     CirRandom                                                        *
!*                                                                      *
!************************************************************************

! ... return two random numbers x,y on a unit circle

subroutine CirRandom(iseed,x,y)
   implicit none
   real(8), parameter ::  pi2 = 6.283185308d0
   integer(4), intent(inout) :: iseed  ! seed to random number generator
   real(8),    intent(out)   :: x, y   ! two random numbers on a unit circle
   real(8) :: theta, Random
   theta = pi2*Random(iseed)
   x = cos(theta)
   y = sin(theta)
end subroutine CirRandom

!************************************************************************
!*                                                                      *
!*     SphRandom                                                        *
!*                                                                      *
!************************************************************************

! ... return three random numbers x,y,z on a unit sphere

subroutine SphRandom(iseed,x,y,z)
   implicit none
   real(8), parameter ::  pi2 = 6.283185308d0
   integer(4), intent(inout) :: iseed   ! seed to random number generator
   real(8),    intent(out)   :: x, y, z ! three random numbers on a unit sphere
   real(8) :: ct, st, psi, Random
   ct = 2.0d0*(Random(iseed)-0.5d0)
   st = sqrt(1.0d0-ct*ct)
   psi = pi2*Random(iseed)
   x = st*cos(psi)
   y = st*sin(psi)
   z = ct
end subroutine SphRandom

!************************************************************************
!*                                                                      *
!*     InvFlt                                                           *
!*                                                                      *
!************************************************************************

! ... return the inverse of a real variable

real(8) function InvFlt(real)
   implicit none
   real(8), intent(in) :: real
   InvFlt = 0.0d0
   if (real > 0.0d0) InvFlt = 1.0d0/real
end function InvFlt

!************************************************************************
!*                                                                      *
!*     GetRelDiff                                                       *
!*                                                                      *
!************************************************************************

! ... calculate relative differences

real(8) function GetRelDiff(a, b)
   implicit none
   real(8), intent(in) :: a, b
   GetRelDiff = 0.0d0
!   if (abs(a) /= 0) GetRelDiff = (b-a)/abs(a)
   if ((abs(a) > 1.0d-15) .and. (abs(b) > 1.0d-15)) GetRelDiff = (b-a)/abs(a)
end function GetRelDiff

!************************************************************************
!*                                                                      *
!*     CpuAdd                                                           *
!*                                                                      *
!************************************************************************

! ... add and write total cpu time elapsed since start

subroutine CpuAdd(txwhattodo,label,level,unit)

   implicit none
   integer(4), parameter :: mnlabel = 100
   character(*), intent(in) :: txwhattodo !'start'
                                           !'stop'
                                           !'write'
                                           !'interrupt' disable the effect of 'start', 'stop', and 'write'
                                           !'resume' enable them again
   character(*), intent(in)    :: label    ! name of section
   integer(4),   intent(in)    :: level    ! size of indentation (0-
   integer(4),   intent(in)    :: unit     ! output unit
   integer(4),    save :: nlabel = 0       ! number of sections in the list
   real(8),       save :: cpustart(mnlabel)
   real(8),       save :: cputot(mnlabel) = 0.0d0
   character(40), save :: xlabel(mnlabel), s = '                          '
   integer(4),    save :: xlevel(mnlabel)
   logical,       save :: interrupt = .false.
   character(2)  :: str
   integer(4) :: i
   character(LEN=len(txwhattodo))   :: whattodo
   real(8) :: SecondsSinceStart, tot

   whattodo = txwhattodo

   call LowerCase(whattodo)

! ... handel interrupt tasks

   if (whattodo == 'interrupt') then
      interrupt =.true.
   else if (whattodo == 'resume' ) then
      interrupt = .false.
      return
   end if

! ... test for interrupt

   if (interrupt) return

! ... get the order of the incoming label, possibly extending the list

   if (whattodo /= 'write') then
      do i = 1,nlabel
         if (trim(label) == xlabel(i)) exit
      end do
      if (i > nlabel) then
         if (i > mnlabel) call stop('CpuAdd','i>mnlabel',unit)
         if (whattodo == 'start') then
            nlabel = i
            xlabel(nlabel) = trim(label)
            xlevel(nlabel) = level
         else if (whattodo == 'stop') then
            write(unit,'(3a)') 'label',trim(label),'label'
            call stop('CpuAdd','no matching start',unit)
         else
            call stop('CpuAdd','error in whattodo',unit)
         end if
      end if
   end if

! ... do the work

   if (whattodo == 'start') then
      cpustart(i) = SecondsSinceStart()
   else if (whattodo == 'stop') then
      cputot(i) = cputot(i)+SecondsSinceStart()-cpustart(i)
   else if (whattodo == 'write') then
      write(str,'(i2)') 15+max(0,maxval(len_trim(xlabel(1:nlabel))))
      tot = sum(cputot(1:nlabel), 1, xlevel(1:nlabel)==0)  ! total cpu time by summing top level
      if (tot == 0.0d0) tot = 1.0d0                        ! avoid overflow
      write(unit,*)
      if (tot < 3600.0) then
         write(unit,'(t10,a,t'//str//',a,5x,a)') 'section','cpu time (s)','fraction'
         write(unit,'(t10,a,t'//str//',a,5x,a)') '-------','------------','--------'
         write(unit,'(t10,a,t'//str//',f10.2,3x,f10.2)') (s(1:xlevel(i))//xlabel(i),cputot(i),cputot(i)/tot,i = 1,nlabel)
      else
         write(unit,'(t10,a,t'//str//',a,5x,a)') 'section','cpu time (h)','fraction'
         write(unit,'(t10,a,t'//str//',a,5x,a)') '-------','------------','--------'
         write(unit,'(t10,a,t'//str//',f10.2,3x,f10.2)') (s(1:xlevel(i))//xlabel(i),cputot(i)/3600.0,cputot(i)/tot,i = 1,nlabel)
      end if
   else
      call Warn('CpuAdd','error in whattodo',unit)
   end if

end subroutine CpuAdd

!************************************************************************
!*                                                                      *
!*     CpuLeft                                                          *
!*                                                                      *
!************************************************************************

! ... check to see if the execution should be stopped

!     if tcycle > tleft then the program is stopped, where
!        tcycle is the cpu time between the last two calls of CpuLeft
!        tleft  is maxcpu-cpu, where cpu is the used cpu time

!     the routine should be called at the beginning of the loop

subroutine CpuLeft(maxcpu,unit)
   implicit none
   integer(4), intent(in) :: maxcpu   ! max cpu in seconds
   integer(4), intent(in) :: unit     ! output unit
   logical, save :: first = .true.
   real(8), save :: t1
   real(8) :: t2, SecondsSinceStart, tcycle, tleft

   t2 = SecondsSinceStart()
   if (first) then
      first = .false.
   else
      tcycle = t2-t1
      tleft = real(maxcpu-t2)
      if (tcycle>tleft) then
         call WriteHead(2,'abort',unit)
         write(unit,'(a,f8.2)') 'cpu time left       (h) =',tleft/3600.0
         write(unit,'(a,f8.2)') 'cpu time last cycle (h) =',tcycle/3600.0
         write(unit,'(a,f8.2)') 'total cpu time used (h) =',t2/3600.0
         write(unit,*)
         write(unit,*)'program stopped due to too little cpu time left'
         stop 1
      end if
   end if
   t1 = t2
end subroutine CpuLeft

!************************************************************************
!*                                                                      *
!*     CpuTot                                                           *
!*                                                                      *
!************************************************************************

! ... write total cpu time elapsed since start

subroutine CpuTot(unit)
   implicit none
   integer(4), intent(in) :: unit
   real(8) :: SecondsSinceStart, t
   t = SecondsSinceStart()
   write(unit,*)
   if (t < 3600.0) then
      write(unit,'(a,2f12.2)') 'total cpu time since start (s)',t
   else
      write(unit,'(a,2f12.2)') 'total cpu time since start (h)',t/3600.0
   end if
end subroutine CpuTot

!************************************************************************
!*                                                                      *
!*     SecondsSinceStart                                                *
!*                                                                      *
!************************************************************************

! ... return cpu time used since start in seconds

real(8) function SecondsSinceStart()
   implicit none
!  SecondsSinceStart = 1.0D-3 * mscpu()        ! Univac
!  integer(4), tused                ! Nord
!  SecondsSinceStart = 20.D-3 * tused(0)       ! Nord
!  SecondsSinceStart = secvax()                ! Vax
!  istat = sys&gettime(SecondsSinceStart,wait) ! Fps164
!  SecondsSinceStart = 1.0d-2 * mclock()       ! Ibm Aix / Risc 6000
!  SecondsSinceStart = 1.0d-3 * mclock()       ! Intel
   call cpu_time(SecondsSinceStart)            ! fortran 95
end function SecondsSinceStart

!************************************************************************
!*                                                                      *
!*     WriteVec                                                         *
!*                                                                      *
!************************************************************************

! ... write a vector

subroutine WriteVec(n,vec,label,ilow,iupp,ifreq,unit)
   implicit none
   integer(4),   intent(in) :: n
   real(8),      intent(in) :: vec(n)
   character(*), intent(in) :: label
   integer(4),   intent(in) :: ilow, iupp, ifreq, unit
   integer(4) :: i
   write(unit,*)
   write(unit,'(a,a,5x,a,i4,a,i4,a,i4,a,/)') 'vector: ',label,'(',ilow,',',iupp,',',ifreq,')'
   write(unit,'(5g16.8)') (vec(i),i = ilow,iupp,ifreq)
end subroutine WriteVec

!************************************************************************
!*                                                                      *
!*     WriteMat                                                         *
!*                                                                      *
!************************************************************************

! ... write a matrix

subroutine WriteMat(nrow,ncol,amat,label,irl,iru,irf,icl,icu,icf,unit)
   implicit none
   integer(4),   intent(in) :: nrow, ncol
   real(8),      intent(in) :: amat(nrow,ncol)
   character(*), intent(in) :: label
   integer(4),   intent(in) :: irl, iru, irf, icl, icu, icf, unit
   integer(4) :: ncolum, nparts, m, i, j, ilow, iupp
   write(unit,*)
   write(unit,'(a,a,5x,2(a,i3,a,i3,a,i3,a))') 'matrix:   ',label, &
   'rows(',irl,',',iru,',',irf,'), ','columns(',icl,',',icu,',',icf,')'
   ncolum = 5
   nparts = 1+(icu-icl)/(icf*ncolum)
   do m = 1,nparts
      ilow = icl+icf*ncolum*(m-1)
      iupp = min0(icu,ilow-1+icf*ncolum)
      write(unit,*)
      write(unit,'(5i15)') (i,i = ilow,iupp,icf)
      write(unit,*)
      do j = irl,iru,irf
         write(unit,'(i3,2x,5g15.7)') j,(amat(j,i),i = ilow,iupp,icf)
      end do
   end do
end subroutine WriteMat

!************************************************************************
!*                                                                      *
!*     WriteStd                                                         *
!*                                                                      *
!************************************************************************

! ... write a variable of standard form

subroutine WriteStd(mlmax,lmax,q,label,iopt,unit)
   implicit none
   integer(4), parameter :: tcenter(0:1) = [57, 20]
   real(8), parameter :: Small = 1.0d-15
   integer(4),   intent(in)  :: mlmax
   integer(4),   intent(in)  :: lmax
   complex(8)  , intent(in)  :: q(0:mlmax,-mlmax:mlmax)
   character(*), intent(in)  :: label
   integer(4),   intent(in)  :: iopt
   integer(4),   intent(in)  :: unit
   integer(4) :: l, m, nchar
   character(2) :: string
   nchar = len_trim(label)
   write(string,'(i2)') tcenter(iopt)-nchar/2
   write(unit,*)
   if (iopt == 0) then                                              ! short form
      write(unit,'(a,'//string//'x,a)') 'l', trim(label)
      write(unit,'(a,46x,3(2f10.5,3x))') '0', q(0, 0:0)
      write(unit,'(a,23x,3(2f10.5,3x))') '1', q(1,-1:1)
      write(unit,'(a,5(2f10.5,3x))')     '2', q(2,-2:2)
   else if (iopt == 1) then                                         ! long form
      write(unit,'(a,'//string//'x,a)') '  l   m', trim(label)
      do l = 0, lmax
         do m = -l, l
            if (real(q(l,m)*conjg(q(l,m))) > Small)  write(unit,'(2i4,2g20.10)') l, m, q(l,m)
         end do
      end do
   else
      call Stop('WriteStd', 'invalid value of iopt', 6)
   end if
end subroutine WriteStd

!************************************************************************
!*                                                                      *
!*     WriteFront                                                       *
!*                                                                      *
!************************************************************************

! ... write a front

subroutine WriteFront(txProgName, txShortInfo, txVerDate, txAuthorMain, txAuthorCont, unit)

   use MollibModule, only: Center
   implicit none
   integer(4), parameter :: nw = 110           ! width of frame
   integer(4), parameter :: nlen = nw - 4      ! width for text
   integer(4), parameter :: mnsubstr = 100     ! maximum number of substrings
   character(*), intent(in) :: txProgName      ! name of program
   character(*), intent(in) :: txShortInfo     ! short description (\ signals line break)
   character(*), intent(in) :: txVerDate       ! version and ate
   character(*), intent(in) :: txAuthorMain    ! main author
   character(*), intent(in) :: txAuthorCont    ! contributing authors (\ signal line break)
   integer(4),   intent(in) :: unit            ! unit
   integer(4)      :: i, nsubstr, ilow(mnsubstr), iupp(mnsubstr)
   character(3) :: fmt

   write(fmt,'(i3)') nw
   write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
   write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'

   write(unit,'('//fmt//'a)') '**',Center(nlen,txProgName),'**'
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'

   if (len(txShortInfo) > 0) then
      call SubStr(txShortInfo, '\', nsubstr, ilow, iupp)
      if (nsubstr > mnsubstr) call Stop('WriteFront', 'nsubstr > mnsubstr', unit)
      do i = 1, nsubstr
         write(unit,'('//fmt//'a)') '**',Center(nlen,txShortInfo(ilow(i):iupp(i))),'**'
      end do
      write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
   end if

   write(unit,'('//fmt//'a)') '**',Center(nlen,txVerDate),'**'
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'

   write(unit,'('//fmt//'a)') '**',Center(nlen,txAuthorMain),'**'
   write(unit,'('//fmt//'a)') '**',Center(nlen,'Physical Chemistry     '),'**'
   write(unit,'('//fmt//'a)') '**',Center(nlen,'Department of Chemistry'),'**'
   write(unit,'('//fmt//'a)') '**',Center(nlen,'Lund University        '),'**'
   write(unit,'('//fmt//'a)') '**',Center(nlen,'P.O. Box 124           '),'**'
   write(unit,'('//fmt//'a)') '**',Center(nlen,'S-221 00 Lund, SWEDEN  '),'**'

   if (len(txAuthorCont) > 0) then
      write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
      write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
      write(unit,'('//fmt//'a)') '**',Center(nlen,'with contributions from'),'**'
      write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
      call SubStr(txAuthorCont, '\', nsubstr, ilow, iupp)
      if (nsubstr > mnsubstr) call Stop('WriteFront', 'nsubstr > mnsubstr', unit)
      do i = 1, nsubstr
         write(unit,'('//fmt//'a)') '**',Center(nlen,txAuthorCont(ilow(i):iupp(i))),'**'
      end do
   end if

   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
   write(unit,'('//fmt//'a)') '**',(' ',i = 1,nlen),'**'
   write(unit,'('//fmt//'a)') ('*',i = 1,nw)
   write(unit,'('//fmt//'a)') ('*',i = 1,nw)
   write(unit,*)

end subroutine WriteFront

!************************************************************************
!*                                                                      *
!*     WriteHead                                                        *
!*                                                                      *
!************************************************************************

! ... write a heading

subroutine WriteHead(ilevel,string,unit)
   use MollibModule, only: Center, SpaceOut
   implicit none
   integer(4), parameter ::  nw = 110      ! width of frame
   integer(4), parameter ::  nlen = nw - 2 ! maximal width of text
   integer(4),   intent(in) :: ilevel      ! 1, 2, or 3
   character(*), intent(in) :: string      ! text
   integer(4),   intent(in) :: unit        ! unit
   integer(4) :: i
   character(3) :: fmt

   write(fmt,'(i3)') nw
   write(unit,'()')
   if (ilevel == 1) then
      write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
      write(unit,'(a,t'//fmt//',a)') '*','*'
      write(unit,'('//fmt//'a)')     '*',Center(nlen,SpaceOut(trim(string))),'*'
      write(unit,'(a,t'//fmt//',a)') '*','*'
      write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
   else if (ilevel == 2) then
      write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
      write(unit,'('//fmt//'a)')     '*',Center(nlen,trim(string)),'*'
      write(unit,'('//fmt//'a)')    ('*',i = 1,nw)
   else if (ilevel == 3) then
      write(unit,'('//fmt//'a)')    ('.',i = 1,nw)
      write(unit,'('//fmt//'a)')     '.',Center(nlen,trim(string)),'.'
      write(unit,'('//fmt//'a)')    ('.',i = 1,nw)
   end if
   write(unit,'()')
end subroutine WriteHead

!************************************************************************
!*                                                                      *
!*     WriteDateTime                                                    *
!*                                                                      *
!************************************************************************

! ... write date and time

subroutine WriteDateTime(unit)
   implicit none
   integer(4), intent(in) :: unit   ! output unit
   character(8)  :: date
   character(10) :: time
   call date_and_time(date,time)
   write(unit,'(a,5x,a)') date(1:4)//'-'//date(5:6)//'-'//date(7:8),time(1:2)//':'//time(3:4)//':'//time(5:6)
end subroutine WriteDateTime

!************************************************************************
!*                                                                      *
!*     WriteIOStat                                                      *
!*                                                                      *
!************************************************************************

! ... write value of iostat and take appropriate action

subroutine WriteIOStat(name,text,iostat,iopt,unit)
   use MollibModule, only: Center, SpaceOut
   implicit none
   integer(4), parameter ::  nw = 110      ! width of frame
   integer(4), parameter ::  nlen = nw - 2 ! maximal width of text
   character(*), intent(in) :: name        ! name of routine
   character(*), intent(in) :: text        ! text to be displayed
   integer(4),   intent(in) :: iostat      ! value of iostat
   integer(4),   intent(in) :: iopt       ! determine the action
                                           !=0 no action
                                           !=1 stop if iostat>0 (an error)
                                           !=2 stop if iostat/=0 (an error or end-of-file)
   integer(4),   intent(in) :: unit        ! output unit
   integer(4) :: i
   character(3) :: fmt
   character(nlen) :: string, txiostat

   write(fmt,'(i3)') nw
   string = 'message from '//trim(name)

   write(txiostat,'(i3)') iostat
   txiostat='iostat = '//txiostat

   write(unit,'('//fmt//'a)')    (':',i = 1,nw)
   write(unit,'(a,t'//fmt//',a)') ':',':'
   write(unit,'('//fmt//'a)')     ':',Center(nlen,SpaceOut(trim(string))),':'
   write(unit,'(a,t'//fmt//',a)') ':',':'
   write(unit,'('//fmt//'a)')     ':',Center(nlen,text),':'
   write(unit,'(a,t'//fmt//',a)') ':',':'
   write(unit,'('//fmt//'a)')     ':',Center(nlen,txiostat),':'
   write(unit,'(a,t'//fmt//',a)') ':',':'
   write(unit,'('//fmt//'a)')    (':',i = 1,nw)

   if ((iopt==1 .and. iostat>=1) .or. (iopt==2 .and. iostat/=0)) stop 1
end subroutine WriteIOStat

!************************************************************************
!*                                                                      *
!*     Warn                                                             *
!*                                                                      *
!************************************************************************

! ... write a warning message

subroutine Warn(name,text,unit)
   use MollibModule, only: Center, SpaceOut
   implicit none
   integer(4), parameter ::  nw = 110      ! width of frame
   integer(4), parameter ::  nlen = nw - 2 ! maximal width of text
   character(*), intent(in) :: name        ! name of routine
   character(*), intent(in) :: text        ! text to be displayed
   integer(4),   intent(in) :: unit        ! output unit
   integer(4) :: i
   character(3) :: fmt
   character(nlen) :: string

   write(fmt,'(i3)') nw
   string = 'warning from '//trim(name)

   write(unit,'('//fmt//'a)')    ('%',i = 1,nw)
   write(unit,'(a,t'//fmt//',a)') '%','%'
   write(unit,'('//fmt//'a)')     '%',Center(nlen,SpaceOut(trim(string))),'%'
   write(unit,'(a,t'//fmt//',a)') '%','%'
   write(unit,'('//fmt//'a)')     '%',Center(nlen,text),'%'
   write(unit,'(a,t'//fmt//',a)') '%','%'
   write(unit,'('//fmt//'a)')    ('%',i = 1,nw)
end subroutine Warn

!************************************************************************
!*                                                                      *
!*     Stop                                                             *
!*                                                                      *
!************************************************************************

! ... write a stop message and stop process

subroutine Stop(name,text,unit)
   use, intrinsic :: iso_fortran_env, only : ustdout=>output_unit
   implicit none
   integer(4), parameter ::  nw = 110      ! width of frame of box
   integer(4), parameter ::  nlen = nw - 2 ! maximal width of text
   character(*), intent(in) :: name        ! name of routine
   character(*), intent(in) :: text        ! text to be displayed
   integer(4),   intent(in) :: unit        ! output unit
   integer(4) :: i
   character(3) :: fmt

   write(fmt,'(i3)') nw
   call StopUnit(ustdout)
   if (unit /= ustdout) call StopUnit(unit)
#if defined (_PAR_)
   call par_finalize
#endif
   stop 1

contains

!........................................................................

subroutine StopUnit(unit)
   use MollibModule, only: Center, SpaceOut
   integer(4),   intent(in) :: unit    ! output unit
   character(nlen) :: string

   string = 'stop in '//trim(name)
   write(unit,'('//fmt//'a)')    ('#',i = 1,nw)
   write(unit,'(a,t'//fmt//',a)') '#','#'
   write(unit,'('//fmt//'a)')     '#',Center(nlen,SpaceOut(trim(string))),'#'
   write(unit,'(a,t'//fmt//',a)') '#','#'
   write(unit,'('//fmt//'a)')     '#',Center(nlen,text),'#'
   write(unit,'(a,t'//fmt//',a)') '#','#'
   write(unit,'('//fmt//'a)')    ('#',i = 1,nw)
end subroutine StopUnit

!........................................................................

end subroutine Stop

!************************************************************************
!*                                                                      *
!*     LowerCase                                                        *
!*                                                                      *
!************************************************************************

! ... change a string to lower case

subroutine LowerCase(string)
   implicit none
   character(*), intent(inout) :: string   ! string to be converted
   integer(4) :: i, ldiff
   ldiff = ichar('a')-ichar('A')
   do i = 1,len(string)
      if (string(i:i) >= 'A' .and. string(i:i) <= 'Z') string(i:i) = char(ichar(string(i:i))+ldiff)
   end do
end subroutine LowerCase

!************************************************************************
!*                                                                      *
!*     UpperCase                                                        *
!*                                                                      *
!************************************************************************

! ... change a string to upper case

subroutine UpperCase(string)
   implicit none
   character(*), intent(inout) :: string   ! string to be converted
   integer(4) :: i, ldiff
   ldiff = ichar('A')-ichar('a')
   do i = 1,len(string)
      if (string(i:i) >= 'a' .and. string(i:i) <= 'z') string(i:i) = char(ichar(string(i:i))+ldiff)
   end do
end subroutine UpperCase

!************************************************************************
!*                                                                      *
!*     SubStr                                                           *
!*                                                                      *
!************************************************************************
!
! ... determine the number of substrings of a string

subroutine SubStr(string,sep,no,ilow,iupp)

   implicit none
   character(*), intent(in)  :: string  ! string to be examined
   character(1), intent(in)  :: sep     ! separation character
   integer(4),   intent(out) :: no      ! number of substrings
   integer(4),   intent(out) :: ilow(*) ! lower index of a substring
   integer(4),   intent(out) :: iupp(*) ! upper index of a substring
   integer(4) :: length, i

   length = len_trim(string)
   no = 0
   if (string(1:1) /= sep) then
      no = no+1
      ilow(no) = 1
   end if
   do i = 2,length
      if (string(i:i) == sep .and. string(i-1:i-1) /= sep) then
         iupp(no) = i-1
      end if
      if (string(i:i) /= sep .and. string(i-1:i-1) == sep) then
         no = no+1
         ilow(no) = i
      end if
   end do
   if (string(length:length) /= sep) iupp(no) = length

end subroutine SubStr

!************************************************************************
!*                                                                      *
!*     SignMagn                                                         *
!*                                                                      *
!************************************************************************

! ... get sign and magnitude of a real number

subroutine SignMagn(xin, isign, imagn)
   implicit none
   real(8), parameter :: zero = 0.0d0 , one = 1.0d0 , ten = 10.0d0
   real(8), intent(in)  :: xin      ! real number
   integer(4), intent(out) :: isign  ! sign
   integer(4), intent(out) :: imagn  ! magnitude
   real(8) :: x
   isign = 1                         ! initialize
   imagn = 0                         ! initialize
   x = xin
   if (x < zero) then                ! x < 0 (invert x)
      isign = -1
      x = -x
   endif
   if (x == zero) then                    ! x = 0 (requires special care)
      return
   end if
   if ((zero < x) .and. (x < one)) then   ! 0 < x < 1
      do
         imagn = imagn - 1
         x = x*ten
         if (x >= one) exit
      end do
   else if (x >= ten) then                ! x >= 10
      do
         imagn = imagn + 1
         x = x/ten
         if (x < ten) exit
      end do
   else                                   ! 1 =< x < 10
      imagn = 0
   endif
end subroutine SignMagn

!************************************************************************
!*                                                                      *
!*     Advance                                                          *
!*                                                                      *
!************************************************************************

! ... read from an external unit until a string match is found

subroutine Advance(targetstring,unit,str,ok)

   implicit none
   character(*), intent(in)  :: targetstring
   integer,      intent(in)  :: unit
   character(*), intent(out) :: str
   logical,      intent(out) :: ok
   integer :: ios

   ok = .false.
   do
      read(unit,'(a)',iostat = ios) str             ! fix  1 -> unit   2010-11-29
      if (ios /= 0) then
!        write(*,*) 'ios = ',ios
         return
      endif
      if (index(str,targetstring) > 0) exit
   end do
   ok = .true.
end subroutine Advance

!************************************************************************
!*                                                                      *
!*     Plot                                                             *
!*                                                                      *
!************************************************************************

! ... plot a function

subroutine Plot(txsa,nin,sain,txyscale,x0,x1,y0,y1,unit)

   implicit none
   character(*), intent(in) :: txsa      ! heading of the plot
   integer(4),   intent(in) :: nin       ! number of values to be plotted
   real(8),      intent(in) :: sain(nin) ! function to be plotted
   character(*), intent(in) :: txyscale  ! 'both' -> y-scale is taken from y0 and y1
                                         ! 'low' -> lower y-scale is taken from y0, the other automacic
                                         ! 'upp' -> upper y-scale is taken from y1, the other automatic
                                         ! 'none' -> automatic scaling
   real(8),      intent(in) :: x0, x1    ! lower and upper limit of x-axis
   real(8),      intent(in) :: y0, y1    ! lower and upper limit of y-axis
   integer(4),   intent(in) :: unit      ! output unit

   integer(4), parameter :: nout = 100 , ny = 40
   real(8) :: saout(nout)
   character(1) :: char(nout)
   integer(4) :: i, imax, imin, j, m
   real(8) :: amax, amin, ymin, ymax, dy, yymin, yymax, div

! ... transform sain to saout

   div = dfloat(nin)/dfloat(nout)
   do i = 1,nout
      amax = i*div
      amin = (i-1)*div
      imax = int(amax-1.0d-6)+1
      imin = int(amin-1.0d-6)+1
      if (imax == imin) then
          saout(i) = sain(imax)
      else
         saout(i) = ((imin)-amin)*sain(imin)+(amax-(imax-1))*sain(imax)
         do j = imin+1,imax-1
            saout(i) = saout(i)+sain(j)
         end do
         saout(i) = saout(i)/(amax-amin)
      end if
   end do

! ... determine ymin, ymax and dy

   if (txyscale == 'both') then
      ymin = y0
      ymax = y1
   else if (txyscale == 'low') then
      ymin = y0
      ymax = maxval(saout(1:nout))
   else if (txyscale == 'upp') then
      ymin = minval(saout(1:nout))
      ymax = y1
   else if (txyscale == 'none') then
      ymin = minval(saout(1:nout))
      ymax = maxval(saout(1:nout))
   end if
   ymax = ymax*(1+sign(1.0d-6,ymax))
   ymin = ymin*(1-sign(1.0d-6,ymin))
   dy = (ymax-ymin)/ny

! ... plot

   write(unit,*)
   write(unit,*)
   write(unit,'(6x,a,5x,a,e10.3,5x,a,e10.3)') txsa,'y-min = ',ymin,'y-max = ',ymax
   write(unit,*)
   do m = ny,1,-1
      yymin = ymin + (m-1)*dy
      yymax = ymin + m*dy
      imax = 0
      do i = 1,nout
         if (saout(i) >= yymin .and. saout(i) <= yymax) then
            char(i) = '*'
            imax = i
         else
            char(i) = ' '
         end if
      end do
      write(unit,'(3x,a,100a)') '!',(char(i),i = 1,imax)
   end do
   write(unit,"(3x,'!',20('----+'),/' ','  !',10(9x,'!'))")
   write(unit,'(f10.3,t95,f10.3)') x0,x1

end subroutine Plot

!************************************************************************
!*                                                                      *
!*     BrentMod                                                         *
!*                                                                      *
!************************************************************************

! ... get the minimum of a function
!     "Numerical recipes" by Press, Flannery, Teukolsky, and Vetterling, Cambridge, 1986.
!     modified

real(8) function BrentMod(ax,bx,cx,f,tol,xthr,xmin)
   implicit none
   integer(4), parameter :: itmax = 100
   real(8), parameter :: cgold = .3819660d0 , zeps = 1.0d-10
   real(8), intent(in) :: ax, bx, cx   ! bx between ax and cx, f(bx) smallest
   real(8), external :: F              ! function to use
   real(8), intent(in) :: tol          ! dersired fractional tolerance
   real(8), intent(in) :: xthr         ! threshold for minimum searching
   real(8), intent(out) :: xmin        ! x value of the minimum
   integer(4) :: iter
   real(8) :: a, b, d, e, etemp, p, q, r, u, v, w, x, xm, fu, fv, fw, fx
   real(8) :: tol1, tol2

   a = min(ax,cx)
   b = max(ax,cx)
   v = bx
   w = v
   x = v
   e = 0.0d0
   fx = F(x)
   if (fx <= xthr) then  ! exit if f <= xthr
!          write(*,'(a,f12.5)') 'max exit 1', -fx
          goto 3
   end if
   fv = fx
   fw = fx
   do iter = 1,itmax
      xm = 0.5d0*(a+b)
      tol1 = tol*abs(x)+zeps
      tol2 = 2.*tol1
      if (abs(x-xm) <= (tol2-0.5d0*(b-a))) goto 3
      if (abs(e) > tol1) then
         r = (x-w)*(fx-fv)
         q = (x-v)*(fx-fw)
         p = (x-v)*q-(x-w)*r
         q = 2.0d0*(q-r)
         if (q > 0.0d0) p = -p
         q = abs(q)
         etemp = e
         e = d
         if (abs(p) >= abs(0.5d0*q*etemp).or.p <= q*(a-x).or. p >= q*(b-x)) goto 1
         d = p/q
         u = x+d
         if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
         goto 2
      end if
1     if (x >= xm) then
         e = a-x
      else
         e = b-x
      end if
      d = cgold*e
2     if (abs(d) >= tol1) then
         u = x+d
      else
         u = x+sign(tol1,d)
      end if
      fu = F(u)
      if (fu <= xthr) then    ! exit if f <= xthr
          x = u
          fx = fu
!          write(*,'(a,f12.5)') 'max exit 2', -fx
          goto 3
      end if
      if (fu <= fx) then
         if (u >= x) then
            a = x
         else
            b = x
         end if
         v = w
         fv = fw
         w = x
         fw = fx
         x = u
         fx = fu
      else
         if (u < x) then
            a = u
         else
            b = u
         end if
         if (fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
         else if (fu <= fv .or. v == x .or. v == w) then
            v = u
            fv = fu
          end if
      end if
!    write(*,'(a,i5,3f15.8)') 'iter, fx', iter, x, fx, fu
   end do
   call stop('BrentMod','exceeded maximum iterations',6)
3   xmin = x
   BrentMod = fx

end function BrentMod

!************************************************************************
!*                                                                      *
!*     KnuthShuffle                                                     *
!*                                                                      *
!************************************************************************

! ... Shuffle 1-Dimensional List if Integers   !Pascal Hebbeker
!     "The Art of Computer Programming" Second Edition Donald E. Knuth 1981
!     modified from http://rosettacode.org/wiki/Knuth_shuffle#Fortran


subroutine KnuthShuffle(a,asize, iseed)

   implicit none

   integer(4), intent(in) :: asize !size needed to cope with allocatable arrays
   integer(4), intent(inout) :: a(asize)
   integer(4), intent(inout) :: iseed   ! seed to random number generator
   integer :: i, randpos, itmp
   real(8) :: Random

   do i = size(a), 2, -1
      randpos = int(Random(iseed) * i) + 1
      itmp = a(randpos)
      a(randpos) = a(i)
      a(i) = itmp
   end do

end subroutine KnuthShuffle
