! ... 'version 6.4.7, Sep 18, 2015',

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
!*     MixedDriver                                                      *
!*                                                                      *
!************************************************************************

! ... mixed task driver

subroutine MixedDriver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='MixedDriver'

   integer(4), parameter :: mnmode = 6
   integer(4), parameter :: mnchoice = 10
   character(70), save :: txmixed(mnmode) = [ 'generate a potential energy curve            ',   &
                                              'calculate potential energies on a lattice    ',   &
                                              'calculate the global potential energy minimum',   &
                                              'generate random coordinates                  ',   &
                                              'calculate the second virial coefficient      ',   &
                                              'calculate orient. averaged potential energy  ' ]
   character(10) :: txpt1, txpt2
   integer(4), save ::  mode(mnchoice), nchoice, iseedst, ip
   real(8), save :: coord1(12), coord2(12)
   integer(4) :: ipt, iopt

   namelist /nmlMixed/ mode, txpt1, txpt2, coord1, coord2, ip

#if defined (_PAR_)
   call Stop(txroutine, '_PAR_ not supported', uout)
#endif

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

! ...............  read input data   ................

      mode = 0

      rewind(uin)
      read(uin,nmlMixed)

      nchoice = count( (mode >= 1) .and. (mode <= mnmode) )   ! determine nchoice
      iseedst = iseed                                         ! start different tasks with the same seed

! ... change to two particles

      nppt = 0
      do ipt = 1, npt
         if (txpt1 == txpt(ipt)) nppt(ipt) = nppt(ipt)+1
         if (txpt2 == txpt(ipt)) nppt(ipt) = nppt(ipt)+1
      end do
      call SetObjectParam1
      if (np /= 2) call Stop(txroutine, 'np /= 2', uout)

! ... set box parameters appropriate to the different tasks

       txbc = 'sph'
       sphrad = 10.0
       call SetBoxParam

   case (iWriteInput)

! ...............  write input data   ................

      call WriteHead(1, 'mixed data', uout)
      write(uout,'(a,a)') ('mode : ', txmixed(mode(iopt)), iopt = 1, nchoice)
      write(uout,'()')
      write(uout,'(a,t35,6i10)')   'particle to be moved           = ', ip
      call CopyCoord(np, coord1, coord2, ro, ori)
      call OrthoOri(np, 1, np, ori, 1.0d-4, uout)
      call SetAtomProp(1, np, .false.)
      call WritePartAtomProp
      call CpuTot(uout)

! ...............   perform the mixed tasks  ................

      do iopt = 1, nchoice
         call CopyCoord(np, coord1, coord2, ro, ori)
         call OrthoOri(np, 1, np, ori, 1.0d-4, uout)
         call SetAtomProp(1, np, .false.)
         iseed = iseedst
         if (mode(iopt) == 1) call Mixed1(ip)
         if (mode(iopt) == 2) call Mixed2(ip)
         if (mode(iopt) == 3) call Mixed3(ip)
         if (mode(iopt) == 4) call Mixed4(ip)
         if (mode(iopt) == 5) call Mixed5(ip)
         if (mode(iopt) == 6) call Mixed6(ip)
         call CpuTot(uout)
      end do

   case (iBeforeSimulation)

      call ImageDriver(iAfterSimulation)
      stop

   end select

end subroutine MixedDriver

!************************************************************************
!*                                                                      *
!*     Mixed1                                                           *
!*                                                                      *
!************************************************************************

! ... generate a potential energy curve

subroutine Mixed1(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='potential energy curve'
   integer(4) :: iaxis, mstep, m
   real(8)    :: dxx, dyy, dzz, dr, umax, dum(3), utwob

   namelist /nmlMixed1/ dxx, dyy, dzz, iaxis, dr, mstep, umax

   !rewind(uin)          ! rewind does not yet work
   read(uin,nmlMixed1)
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6f10.3)') 'translation steps              = ', dxx, dyy, dzz
   write(uout,'(a,6i10)')   'rotation axis                  = ', iaxis
   write(uout,'(a,6f10.3)') 'rotation angle                 = ', dr
   write(uout,'(a,6i10)')   'maximum numbe of steps         = ', mstep
   write(uout,'(a,6f10.3)') 'upper energy boundary          = ', umax

   dr = dr*sclang

   call WriteHead(2, 'results', uout)
   write(uout,'(3x,3(5x,a,6x),5x,a,10x,a)') 'x', 'y', 'z', 'dr', 'potential energy'
   write(uout,'(3x,3(5x,a,6x),5x,a,10x,a)') '-', '-', '-', '--', '----------------'
   call SetAtomProp(1, 2, .false.)
   call UTwoBodyPair(1, 2, utwob, dum)
!   utwob = utwob/AuTokJ/Bohr                                 !!! temporary
   write(uout,'(4f12.3,f20.4)') ro(1:3,ip), Zero, utwob
   if (ilist/= 0) then
      write(ulist,*) 'potential energy curve'
      write(ulist,'(4(f12.3,a),f20.4)') ro(1,ip), tab, ro(2,ip), tab, ro(3,ip), tab, Zero, tab, utwob
   end if
   do m = 1, mstep
      call NewCoord(ip, dxx, dyy, dzz, iaxis, dr)
      call UTwoBodyPair(1, 2, utwob, dum)
!      utwob = utwob/AuTokJ/Bohr                                 !!! temporary
      if (utwob > umax) return
      write(uout,'(4f12.3,f20.4)') ro(1:3,ip), m*dr/sclang, utwob
      if (ilist/= 0) write(ulist,'(4(f12.3,a),2f20.4)') ro(1,ip), tab, ro(2,ip), tab, ro(3,ip), tab, m*dr/sclang, tab, utwob
   end do

end subroutine Mixed1

!************************************************************************
!*                                                                      *
!*     Mixed2                                                           *
!*                                                                      *
!************************************************************************

! ... calculate potential energies on a 2d lattice

subroutine Mixed2(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='potential energies on a 2d lattice'
   integer(4) :: itmax, iaxis, nrdiv1, nrdiv2, m1, m2, it, imove
   real(8)    :: udel, dr, rdist, rlow1, rupp1, rlow2, rupp2
   real(8)    :: dr1, dr2, ulim, utwob, dum(3)
   logical    :: lmove

   namelist /nmlMixed2/ itmax, udel, dr, iaxis, rdist, rlow1, rupp1, nrdiv1, rlow2, rupp2, nrdiv2

   rewind(uin)
   read(uin,nmlMixed2)
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6i10)')      'maximum number of iterations   = ', itmax
   write(uout,'(a,6f10.3)')    'potential energy tolerance     = ', udel
   write(uout,'(a,6f10.3)')    'rotation angle                 = ', dr
   write(uout,'(a,i10,f10.3)') 'coordinate along axis no       = ', iaxis, rdist
   write(uout,'(a,2f10.3,i10)')'rlow1, rupp1, nrdiv1           = ', rlow1, rupp1, nrdiv1
   write(uout,'(a,2f10.3,i10)')'rlow2, rupp2, nrdiv2           = ', rlow2, rupp2, nrdiv2

   if (iaxis == 1) ro(1,ip) = rdist
   if (iaxis == 2) ro(2,ip) = rdist
   if (iaxis == 3) ro(3,ip) = rdist
   dr1 = (rupp1-rlow1)/nrdiv1
   dr2 = (rupp2-rlow2)/nrdiv2
   imove = 0

   dr = dr*sclang

   call WriteHead(2, 'results', uout)
   write(uout,'(3x,3(5x,a,6x),5x,a,5x,a)') 'x', 'y', 'z', 'potential energy'
   write(uout,'(3x,3(5x,a,6x),5x,a,5x,a)') '-', '-', '-', '----------------'
   do m1 = 0, nrdiv1
      do m2 = 0, nrdiv2
         if (iaxis == 1) then
            ro(2,ip) = rlow1+m1*dr1
            ro(3,ip) = rlow2+m2*dr2
         else if (iaxis == 2) then
            ro(3,ip) = rlow1+m1*dr1
            ro(1,ip) = rlow2+m2*dr2
         else if (iaxis == 3) then
            ro(1,ip) = rlow1+m1*dr1
            ro(2,ip) = rlow2+m2*dr2
         end if
         call SetAtomProp(ip, ip, .false.)
         call UTwoBodyPair(1, 2, utwob, dum)
         do it = 1, itmax
            lmove = .false.
            do imove = 1, 3
               call LineMove(ip, imove, Zero, Zero, Zero, dr, udel, utwob, lmove)
            end do
            if (.not.lmove) exit
         end do
         ulim = 10000.0
         write(uout,'(4f12.3,f20.4)') ro(1:3,ip), min(utwob,ulim)
         if (lmove) write(uout,'(a)') 'sorry, this point was not converged'
         if (ilist/= 0) write(ulist,'(4(f12.3,a),f20.4)') ro(1,ip), tab, ro(2,ip), tab, ro(3,ip), tab, min(utwob,ulim)
      end do
   end do

end subroutine Mixed2

!************************************************************************
!*                                                                      *
!*     Mixed3                                                          *
!*                                                                      *
!************************************************************************

! ... calculate the global potential energy minimum

subroutine Mixed3(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='global potential energy minimum'
   integer(4) :: itmax, imove, it, m
   real(8)    :: udel, dxx, dyy, dzz, dr, utwob, dum(3), ufac
   logical    :: lmove

   namelist /nmlMixed3/ itmax, udel, dxx, dyy, dzz, dr

   ufac = one
!   ufac = one/(AuTokJ*Bohr)     !!! temporary

   rewind(uin)
   read(uin,nmlMixed3)
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6i10)')   'maximum number of iterations   = ', itmax
   write(uout,'(a,6f10.3)') 'potential energy tolerance     = ', udel
   write(uout,'(a,6f10.3)') 'translation steps              = ', dxx, dyy, dzz
   write(uout,'(a,6f10.3)') 'rotation angle                 = ', dr

   dr = dr*sclang

   call WriteHead(2, 'results', uout)
   call UTwoBodyPair(1, 2, utwob, dum)
   write(uout,'(a,t15,a,t31,a,t57,a,t80,a,t102,a)')                     &
   'iter. no', 'pot. energy', 'center of mass', 'x''-axis', 'y''-axis', 'h''-axis'
   write(uout,'(a,t15,a,t31,a,t57,a,t80,a,t102,a)')                     &
   '---------', '----------', '--------------', '-------', '-------', '-------'
   write(uout,'(t31,3(a,6x),t53,3(a,6x),t76,3(a,6x),t99,3(a,6x))') ('x','y','z',m = 1,4)
   it = 0
   write(uout,'(i4,t13,f12.4,t26,3f7.3,t48,3f7.3,t71,3f7.3,t94,3f7.3)') it, utwob*ufac, ro(1:3,ip), ori(1:3,1:3,ip)
   do it = 1, itmax
      lmove = .false.
      do imove = 1, 6
         call LineMove(ip, imove, dxx, dyy, dzz, dr, udel, utwob, lmove)
      end do
      write(uout,'(i4,t13,f12.4,t26,3f7.3,t48,3f7.3,t71,3f7.3,t94,3f7.3)') it, utwob*ufac, ro(1:3,ip), ori(1:3,1:3,ip)
      if (.not.lmove) exit
   end do
   if (lmove) write(uout,'(a)') 'sorry, no convergence'
   call SetAtomProp(1, 2, .false.)
   call WritePartAtomProp

end subroutine Mixed3

!************************************************************************
!*                                                                      *
!*     Mixed4                                                           *
!*                                                                      *
!************************************************************************

! ... generate random coordinates

subroutine Mixed4(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='random coordinates'
   integer(4) :: iwr, no, ntry, itry, ino, m, m2
   real(8)    :: rlow, rupp, thlow, thupp, filow, fiupp, ulow, uupp, dum(3), utwob

   namelist /nmlMixed4/ iwr, rlow, rupp, thlow, thupp, filow, fiupp, ulow, uupp, no, ntry

   rewind(uin)
   read(uin,nmlMixed4)
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6i10)')   'weight exponent in r           = ', iwr
   write(uout,'(a,6f10.3)') 'lower and upper radial dist.   = ', rlow, rupp
   write(uout,'(a,6f10.3)') 'lower and upper theta angel    = ', thlow, thupp
   write(uout,'(a,6f10.3)') 'lower and upper fi angel       = ', filow, fiupp
   write(uout,'(a,6f10.3)') 'lower and upper energy level   = ', ulow, uupp
   write(uout,'(a,6i10)')   'number of desired points       = ', no
   write(uout,'(a,6i10)')   'maximum number of attempts     = ', ntry

   thlow = thlow*sclang
   thupp = thupp*sclang
   filow = filow*sclang
   fiupp = fiupp*sclang

   call WriteHead(2, 'results', uout)
   write(uout,'(a,t15,a,t31,a,t57,a,t80,a,t102,a)')                     &
   'iter. no', 'pot. energy', 'center of mass', 'x''-axis', 'y''-axis', 'h''-axis'
   write(uout,'(a,t15,a,t31,a,t57,a,t80,a,t102,a)')                     &
   '---------', '----------', '--------------', '-------', '-------', '-------'
   write(uout,'(t31,3(a,6x),t53,3(a,6x),t76,3(a,6x),t99,3(a,6x))')      &
   ('x','y','z',m = 1,4)
   ino = 0
   do itry = 1, ntry
      call RandomPos(iseed, iwr, rlow, rupp, thlow, thupp, filow, fiupp, ro(1:3,ip))
      call SetPartOriRandom(iseed, ori(1,1,ip))
      call SetAtomProp(ip, ip, .false.)
      call UTwoBodyPair(1, 2, utwob, dum)
      if (utwob > ulow .and. utwob < uupp) then
         ino = ino+1
         write(uout,'(i4,t13,f12.4,t26,3f7.3,t48,3f7.3,t71,3f7.3,t94,3f7.3)' ) ino, utwob, ro(1:3,ip), ori(1:3,1:3,ip)
         if (ilist/= 0)                                                           &
            write(ulist,'((i4,a),(f12.4,a),3(f7.3,a),9(f7.3,a))') &
            ino, tab, utwob, (tab,ro(m,ip),m = 1,3), ((tab,ori(m,m2,ip), m = 1,3), m2 = 1,3)
         if (ino == no) exit
      end if
   end do
   if (ino < no) write(uout,'(a)') 'sorry, no convergence'
   write(uout,'()')
   write(uout,'(a,6i10)') 'number of trials               = ', itry

end subroutine Mixed4

!************************************************************************
!*                                                                      *
!*     Mixed5                                                           *
!*                                                                      *
!************************************************************************

! ... calculate the second virial coefficient

subroutine Mixed5(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='second virial coefficient'
   integer(4) :: nblock
   integer(4) :: iwr, no, m
   real(8)    :: rlow, rupp, thlow, thupp, filow, fiupp, dum(3), utwob
   real(8)    :: sum, sum1, sum2, term, r1, rfac, fac

   namelist /nmlMixed5/ iwr, rlow, rupp, thlow, thupp, filow, fiupp, no, nblock

   nblock = 10

   rewind(uin)
   read(uin,nmlMixed5)
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6i10)')   'weight exponent in r           = ', iwr
   write(uout,'(a,6f10.3)') 'lower and upper radial dist.   = ', rlow, rupp
   write(uout,'(a,6f10.3)') 'lower and upper theta angel    = ', thlow, thupp
   write(uout,'(a,6f10.3)') 'lower and upper fi angel       = ', filow, fiupp
   write(uout,'(a,6i10)')   'number of desired points       = ', no

   thlow = thlow*sclang
   thupp = thupp*sclang
   filow = filow*sclang
   fiupp = fiupp*sclang

   call WriteHead(2, 'results', uout)
   sum1 = Zero
   sum2 = Zero
   do m = 1, no
      if (mod(m-1,no/nblock) == 0) sum = Zero
      call RandomPos(iseed, iwr, rlow, rupp, thlow, thupp, filow, fiupp, ro(1:3,ip))
      call SetPartOriRandom(iseed, ori(1,1,ip))
      call SetAtomProp(ip, ip, .false.)
      call UTwoBodyPair(1, 2, utwob, dum)
      rfac = One
      if (iwr/= 2) then
         r1 = sqrt((ro(1,1)-ro(1,2))**2+(ro(2,1)-ro(2,2))**2+(ro(3,1)-ro(3,2))**2)
         rfac = r1**(2-iwr)
      end if
!      write(54,'(a,f10.4,e15.5)') 'r1,utwob',r1,utwob  ! Wrong computation without this line using intel compiler!!!
!      rewind(54)                                         ! Just to limit the size of unit 54
      if (utwob*beta <-100) then
         term = Zero
      else if (utwob*beta > 100) then
         term =-rfac
      else
         term = (exp(-utwob*beta)-1.)*rfac
      end if
      sum = sum+term
      if (mod(m,no/nblock) == 0) then
         rfac = (rupp**(1+iwr)-rlow**(1+iwr))/(1+iwr)
         sum =-sum*2*Pi*rfac/(no/nblock)
         sum1 = sum1+sum
         sum2 = sum2+sum**2
      end if
   end do
   sum1 = sum1/nblock
   sum2 = sqrt(max(Zero,(sum2/nblock-sum1**2)/(nblock-1)))
   fac = AvNo*(scllen/1.0d-2)**3
   write(uout,'(a,f15.4)') 'second virial coeff.                        = ', sum1
   write(uout,'(a,f15.4)') 'sd of the second virial coeff.              = ', sum2
   write(uout,'(a,f15.4)') 'second virial coeff. (cm**3/mole)           = ', fac*sum1
   write(uout,'(a,f15.4)') 'sd of the second virial coeff. (cm**3/mole) = ', fac*sum2

end subroutine Mixed5

!************************************************************************
!*                                                                      *
!*     Mixed6                                                           *
!*                                                                      *
!************************************************************************

! ... calculate orientational averaged potential energy

subroutine Mixed6(ip)

   use MolModule
   implicit none

   integer(4), intent(in) :: ip

   character(80), parameter :: txheading ='orientational averaged potential energy'
   integer(4) :: no, nbin, ibin, iwr, ino
   real(8), allocatable :: sa(:)
   real(8)    :: bin, bini, rlow, rupp, thlow, thupp, filow, fiupp, ulow, uupp, utwob
   real(8)    :: uboltz, uave, uw, usum, w, norm, dum(3)

   namelist /nmlMixed6/ rlow, rupp, thlow, thupp, filow, fiupp, ulow, uupp, no, nbin

   rewind(uin)
   read(uin,nmlMixed6)
   allocate(sa(-1:nbin))
   sa = 0.0E+00
   call WriteHead(1, txheading, uout)
   write(uout,'(a,6f10.3)') 'lower and upper radial dist.   = ', rlow, rupp
   write(uout,'(a,6f10.3)') 'lower and upper theta angel    = ', thlow, thupp
   write(uout,'(a,6f10.3)') 'lower and upper fi angel       = ', filow, fiupp
   write(uout,'(a,6f10.3)') 'lower and upper energy level   = ', ulow, uupp
   write(uout,'(a,6i10)')   'number of desired points       = ', no
   write(uout,'(a,6i10)')   'number of desired grids        = ', nbin

   if (tempst > Zero) beta   = sclene/(GasConstant*tempst*scltem)
   thlow = thlow*sclang
   thupp = thupp*sclang
   filow = filow*sclang
   fiupp = fiupp*sclang

   call WriteHead(2, 'results', uout)
   bin = (uupp-ulow)/nbin
   bini = One/bin
   sa = Zero
   usum = Zero
   uboltz = Zero
   uw = Zero
   do ino = 1, no
      call RandomPos(iseed, iwr, rlow, rupp, thlow, thupp, filow, fiupp, ro(1:3,ip))
      call SetPartOriRandom(iseed, ori(1,1,ip))
      call SetAtomProp(ip, ip, .false.)
      call UTwoBodyPair(1, 2, utwob, dum)
      ibin = max(-1,min(floor(bini*(utwob-ulow)),int(nbin)))
      sa(ibin) = sa(ibin)+One
      usum = usum+utwob
! ... temporary fix for hard-core overlap of charged atoms
      if (utwob <-15.0d0) utwob = 10000.0d0
      if (utwob*beta < 80.0d0) then
         w = exp(-utwob*beta)
         uboltz = uboltz+utwob*w
         uw = uw+w
      end if
   end do
   norm = One/no
   usum = usum*norm
   uboltz = uboltz/uw
   sa = sa*norm
   uave = Zero
   do ibin = 0, nbin
      uave = uave+(ulow+(0.5+ibin)*bin)*sa(ibin)
   end do
   write(uout,'(a,f10.4)')  'potential energy, direct evaluation     = ', usum
   write(uout, '(a,f10.4)') 'potential energy, from distr. function  = ', uave
   write(uout,'(a,f10.4)')  'potential energy, boltzmann weighted    = ', uboltz
   call Plot('potential energy', nbin, sa(0), 'none', ulow, uupp, dum(1), dum(1), uout)
   if (ilist/= 0) then
      write(ulist,*) txheading
      write(ulist,'(2g15.5)') (ulow+(0.5+ibin)*bin,sa(ibin),ibin = 0,nbin-1)
   end if
   deallocate(sa)
end subroutine Mixed6

!************************************************************************
!*                                                                      *
!*     LineMove                                                         *
!*                                                                      *
!************************************************************************

! ... perform one movement, rotation or translation, to a new minimum

subroutine LineMove(ip, imove, dxx, dyy, dzz, dr, udel, utwob, lmove)

   use MolModule
   implicit none

   integer(4), intent(in)    :: ip            ! number of the particle to be moved
   integer(4), intent(in)    :: imove         != 1-3 rotation, 4-6 translation
                                              != 1, 4 x-axis, 2, 5 y-axis, 3, 6 z-axis
   real(8),    intent(in)    :: dxx, dyy, dzz ! translational displacement step
   real(8),    intent(in)    :: dr            ! rotational displacement step
   real(8),    intent(in)    :: udel          ! convergency level for pot energy
   real(8),    intent(inout) :: utwob         ! pot energy before and after move
   logical,    intent(out)   :: lmove         !.true. if the particle is moved

   real(8) :: roold(3), oriold(3,3), dx, dy, dz, ddr, dum(3), utwobold
   integer(4) :: i

   dx = Zero
   dy = Zero
   dz = Zero
   ddr = Zero
   if (imove >= 1 .and. imove <= 3) ddr = dr
   if (imove == 4) dx = dxx
   if (imove == 5) dy = dyy
   if (imove == 6) dz = dzz

   do i = 1, -1, -2
      do
         utwobold       = utwob
         roold(1:3)      = ro(1:3,ip)
         oriold(1:3,1:3) = ori(1:3,1:3,ip)
         call NewCoord(ip, i*dx, i*dy, i*dz, imove, i*ddr)
         call UTwoBodyPair(1, 2, utwob, dum)
!            write(*,'(a,2e15.5)') '- old',utwob, utwobold
         if (utwob-utwobold <-udel) then
            lmove = .true.
         else
            utwob           = utwobold
            ro(1:3,ip)      = roold(1:3)
            ori(1:3,1:3,ip) = oriold(1:3,1:3)
            exit
         end if
      end do
   end do

end subroutine LineMove

!************************************************************************
!*                                                                      *
!*     NewCoord                                                         *
!*                                                                      *
!************************************************************************

! ... translate and rotate a particle and update its atom coordinates

subroutine NewCoord(ip, dxx, dyy, dzz, iaxis, dr)

   use MolModule
   implicit none

   integer(4), intent(in)  :: ip
   real(8),    intent(in)  :: dxx, dyy, dzz
   integer(4), intent(in)  :: iaxis
   real(8),    intent(in)  :: dr

   real(8) :: Pih

   Pih = Half*Pi
   ro(1,ip) = ro(1,ip)+dxx
   ro(2,ip) = ro(2,ip)+dyy
   ro(3,ip) = ro(3,ip)+dzz
   if (iaxis == 1) then
      call EulerRot('func/rot', 'rad', -Pih, dr, Pih, ori(1,1,ip), ori(2,1,ip), ori(3,1,ip))
      call EulerRot('func/rot', 'rad', -Pih, dr, Pih, ori(1,2,ip), ori(2,2,ip), ori(3,2,ip))
      call EulerRot('func/rot', 'rad', -Pih, dr, Pih, ori(1,3,ip), ori(2,3,ip), ori(3,3,ip))
   else if (iaxis == 2) then
      call EulerRot('func/rot', 'rad', Zero, dr, Zero, ori(1,1,ip), ori(2,1,ip), ori(3,1,ip))
      call EulerRot('func/rot', 'rad', Zero, dr, Zero, ori(1,2,ip), ori(2,2,ip), ori(3,2,ip))
      call EulerRot('func/rot', 'rad', Zero, dr, Zero, ori(1,3,ip), ori(2,3,ip), ori(3,3,ip))
   else if (iaxis == 3) then
      call EulerRot('func/rot', 'rad', dr, Zero, Zero, ori(1,1,ip), ori(2,1,ip), ori(3,1,ip))
      call EulerRot('func/rot', 'rad', dr, Zero, Zero, ori(1,2,ip), ori(2,2,ip), ori(3,2,ip))
      call EulerRot('func/rot', 'rad', dr, Zero, Zero, ori(1,3,ip), ori(2,3,ip), ori(3,3,ip))
   end if
   call SetAtomProp(ip, ip, .false.)

end subroutine NewCoord

!************************************************************************
!*                                                                      *
!*     RandomPos                                                        *
!*                                                                      *
!************************************************************************

! ... generate a random position

subroutine RandomPos(iseed, iwr, rlow, rupp, thlow, thupp, filow, fiupp, r)

   implicit none

   integer(4), intent(inout)  :: iseed         ! seed to random number generator
   integer(4), intent(in)  :: iwr           ! radial weigth
   real(8),    intent(in)  :: rlow, rupp    ! lower and upper radial distance
   real(8),    intent(in)  :: thlow, thupp  ! lower and upper theta angle
   real(8),    intent(in)  :: filow, fiupp  ! lower and upper phi angle
   real(8),    intent(out) :: r(3)          ! coordinate of random position

   real(8) :: r1, cth, sth, fi, random

!  ... position is delimited by rlow, rupp, thlow, thupp, filow, and fiupp

   r1 = rlow+Random(iseed)**(1.0d0/(1+iwr))*(rupp-rlow)
   cth = cos(thlow)+Random(iseed)*(cos(thupp)-cos(thlow))
   sth = sqrt(1.0d0-cth*cth)
   fi = filow+Random(iseed)*(fiupp-filow)
   r(1) = r1*sth*cos(fi)
   r(2) = r1*sth*sin(fi)
   r(3) = r1*cth

end subroutine RandomPos

!************************************************************************
!*                                                                      *
!*     CopyCoord                                                        *
!*                                                                      *
!************************************************************************

! ... copy coordinates

subroutine CopyCoord(np, coord1, coord2, ro, ori)
   implicit none
   integer(4), intent(in)  :: np
   real(8),    intent(in)  :: coord1(12), coord2(12)
   real(8),    intent(out) :: ro(3,np)
   real(8),    intent(out) :: ori(3,3,np)
   ro(1:3,1)    = coord1(1:3)
   ori(1:3,1,1) = coord1(4:6)
   ori(1:3,2,1) = coord1(7:9)
   ori(1:3,3,1) = coord1(10:12)
   ro(1:3,2)    = coord2(1:3)
   ori(1:3,1,2) = coord2(4:6)
   ori(1:3,2,2) = coord2(7:9)
   ori(1:3,3,2) = coord2(10:12)
end subroutine CopyCoord
