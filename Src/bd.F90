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
!> \page bd bd.F90
!! **BDDriver**
!! *Brownian dynamics driver*
!************************************************************************

!> \page nmlBD
!! The namelist \ref nmlBD contains variables that control the BD simulation. The simulation is performed in the configurational space according to Ermark, 1975.
!! * Variables:
!!  * \subpage nmlBDtstep
!!  * \subpage dcoeff
subroutine BDDriver(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='BDDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IOBD(iStage)

   case (iWriteInput)

      call IOBD(iStage)

   case (iSimulationStep)

      if (ltime) call CpuAdd('start', txroutine, 0, uout)
      call BDStep
      call UTotal(iStage)
      if (itest == 1) call TestSimulation
      if (ltime) call CpuAdd('stop', txroutine,  0, uout)

   end select

end subroutine BDDriver

!************************************************************************
!> \page bd bd.F90
!! **IOBD**
!! *perform i/o on Brownian dynamics variables*
!************************************************************************


subroutine IOBD(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOBD'
   character(80), parameter :: txheading ='brownian dynamcis data'

   namelist /nmlBD/  tstep, dcoeff

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(dcoeff)) then
         allocate(dcoeff(npt))
         dcoeff = 0.0E+00
      end if

      rewind(uin)
      read(uin,nmlBD)

   case (iWriteInput)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,4(f10.3,2x))') 'time step                   = ', tstep
         write(uout,'(a,t35,4(f10.3,2x))') 'self diffusion coefficient  = ', dcoeff(1:npt)
         write(uout,'(a,t35,4(f10.3,2x))') 'friction coefficient        = ',        &
          GasConstant*temp*scltem/(masspt(1:npt)*sclmas*dcoeff(1:npt)*scldif)*scltim
         write(uout,'(a,t35,4(f10.3,2x))') 'friction coeff. * time step = ',        &
          tstep*scltim*GasConstant*temp*scltem/(masspt(1:npt)*sclmas*dcoeff(1:npt)*scldif)
      end if

   end select

end subroutine IOBD

!************************************************************************
!> \page bd bd.F90
!! **BDStep**
!! *perform one bd step (configuration space)*
!************************************************************************


subroutine BDStep

   use MDModule
   implicit none

   character(40), parameter :: txroutine ='BDStep'

   real(8), allocatable, save :: fac1(:), fac2(:)
   logical, save :: first = .true.
   integer(4) :: ip, ipt
   real(8)    :: GauRandom

   if (ltrace) call WriteTrace(2, txroutine, iSimulationStep)

   if (first) then
      allocate(fac1(npt), fac2(npt))
      fac1 = 0.0E+00
      fac2 = 0.0E+00
      fac1(1:npt) = tstep*scltim*dcoeff(1:npt)*scldif*sclfor/(GasConstant*temp*scltem)/scllen
      fac2(1:npt) = sqrt(Two*dcoeff(1:npt)*scldif*tstep*scltim)/scllen
      first = .false.
   end if

   do ip = 1, np
      ipt = iptpn(ip)
      drostep(1,ip) = fac1(ipt)*forceo(1,ip) + fac2(ipt)*GauRandom(iseed)
      drostep(2,ip) = fac1(ipt)*forceo(2,ip) + fac2(ipt)*GauRandom(iseed)
      drostep(3,ip) = fac1(ipt)*forceo(3,ip) + fac2(ipt)*GauRandom(iseed)
      ro(1,ip) = ro(1,ip) + drostep(1,ip)
      ro(2,ip) = ro(2,ip) + drostep(2,ip)
      ro(3,ip) = ro(3,ip) + drostep(3,ip)
      call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
   end do
   call SetAtomPos(1, np, .false.)

end subroutine BDStep

