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
!> \page md md.F90
!! **MDModule**
!! *module for md*
!************************************************************************

!> \page nmlMD
!! The namelist \ref nmlMD contains variables that control the MD simulation.
!! * Variables:
!!  * \subpage integ
!!  * \subpage nmlMDtstep
!!  * \subpage tvvite
!!  * \subpage nvvite
!!  * \subpage lsetvel
!!  * \subpage lzeromom
!!  * \subpage tvscl
!!  * \subpage tlscl
!!  * \subpage compre
module MDModule

   use MolModule
!> \page integ
!! `character(6)`
!! * `verlver` Integration according to the velocity form of the Verlet algorithm.
!! * `gear3` Integration according to a third-order Gear algorithm.
!! * `gear4` Integration according to a fourth-order Gear algorithm.
   character(6)  :: integ
!> \page nmlMDtstep
!! `real`
!! * Time step of the MD integration.

!> \page nmlBDtstep
!! `real`
!! * Time step of the BD integration.
   real(8)       :: tstep
!> \page tvvite
!! `real`
!! **default:** `0.0`
!! * Factor that determines the initial quaternion velocity for the iteration of the quaternion velocity (only \ref integ='velver'; 0.0 is fine).
   real(8)       :: tvvite
!> \page nvvite
!! `integer`
!! **default:** 2
!! * Number of iterations for the quaternion velocities (only \ref integ='velver'; 2 is preferred).
   integer(4)    :: nvvite
!> \page lsetvel
!! `logical`
!! **default:** `.true.` (only \ref txstart = 'setconf' .or. 'zero'); `.false.`  (only \ref txstart = 'continue')
!! * `.true.`: Linear and angular velocities are set according to a Maxwell distribution using the temperature \ref temp.
!! * `.false.`: No set of linear and angular velocities.
   logical       :: lsetvel
!> \page lzeromom
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Modify linear and angular velocities (only possible) to get zero linear and angular moments while preserving the translational and rotational temperature (only \ref lsetvel=.true.).
!! * `.false.`: Nothing.
   logical       :: lzeromom
!> \page tvscl
!! `real`
!! **default:** `0.0`
!! * Time constant for the velocity scaling. The scaling is applied if \ref tvscl>0.0 and is then performed every time step according to
!!  velocity(new)=velocity(old)* sqrt(1+(tstep/\ref tvscl)*(\ref temp /t-1)), where t is the instantaneous temperature (Berendsen et al., 1984).
!!  \ref tvscl=tstep scales the velocities to give t=\ref temp.
   real(8)       :: tvscl
!> \page tlscl
!! `real`
!! **default:** `0.0`
!! * Time constant for the length scaling. The scaling is applied if \ref tlscl >0.0 and is then performed every time step according to
!!  length(new) = length(old)*(1+x)**(-1/3), x=(tstep/\ref tlscl)*\ref compre*(p-\ref prsr), where p is the instantaneously pressure (Berendsen et
!!  al., 1984).
   real(8)       :: tlscl
!> \page compre
!! `real`
!! * Compressibility used for the length scaling.
   real(8)       :: compre

   real(8)       :: linmom(3)              ! total linear moment
   real(8)       :: angmom(3)              ! total angular moment

   real(8), allocatable :: rodd(:,:)       ! particle acceleration
   real(8), allocatable :: roddo(:,:)      ! particle accelaration, old
   real(8), allocatable :: roddd(:,:)      ! time derivative of rodd
   real(8), allocatable :: quadd(:,:)      ! quaterion accleration
   real(8), allocatable :: quaddo(:,:)     ! quaterion acceleration, old
   real(8), allocatable :: quaddd(:,:)     ! time derivative of quadd
end module MDModule

!************************************************************************
!> \page md md.F90
!! **MDDriver**
!! *molecular dynamics driver*
!************************************************************************


subroutine MDDriver(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MDDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IOMD(iStage)
      allocate(rodd(3,np_alloc))
      rodd = 0.0E+00
      allocate(roddo(3,np_alloc))
      roddo = 0.0E+00
      allocate(roddd(3,np_alloc))
      roddd = 0.0E+00
      allocate(quadd(0:3,np_alloc))
      quadd = 0.0E+00
      allocate(quaddo(0:3,np_alloc))
      quaddo = 0.0E+00
      allocate(quaddd(0:3,np_alloc))
      quaddd = 0.0E+00
      rodd = Zero
      roddo = Zero
      roddd = Zero
      quadd = Zero
      quaddo = Zero
      quaddd = Zero

   case (iWriteInput)

      call IOMD(iStage)

   case (iSimulationStep)

      if (ltime) call CpuAdd('start', txroutine, 0, uout)
      if (integ == 'velver') then
         call VelVer(iStage)
      else if (integ(1:4) == 'gear') then
         call Gear(iStage)
      else
         call Stop(txroutine, 'unsupported value of integ', uout)
      end if
      call GetKinEnergy
      if (tvscl > Zero) call ScaleVel
      if (tlscl > Zero) call ScaleLength
      if (itest == 1) call TestSimulation
      if (ltime) call CpuAdd('stop', txroutine, 0, uout)

   case (iAfterSimulation)

      deallocate(rodd)
      deallocate(roddo)
      deallocate(roddd)
      deallocate(quadd)
      deallocate(quaddo)
      deallocate(quaddd)

   end select

   if (ltrace) call WriteTrace(1, trim(txroutine)//'_exit', iStage)

end subroutine MDDriver

!************************************************************************
!> \page md md.F90
!! **IOMD**
!! *perform i/o on molecular dynamics variables*
!************************************************************************


subroutine IOMD(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOMD'

   namelist /nmlMD/ integ, tstep, tvvite, nvvite, lsetvel, lzeromom, tvscl, tlscl, compre

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      tvvite   = Zero
      nvvite   = 2
      if ((txstart == 'setconf') .or. (txstart == 'zero')) lsetvel = .true.
      if (txstart == 'continue') lsetvel = .false.
      lzeromom =.false.
      tvscl    = Zero
      tlscl    = Zero

      rewind(uin)
      read(uin,nmlMD)

      call LowerCase(integ)

   case (iWriteInput)

      if (master) then
         call WriteHead(2, 'md data', uout)
         write(uout,'(a,t35,a)')     'integration method             = ', integ
         write(uout,'(a,t35,f10.3)') 'total simulation time          = ', tstep*nstep
         write(uout,'(a,t35,f10.3)') 'time step                      = ', tstep
         write(uout,'(a,t35,f10.3)') 'factor in iteration (velver)   = ', tvvite
         write(uout,'(a,t35,i10)')   'no of iterations (velver)      = ', nvvite
         write(uout,'(a,t35,l10)')   'set random velocities          = ', lsetvel
         write(uout,'(a,t35,l10)')   'zero total moments             = ', lzeromom
         write(uout,'(a,t35,f10.3)') 'vel. scaling time const.       = ', tvscl
         write(uout,'(a,t35,f10.3)') 'length scaling time const.     = ', tlscl
         write(uout,'(a,t35,f10.3)') 'used isotherm. compress.       = ', compre
         if ((tvscl > Zero) .and. (tvscl < tstep)) call Stop(txroutine, 'Zero < tvscl < tstep is not supported', uout)
      end if

   end select

end subroutine IOMD

!************************************************************************
!> \page md md.F90
!! **VelVer**
!! *perform one step according to the velocity form of the Verlet algorithm*
!************************************************************************


!     ref: swope et al. JCP 76, 637 (1982)
!     the combination of the velocity form of the Verlet algorithm and quaterions is new

subroutine VelVer(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='VelVer'
   real(8), allocatable, save  :: quado(:,:)
   real(8)    :: fac
   integer(4) :: iplow, ipupp, ip, m

   if (ltrace) call WriteTrace(2, txroutine, iSimulationStep)

   if(.not.allocated(quado)) then
      allocate(quado(0:3,np))
      quado = 0.0E+00
   end if

#if defined (_XXXPAR_)
   load = ceiling(np/float(nproc))    ! mpi_gather requires eqaual load (ipmyid(1:2) may provide unequal load)
   imyid(1) = 1 + load*myid
   imyid(2) = load*(myid+1)
   par_len = 3*(imyid(2)-imyid(1)+1)
   iplow = imyid(1)
   ipupp = min(imyid(2),np)           ! to assure that ipupp not exceed np
!  if(master) write(*,'(i3,a,4i5)') myid,'master: iplow, ipupp, par_len',iplow, ipupp, par_len
!  if(slave ) write(*,'(i3,a,4i5)') myid,'slave : iplow, ipupp, par_len',iplow, ipupp, par_len
#else
   iplow = 1
   ipupp = np
#endif

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

! ............... new positions and quaternions ...............

   fac = Half*tstep**2
   do ip = iplow, ipupp
      drostep(1,ip) = tstep*rod(1,ip) + fac*rodd(1,ip)
      drostep(2,ip) = tstep*rod(2,ip) + fac*rodd(2,ip)
      drostep(3,ip) = tstep*rod(3,ip) + fac*rodd(3,ip)
      ro(1,ip)  = ro(1,ip) + drostep(1,ip)
      ro(2,ip)  = ro(2,ip) + drostep(2,ip)
      ro(3,ip)  = ro(3,ip) + drostep(3,ip)
      call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
   end do

   if (lpolyatom) then
      do ip = iplow, ipupp
         qua(0,ip) = qua(0,ip) + tstep*quad(0,ip) + fac*quadd(0,ip)
         qua(1,ip) = qua(1,ip) + tstep*quad(1,ip) + fac*quadd(1,ip)
         qua(2,ip) = qua(2,ip) + tstep*quad(2,ip) + fac*quadd(2,ip)
         qua(3,ip) = qua(3,ip) + tstep*quad(3,ip) + fac*quadd(3,ip)
      end do
      call QuaNorm(np, iplow, ipupp, qua)
      call QuaToOri(np, iplow, ipupp, qua, ori)
   end if

   call SetAtomProp(iplow, ipupp, lintsite)

#if defined (_XXXPAR_)
!   if (master) write(uout,'(a, 48f8.3)') 'master, before: ro  ',ro(1,1:np)
!   if (master) write(uout,'(a, 48f8.3)') 'master, before: r   ',r(1,1:np)
!   if (master) write(uout,'(a, 48f8.3)') 'master, before: dr  ',drostep(1,1:np)
!   if (slave ) write(uout,'(a, 48f8.3)') 'slave, before: ro   ',ro(1,1:np)
!   if (slave ) write(uout,'(a, 48f8.3)') 'slave, before: r    ',r(1,1:np)
   call par_gather_reals(ro(1,iplow),par_len,ro(1,1),par_len, vaux)
   call par_gather_reals(r(1,iplow),par_len,r(1,1),par_len, vaux)
   call par_gather_reals(drostep(1,iplow),par_len,drostep(1,1),par_len, vaux)
   call par_bc_reals(ro,3*np)
   call par_bc_reals(r,3*np)
   call par_bc_reals(drostep,3*np)
   if (lpolyatom) then
      call par_gather_reals(ori(1,1,iplow),3*par_len,ori(1,1,1),3*par_len, vaux)
      call par_bc_reals(ori,9*np)
   end if
!   if (master) write(uout,'(a, 48f8.3)') 'master, after: ro   ',ro(1,1:np)
!   if (master) write(uout,'(a, 48f8.3)') 'master, after: r    ',r(1,1:np)
!   if (master) write(uout,'(a, 48f8.3)') 'master, after: dr   ',drostep(1,1:np)
!   if (slave ) write(uout,'(a, 48f8.3)') 'slave, after: ro    ',ro(1,1:np)
!   if (slave ) write(uout,'(a, 48f8.3)') 'slave, after: r     ',r(1,1:np)
#endif

! ............... energy and force evaulation ...............

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   call UTotal(iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

! ............... new linear and quaternion acceleration and velocity ...............

! ... linear

   do ip = iplow, ipupp
      roddo(1,ip) = rodd(1,ip)
      roddo(2,ip) = rodd(2,ip)
      roddo(3,ip) = rodd(3,ip)
   end do
   call GetLinAcc(iplow, ipupp)
   do ip = iplow, ipupp
      rod(1,ip) = rod(1,ip) + Half*tstep*(rodd(1,ip)+roddo(1,ip))
      rod(2,ip) = rod(2,ip) + Half*tstep*(rodd(2,ip)+roddo(2,ip))
      rod(3,ip) = rod(3,ip) + Half*tstep*(rodd(3,ip)+roddo(3,ip))
   end do
#if defined (_XXXPAR_)
   call par_gather_reals(rod(1,iplow),par_len,rod(1,1),par_len, vaux)
   call par_bc_reals(rod,3*np)
#endif

! ... quaternions

   if (lpolyatom) then
      do ip = iplow, ipupp
         quaddo(0,ip) = quadd(0,ip)
         quaddo(1,ip) = quadd(1,ip)
         quaddo(2,ip) = quadd(2,ip)
         quaddo(3,ip) = quadd(3,ip)
      end do

! ... save old and set trial quaternion velocities

      do ip = iplow, ipupp
         quado(0,ip) = quad(0,ip)
         quado(1,ip) = quad(1,ip)
         quado(2,ip) = quad(2,ip)
         quado(3,ip) = quad(3,ip)
         quad(0,ip)  = quado(0,ip) + (tvvite*tstep)*quaddo(0,ip)
         quad(1,ip)  = quado(1,ip) + (tvvite*tstep)*quaddo(1,ip)
         quad(2,ip)  = quado(2,ip) + (tvvite*tstep)*quaddo(2,ip)
         quad(3,ip)  = quado(3,ip) + (tvvite*tstep)*quaddo(3,ip)
      end do

      do m = 1, nvvite
         call QuaVelToAngVel(np, iplow, ipupp, qua, quad, angvelo)
         call GetAngAcc(iplow, ipupp)
         do ip = iplow, ipupp
            quad(0,ip) = quado(0,ip) + (Half*tstep)*(quadd(0,ip)+quaddo(0,ip))
            quad(1,ip) = quado(1,ip) + (Half*tstep)*(quadd(1,ip)+quaddo(1,ip))
            quad(2,ip) = quado(2,ip) + (Half*tstep)*(quadd(2,ip)+quaddo(2,ip))
            quad(3,ip) = quado(3,ip) + (Half*tstep)*(quadd(3,ip)+quaddo(3,ip))
         end do
      end do
#if defined (_XXXPAR_)
   call par_gather_reals(angvelo(1,iplow),par_len,angvelo(1,1),par_len, vaux)
   call par_bc_reals(angvelo,3*np)
#endif

   end if

! ... set back to geometric sites

   if (lintsite) call SetAtomPos(iplow, ipupp, .false.)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine VelVer

!************************************************************************
!> \page md md.F90
!! **Gear**
!! *perform one step according to the Gear algorithm*
!************************************************************************


!     ref: Gear "numerical initial value problems in ordinary differential equations"
!     prentice-hall, englewood cliffs, n.j. 1971 and sonnenschein, j. comp. phys. 59, 347 (1985)

subroutine Gear(iStage)

   use MDModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Gear'
   real(8), save :: tay(1:3), gearc(0:3)
   real(8) :: fac, drodd(3), dquadd(0:3)
   integer(4) :: ip, m, iplow, ipupp

   iplow = 1
   ipupp = np

   if (ltrace) call WriteTrace(2, txroutine, iSimulationStep)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

! ... set taylor and gear coefficients  (need only to be done once)

   if (integ(5:5) == '3') then
      tay(1) = tstep
      tay(2) = Half*tay(1)*tstep
      tay(3) = Zero
      gearc(0) = Zero
      gearc(1) = One
      gearc(2) = One
      gearc(3) = Zero
   else if (integ(5:5) == '4') then
      tay(1) = tstep
      tay(2) = Half*tay(1)*tstep
      tay(3) = Third*tay(2)*tstep
      gearc(0) = Sixth
      gearc(1) = Five*Sixth
      gearc(2) = One
      gearc(3) = Third
   else
      call Stop(txroutine, 'unsupported order of gear', uout)
   end if

   fac = Half*tstep**2
   do m = 0, 3
      gearc(m) = fac*gearc(m)
      fac = fac*real(m+1)/tstep
   end do

! ..............   predict   ............

! ... linear motion

   do ip = iplow, ipupp
      drostep(1,ip) = tay(1)*rod(1,ip) + tay(2)*rodd(1,ip) + tay(3)*roddd(1,ip)
      drostep(2,ip) = tay(1)*rod(2,ip) + tay(2)*rodd(2,ip) + tay(3)*roddd(2,ip)
      drostep(3,ip) = tay(1)*rod(3,ip) + tay(2)*rodd(3,ip) + tay(3)*roddd(3,ip)
      ro(1,ip) = ro(1,ip) + drostep(1,ip)
      ro(2,ip) = ro(2,ip) + drostep(2,ip)
      ro(3,ip) = ro(3,ip) + drostep(3,ip)
      rod(1,ip) = rod(1,ip) + tay(1)*rodd(1,ip) + tay(2)*roddd(1,ip)
      rod(2,ip) = rod(2,ip) + tay(1)*rodd(2,ip) + tay(2)*roddd(2,ip)
      rod(3,ip) = rod(3,ip) + tay(1)*rodd(3,ip) + tay(2)*roddd(3,ip)
      rodd(1,ip) = rodd(1,ip) + tay(1)*roddd(1,ip)
      rodd(2,ip) = rodd(2,ip) + tay(1)*roddd(2,ip)
      rodd(3,ip) = rodd(3,ip) + tay(1)*roddd(3,ip)
      call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
   end do

! ... angular motion

   if (lpolyatom) then
      do ip = iplow, ipupp
         qua(0,ip)  = qua(0,ip)  +tay(1)*quad(0,ip) +tay(2)*quadd(0,ip)+tay(3)*quaddd(0,ip)
         qua(1,ip)  = qua(1,ip)  +tay(1)*quad(1,ip) +tay(2)*quadd(1,ip)+tay(3)*quaddd(1,ip)
         qua(2,ip)  = qua(2,ip)  +tay(1)*quad(2,ip) +tay(2)*quadd(2,ip)+tay(3)*quaddd(2,ip)
         qua(3,ip)  = qua(3,ip)  +tay(1)*quad(3,ip) +tay(2)*quadd(3,ip)+tay(3)*quaddd(3,ip)
         quad(0,ip) = quad(0,ip) +tay(1)*quadd(0,ip)+tay(2)*quaddd(0,ip)
         quad(1,ip) = quad(1,ip) +tay(1)*quadd(1,ip)+tay(2)*quaddd(1,ip)
         quad(2,ip) = quad(2,ip) +tay(1)*quadd(2,ip)+tay(2)*quaddd(2,ip)
         quad(3,ip) = quad(3,ip) +tay(1)*quadd(3,ip)+tay(2)*quaddd(3,ip)
         quadd(0,ip) = quadd(0,ip)+tay(1)*quaddd(0,ip)
         quadd(1,ip) = quadd(1,ip)+tay(1)*quaddd(1,ip)
         quadd(2,ip) = quadd(2,ip)+tay(1)*quaddd(2,ip)
         quadd(3,ip) = quadd(3,ip)+tay(1)*quaddd(3,ip)
      end do
      call QuaNorm(np, iplow, ipupp, qua)
      call QuaToOri(np, iplow, ipupp, qua, ori)
   end if

   call SetAtomProp(iplow, ipupp, lintsite)

! ............... energy and force evaulation ...............

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   call UTotal(iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

! ............   correct   .............

! ... linear motion

   do ip = iplow, ipupp
      roddo(1,ip) = rodd(1,ip)
      roddo(2,ip) = rodd(2,ip)
      roddo(3,ip) = rodd(3,ip)
   end do
   call GetLinAcc(iplow, ipupp)
   do ip = iplow, ipupp
      drodd(1)  = rodd(1,ip) - roddo(1,ip)
      drodd(2)  = rodd(2,ip) - roddo(2,ip)
      drodd(3)  = rodd(3,ip) - roddo(3,ip)
      ro(1,ip)   = ro(1,ip)   + gearc(0)*drodd(1)
      ro(2,ip)   = ro(2,ip)   + gearc(0)*drodd(2)
      ro(3,ip)   = ro(3,ip)   + gearc(0)*drodd(3)
      rod(1,ip)  = rod(1,ip)  + gearc(1)*drodd(1)
      rod(2,ip)  = rod(2,ip)  + gearc(1)*drodd(2)
      rod(3,ip)  = rod(3,ip)  + gearc(1)*drodd(3)
      rodd(1,ip) = roddo(1,ip)+ gearc(2)*drodd(1)
      rodd(2,ip) = roddo(2,ip)+ gearc(2)*drodd(2)
      rodd(3,ip) = roddo(3,ip)+ gearc(2)*drodd(3)
      roddd(1,ip) = roddd(1,ip)+ gearc(3)*drodd(1)
      roddd(2,ip) = roddd(2,ip)+ gearc(3)*drodd(2)
      roddd(3,ip) = roddd(3,ip)+ gearc(3)*drodd(3)
      call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
   end do

! ... angular motion

   if (lpolyatom) then
      do ip = iplow, ipupp
         quaddo(0,ip) = quadd(0,ip)
         quaddo(1,ip) = quadd(1,ip)
         quaddo(2,ip) = quadd(2,ip)
         quaddo(3,ip) = quadd(3,ip)
      end do
      call QuaVelToAngVel(np, iplow, ipupp, qua, quad, angvelo)
      call GetAngAcc(iplow, ipupp)
      do ip = iplow, ipupp
         dquadd(0)   = quadd(0,ip) - quaddo(0,ip)
         dquadd(1)   = quadd(1,ip) - quaddo(1,ip)
         dquadd(2)   = quadd(2,ip) - quaddo(2,ip)
         dquadd(3)   = quadd(3,ip) - quaddo(3,ip)
         qua(0,ip)   = qua(0,ip)   + gearc(0)*dquadd(0)
         qua(1,ip)   = qua(1,ip)   + gearc(0)*dquadd(1)
         qua(2,ip)   = qua(2,ip)   + gearc(0)*dquadd(2)
         qua(3,ip)   = qua(3,ip)   + gearc(0)*dquadd(3)
         quad(0,ip)  = quad(0,ip)  + gearc(1)*dquadd(0)
         quad(1,ip)  = quad(1,ip)  + gearc(1)*dquadd(1)
         quad(2,ip)  = quad(2,ip)  + gearc(1)*dquadd(2)
         quad(3,ip)  = quad(3,ip)  + gearc(1)*dquadd(3)
         quadd(0,ip) = quaddo(0,ip)+ gearc(2)*dquadd(0)
         quadd(1,ip) = quaddo(1,ip)+ gearc(2)*dquadd(1)
         quadd(2,ip) = quaddo(2,ip)+ gearc(2)*dquadd(2)
         quadd(3,ip) = quaddo(3,ip)+ gearc(2)*dquadd(3)
         quaddd(0,ip) = quaddd(0,ip)+ gearc(3)*dquadd(0)
         quaddd(1,ip) = quaddd(1,ip)+ gearc(3)*dquadd(1)
         quaddd(2,ip) = quaddd(2,ip)+ gearc(3)*dquadd(2)
         quaddd(3,ip) = quaddd(3,ip)+ gearc(3)*dquadd(3)
      end do
      call QuaNorm(np, iplow, ipupp, qua)
      call QuaToOri(np, iplow, ipupp, qua, ori)
   end if

   call SetAtomProp(iplow, ipupp, .false.)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine Gear

!************************************************************************
!> \page md md.F90
!! **GetLinAcc**
!! *calculate new linear acceleration*
!************************************************************************


!     uses forces(lab)

subroutine GetLinAcc(iplow, ipupp)

   use MDModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp

   integer(4)   :: ip

   do ip = iplow, ipupp
      rodd(1,ip) = massip(ip)*forceo(1,ip)
      rodd(2,ip) = massip(ip)*forceo(2,ip)
      rodd(3,ip) = massip(ip)*forceo(3,ip)
   end do

end subroutine GetLinAcc

!************************************************************************
!> \page md md.F90
!! **GetAngAcc**
!! *calculate new quaternion accelerations*
!************************************************************************


!     uses torques(lab), orientation, angular velocity, quaternion, and quaternion velocities

subroutine GetAngAcc(iplow, ipupp)

   use MDModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp

   integer(4) :: ip
   real(8)    :: tor1, tor2, tor3, angvelod1, angvelod2, angvelod3, fac

   do ip = iplow, ipupp

! ... transform torque from lab to principal axes

      tor1 = torqueo(1,ip)*ori(1,1,ip)+torqueo(2,ip)*ori(2,1,ip)+torqueo(3,ip)*ori(3,1,ip)
      tor2 = torqueo(1,ip)*ori(1,2,ip)+torqueo(2,ip)*ori(2,2,ip)+torqueo(3,ip)*ori(3,2,ip)
      tor3 = torqueo(1,ip)*ori(1,3,ip)+torqueo(2,ip)*ori(2,3,ip)+torqueo(3,ip)*ori(3,3,ip)

! ... calculate the angular acceleration (principal frame)

      angvelod1 = (tor1+angvelo(2,ip)*angvelo(3,ip)*(momp(2,ip)-momp(3,ip)))*momip(1,ip)
      angvelod2 = (tor2+angvelo(3,ip)*angvelo(1,ip)*(momp(3,ip)-momp(1,ip)))*momip(2,ip)
      angvelod3 = (tor3+angvelo(1,ip)*angvelo(2,ip)*(momp(1,ip)-momp(2,ip)))*momip(3,ip)

! ... calculate quaternion acceleration

      fac = -Two*(quad(2,ip)**2+quad(1,ip)**2+quad(3,ip)**2+quad(0,ip)**2)
      quadd(0,ip) = Half*(-qua(1,ip)*angvelod1+qua(2,ip)*angvelod2-qua(3,ip)*angvelod3+qua(0,ip)*fac)
      quadd(1,ip) = Half*(+qua(0,ip)*angvelod1-qua(3,ip)*angvelod2-qua(2,ip)*angvelod3+qua(1,ip)*fac)
      quadd(2,ip) = Half*(-qua(3,ip)*angvelod1-qua(0,ip)*angvelod2+qua(1,ip)*angvelod3+qua(2,ip)*fac)
      quadd(3,ip) = Half*(+qua(2,ip)*angvelod1+qua(1,ip)*angvelod2+qua(0,ip)*angvelod3+qua(3,ip)*fac)

   end do

end subroutine GetAngAcc

!************************************************************************
!> \page md md.F90
!! **ScaleLength**
!! *scale lengths through a pressure coupling to an external bath*
!************************************************************************


!     ref. Berendsen, JCP 81, 3684 (1984)

subroutine ScaleLength

   use MDModule
   implicit none

   character(40), parameter :: txroutine ='ScaleLength'
   integer(4) :: ip
   real(8)    :: fac, fac1

   if (lmd) then
      fac1 = tstep/tlscl
   else
      fac1 = 1.0d0
   end if
   fac = fac1*compre*(prsr-prsrst)

   if (fac <-One) call Stop(txroutine, 'fac <-1.0', uout)
   fac = (One+fac)**(0.3333333d0)
   if (lbcbox) then
      boxlen(1:3) = boxlen(1:3)*fac
   else if (lbcrd .or. lbcto) then
      cellside = cellside*fac
   else
      call Stop(txroutine,'tlscl > 0.0 and non-periodic boundary conditions',uout)
   endif
   call SetBoxParam
   if (lewald) call EwaldSetup
   do ip = 1, np
      ro(1,ip) = ro(1,ip) * fac
      ro(2,ip) = ro(2,ip) * fac
      ro(3,ip) = ro(3,ip) * fac
   end do
   call SetAtomPos(1, np, .false.)

end subroutine ScaleLength

!************************************************************************
!> \page md md.F90
!! **ScaleVel**
!! *scale velocities and calculate new kinetic energies*
!************************************************************************


subroutine ScaleVel
   call ScaleLinVel
   call ScaleAngVel
   call GetKinEnergy
end subroutine ScaleVel

!************************************************************************
!> \page md md.F90
!! **ScaleLinVel**
!! *scale linear velocities through a thermal coupling to an external bath*
!************************************************************************


!     ref: Berendsen, JCP 81, 3684 (1984).

subroutine ScaleLinVel

   use MDModule
   implicit none

   integer(4) :: ip
   real(8)    :: fac

   fac = Zero
   if (temptra > Zero) fac = sqrt(1+(tstep/tvscl)*(tempst/temptra-One))
   do ip = 1, np
      rod(1,ip) = rod(1,ip) * fac
      rod(2,ip) = rod(2,ip) * fac
      rod(3,ip) = rod(3,ip) * fac
   end do

end subroutine ScaleLinVel

!************************************************************************
!> \page md md.F90
!! **ScaleAngVel**
!! *scale angular velocities through a thermal coupling to an external bath*
!************************************************************************


!     ref: Berendsen, JCP 81, 3684 (1984).

subroutine ScaleAngVel

   use MDModule
   implicit none

   integer(4) :: ip
   real(8)    :: fac

   fac = Zero
   if (temprot > Zero) fac = sqrt(1+(tstep/tvscl)*(tempst/temprot-One))
   do ip = 1, np
      angvelo(1,ip) = angvelo(1,ip) * fac
      angvelo(2,ip) = angvelo(2,ip) * fac
      angvelo(3,ip) = angvelo(3,ip) * fac
   end do
   call AngVelToQuaVel(np, 1, np, qua, angvelo, quad)

end subroutine ScaleAngVel

!************************************************************************
!> \page md md.F90
!! **SetVel**
!! *set linear and angular velocities*
!************************************************************************


subroutine SetVel

   use MDModule
   implicit none

   call SetLinVel
   call SetAngVel

end subroutine SetVel

!************************************************************************
!> \page md md.F90
!! **SetLinVel**
!! *set linear velocities according to the maxwell distribution*
!************************************************************************


subroutine SetLinVel

   use MDModule
   implicit none

   integer(4) :: ip, ipt
   real(8)    :: fac, fact, GauRandom

   fac = GasConstant*tempst*scltem/(sclmas*sclvel**2)
   do ipt = 1, npt
      fact = sqrt(fac*massipt(ipt))
      do ip = ipnpt(ipt), ipnpt(ipt)+nppt(ipt)-1
         rod(1,ip) = fact * GauRandom(iseed)
         rod(2,ip) = fact * GauRandom(iseed)
         rod(3,ip) = fact * GauRandom(iseed)
      end do
   end do

end subroutine SetLinVel

!************************************************************************
!> \page md md.F90
!! **SetAngVel**
!! *set angular velocities according to the maxwell distribution*
!************************************************************************


subroutine SetAngVel

   use MDModule
   implicit none

   integer(4) :: ip, ipt
   real(8)    :: fac, facr(3), GauRandom

   fac = GasConstant*tempst*scltem/(sclmas*sclvel**2)
   do ipt = 1, npt
      facr(1:3) = sqrt(fac*momipt(1:3,ipt))
      do ip = ipnpt(ipt), ipnpt(ipt)+nppt(ipt)-1
         angvelo(1,ip) = facr(1) * GauRandom(iseed)
         angvelo(2,ip) = facr(2) * GauRandom(iseed)
         angvelo(3,ip) = facr(3) * GauRandom(iseed)
      end do
   end do
   call AngVelToQuaVel(np, 1, np, qua, angvelo, quad)

end subroutine SetAngVel

!************************************************************************
!> \page md md.F90
!! **SetZeroMom**
!! *set zero linear and angular moments*
!************************************************************************


subroutine SetZeroMom

   use MDModule
   implicit none

                  call SetZeroLinMom
   if (lpolyatom) call SetZeroAngMom

end subroutine SetZeroMom

!************************************************************************
!> \page md md.F90
!! **SetZeroLinMom**
!! *set zero linear moments*
!************************************************************************

!     if possible keep the linear moments consistent with temperature = tempst

subroutine SetZeroLinMom

   use MDModule
   implicit none

   integer(4), save :: mmax = 5
   integer(4) :: m
   real(8)    :: tollinmom(3), fac
   logical    :: ltest =.false.

   if (np < 2) return

! ... calculate tollinmom

   tollinmom(1) = 1.0d-4*sqrt(sum((massp(1:np)*rod(1,1:np))**2))
   tollinmom(2) = 1.0d-4*sqrt(sum((massp(1:np)*rod(2,1:np))**2))
   tollinmom(3) = 1.0d-4*sqrt(sum((massp(1:np)*rod(3,1:np))**2))

   if (ltest) write(uout,'(a,3e15.5)')       'tol  : linear moments         ', tollinmom

! ... calculate initial linear moments

   call GetLinMom

   if (ltest) call GetKinEnergy
   if (ltest) write(uout,'(a,3e15.5,f15.5)') 'init : linear moments and temp', linmom(1:3), temptra

   do m = 1, mmax

! ... adjust linear velocity so linear moments becomes Zero

      linmom(1:3) = linmom(1:3)/np
      rod(1,1:np) = rod(1,1:np)-linmom(1)*massip(1:np)
      rod(2,1:np) = rod(2,1:np)-linmom(2)*massip(1:np)
      rod(3,1:np) = rod(3,1:np)-linmom(3)*massip(1:np)

      if (ltest) call GetLinMom
      if (ltest) call GetKinEnergy
      if (ltest) write(uout,'(a,3e15.5,f15.5)') 'Zero : linear moments and temp', linmom(1:3), temptra

! ... calculate temptra and scale linear velocities to tempst

      call GetKinEnergy
      fac = sqrt(tempst/temptra)
      rod(1:3,1:np) = rod(1:3,1:np)*fac

! ... check if the linear moments remain sufficiently small

      call GetLinMom
      call GetKinEnergy
      if (ltest) write(uout,'(a,3e15.5,f15.5)') 'Set T: linear moments and temp', linmom(1:3), temptra
      if (count( abs(linmom(1:3)) < tollinmom(1:3) ) == 3) exit

   end do

   if (m > mmax) call Warn('SetZeroLinMom', 'not able to obtain zero linear moments', uout)
   if (ltest) write(uout,'(a,3e15.5,f15.5)') 'final: linear moments and temp', linmom(1:3), temptra

end subroutine SetZeroLinMom

!************************************************************************
!> \page md md.F90
!! **SetZeroAngMom**
!! *set zero angular moments*
!************************************************************************

!      if possible, keep the angular moments consistent with temperature = tempst

subroutine SetZeroAngMom

   use MDModule
   implicit none

   integer(4), save :: mmax = 5
   integer(4) :: m
   real(8)    :: tolangmom(3), fac
   logical    :: ltest =.false.

   if (np < 2) return

! ... calculate tolangmom

   tolangmom(1) = 1.0d-4*sqrt(sum( (massp(1:np)*(ro(2,1:np)*rod(3,1:np)-ro(3,1:np)*rod(2,1:np)))**2   &
                       +(ori(1,1,1:np)*momp(1,1:np)*angvelo(1,1:np)   &
                       + ori(1,2,1:np)*momp(2,1:np)*angvelo(2,1:np)   &
                       + ori(1,3,1:np)*momp(3,1:np)*angvelo(3,1:np))**2 ))
   tolangmom(2) = 1.0d-4*sqrt(sum( (massp(1:np)*(ro(3,1:np)*rod(1,1:np)-ro(1,1:np)*rod(3,1:np)))**2   &
                       +(ori(2,1,1:np)*momp(1,1:np)*angvelo(1,1:np)   &
                       + ori(2,2,1:np)*momp(2,1:np)*angvelo(2,1:np)   &
                       + ori(2,3,1:np)*momp(3,1:np)*angvelo(3,1:np))**2 ))
   tolangmom(3) = 1.0d-4*sqrt(sum( (massp(1:np)*(ro(1,1:np)*rod(2,1:np)-ro(2,1:np)*rod(1,1:np)))**2   &
                       +(ori(3,1,1:np)*momp(1,1:np)*angvelo(1,1:np)   &
                       + ori(3,2,1:np)*momp(2,1:np)*angvelo(2,1:np)   &
                       + ori(3,3,1:np)*momp(3,1:np)*angvelo(3,1:np))**2 ))

   if (ltest) call GetKinEnergy
   if (ltest) write(uout,'(a,3e15.5)')       'tol  : angular moments         ', tolangmom

! ... calculate initial angular moments

   call GetAngMom

   if (ltest) write(uout,'(a,3e15.5,f15.5)') 'init : angular moments and temp', linmom(1:3), temptra

   do m = 1, mmax

! ... adjust angular velocity so angular moments becomes Zero

      angmom(1:3) = angmom(1:3)/np
      angvelo(1,1:np) = angvelo(1,1:np)         &
         -(angmom(1)*ori(1,1,1:np)+angmom(2)*ori(2,1,1:np)+angmom(3)*ori(3,1,1:np))*momip(1,1:np)
      angvelo(2,1:np) = angvelo(2,1:np)         &
         -(angmom(1)*ori(1,2,1:np)+angmom(2)*ori(2,2,1:np)+angmom(3)*ori(3,2,1:np))*momip(2,1:np)
      angvelo(3,1:np) = angvelo(3,1:np)         &
         -(angmom(1)*ori(1,3,1:np)+angmom(2)*ori(2,3,1:np)+angmom(3)*ori(3,3,1:np))*momip(3,1:np)

      if (ltest) call GetAngMom
      if (ltest) call GetKinEnergy
      if (ltest) write(uout,'(a,3e15.5,f15.5)') 'Zero : angular moments and temp', angmom(1:3), temprot

! ... calculate temprot and scale angular velocities to tempst

      call GetKinEnergy
      fac = sqrt(tempst/temprot)
      angvelo(1:3,1:np) = angvelo(1:3,1:np)*fac

! ... check if the angular moments remain sufficiently small

      call GetAngMom
      call GetKinEnergy
      if (ltest) write(uout,'(a,3e15.5,f15.5)') 'Set T: angular moments and temp', angmom(1:3), temprot
      if (count( abs(angmom(1:3)) < tolangmom(1:3) ) == 3) exit

   end do

   if (m > mmax) call Warn('SetZeroAngMom', 'not able to obtain zero angular moments', uout)
   if (ltest) write(uout,'(a,3e15.5,f15.5)') 'final: angular moments and temp', angmom(1:3), temprot

   call AngVelToQuaVel(np, 1, np, qua, angvelo, quad)

end subroutine SetZeroAngMom

!************************************************************************
!> \page md md.F90
!! **GetLinMom**
!! *calculate total linear moment*
!************************************************************************


subroutine GetLinMom
   use MDModule
   implicit none
   linmom(1) = sum(massp(1:np)*rod(1,1:np))
   linmom(2) = sum(massp(1:np)*rod(2,1:np))
   linmom(3) = sum(massp(1:np)*rod(3,1:np))
end subroutine GetLinMom

!************************************************************************
!> \page md md.F90
!! **GetAngMom**
!! *calculate total angular moment*
!************************************************************************


subroutine GetAngMom
   use MDModule
   implicit none
   angmom(1) = sum( massp(1:np)*(ro(2,1:np)*rod(3,1:np)-ro(3,1:np)*rod(2,1:np))  &
                  +ori(1,1,1:np)*momp(1,1:np)*angvelo(1,1:np)                    &
                  +ori(1,2,1:np)*momp(2,1:np)*angvelo(2,1:np)                    &
                  +ori(1,3,1:np)*momp(3,1:np)*angvelo(3,1:np) )
   angmom(2) = sum( massp(1:np)*(ro(3,1:np)*rod(1,1:np)-ro(1,1:np)*rod(3,1:np))  &
                  +ori(2,1,1:np)*momp(1,1:np)*angvelo(1,1:np)                    &
                  +ori(2,2,1:np)*momp(2,1:np)*angvelo(2,1:np)                    &
                  +ori(2,3,1:np)*momp(3,1:np)*angvelo(3,1:np) )
   angmom(3) = sum( massp(1:np)*(ro(1,1:np)*rod(2,1:np)-ro(2,1:np)*rod(1,1:np))  &
                  +ori(3,1,1:np)*momp(1,1:np)*angvelo(1,1:np)                    &
                  +ori(3,2,1:np)*momp(2,1:np)*angvelo(2,1:np)                    &
                  +ori(3,3,1:np)*momp(3,1:np)*angvelo(3,1:np) )
end subroutine GetAngMom

!************************************************************************
!> \page md md.F90
!! **WriteLinMom**
!! *write total linear moment*
!************************************************************************


subroutine WriteLinMom
   use MDModule
   implicit none
   write(uout,'(a,t35,3g12.3)') 'linear moments (x,y,z)         = ', linmom(1:3)
end subroutine WriteLinMom

!************************************************************************
!> \page md md.F90
!! **WriteAngMom**
!! *write total angular moment*
!************************************************************************


subroutine WriteAngMom
   use MDModule
   implicit none
   write(uout,'(a,t35,3g12.3)') 'angular moments (x,y,z)        = ', angmom(1:3)
end subroutine WriteAngMom

!************************************************************************
!> \page md md.F90
!! **GetKinEnergy**
!! *calculate total kinetic energies and temperatures*
!************************************************************************


subroutine GetKinEnergy

   use MDModule
   implicit none

   real(8) :: ekintra, ekinrot, itradegfreetot, irotdegfreetot

   ekintra = sum(Half*massp(1:np)*(rod(1,1:np)**2+rod(2,1:np)**2+rod(3,1:np)**2))
   ekinrot = sum(Half*(momp(1,1:np)*angvelo(1,1:np)**2    &
                      +momp(2,1:np)*angvelo(2,1:np)**2    &
                      +momp(3,1:np)*angvelo(3,1:np)**2))
   ekintra = ekintra*sclmas*sclvel**2/sclene
   ekinrot = ekinrot*sclmas*sclvel**2/sclene
   ekin    = ekintra + ekinrot
   temptra = Zero
   temprot = Zero
   temp    = Zero
   itradegfreetot = sum(itradegfree)
   irotdegfreetot = sum(irotdegfree)
   if (itradegfreetot > 0)           temptra = ekintra*sclene/(Half*(itradegfreetot)*GasConstant)/scltem
   if (irotdegfreetot > 0)           temprot = ekinrot*sclene/(Half*(irotdegfreetot)*GasConstant)/scltem
   if (itradegfreetot+irotdegfreetot > 0) temp = ekin*sclene/(Half*(itradegfreetot+irotdegfreetot)*GasConstant)/scltem

end subroutine GetKinEnergy

!************************************************************************
!> \page md md.F90
!! **GetTStep**
!! *return the time step*
!************************************************************************


function GetTStep()
   use MDModule
   implicit none
   real(8) :: GetTStep
   GetTStep = tstep
end function GetTStep

!************************************************************************
!> \page md md.F90
!! **GetTime**
!! *return the time*
!************************************************************************


function GetTime()
   use MDModule
   implicit none
   real(8) :: GetTime
   if (lsim)  GetTime = ((min(istep1,nstep1)-1)*float(nstep2) + float(min(istep2,nstep2)))*tstep
   if (lana) GetTime = ((min(istep1,nstep1)-1)*float(nstep2) + float(min(istep2,nstep2/idump))*idump)*tstep
end function GetTime

!************************************************************************
!> \page md md.F90
!! **GetlSetVel**
!! *return the value of \ref lsetvel*
!************************************************************************


function GetlSetVel()
   use MDModule
   implicit none
   logical(4) :: GetlSetVel
   GetlSetVel = lsetvel
end function GetlSetVel

!************************************************************************
!> \page md md.F90
!! **GetlZeroMom**
!! *return the value of \ref lzeromom*
!************************************************************************


function GetlZeroMom()
   use MDModule
   implicit none
   logical(4) :: GetlZeroMom
   GetlZeroMom = lzeromom
end function GetlZeroMom

!************************************************************************
!> \page md md.F90
!! **TestSimulationMD**
!! *write MD test output*
!************************************************************************


subroutine TestSimulationMD(unit)

   use MDModule
   implicit none

   character(33), parameter :: fmt33  = '(i4,t10,a,t23,3es12.4,2x,3es12.4)'
   character(33), parameter :: fmt44  = '(i4,t10,a,t23,4es12.4,2x,4es12.4)'
   integer(4),    intent(in) :: unit
   integer(4) :: ip

   write(unit,'(a)') 'forces and torques (lab frame)'
   write(unit,fmt33) (ip,txpt(iptpn(ip)),forceo(1:3,ip),torqueo(1:3,ip),ip = 1,np)
   write(unit,'(a)') 'linear acceleration, (old and new)'
   write(unit,fmt33) (ip,txpt(iptpn(ip)),roddo(1:3,ip),rodd(1:3,ip),ip = 1,np)
   write(unit,'(a)') 'quaternion acceleration (old and new)'
   write(unit,fmt44) (ip,txpt(iptpn(ip)),quaddo(0:3,ip),quadd(0:3,ip),ip = 1,np)
   write(unit,'(a)') 'linear and angular velocities'
   write(unit,fmt33) (ip,txpt(iptpn(ip)),rod(1:3,ip),angvelo(1:3,ip),ip = 1,np)
   write(unit,'(a)') 'quaternions and quaternion velocities'
   write(unit,fmt44) (ip,txpt(iptpn(ip)),qua(0:3,ip),quad(0:3,ip),ip = 1,np)

end subroutine TestSimulationMD

