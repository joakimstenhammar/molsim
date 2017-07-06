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
!*     Group                                                            *
!*                                                                      *
!************************************************************************

! ... classify reference and field particles into groups and make averages

subroutine Group(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Group'
   character(20) :: ref           ! text label for selecting procedure of how to divide particles into rerernce groups
   character(20) :: field         ! text lable for selecting procedure of how to divide particles into field groups
   logical       :: lwref         !.true. => reference group data are written on file FGROUP
   character(20) :: txtype(2)     ! holds ref and field
   character(10), allocatable :: no(:)    ! for internal reading
   logical       :: lsetconf
   integer(4)    :: ip, ipt, igr, jgr,  m, ivar
   character(80) :: str
   logical       :: ok
   save

   namelist /nmlGroup/ ref, field, lwref

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

      ref   = 'type=all'
      field = 'type=all'
      lwref =.false.

      rewind(uin)
      call Advance('nmlGroup',uin,str,ok)
      rewind(uin)
      if (ok) read(uin,nmlGroup)                ! an omission of nmlGroup is technically allowed

      call LowerCase(ref)
      call LowerCase(field)

      txtype(1) = ref
      txtype(2) = field

! ... set ngr

      do m = 1, 2
         ngr(m) = 0
         if (txtype(m)(1:8) == 'type=all') then
            call GroupAll(iStage, m)
         else if (txtype(m)(1:5) == 'type=') then
            call GroupType(iStage, m, txtype)
         else
            call GroupUser(iStage, m, txtype, lsetconf)
            if (.not.lsetconf) call Stop(txroutine, 'unsupported value of ref or field', uout)
         end if
         if (ngr(m) <= 0)   call Stop(txroutine, 'ngr(m) <= 0', uout)
      end do

! ... set maxngr, ngrgr, and ngrvar

      maxngr = maxval(ngr(1:2))
      ngrgr = ngr(1)*ngr(2)
      ngrvar = 2+sum(ngr(1:2))

   case (iWriteInput)

! ... allocate memory

      allocate(grvar(ngrvar))
      allocate(igrpn(np_alloc,2))
      igrpn = 0
      allocate(iptgr(maxngr,2))
      iptgr = 0
      allocate(iatgr(maxngr,2))
      iatgr = 0
      allocate(natgr(maxngr,2))
      natgr = 0
      allocate(igrgr(maxngr,maxngr))
      igrgr = 0
      allocate(igrpnt(2,0:maxngr))
      igrpnt = 0
      allocate(txgr(maxngr))
      txgr = ""
      allocate(txgrgr(ngrgr))
      txgrgr = ""
      allocate(no(maxngr))
      no = ""

! ... set igrpnt

      ivar = 0
      do m = 1, 2
         do igr = 0, ngr(m)
            ivar = ivar + 1
            igrpnt(m,igr) = ivar
         end do
      end do

! ... set iptgr and label

      do m = 1, 2
         grvar(igrpnt(m,0))%label = 'not in any group'
         if (txtype(m)(1:8) == 'type=all') then
            call GroupAll(iStage, m)
         else if (txtype(m)(1:5) == 'type=') then
            call GroupType(iStage, m, txtype)
         else
            call GroupUser(iStage, m, txtype, lsetconf)
            if (.not.lsetconf) call Stop(txroutine, 'unsupported value of ref or field', uout)
         end if
      end do

! ... set ngrgr, igrgr, txgr, and txgrgr

      do igr = 1, maxngr
         write(no(igr),'(i10)') igr
      end do
      ngrgr = 0
      do igr = 1, ngr(1)
         txgr(igr) = 'gr:'//trim(adjustl(no(igr)))
         do jgr = 1, ngr(2)
            ngrgr = ngrgr+1
            igrgr(igr,jgr) = ngrgr
            txgrgr(ngrgr) = 'gr:'//trim(adjustl(no(igr)))//','//trim(adjustl(no(jgr)))
         end do
      end do

! ... set natgr and iatgr with help of iptgr

      do m = 1, 2
         do igr = 1, ngr(m)
            ipt = iptgr(igr,m)
            natgr(igr,m) = natpt(ipt)
            iatgr(igr,m) = iatpt(ipt)
         end do
      end do

      if (lwref .and. master) call FileOpen(ugroup, fgroup, 'form/noread')   ! open FGROUP

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, ngrvar, grvar)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) grvar

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, ngrvar, grvar)

   case (iSimulationStep)

! ... calculate igrpn and grvar()%value

      do m = 1, 2
         grvar(igrpnt(m,0:ngr(m)))%value = 0
         if (txtype(m)(1:8) == 'type=all') then
            call GroupAll(iStage, m)
         else if (txtype(m)(1:5) == 'type=') then
            call GroupType(iStage, m, txtype)
         else
            call GroupUser(iStage, m, txtype, lsetconf)
            if (.not.lsetconf) call Stop(txroutine, 'unsupported value of ref or field', uout)
         end if

! ... check that particles divided into groups are of right type


         do ip = 1, np
            igr = igrpn(ip,m)
            if (igr >= 1) then
               if (iptpn(ip) /= iptgr(igr,m)) call Stop(txroutine, 'iptpn(ip) /= iptgr(igr,m)', uout)
            end if
         end do

      end do

! ... sample

      call ScalarSample(iStage, 1, ngrvar, grvar)

      if (lwref .and. master) write(ugroup,'(200i1)') igrpn(1:np,1)   ! write group data

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, ngrvar, grvar)
      if (lsim .and. master) write(ucnf) grvar

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, ngrvar, grvar)
         call WriteHead(2, 'group data', uout)
         if (lwref) write(uout,'(a)') 'ref. group data written'
         do m = 1, 2
            write(uout,'()')
            if (m == 1) write(uout,'(a,a)') 'ref. particles: ', txtype(m)
            if (m == 2) write(uout,'(a,a)') 'field particles: ', txtype(m)
            write(uout,'(a,t7,a,t30,a,5x,a)') 'no', 'group label', 'no of particles', 'precision'
            write(uout,'(a,t7,a,t30,a,5x,a)') '--', '-----------', '---------------', '---------'
            write(uout,'()')
            do igr = 0, ngr(m)
               ivar = igrpnt(m,igr)
               write(uout,'(i2,t7,a,t30,f10.3,7x,f10.3)') igr, trim(grvar(ivar)%label), grvar(ivar)%avs1, grvar(ivar)%avsd
            end do
         end do
      end if

      deallocate(grvar)
!     deallocate(igrpn)        ! needen in dynamic
      deallocate(iptgr)
!     deallocate(iatgr)        ! needed in static (for pressure evaluation)
!     deallocate(natgr)        ! needed in static (for pressure evaluation)
!     deallocate(igrgr)        ! needed in static (for pressure evaluation)
      deallocate(igrpnt)
!     deallocate(txgr)         ! needed in statistics at case (iAfterSimulation)
      deallocate(txgrgr)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine Group

!************************************************************************
!*                                                                      *
!*     GroupAll                                                         *
!*                                                                      *
!************************************************************************

! ... group division: each particle type forms one group

subroutine GroupAll(iStage, m)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m

   integer(4) :: ip, igr

   select case (iStage)
   case (iReadInput)
      ngr(m) = npt
   case (iWriteInput)
      do igr = 1, ngr(m)
         iptgr(igr,m) = igr
         grvar(igrpnt(m,igr))%label = txpt(igr)
      end do
   case (iSimulationStep)
      do ip = 1, np
         igr = iptpn(ip)
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupAll

!************************************************************************
!*                                                                      *
!*     GroupType                                                        *
!*                                                                      *
!************************************************************************

! ... group division: select particle type forms one group

subroutine GroupType(iStage, m, txtype)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   integer(4),    intent(in) :: m
   character(20), intent(in) :: txtype(2)

   integer(4) :: ip, ipt, igr

   select case (iStage)
   case (iReadInput)
      do ipt = 1, npt
         if (txtype(m) == 'type='//Trim(txpt(ipt))) ngr(m) = 1
      end do
   case (iWriteInput)
      do ipt = 1, npt
         if (txtype(m) == 'type='//Trim(txpt(ipt))) then
            iptgr(1,m) = ipt
            grvar(igrpnt(m,ipt))%label = txpt(ipt)
         end if
      end do
   case (iSimulationStep)
      do ip = 1, np
         igr = 0
         if (iptpn(ip) == iptgr(1,m)) igr = 1
         igrpn(ip,m) = igr
         grvar(igrpnt(m,igr))%value = grvar(igrpnt(m,igr))%value + 1
      end do
   end select

end subroutine GroupType
