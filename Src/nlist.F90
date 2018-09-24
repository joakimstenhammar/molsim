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

! relation among interaction-list routines

!     NListDriver
!       !
!       !---------- IONList
!       !
!       !---------- LoadBalance
!       !
!       !---------- NList
!                      !
!                      ! lvlist
!                      !-------- VListAver
!                      !
!                      ! lllist
!                      !-------- LListAver
!                      !
!                      ! lCellList
!                      !-------- CellListAver
!                      !
!                      ! lvlist
!                      !-------- SetVList
!                      !            !
!                      !            !-------- SetVListMD
!                      !            !
!                      !            !-------- SetVListMDLList
!                      !            !
!                      !            !-------- SetVListMC
!                      !            !
!                      !            !-------- SetVListMCLList
!                      !
!                      ! lllist
!                      !-------- SetLList
!                      !
!                      ! lCellList
!                      !-------- SetCellList
!                      !
!                      ! lvlist
!                      !-------- TestVList
!                      !
!                      ! lllist
!                      !-------- TestLList
!                      !
!                      ! lCellList
!                      !-------- TestCellList

!************************************************************************
!> \page nlist nlist.F90
!! **NListModule**
!! *module for neighbour (verlet and liked) list*
!************************************************************************

!> \page nmlIntList
!! The namelist  \ref nmlIntList contains variables that control the calculation of lists of nonbonded pairs to be considered in two-body energy evaluations.
!! * Variables:
!!  * \subpage txintlist
!!  * \subpage lvlistllist
!!  * \subpage inlist
!!  * \subpage drnlist
!!  * \subpage facnneigh

module NListModule

   use MolModule


   logical    :: lnolist      ! flag for no interaction
!> \page inlist
!! `integer`
!! **default:** `0`
!! * `0`: Automatic check whether a new neighbour list should be calculated or not. If the sum of the displacements of any two
!!   molecules, since last calculation, is larger than \ref drnlist a new list is generated. The number of neighbour list calculations is
!!   written after each macrostep.
!! * >`0`: Interval of calculating a new neighbour list.
   integer(4) :: inlist
!> \page drnlist
!! `real`
!! **default:** `0.0`
!! * Distance added to the potential cutoff for calculating the neighbour list.
   real(8)    :: drnlist
!> \page lvlistllist
!! `logical`
!! **default:** `.false.`
!> * `.true.`: Use of liked lists to crease Verlet neighbour lists (only txintlist='nlist').
   logical    :: lvlistllist
   integer(4) :: ndivllist(3) ! number of liked-list cells in one direction
   integer(4) :: ncellllist   ! number of liked-list cells
   real(8)    :: celli(3)     ! 1 / cell lengths
   integer(4) :: npartperproc ! number of particles per processor
   integer(4) :: npmyid       ! number of particles handled by processor myid
   integer(4) :: maxnneigh    ! maximum number of neighbours; used in memory allocation

   real(8), allocatable :: drosum(:,:)  ! sum of displacements

end module NListModule

!************************************************************************
!> \page nlist nlist.F90
!! **NListDriver**
!! *driver of calls for making lists of interacting particles*
!************************************************************************


subroutine NListDriver(iStage)

   use NListModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NListDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IONList(iStage)

   case (iWriteInput)

      call IONList(iStage)
      call LoadBalance !(iStage) iStage is not needed
      call NList(iStage)

   case (iBeforeMacrostep)

      call NList(iStage)

   case (iSimulationStep)

      call NList(iStage)

   case (iAfterMacrostep)

      call NList(iStage)

   case (iAfterSimulation)

      call IONList(iStage)
      call NList(iStage)

   end select

end subroutine NListDriver

!************************************************************************
!> \page nlist nlist.F90
!! **IONList**
!! *perform i/o on neighbour list variables*
!************************************************************************

!> \page txintlist
!! `character(8)`
!! **default:** '`vlist`'
!! * '`vlist`': Verlet neighbour lists. Further specification is given by \ref inlist and \ref drnlist.
!! * '`llist`': Linked lists.

!> \page facnneigh
!! `real`
!! **default:** `2.0`
!! * `0`: Multiplicative constant of radius for calculation number of interacting particles.

subroutine IONList(iStage)

   use NListModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IONList'
   character(80), parameter :: txheading ='neighbour list'

   character(8), save :: txintlist
   real(8)     , save :: facnneigh
   integer(4) :: ierr

   namelist /nmlIntList/ txintlist, lvlistllist, inlist, drnlist, facnneigh

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      txintlist     ='vlist'
      lvlistllist   =.false.
      inlist        = 0
      drnlist       = Zero
      facnneigh     = 2.0D0

      rewind(uin)
      read(uin,nmlIntList)

      call LowerCase(txintlist)

      if(txintlist == 'nlist') txintlist = 'vlist'       ! handle old input'

      if (ldipolesph) txintlist = 'nolist'                 ! force no interaction lists

! ... set lnolist, lvlist, and lllist

      lnolist = .false.
      lvlist = .false.
      lllist = .false.
      lCellList = .false.
      if (txintlist == 'nolist') then
         lnolist = .true.
      else if (txintlist == 'vlist') then
         lvlist = .true.
      else if (txintlist == 'llist') then
         lllist = .true.
      else if (txintlist == 'celllist') then
         lCellList = .true.
      else
         call Stop(txroutine, 'unsupported txintlist', uout)
      end if

! ... determine maximum number of neighbours

      maxnneigh = min(int(facnneigh*(4*pi/3*(rcut+drnlist)**3/vol*np)), np-1)
      maxnneigh = Two*maxnneigh
      maxnneigh = max(250, maxnneigh)                                        ! 2014-10-29
      maxnneigh = max(500, maxnneigh)                                        ! 2014-11-30
      maxnneigh = max(1000, maxnneigh)                                       ! 2015-06-09
      maxnneigh = max(2000, maxnneigh)                                       ! 2015-08-06
      maxnneigh = max(np-1, maxnneigh)                                         ! 2016-08-15
      if(maxnneigh < 0) call Stop('IONlist', 'maxnneigh < 0', uout)

! ... set npartperproc

      if (lmd .or. lmcall .or. lbd) then
         if (lvlistllist) then
            npartperproc = np
         else
            npartperproc = np/nproc + 1
         end if
      else if (lmc) then
         npartperproc = np
      end if
      if(lmvt) npartperproc = facmvt*npartperproc      ! prepare for memory

! ... check unsupported conditions

#if defined (_PAR_)
      if (lllist) call Stop(txroutine, 'lllist .and. _PAR_ not supported ', uout)
      if (lvlistllist) call Stop(txroutine, 'lvlistllist .and. _PAR_ not supported', uout)
#endif
      if (lllist .and. (.not.lbcbox)) call Stop(txroutine, 'lllist .and. non-cubic BCs', uout)

! ... allocate memory

      if (.not.allocated(ipnploc)) then
         allocate(ipnploc(npartperproc), nneighpn(npartperproc), jpnlist(maxnneigh,npartperproc), stat = ierr)
         ipnploc = 0
         nneighpn = 0
         jpnlist = 0
         if(ierr /= 0) then
            write(*,'(a,i10)') 'maxnneigh   = ', maxnneigh
            call WriteIOStat(txroutine, 'memory allocation failed', ierr, 2, 6)
         end if
      end if

   case (iWriteInput)

      if (master) then
         call WriteHead(2, txheading, uout)
         call WriteSub
      end if

   case (iAfterSimulation)

      if (master) then
         call WriteHead(2, txheading, uout)
         call WriteSub
      end if

     ! deallocate(ipnploc, nneighpn, jpnlist)

   end select

contains

!........................................................................

subroutine WriteSub
      if (txintlist == 'nolist') then
         write(uout,'(a)') 'no neighbour list'
      else if (txintlist == 'vlist') then
         write(uout,'(a)') 'verlet list'
         if (lvlistllist) write(uout,'(a)') 'linked lists used to make verlet list'
      else if (txintlist == 'llist') then
         write(uout,'(a)') 'linked list'
      else if (txintlist == 'celllist') then
         write(uout,'(a)') 'cell list'
      end if

      if(lvlist .or. lllist) then
         if (inlist == 0) then
            write(uout,'(a)') 'automatic control of update frequency'
         else if (inlist > 0) then
            write(uout,'(a,t35,i10)') 'update interval                = ', inlist
         end if
      end if
      if(lvlist) then
         write(uout,'(a,t35,f12.3)') 'layer thickness                = ', drnlist
         write(uout,'(a,t35,f12.3)') 'factor used to calculate nneigh= ', facnneigh
         write(uout,'(a,t35,i12  )') 'maximal number of neighbours   = ', maxnneigh
      end if
end subroutine WriteSub

!........................................................................

end subroutine IONList

!いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
!いいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいいい
!************************************************************************
!> \page nlist nlist.F90
!! **LoadBalance**
!! *control the calls of LoadBalanceRealSpace and LoadBalanceRecSpace*
!************************************************************************


subroutine LoadBalance !(iStage) iStage is not needed

   use NListModule
   implicit none

   !integer(4), intent(in) :: iStage

   if (lchain) call LoadBalanceRealSpace(myid, master, nproc, 1, nc, icmyid, itest, 'icmyid', uout)  ! set icmyid
               call LoadBalanceRealSpace(myid, master, nproc, 1, np, ipmyid, itest, 'ipmyid', uout)  ! set ipmyid
               call LoadBalanceRealSpace(myid, master, nproc, 1, na, iamyid, itest, 'iamyid', uout)  ! set iamyid
   if (lewald) call LoadBalanceRecSpace

end subroutine LoadBalance

!************************************************************************
!> \page nlist nlist.F90
!! **LoadBalanceRealSpace**
!! *set iobjmyid for load balancing of single loop with equal load*
!************************************************************************


!     gives nobj/nproc or nobj/nproc+1 objects per process, where nobj = iobjupp - iobjlow + 1

!     use:     do i = iobjmyid(1) , iobjmyid(2)

!     example: nproc = 3, iobjlow = 1, iobjupp = 2       myid   iobjmyid(1)  iobjmyid(2)
!                                                         0        1            0
!                                                         1        1            1
!                                                         2        2            2

!     example: nproc = 3, iobjlow = 1, iobjupp = 5       myid   iobjmyid(1)  iobjmyid(2)
!                                                         0        1            1
!                                                         1        2            3
!                                                         2        4            5

subroutine LoadBalanceRealSpace(myid, master, nproc, iobjlow, iobjupp, iobjmyid, itest, str, unit)

   implicit none

   integer(4), intent(in)  :: myid         ! identity of process
   logical(4), intent(in)  :: master       ! true if master
   integer(4), intent(in)  :: nproc        ! number of processes
   integer(4), intent(in)  :: iobjlow      ! number of the lower object
   integer(4), intent(in)  :: iobjupp      ! number of the upper object
   integer(4), intent(out) :: iobjmyid(2)  ! distribution of objects among the processes
   integer(4), intent(in)  :: itest        ! test flag
   character(*), intent(in):: str          ! variable name
   integer(4), intent(in)  :: unit         ! output unit

   character(40), parameter :: txroutine='LoadBalanceRealSpace'
   integer(4) :: nobj, nobjmyid
#if defined (_PAR_)
   integer(4) :: ivaux(1000), iaux
#endif

   nobj = iobjupp - iobjlow + 1

! ... check conditions

   if (myid < 0) call Stop(txroutine, 'myid < 0', unit)
   if (nproc < 1) call Stop(txroutine, 'nproc < 1', unit)
   if (nobj < 1) call Stop(txroutine, 'nobj < 1', unit)

! ... assign

   iobjmyid(1) = (myid*nobj)/nproc + iobjlow
   iobjmyid(2) = ((myid+1)*nobj)/nproc + iobjlow - 1

! ... sum up

   nobjmyid = iobjmyid(2) - iobjmyid(1) + 1
#if defined (_PAR_)
   call par_allreduce_int(nobjmyid, iaux)
#endif
   if (nobjmyid /= nobj) call Stop(txroutine, 'error in values of iobjmyid', unit)

   if (itest == 4) call TestLoadBalanceRealSpace(unit)

contains

!........................................................................

subroutine TestLoadBalanceRealSpace(unit)
   integer(4), intent(in) :: unit
   integer(4) :: iobjmyidx(2,0:nproc-1), iproc
#if defined (_PAR_)
   iobjmyidx(1:2,0:nproc-1) = 0
#endif
   iobjmyidx(1:2,myid) = iobjmyid(1:2)
#if defined (_PAR_)
   call par_allreduce_ints(iobjmyidx, ivaux, 2*nproc)
#endif
   if (master) then
      call WriteHead(3, 'Test'//trim(txroutine), unit)
      write(unit,'(a)') 'process      '//trim(str)//'(1:2)'
      write(unit,'(a)') '-------      -----------'
      write(unit,'(i4,t13,i6,t20,i6)') (iproc, iobjmyidx(1:2,iproc),  iproc = 0, nproc-1)
   end if
end subroutine TestLoadBalanceRealSpace

!........................................................................

end subroutine LoadBalanceRealSpace

!************************************************************************
!> \page nlist nlist.F90
!! **LoadBalanceRecSpace**
!! *set vbounds for load balancing of reciprocal space loop over nz and ny*
!************************************************************************


!     use:
!
!       ikvec2 = 0
!       do nz = 0, ncut
!          do ny = 0, ncut
!             if (ny**2+nz**2 > ncut2) cycle
!             ikvec2 = ikvec2+1
!             if (ikvec2 < kvecmyid(1) .or. ikvec2 > kvecmyid(2)) cycle

subroutine LoadBalanceRecSpace

   use NListModule
   implicit none

   character(40), parameter :: txroutine='LoadBalanceRecSpace'
   real(8), parameter    :: ddelta = 1.0d-10  ! to ensure correct roundoff
   real(8)    :: aver, sum
   integer(4), allocatable :: load(:,:)
   integer(4) :: vbounds(2,0:mnproc-1)  ! k-vector boundary
   integer(4) :: kvecoffs(0:mnproc-1)   ! k-vector offset
   integer(4) :: iproc, ikvec2, nz, ny, nx, Getnkvec

   allocate(load(0:ncut,0:ncut))
   load = 0

! ... determine load for each (nz,ny)

   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         load(nz,ny) = 0
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! no contribution for odd nx+ny+nz for RD and TO b.c.
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            load(nz,ny) = load(nz,ny)+1
         end do
      end do
   end do

! ... set vbounds(2,0:nproc-1)

   aver = real(Getnkvec())/real(nproc)
   iproc = 0
   sum = Zero
   ikvec2 = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         ikvec2 = ikvec2+1
         sum = sum+load(nz,ny)
         if (sum-(iproc+1)*aver > -ddelta) then
            if (iproc < nproc) kvecoffs(iproc+1) = int(sum)
            vbounds(2,iproc) = ikvec2
            iproc = iproc+1
         end if
      end do
   end do
   kvecoffs(0) = 0

   deallocate(load)

! ... set vbounds(1,0:nproc-1)

   vbounds(1,0) = 1
   vbounds(1,1:nproc-1) = vbounds(2,0:nproc-2)+1

! ... set kvecmyid and kvecoffmyid

   kvecmyid(1:2) = vbounds(1:2,myid)
   kvecoffmyid = kvecoffs(myid)

   if (master .and. itest == 4) call TestLoadBalanceRecSpace(uout)

contains

!........................................................................

subroutine TestLoadBalanceRecSpace(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a)') 'process   kvecmyid(1)   kvecmyid(2)      kvecoffmyid'
   write(unit,'(a)') '-------   -----------   -----------      -----------'
   write(unit,'(i4,t15,i4,t30,i5,t45,i5)') &
        (iproc, vbounds(1,iproc), vbounds(2,iproc), kvecoffs(iproc), iproc = 0, nproc-1)
end subroutine TestLoadBalanceRecSpace

!........................................................................

end subroutine LoadBalanceRecSpace

!************************************************************************
!> \page nlist nlist.F90
!! **NList**
!! *control the calls of SetVList, VListAver, SetLList, and LListAver*
!************************************************************************


subroutine NList(iStage)

   use NListModule
   use CellListModule, only: InitCellList, SetCellList, CellListAver, TestCellList
   implicit none
   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NList'
   integer(4) :: ip, ntotpppair, ntotaapair
   real(8)    :: r2, r2max1, r2max2, naux

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iWriteInput)

      if (lvlist) then
         call SetVList
         call VListAver(iStage)
         if (itest == 4) call TestVList(uout)
      end if
      if (lllist) then
         call SetLList(rcut+drnlist)
         write(uout,'(a,t35,3i10)')  'number of cells in x, y, and z = ', ndivllist(1:3)
         write(uout,'(a,t35,i10)')   'total number of cells          = ', ncellllist
         call LListAver(iStage)
         if (itest == 4) call TestLList(uout)
      end if
      if (lCellList) then
         call InitCellList(rcut + drnlist, iStage)
         call SetCellList()
         call CellListAver(iStage)
         if (itest == 4) call TestCellList(uout)
      end if

      if (.not.allocated(drosum)) then
         allocate(drosum(3,np_alloc))
         drosum = 0.0E+00
      end if

   case (iBeforeMacrostep)

      drosum = Zero

   case (iSimulationStep)

      if(lvlist .or. lllist) then   !do not update the celllist regulary, instead update the celllist in mc.F90 after an successful step
         if (inlist == 0) then                             ! adjusted update interval
            r2max1 = Zero
            r2max2 = Zero
            do ip = 1, np
               drosum(1:3,ip) = drosum(1:3,ip)+drostep(1:3,ip)
               r2 = sum(drosum(1:3,ip)**2)
               if (r2 > r2max1) then
                  r2max1 = r2
               else if (r2 > r2max2) then
                  r2max2 = r2
               end if
            end do
            if (sqrt(r2max1)+sqrt(r2max2) > drnlist) then
               drosum = Zero
               call NListSub
            end if
         else if (inlist > 0) then                        ! regular update interval
            if (mod(istep2,inlist) == 0) call NListSub
         else
            call Stop(txroutine, 'inlist negative', uout)
         end if
      end if

   case (iAfterMacrostep)

      if (lvlist) call VListAver(iStage)
      if (lCellList) call CellListAver(iStage)

   case (iAfterSimulation)

      if (lvlist) then
         call VListAver(iStage)
         call SetVList
         call CalcTotNeighbourPair(ntotpppair, ntotaapair)
         naux = real(ntotaapair)*real((nstep1-(nstep1beg-1))*nstep2)
         if ((naux > Zero) .and. master) then
            write(uout,'(i10,a)')  ntotaapair, ' site-site interactions / step'
            write(uout,*)
         end if
      end if
      if (lllist) then
         write(uout,'(a,t35,3i10)')  'number of cells in x, y, and z = ', ndivllist(1:3)
         write(uout,'(a,t35,i10)')   'total number of cells          = ', ncellllist
         call LListAver(iStage)
      end if
      if (lCellList) call CellListAver(iStage)

      if (allocated(drosum)) deallocate(drosum)
      if (allocated(lcellllist)) deallocate(lcellllist, headllist, jpllist)

   end select

contains

!........................................................................

subroutine NListSub
   if (lvlist) then
      call SetVList
      call VListAver(iStage)
   end if
   if (lllist) then
      call SetLList(rcut+drnlist)
      call LListAver(iStage)
   end if
end subroutine NListSub

!........................................................................

end subroutine NList

!************************************************************************
!> \page nlist nlist.F90
!! **SetVList**
!! *set neighbour (verlet och liked) lists*
!************************************************************************


!     variable                         serial               parallel md         parallel mc
!     --------                         -------              ----------          -----------
!     iploc                            = ip                 local ip            = ip
!     npmyid                           = np                 < np                = np
!     ipnploc(iploc)                   = ip                 = ip                = ip
!     nneighpn(iploc)                  number of neigh      number of neigh     number of neigh
!     jpnlist(nneighpn(iploc),iploc)   = jp                 = jp                = jp

subroutine SetVList

   use NListModule
   implicit none

   character(40), parameter :: txroutine ='SetVList'

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

! ... generate npmyid, ipnploc, nneighpn, and jpnlist

   if (lmd .or. lmcall .or. lbd) then
      if (lvlistllist) then            ! use linked lists to generate verlet list
         call SetLList(rcut+drnlist)
         call SetVListMDLList
      else
         call SetVListMD
      end if
   else if (lmc) then
      if (lvlistllist) then           ! use linked lists to generate verlet list
         call SetLList(rcut+drnlist)
         call SetVListMCLList
      else
         call SetVListMC
      end if
   end if

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine SetVList

!************************************************************************
!> \page nlist nlist.F90
!! **SetVListMD**
!! *set verlet lists for md*
!************************************************************************


subroutine SetVListMD

   use NListModule
   implicit none

   character(10), parameter :: txroutine ='SetVListMD'
   integer(4) :: ip, jp, iploc, itemp
   real(8)    :: rs2cut, r2, dx, dy, dz

   rs2cut = (rcut+drnlist)**2

   iploc = 0
   do ip = ipmyid(1), ipmyid(2)
      iploc = iploc + 1
      if (iploc > npartperproc) call Stop(txroutine, 'iploc > npartperproc', uout)
      nneighpn(iploc) = 0
      ipnploc(iploc) = ip
      do jp = 1, np
         itemp = mod(int(ip+jp),2)
         if (ip < jp .and. itemp == 0) cycle   ! (ip, jp) = (odd, odd) or (even, even)
         if (ip > jp .and. itemp /= 0) cycle   ! (ip, jp) = (odd, even) or (even, odd)
         if (ip == jp) cycle
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rs2cut) then
            nneighpn(iploc) = nneighpn(iploc)+1
            if (nneighpn(iploc) > maxnneigh) call Stop(txroutine, 'nneighpn(iploc) > maxnneigh', uout)
            jpnlist(nneighpn(iploc),iploc) = jp
         end if
      end do
   end do
   npmyid = iploc

end subroutine SetVListMD

!************************************************************************
!> \page nlist nlist.F90
!! **SetVListMDLList**
!! *set verlet lists for md using linked lists*
!************************************************************************

!     serial version

subroutine SetVListMDLList

   use NListModule
   implicit none

   character(15), parameter :: txroutine ='SetVListMDLList'
   integer(4) :: ip, jp
   real(8)    :: rs2cut, r2, dx, dy, dz
   integer    :: icell

   rs2cut = (rcut+drnlist)**2

   do ip = 1, np
      nneighpn(ip) = 0
      ipnploc(ip) = ip

      call SetNCell(ro(1:3,ip), rcut, lcellllist)
      do icell = 1, ncellllist
         if (.not.lcellllist(icell)) cycle
         jp = headllist(icell)

         do
            if (jp == 0) exit
            if (ip > jp .and. mod(int(ip+jp),2) == 0) goto 199
            if (ip < jp .and. mod(int(ip+jp),2) /= 0) goto 199
            if (ip == jp) goto 199
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < rs2cut) then
               nneighpn(ip) = nneighpn(ip)+1
               if (nneighpn(ip) > maxnneigh) call Stop(txroutine, 'nneighpn(ip) > maxnneigh', uout)
               jpnlist(nneighpn(ip),ip) = jp
            end if
199         continue
            jp = jpllist(jp)
         end do
      end do
   end do
   npmyid = np

end subroutine SetVListMDLList

!************************************************************************
!> \page nlist nlist.F90
!! **SetVListMC**
!! *set verlet lists for mc*
!************************************************************************


subroutine SetVListMC

   use NListModule
   implicit none

   character(10), parameter :: txroutine ='SetVListMC'
   integer(4) :: ip, jp
   real(8)    :: rs2cut, r2, dx, dy, dz
#if defined(_PAR_)
   integer(4) :: iploc
#endif

   rs2cut = (rcut+drnlist)**2

#if defined(_PAR_)

   iploc = 0
   do ip = 1, np
      iploc = iploc+1
      nneighpn(iploc) = 0
      ipnploc(iploc) = ip
!      do jp = ipmyid(1), ipmyid(2)
       do jp = myid+1, np, nproc
         if (jp == ip) cycle
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rs2cut) then
            nneighpn(iploc) = nneighpn(iploc)+1
            if (nneighpn(iploc) > maxnneigh) call Stop(txroutine, 'nneighpn(iploc) > maxnneigh', uout)
            jpnlist(nneighpn(iploc),iploc) = jp
         end if
      end do
   end do
   npmyid = iploc

#else

   nneighpn(1:np) = 0
   do ip = 1, np
      ipnploc(ip) = ip
      do jp = ip+1, np
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < rs2cut) then
            nneighpn(ip) = nneighpn(ip)+1
            nneighpn(jp) = nneighpn(jp)+1
            if (nneighpn(ip) > maxnneigh) write(*,*) 'ip, nneighpn(ip), maxnneigh',ip, nneighpn(ip), maxnneigh
            if (nneighpn(jp) > maxnneigh) write(*,*) 'jp, nneighpn(jp), maxnneigh',jp, nneighpn(jp), maxnneigh
            if (nneighpn(jp) > maxnneigh) call Stop(txroutine, 'nneighpn(jp) > maxnneigh', uout)
            if (nneighpn(ip) > maxnneigh) call Stop(txroutine, 'nneighpn(ip) > maxnneigh', uout)
            if (nneighpn(jp) > maxnneigh) call Stop(txroutine, 'nneighpn(jp) > maxnneigh', uout)
            jpnlist(nneighpn(ip),ip) = jp
            jpnlist(nneighpn(jp),jp) = ip
         end if
      end do
   end do
   npmyid = np

#endif

end subroutine SetVListMC

!************************************************************************
!> \page nlist nlist.F90
!! **SetVListMCLList**
!! *set verlet lists for mc using linked lists*
!************************************************************************


subroutine SetVListMCLList

   use NListModule
   implicit none

   character(40), parameter :: txroutine ='SetVListMCLList'
   integer(4) :: ip, jp
   real(8)    :: rs2cut, r2, dx, dy, dz
   integer    :: icell

   rs2cut = (rcut+drnlist)**2

   nneighpn(1:np) = 0
   do ip = 1, np
      ipnploc(ip) = ip

      call SetNCell(ro(1:3,ip), rcut, lcellllist)
      do icell = 1, ncellllist
         if (.not.lcellllist(icell)) cycle
         jp = headllist(icell)
         do
            if (jp == 0) exit
            if (jp <= ip) goto 199

            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < rs2cut) then
               nneighpn(ip) = nneighpn(ip)+1
               nneighpn(jp) = nneighpn(jp)+1
               if (nneighpn(ip) > maxnneigh) call Stop(txroutine, 'nneighpn(ip) > maxnneigh', uout)
               if (nneighpn(jp) > maxnneigh) call Stop(txroutine, 'nneighpn(jp) > maxnneigh', uout)
               jpnlist(nneighpn(ip),ip) = jp
               jpnlist(nneighpn(jp),jp) = ip
            end if
199         continue
            jp = jpllist(jp)
         end do
      end do
   end do
   npmyid = np

end subroutine SetVListMCLList

!************************************************************************
!> \page nlist nlist.F90
!! **VListAver**
!! *calculate statistics involving verlet lists and particle partitioning*
!************************************************************************


subroutine VListAver(iStage)

   use NListModule
   implicit none

   integer(4), intent(in) :: iStage

   character(80), parameter :: txheading ='verlet list statistics'
   integer(4), save :: ncall           ! sum of number of calls of SetVList
   real(8), save    :: nppp(4)         ! number of particles per processor
   real(8), save    :: nneigh(4)       ! number of neigbours per particle and processor
   real(8), save    :: nload(4)        ! load (sum(particle*neighbour)) per processor
   real(8), save    :: dnload(4)       ! largest load of a processor - average load per processor
   real(8)    :: huge, dload, fac
   real(8)    :: npppmyid(0:mnproc-1)
   real(8)    :: nnmin(0:mnproc-1), nnmax(0:mnproc-1)
   real(8)    :: nneighmyid(0:mnproc-1), nneigh2myid(0:mnproc-1)

   select case (iStage)
   case (iWriteInput)

      ncall      = 1         ! 1 for the initial call before the simulation
      nppp(1)    = Zero
      nppp(2)    = Zero
      nppp(3)    =+huge(One)
      nppp(4)    =-huge(One)
      nneigh(1)  = Zero
      nneigh(2)  = Zero
      nneigh(3)  =+huge(One)
      nneigh(4)  =-huge(One)
      nload(1)   = Zero
      nload(2)   = Zero
      nload(3)   =+huge(One)
      nload(4)   =-huge(One)
      dnload(3)  =+huge(One)
      dnload(4)  =-huge(One)

   case (iSimulationStep)

      ncall = ncall+1

   case (iAfterMacrostep)

#if defined (_PAR_)
      npppmyid(0:nproc-1)    = 0
      nneighmyid(0:nproc-1)  = 0
      nneigh2myid(0:nproc-1) = 0
      nnmin(0:nproc-1)       = 0
      nnmax(0:nproc-1)       = 0
#endif

      npppmyid(myid) = npmyid                              ! number of particles on myid
      nneighmyid(myid) = sum(nneighpn(1:npmyid))           ! sum of neighbours on myid
      nneigh2myid(myid) = sum(dble(nneighpn(1:npmyid))**2) ! sum of neighbours squared on myid
      nnmin(myid) = minval(nneighpn(1:npmyid))             ! smallest number of neighbours on myid
      nnmax(myid) = maxval(nneighpn(1:npmyid))             ! largest number of neighbours on myid

#if defined (_PAR_)
      call par_allreduce_reals(npppmyid,    vaux, nproc)
      call par_allreduce_reals(nneighmyid,  vaux, nproc)
      call par_allreduce_reals(nneigh2myid, vaux, nproc)
      call par_allreduce_reals(nnmin,       vaux, nproc)
      call par_allreduce_reals(nnmax,       vaux, nproc)
#endif

!     if (master) then
!        write(uout,'(i5)') nproc
!        write(uout,'(a,10f12.1)') 'npppmyid(0:nproc-1)',   npppmyid(0:nproc-1)
!        write(uout,'(a,10f12.1)') 'nneighmyid(0:nproc-1) ', nneighmyid(0:nproc-1)
!        write(uout,'(a,10f12.1)') 'nneigh2myid(0:nproc-1)', nneigh2myid(0:nproc-1)
!        write(uout,'(a,10f12.1)') 'nnmin(0:nproc-1)     ', nnmin(0:nproc-1)
!        write(uout,'(a,10f12.1)') 'nnmax(0:nproc-1)     ', nnmax(0:nproc-1)
!     end if

      nppp(1) = nppp(1) + sum(npppmyid(0:nproc-1))
      nppp(2) = nppp(2) + sum(npppmyid(0:nproc-1)**2)
      nppp(3) = min(nppp(3),minval(npppmyid(0:nproc-1)))
      nppp(4) = max(nppp(4),maxval(npppmyid(0:nproc-1)))

      nneigh(1) = nneigh(1) + sum(nneighmyid(0:nproc-1))
      nneigh(2) = nneigh(2) + sum(nneigh2myid(0:nproc-1))
      nneigh(3) = min(nneigh(3),minval(nnmin(0:nproc-1)))
      nneigh(4) = max(nneigh(4),maxval(nnmax(0:nproc-1)))

      nload(1) = nload(1) + sum(nneighmyid(0:nproc-1))
      nload(2) = nload(2) + sum(nneighmyid(0:nproc-1)**2)
      nload(3) = min(nload(3),minval(nneighmyid(0:nproc-1)))
      nload(4) = max(nload(4),maxval(nneighmyid(0:nproc-1)))

      dload = maxval(nneighmyid(0:nproc-1))-sum(nneighmyid(0:nproc-1))/nproc
      dnload(1) = dnload(1) + dload
      dnload(2) = dnload(2) + dload**2
      dnload(3) = min(dnload(3),dload)
      dnload(4) = max(dnload(4),dload)

   case (iAfterSimulation)

      if (master) then
         if (lmd .or. lmcall .or. lbd) then
            fac = One
         else if (lmc) then
            fac = nproc
         end if

          nppp(1)   = nppp(1)/(nstep1*nproc)
          nppp(2)   = sqrt(nppp(2)/(nstep1*nproc)-nppp(1)**2)
          nneigh(1) = nneigh(1)/(fac*nstep1*np)
          nneigh(2) = sqrt(nneigh(2)/(fac*nstep1*np)-nneigh(1)**2)
          nload(1)  = nload(1)/(nstep1*nproc)
          nload(2)  = sqrt(nload(2)/(nstep1*nproc)-nload(1)**2)
          dnload(1) = dnload(1)/(nstep1)
          dnload(2) = sqrt(dnload(2)/(nstep1)-dnload(1)**2)

         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,i10)')   'no of neighbour list calls     = ', ncall
         write(uout,'(a,t35,f10.1)') 'no of steps/passes per call    = ', real((nstep1-(nstep1beg-1))*nstep2)/ncall
         write(uout,'()')
         write(uout,'(t12,a,t34,a,t70,a,t85,a)') &
         'no of part per proc', 'no of neighbours per part and proc', 'load of proc', 'load(max)-load(aver)'
         write(uout,'(t12,a,t34,a,t70,a,t85,a)') &
         '-------------------', '----------------------------------', '------------', '--------------------'
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'average', nppp(1), nneigh(1), nload(1), dnload(1)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'one sd ', nppp(2), nneigh(2), nload(2), dnload(2)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'minimum', nppp(3), nneigh(3), nload(3), dnload(3)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'maximum', nppp(4), nneigh(4), nload(4), dnload(4)
      end if

   end select

end subroutine VListAver

!************************************************************************
!> \page nlist nlist.F90
!! **TestVList**
!! *write verlet list and node decomposition test output*
!************************************************************************


subroutine TestVList(unit)

   use NListModule
   implicit none

   integer(4), intent(in) :: unit

   character(40), parameter :: txroutine ='TestVList'
   integer(4) :: iploc, m, iproc, npmyidproc(0:mnproc-1), sumnneighpnproc(0:mnproc-1)
   integer(4), allocatable :: nneighbour(:)

! ... allocate memory

   if (.not.allocated(nneighbour)) then
      allocate(nneighbour(maxnneigh))
      nneighbour = 0
   end if

! ... write process, number of particles, and total number of neighbours

#if defined (_PAR_)
   npmyidproc(0:mnproc-1)      = 0
   sumnneighpnproc(0:mnproc-1) = 0
#endif

   npmyidproc(myid)            = npmyid
   sumnneighpnproc(myid)       = sum(nneighpn(1:npmyid))

#if defined (_PAR_)
   call par_allreduce_ints(npmyidproc,      ivaux, nproc)
   call par_allreduce_ints(sumnneighpnproc, ivaux, nproc)
#endif

   if (master) then
      call WriteHead(3, txroutine, uout)
      write(unit,'(a)') 'process   no of particles      total no of neighbours'
      write(unit,'(a)') '-------   ---------------      ----------------------'
      do iproc = 0, nproc-1
         write(unit,'(i4,t15,i4,t40,i10)') iproc, npmyidproc(iproc), sumnneighpnproc(iproc)
      end do
   end if

! ... write neigbouring particles for each particle on master

   if (master) then
      write(unit,'()')
      write(unit,'(a)') 'neighbour lists (master only)'
      write(unit,'(a)') '-----------------------------'
      write(unit,'(a)') 'no of nb      frequence'
      write(unit,'(a)') '--------      ---------'
      nneighbour(1:maxnneigh) = 0
      do iploc = 1, npmyid
         nneighbour(nneighpn(iploc)) = nneighbour(nneighpn(iploc))+1
      end do
      write(unit,'(i5,i15)') (m,nneighbour(m),m = minval(nneighpn(1:npmyid)),maxval(nneighpn(1:npmyid)))
      write(unit,'(a)') 'particle   no of nb     neighbours'
      write(unit,'(a)') '--------   --------     ----------'
      do iploc = 1, npmyid
         write(unit,'(i5,i10)')     ipnploc(iploc), nneighpn(iploc)
         write(unit,'(t25,20(i5))') jpnlist(1:nneighpn(iploc),iploc)
      end do
   end if

   write(unit,'()')
   write(unit,'(a)') 'neighbour lists (all nodes)'
   write(unit,'(a)') '---------------------------'
   write(unit,'(a)') 'node   particle   no of nb     neighbours'
   write(unit,'(a)') '----   --------   --------     ----------'
   do iploc = 1, npmyid
      write(unit,'(i2,t7,i5,i10)') myid, ipnploc(iploc), nneighpn(iploc)
      write(unit,'(t31,20(i5))') jpnlist(1:nneighpn(iploc),iploc)
   end do

   deallocate(nneighbour)

end subroutine TestVList

!************************************************************************
!> \page nlist nlist.F90
!! **SetLList**
!! *set linked lists*
!************************************************************************

!     serial version

subroutine SetLList(distance)

   use NListModule
   implicit none

   real(8), intent(in) :: distance

   character(40), parameter :: txroutine ='SetLList'
   integer(4) ::    ip, ix, iy, iz, icell, ierr

   if (txbc /= 'xyz') call Stop(txroutine, 'txbc /= ''xyz''', uout)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   do ip = 1, np
     ipnploc(ip) = ip
   end do
   npmyid = np

! ... set ndivllist, ncellllist, and celli

   ndivllist(1:3) = int(boxlen(1:3)/distance)
   ncellllist = product(ndivllist(1:3))
   celli(1:3) = real(ndivllist(1:3))*boxleni(1:3)

! ... allocate memory

      ncellllist = min(ncellllist, 1000)
      if (.not.allocated(lcellllist)) then
         allocate(lcellllist(ncellllist), headllist(ncellllist), jpllist(np_alloc), stat = ierr)
         lcellllist = .false.
         headllist = 0
         jpllist = 0
         if(ierr /= 0) then
            write(*,'(a,i10)') 'ncellllist   = ', ncellllist
            call WriteIOStat(txroutine, 'memory allocation failed', ierr, 2, 6)
         end if
      end if

! ... set headllist and jpllist

   headllist(1:ncellllist) = 0
   do ip = np, 1, -1       ! reverse order to get lists identical to previous code (no problem with ip = 1, np)
      ix = aint((ro(1,ip)+boxlen2(1))*celli(1))
      iy = aint((ro(2,ip)+boxlen2(2))*celli(2))
      iz = aint((ro(3,ip)+boxlen2(3))*celli(3))
      icell = 1 + ix + iy*ndivllist(1) + iz*ndivllist(1)*ndivllist(2)
      jpllist(ip) = headllist(icell)
      headllist(icell) = ip
   end do

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine SetLList

!************************************************************************
!> \page nlist nlist.F90
!! **SetNCell**
!! *set list of cells to be used for a position*
!************************************************************************

!     note, sum(lcellllist(1:ncellllist)) < 27 is possible and may apear when distance < cell length

subroutine SetNCell(rr, distance, lcellllistloc)

   use NListModule
   implicit none

   real(8), intent(in)    :: rr(3)                       ! position
   real(8), intent(in)    :: distance                    ! radius
   logical, intent(out)   :: lcellllistloc(*)            ! list of cells within rr+-distance

   real(8)                :: dx, dy, dz
   integer(4)             :: ix, iy, iz, icell, idx, idy, idz

   lcellllistloc(1:ncellllist) = .false.
   do idx = -1, 1
      do idy = -1, 1
         do idz = -1, 1
            dx = rr(1)+idx*distance
            dy = rr(2)+idy*distance
            dz = rr(3)+idz*distance
            call PBC(dx,dy,dz)
            ix = aint((dx+boxlen2(1))*celli(1))
            iy = aint((dy+boxlen2(2))*celli(2))
            iz = aint((dz+boxlen2(3))*celli(3))
            icell = 1 + ix + iy*ndivllist(1) + iz*ndivllist(1)*ndivllist(2)
            lcellllistloc(icell) = .true.
         end do
      end do
   end do

end subroutine SetNCell

!************************************************************************
!> \page nlist nlist.F90
!! **LListAver**
!! *calculate statistics involving linked cell lists and particle partitioning*
!************************************************************************


subroutine LListAver(iStage)

   use NListModule
   implicit none

   integer(4), intent(in) :: iStage

   character(80), parameter :: txheading ='linked list statistics'
   integer(4), save :: ncall           ! sum of number of calls of SetLList

   select case (iStage)
   case (iWriteInput)

      ncall = 1                        ! 1 for the inital call before the simulation

   case (iSimulationStep)

      ncall = ncall+1

   case (iAfterSimulation)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,i10)')   'no of linked cell list calls   = ', ncall
         write(uout,'(a,t35,f10.1)') 'no of steps/passes per call    = ', real((nstep1-(nstep1beg-1))*nstep2)/ncall
      end if

   end select

end subroutine LListAver

!************************************************************************
!> \page nlist nlist.F90
!! **TestLList**
!! *write linked list*
!************************************************************************


subroutine TestLList(unit)

   use NListModule
   implicit none
   integer(4), intent(in) :: unit
   integer(4) :: icell, jp

   call WriteHead(3, 'TestLList', uout)
   write(unit,'(a)') 'particle     list (next particle)'
   write(unit,'(i5,i10)') (jp,jpllist(jp),jp=1,np)
   write(unit,'(a)') ' cell      particles (first is head)'
   do icell = 1, ncellllist
      write(unit,'(i5)') icell
      jp = headllist(icell)
      do
         if (jp == 0) exit
         write(unit,'(i15)') jp
         jp = jpllist(jp)
      end do
   end do

end subroutine TestLList

!************************************************************************
!> \page nlist nlist.F90
!! **CalcTotNeighbourPair**
!! *calculate number of neigbouring pairs*
!************************************************************************


subroutine CalcTotNeighbourPair(ntotpppair, ntotaapair)

   use NListModule
   implicit none

   integer(4), intent(out) :: ntotpppair  ! total number of particle-particle pairs
   integer(4), intent(out) :: ntotaapair  ! total number of atom-atom pairs

   integer(4) :: ip, iploc, ipt, jp, jploc, jpt

   ntotpppair = 0
   ntotaapair = 0
   do iploc = 1, npmyid
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         jpt = iptpn(jp)
         ntotpppair = ntotpppair+1
         ntotaapair = ntotaapair+napt(ipt)*napt(jpt)
      end do
   end do

#if defined (_PAR_)
   call par_allreduce_int(ntotpppair, iaux)      ! allreduce of ntotpppair
   call par_allreduce_int(ntotaapair, iaux)      ! allreduce of ntotaapair
#endif

end subroutine CalcTotNeighbourPair

!************************************************************************
!> \page nlist nlist.F90
!! **Getdrnlist**
!! *return the value of \ref drnlist*
!************************************************************************


function Getdrnlist()
   use NListModule
   implicit none
   real(8) :: Getdrnlist
   Getdrnlist = drnlist
end function Getdrnlist

!************************************************************************
!> \page nlist nlist.F90
!! **Getnpmyid**
!! *return the value of npmyid*
!************************************************************************


function Getnpmyid()
   use NListModule
   implicit none
   integer(4) :: Getnpmyid
   Getnpmyid = npmyid
end function Getnpmyid

!************************************************************************
!> \page nlist nlist.F90
!! **Getdrnlist**
!! *return the value of ncellllist*
!************************************************************************


function Getncellllist()
   use NListModule
   implicit none
   integer(4) :: Getncellllist
   Getncellllist = ncellllist
end function Getncellllist
