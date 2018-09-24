! ... 'version 6.4.7, Sep 18, 2015'

!************************************************************************
!************************************************************************
!**                                                                    **
!**  Copyright 1997                                                    **
!**                                                                    **
!**  Per Linse                                                         **
!**  Physical Chemistry 1                                              **
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

#if defined (_PAR_)

!************************************************************************
!************************************************************************
!**                                                                    **
!**  this file contains the parallel interface to which molsim         **
!**  is connected. the precise implementation of the parallelism       **
!**  is thus relied upon here. in order to keep everything as          **
!**  simple as possible, all error handling is supposed to take        **
!**  place in the routines in this file.                               **
!**                                                                    **
!************************************************************************
!************************************************************************

!************************************************************************
!> \page parallel parallel.F90
!! **ParallelModule**
!! *module to be used in routines invoking mpi*
!************************************************************************


module ParallelModule

   include 'mpif.h'

   integer(4), parameter :: rootid = 0             ! id of root process

end module ParallelModule

!************************************************************************
!> \page parallel parallel.F90
!! **par_initialize**
!! *initialize parallel matter*
!************************************************************************


subroutine par_initialize

# ifdef F03_CBIND
   use, intrinsic :: iso_c_binding       ! for FFTW
# endif
   use ParallelModule

   implicit none
   integer(4) :: ierr

! ... we need to do some checking for compatibility of different real and complex types

  real(8) :: std_real
  complex(8) :: std_complex
  real(C_DOUBLE) :: fftw_real
  complex(C_DOUBLE_COMPLEX) :: fftw_complex

! ... initialize mpi

   call mpi_init(ierr)
   if (ierr/=mpi_success) call par_error('par_init',ierr)

! ... testing MOLSIM and FFTW compatibility

    if (kind(std_real) .ne. kind(fftw_real) ) then
      call stop('par_initialize', 'REAL types from FFTW incompatible with Molsim REAL types',6)
   endif
   if (kind(std_complex) .ne. kind(fftw_complex) ) then
      call stop('par_initialize', 'COMPLEX types from FFTW incompatible with Molsim COMPLEX types',6)
   endif

end subroutine par_initialize

!************************************************************************
!> \page parallel parallel.F90
!! **par_finalize**
!! *finalize parallel matter*
!************************************************************************


subroutine par_finalize

   use ParallelModule
   implicit none

   integer(4) :: ierr

! ... finalize mpi

   call mpi_finalize(ierr)
   if (ierr/=mpi_success) call par_error('par_init',ierr)

end subroutine par_finalize

!************************************************************************
!> \page parallel parallel.F90
!! **par_comm_size**
!! *get number of processes*
!************************************************************************


subroutine par_comm_size(nproc)

   use ParallelModule
   implicit none

   integer(4), intent(out) :: nproc
   integer(4) :: ierr

! ... get number of processes

   call mpi_comm_size(mpi_comm_world,nproc,ierr)
   if (ierr/=mpi_success) call par_error('par_comm_size',ierr)

end subroutine par_comm_size

!************************************************************************
!> \page parallel parallel.F90
!! **par_comm_rank**
!! *get my id and set master and slave*
!************************************************************************


subroutine par_comm_rank(myid,master,slave)

   use ParallelModule
   implicit none

   integer(4), intent(in)  :: myid
   logical,    intent(out) :: master
   logical,    intent(out) :: slave

   integer(4) :: ierr

! ... get my id

   call mpi_comm_rank(mpi_comm_world,myid,ierr)

! ... set master and slave

   master=.false.
   if (myid==rootid) master = .true.
   slave = .not.master

end subroutine par_comm_rank

!************************************************************************
!> \page parallel parallel.F90
!! **par_barrier**
!! *barrier synchronisation*
!************************************************************************


subroutine par_barrier

   use ParallelModule
   implicit none

   integer(4)               :: ierr

   call mpi_barrier(mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_barrier',ierr)

end subroutine par_barrier

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_characters**
!! *broadcast character variables*
!************************************************************************


subroutine par_bc_characters(buff,icount)

   use ParallelModule
   character(*), intent(in) :: buff
   integer(4)  , intent(in) :: icount

   integer(4)               :: ierr, ilen

   ilen = len(buff)
   call mpi_bcast(buff,ilen*icount,mpi_character,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_characters',ierr)

end subroutine par_bc_characters

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_logicals**
!! *broadcast logicals*
!************************************************************************


subroutine par_bc_logicals(buff,icount)

   use ParallelModule
   implicit none

   logical   , intent(in)   :: buff(*)
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr

   call mpi_bcast(buff,icount,mpi_logical,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_logicals',ierr)

end subroutine par_bc_logicals

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_logical**
!! *broadcast logical*
!************************************************************************


subroutine par_bc_logical(buff)

   use ParallelModule
   implicit none

   logical   , intent(inout)   :: buff

   integer(4)               :: ierr

   call mpi_bcast(buff,1,mpi_logical,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_logicals',ierr)

end subroutine par_bc_logical

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_ints**
!! *broadcast integers*
!************************************************************************


subroutine par_bc_ints(buff,icount)

   use ParallelModule
   implicit none

   integer(4), intent(in)   :: buff(*)
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr

   call mpi_bcast(buff,icount,mpi_integer4,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_ints',ierr)

end subroutine par_bc_ints

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_int**
!! *broadcast of scalar integer*
!************************************************************************


subroutine par_bc_int(buff)

   use ParallelModule
   implicit none

   integer(4), intent(inout)   :: buff

   integer(4)               :: ierr

   call mpi_bcast(buff,1,mpi_integer4,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_int',ierr)

end subroutine par_bc_int

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_ints8**
!! *broadcast integers of kind 8*
!************************************************************************


subroutine par_bc_ints8(buff,icount)

   use ParallelModule
   implicit none

   integer(8), intent(inout)   :: buff(*)
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr

   call mpi_bcast(buff,icount,mpi_integer8,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_ints8',ierr)

end subroutine par_bc_ints8

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_reals**
!! *broadcast double precision reals*
!************************************************************************


subroutine par_bc_reals(buff,icount)

   use ParallelModule
   implicit none

   real(8)   , intent(in)   :: buff(*)
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr

   call mpi_bcast(buff,icount,mpi_real8,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_reals',ierr)

end subroutine par_bc_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_real**
!! *broadcast double precision real*
!************************************************************************


subroutine par_bc_real(buff)

   use ParallelModule
   implicit none

   real(8)   , intent(inout)   :: buff

   integer(4)               :: ierr

   call mpi_bcast(buff,1,mpi_real8,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_reals',ierr)

end subroutine par_bc_real

!************************************************************************
!> \page parallel parallel.F90
!! **par_bc_comps**
!! *broadcast double precision complex variables*
!************************************************************************


subroutine par_bc_comps(buff,icount)

   use ParallelModule
   implicit none

   complex(8), intent(in)   :: buff(*)
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr

   call mpi_bcast(buff,icount,mpi_double_complex,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_bc_comps',ierr)

end subroutine par_bc_comps


!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_logicals**
!! *perform global logical or and redistribution of logical variables*
!************************************************************************


subroutine par_allreduce_logicals(buff,temp,icount)

   use ParallelModule
   implicit none
   logical, intent(out)     :: buff(*)
   logical, intent(in)      :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_logical,mpi_lor,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_logicals',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_allreduce_logicals

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_logical**
!! *perform global logical or and redistribution of logical variables*
!************************************************************************


subroutine par_allreduce_logical(buff,temp)

   use ParallelModule
   implicit none
   logical, intent(inout)   :: buff
   logical, intent(in)      :: temp  ! temporary variable, should be able to hold buff
   integer(4), parameter    :: icount = 1

   integer(4)               :: ierr

   call mpi_allreduce(buff,temp,icount,mpi_logical,mpi_lor,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_logicals',ierr)
   buff = temp

end subroutine par_allreduce_logical

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_ints**
!! *perform global sum and redistribution of integer variables*
!************************************************************************



subroutine par_allreduce_ints(buff,temp,icount)

   use ParallelModule
   implicit none
   integer(4), intent(out)  :: buff(*)
   integer(4), intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_integer4,mpi_sum,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_ints',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_allreduce_ints

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_int**
!! *perform global sum and redistribution of integer variable*
!************************************************************************



subroutine par_allreduce_int(buff,temp)

   use ParallelModule
   implicit none
   integer(4), intent(inout)  :: buff
   integer(4), intent(in)     :: temp  ! temporary variable, should be able to hold buff
   integer(4), parameter      :: icount = 1

   integer(4)               :: ierr

   call mpi_allreduce(buff,temp,icount,mpi_integer4,mpi_sum,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_ints',ierr)
   buff = temp

end subroutine par_allreduce_int

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_reals**
!! *perform global sum and redistribution of double precision variables*
!************************************************************************


subroutine par_allreduce_reals(buff,temp,icount)

   use ParallelModule
   implicit none
   real(8)   , intent(out)  :: buff(*)
   real(8)   , intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_real8,mpi_sum,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_reals',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_allreduce_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_real**
!! *perform global sum and redistribution of double precision variable*
!************************************************************************


subroutine par_allreduce_real(buff,temp)

   use ParallelModule
   implicit none
   real(8)   , intent(inout)  :: buff
   real(8)   , intent(in)     :: temp  ! temporary variable, should be able to hold buff
   integer(4), parameter      :: icount = 1

   integer(4)               :: ierr

   call mpi_allreduce(buff,temp,icount,mpi_real8,mpi_sum,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_reals',ierr)
   buff = temp

end subroutine par_allreduce_real

!************************************************************************
!> \page parallel parallel.F90
!! **par_allreduce_comps**
!! *perform global sum and redistribution of double precision complex variables*
!************************************************************************


subroutine par_allreduce_comps(buff,temp,icount)

   use ParallelModule
   implicit none
   complex(8), intent(out)  :: buff(*)
   complex(8), intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allreduce_comps',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_allreduce_comps

!************************************************************************
!> \page parallel parallel.F90
!! **par_max_ints**
!! *find maximum value and redistribute integer variables*
!************************************************************************



subroutine par_max_ints(buff,temp,icount)

   use ParallelModule
   implicit none
   integer(4), intent(out)  :: buff(*)
   integer(4), intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_integer4,mpi_max,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_max_ints',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_max_ints

!************************************************************************
!> \page parallel parallel.F90
!! **par_min_reals**
!! *find minimum value and redistribute double precision variables*
!************************************************************************


subroutine par_min_reals(buff,temp,icount)

   use ParallelModule
   implicit none

   real(8)   , intent(out)  :: buff(*)
   real(8)   , intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_real8,mpi_min,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_min_reals',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_min_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_max_reals**
!! *find maximum value and redistribute double precision variables*
!************************************************************************


subroutine par_max_reals(buff,temp,icount)

   use ParallelModule
   implicit none

   real(8)   , intent(out)  :: buff(*)
   real(8)   , intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_allreduce(buff,temp,icount,mpi_real8,mpi_max,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_max_reals',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_max_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_reduce_reals**
!! *perform global sum and redistribution of double precision variables*
!************************************************************************


subroutine par_reduce_reals(buff,temp,icount)

   use ParallelModule
   implicit none

   real(8)   , intent(out)  :: buff(*)
   real(8)   , intent(in)   :: temp(*)  ! temporary array, should be able to hold buff
   integer(4), intent(in)   :: icount

   integer(4)               :: ierr, i

   call mpi_reduce(buff,temp,icount,mpi_real8,mpi_sum,rootid,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_reduce_reals',ierr)
   do i = 1,icount
      buff(i) = temp(i)
   end do

end subroutine par_reduce_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_allgather_shortints**
!! *perform a distribution of array of integer(2)*
!************************************************************************


subroutine par_allgather_shortints(buff,temp,isendcount,isendoffs,ireceivecounts,iptrs)

   use ParallelModule
   implicit none

   integer(2)               :: buff(*)
   integer(2)               :: temp(*)
   integer(4), intent(in)   :: isendcount,isendoffs,ireceivecounts(*)
   integer(4), intent(in)   :: iptrs(*)

   integer(4)               :: ierr,i

   do i = 1,isendcount
      temp(i)=buff(i+isendoffs)
   end do
   call mpi_allgatherv(temp,isendcount,mpi_integer2,buff,ireceivecounts,iptrs,mpi_integer2,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allgather_shortints',ierr)

end subroutine par_allgather_shortints

!************************************************************************
!> \page parallel parallel.F90
!! **par_allgather_ints**
!! *perform a distribution of array of integer(4)*
!************************************************************************


subroutine par_allgather_ints(buff,temp,isendcount,isendoffs,ireceivecounts,iptrs)

   use ParallelModule
   implicit none

   integer(4)               :: buff(*)
   integer(4)               :: temp(*)
   integer(4), intent(in)   :: isendcount,isendoffs,ireceivecounts(*)
   integer(4), intent(in)   :: iptrs(*)
   integer(4)               :: ierr,i

   do i = 1,isendcount
      temp(i)=buff(i+isendoffs)
   end do
   call mpi_allgatherv(temp,isendcount,mpi_integer4,buff,ireceivecounts,iptrs,mpi_integer4,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_allgather_ints',ierr)

end subroutine par_allgather_ints

!************************************************************************
!> \page parallel parallel.F90
!! **par_scatter_reals**
!! *scatter double precision data from root to the other processes*
!************************************************************************


subroutine par_scatter_reals(sbuff,sicount,rbuff,ricount,temp)

   use ParallelModule
   implicit none
   real(8)   , intent(in)   :: sbuff(*)
   integer(4), intent(in)   :: sicount
   real(8)   , intent(out)  :: rbuff(*)
   integer(4), intent(out)  :: ricount
   real(8)   , intent(in)   :: temp(*)  ! temporary array, should be able to hold rbuff

   integer(4)               :: ierr,i

   call mpi_scatter(sbuff,sicount,mpi_real8,rbuff,ricount,mpi_real8,0,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_scatter_reals',ierr)
   do i = 1,ricount
      rbuff(i) = temp(i)
   end do

end subroutine par_scatter_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_gather_reals**
!! *gather double precision data from all processes except root to root*
!************************************************************************


subroutine par_gather_reals(sbuff,sicount,rbuff,ricount,temp)

   use ParallelModule
   implicit none
   real(8)   , intent(in)   :: sbuff(*)
   integer(4), intent(in)   :: sicount
   real(8)   , intent(out)  :: rbuff(*)
   integer(4), intent(out)  :: ricount
   real(8)   , intent(inout):: temp(*)  ! temporary array, should be able to hold rbuff

   integer(4)               :: ierr,i

   do i = 1,sicount
      temp(i) = sbuff(i)
   end do
   call mpi_gather(temp,sicount,mpi_real8,rbuff,ricount,mpi_real8,0,mpi_comm_world,ierr)
   if (ierr/=mpi_success) call par_error('par_gather_reals',ierr)

end subroutine par_gather_reals

!************************************************************************
!> \page parallel parallel.F90
!! **par_handshake**
!! *hand shaking between master and slaves*
!************************************************************************


subroutine par_handshake(myid,master,slave,nproc,unit)

   use ParallelModule
   implicit none

   integer(4), intent(in) :: myid
   logical,    intent(in) :: master
   logical,    intent(in) :: slave
   integer(4), intent(in) :: nproc
   integer(4), intent(in) :: unit

   integer(4)      :: ierr, procnamelen, n
   integer(4)      :: istatus(mpi_status_size), sid
   character(50), allocatable :: procname(:)

   allocate(procname(0:nproc))
   procname = ""

! ... get hostnames where the processes are residing

   procname(myid)='                                                  '
   call mpi_get_processor_name(procname(myid),procnamelen,ierr)

! ... slaves send their id to master

   if (slave) call mpi_send(myid,1,mpi_integer4,0,myid,mpi_comm_world,ierr)
   if (slave) call mpi_send(procname(myid),50,mpi_character,0,myid,mpi_comm_world,ierr)

! ... and master receive the id of the slaves

   if (master) then
      write(unit,'(a,10x,a)') 'master: ', procname(myid)
      do n=1,nproc-1
         call mpi_recv(sid,1,mpi_integer4,n,n,mpi_comm_world,istatus,ierr)
         call mpi_recv(procname(myid),50,mpi_character,n,n,mpi_comm_world,istatus,ierr)
         write(unit,'(a,i4,a,a)') 'slave ', sid, ' is up: ', procname(myid)
      end do
      write(unit,*)
   end if
   if (ierr/=mpi_success) call par_error('par_handshake',ierr)
   if (master) call FileFlush(unit)
   deallocate(procname)

end subroutine par_handshake

!************************************************************************
!> \page parallel parallel.F90
!! **par_timing**
!! *parallel timing*
!************************************************************************


subroutine par_timing(whattodo,master,nproc,unit)

   use ParallelModule
   implicit none

   character(*), intent(in) :: whattodo           ! 'start' and 'stop'
   logical,      intent(in) :: master
   integer(4),   intent(in) :: nproc
   integer(4),   intent(in) :: unit

   integer(4), save         :: mpi_t0, mpi_t1     ! parallel wall-clock timing

   if(master) then
      if (whattodo=='start') then
         mpi_t0=mpi_wtime()
      else if (whattodo=='stop') then
         mpi_t1=mpi_wtime() - mpi_t0
         call WriteHead(2,'mpi timing',unit)
         write(unit,'(a,f12.3,i12)') 'total wall time since start (h/s)  =',mpi_t1/3600.0, mpi_t1
         write(unit,'(a,i12)')       'number of processes                =',nproc
      end if
   end if

end subroutine par_timing

!************************************************************************
!> \page parallel parallel.F90
!! **par_error**
!! *error handler for par routines*
!************************************************************************


subroutine par_error(caller,ierr)

   use ParallelModule
   implicit none

   character(*), intent(in) :: caller
   integer(4)  , intent(in) :: ierr

   if (ierr==mpi_err_buffer   ) call Stop(caller,'invalid buffer !'           ,6)
   if (ierr==mpi_err_count    ) call Stop(caller,'invalid icount !'           ,6)
   if (ierr==mpi_err_type     ) call Stop(caller,'invalid type !'             ,6)
   if (ierr==mpi_err_tag      ) call Stop(caller,'invalid tag !'              ,6)
   if (ierr==mpi_err_comm     ) call Stop(caller,'invalid communicator !'     ,6)
   if (ierr==mpi_err_rank     ) call Stop(caller,'invalid rank !'             ,6)
   if (ierr==mpi_err_root     ) call Stop(caller,'invalid root !'             ,6)
   if (ierr==mpi_err_group    ) call Stop(caller,'invalid group !'            ,6)
   if (ierr==mpi_err_op       ) call Stop(caller,'invalid operation !'        ,6)
   if (ierr==mpi_err_topology ) call Stop(caller,'invalid topology !'         ,6)
   if (ierr==mpi_err_dims     ) call Stop(caller,'invalid dimen. arg. !'      ,6)
   if (ierr==mpi_err_arg      ) call Stop(caller,'invalid argument !'         ,6)
   if (ierr==mpi_err_unknown  ) call Stop(caller,'unknown error !'            ,6)
   if (ierr==mpi_err_truncate ) call Stop(caller,'message truncated on recv !',6)
   if (ierr==mpi_err_other    ) call Stop(caller,'other error !'              ,6)
   if (ierr==mpi_err_intern   ) call Stop(caller,'internal error code !'      ,6)
   if (ierr==mpi_err_in_status) call Stop(caller,'error value in status !'    ,6)
   if (ierr==mpi_err_request  ) call Stop(caller,'illegal mpi_request hand. !',6)
   if (ierr==mpi_err_lastcode ) call Stop(caller,'lastcode'                   ,6)
                                call Stop(caller,'no error found !'           ,6)

end subroutine par_error

#endif
