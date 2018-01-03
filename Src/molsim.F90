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
!*     MolsimDriver                                                     *
!*                                                                      *
!************************************************************************

! ... driver of the molsim program

program MolsimDriver

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='MolsimDriver'

! ... set nproc, myid, master, and slave

   nproc  = 1
   myid   = 0
   master = .true.
   slave  = .false.
#if defined (_PAR_)
   call par_initialize
   call par_comm_size(nproc)
   if (master .and. (nproc > mnproc)) call Stop(txroutine, 'nproc > mnproc',uout)
   call par_comm_rank(myid, master, slave)
   if (master) write(*,'(a,i5,/)') 'Molsim: start of parallel run, number of processes =' ,nproc
   call par_handshake(myid, master, slave, nproc, ustdout)
   call par_timing('start', master, nproc, uout)
#endif

! ............... initiate and read input data ..............

   call IOMolsim(iReadInput)
   call IOSystem(iReadInput)
   call IOScale(iReadInput)
   call Particle(iReadInput)
   call PotentialDriver(iReadInput)
   if (.not.lmix) call Coordinate(iReadInput)
   call MolsimDriverSub(iReadInput)

! ............... process and write input data ..............

   call IOMolsim(iWriteInput)
   call IOSystem(iWriteInput)
   call IOScale(iWriteInput)
   call Particle(iWriteInput)
   call PotentialDriver(iWriteInput)
   call IOCnf('open')
   if (.not.lmix) call Coordinate(iWriteInput)
   call MolsimDriverSub(iWriteInput)

! ............... initiate simulation ................

   call IOMolsim(iBeforeSimulation)
   call IOSystem(iBeforeSimulation)
   call PotentialDriver(iBeforeSimulation)
   call MolsimDriverSub(iBeforeSimulation)
   call IOCnf('close')

   do istep1 = nstep1beg, nstep1

! ............... initiate one macrostep  ................

      if (maxcpu > 0) call CpuLeft(maxcpu, uout)
      call IOMolsim(iBeforeMacrostep)
      call IOSystem(iBeforeMacrostep)
      call MolsimDriverSub(iBeforeMacrostep)

      do istep2 = 1, nstep2
         istep = nstep2*(istep1-1)+istep2

! .............. perform one time step/mass ............

         call IOMolsim(iSimulationStep)
         call MolsimDriverSub(iSimulationStep)

      end do

! .............. average over one macrostep ...............

      call IOCnf('open')
      call IOMolsim(iAfterMacrostep)
      call PotentialDriver(iAfterMacrostep)
      call Coordinate(iAfterMacrostep)
      call MolsimDriverSub(iAfterMacrostep)
      call IOCnf('close')

   end do

! ............... make grand averages ................

   call IOMolsim(iAfterSimulation)
   call IOSystem(iAfterSimulation)
   call PotentialDriver(iAfterSimulation)
   call Coordinate(iAfterSimulation)
   call MolsimDriverSub(iAfterSimulation)
   call Particle(iAfterSimulation)          ! has to be after call of Image

   if (ilist/= 0 .and. master) close(ulist)

   if (master) call WriteHead(2, 'timing statistics', uout)
   if (ltime .and. master) call CpuAdd('write', ' ', 0, uout)
   if (master) call CpuTot(uout)

#if defined (_PAR_)
   call par_timing('stop', master, nproc, uout)
#endif

   if (master) close(uout)

#if defined (_PAR_)
      call par_finalize
#endif

contains

!........................................................................

subroutine MolsimDriverSub(iStage)
   integer(4), intent(in) :: iStage
   if (lsim) then
      if (lmd)    call MDDriver(iStage)
      if (lmc)    call MCDriver(iStage)
      if (lmcall) call MCAllDriver(iStage)
      if (lbd)    call BDDriver(iStage)
   else if (lmix) then
      call MixedDriver(iStage)
   end if
   call NListDriver(iStage)
   if (lsim) then
      if (lcont)  call ControlAver(iStage)
      if (laver)  call MainAver(iStage)
      if (lti)    call ThermoInteg(iStage)
      if (ldist)  call DistFunc(iStage)
   end if
   if (ldump)    call DumpDriver(iStage)
   if (lgroup)   call Group(iStage)
   if (lstatic)  call StaticDriver(iStage)
   if (ldynamic) call DynamicDriver(iStage)
   if (limage)   call ImageDriver(iStage)
end subroutine MolsimDriverSub

!........................................................................

end program MolsimDriver

!************************************************************************
!*                                                                      *
!*     IOMolsim                                                         *
!*                                                                      *
!************************************************************************

! ... perform mixed tasks, mainly i/o

subroutine IOMolsim(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOMolsim'
   character(4) :: txistep1

   character(len=128) :: arg
   integer(4) :: i

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

      !parse command line arguments
      if (command_argument_count() == 0) then
         call Stop(txroutine, 'No arguments provided', ustdout)
      end if

      do i = 1, command_argument_count()
         call get_command_argument(i, arg)

         select case (arg)
         case ('-v', '--version', '-V')
            if(master) then
               write(*,'(a)') txVersionDate
# ifdef _TEST_
               write(*,'(a)') "mode = test"
# elif _NORMAL_
               write(*,'(a)') "mode = normal"
# elif _WARN_
               write(*,'(a)') "mode = warn"
# elif _DEBUG_
               write(*,'(a)') "mode = debug"
# elif _QUICK_
               write(*,'(a)') "mode = quick"
# elif _GRPOF_
               write(*,'(a)') "mode = gprof"
# endif
            end if
            stop 0
         end select
      end do

      !set filenames from the last argument
      call get_command_argument(command_argument_count(), project)
      fin   = trim(adjustl(project))//'.in'
      fout  = trim(adjustl(project))//'.out'
      fcnf  = trim(adjustl(project))//'.cnf'
      flist = trim(adjustl(project))//'.list'
      fuser = trim(adjustl(project))//'.user'
      fwrl  = trim(adjustl(project))//'.wrl'
      fvtf  = trim(adjustl(project))//'.vtf'
      ftcl  = trim(adjustl(project))//'.tcl'
      fgroup= trim(adjustl(project))//'.group'
      fpos  = trim(adjustl(project))//'.pos'
      fori  = trim(adjustl(project))//'.ori'
      fliv  = trim(adjustl(project))//'.liv'
      fanv  = trim(adjustl(project))//'.anv'
      ffor  = trim(adjustl(project))//'.for'
      ftor  = trim(adjustl(project))//'.tor'
      fidm  = trim(adjustl(project))//'.idm'
      flaz  = trim(adjustl(project))//'.laz'
      futot = trim(adjustl(project))//'.utot'

! ... open FIN and FOUT

      call FileOpen(uin, fin , 'form/noread')
      if (master) call FileOpen(uout, fout, 'form/noread')

   case (iWriteInput)

      beta = sclene/(GasConstant*tempst*scltem)

! ... open FLIST

      if (master .and. ilist/= 0) call FileOpen(ulist, flist, 'form/noread')

! ... check unsupported conditions

      if (lbd .and. lpolyatom) call Stop(txroutine, 'bd .and. lpolyatom is not allowed', uout)
      if (lmc .and. lpolarization) call Stop(txroutine, 'lmc .and. lpolarization is not allowed', uout)

! ............... write input data ................

      if (master) then
         call WriteFront('Molsim', 'integrated md/mc/bd simulation program belonging to the molsim package', &
                         txVersionDate, txAuthor,'&
                         & Anna Akinchina, Fredrik Carlsson, Samuel Edgecombe, Yoshikatsu Hayashi \&
                         & Pascal Hebbeker, Cornelius Hofzumahaus, Niklas Källrot, Björn Linse, Vladimir Lobaskin \&
                         & Thomas M. Nymand, Alberto Pais, Jurij Rescic, Stefanie Schneider, Marie Skepö \&
                         & Joakim Stenhammar, Anders Wallqvist, Jos van Rijssel, Erik Wernersson, Per-Olof Åstrand &
                         &', uout)
         call WriteDateTime(uout)
         call WriteHead(1, 'input data', uout)
         call WriteHead(2, 'general information', uout)
      end if

   case (iBeforeSimulation)

      close(uin)

      call WarnPartOutsideBox(1, np)
      if (lradatbox) Call WarnAtomOutsideBox(1, na)
#if defined (ALARIK_INTEL)
#else
      if (np < 10000) call WarnHCOverlap(1, np)
#endif
      if (itest == 1)   call TestSimulation

      if (master) call FileFlush(uout)

   case (iBeforeMacrostep)

      if (nstep1 == 0) call Stop(txroutine, 'nstep1 == 0', uout)
      if (lsim .and. master) write(uout, *)
      if (lsim .and. master) call WriteDateTime(uout)
      if (lsim .and. master) call CpuTot(uout)
      if (lsim .and. master) then
         write(txistep1,'(i4)') istep1
         call WriteHead(1, 'result of macrostep '//txistep1, uout)
      end if

   case (iSimulationStep)

      if (master) call FileFlush(uout)

   case (iAfterMacrostep)

      call WarnPartOutsideBox(1, np)
      if (lradatbox) Call WarnAtomOutsideBox(1, na)
#if defined (ALARIK_INTEL)
#else
      if (np < 10000) call WarnHCOverlap(1, np)
#endif
      if (master) call FileFlush(uout)
      if (lsim .and. master) call WriteDateTime(ustdout)
      if (lsim .and. master) write(ustdout,'(a,i4,a,/)') 'macrostep', istep1 , ' is completed'
      if (lsim .and. master) call FileFlush(ustdout)

   case (iAfterSimulation)

      if (master) call WriteHead(1, 'final result', uout)
      if (master) call WriteDateTime(uout)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine IOMolsim

!************************************************************************
!*                                                                      *
!*     IOSystem                                                         *
!*                                                                      *
!************************************************************************

! ... perform i/o on general system variables

subroutine IOSystem(iStage)

   use MolModule
   use Random_Module, only: ix, iy, k4b
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOSystem'
   real(8) :: SecondsSinceStart
   integer(k4b)  :: ixseed, iyseed, iseedtmp
   real(8)       :: Random
   logical, save :: luseXYseed

   namelist /nmlSystem/ txtitle,                                                               &
                        txmode,                                                                &
                        txmethod, txensemb, txbc, txstart, txuser,                             &
                        boxlen, cellside, sphrad, ellrad, cylrad, cyllen, lenscl,              &
                        temp, prsr,                                                            &
                        nstep1, nstep2,                                                        &
                        iseed, ixseed, iyseed, luseXYseed, maxcpu,                             &
                        lcont, laver, lti,   ldist, ldump, lgroup, lstatic, ldynamic, limage,  &
                        itest, ipart, iatom, iaver, ishow, iplot,  ilist,                      &
                        ltrace, lblockaver,                                                    &
                        ltime


   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      txmode     = 'simulation'
      txuser     = ' '
      lenscl     = One
      maxcpu     = 0
      lcont      = .false.
      laver      = .false.
      lti        = .false.
      ldist      = .false.
      ldump      = .false.
      lgroup     = .false.
      lstatic    = .false.
      ldynamic   = .false.
      limage     = .false.
      itest      = 0
      ipart      = 0
      iatom      = 0
      iaver      = 0
      ishow      = 0
      iplot      = 0
      ilist      = 0
      ltrace     = .false.
      lblockaver = .false.
      ltime      = .true.
      ixseed     = -1
      iyseed     = -1
      luseXYseed = .false.

      rewind(uin)
      read(uin,nmlSystem)

      call LowerCase(txmethod)
      call LowerCase(txensemb)
      call LowerCase(txbc)
      call LowerCase(txstart)
      call LowerCase(txuser)


! ... set lsim and lana

      lsim =.false.
      lana =.false.
      lmix =.false.
      if (txmode == 'simulation') then
         lsim =.true.            ! simulation mode
      else if (txmode == 'analysis') then
         lana =.true.            ! analysis mode
      else if (txmode == 'mixed') then
         lmix =.true.            ! mixed mode
      else
         call stop(txroutine, 'error in txmode', uout)
      end if

! ... set method and ensemble variables

      lnve    = .false.
      lmd     = .false.
      lmc     = .false.
      lmcall  = .false.
      lbd     = .false.
      lnvt    = .false.
      lntp    = .false.
      lmvt    = .false.
      if (txmethod == 'md') then
         lmd =.true.
         if (txensemb == 'nve') then
            lnve =.true.
         else
            call Stop(txroutine, 'md: unsupported ensemble', uout)
         end if
      else if (txmethod == 'mc') then
         lmc =.true.
         if ((txensemb == 'nvt') .or. (txensemb == 'ntv')) then
            lnvt =.true.
         else if ((txensemb == 'ntp') .or. (txensemb == 'npt')) then
            lntp =.true.
         else if ((txensemb == 'mvt') .or. (txensemb == 'mtv')) then
            lmvt =.true.
         else
            call Stop(txroutine, 'mc: unsupported txensemb', uout)
         end if
      else if (txmethod == 'mcall') then
         lmcall =.true.
         if ((txensemb == 'nvt') .or. (txensemb == 'ntv')) then
            lnvt =.true.
         else if ((txensemb == 'npt') .or. (txensemb == 'ntp')) then
            lntp =.true.
         else if ((txensemb == 'mvt') .or. (txensemb == 'mtv')) then
            lmvt =.true.
         else
            call Stop(txroutine, 'mcall: unsupported txensemb', uout)
         end if
      else if (txmethod == 'bd') then
         lbd =.true.
         if ((txensemb == 'nvt') .or. (txensemb == 'ntv')) then
            lnvt =.true.
         else
            call Stop(txroutine, 'bd: unsupported txensemb', uout)
         end if
      else
         write(*,*) 'method',txmethod
         call Stop(txroutine, 'unsupported txmethod', uout)
      end if

! ... specific for lana: force ldump and lgroup to be true

      if (lana) then
         ldump  =.true.
         lgroup =.true.
      end if

! ... specific for lmix

      if (lmix) then
         lmd     = .false.
         lmc     = .false.
         lmcall  = .false.
         lbd     = .false.
         iatom = 1
         ipart = 1
      end if

! ... set box related variables

      call SetBoxParam

      if (lntp .and. (.not.lbcbox)) call Stop(txroutine,'ntp and non-cubic geometry not supported',uout)

! ... set some other variables

      nstep  = nstep1*nstep2
      tempst = temp
      prsrst = prsr
      volst  = vol

      if (luseXYseed) then ! the ix and iy are to be used
         ! initialize with a negative numberthe random number generator to set the value of am (see Random function)
         iseedtmp = -abs(iseed)
         raux = Random(iseedtmp)
         ! set the ix and iy values as requested in the input file
         ix = ixseed
         iy = iyseed
      else
         if (iseed <= 0) iseed = int(1.0e+6*SecondsSinceStart())
         ! the first seed should be negative, as the random seed generator is only initialized properly when the passed seed is
         ! negative. Note that the seed returned is positive an not changing on subsequent calls of Random(iseed) when iseed > 0:
         iseed = -abs(iseed)
      end if


   case (iWriteInput)

      if (txuser == 'jos') call JosUser(2)

      if (master) then

!        if (ldipole.and.lmc) call Stop(txroutine, 'ldipole and mc is not supported', uout)
         if (ldipole.and.lbd) call Stop(txroutine, 'iadip and bd is not supported', uout)
         if (lpolarization.and.lmc) call Stop(txroutine, 'lpolarization and mc is not supported', uout)
         if (lpolarization.and.lbd) call Stop(txroutine, 'lpolarization and bd is not supported', uout)
         if (ltrace .and. nstep2 > 100) call Stop(txroutine, 'ltrace .and. nstep2 > 100' , uout)

         write(uout,'()')
         write(uout,'(a,a)') 'title   : ', txtitle
         write(uout,'()')
         write(uout,'(a,a)') 'key data: mode           : ', txmode
         write(uout,'(a,a)') '          method         : ', txmethod
         write(uout,'(a,a)') '          ensemble       : ', txensemb
         write(uout,'(a,a)') '          boundary cond. : ', txbc
         write(uout,'(a,a)') '          start          : ', txstart
         if (txuser(1:1) /= ' ') &
         write(uout,'(a,a)') '          txuser         : ', txuser
         write(uout,'()')
         write(uout,'()')
         write(uout,'(a,t30,i10,a,i10,a,i10)') 'number of time steps/passes', nstep1, '  *', nstep2, '  = ', nstep
         write(uout,'()')
         write(uout,'(a,t35,f10.3)') 'temperature                    = ', temp
         write(uout,'(a,t35,f10.3)') 'pressure                       = ', prsr
         if (lbcbox) then
            write(uout,'(a,t35,3(f10.3,2x))') 'box lengths (x,y,z)            = ', boxlen
         else if (lbcrd) then
            write(uout,'(a,t35,f10.3)') 'side length rhombic dodecaheron= ', cellside
         else if (lbcto) then
            write(uout,'(a,t35,f10.3)') 'side length truncated octaheron= ', cellside
         else if (lbcsph) then
            write(uout,'(a,t35,f10.3)') 'spherical cell radius          = ', sphrad
         else if (lbcell) then
            write(uout,'(a,t35,3(f10.3,2x))') 'ellipsoidal radii (x,y,z)      =  ', ellrad
         else if (lbccyl) then
            write(uout,'(a,t35,f10.3)') 'cylindrical cell radius        = ', cylrad
            write(uout,'(a,t35,f10.3)') 'cylindrical cell length        = ', cyllen
         end if
         write(uout,'(a,t35,e10.3)') 'volume                         = ', vol
         write(uout,'()')
         if(luseXYseed) then
            write(uout,'(a)')         'using iseed, ixseed and iyseed from the input'
            write(uout,'(a,t35,i12)') 'iseed of random generator      = ', iseed
            write(uout,'(a,t35,i12)') 'ixseed of random generator     = ', ix
            write(uout,'(a,t35,i12)') 'iyseed of random generator     = ', iy
         else
            write(uout,'(a,t35,i10)') 'seed of random generator       = ', abs(iseed)
         end if
         write(uout,'()')
         if (maxcpu > 0) then
            write(uout,'(a,t35,i10)') 'maximum cpu time (seconds)     = ', maxcpu
            write(uout,'()')
         end if

         write(uout,'(a)') 'optional subroutine used'
         write(uout,'(a)') '------------------------'
         if (lcont)    write(uout,'(a)') '   ControlAver'
         if (laver)    write(uout,'(a)') '   MainAver'
         if (lti)      write(uout,'(a)') '   ThermoInteg'
         if (ldist)    write(uout,'(a)') '   DistFunc'
         if (ldump)    write(uout,'(a)') '   Dump'
         if (lgroup)   write(uout,'(a)') '   Group'
         if (lstatic)  write(uout,'(a)') '   Static'
         if (ldynamic) write(uout,'(a)') '   Dynamic'
         if (limage)   write(uout,'(a)') '   Image'

         write(uout,'()')
         write(uout,'(a)') 'parameters controlling the output interval (0 = no)'
         write(uout,'(a)') '-------------------------------------------------'
         write(uout,'(a,i10)') '   itest     ', itest
         write(uout,'(a,i10)') '   ipart     ', ipart
         write(uout,'(a,i10)') '   iatom     ', iatom
         write(uout,'(a,i10)') '   iaver     ', iaver
         write(uout,'(a,i10)') '   ishow     ', ishow
         write(uout,'(a,i10)') '   iplot     ', iplot
         write(uout,'(a,i10)') '   ilist     ', ilist
         write(uout,'()')

         if (ltrace) then
            write(uout,'(a)') 'trace data is written on trace.master.data'
            write(uout,'()')
         end if
         if (lblockaver) then
            call SetlBlockAver(lblockaver)
            write(uout,'(a)') 'blockavering data is written on blockaver.data'
            write(uout,'()')
         end if

         write(uout,'(a)') 'external units'
         write(uout,'(a)') '--------------'
         write(uout,'(a,t15,a,t35,a)')  'number', 'generic name', 'file name'
         write(uout,'(a,t15,a,t35,a)')  '------', '------------', '---------'
         write(uout,'(i4,t20,a,t35,a)')                 uin,    'fin  ', trim(fin)
         write(uout,'(i4,t20,a,t35,a)')                 uout,   'fout ', trim(fout)
         if (lsim) write(uout,'(i4,t20,a,t35,a)')       ucnf,   'fcnf ', trim(fcnf)
         if (ilist/= 0) write(uout,'(i4,t20,a,t35,a)')  ulist,  'flist', trim(flist)
         if (limage) write(uout,'(i4,t20,a,t35,a)')     uwrl,   'fwrl ', trim(fwrl)
         if (limage) write(uout,'(i4,t20,a,t35,a)')     uvtf,   'fvtf ', trim(fvtf)
         if (limage) write(uout,'(i4,t20,a,t35,a)')     utcl,   'ftcl ', trim(ftcl)    ! TO BE FIXED
      end if

   case (iAfterSimulation)

      if (master) then

         write(uout,'()')
         write(uout,'(a,a)') 'title   : ', txtitle
         write(uout,'()')
         write(uout,'(a,a)') 'key data: mode           : ', txmode
         write(uout,'(a,a)') '          method         : ', txmethod
         write(uout,'(a,a)') '          ensemble       : ', txensemb
         write(uout,'(a,a)') '          boundary cond. : ', txbc
         write(uout,'(a,a)') '          start          : ', txstart
         if (txuser(1:1) /= ' ') &
         write(uout,'(a,a)') '          txuser         : ', txuser
         write(uout,'()')
         write(uout,'()')
         write(uout,'(a,t30,i10,a,i10,a,i10)') 'number of time steps/passes', nstep1, '  *', nstep2, '  = ', nstep
         write(uout,'()')
         write(uout,'(a,t35,f10.3)') 'temperature                    = ', temp
         write(uout,'(a,t35,f10.3)') 'pressure                       = ', prsr
         if (lbcbox) then
            write(uout,'(a,t35,3(f10.3,2x))') 'box lengths (x,y,z)            = ', boxlen
         else if (lbcrd) then
            write(uout,'(a,t35,f10.3)') 'side length rhombic dodecaheron= ', cellside
         else if (lbcto) then
            write(uout,'(a,t35,f10.3)') 'side length truncated octaheron= ', cellside
         else if (lbcsph) then
            write(uout,'(a,t35,f10.3)') 'spherical cell radius          = ', sphrad
         else if (lbcell) then
            write(uout,'(a,t35,3(f10.3,2x))') 'ellipsoidal radii (x,y,z)      =  ', ellrad
         else if (lbccyl) then
            write(uout,'(a,t35,f10.3)') 'cylindrical cell radius        = ', cylrad
            write(uout,'(a,t35,f10.3)') 'cylindrical cell length        = ', cyllen
         end if
         write(uout,'(a,t35,e10.3)') 'volume                         = ', vol
         write(uout,'()')
         write(uout,'(a,t35,i12)') 'iseed of random generator      = ', iseed
         write(uout,'(a,t35,i12)') 'ixseed of random generator     = ', ix
         write(uout,'(a,t35,i12)') 'iyseed of random generator     = ', iy
         write(uout,'()')

      end if

   case (iBeforeMacrostep)

      if (txuser == 'sim_annealing') then
          call WriteHead(2, 'simulated annealing', uout)
          temp = tempst*0.5**(istep1-1)
          temp = tempst*0.8**(istep1-1)
          beta = sclene/(GasConstant*temp*scltem)
          write(uout,'(a,t35,f10.3)') 'temperature                    = ', temp
          write(uout,'(a,t35,f10.3)') 'beta                           = ', beta
      end if

   end select

end subroutine IOSystem

!************************************************************************
!*                                                                      *
!*     IOScale                                                          *
!*                                                                      *
!************************************************************************

! ... perform i/o on scaling variables

!     q(si) = q(internal)*(scaling factor)

subroutine IOScale(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IOScale'

   namelist /nmlScale/ scllen, sclmas, scltem, scltim, sclene, sclhca, sclpre, scldif, sclang

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      scllen = 1.0d-10
      sclmas = 1.0d-03
      scltem = 1.0d0
      scltim = 1.0d-12
      sclene = 1.0d+03
      sclhca = 1.0d0
      sclpre = 1.0d+06
      scldif = 1.0d-09
      sclang = DegToRad

      rewind(uin)
      read(uin,nmlScale)

      sclvol = scllen**3
      sclvel = scllen/scltim
      sclacc = sclvel/scltim
      sclfor = sclacc*sclmas
      sclmin = sclmas*scllen**2

      Epsi0FourPi = ECh**2/(FourPi*Epsi0*scllen)*AvNo/sclene

   case (iWriteInput)

      if (master) then
         call WriteHead(2, 'units', uout)
         write(uout,'(a)') 'quantity                unit'
         write(uout,'(a)') '--------                ----'
         write(uout,'(a,es12.3,a)') 'length           ', scllen, ' m'
         write(uout,'(a,es12.3,a)') 'mass             ', sclmas, ' kg/mol'
         write(uout,'(a,es12.3,a)') 'temperature      ', scltem, ' k'
         write(uout,'(a,es12.3,a)') 'time             ', scltim, ' s'
         write(uout,'(a,es12.3,a)') 'energy           ', sclene, ' j/mol'
         write(uout,'(a,es12.3,a)') 'heat capacity    ', sclhca, ' j/mol*k'
         write(uout,'(a,es12.3,a)') 'diff. coeff.     ', scldif, ' m**2/s'
         write(uout,'(a,es12.3,a)') 'pressure         ', sclpre, ' pa'
         write(uout,'(a,es12.3,a)') 'volume           ', sclvol, ' m**3'
         write(uout,'(a,es12.3,a)') 'linear vel.      ', sclvel, ' m/s'
         write(uout,'(a,es12.3,a)') 'linear acc.      ', sclacc, ' m/s**2'
         write(uout,'(a,es12.3,a)') 'force            ', sclfor, ' n/mol'
         write(uout,'(a,es12.3,a)') 'moment of inertia', sclmin, ' kg*m**2/mol'
         write(uout,'(a,es12.3,a)') 'angle            ', sclang, ' rad'
      end if

   end select

end subroutine IOScale

!************************************************************************
!*                                                                      *
!*     IOCnf                                                            *
!*                                                                      *
!************************************************************************

! ... perform i/o on the configuration file

subroutine IOCnf(str)

   use MolModule
   use Random_Module, only: am, ix, iy, k4b
   use ISO_FORTRAN_ENV
   implicit none

   character(*) :: str

   character(40), parameter :: txroutine ='IOCnf'
   integer(4)   :: ip, npread
   logical      :: GetlSetVel
   integer(4)   :: ierr
   integer(k4b) :: ik4bAux

!   write(*,*) 'IOCNF start: str= ',str

   if (str == 'open') then

      if (master) call FileOpen(ucnf, fcnf, 'unform/noread')

   else if (str == 'read') then

      if (master) then
         rewind(ucnf)
         if (txstart == 'zero') then
!           read(ucnf) nstep1done, iseed, boxlen               ! -1999-12-21, 2003-09-30--2005-09-02
!           read(ucnf) nstep1done, iseed, raux, raux, raux     !  1999-12-21--2003-09-30
!           read(ucnf) nstep1done, iaux, boxlen                !  2005-09-02--2007-03-16
            if (lntp) then
               read(ucnf, iostat=ierr) nstep1done, iaux, raux, ik4bAux, ik4bAux, boxlen !  2017-08-22-
               if(ierr .ne. 0) then ! try old format
                  backspace(ucnf)
                  read(ucnf) nstep1done, iaux, boxlen          !  2007-03-16--2017-08-22
               end if
            else
               ! read(ucnf) nstep1done, iaux, raux, raux, raux !  2007-03-16--2017-08-22
               read(ucnf, iostat=ierr) nstep1done              !  2017-08-22-
            end if
         else if (txstart == 'continue') then
            read(ucnf) nstep1done, iseed, am, ix, iy, boxlen   ! 2017-01-05-

            !read(ucnf) nstep1done, iseed, boxlen              ! -1999-12-21, 2003-09-30-2017-01-05
!           read(ucnf) nstep1done, iseed, raux, raux, raux     !  1999-12-21--2003-09-30
         end if
         read(ucnf) npread, iaux, raux, raux, laux, laux
         if (lmvt) np = npread
         read(ucnf) (ro(1:3,ip),qua(0:3,ip),ip = 1,np)
         if (lclink) read(ucnf) nbondcl(1:np), bondcl(1:maxvalnbondcl,1:np)
         if (lweakcharge) then
            read(ucnf) laz(1:np)
            where (.not. laz(1:np)) az(1:np) = Zero
         end if
         if (lmd .and. .not.GetlSetVel()) then
            read(ucnf,end = 998) (rod(1:3,ip),quad(0:3,ip),ip = 1,np)
 998        continue
            call QuaVelToAngVel(np, 1, np, qua, quad, angvelo)
          end if
      end if

#if defined (_PAR_)
      call par_bc_int (nstep1done)
      call par_bc_int (iseed)
      call par_bc_real (am)
      call par_bc_int (ix)
      call par_bc_int (iy)
      call par_bc_reals(boxlen    , 3   )
      call par_bc_reals(ro        , 3*np)
      call par_bc_reals(qua       , 4*np)
      if (lclink) then
         call par_bc_ints(nbondcl ,   np)
         call par_bc_ints(bondcl  , maxvalnbondcl*np)
      end if
      if (lweakcharge) then
         call par_bc_logicals(laz ,   np)
         call par_bc_reals(az     ,   np)
      end if
      if (lmd .and. .not.GetlSetVel()) then
         call par_bc_reals(rod    , 3*np)
         call par_bc_reals(quad   , 4*np)
         call par_bc_reals(angvelo, 3*np)
      end if
#endif

      call QuaToOri(np, 1, np, qua, ori)
      call SetBoxParam
      if (lewald) call EwaldSetup
      if (lmvt) then
         if (npt > 1) call Stop(txroutine, 'mvt .and. npt > 1', uout)
         nppt(1) = np
         call SetObjectParam1
         call SetObjectParam2
                                   ! neighbour list is generated later
      end if

   else if (str == 'write' .and. master) then

      rewind(ucnf)
      if (lmc .or. lmcall) call OriToQua(np, 1, np, ori, qua)
      write(ucnf) istep1, iseed, am, ix, iy, boxlen
      write(ucnf) np, iaux, raux, raux, laux, laux
      write(ucnf) (ro(1:3,ip),qua(0:3,ip),ip = 1,np)
      if (lclink) write(ucnf) nbondcl(1:np), bondcl(1:maxvalnbondcl,1:np)
      if (lweakcharge) write(ucnf) laz(1:np)
      if (lmd) write(ucnf) (rod(1:3,ip),quad(0:3,ip),ip = 1,np)
      call FileFlush(ucnf)

   else if (str == 'close' .and. master) then

      close(ucnf, iostat = ierr)

   end if

end subroutine IOCnf

!************************************************************************
!*                                                                      *
!*     ControlAver                                                      *
!*                                                                      *
!************************************************************************

! ... calculate control quantities

subroutine ControlAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ControlAver'
   integer(4) :: ipt

   if (slave) return    ! master only

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iAfterMacrostep, iAfterSimulation)
      call WriteHead(2, 'control statistics', uout)
      write(uout,'(a,t30,a,6(10x,a))') 'quantity', '     ', txpt(1:npt)
      write(uout,'(a,t30,a,6(10x,a))') '--------', '     ', ('----------',ipt = 1,npt)
   end select
   if (lmd)    call MDAver(iStage)
   if (lmcall) call MDAver(iStage)
   call DispAver(iStage)
   call OriOrderAver(iStage)
   call PosOriAver(iStage)

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine ControlAver

!************************************************************************
!*                                                                      *
!*     MDAver                                                           *
!*                                                                      *
!************************************************************************

! ... calculate averages of md specific variables

subroutine MDAver(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MDAver'
   integer(4), save :: nsamp1, nsamp2
   real(8), allocatable, save :: for2s2(:),   for2s1(:)   ! force
   real(8), allocatable, save :: tor2s2(:,:), tor2s1(:,:) ! torque
   real(8), allocatable, save :: lims2(:,:),  lims1(:,:)  ! linear moment
   real(8), allocatable, save :: anms2(:,:),  anms1(:,:)  ! angular moment
   real(8), allocatable, save :: ttras2(:),   ttras1(:)   ! translational temp
   real(8), allocatable, save :: trots2(:),   trots1(:)   ! rotational temp
   integer(4) :: ip, ipt
   real(8)    :: norm, fac

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      if(.not.allocated(for2s2)) then
         allocate( for2s2(npt), for2s1(npt), tor2s2(3,npt), tor2s1(3,npt), &
         lims2(3,npt), lims1(3,npt), anms2(3,npt), anms1(3,npt), ttras2(npt), ttras1(npt), trots2(npt), trots1(npt) )
         for2s2 = 0.0E+00
         for2s1 = 0.0E+00
         tor2s2 = 0.0E+00
         tor2s1 = 0.0E+00
         lims2 = 0.0E+00
         lims1 = 0.0E+00
         anms2 = 0.0E+00
         anms1 = 0.0E+00
         ttras2 = 0.0E+00
         ttras1 = 0.0E+00
         trots2 = 0.0E+00
         trots1 = 0.0E+00
      end if

      nsamp1 = 0
      for2s1 = Zero
      tor2s1 = Zero
      lims1  = Zero
      anms1  = Zero
      ttras1 = Zero
      trots1 = Zero
      if (txstart == 'continue') read(ucnf) nsamp1, for2s1, tor2s1, lims1, anms1, ttras1, trots1

   case (iBeforeMacrostep)

      nsamp2 = 0
      for2s2 = Zero
      tor2s2 = Zero
      lims2  = Zero
      anms2  = Zero
      ttras2 = Zero
      trots2 = Zero

   case (iSimulationStep)

      nsamp2 = nsamp2 + 1
      do ip = 1, np
         ipt = iptpn(ip)
         for2s2(ipt)   = for2s2(ipt)  + sum(forceo(1:3,ip)**2)
         tor2s2(1,ipt) = tor2s2(1,ipt)+ (dot_product(torqueo(1:3,ip),ori(1:3,1,ip)))**2
         tor2s2(2,ipt) = tor2s2(2,ipt)+ (dot_product(torqueo(1:3,ip),ori(1:3,2,ip)))**2
         tor2s2(3,ipt) = tor2s2(3,ipt)+ (dot_product(torqueo(1:3,ip),ori(1:3,3,ip)))**2
         lims2(1,ipt)  = lims2(1,ipt) + massp(ip)*rod(1,ip)
         lims2(2,ipt)  = lims2(2,ipt) + massp(ip)*rod(2,ip)
         lims2(3,ipt)  = lims2(3,ipt) + massp(ip)*rod(3,ip)
         anms2(1,ipt)  = anms2(1,ipt) + massp(ip)*(ro(2,ip)*rod(3,ip)-ro(3,ip)*rod(2,ip)) &
                                      + dot_product(ori(1,1:3,ip)*momp(1:3,ip),angvelo(1:3,ip))
         anms2(2,ipt)  = anms2(2,ipt) + massp(ip)*(ro(3,ip)*rod(1,ip)-ro(1,ip)*rod(3,ip)) &
                                      + dot_product(ori(2,1:3,ip)*momp(1:3,ip),angvelo(1:3,ip))
         anms2(3,ipt)  = anms2(3,ipt) + massp(ip)*(ro(1,ip)*rod(2,ip)-ro(2,ip)*rod(1,ip)) &
                                      + dot_product(ori(3,1:3,ip)*momp(1:3,ip),angvelo(1:3,ip))
         ttras2(ipt)   = ttras2(ipt)  + Half*massp(ip)*dot_product(rod(1:3,ip),rod(1:3,ip))
         trots2(ipt)   = trots2(ipt)  + Half*dot_product(momp(1:3,ip),angvelo(1:3,ip)**2)
      end do

   case (iAfterMacrostep)

      norm = InvInt(nsamp2)
      fac = sclmas*sclvel**2/GasConstant/scltem
      nsamp1 = nsamp1 + 1
      for2s2 = for2s2*norm
      tor2s2 = tor2s2*norm
      lims2  = lims2 *norm
      anms2  = anms2 *norm
      ttras2 = ttras2/(Half*itradegfree*nstep2)*fac
      where (irotdegfree /= 0) trots2 = trots2/(Half*irotdegfree*nstep2)*fac
      for2s1 = for2s1 + for2s2
      tor2s1 = tor2s1 + tor2s2
      lims1  = lims1  + lims2
      anms1  = anms1  + anms2
      ttras1 = ttras1 +ttras2
      trots1 = trots1 +trots2

      call  MDAverWrite(for2s2, tor2s2, lims2, anms2, ttras2, trots2)

      write(ucnf) nsamp1, for2s1, tor2s1, lims1, anms1, ttras1, trots1

   case (iAfterSimulation)

      norm   = InvInt(nsamp1)
      for2s1 = for2s1*norm
      tor2s1 = tor2s1*norm
      lims1  = lims1*norm
      anms1  = anms1*norm
      ttras1 = ttras1*norm
      trots1 = trots1*norm

      call  MDAverWrite(for2s1, tor2s1, lims1, anms1, ttras1, trots1)

      deallocate(for2s2, for2s1, tor2s2, tor2s1, lims2, lims1, anms2, anms1, ttras2, ttras1, trots2, trots1)

   end select

contains

!........................................................................

subroutine MDAverWrite(for2, tor2, lim, anm, ttra, trot)
   real(8), intent(in) :: for2(npt), tor2(3,npt), lim(3,npt), anm(3,npt), ttra(npt), trot(npt)
   write(uout,'()')
   write(uout,'(a,t35,6e20.4)') '<force**2>                     = ', for2(1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<torque**2>x''                  = ', tor2(1,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<torque**2>y''                  = ', tor2(2,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<torque**2>z''                  = ', tor2(3,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<linear mom.>x                 = ', lim(1,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<linear mom.>y                 = ', lim(2,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<linear mom.>z                 = ', lim(3,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<angular mom.>x                = ', anm(1,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<angular mom.>y                = ', anm(2,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6e20.4)') '<angular mom.>z                = ', anm(3,1:npt)/nppt(1:npt)
   write(uout,'(a,t35,6f20.4)') '<t(tran)>                      = ', ttra
   write(uout,'(a,t35,6f20.4)') '<t(rot)>                       = ', trot
end subroutine MDAverWrite

!........................................................................

end subroutine MDAver

!************************************************************************
!*                                                                      *
!*     DispAver                                                         *
!*                                                                      *
!************************************************************************

! ... calculate averages of displacements

subroutine DispAver(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DispAver'
   integer(4), save :: nsamp1, nsamp2
   real(8), allocatable, save :: dros2(:,:), dros1(:,:) ! sum of displacements for particle ip
   real(8), allocatable, save :: rr2(:), rr2sd(:)       ! mean square displacement/step for particle type ipt
   real(8), allocatable, save :: rr3(:), rr3sd(:)       ! mean square displacement/step for particle type ipt
   integer(4) :: ip, ipt
   real(8)    :: norm, term

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      if(.not.allocated(dros2)) then
         allocate(dros2(3,np_alloc), dros1(3,np_alloc), rr2(npt), rr2sd(npt), rr3(npt), rr3sd(npt))
         dros2 = 0.0E+00
         dros1 = 0.0E+00
         rr2 = 0.0E+00
         rr2sd = 0.0E+00
         rr3 = 0.0E+00
         rr3sd = 0.0E+00
      end if

      nsamp1 = 0
      dros1  = Zero
      if (txstart == 'continue') read(ucnf) nsamp1, dros1

   case (iBeforeMacrostep)

      nsamp2 = 0
      dros2  = Zero

   case (iSimulationStep)

      nsamp2 = nsamp2 + 1
      dros2  = dros2 + drostep

   case (iAfterMacrostep)

      nsamp1 = nsamp1 + 1
      dros1 = dros1 + dros2
      rr2   = Zero
      rr2sd = Zero
      rr3   = Zero
      rr3sd = Zero
      do ip = 1, np
         ipt = iptpn(ip)
         term       = sum(dros2(1:3,ip)**2)
         rr2(ipt)   = rr2(ipt) + term
         rr2sd(ipt) = rr2sd(ipt) + term**2
         term       = sum(dros1(1:3,ip)**2)
         rr3(ipt)   = rr3(ipt) + term
         rr3sd(ipt) = rr3sd(ipt) + term**2
      end do
      rr2   = rr2/nppt(1:npt)
      rr3   = rr3/nppt(1:npt)
      do ipt = 1, npt
         rr2sd(ipt) = sqrt(max(Zero,(rr2sd(ipt)/nppt(ipt)-rr2(ipt)**2)*InvInt(nppt(ipt)-1)))
         rr3sd(ipt) = sqrt(max(Zero,(rr3sd(ipt)/nppt(ipt)-rr3(ipt)**2)*InvInt(nppt(ipt)-1)))
      end do

      write(uout,'()')
      norm = InvInt(nsamp2)
      write(uout,'(a,t35,6e20.4)') '<dr**2>/step (macrostep)       = ', rr2*norm
      write(uout,'(a,t35,6e20.4)') 'precision                      = ', rr2sd*norm
      norm = InvInt(istep1*nsamp2)
      write(uout,'(a,t35,6e20.4)') '<dr**2>/step (from start)      = ', rr3*norm
      write(uout,'(a,t35,6e20.4)') 'precision                      = ', rr3sd*norm

      write(ucnf) nsamp1, dros1

   case (iAfterSimulation)

      rr3   = Zero
      rr3sd = Zero
      do ip = 1, np
         ipt = iptpn(ip)
         term       = sum(dros1(1:3,ip)**2)
         rr3(ipt)   = rr3(ipt) + term
         rr3sd(ipt) = rr3sd(ipt) + term**2
      end do
      rr3  = rr3/nppt(1:npt)
      do ipt = 1, npt
         rr3sd(ipt) = sqrt(max(Zero,(rr3sd(ipt)/nppt(ipt)-rr3(ipt)**2)*InvInt(nppt(ipt)-1)))
      end do

      norm = InvInt(nsamp1*nsamp2)
      write(uout,'()')
      write(uout,'(a,t35,6e20.4)') '<dr**2>/step (from start)      = ', rr3*norm
      write(uout,'(a,t35,6e20.4)') 'precision                      = ', rr3sd*norm
      write(uout,'(a,t35,6f20.4)') 'sqrt(<dr**2>)(from start)      = ', sqrt(rr3)
      write(uout,'(a,t35,6f20.4)') 'precision                      = ', sqrt(rr3+rr3sd) - sqrt(rr3)

      deallocate(dros1, dros2, rr2, rr2sd, rr3, rr3sd)

!     call DispAverData
!     call DispPartData

   end select

contains

!........................................................................

subroutine DispAverData
   open(90, file = 'dispaver.data', position = 'append')
   write(90,'(5g12.3)') sqrt(rr3(1:npt)*norm)
   close(90)
end subroutine DispAverData

!........................................................................

subroutine DispPartData
   open(90, file = 'disppart.data', position = 'append')
   write(90,'(i5,g12.3)') (ip,sqrt(sum(dros1(1:3,ip)**2)),ip = 1, np)
   close(90)
end subroutine DispPartData

!........................................................................

end subroutine DispAver

!************************************************************************
!*                                                                      *
!*     OriOrderAver                                                     *
!*                                                                      *
!************************************************************************

! ... calculate averages of orientational order

subroutine OriOrderAver(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='OriOrderAver'
   integer(4), save :: nsamp1
   real(8), allocatable, save :: orist(:,:,:)               ! initial orientation of particle ip
   real(8), allocatable, save :: oris2(:,:), oris1(:,:)     ! orientational order parameter (0 to 1)
   integer(4) :: ip, ipt
   real(8)    :: norm

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      if(.not.allocated(orist)) then
         allocate(orist(3,3,np_alloc), oris2(3,npt), oris1(3,npt))
         orist = 0.0E+00
         oris2 = 0.0E+00
         oris1 = 0.0E+00
      end if

      nsamp1 = 0
      orist  = ori
      oris1  = Zero
      if (txstart == 'continue') read(ucnf) nsamp1, orist, oris1

   case (iBeforeMacrostep)

      oris2 = Zero

   case (iAfterMacrostep)

      nsamp1 = nsamp1 + 1
      do ip = 1, np
         ipt = iptpn(ip)
         oris2(1,ipt) = oris2(1,ipt)+Half*(Three*(dot_product(ori(1:3,1,ip),orist(1:3,1,ip)))**2-1)
         oris2(2,ipt) = oris2(2,ipt)+Half*(Three*(dot_product(ori(1:3,2,ip),orist(1:3,2,ip)))**2-1)
         oris2(3,ipt) = oris2(3,ipt)+Half*(Three*(dot_product(ori(1:3,3,ip),orist(1:3,3,ip)))**2-1)
      end do

      oris2(1,:) = oris2(1,:)/nppt(1:npt)
      oris2(2,:) = oris2(2,:)/nppt(1:npt)
      oris2(3,:) = oris2(3,:)/nppt(1:npt)
      oris1 = oris1 + oris2

      call OriOrderAverWrite(oris2)

      write(ucnf) nsamp1, orist, oris1

   case (iAfterSimulation)

      norm = InvInt(nsamp1)
      oris1 = oris1*norm

      call OriOrderAverWrite(oris1)

      deallocate(orist, oris2, oris1)

   end select

contains

!........................................................................

subroutine OriOrderAverWrite(oriori)
   real(8), intent(in) :: oriori(3,npt)
      write(uout,'()')
      write(uout,'(a,t35,6f20.4)') 'orientation order, x''-axis     = ', oriori(1,1:npt)
      write(uout,'(a,t35,6f20.4)') 'orientation order, y''-axis     = ', oriori(2,1:npt)
      write(uout,'(a,t35,6f20.4)') 'orientation order, z''-axis     = ', oriori(3,1:npt)
end subroutine OriOrderAverWrite

!........................................................................

end subroutine OriOrderAver

!************************************************************************
!*                                                                      *
!*     PosOriAver                                                       *
!*                                                                      *
!************************************************************************

! ... calculate averages of positions and orientations

subroutine PosOriAver(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='PosOriAver'
   integer(4), save :: nsamp1, nsamp2
   real(8), allocatable, save :: roavs2(:,:), roavs1(:,:)       ! average of ro for particles of type ipt
   real(8), allocatable, save :: oriavs2(:,:,:), oriavs1(:,:,:) ! average of ori for particles of type ipt
   integer(4) :: ip, ipt
   real(8)    :: norm

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iBeforeSimulation)

      if(.not.allocated(roavs2)) then
         allocate(roavs2(3,npt), roavs1(3,npt), oriavs2(3,3,npt), oriavs1(3,3,npt))
         roavs2 = 0.0E+00
         roavs1 = 0.0E+00
         oriavs2 = 0.0E+00
         oriavs1 = 0.0E+00
      end if

      nsamp1  = 0
      roavs1  = Zero
      oriavs1 = Zero
      if (txstart == 'continue') read(ucnf) nsamp1, roavs1, oriavs1

   case (iBeforeMacrostep)

      nsamp2  = 0
      roavs2  = Zero
      oriavs2 = Zero

   case (iSimulationStep)

      nsamp2 = nsamp2 + 1
      do ip = 1, np
         ipt = iptpn(ip)
         roavs2(1:3,ipt)      = roavs2(1:3,ipt)     +ro(1:3,ip)
         oriavs2(1:3,1:3,ipt) = oriavs2(1:3,1:3,ipt)+ori(1:3,1:3,ip)
      end do

   case (iAfterMacrostep)

      norm = InvInt(nsamp2)
      nsamp1 = nsamp1 + 1
      do ipt = 1, npt
         roavs2(1:3,ipt)      = roavs2(1:3,ipt)     *norm/nppt(ipt)
         oriavs2(1:3,1:3,ipt) = oriavs2(1:3,1:3,ipt)*norm/nppt(ipt)
      end do
      roavs1  = roavs1  + roavs2
      oriavs1 = oriavs1 + oriavs2

      call PosOriAverWrite(roavs2, oriavs2)

      write(ucnf) nsamp1, roavs1, oriavs1

   case (iAfterSimulation)

      norm = InvInt(nsamp1)
      roavs1  = roavs1*norm
      oriavs1 = oriavs1*norm

      call PosOriAverWrite(roavs1, oriavs1)

      deallocate(roavs2, roavs1, oriavs2, oriavs1)

   end select

contains

!........................................................................

subroutine PosOriAverWrite(rotav, oriav)
   real(8), intent(in) :: rotav(3,npt), oriav(3,3,npt)
      integer(4) :: isign, imagn
      character  :: char*1, fmt*20
      call SignMagn(maxval(abs(rotav)), isign, imagn)
      write(char,'(i1)') 3 - min(3,max(0,imagn))
      fmt = '(a,t37,6(3f6.'//char//',2x))'
      write(uout,'()')
      write(uout,fmt)                   '<ro> (x,y,z)                   = ', rotav(1:3,1:npt)
      write(uout,'(a,t37,6(3f6.2,2x))') '<x''> (x,y,z)                   = ', oriav(1:3,1,1:npt)
      write(uout,'(a,t37,6(3f6.2,2x))') '<y''> (x,y,z)                   = ', oriav(1:3,2,1:npt)
      write(uout,'(a,t37,6(3f6.2,2x))') '<z''> (x,y,z)                   = ', oriav(1:3,3,1:npt)
end subroutine PosOriAverWrite

!........................................................................

end subroutine PosOriAver

!************************************************************************
!*                                                                      *
!*     MainAver                                                         *
!*                                                                      *
!************************************************************************

! ... calculate various quantities

subroutine MainAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MainAver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)
                      call ThermoAver(iStage)
   if (lweakcharge)   call ChargeAver(iStage)
   if (lpolarization) call IndDipMomAver(iStage)
   if (lchain)        call ChainAver(iStage)
   if (lnetwork)      call NetworkAver(iStage)
   if (lhierarchical) call HierarchicalAver(iStage)
   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

end subroutine MainAver

!************************************************************************
!*                                                                      *
!*     ThermoAver                                                       *
!*                                                                      *
!************************************************************************

! ... calculate averages of thermodynamic quantities

!          no            quantity
!          --            ---------------------
!          1             total energy
!          2             kinetic energy
!          3             total potential energy
!          4             total two-body potential energy  u%twob(iptjpt = 0)
!          5             u%twob(iptjpt = 1)
!          6             u%twob(iptjpt = 2)
!          ...           ...
!          nptpt+5       potential energy from reciprocal space  u%rec
!          nptpt+6       electrostatic energy  u%stat
!          nptpt+7       polarization energy  u%pol
!          nptpt+8       total one-particle energy u%oneb(ipt = 0)
!          nptpt+9       u%oneb(ipt = 1)
!          nptpt+10      u%oneb(ipt = 2)
!          ...           ...
!          nptpt+npt+9   bond energy  u%bond
!          nptpt+npt+10  angle energy  u%angle
!          nptpt+npt+11  crosslink energy u%crosslink
!          nptpt+npt+12  external potential u%external
!          nptpt+npt+13  enthalpy
!          nptpt+npt+14  temperature
!          nptpt+npt+15  pressure
!          nptpt+npt+16  volume
!          nptpt+npt+17  number of particles

subroutine ThermoAver(iStage)

   use MolModule
#if !defined (_NOIEEE_)
   use, intrinsic :: IEEE_ARITHMETIC
#endif
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ThermoAver'
   character(40), parameter :: txroutine_aux ='MainAver'
   character(80), parameter :: txheading ='thermodynamic quantities'
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4),       save :: ietot, iekin, iutot, iutwob, iurec, iustat, iupol, iuoneb, &
                             iubond, iuangle, iuclink, iuext, ihtot, itemp, iprsr, ivol, inpart
   character(28),    save :: fmt1  = '(a,t35,f15.5,f15.5,e15.3)'
   character(29),    save :: fmt1b = '(a,t35,es15.5,es15.5,e15.3)'
   character(28),    save :: fmt1c = '(a,t35,f15.5)'
   character(60),    save :: fmt2  = '(a,t35,2(f15.5,f15.5),f15.0)'
   character(60),    save :: fmt2b = '(a,t35,2(es15.5,es15.5),f15.0)'
   type(potenergy_var) :: ucheck
   integer(4) :: ipt, iptjpt, ivar
   real(8)    :: rtemp, rtemp2, cv, cp, InvFlt, GetRelDiff, fac

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      ietot   =           1
      iekin   = ietot   + 1
      iutot   = iekin   + 1
      iutwob  = iutot   + 1
      iurec   = iutwob  + 1 + nptpt
      iustat  = iurec   + 1
      iupol   = iustat  + 1
      iuoneb  = iupol   + 1
      iubond  = iuoneb  + 1 + npt
      iuangle = iubond  + 1
      iuclink = iuangle + 1
      iuext   = iuclink + 1
      ihtot   = iuext   + 1
      itemp   = ihtot   + 1
      iprsr   = itemp   + 1
      ivol    = iprsr   + 1
      inpart  = ivol    + 1
      nvar    = inpart
      allocate(var(nvar))

      var(ietot                  )%label = 'total energy                   = '
      var(iekin                  )%label = 'kinetic energy                 = '
      var(iutot                  )%label = 'total potential energy         = '
      var(iutwob                 )%label = ' total two-body energy         = '
      var(iutwob+1:iutwob+nptpt)%label = '  '//txptpt(1:nptpt)//'        = '
      var(iurec                  )%label = ' pot. energy from rec. space   = '
      var(iustat                 )%label = ' electrostatic energy          = '
      var(iupol                  )%label = ' polarization energy           = '
      var(iuoneb                 )%label = ' total one-body energy         = '
      var(iuoneb+1:iuoneb+npt    )%label = '  '//txpt(1:npt)// '                   = '
      var(iubond                 )%label = ' bond energy                   = '
      var(iuangle                )%label = ' angle energy                  = '
      var(iuclink                )%label = ' crosslink energy              = '
      var(iuext                  )%label = ' external potential            = '
      var(ihtot                  )%label = 'enthalpy                       = '
      var(itemp                  )%label = 'temperature                    = '
      var(iprsr                  )%label = 'pressure (excl. cont. contr.)  = '
      var(ivol                   )%label = 'volume                         = '
      var(inpart                 )%label = 'number of particles            = '

      var(1:itemp-1)%norm  = One/np      ! these quantities will be given per particle
      var(itemp:nvar)%norm = One

   case (iWriteInput)

! ... initiate energy, forces, and pressure

      if (ltime) call CpuAdd('stop', txroutine_aux, 0, uout)
      if (ltime) call CpuAdd('interrupt', ' ', 0, uout)

      call SetAtomProp(1, np, lintsite)
      call UTotal(iStage)
      if (lintsite) call SetAtomPos(1, np, .false.)

      if (ltime) call CpuAdd('resume', ' ', 0, uout)
      if (ltime) call CpuAdd('start', txroutine_aux, 0, uout)

      if (lmd)                 call GetLinAcc(1, np)
      if (lmd .and. lpolyatom) call GetAngAcc(1, np)
      if (lmd)                 call GetKinEnergy
      prsr = (np*Boltz*temp*scltem - virial*sclene/(Three*AvNo))/(vol*sclvol)/sclpre
      htot = (u%tot+prsr*sclpre*vol*sclvol*AvNo/sclene)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t37,a)') 'quantity', 'initial value'
         write(uout,'(a,t37,a)') '&-------', '-------------'
         if (lmd)                     write(uout,fmt1) var(ietot         )%label, (ekin+u%tot)/np
         if (lmd)                     write(uout,fmt1) var(iekin         )%label, ekin/np
                                      write(uout,fmt1) var(iutot         )%label, u%tot/np
                                      write(uout,fmt1c)(var(iutwob+iptjpt)%label, u%twob(iptjpt)/np,iptjpt = 0,nptpt)
         if (lcharge .and. lewald)    write(uout,fmt1) var(iurec         )%label, u%rec/np
         if (lweakcharge .and. lewald) write(uout,fmt1) var(iurec         )%label, u%rec/np
         if (ldipole .or. ldipolesph) write(uout,fmt1) var(iustat        )%label, u%stat/np
         if (lpolarization)           write(uout,fmt1) var(iustat        )%label, u%stat/np
         if (lpolarization)           write(uout,fmt1) var(iupol         )%label, u%pol/np
         if (ldieldis)                write(uout,fmt1c)(var(iuoneb+ipt   )%label, u%oneb(ipt)/np,ipt = 0,npt)
         if (lchain)                  write(uout,fmt1) var(iubond        )%label, u%bond/np
         if (lchain)                  write(uout,fmt1) var(iuangle       )%label, u%angle/np
         if (lclink)                  write(uout,fmt1) var(iuclink       )%label, u%crosslink/np
         if (luext)                   write(uout,fmt1) var(iuext         )%label, u%external/np
                                      write(uout,fmt1) var(ihtot         )%label, htot/np
                                      write(uout,fmt1) var(itemp         )%label, temp
                                      write(uout,fmt1) var(iprsr         )%label, prsr
                                      write(uout,fmt1) var(inpart        )%label, real(np)
                                      write(uout,fmt1b)var(ivol          )%label, vol
                                      write(uout,'()')

         if (temp > Zero) then
            write(uout,fmt1) 'total potential energy (NkT)   = ', u%tot*sclene*InvFlt(np*GasConstant*temp*scltem)
            write(uout,fmt1) 'pressure (NkT/V)(excl. contact)= ', prsr*sclpre*vol*sclvol*InvFlt(np*Boltz*temp*scltem)
            write(uout,fmt1b) 'pressure (kT) (excl. contact)  = ', prsr*sclpre*sclvol*InvFlt(Boltz*temp*scltem)
         end if

      end if

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

      if (iaver/= 0 .and. master) then
         write(uout,*)
         write(uout,'(2x,a,t15,a,t30,a,t45,a,t60,a)') 'step/pass', 'pot. energy', 'temperature', 'pressure', 'volume'
         call FileFlush(uout)
      end if

   case (iSimulationStep)

      htot = u%tot + prsr*sclpre*vol*sclvol*AvNo/sclene

      var(ietot                )%value = ekin+u%tot
      var(iekin                )%value = ekin
      var(iutot                )%value = u%tot
      var(iutwob:iutwob+nptpt)%value = u%twob(0:nptpt)
      var(iurec                )%value = u%rec
      var(iustat               )%value = u%stat
      var(iupol                )%value = u%pol
      var(iuoneb:iuoneb+npt    )%value = u%oneb(0:npt)
      var(iubond               )%value = u%bond
      var(iuangle              )%value = u%angle
      var(iuclink              )%value = u%crosslink
      var(iuext                )%value = u%external
      var(ihtot                )%value = htot
      var(itemp                )%value = temp
      var(iprsr                )%value = prsr
      var(ivol                 )%value = vol
      var(inpart               )%value = np

      call ScalarSample(iStage, 1, nvar, var)

      if (iaver/= 0 .and. master) then
         if (mod(istep2,iaver) == 0) then
            write(uout,'(i6,3f15.3,f15.3)')  istep2, &
            var(iutot)%avs2/(real(np)*istep2), var(itemp)%avs2/istep2, var(iprsr)%avs2/istep2, var(ivol)%avs2/istep2
            call FileFlush(uout)
         end if
      end if

     continue

   case (iAfterMacrostep)

      if(.not.allocated(ucheck%twob)) then
         allocate(ucheck%twob(0:nptpt))
      end if
      if (lmc) then
         if (ltime) call CpuAdd('stop', txroutine_aux, 0, uout)
         if (ltime) call CpuAdd('interrupt', ' ', 0, uout)
         ucheck = u                                          ! save current u
         call UTotal(iStage)                                 ! calcuate potential energies from scratch
         ucheck%crosslink     = GetRelDiff(ucheck%crosslink,u%crosslink)
!         call TestThermoAver(uout)
         ucheck%tot           = GetRelDiff(ucheck%tot,u%tot) ! calculate relative differences
         do iptjpt = 0, nptpt
            ucheck%twob(iptjpt) = GetRelDiff(ucheck%twob(iptjpt),u%twob(iptjpt))
         end do
         ucheck%rec           = GetRelDiff(ucheck%rec,u%rec)
         ucheck%stat          = GetRelDiff(ucheck%stat,u%stat)
         do ipt = 0, npt
            ucheck%oneb(ipt)  = GetRelDiff(ucheck%oneb(ipt),u%oneb(ipt))
         end do
         ucheck%bond          = GetRelDiff(ucheck%bond,u%bond)
         ucheck%angle         = GetRelDiff(ucheck%angle,u%angle)
         ucheck%external      = GetRelDiff(ucheck%external,u%external)
         if (ltime) call CpuAdd('resume', ' ', 0, uout)
         if (ltime) call CpuAdd('start', txroutine_aux, 0, uout)
      end if

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var

      rtemp = GasConstant*(var(itemp)%avs2*scltem)
      rtemp2 = GasConstant*(var(itemp)%avs2*scltem)**2
      if (lnve) then
         if(rtemp == 0.0d0) then !prevent division by zero
#if !defined (_NOIEEE_)
            cv = IEEE_VALUE(cv,IEEE_QUIET_NAN)
#else
            cv = huge(cv)
#endif
         else
            cv = 1.5*GasConstant/(one-(var(iekin)%fls2*sclene)**2/(1.5*np*rtemp**2))/sclhca
         end if
      end if
      if (lnvt) cv = (var(iutot)%fls2*sclene)**2/(rtemp2*sclhca*np)
      if (lntp) cp = (var(ihtot)%fls2*sclene)**2/(rtemp2*sclhca*np)

      call ScalarNorm(iStage, 1, nvar, var, 0)

      tempaver = var(itemp)%avs2
      prsraver = var(iprsr)%avs2
      volaver = var(ivol)%avs2
      npartaver = var(inpart)%avs2

      if (master) then
         call WriteHead(2, txheading, uout)
         if (lmc) then
            write(uout,900) 'quantity', 'average', 'fluctuation', '      check'
            write(uout,900) '&&------', '-------', '-----------', '      -----'
         else
            write(uout,900) 'quantity', 'average', 'fluctuation'
            write(uout,900) '&&------', '-------', '-----------'
         end if
  900    format(a,t43,a,t54,a,t69,a)
         if (lmd) write(uout,fmt1) var(ietot )%label, var(ietot )%avs2, var(ietot )%fls2
         if (lmd) write(uout,fmt1) var(iekin )%label, var(iekin )%avs2, var(iekin )%fls2
         call ThermoAverSub(lmc, fmt1, var(iutot)%label, var(iutot)%avs2, var(iutot)%fls2, ucheck%tot)
         do ivar = iutwob, iutwob+nptpt
         call ThermoAverSub(lmc, fmt1, var(ivar )%label, var(ivar )%avs2, var(ivar )%fls2, ucheck%twob(ivar-iutwob))
         end do
         if (lcharge .and. lewald) call ThermoAverSub(lmc, fmt1, var(iurec)%label, var(iurec)%avs2, var(iurec)%fls2, ucheck%rec)
         if (lweakcharge .and. lewald) call ThermoAverSub(lmc, fmt1, var(iurec)%label, var(iurec)%avs2, var(iurec)%fls2, ucheck%rec)
         if (ldipole .or. ldipolesph) call ThermoAverSub(lmc, fmt1, var(iustat)%label, var(iustat)%avs2, var(iustat)%fls2, ucheck%stat)
         if (lpolarization) call ThermoAverSub(lmc, fmt1, var(iustat)%label, var(iustat)%avs2, var(iustat)%fls2, ucheck%stat)
         if (lpolarization) call ThermoAverSub(lmc, fmt1, var(iupol )%label, var(iupol )%avs2, var(iupol )%fls2, ucheck%pol )
         if (ldieldis) then
            do ivar = iuoneb, iuoneb+npt
            call ThermoAverSub(lmc, fmt1, var(ivar )%label, var(ivar )%avs2, var(ivar )%fls2, ucheck%oneb(ivar-iuoneb))
            end do
         end if
         if (lchain) call ThermoAverSub(lmc, fmt1, var(iubond )%label, var(iubond )%avs2, var(iubond )%fls2, ucheck%bond)
         if (lchain) call ThermoAverSub(lmc, fmt1, var(iuangle)%label, var(iuangle)%avs2, var(iuangle)%fls2, ucheck%angle)
         if (lclink) call ThermoAverSub(lmc, fmt1, var(iuclink)%label, var(iuclink)%avs2, var(iuclink)%fls2, ucheck%crosslink)
         if (luext ) call ThermoAverSub(lmc, fmt1, var(iuext  )%label, var(iuext  )%avs2, var(iuext  )%fls2, ucheck%external)
         write(uout,fmt1) var(ihtot)%label, var(ihtot)%avs2, var(ihtot)%fls2
         if (lnve .or. lnvt) write(uout,fmt1) 'cv                             = ', cv
         if (lntp)           write(uout,fmt1) 'cp                             = ', cp
         write(uout,fmt1) var(itemp )%label, var(itemp )%avs2, var(itemp )%fls2
         write(uout,fmt1) var(iprsr )%label, var(iprsr )%avs2, var(iprsr )%fls2
         write(uout,fmt1) var(inpart)%label, var(inpart)%avs2, var(inpart)%fls2
         write(uout,fmt1b)var(ivol  )%label, var(ivol  )%avs2, var(ivol  )%fls2
         write(uout,'()')
         write(uout,fmt1) 'total potential energy (NkT)   = ', var(iutot)%avs2*sclene*InvFlt(rtemp)
         write(uout,fmt1) 'pressure (NkT/V)(excl. contact)= ', var(iprsr)%avs2*sclpre*var(ivol)%avs2*sclvol*InvFlt(np*rtemp)*AvNo
         write(uout,fmt1b) 'pressure (kT) (excl. contact)  = ', var(iprsr)%avs2*sclpre*sclvol*InvFlt(rtemp)*AvNo
         write(uout,fmt1b) 'number density                 = ', var(inpart)%avs2/var(ivol)%avs2
         if (txuser == 'sim_annealing') then  ! output from the last macrostep
            call FileOpen(uuser, fuser, 'form/noread')
            if (istep1 == nstep1) write(uuser,'(f10.1,f12.4)')  qsuperball, np*var(iutot)%avs2/(AutokJ*Bohr)
         end if

      end if

   case (iAfterSimulation)

   if (ltrace) call WriteTrace(2, trim(txroutine)// '1', iStage)
      if (master) then

         call ScalarSample(iStage, 1, nvar, var)
         rtemp = GasConstant*(var(itemp)%avs1*scltem)
         rtemp2 = GasConstant*(var(itemp)%avs1*scltem)**2
         if (lnve) then
            if(rtemp == 0.0d0) then !prevent division by zero
#if !defined (_NOIEEE_)
               cv = IEEE_VALUE(cv,IEEE_QUIET_NAN)
#else
               cv = huge(cv)
#endif
            else
               cv = 1.5*GasConstant/(one-(var(iekin)%fls1*sclene)**2/(1.5*np*rtemp**2))/sclhca
            end if
         end if
         if (lnvt) cv = (var(iutot)%fls1*sclene)**2/(rtemp2*sclhca*np)
         if (lntp) cp = (var(ihtot)%fls1*sclene)**2/(rtemp2*sclhca*np)

         call ScalarNorm(iStage, 1, nvar, var, 0)
         tempaver = var(itemp)%avs1
         prsraver = var(iprsr)%avs1
         volaver = var(ivol )%avs1
         npartaver = var(inpart)%avs1

         call WriteHead(2, txheading, uout)
         write(uout,910) 'quantity', 'average', ' precision', 'fluctuation', 'precision','stat eff'
         write(uout,910) '&&&-----', '-------', ' ---------', '-----------', '---------','--------'
  910    format(a,t43,a,t55,a,t69,a,t86,a,t102,a)
         if (lmd) call TempWrite(ietot  ,fmt2,uout)
         if (lmd) call TempWrite(iekin  ,fmt2,uout)
                  call TempWrite(iutot  ,fmt2,uout)
         do ivar = iutwob, iutwob+nptpt
                  call TempWrite(ivar   ,fmt2,uout)
         end do
         if (lcharge .and. lewald) call TempWrite(iurec ,fmt2,uout)
         if (lweakcharge .and. lewald) call TempWrite(iurec ,fmt2,uout)
         if (ldipole .or. ldipolesph) call TempWrite(iustat ,fmt2,uout)
         if (lpolarization) call TempWrite(iustat ,fmt2,uout)
         if (lpolarization) call TempWrite(iupol  ,fmt2,uout)
         if (ldieldis) then
            do ivar = iuoneb, iuoneb+npt
               call TempWrite(ivar   ,fmt2,uout)
            end do
         end if
         if (lchain) call TempWrite(iubond ,fmt2,uout)
         if (lchain) call TempWrite(iuangle,fmt2,uout)
         if (lclink) call TempWrite(iuclink,fmt2,uout)
         if (luext ) call TempWrite(iuext  ,fmt2,uout)
                     call TempWrite(ihtot  ,fmt2,uout)
         if (lnve .or. lnvt) write(uout,fmt1) 'cv                             = ', cv
         if (lntp)           write(uout,fmt1) 'cp                             = ', cp
                     call TempWrite(itemp  ,fmt2,uout)
                     call TempWrite(iprsr  ,fmt2,uout)
                     call TempWrite(inpart ,fmt2,uout)
                     call TempWrite(ivol   ,fmt2b,uout)
         write(uout,'()')
         fac = sclene*InvFlt(rtemp)
         write(uout,fmt2) 'total potential energy (NkT)   = ', var(iutot)%avs1*fac, var(iutot)%avsd*fac
         fac = sclpre*var(ivol)%avs1*sclvol*InvFlt(np*rtemp)*AvNo
         write(uout,fmt2) 'pressure (NkT/V)(excl. contact)= ', var(iprsr)%avs1*fac, var(iprsr)%avsd*fac
         fac = sclpre*sclvol*InvFlt(rtemp)*AvNo
         write(uout,fmt2b) 'pressure (kT) (excl. contact)  = ', var(iprsr)%avs1*fac, var(iprsr)%avsd*fac
         write(uout,fmt2b) 'number density                 = ', var(inpart)%avs1/var(ivol)%avs1, &
           var(inpart)%avs1/var(ivol)%avs1*(var(inpart)%avsd/var(inpart)%avs1+var(ivol)%avsd/var(ivol)%avs1)

         fac = sclpre*var(ivol)%avs1*sclvol*InvFlt(np*rtemp)*AvNo
         prsrreds3 = var(iprsr)%avs1*fac
         prsrredsd = var(iprsr)%avsd*fac

     !   if (allocated(var)) deallocate(var) ! don't allways work with lmpt

      end if

   end select

contains

!........................................................................

subroutine ThermoAverSub(lmc, fmt, txv, data1, data2, data3)
   logical,      intent(in) :: lmc
   character(*), intent(in) :: txv
   character(*), intent(in) :: fmt
   real(8),      intent(in) :: data1, data2, data3
   if (lmc) then
      write(uout,fmt) txv, data1, data2, data3
   else
      write(uout,fmt) txv, data1, data2
   end if
end subroutine ThermoAverSub

!........................................................................

subroutine TempWrite(i,fmt,unit)
   integer(4),   intent(in) :: i
   character(*), intent(in) :: fmt
   integer(4),   intent(in) :: unit
   integer(4) :: nblocklen              ! number of block lengths
   nblocklen = var(i)%nblocklen
!  write(unit,fmt) trim(var(i)%label), var(i)%avs1, var(i)%avsd, var(i)%fls1, var(i)%flsd  ! fixed blocklen
   write(unit,fmt) trim(var(i)%label), var(i)%av_s1(1), var(i)%av_sd_extrap, var(i)%fl_s1(nblocklen), &
                                       var(i)%fl_sd_extrap, var(i)%av_sd_stateff           ! variable blocklen
end subroutine TempWrite

!........................................................................

subroutine TestThermoAver(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,i5)') 'istep1 ', istep1
   write(unit,'(a,g15.5)')  'ucheck%tot(old) ', ucheck%tot
   write(unit,'(a,g15.5)')  'u%tot(upd)      ', u%tot
   write(unit,'(a,9g15.5)') 'ucheck%twob(old)', ucheck%twob
   write(unit,'(a,9g15.5)') 'u%twob(upd)     ', u%twob
end subroutine TestThermoAver

!........................................................................

end subroutine ThermoAver

!************************************************************************
!*                                                                      *
!*     ChargeAver                                                       *
!*                                                                      *
!************************************************************************

! ... calculate averages of related to weak charges

!          no            quantity
!          --            ---------------------
!          1             total charge
!          2             fractional charge: atom (iat = 1)
!          3             fractional charge: atom (iat = 2)
!          ...           ...

subroutine ChargeAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChargeAver'
   character(80), parameter :: txheading ='charge averages'
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4) :: ia, iat, ivar
   real(8) :: InvFlt

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1+nat
      allocate(var(nvar))

      var(1              )%label = 'total charge                   = '
      var(2:1+nat        )%label = 'fractional charge: '//txat(1:nat)//'  = '
      var(1              )%norm = One
      var(2:1+nat        )%norm = One/(nppt(iptat(1:nat))*naat(1:nat))
      do iat = 1, nat
         var(1+iat)%norm = var(1+iat)%norm*InvFlt(abs(zat(iat)))
      end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var%value = Zero
      do ia = 1, na
         ivar = 1
         if (laz(ia)) var(ivar)%value = var(ivar)%value + az(ia)
         iat = iatan(ia)
         ivar = 1+iat
         if (laz(ia)) var(ivar)%value = var(ivar)%value + abs(az(ia))
      end do
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.3,f15.0)', uout)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.3,f15.0)', uout)

      deallocate(var)

   end select

end subroutine ChargeAver

!************************************************************************
!*                                                                      *
!*     IndDipMomAver                                                    *
!*                                                                      *
!************************************************************************

! ... calculate averages of induced dipole moment

!          no            quantity
!          --            ---------------------
!          1             induced dipole moment: total
!          2             induced dipole moment: particle (ipt = 1)
!          3             induced dipole moment: particle (ipt = 2)
!          ...           ...
!          npt+2         induced dipole moment: atom (kat = 1)
!          npt+3         induced dipole moment: atom (kat = 2)
!          ...           ...

subroutine IndDipMomAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IndDipMovAver'
   character(80), parameter :: txheading ='induced dipole moments'
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4) :: ip, ipt, ia, iat, ivar

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = npt+nat+1
      allocate(var(nvar))

      var(1              )%label = 'ind. dip. mom.: total          = '
      var(2:1+npt        )%label = 'ind. dip. mom.: '//txpt(1:npt)//'     = '
      var(2+npt:1+npt+nat)%label = 'ind. dip. mom.: '//txat(1:nat)//'     = '
      var(1              )%norm = One
      var(2:1+npt        )%norm = One/nppt(1:npt)
      var(2+npt:1+npt+nat)%norm = One/(naat(1:nat)*nppt(iptat(1:nat)))

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var%value = Zero
      ivar = 1
      var(ivar)%value = sqrt(idmsys(1)**2+idmsys(2)**2+idmsys(3)**2)
      do ip = 1, np
         ipt = iptpn(ip)
         ivar = 1+ipt
         var(ivar)%value = var(ivar)%value + sqrt(idmo(1,ip)**2+idmo(2,ip)**2+idmo(3,ip)**2)
      end do
      do ia = 1, na
         iat = iatan(ia)
         ivar = 1+npt+iat
         var(ivar)%value = var(ivar)%value + sqrt(idm(1,ia)**2+idm(2,ia)**2+idm(3,ia)**2)
      end do
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

end subroutine IndDipMomAver

!************************************************************************
!*                                                                      *
!*     ChainAver                                                        *
!*                                                                      *
!************************************************************************

! ... calculate averages of chain quantities

subroutine ChainAver(iStage)

   use MolModule
#if !defined (_NOIEEE_)
   use, intrinsic :: IEEE_ARITHMETIC
#endif
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainAver'
   character(80), parameter :: txheading ='chain quantities'
   integer(4)   , parameter :: ntype = 13          ! number of types of properties
   integer(4)      , save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   type(chainprop_var) :: ChainProperty
   integer(4) :: ic, ict, ioffset
   real(8) :: PerLengthRg, Asphericity, InvFlt

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = nct*ntype
      allocate(var(nvar))

      do ict = 1, nct
         ioffset = ntype*(ict-1)
         var( 1+ioffset)%label = '<r(bb)**2>**0.5                = ' ! rms bead-to-bead separation
         var( 2+ioffset)%label = '<r(ee)**2>**0.5                = ' ! rms end-to-end separation
         var( 3+ioffset)%label = '<r(g)**2>**0.5                 = ' ! rms radius of gyration
         var( 4+ioffset)%label = 'smallest rms mom. p.a.         = ' ! smallest rms moment along a prinical axis
         var( 5+ioffset)%label = 'intermediate rms mom. p.a.     = ' ! intermediate rms moment along a prinical axis
         var( 6+ioffset)%label = 'largest rms mom. p.a.          = ' ! largest rms moment along a prinical axis
         var( 7+ioffset)%label = '<persist. l.>  (r(ee))         = ' ! persistence length (ref JCP 107, 1279 (1997))
         var( 8+ioffset)%label = '<persist. l.>  (r(g))          = ' ! persistence length (ref JCP 107, 1279 (1997))
         var( 9+ioffset)%label = '<angle>                        = ' ! angle between consecutive beads
         var(10+ioffset)%label = '<cos(180-angle)>               = ' ! cos(180-angle)
         var(11+ioffset)%label = '<r(ee)**2/r(g)**2>             = ' ! square end-to-end separation / square radius of gyration
         var(12+ioffset)%label = '<asphericity>                  = ' ! asphericity
         var(13+ioffset)%label = '<toroid param>                 = ' ! toroid parameter
         var(1+ioffset:13+ioffset)%norm = One/ncct(ict)
     end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      if (lbd .and. txuser == 'chainads' .and. mod(istep2,10) /= 0) goto 10

      var%value = Zero
      do ic = 1, nc
         ict = ictcn(ic)
         ioffset = ntype*(ict-1)
         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)
         var( 1+ioffset)%value = var( 1+ioffset)%value + ChainProperty%rbb2
         var( 2+ioffset)%value = var( 2+ioffset)%value + ChainProperty%ree2
         var( 3+ioffset)%value = var( 3+ioffset)%value + ChainProperty%rg2
         var( 4+ioffset)%value = var( 4+ioffset)%value + ChainProperty%rg2s
         var( 5+ioffset)%value = var( 5+ioffset)%value + ChainProperty%rg2m
         var( 6+ioffset)%value = var( 6+ioffset)%value + ChainProperty%rg2l
         var( 7+ioffset)%value = var( 7+ioffset)%value + ChainProperty%lpree
         var( 8+ioffset)%value = var( 8+ioffset)%value + ChainProperty%lprg
         var( 9+ioffset)%value = var( 9+ioffset)%value + ChainProperty%angle
         var(10+ioffset)%value = var(10+ioffset)%value + ChainProperty%cos
         var(11+ioffset)%value = var(11+ioffset)%value + ChainProperty%shape
         var(12+ioffset)%value = var(12+ioffset)%value + ChainProperty%asph
         var(13+ioffset)%value = var(13+ioffset)%value + ChainProperty%torp
      end do
!      write(*,*) 'ree2',chainproperty%ree2
      call ScalarSample(iStage, 1, nvar, var)

  10 continue

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call WriteHead(2, txheading, uout)
      do ict = 1, nct
         ioffset = ntype*(ict-1)
         call ScalarNorm(iStage, 1+ioffset, 6+ioffset, var, 2)
         call ScalarNorm(iStage, 7+ioffset, ntype+ioffset, var, 0)
         if (nct > 1) write(uout,'(2a)') 'chain: ', txct(ict)
         if (nct > 1) write(uout,'(a)')  '------------------'
         call ScalarWrite(iStage, 1+ioffset, ntype+ioffset, var, 1, '(a,t35,4f15.5,f15.0)', uout)
         write(uout,'()')
         write(uout,'(a,t35,2f15.5)') 'persist. l. (<r(ee)**2>)       = ', &
           var(2+ioffset)%avs2**2*InvFlt(Two*(npct(ict)-1)*(var(1+ioffset)%avs2))+Half*(var(1+ioffset)%avs2)
         write(uout,'(a,t35,2f15.5)') 'persist. l. (<r(g)**2>)        = ', &
           PerLengthRg(var(3+ioffset)%avs2**2, (npct(ict)-1)*(var(1+ioffset)%avs2))
         if (var(10+ioffset)%avs2 == One) then
            write(uout,'(a,t35,2f15.5)') 'persist. l. (<cos(180-angle)>) = ', &
#if !defined (_NOIEEE_)
            IEEE_VALUE((var(10+ioffset)%avs2),IEEE_QUIET_NAN)
#else
            huge(var(10+ioffset)%avs2)
#endif
         else
            write(uout,'(a,t35,2f15.5)') 'persist. l. (<cos(180-angle)>) = ', &
            (var(1+ioffset)%avs2)/(One-var(10+ioffset)%avs2)
         end if
         write(uout,'(a,t35,2f15.5)') '<r(ee)**2>/<r(g)**2>           = ', &
           var(2+ioffset)%avs2**2*InvFlt(var(3+ioffset)%avs2**2)
         write(uout,'(a,t35,2f15.5)') 'asphericity (<2:nd moments>)   = ', &
           Asphericity(var(4+ioffset)%avs2**2,var(5+ioffset)%avs2**2,var(6+ioffset)%avs2**2)
         write(uout,'()')
      end do

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call WriteHead(2, txheading, uout)
      do ict = 1, nct
         ioffset = ntype*(ict-1)
         call ScalarNorm(iStage, 1+ioffset, 6+ioffset, var, 2)
         call ScalarNorm(iStage, 7+ioffset, ntype+ioffset, var, 0)
         if (nct > 1) write(uout,'(2a)') 'chain: ', txct(ict)
         if (nct > 1) write(uout,'(a)')  '------------------'
         call ScalarWrite(iStage, 1+ioffset, ntype+ioffset, var, 1, '(a,t35,4f15.5,f15.0)', uout)
         write(uout,'()')
         write(uout,'(a,t35,2f15.5)') 'persist. l. (<r(ee)**2>)       = ', &
           var(2+ioffset)%avs1**2*InvFlt(Two*(npct(ict)-1)*(var(1+ioffset)%avs1))+Half*(var(1+ioffset)%avs1)
         write(uout,'(a,t35,2f15.5)') 'persist. l. (<r(g)**2>)        = ', &
           PerLengthRg(var(3+ioffset)%avs1**2, (npct(ict)-1)*(var(1+ioffset)%avs1))
         if (var(10+ioffset)%avs1 == One) then
            write(uout,'(a,t35,2f15.5)') 'persist. l. (<cos(180-angle)>) = ', &
#if !defined (_NOIEEE_)
            IEEE_VALUE((var(10+ioffset)%avs1),IEEE_QUIET_NAN)
#else
            huge(var(10+ioffset)%avs1)
#endif
         else
            write(uout,'(a,t35,2f15.5)') 'persist. l. (<cos(180-angle)>) = ', &
            (var(1+ioffset)%avs1)/(One-var(10+ioffset)%avs1)
         end if
         write(uout,'(a,t35,2f15.5)') '<r(ee)**2>/<r(g)**2>           = ', &
           var(2+ioffset)%avs1**2*InvFlt(var(3+ioffset)%avs1**2)
         write(uout,'(a,t35,2f15.5)') 'asphericity (<2:nd moments>)   = ', &
           Asphericity(var(4+ioffset)%avs1**2,var(5+ioffset)%avs1**2,var(6+ioffset)%avs1**2)
         write(uout,'()')
      end do

      deallocate(var)

   end select

end subroutine ChainAver

!************************************************************************
!*                                                                      *
!*     NetworkAver                                                      *
!*                                                                      *
!************************************************************************

! ... calculate averages of network quantities

subroutine NetworkAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NetworkAver'
   character(80), parameter :: txheading ='network quantities'
   integer(4)   , save      :: ntype ! number of types of properties is being set below
   integer(4)   , save      :: nvar
   type(scalar_var), allocatable, save :: var(:)
   type(networkprop_var) :: NetworkProperty
   integer(4) :: inw, inwt, ioffset
   real(8) :: Asphericity

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      ntype = merge(12,11,lweakcharge) ! ntype equals 12 only if lweakcharge, else ntype = 11

      nvar = nnwt*ntype
      allocate(var(nvar))

      ! ... please note, that properties of type 1-7 are "rms" quantities. They have to be normalized by "ScalarNorm" with
      ! ... iopt = 2. If new properties are wished to be incorporated, please adjust the call of "ScalarNorm" below for
      ! ... both stages, "IAfterMacrostep" and "IAfterSimulation".

      do inwt = 1, nnwt
         ioffset = ntype*(inwt-1)
         var(1+ioffset)%label  = '<r(g)**2>**0.5                 = ' ! rms radius of gyration
         var(2+ioffset)%label  = '<r(g)**2_x>**0.5               = ' ! rms radius of gyration projection on the x-axis
         var(3+ioffset)%label  = '<r(g)**2_y>**0.5               = ' ! rms radius of gyration projection on the y-axis
         var(4+ioffset)%label  = '<r(g)**2_z>**0.5               = ' ! rms radius of gyration projection on the z-axis
         var(5+ioffset)%label  = 'smallest rms mom. p.a.         = ' ! smallest rms moment along a prinical axis
         var(6+ioffset)%label  = 'intermediate rms mom. p.a.     = ' ! intermediate rms moment along a prinical axis
         var(7+ioffset)%label  = 'largest rms mom. p.a.          = ' ! largest rms moment along a prinical axis
         var(8+ioffset)%label  = '<asphericity>                  = ' ! asphericity
         var(9+ioffset)%label  = '<xtheta>                       = ' ! angle of axes of largest extension and x-axes of main frame
         var(10+ioffset)%label = '<ytheta>                       = ' ! angle of axes of largest extension and y-axes of main frame
         var(11+ioffset)%label = '<ztheta>                       = ' ! angle of axes of largest extension and z-axes of main frame
         if (lweakcharge) &
         var(12+ioffset)%label = '<alpha>                        = ' ! degree of ionization
         var(1+ioffset:ntype+ioffset)%norm = One/nnwnwt(inwt)
     end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var%value = Zero
      do inw = 1, nnw
         inwt = inwtnwn(inw)
         ioffset = ntype*(inwt-1)
         call CalcNetworkProperty(inw, NetworkProperty)
         var(1+ioffset)%value  = var(1+ioffset)%value  + NetworkProperty%rg2
         var(2+ioffset)%value  = var(2+ioffset)%value  + NetworkProperty%rg2x
         var(3+ioffset)%value  = var(3+ioffset)%value  + NetworkProperty%rg2y
         var(4+ioffset)%value  = var(4+ioffset)%value  + NetworkProperty%rg2z
         var(5+ioffset)%value  = var(5+ioffset)%value  + NetworkProperty%rg2s
         var(6+ioffset)%value  = var(6+ioffset)%value  + NetworkProperty%rg2m
         var(7+ioffset)%value  = var(7+ioffset)%value  + NetworkProperty%rg2l
         var(8+ioffset)%value  = var(8+ioffset)%value  + NetworkProperty%asph
         var(9+ioffset)%value  = var(9+ioffset)%value  + NetworkProperty%theta(1)
         var(10+ioffset)%value = var(10+ioffset)%value + NetworkProperty%theta(2)
         var(11+ioffset)%value = var(11+ioffset)%value + NetworkProperty%theta(3)
         if (lweakcharge) &
         var(12+ioffset)%value = var(12+ioffset)%value + NetworkProperty%alpha
      end do
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call WriteHead(2, txheading, uout)
      do inwt = 1, nnwt
         ioffset = ntype*(inwt-1)
         call ScalarNorm(iStage, 1+ioffset, 7+ioffset, var, 2)
         call ScalarNorm(iStage, 8+ioffset, ntype+ioffset, var, 0)
         if (nnwt > 1) write(uout,'(2a)') 'network: ', txnwt(inwt)
         if (nnwt > 1) write(uout,'(a)')  '------------------'
         call ScalarWrite(iStage, 1+ioffset, ntype+ioffset, var, 1, '(a,t35,4f15.5,f15.0)', uout)
         write(uout,'()')
         write(uout,'(a,t35,2f15.5)') 'asphericity (<2:nd moments>)   = ', &
           Asphericity(var(5+ioffset)%avs2**2,var(6+ioffset)%avs2**2,var(7+ioffset)%avs2**2)
         write(uout,'()')
      end do

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call WriteHead(2, txheading, uout)
      do inwt = 1, nnwt
         ioffset = ntype*(inwt-1)
         call ScalarNorm(iStage, 1+ioffset, 7+ioffset, var, 2)
         call ScalarNorm(iStage, 8+ioffset, ntype+ioffset, var, 0)
         if (nnwt > 1) write(uout,'(2a)') 'network: ', txnwt(inwt)
         if (nnwt > 1) write(uout,'(a)')  '------------------'
         call ScalarWrite(iStage, 1+ioffset, ntype+ioffset, var, 1, '(a,t35,4f15.5,f15.0)', uout)
         write(uout,'()')
         write(uout,'(a,t35,2f15.5)') 'asphericity (<2:nd moments>)   = ', &
           Asphericity(var(5+ioffset)%avs1**2,var(6+ioffset)%avs1**2,var(7+ioffset)%avs1**2)
         write(uout,'()')
      end do

      deallocate(var)

   end select

end subroutine NetworkAver

!***********************************************************************
!*                                                                      *
!*     HierarchicalAver                                                 *
!*                                                                      *
!************************************************************************

! ... calculate averages of chain quantities

subroutine HierarchicalAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='HierarchicalAver'
   character(80), parameter :: txheading ='hierarchical quantities'
   integer(4)      , save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   real(8) :: rg

   if (slave) return   ! master only

   if (nh > 1) return  ! routine only adapted for nh = 1 (one hierarchical structure)

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      nvar = 1
      allocate(var(nvar))

      var(1)%label = '<r(g)**2>**0.5                 = ' ! rms radius of gyration
      var(1)%norm = One

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      call HierarchicalRg(rg)
      var(1)%value = rg
      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,f15.5,f15.0)', uout)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

      deallocate(var)

   end select

contains

!........................................................................

subroutine HierarchicalRg(vsumr)

! ... assumes that ipnhn and nphn are apporiate (contigeously ip)
! ... size of hierarchical structure should be smaller than half box size

   real(8), intent(out) :: vsumr    ! rms radius of gyration
   real(8) :: xcom, ycom, zcom, dx, dy, dz, r2, rref(1:3)
   integer(4) :: ip

   rref(1:3) = ro(1:3,ipnhn)
   do ip = ipnhn, ipnhn + nphn - 1
      dx = ro(1,ip) - rref(1)
      dy = ro(2,ip) - rref(2)
      dz = ro(3,ip) - rref(3)
      call PBC(dx,dy,dz)
      vaux(1,ip) = rref(1) + dx
      vaux(2,ip) = rref(2) + dy
      vaux(3,ip) = rref(3) + dz
   end do

   xcom = sum(vaux(1,ipnhn:ipnhn+nphn-1))/nphn
   ycom = sum(vaux(2,ipnhn:ipnhn+nphn-1))/nphn
   zcom = sum(vaux(3,ipnhn:ipnhn+nphn-1))/nphn
   vsumr = Zero
   do ip = ipnhn, ipnhn + nphn - 1
      dx = vaux(1,ip)-xcom
      dy = vaux(2,ip)-ycom
      dz = vaux(3,ip)-zcom
      r2 = dx**2+dy**2+dz**2
      vsumr = vsumr + r2
   end do
   vsumr = sqrt(vsumr/nphn)

end subroutine HierarchicalRg

!........................................................................

end subroutine HierarchicalAver

!************************************************************************
!*                                                                      *
!*     ThermoInteg                                                      *
!*                                                                      *
!************************************************************************

! ... handle themodynamic integration

subroutine ThermoInteg(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ThermoInteg'
   character(80), parameter :: txheading ='thermodynamic integration'
   real(8),          save :: lambda
   real(8),          save :: powercharge
   real(8),          save :: powerkbond
   real(8),          save :: powerkangle
   integer(4),       save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   real(8) :: data(9)

   namelist /nmlThermoInteg/ lambda, powercharge, powerkbond, powerkangle

   if (slave) return    ! master only

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      powercharge = One
      powerkbond = 9.0d0
      powerkangle = 6.0d0

      rewind(uin)
      read(uin,nmlThermoInteg)

      zat(1:nat) = lambda**powercharge*zat(1:nat)
      az(1:na) = lambda**powercharge*az(1:na)
      bond(1:nct)%k = lambda**powerkbond*bond(1:nct)%k
      angle(1:nct)%k = lambda**powerkangle*angle(1:nct)%k

      nvar = 4
      allocate(var(nvar))

      var(1)%label = '<dUtot/dlambda>                = '
      var(2)%label = '<dUtwobody/dlambda>            = '
      var(3)%label = '<dUbond/dlambda>               = '
      var(4)%label = '<dUangle/dlambda>              = '
      var(1)%norm = One/np
      var(2)%norm = One/np
      var(3)%norm = One/np
      var(4)%norm = One/np

   case (iWriteInput)

      call WriteHead(2, txheading, uout)
      write(uout,'(a)') 'Only applicable for systems with (i) hard-core, (ii) Coulomb, (iii) bond, and/or (iv) angle potential'
      write(uout,*)
      write(uout,'(a,t35,f10.2,a)') 'lambda                         = ', lambda
      write(uout,'(a,t35,f10.2,a)') 'power of lambda (charge)       = ', powercharge , '   (e -> e*lambda**power)'
      write(uout,'(a,t35,f10.2,a)') 'power of lambda (bond%p)       = ', powerkbond  , '   (bond%k -> bond%k*lambda**power)'
      write(uout,'(a,t35,f10.2,a)') 'power of lambda (angle%p)      = ', powerkangle , '   (angle%k -> angle%k*lambda**power)'

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      if (lambda > Zero) then
         var(2)%value = Two*powercharge*(u%twob(0)+u%rec)/lambda
         var(3)%value = powerkbond*u%bond/lambda
         var(4)%value = powerkangle*u%angle/lambda
         var(1)%value = sum(var(2:4)%value)
      else
         var(1:4)%value = Zero
      end if
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

      data(1)     = lambda                                       ! lambda
      data(2:3)   = [ np*beta*var(1)%avs1, np*beta*var(1)%avsd ] ! <dUtot/dlambda>/kT
      data(4:5)   = [ np*beta*var(2)%avs1, np*beta*var(2)%avsd ] ! <dUtwobody/dlambda>/kT
      data(6:7)   = [ np*beta*var(3)%avs1, np*beta*var(3)%avsd ] ! <dUbond/dlambda>/kT
      data(8:9)   = [ np*beta*var(4)%avs1, np*beta*var(4)%avsd ] ! <dUangle/dlambda>/kT
      call WriteVecAppend(9, data, 'thermointeg.data')

      deallocate(var)

   end select

end subroutine ThermoInteg

!************************************************************************
!*                                                                      *
!*     DistFunc                                                         *
!*                                                                      *
!************************************************************************

! ... calculate distribution functions

!     type  label  distribution functions
!     ----  -----  ----------------------
!     1     totu   total potential energy, u%tot
!     2     paru   partial potential energy, u
!     3     bindu  binding energy, ui(ip)
!     4     pairu  pair energy   , upp(i,j)
!     5     rdf    radial d.f., particle-particle (center of mass)
!     6     rdf    radial d.f., atom(mass > maslim)-atom(mass > maslim)
!     7     idm    induced dipole moment: total
!     8     idm    induced dipole moment: particle
!     9     idm    induced dipole moment: atom (polarization > pollim)
!    10     zdens  z-density distribution function: particle

subroutine DistFunc(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DistFunc'
   character(80), parameter :: txheading ='distribution functions'
   integer(4)   , parameter :: ntype = 10
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),                    save :: maslim, pollim
   real(8),                    save :: rcutdist, rcutdist2
   integer(4),                 save :: idist, itestdist
   integer(4),    allocatable, save :: ixat(:), ixatat(:,:)
   type(potenergy_var) :: usave
   integer(4) :: itype, ivar, ibin, idi, ibuf
   integer(4) :: ip, ipt, jp, jpt, iptjpt, ia, ialoc, ialow, iat, iatloc, ja, jaloc, jalow, jat, iatjat
   integer(4) :: iatx, natx, iatjatx, natatx
   real(8), allocatable, save :: ubind(:), gcont(:)
   real(8)    :: dx, dy, dz, dropbc(3), d, usum, fsum, virtwob, virmolecule
   real(8)    :: raa2, raa1, rpp2, rpp1, dipm, norm, dvol, vols2, utobdypp

   integer(4), parameter :: nvar_s = 2
   type(scalar_var), save :: var_s(nvar_s)

   namelist /nmlDist/ vtype, idist, rcutdist, maslim, pollim, itestdist

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

#if defined (_PAR_)
! ... to get right position of ' comm'
      if (ltime) call CpuAdd('start', 'comm', 1, uout)
      if (ltime) call CpuAdd('stop', 'comm', 1, uout)
#endif

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
!     vtype%lable= ['totu ','paru ','bindu','pairu','rdf  ','rdf  ','idm  ','idm  ','idm  ','zdens']
      vtype%min  = [-100.d0,-100.d0,-100.d0, -25.d0, Zero  , Zero  , Zero , Zero  , Zero  , -boxlen(3)]
      vtype%max  = [  Zero ,  Zero ,  Zero ,  25.d0, 10.d0 , 10.d0 , 0.5d0, 0.25d0, 0.25d0,  boxlen(3)]
      vtype%nbin = 100
      idist = 10
      rcutdist = rcut + Two*maxval(racom)
      maslim = 1.0d-4
      pollim = 1.0d-4
      itestdist = 0

      rewind(uin)
      read(uin,nmlDist)

! ... check condition

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)
      if (lmc .and. (rcutdist < rcut)) call Warn(txroutine, 'lmc .and. (rcutdist < rcut): distfunc considers smaller region then energy evaluation',uout)

   case (iWriteInput)

! ... set rcutdist2

      rcutdist2 = rcutdist**2

      if (.not.allocated(ixat)) then
         allocate(ixat(nat), ixatat(nat,nat))
         ixat = 0
         ixatat = 0
      end if

! ... set ixatat (for selecting atom-atom pairs with masses > maslim)

      natatx = 0
      do iat = 1, nat
         do jat = iat, nat
            if (massat(iat) >= maslim .and. massat(jat) >= maslim) then
               natatx = natatx+1
               ixatat(iat,jat) = natatx
               ixatat(jat,iat) = natatx
            else
               ixatat(iat,jat) = 0
               ixatat(jat,iat) = 0
            end if
         end do
      end do

! ... set ixat (for selecting atoms with polarizability > pollim)

      natx = 0
      iat = 0
      do ipt = 1, npt
         do iatloc = 1, natpt(ipt)
            iat = iat+1
            if (maxval(abs(poltensa(1:6,iatloc,ipt))) > pollim) then
               natx = natx+1
               ixat(iat) = natx
            else
               ixat(iat) = 0
            end if
         end do
      end do

! ... set remaining elements of vtype

      vtype%label = ['totu ','paru ','bindu','pairu','rdf  ','rdf  ','idm  ','idm  ','idm  ','zdens']
      vtype%nvar = [   1   , nptpt,     npt,  nptpt,  nptpt, natatx,      1, npt   , natx  ,    npt]

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(natat,ntype), ubind(np_alloc))
      ipnt = 0
      ubind = 0.0E+00

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            if (itype == 1) then
               ivar = ivar+1
               ipnt(1,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            else if (itype == 2) then
               do ipt = 1, npt
                  do jpt = ipt, npt
                     ivar = ivar+1
                     ipnt(iptpt(ipt,jpt),itype) = ivar
                     var(ivar)%label = trim(vtype(itype)%label)//' '//txptpt(iptpt(ipt,jpt))
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
                  end do
               end do
            else if (itype == 3) then
               do ipt = 1, npt
                  ivar = ivar+1
                  ipnt(ipt,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txpt(ipt)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            else if (itype == 4) then
               do ipt = 1, npt
                  do jpt = ipt, npt
                     ivar = ivar+1
                     ipnt(iptpt(ipt,jpt),itype) = ivar
                     var(ivar)%label = trim(vtype(itype)%label)//' '//txptpt(iptpt(ipt,jpt))
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
                  end do
               end do
            else if (itype == 5) then
               do ipt = 1, npt
                  do jpt = ipt, npt
                     ivar = ivar+1
                     ipnt(iptpt(ipt,jpt),itype) = ivar
                     var(ivar)%label = trim(vtype(itype)%label)//' '//txptpt(iptpt(ipt,jpt))
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
                   end do
               end do
            else if (itype == 6) then
               do iat = 1, nat
                  do jat = iat, nat
                     iatjatx = ixatat(iat,jat)
                     if (iatjatx > 0) then
                     ivar = ivar+1
                     ipnt(iatjatx,itype) = ivar
                     var(ivar)%label = trim(vtype(itype)%label)//' '//txatat(iatat(iat,jat))
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
                     end if
                  end do
               end do
            else if (itype == 7) then
               ivar = ivar+1
               ipnt(1,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' total'
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            else if (itype == 8) then
               do ipt = 1, npt
                  ivar = ivar+1
                  ipnt(ipt,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txpt(ipt)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            else if (itype == 9) then
               do iat = 1, nat
                  iatx = ixat(iat)
                  if (iatx > 0) then
                     ivar = ivar+1
                     ipnt(iatx,itype) = ivar
                     var(ivar)%label = trim(vtype(itype)%label)//' '//txat(iat)
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
                  end if
               end do
            else if (itype == 10) then
               do ipt = 1, npt
                  ivar = ivar+1
                  ipnt(ipt,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txpt(ipt)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            end if
         end if
      end do
      if (ivar /= nvar) call Stop(txroutine, 'ivar /= nvar', uout)

      var_s(1)%label = 'volume'
      var_s(2)%label = 'pressure (hs)'
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar_s, var_s)
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) var_s

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar_s, var_s)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1

      var_s(1)%value = vol
      call ScalarSample(iStage,1,1,var_s)              ! volume

! ... sampling of type 1

      itype = 1
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         ibin = max(-1,min(floor(var(ivar)%bini*(u%tot/np-var(ivar)%min)),var(ivar)%nbin))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
      end if

! ... sampling of type 2

      itype = 2
      if (vtype(itype)%l) then
         do idi = 1, vtype(itype)%nvar
            ivar = ipnt(idi,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(u%twob(idi)/np-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end if

      if (mod(istep2,idist) == 0) then

      force(1:3,1:na) = Zero
      virial = Zero

! ............... loop over all particle pairs ..........

      ubind(1:np) = Zero
      virtwob = Zero
      do ip = ipmyid(1), ipmyid(2)
         ipt = iptpn(ip)
         ialow = ianpn(ip)

! ... sampling of type 10

         itype = 10
         if (vtype(itype)%l) then
            ivar = ipnt(ipt,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(r(3,ip)-var(ivar)%min)),var(ivar)%nbin))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

#if defined (_PAR_)
         do jp = 1, np                                     ! primitive loadbalance
            if (ip > jp .and. mod(ip+jp,2) == 0) cycle
            if (ip < jp .and. mod(ip+jp,2) /= 0) cycle
            if (ip == jp) cycle
#else
         do jp = ip+1, np
#endif
            utobdypp = Zero
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            ia = ialow-1
            jalow = ianpn(jp)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC2(dx,dy,dz,dropbc(1),dropbc(2),dropbc(3))
            dx = dx-dropbc(1)
            dy = dy-dropbc(2)
            dz = dz-dropbc(3)
            rpp2 = dx**2+dy**2+dz**2
            if (rpp2 > rcutdist2) cycle

! ... sampling of type 5

            itype = 5
            if (vtype(itype)%l) then
               rpp1 = sqrt(rpp2)
               ivar = ipnt(iptpt(ipt,jpt),itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(rpp1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ............... loop over all atom pairs ..........

            do ialoc = 1, napt(ipt)
               ia = ia+1
               iat = iatan(ia)
               ja = jalow-1
               do jaloc = 1, napt(jpt)
                  ja = ja+1
                  jat = iatan(ja)
                  iatjat = iatat(iat,jat)
                  dx = r(1,ia)-r(1,ja)-dropbc(1)
                  dy = r(2,ia)-r(2,ja)-dropbc(2)
                  dz = r(3,ia)-r(3,ja)-dropbc(3)
                  raa2 = dx**2+dy**2+dz**2
                  raa1 = sqrt(raa2)
                  if (rpp2 > rcut2) goto 455

                  ibuf = iubuflow(iatjat)
                  do
                     if (raa2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                     if (ibuf > nbuf) call Stop(txroutine, 'ibuf > nbuf', uout)
                  end do
                  d = raa2-ubuf(ibuf)
                  usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                                    d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
                  fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                                   d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))

                  force(1,ia) = force(1,ia) + (fsum * dx)
                  force(2,ia) = force(2,ia) + (fsum * dy)
                  force(3,ia) = force(3,ia) + (fsum * dz)
                  force(1,ja) = force(1,ja) - (fsum * dx)
                  force(2,ja) = force(2,ja) - (fsum * dy)
                  force(3,ja) = force(3,ja) - (fsum * dz)
                  virtwob     = virtwob     - (fsum*raa2)
                  utobdypp = utobdypp +usum
                  ubind(ip) = ubind(ip)+usum
                  ubind(jp) = ubind(jp)+usum

                  if (master .and. itestdist == 1) call TestDistFunc1(uout)

  455             continue

! ... sampling of type 6

                  itype = 6
                  if (vtype(itype)%l) then
                     iatjat = ixatat(iat,jat)
                     if (iatjat > 0) then
                        ivar = ipnt(iatjat,itype)
                        ibin = max(-1,min(floor(var(ivar)%bini*(raa1-var(ivar)%min)),int(var(ivar)%nbin)))
                        var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
                     end if
                  end if

               end do
            end do

! ... sampling of type 4

            itype = 4
            if (vtype(itype)%l) then
               ivar = ipnt(iptpt(ipt,jpt),itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(utobdypp-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

         end do

      end do

      if (lmc .or. lmcall .or. lbd) virial = virial + virtwob

! ............... end of looping over all atom pairs ...............

! ... add contributions to forces and the virial from other interactions

      if (lmc .or. lmcall .or. lbd) then
         if (ltime) call CpuAdd('stop', txroutine, 0, uout)
         if (ltime) call CpuAdd('interrupt', ' ', 0, uout)
         usave = u                   ! save potential energies
         if (lcharge .and. lewald) call UEwald
         if (lweakcharge .and. lewald) call UEwald
         if (ldipole .and. lewald) call UDipoleEwald
         if (lchain) call UBond
         if (lchain) call UAngle
         if (lclink) call UCrossLink
!        if (ldieldis) call UDielDis !  no forces are included in UDielDis
         if (luext)  call UExternal
         u = usave                   ! restore potential energies
         if (ltime) call CpuAdd('resume', ' ', 0, uout)
         if (ltime) call CpuAdd('start', txroutine, 0, uout)
      end if

#if defined (_PAR_)
! ... allreduce of ubind, force, and virial
      if (ltime) call CpuAdd('start', 'comm', 0, uout)
      call par_allreduce_reals(ubind, vaux, np  )
      call par_allreduce_reals(force, vaux, 3*na)
      call par_allreduce_real(virial, raux)
      if (ltime) call CpuAdd('stop', 'comm', 0, uout)
#endif

      if (master) then   ! could be parallelized later

! ... sampling of type 3

         itype = 3
         if (vtype(itype)%l) then
            do ip = 1, np
               ipt = iptpn(ip)
               ivar = ipnt(ipt,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(ubind(ip)-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end do
         end if

! ... sampling of type 7

         itype = 7
         if (vtype(itype)%l) then
            ivar = ipnt(1,itype)
            dipm = sqrt(idmsys(1)**2+idmsys(2)**2+idmsys(3)**2)
            ibin = max(-1,min(floor(var(ivar)%bini*(dipm-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sampling of type 8

         itype = 8
         if (vtype(itype)%l) then
            do ip = 1, np
               ipt = iptpn(ip)
               ivar = ipnt(ipt,itype)
               dipm = sqrt(idmo(1,ip)**2+idmo(2,ip)**2+idmo(3,ip)**2)
               ibin = max(-1,min(floor(var(ivar)%bini*(dipm-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end do
         end if

! ... sampling of type 9

         if (vtype(9)%l) then
            do ia = 1, na
               iat = iatan(ia)
               iatx = ixat(iat)
               if (iatx > 0) then
                  ivar = ipnt(iatx,9)
                  dipm = sqrt(idm(1,ia)**2+idm(2,ia)**2+idm(3,ia)**2)
                  ibin = max(-1,min(floor(var(ivar)%bini*(dipm-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if
            end do
         end if

      end if

! ... calculate forces, torques, virial, and  pressure

      virmolecule      = Zero
!     forceo(1:3,1:np) = Zero
!     torqueo(1:3,1:np) = Zero
!     fac = (sclene/scllen)/sfor
      do ia = 1, na
         ip = ipnan(ia)
         dx = r(1,ia)-ro(1,ip)
         dy = r(2,ia)-ro(2,ip)
         dz = r(3,ia)-ro(3,ip)
         virmolecule = virmolecule + dx*force(1,ia) + dy*force(2,ia) + dz*force(3,ia)
!        force(1,ia) = force(1,ia) * fac
!        force(2,ia) = force(2,ia) * fac
!        force(3,ia) = force(3,ia) * fac
!        forceo(1,ip)  = forceo(1,ip)  + force(1,ia)
!        forceo(2,ip)  = forceo(2,ip)  + force(2,ia)
!        forceo(3,ip)  = forceo(3,ip)  + force(3,ia)
!        torqueo(1,ip) = torqueo(1,ip) + dy*force(3,ia)-dz*force(2,ia) + torque(1,ia)*fac
!        torqueo(2,ip) = torqueo(2,ip) + dz*force(1,ia)-dx*force(3,ia) + torque(2,ia)*fac
!        torqueo(3,ip) = torqueo(3,ip) + dx*force(2,ia)-dy*force(1,ia) + torque(3,ia)*fac
      end do
      virial = virial + virmolecule
      if (lmc .or. lmcall .or. lbd) prsr = (np*Boltz*temp*scltem-virial*sclene/(Three*AvNo))/(vol*sclvol)/sclpre

      if (master .and. itestdist == 1) call TestDistFunc2(uout)

      end if

   case (iAfterMacrostep)

! ... reduce var%avs2 to master

#if defined (_PAR_)
         if (ltime) call CpuAdd('start', 'comm', 0, uout)
         do ivar = 1, nvar
            call par_reduce_reals(var(ivar)%avs2(-1), vaux, mnbin_df+2)
         end do
         if (ltime) call CpuAdd('stop', 'comm', 0, uout)
#endif

      if (master) then

         call ScalarSample(iStage,1,1,var_s)           ! volume
         vols2 = var_s(1)%avs2                        ! extract volume averaged over a macrostep

! ............... normalize and sum distribution functions .............

! ... normalize

         do itype = 1, 4
            if (vtype(itype)%l) then
               call DistFuncNorm(ipnt(1,itype), ipnt(vtype(itype)%nvar,itype), var)
            end if
         end do
         if (vtype(5)%l) then
            do ipt = 1, npt
               do jpt = ipt, npt
                  ivar = ipnt(iptpt(ipt,jpt),5)
                  norm = idist*vols2/(FourPiThird*nppt(ipt)*nppt(jpt))
                  if (ipt == jpt) norm = norm*Two
                  do ibin = 0, var(ivar)%nbin
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
                  end do
               end do
            end do
         end if
         if (vtype(6)%l) then
            do iat = 1, nat
               ipt = iptat(iat)
               do jat = iat, nat
                  jpt = iptat(jat)
                  iatjat = ixatat(iat,jat)
                  if (iatjat > 0) then
                     ivar = ipnt(iatjat,6)
                     norm = idist*vols2/(FourPiThird*nppt(ipt)*naat(iat)*nppt(jpt)*naat(jat))
                     if (ipt == jpt .and. iat == jat) norm = norm*Two
                     do ibin = 0, var(ivar)%nbin
                        var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
                     end do
                  end if
               end do
            end do
         end if
         do itype = 7, 9
            if (vtype(itype)%l) then
               call DistFuncNorm(ipnt(1,itype), ipnt(vtype(itype)%nvar,itype), var)
            end if
         end do
         do itype = 10, 10
            if (vtype(itype)%l) then
               do ipt = 1, npt
                  ivar = ipnt(ipt,10)
                  norm = idist/(boxlen(1)*boxlen(2)*var(ivar)%bin)
                  do ibin = 0, var(ivar)%nbin
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm
                  end do
               end do
            end if
         end do

         call DistFuncSample(iStage, nvar, var)

! ... calculate contract contribution to the reduced pressure

         itype = 5
         if (lmonoatom .and. vtype(itype)%l) then
            call PressureContact(ipnt(1,itype), var, vols2, var_s(2)%value)
            call ScalarSample(iSimulationStep, 2, 2, var_s)     ! sample
            call ScalarSample(iStage, 2, 2, var_s)
         end if

         if (lsim .and. master) write(ucnf) var
         if (lsim .and. master) write(ucnf) var_s

      end if

   case (iAfterSimulation)

      if (master) then

         call DistFuncSample(iStage, nvar, var)
         call ScalarSample(iStage, 1, 1, var_s)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,i15)')    'sampling interval              = ', idist
         write(uout,'(a,t35,i15)')    'number of steps/passes used    = ', nstep1*nstep2/idist
         write(uout,'(a,t35,2f15.2)') 'cutoff distance                = ', rcutdist
         write(uout,'(a,t35,2f15.5)') 'mass limit                     = ', maslim
         write(uout,'(a,t35,2f15.5)') 'polarization limit             = ', pollim
         write(uout,'()')
         write(uout,'(a,t35,2es15.5)')'average volume                 = ', var_s(1)%avs1, var_s(1)%avsd

         if (lmonoatom .and. vtype(5)%l) then     ! get contact contribution to reduced pressure
            allocate(gcont(nptpt))
            gcont = 0.0E+00
            do iptjpt = 1, nptpt
               ivar = ipnt(iptjpt,5)
               call gcontact(var(ivar)%min, var(ivar)%bin, var(ivar)%avs1(-1), r1atat(iptjpt), gcont(iptjpt))
            end do
            call ScalarSample(iStage,2,2,var_s)
            write(uout,'()')
            write(uout,'(a,t35,10f15.2)') 'contact values                 = ', gcont
            write(uout,'()')
            write(uout,'(a,t35,2f15.5)') 'pressure (NkT/V)(contact cont.)= ', var_s(2)%avs1, var_s(2)%avsd
            write(uout,'(a,t35,2f15.5)') 'pressure (NkT/V)(total)        = ', &
             prsrreds3+var_s(2)%avs1, sqrt(prsrredsd**2+var_s(2)%avsd**2)
            write(uout,'(a,t35,2es15.5)')'pressure (kT) (total)          = ', &
            np/var_s(1)%avs1*(prsrreds3+var_s(2)%avs1), np/var_s(1)%avs1*sqrt(prsrredsd**2+var_s(2)%avsd**2)
         end if

         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

! ... calculate excess ammount

         itype = 10
         if (vtype(itype)%l) call ExcessAmount(vtype(itype)%nbin, npt, txpt, ipnt(1,itype), var, boxlen, uout)

      end if

      deallocate(var, ipnt, ixat, ixatat, ubind)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine TestDistFunc1(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1',unit)
   write(unit,'(a,2i5)')    'ia, ja       ', ia, ja
   write(unit,'(a,3g15.7)') 'force(1:3,ia)', force(1:3,ia)
   write(unit,'(a,3g15.7)') 'force(1:3,ja)', force(1:3,ja)
   write(unit,'(a,1g15.7)') 'accum virtwob', virtwob
end subroutine TestDistFunc1

!........................................................................

subroutine TestDistFunc2(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a,g15.7)') 'virial(molecular correction)       ', virmolecule
   write(unit,'(a,g15.7)') 'virial(after molecular correction) ', virial
end subroutine TestDistFunc2

!........................................................................

end subroutine DistFunc

!************************************************************************
!*                                                                      *
!*     PressureContact                                                  *
!*                                                                      *
!************************************************************************

! ... calculate contract contribution to the reduced pressure

subroutine PressureContact(ipnt, var, vols2, prsrhs)

   use MolModule
   implicit none

   integer(4), parameter :: npol = 2, ndp = npol+1

   integer(4),   intent(in)  :: ipnt(*)
   type(df_var), intent(in)  :: var(*)
   real(8),      intent(in)  :: vols2
   real(8),      intent(out) :: prsrhs

   integer(4) :: ipt, jpt, iptjpt, ivar
   real(8)    :: denst, densipt, densjpt, rcont, gcont, term

   prsrhs = Zero
   do ipt = 1, npt
      densipt = nppt(ipt)/vols2
      do jpt = 1, npt
         densjpt = nppt(jpt)/vols2
         iptjpt = iptpt(ipt,jpt)
         ivar = ipnt(iptjpt)
         rcont = radat(ipt)+radat(jpt)
         call gcontact(var(ivar)%min, var(ivar)%bin, var(ivar)%avs2(-1), rcont, gcont)
         term = rcont**3*densipt*densjpt*gcont
         prsrhs = prsrhs+term
      end do
   end do
   denst = sum(nppt(1:npt))/vols2
   prsrhs = TwoPi/(Three*denst)*prsrhs

end subroutine PressureContact

!************************************************************************
!*                                                                      *
!*     GContact                                                         *
!*                                                                      *
!************************************************************************

! ... calculate rdf contact value

subroutine GContact(vlow, bin, vs2, rcont, gcont)

   implicit none

   integer(4), parameter :: npol = 2, ndp = npol+1

   real(8), intent(in)    :: vlow
   real(8), intent(in)    :: bin
   real(8), intent(in)    :: vs2(-1:*)
   real(8), intent(in)    :: rcont
   real(8), intent(out)   :: gcont

   integer(4) :: i, ihs
   real(8)    :: xx(ndp), yy(ndp), ww(ndp), a(0:npol), dum1, dum2
   real(8), external :: PolVal

   ihs = int((rcont-vlow-1.0d-10)/bin)
   xx(1:ndp) = vlow+bin*(0.5d0+[ (ihs+i,i = 1,ndp) ] )
   yy(1:ndp) = [ (vs2(ihs+i),i = 1,ndp) ]
!  write(*,'(a,4f10.5)') 'vs2(ihs+1:ihs+ndp)', vs2(ihs+1:ihs+ndp)
!  write(*,'(a,4f10.5)') 'yy(1:3)', yy(1:3)
   ww(1:ndp) = 1.0d0
   if (count(yy(1:ndp) < 1.0d-10) > 0) then     ! Zero contribution or ruff rdf
      gcont = 0.0d0
    else if (yy(1) > 0.5d0) then                 ! logarithmic transformation before fit
      yy(1:ndp) = log(yy(1:ndp))
      call PolFit(npol, ndp, xx, yy, ww, 0, 6, a, dum1, dum2)
      gcont = exp(PolVal(npol,a,rcont))
   else
      call PolFit(npol, ndp, xx, yy, ww, 0, 6, a, dum1, dum2)
      gcont = PolVal(npol,a,rcont)
   end if
   gcont = max(0.0d0,gcont)
!  write(*,'(a,4f10.5)') 'vs2(ihs+1:ihs+ndp)', vs2(ihs+1:ihs+ndp)
!  write(*,'(a,4f10.5)') 'yy(1:3), gcont', yy(1:3), gcont

end subroutine GContact

!************************************************************************
!*                                                                      *
!*     TestSimulation                                                   *
!*                                                                      *
!************************************************************************

! ... write test output

subroutine TestSimulation

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='TestSimulation'
   character(34), parameter :: fmt_3  = '(i4,t10,a,t25,i6,t35,a,t50,3f12.6)'
   character(34), parameter :: fmt_6  = '(i4,t10,a,t25,i6,t35,a,t50,6f12.6)'
   character(22), parameter :: fmt3   = '(i4,t10,a,t23,3es12.4)'
   character(33), parameter :: fmt33  = '(i4,t10,a,t23,3es12.4,2x,3es12.4)'
   character(33), parameter :: fmt44  = '(i4,t10,a,t23,4es12.4,2x,4es12.4)'
   character(44), parameter :: fmt333 = '(i4,t10,a,t23,3es12.4,2x,3es12.4,2x,3es12.4)'
   integer(4) :: ip, jp, jplow, jpupp, ia, m, ncolum, nparts
   real(8) :: dx, dy, dz

   if (slave) return

   call WriteHead(3, txroutine, uout)

! ... energies, pressure, and volume

   if (lmd)    write(uout,'(a,t25,i12)') 'time step     ', istep2
   if (lmc)    write(uout,'(a,t25,i12)') 'configuration ', ipass+np*(istep2-1)
   if (lmcall) write(uout,'(a,t25,i12)') 'configuration ', istep2
   write(uout,'(a,t25,i12)') 'seed of random generator ', iseed
   if (lmd)    write(uout,'(a,t25,f12.5)') 'total energy', (ekin+u%tot)/np
   if (lmd)    write(uout,'(a,t25,f12.5)') 'kinetic energy', ekin/np
   write(uout,'(a,t25,f12.5)')  'tot. pot. energy', u%tot/np
   write(uout,'(a,t25,f12.5)')  'tot. two-body pot. energy', u%twob(0)/np
   write(uout,'(a,t25,6f12.5)') 'two-body pot. energy', u%twob(1:nptpt)/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%rec', u%rec/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%stat', u%stat/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%pol', u%pol/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%bond', u%bond/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%angle', u%angle/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%crosslink', u%crosslink/np
   write(uout,'(a,t25,f12.5)')  'pot. energy, u%external', u%external/np
   write(uout,'(a,t25,f12.5)')  'pressure', prsr
   write(uout,'(a,t25,f12.5)')  'volume', vol

! ... atom data

   write(uout,*)
   write(uout,'(a)') 'atom data'
   write(uout,'(a)') 'atom coordinates'
   write(uout,fmt_3 ) (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),r(1:3,ia),ia = 1,na)
   if (lmd) then
      write(uout,'(a)') 'forces and torques (lab frame)'
      write(uout,fmt33 ) (ia,txat(iatan(ia)),force(1:3,ia),torque(1:3,ia),ia = 1,na)
   end if
   if (ldipole) then
      write(uout,'(a)') 'atom dipole moments'
      write(uout,fmt_3 ) (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),dip(1:3,ia),ia = 1,na)
   end if
   if (lpolarization) then
      write(uout,'(a)') 'atom dipole moments'
      write(uout,fmt_3 ) (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),dip(1:3,ia),ia = 1,na)
      write(uout,'(a)') 'atom polarizabilities (xx,yy,zz,xy,xz,yz)'
      write(uout,fmt_6 ) (ia,txat(iatan(ia)),ipnan(ia),txpt(iptan(ia)),poltens(1:6,ia),ia = 1,na)
   end if

! ... particle data

   write(uout,*)
   write(uout,'(a)') 'particle data'
   write(uout,'(a)') 'particle position'
   write(uout,fmt3  ) (ip,txpt(iptpn(ip)),ro(1:3,ip),ip = 1,np)
   write(uout,'(a)') 'particle orientation'
   write(uout,fmt333) (ip,txpt(iptpn(ip)),ori(1:3,1:3,ip),ip = 1,np)
   if (lmd) call TestSimulationMD(uout)

   write(uout,'(a)') 'particle distance matrix'
   ncolum = 10
   nparts = 1+(np-1)/ncolum
   do m = 1, nparts
      jplow = 1+ncolum*(m-1)
      jpupp = min(np,jplow-1+ncolum)
      write(uout,'(4x,10i10)') (jp,jp = jplow,jpupp)
      do ip = 1, np
         do jp = jplow, jpupp
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC(dx,dy,dz)
            vaux(1,jp) = sqrt(dx**2+dy**2+dz**2)
         end do
         write(uout,'(i4,10f10.4)') ip, vaux(1,jplow:jpupp)
      end do
   end do

end subroutine TestSimulation

