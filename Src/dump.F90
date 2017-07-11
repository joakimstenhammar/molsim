! ... 'version 6.3.4, Sep 18, 2015'

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
!*     DumpModule                                                       *
!*                                                                      *
!************************************************************************

! ... module for dumping

module DumpModule

   use MolModule

   logical       :: ldpos, ldori, ldliv, ldanv, ldfor, ldtor, ldidm, ldlaz, ldutot, ldumpuser ! logical flag for dumping
   integer(4)    :: iplow                                                                     ! lower particle to be dumped
   integer(4)    :: ipupp                                                                     ! upper particle to be dumped
   integer(4)    :: ialow                                                                     ! lower atom to be dumped
   integer(4)    :: iaupp                                                                     ! upper atom to be dumped

! ... external units


end module DumpModule

!************************************************************************
!*                                                                      *
!*     DumpDriver                                                       *
!*                                                                      *
!************************************************************************

! ... dump driver

subroutine DumpDriver(iStage)

   use DumpModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DumpDriver'

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      call IODump(iStage)
      call DoDump('open')
      if (ldumpuser) call DumpUser(iStage)

   case (iWriteInput)

      call IODump(iStage)
      if (ldumpuser) call DumpUser(iStage)

   case (iBeforeSimulation)

      if (lsim) then
         if (txstart /= 'continue') call DoDump('write')
         if (txstart == 'continue') call DoDump('advance')
      end if
      if (lana) call DoDump('read')
      if (ldumpuser) call DumpUser(iStage)

   case (iSimulationStep)

      if (lsim) then
         if (mod(istep2,idump) == 0) then
            call DoDump('write')
            if (ldumpuser) call DumpUser(iStage)
         end if
      end if
      if (lana) then
         call DoDump('read')
         if (ldumpuser) call DumpUser(iStage)
      end if

   case (iAfterSimulation)

      call DoDump('close')
      if (ldumpuser) call DumpUser(iStage)

   end select

end subroutine DumpDriver

!************************************************************************
!*                                                                      *
!*     IODump                                                           *
!*                                                                      *
!************************************************************************

! ... performing i/o on dump variables

subroutine IODump(iStage)

   use DumpModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='IODump'
   character(20), save :: txptdump
   integer(4)          :: ipt

   namelist /nmlDump/ idump, txptdump, ldpos, ldori, ldliv, ldanv, ldfor, ldtor, ldidm, ldlaz, ldutot, ldumpuser
   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      idump = 10
      txptdump = 'all'
      ldpos  = .false.
      ldori  = .false.
      ldliv  = .false.
      ldanv  = .false.
      ldfor  = .false.
      ldtor  = .false.
      ldidm  = .false.
      ldlaz  = .false.
      ldutot = .false.
      ldumpuser = .false.

      rewind(uin)
      read(uin,nmlDump)

      call LowerCase(txptdump)

   case (iWriteInput)

! ... determine iplow, ipupp, ialow, and iaupp

      if (txptdump == 'all') then
         iplow = 1
         ipupp = np
         ialow = 1
         iaupp = na
      else
         iplow = 0
         do ipt = 1, npt
            if (txptdump == txpt(ipt)) then
               iplow = ipnpt(ipt)
               ipupp = ipnpt(ipt) + nppt(ipt)-1
               ialow = ianpn(iplow) + nppt(ipt)*napt(ipt)
            end if
         end do
         if (iplow == 0) call Stop(txroutine, 'check variable txptdump', uout)
      end if

! ... check some conditions

      if (idump <= 0) call Stop(txroutine, 'idump <= 0', uout)
      if (mod(nstep2,idump) /= 0) call Stop(txroutine, 'mod(nstep2,idump) /= 0', uout)
      if (lana .and. (txptdump /= 'all')) call Stop(txroutine, 'lana .and. (txptdump /= ''all'')', uout)

! ... write input data

      if (master) then
         call WriteHead(2, 'dumping data', uout)
         write(uout,'(a,t35,i5)') 'dumping interval               = ', idump
         write(uout,'(a,t35,a )') 'particle type to be dumped     = ', txptdump
         write(uout,'()')
         write(uout,'(a,t20,a)') 'quantity', 'external unit'
         write(uout,'(a,t20,a)') '--------', '-------------'
         if (ldpos)     write(uout,'(a,a)') 'positions              ', trim(fpos)
         if (ldori)     write(uout,'(a,a)') 'orientations           ', trim(fori)
         if (ldliv)     write(uout,'(a,a)') 'linear velocities      ', trim(fliv)
         if (ldanv)     write(uout,'(a,a)') 'angular velocities     ', trim(fanv)
         if (ldfor)     write(uout,'(a,a)') 'forces                 ', trim(ffor)
         if (ldtor)     write(uout,'(a,a)') 'torques                ', trim(ftor)
         if (ldidm)     write(uout,'(a,a)') 'ind dip mom            ', trim(fidm)
         if (ldlaz)     write(uout,'(a,a)') 'atom charge state      ', trim(flaz)
         if (ldutot)    write(uout,'(a,a)') 'potential energy       ', trim(futot)
         if (ldumpuser) write(uout,'(a,a)') 'user dump              ', trim(fuser)
      end if

   end select

end subroutine IODump

!************************************************************************
!*                                                                      *
!*     DoDump                                                           *
!*                                                                      *
!************************************************************************

! ... performing dumping matters

subroutine DoDump(str)

   use DumpModule
   implicit none

   character(*), intent(in) :: str

   integer(4)   :: ip, ia, m, idum = 0
   real(8)      :: dum
   logical      :: ldum

   if (str(1:4) == 'open' .and. master) then

      if (ldpos) call FileOpen(upos, fpos, 'unform/noread')
      if (ldori) call FileOpen(uori, fori, 'unform/noread')
      if (ldliv) call FileOpen(uliv, fliv, 'unform/noread')
      if (ldanv) call FileOpen(uanv, fanv, 'unform/noread')
      if (ldfor) call FileOpen(ufor, ffor, 'unform/noread')
      if (ldtor) call FileOpen(utor, ftor, 'unform/noread')
      if (ldidm) call FileOpen(uidm, fidm, 'unform/noread')
      if (ldlaz) call FileOpen(ulaz, flaz, 'unform/noread')
      if (ldutot) call FileOpen(uutot, futot, 'form/noread')

   else if (str == 'advance' .and. master) then

      do m = 1, (nstep2/idump)*(nstep1beg-1)+1
         if (ldpos) read(upos) idum, dum, dum, dum
         if (ldpos) read(upos) (dum,dum,dum,ip = iplow,ipupp)
         if (ldori) read(uori) (dum,dum,dum,dum,ip = iplow,ipupp)
         if (ldliv) read(uliv) (dum,dum,dum,ip = iplow,ipupp)
         if (ldanv) read(uanv) (dum,dum,dum,ip = iplow,ipupp)
         if (ldfor) read(ufor) (dum,dum,dum,ip = iplow,ipupp)
         if (ldtor) read(utor) (dum,dum,dum,ip = iplow,ipupp)
         if (ldidm) read(uidm) (dum,dum,dum,ip = iplow,ipupp)
         if (ldlaz) read(ulaz) (ldum,ia = ialow,iaupp)
         if (ldutot) read(uutot,*) dum
      end do

   else if (str == 'read') then

      if (master) then
         if (ldpos) read(upos) idum, boxlen(1:3)
         if (ldpos) read(upos) ro(1:3,iplow:ipupp)
         if (ldori) read(uori) qua(0:3,iplow:ipupp)
         if (ldliv) read(uliv) rod(1:3,iplow:ipupp)
         if (ldanv) read(uanv) angvelo(1:3,iplow:ipupp)
         if (ldfor) read(ufor) forceo(1:3,iplow:ipupp)
         if (ldtor) read(utor) torqueo(1:3,iplow:ipupp)
         if (ldidm) read(uidm) idmo(1:3,iplow:ipupp)
         if (ldlaz) then
            read(ulaz) laz(ialow:iaupp)
            where (laz(ialow:iaupp))
               az(ialow:iaupp) = zat(iatan(ialow:iaupp))
            elsewhere
               az(ialow:iaupp) = Zero
            end where
         end if
         if (ldutot) read(uutot,*) u%tot
         if (ldutot) u%tot = u%tot/(sclene/(np*GasConstant*temp*scltem))

      end if

#if defined (_PAR_)
      if (ldpos) call par_bc_reals(boxlen    , 3   )
      if (ldpos) call par_bc_reals(ro     , 3*(ipupp-iplow+1))
      if (ldori) call par_bc_reals(qua    , 4*(ipupp-iplow+1))
      if (ldliv) call par_bc_reals(rod    , 3*(ipupp-iplow+1))
      if (ldanv) call par_bc_reals(angvelo, 3*(ipupp-iplow+1))
      if (ldfor) call par_bc_reals(forceo , 3*(ipupp-iplow+1))
      if (ldtor) call par_bc_reals(torqueo, 3*(ipupp-iplow+1))
      if (ldidm) call par_bc_reals(idmo   , 3*(ipupp-iplow+1))
      if (ldlaz) call par_bc_reals(az     ,   (ipupp-iplow+1))
      if (ldlaz) call par_bc_logicals(laz ,   (ipupp-iplow+1))
      if (ldutot) call par_bc_real(u%tot)
      ! previously it was
      ! if (ldutot) call par_bc_reals(u%tot*sclene/(np*GasConstant*temp*scltem), 1)
      ! which does not work, as the output of the routine is written into a constant (u%tot*sclene/(...) )
#endif

      call QuaToOri(np, iplow, ipupp, qua, ori)
      call SetAtomProp(iplow, ipupp, .false.)
      call SetBoxParam

   else if (str == 'write' .and. master) then

      if (lmc .or. lmcall) call OriToQua(np, iplow, ipupp, ori, qua)
      if (ldpos) then
         write(upos) idum, boxlen(1:3)
         write(upos) ro(1:3,iplow:ipupp)
         call FileFlush(upos)
      end if
      if (ldori) then
         write(uori) qua(0:3,iplow:ipupp)
         call FileFlush(uori)
      end if
      if (ldliv) then
         write(uliv) rod(1:3,iplow:ipupp)
         call FileFlush(uliv)
      end if
      if (ldanv) then
         write(uanv) angvelo(1:3,iplow:ipupp)
         call FileFlush(uanv)
      end if
      if (ldfor) then
         write(ufor) forceo(1:3,iplow:ipupp)
         call FileFlush(ufor)
      end if
      if (ldtor) then
         write(utor) torqueo(1:3,iplow:ipupp)
         call FileFlush(utor)
      end if
      if (ldidm) then
         write(uidm) idmo(1:3,iplow:ipupp)
         call FileFlush(uidm)
      end if
      if (ldlaz) then
         write(ulaz) laz(ialow:iaupp)
         call FileFlush(ulaz)
      end if
      if (ldutot) then
         write(uutot,*) u%tot*sclene/(np*GasConstant*temp*scltem)
         call FileFlush(uutot)
      end if

   else if (str == 'close' .and. master) then

      if (ldpos) close (upos)
      if (ldori) close (uori)
      if (ldliv) close (uliv)
      if (ldanv) close (uanv)
      if (ldfor) close (ufor)
      if (ldtor) close (utor)
      if (ldidm) close (uidm)
      if (ldlaz) close (ulaz)
      if (ldutot) close (uutot)

   end if

end subroutine DoDump


