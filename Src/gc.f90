!************************************************************************
!> \page gc_f90 gc.f90
!! **NPartChange**
!! *perform one number of particle change trial move*
!************************************************************************


!     note, not check to work with weak charges

subroutine NPartChange(iStage)

   use MCModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NPartChange'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap, lWarnHCOverlap, lDelete
   integer(4) :: ip, ipt, jp, jpt, iptjpt, iploc
   real(8) :: prandom, uuu, fforce(3), weight, Random

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   call CpuAdd('start', txroutine, 1, uout)

   imovetype = inpartmove

! ... deletion or insertion?

   lDelete = .true.
   prandom = Random(iseed)
   if (prandom > half) lDelete = .false.

   du%twob = Zero
   lboxoverlap = .false.
   lhsoverlap = .false.
   lhepoverlap = .false.

! ... calculate energy difference

   if (lDelete) then    ! trial delete of particle ipmove

      do jp = 1, np
         jpt = iptpn(jp)
         iptjpt = iptpt(iptmove, jpt)
         if (jp == ipmove) cycle
         call UTwoBodyPair(ipmove, jp, uuu, fforce)
         du%twob(iptjpt) = du%twob(iptjpt) - uuu    ! note, minus sign
      end do
      du%twob(0) = sum(du%twob(1:nptpt))
      weight = nppt(iptmove)/vol * exp(-beta*chempot(iptmove))

   else                 ! trial insert a particle of type iptmove

      nppt(iptmove) = nppt(iptmove) + 1            ! add one particle of type iptmove

      if (sum(nppt(1:npt)) > np_alloc) call stop(txroutine, 'sum(nppt(1:npt)) > np_alloc', uout)

      call SetObjectParam1
      call SetObjectParam2
      ipmove = np                                  ! place the new particle last
      call SetPartPosRandomMC(ipmove)              ! set trial random particle position
      call SetPartOriRandom(iseed,ori(1,1,ipmove)) ! set trial random particle orientation
      call SetAtomProp(ipmove,ipmove,.false.)      ! trial atom properties

      do jp = 1, np-1
         jpt = iptpn(jp)
         iptjpt = iptpt(iptmove, jpt)
         call UTwoBodyPair(ipmove, jp, uuu, fforce)
         if (uuu >= 1d10) then
            lhsoverlap = .true.
            goto 400
         end if
         du%twob(iptjpt) = du%twob(iptjpt) + uuu    ! note, plus sign
      end do
      du%twob(0) = sum(du%twob(1:nptpt))
      weight = vol/nppt(iptmove)*exp(beta*chempot(iptmove))
400   continue

   end if

   du%tot = du%twob(0)                            ! could be generalized

   if (itestmc == 3) call TestNPartChange1(uout)

!  ............. decide new configuration .............

   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   if (ievent == imcaccept) then        ! accepted change of number of partciles

      if (lDelete) then
         nppt(iptmove) = nppt(iptmove) - 1
         ro(1:3,ipmove) = ro(1:3,np)
         ori(1:3,1:3,ipmove) = ori(1:3,1:3,np)
         call SetObjectParam1
         call SetObjectParam2
         call SetAtomProp(ipmove,ipmove,.false.)      ! trial atom properties
      else
      end if
      call SetVList                                   ! update neighbour list
      u%tot           = u%tot           + du%tot
      u%twob(0:nptpt) = u%twob(0:nptpt) + du%twob(0:nptpt)

   else                                ! rejected change of number of paricles

      if (lDelete) then
      else
         nppt(iptmove) = nppt(iptmove) - 1
         call SetObjectParam1
         call SetObjectParam2
      end if

   end if
   if (itestmc == 3) call TestNPartChange2(uout)

   call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine TestNPartChange1(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1', unit)
   write(unit,'(a,t20,i4    )') 'ipmove              ', ipmove
   write(unit,'(a,t20,i4    )') 'iptmove             ', iptmove
   write(unit,'(a,t20,l     )') 'lDelete             ', lDelete
   write(unit,'(a,t20,l     )') 'lhsoverlap          ', lhsoverlap
   write(unit,'(a,t20,3f15.5)') 'volume              ', vol
   write(unit,'(a,t20,i4    )') 'nppt(iptmove)       ', nppt(iptmove)
   write(unit,'(a,t20,3f15.5)') 'weight              ', weight
   write(unit,'(a,t20,3f15.5)') 'du%tot              ', du%tot
   write(unit,'(a,t20,3f15.5)') 'du%twob(0:npt)     ', du%twob(0:npt)
   if (.not. lDelete) then
   write(unit,'(a,t20,3f15.5)') 'trial coordinates   ', ro(1:3,ipmove)
   end if
end subroutine TestNPartChange1

!........................................................................

subroutine TestNPartChange2(unit)
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a,t20,i4    )') 'istep1              ', istep1
   write(unit,'(a,t20,i4    )') 'istep2              ', istep2
   write(unit,'(a,t20,i4    )') 'ipass               ', ipass
   write(unit,'(a,t20,i4    )') 'ievent              ', ievent
   write(unit,'(a,t20,i4    )') 'np                  ', np
   write(unit,'(a,t20,i4    )') 'nppt(iptmove)       ', nppt(iptmove)
   call WritePartAtomProp
end subroutine TestNPartChange2

!........................................................................

end subroutine NPartChange

