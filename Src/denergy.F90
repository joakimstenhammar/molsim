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

! relation among energy routines
!
!    DUTotal
!       !   lcharge
!       !--------------
!       !             !
!       !             !--------- DUTwoBody
!       !             !              !
!       !             !              !----- DUTwoBody(A/ALList/P)New
!       !             !              !
!       !             !              !------DUTwoBody(A/ALList/P)Old
!       !             !  lewald
!       !             !--------- DUTwoBodyEwald
!       !
!       !   lweakcharge
!       !--------------
!       !             !
!       !             !--------- DUTwoBody
!       !             !              !
!       !             !              !----- DUWeakCharge(A/P)New
!       !             !              !
!       !             !              !------DUWeakCharge(A/P)Old
!       !             !  lewald
!       !             !--------- DUWeakChargeEwald
!       !
!       !   ldipole
!       !--------------
!       !             !
!       !             !--------- DUTwoBody
!       !             !              !
!       !             !              !----- DUTwoBodyPNew
!       !             !              !
!       !             !              !------DUTwoBodyPOld
!       !             !
!       !             !--------- DUDipole
!       !                            !
!       !                            !---------- UDipolePNew
!       !                            !
!       !                            !---------- UDipolePOld
!       !                            !   lewald
!       !                            !---------- DUDipoleEwald
!       !   ldipolesph
!       !--------------
!       !             !
!       !             !--------- DUDipoleSph
!       !                            !
!       !                            !----- DUDipoleSphNew
!       !                            !
!       !                            !------DUDipoleSphOld
!       !   ldieldis
!       !--------------
!       !             !
!       !             !---------- DUDielDis
!       !   lchain
!       !---------- DUBond
!       !
!       !   lchain
!       !---------- DUAngle
!       !
!       !   lclink
!       !---------- DUCrossLink
!       !
!       !   lext
!       !---------- DUExternal

!************************************************************************
!> \page denergy denergy.F90
!! **DUTotal**
!! *calculate energy difference between two configurations*
!************************************************************************


!     old configuration is given by ro, ori, r for all particles
!     new configuration is given by ro, ori, r except for moving particles
!     where rotm, oritm, and rtm are used

subroutine DUTotal(lhsoverlap,lhepoverlap)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap
                                                ! =.false. no hard-core overlap
   logical,    intent(out) :: lhepoverlap       ! =.true. hard-external-potential overlap
                                                ! =.false. no hard-external-potential overlap
   character(40), parameter :: txroutine ='DUTotal'

   external UTwoBodyANew, UTwoBodyAOld
   external UTwoBodyANewLList, UTwoBodyAOldLList
   external UTwoBodyANewCellList, UTwoBodyAOldCellList
   external UTwoBodyPNew, UTwoBodyPOld
   external UWeakChargeANew, UWeakChargeAOld
   external UWeakChargePNew, UWeakChargePOld

   if (ltrace) call WriteTrace(4, trim(txroutine), iSimulationStep)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(.not.allocated(utwobnew)) then
      allocate(utwobnew(0:nptpt), utwobold(0:nptpt))
      utwobnew = 0.0E+00
      utwobold = 0.0E+00
   end if

! ............... initiate ...............

   du%tot   = Zero
   du%twob  = Zero
   du%oneb  = Zero
   du%rec   = Zero
   du%stat  = Zero
   du%pol   = Zero
   du%bond  = Zero
   du%angle = Zero
   du%crosslink   = Zero
   du%external    = Zero
   lhsoverlap  = .true.             ! temporary fix to satisfy sensitive compilers

! .............. select appropiate energy routines ............

   if (lcharge) then                  ! atoms possessing charges

      if (lmonoatom) then
         if (lvlist) call DUTwoBody(lhsoverlap, UTwoBodyANew, UTwoBodyAOld)
         if (lllist) call DUTwoBody(lhsoverlap, UTwoBodyANewLList, UTwoBodyAOldLList)
         if (lCellList) call DUTwoBody(lhsoverlap, UTwoBodyANewCellList, UTwoBodyAOldCellList)
      else
         call DUTwoBody(lhsoverlap, UTwoBodyPNew, UTwoBodyPOld)
      end if
      if (lhsoverlap) goto 400
      if (lewald) call DUTwoBodyEwald

   else if (lweakcharge) then  ! atoms possesing weak charges (titrating system)

      if (lmonoatom) then
         call DUTwoBody(lhsoverlap, UWeakChargeANew, UWeakChargeAOld)
      else
         call DUTwoBody(lhsoverlap, UWeakChargePNew, UWeakChargePOld)
      end if

      if (lhsoverlap) goto 400

      if (lewald) call DUWeakChargeEwald

if (itest == 90) then
      call writehead(3,txroutine, uout)                             !cc
      write(uout,'(a,100l4)') 'laz', laz(1:na)                       !cc
      write(uout,'(a,10f10.5)') 'du%twob(0:nptpt)', du%twob(0:nptpt)!cc
end if

   else if (ldipole) then                     ! atoms possessing charges and dipoles

      call DUTwoBody(lhsoverlap, UTwoBodyPNew, UTwoBodyPOld)
      if (lhsoverlap) goto 400

      call DUDipole

   else if (ldipolesph) then                  ! spherical b.c., atoms possessing charges and dipoles, image charge

      call DUDipoleSph(lhsoverlap)
      if (lhsoverlap) goto 400

   else if (ldieldis) then                  ! atoms possesing charge in a system with dielectric discontinuities

      call DUDielDis(lhsoverlap)
      if (lhsoverlap) goto 400

   end if

   if (lchain) call DUBond
   if (lchain) call DUAngle
   if (lclink) call DUCrossLink

   lhepoverlap = .false.
   if (luext) call DUExternal(lhepoverlap)
   if (lhepoverlap) goto 400

400 continue

#if defined (_PAR_)
   if (ltrace) call WriteTrace(4, trim(txroutine)//'_start_Pack', iSimulationStep)
   call PackReduceU(nptpt+1, npt+1, du%tot, du%twob, du%oneb, du%rec, du%stat, du%pol, du%bond, du%angle, du%crosslink, du%external, uaux)
   call par_allreduce_logical(lhsoverlap, laux)
   call par_allreduce_logical(lhepoverlap, laux)
   if (ltrace) call WriteTrace(4, trim(txroutine)//'_after_Pack', iSimulationStep)
#endif

!  call TestDUTotal

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

!........................................................................

contains

subroutine TestDUTotal
   call WriteHead(3, 'Test'//trim(txroutine), uout)
   write(uout,'(a, l    )') 'lhsoverlap  ' ,lhsoverlap
   write(uout,'(a, l    )') 'lhepoverlap ' ,lhepoverlap
   write(uout,'(a,3f14.4)') 'du%tot      ' ,du%tot
   write(uout,'(a,8f14.4)') 'du%twob     ' ,du%twob(0:nptpt)
   write(uout,'(a,8f14.4)') 'du%oneb     ' ,du%oneb(0:npt)
!   write(uout,'(a,8f14.4)') 'du%rec      ' ,du%rec
!   write(uout,'(a,8f14.4)') 'du%stat     ' ,du%stat
!   write(uout,'(a,8f14.4)') 'du%pol      ' ,du%pol
!   write(uout,'(a,8f14.4)') 'du%bond     ' ,du%bond
!   write(uout,'(a,8f14.4)') 'du%angle    ' ,du%angle
!   write(uout,'(a,8f14.4)') 'du%crosslink' ,du%crosslink
!   write(uout,'(a,8f14.4)') 'du%external ' ,du%external
end subroutine TestDUTotal

!........................................................................

end subroutine DUTotal

!**********************************************************************************************************************

!************************************************************************
!> \page denergy denergy.F90
!! **DUTwoBody**
!! *calculate two-body potential energy difference*
!************************************************************************


subroutine DUTwoBody(lhsoverlap, utwobodynew, twobodyold)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap
                                                ! =.false. no hard-core overlap

   character(40), parameter :: txroutine ='DUTwoBody'
   integer(4) :: jp
   external utwobodynew
   external twobodyold

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   du%twob(0:nptpt) = Zero                      ! initiate

   call utwobodynew(lhsoverlap,jp)                ! calculate new two-body potential energy

#if defined (_PAR_)
   call par_allreduce_logical(lhsoverlap, laux)
#endif

   if (lhsoverlap) goto 400                     ! check hard-core overlap

   call twobodyold                               ! calculate old two-body potential energy

   du%tot = du%tot + du%twob(0)                ! update

400 continue

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

end subroutine DUTwoBody

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyANew**
!! *calculate two-body potential energy for new configuration*
!************************************************************************

!     only monoatomic particles

subroutine UTwoBodyANew(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='UTwoBodyANew'

   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

!   write(uout,*) txroutine

   utwobnew(0:nptpt) = Zero
   lhsoverlap =.true.

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
!      write(uout,'(a,i5,3f10.5)') 'ip,rotm(1:3,iploc)',ip, rotm(1:3,iploc)

      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
         end if
         if (r2 > rcut2) cycle
         if (r2 < r2atat(iptjpt)) goto 400

         if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                           d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         utwobnew(iptjpt) = utwobnew(iptjpt) + usum
!         write(uout,'(a,i5,6f10.5)') 'jp,ro(1:3,jp), r2, usum, utwobnew(iptjpt)',jp, ro(1:3,jp), r2, usum, utwobnew(iptjpt)
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-rotm(1,jploc)
            dy = rotm(2,iploc)-rotm(2,jploc)
            dz = rotm(3,iploc)-rotm(3,jploc)
            call PBCr2(dx,dy,dz,r2)
            if (lellipsoid) Then
                if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc))) goto 400
            end if
            if (r2 > rcut2) cycle
            if (r2 < r2atat(iptjpt)) goto 400

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         end do
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)
   lhsoverlap =.false.

  400 continue

end subroutine UTwoBodyANew

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyAOld**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     only monoatomic particles

subroutine UTwoBodyAOld

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UTwoBodyAOld'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   utwobold(0:nptpt) = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
!      write(uout,'(a,i5,3f10.5)') 'ip,ro(1:3,ip)',ip, ro(1:3,ip)
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 > rcut2) cycle
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
         end if

         if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                           d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         utwobold(iptjpt) = utwobold(iptjpt) + usum
!        write(uout,'(a,i5,6f10.5)') 'jp,ro(1:3,jp), r2, usum, utwobold(iptjpt)',jp, ro(1:3,jp), r2, usum, utwobold(iptjpt)
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if
            if (r2 < r2atat(iptjpt)) goto 400

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobold(iptjpt) = utwobold(iptjpt) + usum
         end do
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

end subroutine UTwoBodyAOld

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyANewLList**
!! *calculate two-body potential energy for new configuration*
!************************************************************************

!     only monoatomic particles
!     linked list version

subroutine UTwoBodyANewLList(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='UTwoBodyANewLList'
   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf, icell, Getncellllist
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   utwobnew(0:nptpt) = Zero
   lhsoverlap =.true.

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)

      call SetNCell(rotm(1:3,iploc), rcut, lcellllist)
      do icell = 1, Getncellllist()
         if (.not.lcellllist(icell)) cycle
         jp = headllist(icell)
         do
            if (jp == 0) exit
            if (lptm(jp) .or. ip == jp) goto 101
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-ro(1,jp)
            dy = rotm(2,iploc)-ro(2,jp)
            dz = rotm(3,iploc)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) goto 101
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if
            if (r2 < r2atat(iptjpt)) goto 400

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
 101        continue
            jp = jpllist(jp)
         end do
      end do
   end do

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-rotm(1,jploc)
            dy = rotm(2,iploc)-rotm(2,jploc)
            dz = rotm(3,iploc)-rotm(3,jploc)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if
            if (r2 < r2atat(iptjpt)) goto 400

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         end do
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)
   lhsoverlap =.false.

  400 continue

end subroutine UTwoBodyANewLList

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyAOldLList**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     only monoatomic particles
!     linked list version

subroutine UTwoBodyAOldLList

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UTwoBodyAOldLList'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf, icell, Getncellllist
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   utwobold(0:nptpt) = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)

      call SetNCell(ro(1:3,ip), rcut, lcellllist)
      do icell = 1, Getncellllist()
         if (.not.lcellllist(icell)) cycle
         jp = headllist(icell)
         do
            if (jp == 0) exit
            if (lptm(jp) .or. ip == jp) goto 102
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) goto 102
            if (lellipsoid) Then   ! Also here?
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                                d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobold(iptjpt) = utwobold(iptjpt) + usum
 102        continue
            jp = jpllist(jp)
         end do
      end do
   end do

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if
            if (r2 < r2atat(iptjpt)) goto 400

            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            utwobold(iptjpt) = utwobold(iptjpt) + usum
         end do
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

end subroutine UTwoBodyAOldLList

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyANewCellList**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     only monoatomic particles
!     cell list version

subroutine UTwoBodyANewCellList(lHsOverlap, ipOverlap)

   use MolModule,      only: ro, rotm, iptpn, iptpt, rcut2, r2atat, Zero
   use MolModule,      only: ipnptm, nptm, lptm, lptmdutwob, nptpt
   use MolModule,      only: nproc, myid
   use MolModule,      only: lellipsoid, lsuperball, oritm, ori, radellipsoid2, aellipsoid
   use MolModule,      only: du, iubuflow, ubuf, r2umin
   use EnergyModule,   only: utwobnew
   use CellListModule, only: pcellro, cell_type, ipnext

   implicit none
   character(40), parameter :: txroutine ='UTwoBodyANewCellList'
   logical, intent(out)     :: lHsOverlap
   integer(4), intent(out)  :: ipOverlap

   integer(4) :: ip, jp, ipt, jpt, iptjpt, iploc, jploc, ibuf
   real(8)    :: dr(3), r2, d, usum

   type(cell_type), pointer :: icell, ncell
   integer(4)               :: incell
   logical    :: EllipsoidOverlap, SuperballOverlap

   utwobnew(0:nptpt) = Zero
   lHsOverlap = .true.
   ipOverlap = 1

! ... contribution from pairs where only one particle is moved
   do iploc = 1, nptm
      ip    =  ipnptm(iploc)
      ipt   =  iptpn(ip)
      icell => pcellro(rotm(1:3,iploc))
      do incell = 1 + myid, icell%nneighcell, nproc ! increment with nproc to have parallel execution
         ncell => icell%neighcell(incell)%p
         jp = ncell%iphead
         do jploc = 1, ncell%npart
            ipOverlap = jp
            if (.not. lptm(jp)) then
               jpt = iptpn(jp)
               iptjpt = iptpt(ipt,jpt)
               dr(1:3) = rotm(1:3,iploc)-ro(1:3,jp)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               if (lellipsoid) Then
                  if (EllipsoidOverlap(r2,dr,oritm(:,:,iploc),ori(:,:,jp),radellipsoid2,aellipsoid)) return
               end if
               if (lsuperball) Then
                  if (SuperballOverlap(r2,dr,oritm(:,:,iploc),ori(:,:,jp))) return
               end if
               if (r2 < rcut2) then
                  if (r2 < r2atat(iptjpt)) return ! Hard Sphere overlap
                  if (r2 < r2umin(iptjpt)) return ! outside lower end of tabulated potential energy table

                  ibuf = iubuflow(iptjpt)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = ubuf(ibuf+1) + d*(ubuf(ibuf+2) + d*(ubuf(ibuf+3) + d*(ubuf(ibuf+4) + d*(ubuf(ibuf+5) + d*ubuf(ibuf+6)))))
                  utwobnew(iptjpt) = utwobnew(iptjpt) + usum
               end if
            end if
            jp = ipnext(jp)
        end do
      end do
   end do

! ... contribution from pairs where both particle is moved

   if (lptmdutwob) then
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc + 1 + myid, nptm, nproc ! increment with nproc for parallel simulations
            jp = ipnptm(jploc)
            ipOverlap = jp
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dr(1:3) = rotm(1:3,iploc)-rotm(1:3,jploc)
            call PBCr2(dr(1), dr(2), dr(3), r2)
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,dr,oritm(:,:,iploc),oritm(:,:,jploc),radellipsoid2,aellipsoid)) return
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,dr,oritm(:,:,iploc),oritm(:,:,jploc))) return
            end if
            if (r2 < rcut2) then
               if (r2 < r2atat(iptjpt)) return ! Hard Sphere overlap
               if (r2 < r2umin(iptjpt)) return ! outside lower end

               ibuf = iubuflow(iptjpt)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                                 d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

               utwobnew(iptjpt) = utwobnew(iptjpt) + usum
            end if
         end do
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)
   lhsoverlap =.false.

end subroutine UTwoBodyANewCellList

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyAOldCellList**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     only monoatomic particles
!     cell list version

subroutine UTwoBodyAOldCellList

   use MolModule,      only: ro, iptpn, iptpt, rcut2, Zero
   use MolModule,      only: ipnptm, nptm, lptm, lptmdutwob, nptpt
   use MolModule,      only: nproc, myid
   use MolModule,      only: du, iubuflow, ubuf
   use EnergyModule,   only: utwobold
   use CellListModule, only: cell_type, cellip, ipnext

   implicit none
   character(40), parameter :: txroutine ='UTwoBodyAOldCellList'

   integer(4) :: ip, jp, ipt, jpt, iptjpt, iploc, jploc, ibuf
   real(8)    :: dr(3), r2, d, usum

   type(cell_type), pointer :: icell, ncell
   integer(4)               :: incell

   utwobold(0:nptpt) = Zero

! ... contribution from pairs where only one particle is moved
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      icell => cellip(ip)%p
      do incell = 1 + myid, icell%nneighcell, nproc ! increment with nproc to have parallel execution
         ncell => icell%neighcell(incell)%p
         jp = ncell%iphead
         do jploc = 1, ncell%npart
            if (.not. lptm(jp)) then
               jpt = iptpn(jp)
               iptjpt = iptpt(ipt,jpt)
               dr(1:3) = ro(1:3,ip)-ro(1:3,jp)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               !as the old configuration should be free of overlaps one can skip the checks
               if (r2 < rcut2) then
                  ibuf = iubuflow(iptjpt)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = ubuf(ibuf+1) + d*(ubuf(ibuf+2) + d*(ubuf(ibuf+3) + d*(ubuf(ibuf+4) + d*(ubuf(ibuf+5) + d*ubuf(ibuf+6)))))
                  utwobold(iptjpt) = utwobold(iptjpt) + usum
               end if
            end if
            jp = ipnext(jp)
         end do
      end do
   end do

! ... contribution from pairs where both particle is moved

   if (lptmdutwob) then
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc + 1 + myid, nptm, nproc ! increment with nproc for parallel simulations
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dr(1:3) = ro(1:3,ip)-ro(1:3,jp)
            call PBCr2(dr(1), dr(2), dr(3), r2)
            !as the old configuration should be free of overlaps one can skip the checks
            if (r2 < rcut2) then
               ibuf = iubuflow(iptjpt)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = ubuf(ibuf+1) + d*(ubuf(ibuf+2) + d*(ubuf(ibuf+3) + d*(ubuf(ibuf+4) + d*(ubuf(ibuf+5) + d*ubuf(ibuf+6)))))
               utwobold(iptjpt) = utwobold(iptjpt) + usum
            end if
         end do
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

end subroutine UTwoBodyAOldCellList

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyPNew**
!! *calculate two-body potential energy for new configuration*
!************************************************************************


subroutine UTwoBodyPNew(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf
   integer(4) :: ia, ialoc, ialow, iaupp, kialow, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum
   integer(4) :: jaloc, kjalow
   logical    :: EllipsoidOverlap, SuperballOverlap

   lhsoverlap =.true.
   utwobnew(0:nptpt) = Zero

   kialow = 0
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
         end if

         usum = Zero
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         ialoc = kialow
         do ia = ialow, iaupp
            iat = iatan(ia)
            ialoc = ialoc+1
            do ja = jalow, jaupp
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = rtm(1,ialoc)-r(1,ja)-dxopbc
               dy = rtm(2,ialoc)-r(2,ja)-dyopbc
               dz = rtm(3,ialoc)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               if (r2 < r2atat(iatjat)) goto 400

               if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                       d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            end do
         end do
         utwobnew(iptjpt) = utwobnew(iptjpt) + usum
      end do
      kialow = kialow+napt(ipt)
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      kialow = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
         kjalow = kialow+napt(ipt)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-rotm(1,jploc)
            dy = rotm(2,iploc)-rotm(2,jploc)
            dz = rotm(3,iploc)-rotm(3,jploc)
            call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
            dx = dx-dxopbc
            dy = dy-dyopbc
            dz = dz-dzopbc
            r2 = dx**2+dy**2+dz**2
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
            end if

            usum = Zero
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
            ialoc = kialow
            do ia = ialow, iaupp
               iat = iatan(ia)
               ialoc = ialoc+1
               jaloc = kjalow
               do ja = jalow, jaupp
                  jaloc = jaloc+1
                  jat = iatan(ja)
                  iatjat = iatat(iat,jat)
                  dx = rtm(1,ialoc)-rtm(1,jaloc)-dxopbc
                  dy = rtm(2,ialoc)-rtm(2,jaloc)-dyopbc
                  dz = rtm(3,ialoc)-rtm(3,jaloc)-dzopbc
                  r2 = dx**2+dy**2+dz**2
                  if (r2 < r2atat(iatjat)) goto 400

                  if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
                  ibuf = iubuflow(iatjat)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                          d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

               end do
            end do
            kjalow = kjalow+napt(jpt)
            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         end do
         kialow = kialow+napt(ipt)
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)

   lhsoverlap =.false.

  400 continue

end subroutine UTwoBodyPNew

!************************************************************************
!> \page denergy denergy.F90
!! **UTwoBodyPOld**
!! *calculate two-body potential energy for old configuration*
!************************************************************************


subroutine UTwoBodyPOld

   use EnergyModule
   implicit none

   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf
   integer(4) :: ia, ialoc, ialow, iaupp, kialow, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum
   integer(4) :: jaloc, kjalow
   logical    :: EllipsoidOverlap, SuperballOverlap

   utwobold(0:nptpt) = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         usum = Zero
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            iat = iatan(ia)
            do ja = jalow, jaupp
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2

               if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                       d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            end do
         end do
         utwobold(iptjpt) = utwobold(iptjpt) + usum
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      kialow = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
         kjalow = kialow+napt(ipt)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
            dx = dx-dxopbc
            dy = dy-dyopbc
            dz = dz-dzopbc
            r2 = dx**2+dy**2+dz**2
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
            end if

            usum = Zero
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
            ialoc = kialow
            do ia = ialow, iaupp
               iat = iatan(ia)
               ialoc = ialoc+1
               jaloc = kjalow
               do ja = jalow, jaupp
                  jaloc = jaloc+1
                  jat = iatan(ja)
                  iatjat = iatat(iat,jat)
                  dx = r(1,ia)-r(1,ja)-dxopbc
                  dy = r(2,ia)-r(2,ja)-dyopbc
                  dz = r(3,ia)-r(3,ja)-dzopbc
                  r2 = dx**2+dy**2+dz**2
                  if (r2 < r2atat(iatjat)) goto 400

                  if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
                  ibuf = iubuflow(iatjat)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                          d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

               end do
               kjalow = kjalow+napt(jpt)
            end do
            utwobold(iptjpt) = utwobold(iptjpt) + usum
         end do
         kialow = kialow+napt(ipt)
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

end subroutine UTwoBodyPOld

!************************************************************************
!> \page denergy denergy.F90
!! **DUTwoBodyEwald**
!! *calculate two-body potential energy difference; k-space*
!************************************************************************


subroutine DUTwoBodyEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUTwoBodyEwald'

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... initiate

   du%rec = Zero

! ... calculate

   if (txewaldrec == 'std') then
       call DUTwoBodyEwaldRecStd
       if (lewald2dlc) call DUTwoBodyEwaldRec2dlc
   else if (txewaldrec == 'spm') then
       call DUTwoBodyEwaldRecSPM
   end if
   call DUTwoBodyEwaldSelf
   if (lsurf .and. master) call DUTwoBodyEwaldSurf

! ... update

   du%rec = EpsiFourPi*du%rec
   du%tot = du%tot + du%rec

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine DUTwoBodyEwaldRecStd
   character(40), parameter :: txroutine ='UEwaldRecStd'
   integer(4) :: kn, nx, ny, nz, ia, ialoc, ikvec2
   real(8)    :: term, termnew, termold

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   call EwaldSetArrayTM

   kn = kvecoffmyid
   ikvec2 = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         ikvec2 = ikvec2+1
         if (ikvec2 < kvecmyid(1) .or. ikvec2 > kvecmyid(2)) cycle  ! parallelize over k-vectors
         do ialoc = 1, natm
            ia = ianatm(ialoc)
            eikyzm(ia)      = conjg(eiky(ia,ny))     *eikz(ia,nz)
            eikyzp(ia)      =       eiky(ia,ny)      *eikz(ia,nz)
            eikyzmtm(ialoc) = conjg(eikytm(ialoc,ny))*eikztm(ialoc,nz)
            eikyzptm(ialoc) =       eikytm(ialoc,ny) *eikztm(ialoc,nz)
         end do

         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            sumeikrtm(kn,1) = sumeikr(kn,1)
            sumeikrtm(kn,2) = sumeikr(kn,2)
            sumeikrtm(kn,3) = sumeikr(kn,3)
            sumeikrtm(kn,4) = sumeikr(kn,4)
            do ialoc = 1, natm
               ia = ianatm(ialoc)
               sumeikrtm(kn,1) = sumeikrtm(kn,1)+az(ia)*  &
                  (conjg(eikxtm(ialoc,nx))*eikyzmtm(ialoc) - conjg(eikx(ia,nx))*eikyzm(ia))
               sumeikrtm(kn,2) = sumeikrtm(kn,2)+az(ia)*  &
                  (conjg(eikxtm(ialoc,nx))*eikyzptm(ialoc) - conjg(eikx(ia,nx))*eikyzp(ia))
               sumeikrtm(kn,3) = sumeikrtm(kn,3)+az(ia)*  &
                        (eikxtm(ialoc,nx) *eikyzmtm(ialoc) -       eikx(ia,nx) *eikyzm(ia))
               sumeikrtm(kn,4) = sumeikrtm(kn,4)+az(ia)*  &
                        (eikxtm(ialoc,nx) *eikyzptm(ialoc) -       eikx(ia,nx) *eikyzp(ia))
            end do

            termnew = real(sumeikrtm(kn,1))**2 + aimag(sumeikrtm(kn,1))**2 + real(sumeikrtm(kn,2))**2 + aimag(sumeikrtm(kn,2))**2 &
                    + real(sumeikrtm(kn,3))**2 + aimag(sumeikrtm(kn,3))**2 + real(sumeikrtm(kn,4))**2 + aimag(sumeikrtm(kn,4))**2
            termold = real(sumeikr(kn,1))**2   + aimag(sumeikr(kn,1))**2   + real(sumeikr(kn,2))**2   + aimag(sumeikr(kn,2))**2 &
                    + real(sumeikr(kn,3))**2   + aimag(sumeikr(kn,3))**2   + real(sumeikr(kn,4))**2   + aimag(sumeikr(kn,4))**2
            term    = kfac(kn)*(termnew - termold)
            du%rec   = du%rec + term

         end do
      end do
   end do

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUTwoBodyEwaldRecStd

!........................................................................

subroutine DUTwoBodyEwaldRec2dlc
   implicit none

   character(40), parameter :: txroutine ='DUTwoBodyEwaldRec2dlc'

   integer(4) :: kn, nx, ny, ia, ialoc
   real(8)    :: term, termnew, termold
   real(8)    :: kx, ky, kp, sinhkpztm, coshkpztm, sinhkpz, coshkpz, termi

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

   call EwaldSetArray2dTM                 ! calculate sinkxt, coskxt, sinkyt, and coskyt for moving particles

   kn = 0
   do ny = 0, ncut2d
      ky = ny*TwoPiBoxi(2)
      do nx = 0, ncut2d
         kx = nx*TwoPiBoxi(1)
         if (nx**2+ny**2 > ncut2d2) cycle
         if (nx == 0 .and. ny == 0) cycle
         kn = kn+1
         kp = sqrt(kx**2+ky**2)            ! k parallel to the surface

         sumtrigtm(kn,1) = sumtrig(kn,1)
         sumtrigtm(kn,2) = sumtrig(kn,2)
         sumtrigtm(kn,3) = sumtrig(kn,3)
         sumtrigtm(kn,4) = sumtrig(kn,4)
         sumtrigtm(kn,5) = sumtrig(kn,5)
         sumtrigtm(kn,6) = sumtrig(kn,6)
         sumtrigtm(kn,7) = sumtrig(kn,7)
         sumtrigtm(kn,8) = sumtrig(kn,8)
         do ialoc = 1, natm
            ia = ianatm(ialoc)
            term  = exp(kp*rtm(3,ialoc))
            termi = One/term
            sinhkpztm = Half*(term - termi)
            coshkpztm = Half*(term + termi)
            term  = exp(kp*r(3,ia))
            termi = One/term
            sinhkpz = Half*(term - termi)
            coshkpz = Half*(term + termi)
            sumtrigtm(kn,1) = sumtrigtm(kn,1) + az(ia) * &
             (sinkxtm(ialoc,nx) * sinkytm(ialoc,ny) * sinhkpztm - sinkx(ia,nx) * sinky(ia,ny) * sinhkpz)
            sumtrigtm(kn,2) = sumtrigtm(kn,2) + az(ia) * &
             (coskxtm(ialoc,nx) * sinkytm(ialoc,ny) * sinhkpztm - coskx(ia,nx) * sinky(ia,ny) * sinhkpz)
            sumtrigtm(kn,3) = sumtrigtm(kn,3) + az(ia) * &
             (sinkxtm(ialoc,nx) * coskytm(ialoc,ny) * sinhkpztm - sinkx(ia,nx) * cosky(ia,ny) * sinhkpz)
            sumtrigtm(kn,4) = sumtrigtm(kn,4) + az(ia) * &
             (coskxtm(ialoc,nx) * coskytm(ialoc,ny) * sinhkpztm - coskx(ia,nx) * cosky(ia,ny) * sinhkpz)
            sumtrigtm(kn,5) = sumtrigtm(kn,5) + az(ia) * &
             (sinkxtm(ialoc,nx) * sinkytm(ialoc,ny) * coshkpztm - sinkx(ia,nx) * sinky(ia,ny) * coshkpz)
            sumtrigtm(kn,6) = sumtrigtm(kn,6) + az(ia) * &
             (coskxtm(ialoc,nx) * sinkytm(ialoc,ny) * coshkpztm - coskx(ia,nx) * sinky(ia,ny) * coshkpz)
            sumtrigtm(kn,7) = sumtrigtm(kn,7) + az(ia) * &
             (sinkxtm(ialoc,nx) * coskytm(ialoc,ny) * coshkpztm - sinkx(ia,nx) * cosky(ia,ny) * coshkpz)
            sumtrigtm(kn,8) = sumtrigtm(kn,8) + az(ia) * &
             (coskxtm(ialoc,nx) * coskytm(ialoc,ny) * coshkpztm - coskx(ia,nx) * cosky(ia,ny) * coshkpz)
         end do
         termnew = -sumtrigtm(kn,1)**2 - sumtrigtm(kn,2)**2 - sumtrigtm(kn,3)**2 - sumtrigtm(kn,4)**2 &
                   +sumtrigtm(kn,5)**2 + sumtrigtm(kn,6)**2 + sumtrigtm(kn,7)**2 + sumtrigtm(kn,8)**2
         termold = -sumtrig(kn,1)**2   - sumtrig(kn,2)**2   - sumtrig(kn,3)**2   - sumtrig(kn,4)**2 &
                   +sumtrig(kn,5)**2   + sumtrig(kn,6)**2   + sumtrig(kn,7)**2   + sumtrig(kn,8)**2
         term    = kfac2d(kn) * (termnew - termold)
         du%rec   = du%rec + term

      end do
   end do

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUTwoBodyEwaldRec2dlc

!........................................................................

subroutine DUTwoBodyEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none
   character(40), parameter :: txroutine ='DUTwoBodyEwaldRecSPM'

   integer(4) :: m, ia, ialoc, ix, iy, iz, nx, nxa, ny, nya, nz, nza, nxtm, nytm, nztm
   real(8)    :: q, qz, qyz, qtm, qztm, qyztm
   real(8)    :: Qsum
   complex(8) :: vtm(4)

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... initiate

   urecnew = Zero

! ... make generalized charge mesh for trial move

   if (ltime) call CpuAdd('start', 'DUMakeGQMesh', 4, uout)
   QMeshTM = QMesh
   do ialoc = 1, natm
      ia = ianatm(ialoc)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*rtm(m,ialoc), meshmaxtm(m,ialoc), &
            splinetm(:,m,ialoc), derivtm(:,m,ialoc), deriv2tm(:,m,ialoc) )
         meshmaxtm(m,ialoc) = mod(meshmaxtm(m,ialoc)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         qz = az(ia)*spline(iz,3,ia)
         nztm = mod(s(3)+meshmaxtm(3,ialoc)-iz, s(3))
         qztm = az(ia)*splinetm(iz,3,ialoc)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            qyz = qz*spline(iy,2,ia)
            nytm = mod(s(2)+meshmaxtm(2,ialoc)-iy, s(2))
            qyztm = qztm*splinetm(iy,2,ialoc)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               q = qyz*spline(ix,1,ia)
               nxtm = mod(s(1)+meshmaxtm(1,ialoc)-ix, s(1))
               qtm = qyztm*splinetm(ix,1,ialoc)
               QMeshTM(nx+1,ny+1,nz+1) = QMeshTM(nx+1,ny+1,nz+1) - q
               QMeshTM(nxtm+1,nytm+1,nztm+1) = QMeshTM(nxtm+1,nytm+1,nztm+1) + qtm
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'DUMakeGQMesh', 4, uout)

! ... Fourier transformation of the generalized charge distribution

   if (ltime) call CpuAdd('start', 'DUMakeFFT', 4, uout)
   call fftw_execute_dft_r2c( plan_fwd, QmeshTM, FQMesh)
   if (ltime) call CpuAdd('stop', 'DUMakeFFT', 4, uout)

! ... calculate the generalized influence function in the reciprocal space

   if (ltime) call CpuAdd('start', 'DUCalcGIF', 4, uout)
   do nz = 0,s(3)/2
      nza = mod(s(3)-nz,s(3))
      do ny = 0,s(2)/2
         nya = mod(s(3)-ny,s(3))
         do nx = 0,s(1)/2
            nxa = mod(s(3)-nx,s(3))
            vtm = [FQMesh(nx+1,ny+1,nz+1), FQMesh(nx+1,ny+1,nza+1), &
                 FQMesh(nx+1,nya+1,nz+1), FQMesh(nx+1,nya+1,nza+1)]
            Qsum = sum(real(vtm(1:4))**2+aimag(vtm(1:4))**2)
            urecnew = urecnew + Qsum*energyfac(nx,ny,nz)
         end do
      end do
   end do
   du%rec = urecnew - urecold
   if (ltime) call CpuAdd('stop', 'DUCalcGIF', 4, uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

# endif

end subroutine DUTwoBodyEwaldRecSPM

!........................................................................

subroutine DUTwoBodyEwaldSelf     ! no contribution to the energy
end subroutine DUTwoBodyEwaldSelf

!........................................................................

subroutine DUTwoBodyEwaldSurf
   integer(4) :: ia, ialoc
   real(8)    :: fac, term, sumqrx, sumqry, sumqrz, sumqrxt, sumqryt, sumqrzt

   if (.not.lewald2dlc) then

      fac = TwoPi/(Three*vol)
      sumqrx = sum(az(1:na)*r(1,1:na))
      sumqry = sum(az(1:na)*r(2,1:na))
      sumqrz = sum(az(1:na)*r(3,1:na))
      sumqrxt = sumqrx
      sumqryt = sumqry
      sumqrzt = sumqrz
      do ialoc = 1, natm
         ia = ianatm(ialoc)
         sumqrxt = sumqrxt + az(ia)*(rtm(1,ialoc)-r(1,ia))
         sumqryt = sumqryt + az(ia)*(rtm(2,ialoc)-r(2,ia))
         sumqrzt = sumqrzt + az(ia)*(rtm(3,ialoc)-r(3,ia))
      end do
      term = fac*((sumqrxt**2+sumqryt**2+sumqrzt**2) - (sumqrx**2+sumqry**2+sumqrz**2))
      du%rec = du%rec + term

   else

      fac = TwoPi/vol
      sumqrz = sum(az(1:na)*r(3,1:na))
      sumqrzt = sumqrz
      do ialoc = 1, natm
         ia = ianatm(ialoc)
         sumqrzt = sumqrzt + az(ia)*(rtm(3,ialoc)-r(3,ia))
      end do
      term = fac*(sumqrzt**2 - sumqrz**2)
      du%rec = du%rec + term

   end if
end subroutine DUTwoBodyEwaldSurf

!........................................................................

end subroutine DUTwoBodyEwald

!************************************************************************
!> \page denergy denergy.F90
!! **DUWeakChargeEwald**
!! *calculate two-body potential energy difference for titrating systems; k-space*
!************************************************************************


subroutine DUWeakChargeEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUWeakChargeEwald'

   if(ltime) call CpuAdd('start', txroutine, 2, uout)

! ... initiate
   du%rec = Zero

! ... calculate
   if (txewaldrec == 'std') then
       call DUWeakChargeEwaldRecStd
       if (lewald2dlc) call Stop(txroutine,'invalid choice of lewald2dlc .and. lweakcharge',uout)
   else
      call Stop(txroutine,'invalid choice of txewaldrec .and. lweakcharge',uout)
   end if

   call DUWeakChargeEwaldSelf

   if (lsurf .and. master) call DUWeakChargeEwaldSurf

! ... update
   du%rec = EpsiFourPi*du%rec
   du%tot = du%tot + du%rec

   if(ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine DUWeakChargeEwaldRecStd

   character(40), parameter :: txroutine ='DUWeakChargeEwaldRecStd'
   integer(4) :: kn, nx, ny, nz, ia, ialoc, ikvec2
   real(8)    :: term, termnew, termold
   !real(8)    :: kx, ky, kp, sinhkpztm, coshkpztm, sinhkpz, coshkpz, termi

   if(ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   call EwaldSetArrayTM

   kn = kvecoffmyid
   ikvec2 = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         ikvec2 = ikvec2+1
         if (ikvec2 < kvecmyid(1) .or. ikvec2 > kvecmyid(2)) cycle  ! parallelize over k-vectors
         do ialoc = 1, natm
            ia = ianatm(ialoc)
            eikyzm(ia)      = conjg(eiky(ia,ny))     *eikz(ia,nz)
            eikyzp(ia)      =       eiky(ia,ny)      *eikz(ia,nz)
            eikyzmtm(ialoc) = conjg(eikytm(ialoc,ny))*eikztm(ialoc,nz)
            eikyzptm(ialoc) =       eikytm(ialoc,ny) *eikztm(ialoc,nz)
         end do

         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            sumeikrtm(kn,1) = sumeikr(kn,1)
            sumeikrtm(kn,2) = sumeikr(kn,2)
            sumeikrtm(kn,3) = sumeikr(kn,3)
            sumeikrtm(kn,4) = sumeikr(kn,4)
            do ialoc = 1, natm
               ia = ianatm(ialoc)
               sumeikrtm(kn,1) = sumeikrtm(kn,1) &
                                 + aztm(ialoc)*conjg(eikxtm(ialoc,nx))*eikyzmtm(ialoc) &
                                 - az(ia)     *conjg(eikx(ia,nx))*eikyzm(ia)
               sumeikrtm(kn,2) = sumeikrtm(kn,2) &
                                 + aztm(ialoc)*conjg(eikxtm(ialoc,nx))*eikyzptm(ialoc) &
                                 - az(ia)     *conjg(eikx(ia,nx))*eikyzp(ia)
               sumeikrtm(kn,3) = sumeikrtm(kn,3) &
                                 + aztm(ialoc)*      eikxtm(ialoc,nx) *eikyzmtm(ialoc) &
                                 - az(ia)     *      eikx(ia,nx) *eikyzm(ia)
               sumeikrtm(kn,4) = sumeikrtm(kn,4) &
                                 + aztm(ialoc)*      eikxtm(ialoc,nx) *eikyzptm(ialoc) &
                                 - az(ia)     *      eikx(ia,nx) *eikyzp(ia)
            end do
            termnew = real(sumeikrtm(kn,1))**2 + aimag(sumeikrtm(kn,1))**2 + real(sumeikrtm(kn,2))**2 + aimag(sumeikrtm(kn,2))**2 &
                    + real(sumeikrtm(kn,3))**2 + aimag(sumeikrtm(kn,3))**2 + real(sumeikrtm(kn,4))**2 + aimag(sumeikrtm(kn,4))**2
            termold = real(sumeikr(kn,1))**2   + aimag(sumeikr(kn,1))**2   + real(sumeikr(kn,2))**2   + aimag(sumeikr(kn,2))**2 &
                    + real(sumeikr(kn,3))**2   + aimag(sumeikr(kn,3))**2   + real(sumeikr(kn,4))**2   + aimag(sumeikr(kn,4))**2
            term    = kfac(kn)*(termnew - termold)
            du%rec  = du%rec + term
         end do
      end do
   end do

   if(ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine DUWeakChargeEwaldRecStd

!........................................................................

subroutine DUWeakChargeEwaldSelf

   real(8)    :: fac
   integer(4) :: ialoc, ia

! ... atomic contribution

   fac = ualpha/sqrt(Pi)
   do ialoc = myid+1, natm, nproc   ! adapted for _PAR_
      ia = ianatm(ialoc)
      du%rec = du%rec + fac * (az(ia)**2 - aztm(ialoc)**2)
   end do

end subroutine DUWeakChargeEwaldSelf

!........................................................................

subroutine DUWeakChargeEwaldSurf

   integer(4) :: ia, ialoc
   real(8)    :: fac, term, sumqrx, sumqry, sumqrz, sumqrxtm, sumqrytm, sumqrztm

   ! ... if lewald2dlc shall be allowed for titrating systems (weak charges)
   ! ... this needs to be differentiated here. Please see DUTwoBodyEwaldSurf

   fac = TwoPi/(Three*vol)
   sumqrx = sum(az(1:na)*r(1,1:na))
   sumqry = sum(az(1:na)*r(2,1:na))
   sumqrz = sum(az(1:na)*r(3,1:na))
   sumqrxtm = sumqrx
   sumqrytm = sumqry
   sumqrztm = sumqrz
   do ialoc = 1, natm
      ia = ianatm(ialoc)
      sumqrxtm = sumqrxtm + aztm(ialoc)*rtm(1,ialoc) - az(ia)*r(1,ia)
      sumqrytm = sumqrytm + aztm(ialoc)*rtm(2,ialoc) - az(ia)*r(2,ia)
      sumqrztm = sumqrztm + aztm(ialoc)*rtm(3,ialoc) - az(ia)*r(3,ia)
   end do
   term = fac*((sumqrxtm**2+sumqrytm**2+sumqrztm**2) - (sumqrx**2+sumqry**2+sumqrz**2))
   du%rec = du%rec + term

end subroutine DUWeakChargeEwaldSurf

!........................................................................

end subroutine DUWeakChargeEwald

!************************************************************************
!> \page denergy denergy.F90
!! **EwaldSetArrayTM**
!! *calculate eikxtm, eikytm, and eikztm arrays for moving particles*
!************************************************************************


subroutine EwaldSetArrayTM

   use EnergyModule
   implicit none

   integer(4) :: ialoc, icut

   do ialoc = 1, natm
      eikxtm(ialoc,0) = cmplx(One,Zero)
      eikytm(ialoc,0) = cmplx(One,Zero)
      eikztm(ialoc,0) = cmplx(One,Zero)
      eikxtm(ialoc,1) = cmplx(cos(TwoPiBoxi(1)*rtm(1,ialoc)),sin(TwoPiBoxi(1)*rtm(1,ialoc)))
      eikytm(ialoc,1) = cmplx(cos(TwoPiBoxi(2)*rtm(2,ialoc)),sin(TwoPiBoxi(2)*rtm(2,ialoc)))
      eikztm(ialoc,1) = cmplx(cos(TwoPiBoxi(3)*rtm(3,ialoc)),sin(TwoPiBoxi(3)*rtm(3,ialoc)))
   end do
   do icut = 2, ncut
      do ialoc = 1, natm
         eikxtm(ialoc,icut) = eikxtm(ialoc,icut-1)*eikxtm(ialoc,1)
         eikytm(ialoc,icut) = eikytm(ialoc,icut-1)*eikytm(ialoc,1)
         eikztm(ialoc,icut) = eikztm(ialoc,icut-1)*eikztm(ialoc,1)
      end do
   end do

end subroutine EwaldSetArrayTM

!************************************************************************
!> \page denergy denergy.F90
!! **EwaldSetArray2dTM**
!! *calculate sinkxtm, coskxtm, sinkytm, and sinkytm arrays for moving particles*
!************************************************************************


subroutine EwaldSetArray2dTM

   use EnergyModule
   implicit none

   integer(4) :: ialoc, icut

   do ialoc = 1, natm
      sinkxtm(ialoc,0) = One             ! ok!!!
      coskxtm(ialoc,0) = One
      sinkytm(ialoc,0) = One             ! ok!!!
      coskytm(ialoc,0) = One
      sinkxtm(ialoc,1) = sin(TwoPiBoxi(1)*rtm(1,ialoc))
      coskxtm(ialoc,1) = cos(TwoPiBoxi(1)*rtm(1,ialoc))
      sinkytm(ialoc,1) = sin(TwoPiBoxi(2)*rtm(2,ialoc))
      coskytm(ialoc,1) = cos(TwoPiBoxi(2)*rtm(2,ialoc))
   end do
   do icut = 2, ncut2d
      do ialoc = 1, natm
         sinkxtm(ialoc,icut) = sinkxtm(ialoc,icut-1)*coskxtm(ialoc,1) + coskxtm(ialoc,icut-1)*sinkxtm(ialoc,1)
         coskxtm(ialoc,icut) = coskxtm(ialoc,icut-1)*coskxtm(ialoc,1) - sinkxtm(ialoc,icut-1)*sinkxtm(ialoc,1)
         sinkytm(ialoc,icut) = sinkytm(ialoc,icut-1)*coskytm(ialoc,1) + coskytm(ialoc,icut-1)*sinkytm(ialoc,1)
         coskytm(ialoc,icut) = coskytm(ialoc,icut-1)*coskytm(ialoc,1) - sinkytm(ialoc,icut-1)*sinkytm(ialoc,1)
      end do
   end do

end subroutine EwaldSetArray2dTM

!************************************************************************
!> \page denergy denergy.F90
!! **EwaldUpdateArray**
!! *update eikx, eiky, eikz, and sumeikr for moving particles*
!************************************************************************

! ... update sinkx, coskx, sinky, cosky, and termsss ... for moving particles

subroutine EwaldUpdateArray

   use EnergyModule
   implicit none

   integer(4) :: ia, ialoc, icut

   if (txewaldrec == 'std') then
      do icut = 0, ncut
         do ialoc = 1, natm
            ia = ianatm(ialoc)
            eikx(ia,icut) = eikxtm(ialoc,icut)
            eiky(ia,icut) = eikytm(ialoc,icut)
            eikz(ia,icut) = eikztm(ialoc,icut)
         end do
      end do
      sumeikr(1:nkvec,1:4) = sumeikrtm(1:nkvec,1:4)
      if (lewald2dlc) then
         do icut = 0, ncut
            do ialoc = 1, natm
               ia = ianatm(ialoc)
               sinkx(ia,icut) = sinkxtm(ialoc,icut)
               coskx(ia,icut) = coskxtm(ialoc,icut)
               sinky(ia,icut) = sinkytm(ialoc,icut)
               cosky(ia,icut) = coskytm(ialoc,icut)
            end do
         end do
         sumtrig(1:nkvec2d,1:8) = sumtrigtm(1:nkvec2d,1:8)
      end if
# ifdef F03_CBIND
   else if (txewaldrec == 'spm') then
      Qmesh = QmeshTM
      urecold = urecnew
      do ialoc = 1, natm
         ia = ianatm(ialoc)
         spline(:,:,ia) = splinetm(:,:,ialoc)
         meshmax(:,ia) = meshmaxtm(:,ialoc)
         deriv(:,:,ia) = derivtm(:,:,ialoc)
         deriv2(:,:,ia) = deriv2tm(:,:,ialoc)
      end do
# endif
   end if

end subroutine EwaldUpdateArray

!***********************************************************************************************************************

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipole**
!! *calculate potential energy difference from charges and dipoles*
!************************************************************************


subroutine DUDipole

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUDipole'

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   du%stat = Zero                  ! initiate

   call DUDipolePNew               ! calculate new two-body potential energy

   call DUDipolePOld               ! calculate old two-body potential energy

   if (lewald) call DUDipoleEwald   ! calculate reciprocal space energy

   du%stat = EpsiFourPi*du%stat
   du%tot  = du%tot + du%stat      ! update

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

end subroutine DUDipole

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipolePNew**
!! *calculate potential energy from charges and dipoles for new configuration*
!************************************************************************


subroutine DUDipolePNew

   use EnergyModule
   implicit none

   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt
   integer(4) :: ia, ialoc, ialow, iaupp, kialow, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r1, r2, r1i, r2i, r3i, r5i
   real(8)    :: ex, pot, fldx, fldy, fldz, usum
   real(8)    :: dotj, dotjr5i, ErfLocal

   usum = Zero

   kialow = 0
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         ialoc = kialow
         do ia = ialow, iaupp
            ialoc = ialoc + 1
            do ja = jalow, jaupp
               dx = rtm(1,ialoc)-r(1,ja)-dxopbc
               dy = rtm(2,ialoc)-r(2,ja)-dyopbc
               dz = rtm(3,ialoc)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2
               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
               else
                  r3i = r2i*r1i
                  r5i = r2i*r3i
               end if

               dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
               pot =  (az(ja)*r1i + dotj*r3i)

               dotjr5i = (Three*r5i) * dotj
               fldx = (+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx
               fldy = (+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy
               fldz = (+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz

               usum = usum + az(ia)*pot - (diptm(1,ialoc)*fldx + diptm(2,ialoc)*fldy + diptm(3,ialoc)*fldz)

            end do

         end do
      end do
      kialow = kialow+napt(ipt)
   end do

   du%stat = du%stat + usum

end subroutine DUDipolePNew

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipolePOld**
!! *calculate potential energy from charges and dipoles for old configuration*
!************************************************************************


subroutine DUDipolePOld

   use EnergyModule
   implicit none

   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt
   integer(4) :: ia, ialow, iaupp, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, ex, dxopbc, dyopbc, dzopbc, r1, r2, r1i, r2i, r3i, r5i
   real(8)    :: pot, fldx, fldy, fldz, usum
   real(8)    :: dotj, dotjr5i, ErfLocal

   usum = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2
               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
               else
                  r3i = r2i*r1i
                  r5i = r2i*r3i
               end if

               dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
               pot =  (az(ja)*r1i + dotj*r3i)

               dotjr5i = (Three*r5i) * dotj
               fldx = (+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx
               fldy = (+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy
               fldz = (+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz

               usum = usum + az(ia)*pot - (dip(1,ia)*fldx + dip(2,ia)*fldy + dip(3,ia)*fldz)

            end do
         end do
      end do
   end do

   du%stat = du%stat - usum

end subroutine DUDipolePOld

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipoleEwald**
!! *calculate potential energy change; k-space*
!************************************************************************


subroutine DUDipoleEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUDipoleEwald'
   real(8)    :: dusum

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... initiate

   dusum = Zero

! ... calculate

   if (txewaldrec == 'std') then
       call DUDipoleEwaldRecStd
       if (lsurf) call DUDipoleEwaldSurf
   else if (txewaldrec == 'spm') then
       call DUDipoleEwaldRecSPM
   end if
   call DUDipoleEwaldSelf

! ... update

   du%stat = du%stat + dusum

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine DUDipoleEwaldRecStd
   character(40), parameter :: txroutine ='DUDipoleEwaldRecStd'
   integer(4) :: kn, nx, ny, nz, ia, ialoc, ikvec2
   real(8)    :: termx, termy, termz, termxtm, termytm, termztm
   real(8)    :: facmm, facmp, facpm, facpp, facmmtm, facmptm, facpmtm, facpptm
   real(8)    :: kx, ky, kz, termnew, termold

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... calculate eikxtm, eikytm, and eikztm for moving particles

   call EwaldSetArrayTM

   kn = kvecoffmyid
   ikvec2 = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         ikvec2 = ikvec2+1
         if (ikvec2 < kvecmyid(1) .or. ikvec2 > kvecmyid(2)) cycle  ! parallelize over k-vectors
         do ialoc = 1, natm
            ia = ianatm(ialoc)
            eikyzm(ia)      = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia)      =       eiky(ia,ny) *eikz(ia,nz)
            eikyzmtm(ialoc) = conjg(eikytm(ialoc,ny))*eikztm(ialoc,nz)
            eikyzptm(ialoc) =       eikytm(ialoc,ny) *eikztm(ialoc,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            sumeikrtm(kn,1) = sumeikr(kn,1)
            sumeikrtm(kn,2) = sumeikr(kn,2)
            sumeikrtm(kn,3) = sumeikr(kn,3)
            sumeikrtm(kn,4) = sumeikr(kn,4)

            kx = TwoPiBoxi(1)*nx
            ky = TwoPiBoxi(2)*ny
            kz = TwoPiBoxi(3)*nz
            do ialoc = 1, natm
               ia = ianatm(ialoc)
               eikr(ia,1)      = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2)      = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3)      =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4)      =       eikx(ia,nx) *eikyzp(ia)
               eikrtm(ialoc,1) = conjg(eikxtm(ialoc,nx))*eikyzmtm(ialoc)
               eikrtm(ialoc,2) = conjg(eikxtm(ialoc,nx))*eikyzptm(ialoc)
               eikrtm(ialoc,3) =       eikxtm(ialoc,nx) *eikyzmtm(ialoc)
               eikrtm(ialoc,4) =       eikxtm(ialoc,nx) *eikyzptm(ialoc)
               termx = kx*dip(1,ia)
               termy = ky*dip(2,ia)
               termz = kz*dip(3,ia)
               termxtm = kx*diptm(1,ialoc)
               termytm = ky*diptm(2,ialoc)
               termztm = kz*diptm(3,ialoc)
               facmm = -termx - termy + termz
               facmp = -termx + termy + termz
               facpm = +termx - termy + termz
               facpp = +termx + termy + termz
               facmmtm = -termxtm - termytm + termztm
               facmptm = -termxtm + termytm + termztm
               facpmtm = +termxtm - termytm + termztm
               facpptm = +termxtm + termytm + termztm
               sumeikrtm(kn,1) = sumeikrtm(kn,1) + cmplx(az(ia), facmmtm) * eikrtm(ialoc,1) - cmplx(az(ia), facmm) * eikr(ia,1)
               sumeikrtm(kn,2) = sumeikrtm(kn,2) + cmplx(az(ia), facmptm) * eikrtm(ialoc,2) - cmplx(az(ia), facmp) * eikr(ia,2)
               sumeikrtm(kn,3) = sumeikrtm(kn,3) + cmplx(az(ia), facpmtm) * eikrtm(ialoc,3) - cmplx(az(ia), facpm) * eikr(ia,3)
               sumeikrtm(kn,4) = sumeikrtm(kn,4) + cmplx(az(ia), facpptm) * eikrtm(ialoc,4) - cmplx(az(ia), facpp) * eikr(ia,4)
            end do

            termnew = real(sumeikrtm(kn,1))**2 + aimag(sumeikrtm(kn,1))**2 &
                    + real(sumeikrtm(kn,2))**2 + aimag(sumeikrtm(kn,2))**2 &
                    + real(sumeikrtm(kn,3))**2 + aimag(sumeikrtm(kn,3))**2 &
                    + real(sumeikrtm(kn,4))**2 + aimag(sumeikrtm(kn,4))**2

            termold = real(sumeikr(kn,1))**2 + aimag(sumeikr(kn,1))**2 &
                    + real(sumeikr(kn,2))**2 + aimag(sumeikr(kn,2))**2 &
                    + real(sumeikr(kn,3))**2 + aimag(sumeikr(kn,3))**2 &
                    + real(sumeikr(kn,4))**2 + aimag(sumeikr(kn,4))**2

            dusum = dusum + kfac(kn)*(termnew - termold)

         end do
      end do
   end do

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

end subroutine DUDipoleEwaldRecStd

subroutine DUDipoleEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none
   character(40), parameter :: txroutine ='DUDipoleEwaldRecSPM'

   integer(4) :: m, ia, ialoc, ix, iy, iz, nx, nxa, ny, nya, nz, nza, nxtm, nytm, nztm
   real(8)    :: q, q1, dipx, dipx1, qtm, q1tm, dipxtm, dipx1tm
   real(8)    :: dipy, dipy1, dipz, dipz1, dipytm, dipy1tm, dipztm, dipz1tm
   real(8)    :: dipx2, dipy2, dipz2, dipx2tm, dipy2tm, dipz2tm
   real(8)    :: splx, sply, splz, splxtm, splytm, splztm
   real(8)    :: Qsum
   real(8)    :: dsdz
   real(8)    :: dsdztm
   complex(8) :: vtm(4)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   urecnew = Zero

! ... make generalized charge mesh for trial move

   if (ltime) call CpuAdd('start', 'DUMakeGQMesh', 4, uout)
   QMeshTM = QMesh
   do ialoc = 1, natm
      ia = ianatm(ialoc)
      dipx2=dip(1,ia)*dkdr(1)  ! dipole moment in mesh units
      dipy2=dip(2,ia)*dkdr(2)  ! dipole moment in mesh units
      dipz2=dip(3,ia)*dkdr(3)  ! dipole moment in mesh units
      dipx2tm =diptm(1,ialoc)*dkdr(1)  ! dipole moment in mesh units
      dipy2tm =diptm(2,ialoc)*dkdr(2)  ! dipole moment in mesh units
      dipz2tm =diptm(3,ialoc)*dkdr(3)  ! dipole moment in mesh units
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*rtm(m,ialoc), meshmaxtm(m,ialoc), &
         splinetm(:,m,ialoc), derivtm(:,m,ialoc), deriv2tm(:,m,ialoc) )
         meshmaxtm(m,ialoc) = mod(meshmaxtm(m,ialoc)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz = spline(iz,3,ia)
         dsdz = deriv(iz,3,ia)
         nztm = mod(s(3)+meshmaxtm(3,ialoc)-iz, s(3))
         splztm = splinetm(iz,3,ialoc)
         dsdztm = derivtm(iz,3,ialoc)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply = spline(iy,2,ia)
            q1 = az(ia)*sply*splz
            dipx1 = dipx2*sply*splz
            dipy1 = dipy2*deriv(iy,2,ia)*splz
            dipz1 = dipz2*sply*dsdz
            nytm = mod(s(2)+meshmaxtm(2,ialoc)-iy, s(2))
            splytm = splinetm(iy,2,ialoc)
            q1tm = az(ia)*splytm*splztm
            dipx1tm = dipx2tm*splytm*splztm
            dipy1tm = dipy2tm*derivtm(iy,2,ialoc)*splztm
            dipz1tm = dipz2tm*splytm*dsdztm
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx = spline(ix,1,ia)
               q = q1*splx
               dipx = dipx1*deriv(ix,1,ia)
               dipy = dipy1*splx
               dipz = dipz1*splx
               nxtm = mod(s(1)+meshmaxtm(1,ialoc)-ix, s(1))
               splxtm = splinetm(ix,1,ialoc)
               qtm = q1tm*splxtm
               dipxtm = dipx1tm*derivtm(ix,1,ialoc)
               dipytm = dipy1tm*splxtm
               dipztm = dipz1tm*splxtm
               QMeshTM(nx+1,ny+1,nz+1) = &
                    QMeshTM(nx+1,ny+1,nz+1) - (q + dipx + dipy + dipz)
               QMeshTM(nxtm+1,nytm+1,nztm+1) = &
                    QMeshTM(nxtm+1,nytm+1,nztm+1) + (qtm + dipxtm + dipytm + dipztm)
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'DUMakeGQMesh', 4, uout)

! ... Fourier transformation of the generalized charge distribution

   if (ltime) call CpuAdd('start', 'DUMakeFFT', 4, uout)
   call fftw_execute_dft_r2c( plan_fwd , QmeshTM, FQMesh)
   if (ltime) call CpuAdd('stop', 'DUMakeFFT', 4, uout)

! ... calculate the generalized influence function in the reciprocal space

   if (ltime) call CpuAdd('start', 'DUCalcGIF', 4, uout)
   do nz = 0,s(3)/2
      nza = mod(s(3)-nz,s(3))
      do ny = 0,s(2)/2
         nya = mod(s(3)-ny,s(3))
         do nx = 0,s(1)/2
            nxa = mod(s(3)-nx,s(3))
            vtm = [FQMesh(nx+1,ny+1,nz+1), FQMesh(nx+1,ny+1,nza+1), &
                 FQMesh(nx+1,nya+1,nz+1), FQMesh(nx+1,nya+1,nza+1)]
            Qsum = sum(real(vtm(1:4))**2+aimag(vtm(1:4))**2)
            urecnew = urecnew + Qsum*energyfac(nx,ny,nz)
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'DUCalcGIF', 4, uout)
   dusum =  urecnew - urecold

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

# endif

end subroutine DUDipoleEwaldRecSPM

!........................................................................

subroutine DUDipoleEwaldSelf      ! no contribution to the energy
end subroutine DUDipoleEwaldSelf

!........................................................................

subroutine DUDipoleEwaldSurf
   integer(4) :: ia, ialoc
   real(8)    :: fac, termnew, termold
   real(8)    :: sumqrx, sumqry, sumqrz, sumdx, sumdy, sumdz
   real(8)    :: sumqrxt, sumqryt, sumqrzt, sumdxt, sumdyt, sumdzt

   fac = TwoPi/(Three*vol)
   sumqrx = sum(az(1:na)*r(1,1:na))
   sumqry = sum(az(1:na)*r(2,1:na))
   sumqrz = sum(az(1:na)*r(3,1:na))
   sumdx = sum(dip(1,1:na))
   sumdy = sum(dip(2,1:na))
   sumdz = sum(dip(3,1:na))
   sumqrxt = sumqrx
   sumqryt = sumqry
   sumqrzt = sumqrz
   sumdxt = sumdx
   sumdyt = sumdy
   sumdzt = sumdz
   do ialoc = 1, natm
      ia = ianatm(ialoc)
      sumqrxt = sumqrxt + az(ia)*(rtm(1,ialoc)-r(1,ia))
      sumqryt = sumqryt + az(ia)*(rtm(2,ialoc)-r(2,ia))
      sumqrzt = sumqrzt + az(ia)*(rtm(3,ialoc)-r(3,ia))
      sumdxt = sumdxt + diptm(1,ialoc)-dip(1,ia)
      sumdyt = sumdyt + diptm(2,ialoc)-dip(2,ia)
      sumdzt = sumdzt + diptm(3,ialoc)-dip(3,ia)
   end do
   termnew = fac*(sumqrxt**2+sumqryt**2+sumqrzt**2 + sumdxt**2+sumdyt**2+sumdzt**2 &
           + Two*(sumqrxt*sumdxt + sumqryt*sumdyt + sumqrzt*sumdzt))
   termold = fac*(sumqrx**2+sumqry**2+sumqrz**2 + sumdx**2+sumdy**2+sumdz**2 &
           + Two*(sumqrx*sumdx + sumqry*sumdy + sumqrz*sumdz))
   dusum = dusum + termnew - termold
end subroutine DUDipoleEwaldSurf

!........................................................................

end subroutine DUDipoleEwald

!***********************************************************************************************************************

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipoleSph**
!! *calculate charge and dipole potential energy difference*
!************************************************************************


subroutine DUDipoleSph(lhsoverlap)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap
                                                ! =.false. no hard-core overlap
   character(40), parameter :: txroutine ='DUDipoleSph'
   integer(4) :: jp

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   du%twob(0:nptpt) = Zero                      ! initiate
   du%stat = Zero                               ! initiate

   call DUDipoleSphNew(lhsoverlap,jp)              ! calculate new potential energy

#if defined (_PAR_)
   call par_allreduce_logical(lhsoverlap, laux)
#endif

   if (lhsoverlap) goto 400                      ! check hard-core overlap

   call DUDipoleSphOld                              ! calculate old potential energy

   du%stat = EpsiFourPi*du%stat
   du%tot = du%tot + du%twob(0) + du%stat      ! update

400 continue

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

end subroutine DUDipoleSph

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipoleSphNew**
!! *calculate charge and dipole potential energy for new configuration*
!************************************************************************


subroutine DUDipoleSphNew(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='DUDipoleSphNew'

   integer(4) :: ip, iploc, ipt, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r1, r2, r1i, r2i, r3i, r5i, d
   real(8)    :: pot, fldx, fldy, fldz, usumtwob, usumstat
   real(8)    :: dotj, dotjr5i

   lhsoverlap =.true.
   utwobnew(0:nptpt) = Zero
   usumstat = Zero
   usumtwob = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         if (jp == ip) cycle
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         r2 = dx**2+dy**2+dz**2

         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i
         if (r2 < r2atat(iptjpt)) goto 400

         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usumtwob = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                 d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         utwobnew(iptjpt) = utwobnew(iptjpt) + usumtwob

         dotj = dip(1,jp)*dx+dip(2,jp)*dy+dip(3,jp)*dz
         pot =  (az(jp)*r1i + dotj*r3i)

         dotjr5i = (Three*r5i) * dotj
         fldx = (+az(jp)*dx - dip(1,jp))*r3i + dotjr5i*dx
         fldy = (+az(jp)*dy - dip(2,jp))*r3i + dotjr5i*dy
         fldz = (+az(jp)*dz - dip(3,jp))*r3i + dotjr5i*dz

         usumstat = usumstat + az(ip)*pot - (diptm(1,iploc)*fldx + diptm(2,iploc)*fldy + diptm(3,iploc)*fldz)

       end do
    end do

    du%stat = du%stat + usumstat

    if (laimage) then

    call SetImageSph(1,np,2)      ! update image particles to trial configuration

    usumstat = Zero            ! add interaction between moved real particle and all image particles
    do iploc = 1, nptm
       ip = ipnptm(iploc)
       ipt = iptpn(ip)
       do jp = ipmyid(1), ipmyid(2)
          dx = rtm(1,iploc)-rimg(1,jp)
          dy = rtm(2,iploc)-rimg(2,jp)
          dz = rtm(3,iploc)-rimg(3,jp)
          r2 = dx**2+dy**2+dz**2
          r1 = sqrt(r2)
          r1i = One/r1
          r2i = r1i**2
          r3i = r2i*r1i
          r5i = r2i*r3i

          dotj = dipimg(1,jp)*dx+dipimg(2,jp)*dy+dipimg(3,jp)*dz
          pot =  (zimg(jp)*r1i + dotj*r3i)

          dotjr5i = (Three*r5i) * dotj
          fldx = (+zimg(jp)*dx - dipimg(1,jp))*r3i + dotjr5i*dx
          fldy = (+zimg(jp)*dy - dipimg(2,jp))*r3i + dotjr5i*dy
          fldz = (+zimg(jp)*dz - dipimg(3,jp))*r3i + dotjr5i*dz

          usumstat = usumstat + Half*az(ip)*pot - Half*(diptm(1,iploc)*fldx + diptm(2,iploc)*fldy + diptm(3,iploc)*fldz)
       end do
    end do

   du%stat = du%stat + usumstat

   usumstat = Zero    ! add interaction between moved image particle and all real particles
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         if (lptm(jp)) cycle  ! avoid double counting
         dx = rimg(1,ip)-ro(1,jp)
         dy = rimg(2,ip)-ro(2,jp)
         dz = rimg(3,ip)-ro(3,jp)
         r2 = dx**2+dy**2+dz**2
         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i

         dotj = dip(1,jp)*dx+dip(2,jp)*dy+dip(3,jp)*dz
         pot =  (az(jp)*r1i + dotj*r3i)

         dotjr5i = (Three*r5i) * dotj
         fldx = (+az(jp)*dx - dip(1,jp))*r3i + dotjr5i*dx
         fldy = (+az(jp)*dy - dip(2,jp))*r3i + dotjr5i*dy
         fldz = (+az(jp)*dz - dip(3,jp))*r3i + dotjr5i*dz

         usumstat = usumstat + Half*zimg(ip)*pot - Half*(dipimg(1,ip)*fldx + dipimg(2,ip)*fldy + dipimg(3,ip)*fldz)
      end do
   end do

   du%stat = du%stat + usumstat

   end if

   if (lptmdutwob) call Stop(txroutine, 'lptmdutwob = .true.', 6)

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)

   lhsoverlap =.false.

  400 continue

end subroutine DUDipoleSphNew

!************************************************************************
!> \page denergy denergy.F90
!! **DUDipoleSphOld**
!! *calculate charge and dipole potential energy for old configuration*
!************************************************************************


subroutine DUDipoleSphOld

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUDIpSphOld'
   integer(4) :: ip, iploc, ipt, jp, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r1, r2, r1i, r2i, r3i, r5i, d
   real(8)    :: pot, fldx, fldy, fldz, dotj, dotjr5i, usumtwob, usumstat

   utwobold(0:nptpt) = Zero
   usumstat = Zero
   usumtwob = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         if (jp == ip) cycle
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         r2 = dx**2+dy**2+dz**2

         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i

         if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usumtwob = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                     d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         utwobold(iptjpt) = utwobold(iptjpt) + usumtwob

         dotj = dip(1,jp)*dx+dip(2,jp)*dy+dip(3,jp)*dz
         pot =  (az(jp)*r1i + dotj*r3i)

         dotjr5i = (Three*r5i) * dotj
         fldx = (+az(jp)*dx - dip(1,jp))*r3i + dotjr5i*dx
         fldy = (+az(jp)*dy - dip(2,jp))*r3i + dotjr5i*dy
         fldz = (+az(jp)*dz - dip(3,jp))*r3i + dotjr5i*dz

         usumstat = usumstat + az(ip)*pot - (dip(1,ip)*fldx + dip(2,ip)*fldy + dip(3,ip)*fldz)
      end do
   end do

   du%stat = du%stat - usumstat

   if (laimage) then

   call SetImageSph(1,np,3)

   usumstat = Zero            ! add interaction between moved real particle and all image particles
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         dx = ro(1,ip)-rimg(1,jp)
         dy = ro(2,ip)-rimg(2,jp)
         dz = ro(3,ip)-rimg(3,jp)
         r2 = dx**2+dy**2+dz**2
         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i

         dotj = dipimg(1,jp)*dx+dipimg(2,jp)*dy+dipimg(3,jp)*dz
         pot =  (zimg(jp)*r1i + dotj*r3i)

         dotjr5i = (Three*r5i) * dotj
         fldx = (+zimg(jp)*dx - dipimg(1,jp))*r3i + dotjr5i*dx
         fldy = (+zimg(jp)*dy - dipimg(2,jp))*r3i + dotjr5i*dy
         fldz = (+zimg(jp)*dz - dipimg(3,jp))*r3i + dotjr5i*dz

         usumstat = usumstat + Half*az(ip)*pot - Half*(dip(1,ip)*fldx + dip(2,ip)*fldy + dip(3,ip)*fldz)
      end do
   end do

   du%stat = du%stat - usumstat

   usumstat = Zero    ! add interaction between moved image particle and all real particles
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         if (lptm(jp)) cycle  ! Avoid double counting
         dx = rimg(1,ip)-ro(1,jp)
         dy = rimg(2,ip)-ro(2,jp)
         dz = rimg(3,ip)-ro(3,jp)
         r2 = dx**2+dy**2+dz**2
         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i

         dotj = dip(1,jp)*dx+dip(2,jp)*dy+dip(3,jp)*dz
         pot =  (az(jp)*r1i + dotj*r3i)

         dotjr5i = (Three*r5i) * dotj
         fldx = (+az(jp)*dx - dip(1,jp))*r3i + dotjr5i*dx
         fldy = (+az(jp)*dy - dip(2,jp))*r3i + dotjr5i*dy
         fldz = (+az(jp)*dz - dip(3,jp))*r3i + dotjr5i*dz

         usumstat = usumstat + Half*(zimg(ip)*pot - (dipimg(1,ip)*fldx + dipimg(2,ip)*fldy + dipimg(3,ip)*fldz))
      end do
   end do

   du%stat = du%stat - usumstat

   end if

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) call Stop(txroutine, 'lptmdutwob = .true.', 6)

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

end subroutine DUDipoleSphOld

!**********************************************************************************************************************

!************************************************************************
!> \page denergy denergy.F90
!! **DUDielDis**
!! *calculate coulomb energy in a system with dielectric discontinuities*
!************************************************************************


subroutine DUDielDis(lhsoverlap)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap
                                                ! =.false. no hard-core overlap

   character(40), parameter :: txroutine ='DUDielDis'

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   du%twob(0:nptpt) = Zero                      ! initiate
   du%oneb(0:npt) = Zero                        ! initiate

   if (txbc == 'xy') Call DUDielDisPlane(lhsoverlap)
   if (txbc == 'sph') Call DUDielDisSph(lhsoverlap)

   if (lhsoverlap) goto 400                     ! check hard-core overlap

   du%oneb(1:npt) = du%oneb(1:npt)/nproc
   du%twob(0) = sum(du%twob(1:nptpt))
   du%oneb(0)  = sum(du%oneb(1:npt))
   du%tot = du%tot + du%twob(0)                 ! update
   du%tot = du%tot + du%oneb(0)                 ! update

400 continue

    if(ipnptm(1) == 6) then
       write(*,*) ipnptm(1)
       write(*,*) 'lhsoverlap',lhsoverlap
       write(*,*) 'du%tot,du%twob(0),du%oneb(0)',du%tot,du%twob(0),du%oneb(0)
    end if

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

end subroutine DUDielDis

!************************************************************************
!> \page denergy denergy.F90
!! **DUDielDisPlane**
!! *calculate coulomb energy in a system with dielectric discontinuity at z = 0*
!************************************************************************

!     restricted to single particle trial move and serial computation

subroutine DUDielDisPlane(lhsoverlap)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap

   character(40), parameter :: txroutine ='DUDielDisPlane'
   integer(4) :: isign, ipt, jpt, ip, jp, iptjpt, iploc, jploc
   real(8) :: dx, dy, dz, r2, ri, rip, rotemp(3), signEpsi1FourPi, signEpsi2FourPi

   lhsoverlap = .true.
   do isign = 0, 1
      signEpsi1FourPi = (-1)**(isign)*Epsi1FourPi ! assign correct sign
      signEpsi2FourPi = (-1)**(isign)*Epsi2FourPi ! assign correct sign
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         if(isign == 0) then                      ! new configuration
            rotemp(1:3) = rotm(1:3,iploc)
         else                                     ! old configuration
            rotemp(1:3) = ro(1:3,ip)
         end if
         if(abs(rotemp(3)) < radat(ipt)) goto 400 ! particle--dielectric discontinuity overlap
      !   do jp = 1, np
      !      if (ip==jp) cycle
         do jploc = 1, nneighpn(ip)
            jp = jpnlist(jploc,ip)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotemp(1)-ro(1,jp)
            dy = rotemp(2)-ro(2,jp)
            dz = rotemp(3)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < r2atat(iptjpt)) goto 400     ! hs overlap
            ri = one/sqrt(r2)
            if ((rotemp(3) < Zero) .and. (ro(3,jp) < Zero)) then  ! ion--ion and ion--image interaction
               dz = rotemp(3)+ro(3,jp)               ! image location
               call PBCr2(dx,dy,dz,r2)
               rip = one/sqrt(r2)
               du%twob(iptjpt) = du%twob(iptjpt) + signEpsi1FourPi*az(ip)*az(jp)*(ri - delta*rip)
            elseif ((rotemp(3) > Zero) .and. (ro(3,jp) > Zero)) then
               dz = rotemp(3)+ro(3,jp)               ! image location
               call PBCr2(dx,dy,dz,r2)
               rip = one/sqrt(r2)
               du%twob(iptjpt) = du%twob(iptjpt) + signEpsi2FourPi*az(ip)*az(jp)*(ri + delta*rip)
            else
               du%twob(iptjpt) = du%twob(iptjpt) + signEpsi1FourPi*az(ip)*az(jp)*(ri - delta*ri)
            end if
         end do
         if (rotemp(3) < Zero) then                               ! ion--self-image interaction and Born energy
           du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi1FourPi*az(ip)**2)*delta/(two*rotemp(3))
           if(epsi1 < epsi2) du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi1FourPi*az(ip)**2)*(one-eta)/radat(ipt)
         else
           du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi2FourPi*az(ip)**2)*delta/(two*rotemp(3))
           if(epsi1 > epsi2) du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi2FourPi*az(ip)**2)*(one-one/eta)/radat(ipt)
         end if
      end do
      lhsoverlap =.false.
!     call TestDUDielDisPlane(isign)
   end do

400 continue

contains

!........................................................................

subroutine TestDUDielDisPlane(isign)
   integer(4), intent(in) :: isign
   character(49), parameter :: str(0:1) = ['after new configuration        ', 'after new and old configuration']
   call WriteHead(2, trim(txroutine)//'  '//str(isign), uout)
   write(uout,*)
   write(uout,'(a,1f12.3)') 'duTotal    =',sum(du%twob(1:nptpt)) + sum(du%oneb(1:npt))
   write(uout,'(a,6f12.3)') 'duTwoBody  =',du%twob(1:nptpt)
   write(uout,'(a,3f12.3)') 'duOneBody  =',du%oneb(1:npt)
end subroutine TestDUDielDisPlane

end subroutine DUDielDisPlane

!************************************************************************
!> \page denergy denergy.F90
!! **DUDielDisSph**
!! *calculate coulomb energy in a system with a radial dielectric discontinuity*
!************************************************************************

!     restricted to single particle trial move and serial computation

subroutine DUDielDisSph(lhsoverlap)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap        ! =.true. hard-core overlap

   character(40), parameter :: txroutine ='DUDielDisSph'
   integer(4) :: isign, ipt, jpt, ip, jp, iptjpt, iploc, jploc
   real(8) :: r1, r2, r12, fac, cosa, ImageIntSph
   real(8) :: signEpsi1FourPi, signEpsi2FourPi, r1temp(3), r2temp(3)

   lhsoverlap = .true.
   do isign = 0, 1                                ! isign = 0: new configuration; isign = 1; old configuration
      signEpsi1FourPi = (-1)**(isign)*Epsi1FourPi ! assign correct sign
      signEpsi2FourPi = (-1)**(isign)*Epsi2FourPi ! assign correct sign
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         if(isign == 0) then                      ! new configuration
            r1temp(1:3) = rtm(1:3,iploc)
         else                                     ! old configuration
            r1temp(1:3) = r(1:3,ip)
         end if
         r1 = sqrt(r1temp(1)**2+r1temp(2)**2+r1temp(3)**2)
         if (abs(r1-boundaryrad) < radat(iptpn(ip))) goto 400  ! particle--dielectric-discontinuity overlap
 !       do jp = 1, na
 !          if (ip == jp) cycle
         do jploc = 1, nneighpn(ip)
            jp = jpnlist(jploc,ip)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            if(isign == 0) then                  ! new configuration
               if(ip == jp) then
                  r2temp(1:3) = rtm(1:3,iploc)
               else                              ! old confiugration
                  r2temp(1:3) = r(1:3,jp)
!                 if (lptm(jp))  call stop(txroutine,' ltpm: not yet adapted', uout)
               end if
            else
               r2temp(1:3) = r(1:3,jp)
            end if
            r2 = sqrt(r2temp(1)**2+r2temp(2)**2+r2temp(3)**2)
            fac = one/(r1*r2)
            cosa = max(-one, min(fac*sum(r1temp(1:3)*r2temp(1:3)), one))
            r12 = sqrt( (r1temp(1)-r2temp(1))**2 + (r1temp(2)-r2temp(2))**2 + (r1temp(3)-r2temp(3))**2 )
            if (r12**2 < r2atat(iptjpt)) goto 400
                                                 ! ion--ion  and ion--self-image interaction
            if ((r1 < boundaryrad) .and. (r2 < boundaryrad)) then
               du%twob(iptjpt) = du%twob(iptjpt) + signepsi2FourPi*(az(ip)*az(jp)) &
                  * (One/(r12*eta) + ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r2,cosa))
            else
               du%twob(iptjpt) = du%twob(iptjpt) + signepsi2FourPi*(az(ip)*az(jp)) &
                  * (One/r12 + ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r2,cosa))
            end if
         end do
                                                 ! ion--self-image interaction and Born energy
         if (r1 < boundaryrad) then
            du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi2FourPi*az(ip)**2)*ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r1,One)
            if (epsi1 < epsi2) du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi1FourPi*az(ip)**2)*(one-eta)/radat(ipt)
         else
            du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi2FourPi*az(ip)**2)*ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r1,One)
            if (epsi1 > epsi2) du%oneb(ipt) = du%oneb(ipt) + (half*signEpsi2FourPi*az(ip)**2)*(one-one/eta)/radat(ipt)
         end if

         if(lbg) then   ! interaction with the uniform volume charge density with the total charge -nppt(1)*az(1)
            if (r1 < boundaryrad) then                                                  ! inside
               du%oneb(ipt) = du%oneb(ipt) + signepsi1FourPi*ubgfac*az(ip)*((half+eta)*boundaryrad**2-half*r1**2)
            else                                                                        ! outside
               du%oneb(ipt) = du%oneb(ipt) + signepsi1FourPi*ubgfac*az(ip)*eta*boundaryrad**3/r1
            end if
         end if

      end do
      lhsoverlap = .false.
!     call TestDUDielDisSph(isign)
   end do

400 continue

contains

!........................................................................

subroutine TestDUDielDisSph(isign)
   integer(4), intent(in) :: isign
   character(49), parameter :: str(0:1) = ['after new configuration        ', 'after new and old configuration']
   call WriteHead(2, trim(txroutine)//'  '//str(isign), uout)
   write(uout,*)
   write(uout,'(a,1f12.3)') 'duTotal    =',sum(du%twob(1:nptpt)) + sum(du%oneb(1:npt))
   write(uout,'(a,6f12.3)') 'duTwoBody  =',du%twob(1:nptpt)
   write(uout,'(a,3f12.3)') 'duOneBody  =',du%oneb(1:npt)
end subroutine TestDUDielDisSph

!........................................................................

end subroutine DUDielDisSph

!**********************************************************************************************************************

!************************************************************************
!> \page denergy denergy.F90
!! **DUBond**
!! *calculate bond potential energy difference*
!************************************************************************


subroutine DUBond

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUBond'
   integer(4) :: ict, ip, iploc, jp_p, jp_m

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   du%bond = Zero
   do iploc = 1, nptm
      ip = ipnptm(iploc)                           ! moving particle
      ict = ictpn(ip)                              ! chain type
      if (ict == 0) cycle                          ! moving particle not a segment
      jp_m = bondnn(1,ip)                          ! lower neighbour of moving segment
      jp_p = bondnn(2,ip)                          ! upper neighbour of moving segment

      if (jp_m/= 0) then                           ! bond between moving segment and its lower neighbour
         if (.not.lptm(jp_m)) then                 ! lower neighbour is not moved
            call DUBondSub(ro(1:3,jp_m),ro(1:3,jp_m))
         else                                      ! lower neighbour is moved
            if(iploc > 1) then
               call DUBondSub(rotm(1:3,iploc-1),ro(1:3,jp_m))          ! orig
            else
               call DUBondSub(rotm(1:3,iptmpn(jp_m)),ro(1:3,jp_m))     ! Jos
            end if
         end if
      end if
      if (jp_p/= 0) then                           ! bond between moving segment and its upper neighbour
         if (.not.lptm(jp_p)) then                 ! upper neighbour is not moved
            call DUBondSub(ro(1:3,jp_p),ro(1:3,jp_p))
         else                                      ! upper neighbour is moved
            continue                               ! already taken account
         end if
      end if
   end do

   du%bond = du%bond/nproc
   du%tot = du%tot + du%bond

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine DUBondSub(ronew,roold)
   real(8), intent(in) :: ronew(3)       ! new position of neighbour
   real(8), intent(in) :: roold(3)       ! old position of neighbour
   real(8) :: dx, dy, dz, r2new, r2old
   dx = rotm(1,iploc)-ronew(1)
   dy = rotm(2,iploc)-ronew(2)
   dz = rotm(3,iploc)-ronew(3)
   call PBCr2(dx,dy,dz,r2new)
   dx = ro(1,ip)-roold(1)
   dy = ro(2,ip)-roold(2)
   dz = ro(3,ip)-roold(3)
   call PBCr2(dx,dy,dz,r2old)
   du%bond = du%bond + bond(ict)%k*((sqrt(r2new)-bond(ict)%eq)**bond(ict)%p - (sqrt(r2old)-bond(ict)%eq)**bond(ict)%p)
end subroutine DUBondSub

!........................................................................

end subroutine DUBond

!************************************************************************
!> \page denergy denergy.F90
!! **DUAngle**
!! *calculate angle potential energy difference*
!************************************************************************

!     adapted to the case where all displaced chain particles form a single string of particles

subroutine DUAngle

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUAngle'
   integer(4) :: ict, ip, iploc, jp_m, jp_p, jp_mm, jp_pp

   du%angle = Zero
   if (count(angle(1:nct)%k /= Zero) == 0) return

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   do iploc = 1, nptm
      ip = ipnptm(iploc)                         ! moving particle
      ict = ictpn(ip)                            ! chain type
      if (ict == 0) cycle                        ! moving particle not a segment
      jp_m = bondnn(1,ip)                        ! lower neighbour of moving segment
      jp_p = bondnn(2,ip)                        ! upper neighbour of moving segment
      if (jp_m /= 0) jp_mm = bondnn(1,jp_m)      ! next lower neighbour of moving segment
      if (jp_p /= 0) jp_pp = bondnn(2,jp_p)      ! next upper neighbour of moving segment
      if ((jp_m /= 0) .and. (jp_p /= 0)) then    ! segment has two neighbours
         if (lptm(jp_m)) then
            if (lptm(jp_p)) then
               if(( iploc > 1 ) .and. (iploc < nptm)) then
                  call DUAngleSub(rotm(1:3,iploc-1), rotm(1:3,iploc), rotm(1:3,iploc+1), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))
               else
                  call DUAngleSub(rotm(1:3,iptmpn(jp_m)), rotm(1:3,iploc), rotm(1:3,iptmpn(jp_p)), &
                     ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))  ! Jos
               end if
            else
               if(iploc > 1) then
                  call DUAngleSub(rotm(1:3,iploc-1), rotm(1:3,iploc), ro(1:3,jp_p), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))
               else
                  call DUAngleSub(rotm(1:3,iptmpn(jp_m)), rotm(1:3,iploc), ro(1:3,jp_p), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))    ! Jos
               end if
            end if
         else
            if (lptm(jp_p)) then
               if(iploc < nptm) then
                  call DUAngleSub(ro(1:3,jp_m), rotm(1:3,iploc), rotm(1:3,iploc+1), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))
               else
                  call DUAngleSub(ro(1:3,jp_m), rotm(1:3,iploc), rotm(1:3,iptmpn(jp_p)), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p)) ! Jos
               endif
            else
               call DUAngleSub(ro(1:3,jp_m), rotm(1:3,iploc), ro(1:3,jp_p), ro(1:3,jp_m), ro(1:3,ip), ro(1:3,jp_p))
            end if
         end if
      end if
      if (jp_m /= 0) then
         if (jp_mm /= 0) then                     ! segment has two lower neighbours
            if (.not.lptm(jp_m)) then             ! neighbour is not moved
               call DUAngleSub(ro(1:3,jp_mm), ro(1:3,jp_m), rotm(1:3,iploc), ro(1:3,jp_mm), ro(1:3,jp_m), ro(1:3,ip))
            end if
         end if
      end if
      if (jp_p /= 0) then
         if (jp_pp /= 0) then                     ! segment has two upper neigbhours
            if (.not.lptm(jp_p)) then             ! neighbour is not moved
               call DUAngleSub(rotm(1:3,iploc), ro(1:3,jp_p), ro(1:3,jp_pp), ro(1:3,ip), ro(1:3,jp_p), ro(1:3,jp_pp))
            end if
         end if
      end if
   end do

   du%angle = du%angle/nproc
   du%tot = du%tot + du%angle

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine DUAngleSub(ronewm, ronew0, ronewp, rooldm, roold0, rooldp)
   real(8), intent(in) :: ronewm(3)       ! new position of lower segment
   real(8), intent(in) :: ronew0(3)       ! new position of middle segment
   real(8), intent(in) :: ronewp(3)       ! new position of upper segment
   real(8), intent(in) :: rooldm(3)       ! old position of lower segment
   real(8), intent(in) :: roold0(3)       ! old position of middle segment
   real(8), intent(in) :: rooldp(3)       ! old position of upper segment
   real(8) :: dxm, dym, dzm, dxp, dyp, dzp, thetanew, thetaold
   real(8) :: angle_rad

   dxm = ronew0(1)-ronewm(1)
   dym = ronew0(2)-ronewm(2)
   dzm = ronew0(3)-ronewm(3)
   call PBC(dxm,dym,dzm)
   dxp = ronewp(1)-ronew0(1)
   dyp = ronewp(2)-ronew0(2)
   dzp = ronewp(3)-ronew0(3)
   call PBC(dxp,dyp,dzp)
   thetanew = angle_rad(-dxm,-dym,-dzm,dxp,dyp,dzp)

   dxm = roold0(1)-rooldm(1)
   dym = roold0(2)-rooldm(2)
   dzm = roold0(3)-rooldm(3)
   call PBC(dxm,dym,dzm)
   dxp = rooldp(1)-roold0(1)
   dyp = rooldp(2)-roold0(2)
   dzp = rooldp(3)-roold0(3)
   call PBC(dxp,dyp,dzp)
   thetaold = angle_rad(-dxm,-dym,-dzm,dxp,dyp,dzp)

   du%angle = du%angle + angle(ict)%k*((thetanew-angle(ict)%eq)**angle(ict)%p - (thetaold-angle(ict)%eq)**angle(ict)%p)

end subroutine DUAngleSub

!........................................................................

end subroutine DUAngle

!************************************************************************
!> \page denergy denergy.F90
!! **DUCrossLink**
!! *calculate crosslink potential energy difference*
!************************************************************************


subroutine DUCrossLink

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUCrossLink'
   integer(4) :: ip, iploc, jp, jploc, icl

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   du%crosslink = Zero
   do iploc = 1, nptm
      ip = ipnptm(iploc)                                ! moving particle
      do icl = 1, nbondcl(ip)
         jp = bondcl(icl,ip)                            ! crosslinked neighbour
         if (.not.lptm(jp)) then                        ! neighbour is not moved
            call DUCrossLinkSub(ro(1:3,jp),ro(1:3,jp))
         else                                           ! neighbour is moved
            if (ip > jp) cycle                          ! avoid double counting
            jploc = iptmpn(jp)                          ! local id of neighbour
            call DUCrossLinkSub(rotm(1:3,jploc),ro(1:3,jp))
         end if
      end do
   end do

   du%crosslink = du%crosslink/nproc                    ! prepare for allreduce
   du%tot = du%tot + du%crosslink

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine DUCrossLinkSub(ronew,roold)
   real(8), intent(in) :: ronew(3)                   ! new position of neighbour
   real(8), intent(in) :: roold(3)                   ! old position of neighbour
   real(8) :: dx, dy, dz, r2new, r2old
   dx = rotm(1,iploc)-ronew(1)
   dy = rotm(2,iploc)-ronew(2)
   dz = rotm(3,iploc)-ronew(3)
   call PBCr2(dx,dy,dz,r2new)
   dx = ro(1,ip)-roold(1)
   dy = ro(2,ip)-roold(2)
   dz = ro(3,ip)-roold(3)
   call PBCr2(dx,dy,dz,r2old)
   du%crosslink = du%crosslink + clink%k*((sqrt(r2new)-clink%eq)**clink%p - (sqrt(r2old)-clink%eq)**clink%p)
end subroutine DUCrossLinkSub

!........................................................................

end subroutine DUCrossLink

!************************************************************************
!> \page denergy denergy.F90
!! **DUExternal**
!! *calculate external potential difference*
!************************************************************************


subroutine DUExternal(lhepoverlap)

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='DUExternal'
   logical, intent(out) :: lhepoverlap
   integer(4) :: ip, iploc, ipt, ia, iat, ialoc
   real(8)    :: r1, r2

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... initiate

   du%external = Zero
   lhepoverlap = .false.
   ialoc = 0

! ... loop over particles subjected to trial moves and select appropriate energy routine

!      if (istep1 == 16 .and. ipass > 160 .and. ipass < 165) write(*,'(a,10i5)') 'DUExternal: nptm, ipnptm', nptm, ipnptm(1:nptm)
!      if (istep1 == 16 .and. ipass > 160 .and. ipass < 165) write(*,'(a,10f10.5)') 'DUExternal: rotm(3,x)', rotm(3,1:nptm)
!      if (istep1 == 16 .and. ipass > 160 .and. ipass < 165) write(*,'(a,10f10.5)') 'DUExternal: ro(3,x)', ro(3,ipnptm(1:nptm))

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      if (txuext(ipt) == 'wall_z') then
         call DUExternalWallZ
      else if (txuext(ipt) == 'superball_wall_z') then
         call DUExternalSuperballWallZ
      else if (txuext(ipt) == 'gravitation_wall_z') then
         call DUExternalGravitationWallZ
      else if (txuext(ipt) == 'laura') then
         call DUExternalSuperballGravitationEstatFieldWallZ
      else if (txuext(ipt) == 'ramp_wall_z') then
         call DUExternalRampWallZ
      else if (txuext(ipt) == 'sw_wall_zlow') then
         call DUExternalSquareWellZlow
      else if (txuext(ipt) == 'lj_wall_z') then
         call DUExternalLJWallZ
      else if (txuext(ipt) == 'lj_wall_z_ts') then
         call DUExternalLJWallZts
      else if (txuext(ipt) == 'lj_wall_z_mod') then
         call DUExternalLJWallZMod
      else if (txuext(ipt) == 'lj_wall_zlow') then
         call DUExternalLJWallZlow
      else if (txuext(ipt) == 'lj_wall_desorb') then
         call DUExternalLJWallZlowDesorb
      else if (txuext(ipt) == 'estat_field') then
         call DUExternalEstatField
      else if (txuext(ipt) == 'hom_charged_walls') then
         call DUExternalHomChargedWall
      else if (txuext(ipt) == 'i_soft_sphere') then
         call DUExternalISoftSphere
      else if (txuext(ipt) == 'Gunnar_soft_sphere') then
         call DUExternalGunnarSoftSphere
      else if (txuext(ipt) == 'out_hard_ellipsoid') then
         call DUExternalOutHardEllipsoid
      else if (txuext(ipt) == 'capsid_shell') then
         call DUExternalCapsidShell
      else if (txuext(ipt) == 'uniform_shell') then
         call DUExternalUniformShell
      else if (txuext(ipt) == 'sphdielboundary_q') then
         call DUExternalSphDielBoundary_q
         du%external = nproc*du%external
      else if (txuext(ipt)(1:17) == 'sphdielboundary_p') then
         call DUExternalSphDielBoundary_p(txuext(ipt)(18:18))
      else if (txuext(ipt) == 'core_shell') then
         call DUExternalCoreShell
      else if (txuext(ipt) == 'insulating_sphere') then
         call DUExternalInsulatingSphere
      else if (txuext(ipt) == 'hollow_sphere') then
         call DUExternalHollowSphere
      else if (txuext(ipt) == 'orienting_field') then
         call DUExternalOrientingField
      else if (txuext(ipt) == 'hard_cylinder') then
         call DUExternalHardCylinder
      else if (txuext(ipt) == 'lekkerkerker-tuinier') then
         call DUExternalLTDepletion
      end if
   end do

   if (lhepoverlap) goto 400

   du%external = du%external/nproc
   du%tot = du%tot + du%external

400 continue

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine DUExternalWallZ                       ! wall at abs(z) = wall_z_ext
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      if (abs(rtm(3,ialoc)) > wall_z_ext) then
         lhepoverlap = .true.
         exit
      end if
   end do
end subroutine DUExternalWallZ

!........................................................................

subroutine DUExternalGravitationWallZ            ! graviatation force and wall at abs(z) = wall_z_ext
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)
      if (abs(rtm(3,ialoc)) > wall_z_ext) then
         lhepoverlap = .true.
         exit
      end if
      du%external = du%external - gravitation_force*(rtm(3,ialoc) - r(3,ia))
   end do
end subroutine DUExternalGravitationWallZ

!........................................................................

subroutine DUExternalSuperballWallZ              ! superball and wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='DUExternalSuperballWallZ'
   real(8) :: dz, oriim(3,3)
   logical :: SuperballOverlap
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)
      dz = two*(rtm(3,ialoc) - sign(wall_z_ext,r(3,ia))) ! distance between atom and its image across the closest z-wall
      oriim(1:3,1:2) = oritm(1:3,1:2,ialoc)
      oriim(1:3,3) = -oritm(1:3,3,ialoc)
      if (SuperballOverlap(dz**2,[zero,zero,dz],oritm(1,1,ialoc),oriim(1,1))) then
         lhepoverlap = .true.
!         write(*,'(a)') 'overlap'
!         write(*,'(a,3i8           )') 'ia, ip', ia, ip
!         write(*,'(a,3f8.3,2x,3f8.3)') 'r, dz', rtm(1:3,ialoc), dz
!         write(*,'(a,9f8.3,2x,9f8.3)') 'ori, oriim', oritm(1:3,1:3,ialoc), oriim(1:3,1:3)
         exit
      end if
   end do
end subroutine DUExternalSuperballWallZ

!........................................................................

subroutine DUExternalSuperballGravitationEstatFieldWallZ              ! superball, gravitation, external field, and wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='DUExternalSuperballGravitationWallZ'
   real(8) ::  dz, oriim(3,3)
   logical :: SuperballOverlap
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)
      dz = two*(rtm(3,ialoc) - sign(wall_z_ext,r(3,ia))) ! distance between atom and its image across the closest z-wall
      oriim(1:3,1:2) = oritm(1:3,1:2,ialoc)
      oriim(1:3,3) = -oritm(1:3,3,ialoc)
      if (SuperballOverlap(dz**2,[zero,zero,dz],oritm(1,1,ialoc),oriim(1,1))) then
         lhepoverlap = .true.
!         write(*,'(a)') 'overlap'
!         write(*,'(a,3i8           )') 'ia, ip', ia, ip
!         write(*,'(a,3f8.3,2x,3f8.3)') 'r, dz', rtm(1:3,ialoc), dz
!         write(*,'(a,9f8.3,2x,9f8.3)') 'ori, oriim', oritm(1:3,1:3,ialoc), oriim(1:3,1:3)
         exit
      end if
      du%external = du%external - gravitation_force*(rtm(3,ialoc) - r(3,ia))
      du%external = du%external + az(ia)*sum((rtm(1:3,ialoc)-r(1:3,ia))*efield_ext(1:3))
      du%external = du%external - sum((diptm(1:3,ialoc)-dip(1:3,ia)) * efield_ext(1:3))
   end do
end subroutine DUExternalSuperballGravitationEstatFieldWallZ

!........................................................................

subroutine DUExternalRampWallZ                     ! walls at z = +-wall_z_ext
   real(8) :: wall_hs, dz, width, slope, uterm
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)

      wall_hs = wall_z_ext - radat(iat)            ! hs interaction at z = wall_hs
      if (abs(rtm(3,ialoc)) > wall_hs) then         ! hs overlap with z-wall
         lhepoverlap = .true.
         exit
      end if
      width = (lambda_ramp_ext-One)*Two*radat(iat)
      dz = abs(rtm(3,ialoc)) - (wall_hs-width)
      if ((dz > Zero) .and. (dz <= width)) then     ! on the slope
         slope = epsilon_ramp_ext/width
         uterm = -slope*dz
         du%external = du%external + uterm
      end if

      dz = abs(r(3,ia)) - (wall_hs-width)
      if ((dz > Zero) .and. (dz <= width)) then     ! on the slope
         slope = epsilon_ramp_ext/width
         uterm = -slope*dz
         du%external = du%external - uterm
      end if

   end do
end subroutine DUExternalRampWallZ

!........................................................................

subroutine DUExternalSquareWellZlow                 ! wall at z = -wall_z-ext
   real(8) :: wall_hs, width, dz, uterm
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)
      wall_hs = wall_z_ext - radat(iat)            ! hs interaction at z = wall_hs
      if (abs(rtm(3,ialoc)) > wall_hs) then
         lhepoverlap = .true.
         exit
      end if
      width = (lambda_sw_ext-One)*Two*radat(iat)
      dz = rtm(3,ialoc) + wall_hs - width          ! dz < 0 => in the well
      uterm = -epsilon_sw_ext*(Half - (One/Pi)*atan(dz/alpha_sw_ext))
      du%external = du%external + uterm
      dz = r(3,ia) + wall_hs - width               ! dz < 0 => in the well
      uterm = -epsilon_sw_ext*(Half - (One/Pi)*atan(dz/alpha_sw_ext))
      du%external = du%external - uterm
   end do
end subroutine DUExternalSquareWellZlow

!........................................................................

subroutine DUExternalLJWallZ
   real(8) :: zi3
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      if (abs(rtm(3,ialoc)) > wall_z_ext) then                  ! across one of the walls
         lhepoverlap = .true.
         exit
      else
         iat = iatan(ia)
         zi3 = One/(wall_z_ext-abs(rtm(3,ialoc)))**3           ! lj-wall at +-wall_z_ext
         du%external = du%external + (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3)
         zi3 = One/(wall_z_ext-abs(r(3,ia)))**3                ! lj-wall at +-wall_z_ext
         du%external = du%external - (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3)
      end if
   end do
end subroutine DUExternalLJWallZ

!........................................................................

subroutine DUExternalLJWallZts
   real(8) :: zi3
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      if (abs(rtm(3,ialoc)) > wall_z_ext) then                  ! across one of the walls
         lhepoverlap = .true.
         exit
      else
         iat = iatan(ia)
         if (wall_z_ext-abs(rtm(3,ialoc)) < zmin_ext(iat)) then
            zi3 = One/(wall_z_ext-abs(rtm(3,ialoc)))**3        ! lj-wall at +-wall_z_ext
            du%external = du%external + (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3 + delta_ext(iat))
         end if
         if (wall_z_ext-abs(r(3,ia)) < zmin_ext(iat)) then
            zi3 = One/(wall_z_ext-abs(r(3,ia)))**3             ! lj-wall at +-wall_z_ext
            du%external = du%external - (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3 + delta_ext(iat))
         end if
      end if
   end do
end subroutine DUExternalLJWallZts

!........................................................................

subroutine DUExternalLJWallZMod
   real(8) :: z1, zi3, ulj, hx, hy, Onephxhy, utemp
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      if (abs(rtm(3,ialoc)) > wall_z_ext) then                  ! across one of the walls
         lhepoverlap = .true.
         exit
      else
         iat = iatan(ia)
         z1 = wall_z_ext-abs(rtm(3,ialoc))           ! lj-wall at +-wall_z_ext
         zi3 = (One/z1)**3
         ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3
         hx = cos(TwoPilxi_ext*rtm(1,ialoc))
         hy = cos(TwoPilyi_ext*rtm(2,ialoc))
         Onephxhy = One + hx*hy
         if (z1 < zmin_ext(iat) ) then
            utemp = ulj+(delta_ext(iat)*c_ext)*Onephxhy
         else
            utemp = ulj*(One-c_ext*Onephxhy)
         end if
         du%external = du%external + utemp

         z1 = wall_z_ext-abs(r(3,ia))                ! lj-wall at +-wall_z_ext
         zi3 = (One/z1)**3
         ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3
         hx = cos(TwoPilxi_ext*r(1,ia))
         hy = cos(TwoPilyi_ext*r(2,ia))
         Onephxhy = One + hx*hy
         if (z1 < zmin_ext(iat) ) then
            utemp = ulj+(delta_ext(iat)*c_ext)*Onephxhy
         else
            utemp = ulj*(One-c_ext*Onephxhy)
         end if
         du%external = du%external - utemp

      end if
   end do
end subroutine DUExternalLJWallZMod

!........................................................................

subroutine DUExternalLJWallZlow
   real(8) :: zi3
      do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
         ialoc = ialoc+1
         if (abs(rtm(3,ialoc)) > wall_z_ext) then      ! across one of the walls
            lhepoverlap = .true.
            exit
         else
            iat = iatan(ia)
            zi3 = One/(rtm(3,ialoc)-(-wall_z_ext))**3 ! lj-wall at -wall_z_ext
            du%external = du%external + (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3)
            zi3 = One/(r(3,ia)-(-wall_z_ext))**3 ! lj-wall at -wall_z_ext
            du%external = du%external - (z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3)
         end if
      end do
end subroutine DUExternalLJWallZlow

!........................................................................

subroutine DUExternalLJWallZlowDesorb                    ! Niklas 2006-12-20 desorb
   character(40), parameter :: txroutine ='DUExternalLJWallZlowDesorb'
   call Stop(txroutine,' mc not valied for this external potential', uout)
end subroutine DUExternalLJWallZlowDesorb

!........................................................................

subroutine DUExternalEstatField
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      du%external = du%external + az(ia)*sum((rtm(1:3,ialoc)-r(1:3,ia))*efield_ext(1:3))
      du%external = du%external - sum((diptm(1:3,ialoc)-dip(1:3,ia)) * efield_ext(1:3))
   end do
end subroutine DUExternalEstatField

!........................................................................

! calculate wall and long-range contribution of electrostatic energy
! (eq. 16 of Jonsson et al, JPC, 84, 2179, 1980)

subroutine DUExternalHomChargedWall
   real(8) :: scd, a, a2, b, z, zsb, zsb2, longrangecontr
   integer(4) :: i
   scd = surfchargeden/(ech*1.d20)
   a = boxlen2(1)
   a2 = a**2
   b = boxlen2(3)

   do ia=ianpn(ip),ianpn(ip)+napt(ipt)-1
      ialoc=ialoc+1
      iat = iatan(ia)

!obsolete    if (abs(rtm(3,ialoc)) > boxlen2(3)-radat(iat)) then
!               lhepoverlap = .true.
!               exit
!            end if

      z = rtm(3,ialoc)    ! trial configuration
      do i = -1,1,2
         zsb = abs(i*b-z)
         zsb2 = zsb**2
         du%external = du%external + EpsiFourPi*zat(iat)*scd*(8.d0*a*log( (sqrt(Two*a2+zsb2)+a)/(sqrt(a2+zsb2)))- &
         Two*zsb*(asin( (a2**2-zsb2**2-Two*a2*zsb2)/(a2+zsb2)**2 )+Half*Pi))
      end do
      if (llongrangecontr) then
         du%external = du%external + EpsiFourPi*zat(iat)*longrangecontr(boxlen2(1), z, scd, mninchden, zdist, chden)
      endif

      z = r(3,ia)         ! old configuration
      do i = -1,1,2
         zsb = abs(i*b-z)
         zsb2 = zsb**2
         du%external = du%external - EpsiFourPi*zat(iat)*scd*(8.d0*a*log( (sqrt(Two*a2+zsb2)+a)/(sqrt(a2+zsb2)))- &
         Two*zsb*(asin( (a2**2-zsb2**2-Two*a2*zsb2)/(a2+zsb2)**2 )+Half*Pi))
      end do
      if (llongrangecontr) then
         du%external = du%external - EpsiFourPi*zat(iat)*longrangecontr(boxlen2(1), z, scd, mninchden, zdist, chden)
       endif

    end do

end subroutine DUExternalHomChargedWall

!........................................................................

subroutine DUExternalISoftSphere
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      r1 = sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      if (r1 > ruext(1,ipt)) du%external = du%external + auext*(r1-ruext(1,ipt))**nuext
      r1 = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if (r1 > ruext(1,ipt)) du%external = du%external - auext*(r1-ruext(1,ipt))**nuext
   end do
end subroutine DUExternalISoftSphere

!........................................................................

subroutine DUExternalGunnarSoftSphere
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      r1 = sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      du%external = du%external + auext*(r1/ruext(1,ipt))**nuext
      r1 = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      du%external = du%external - auext*(r1/ruext(1,ipt))**nuext
   end do
end subroutine DUExternalGunnarSoftSphere

!........................................................................

subroutine DUExternalOutHardEllipsoid
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc=ialoc+1
      r2 = (rtm(1,ialoc)*ruexti(1,ipt))**2+(rtm(2,ialoc)*ruexti(2,ipt))**2+(rtm(3,ialoc)*ruexti(3,ipt))**2
      if (r2 < One) then
         lhepoverlap = .true.
         exit
      end if
   end do
end subroutine DUExternalOutHardEllipsoid

!........................................................................

subroutine DUExternalCapsidShell
   do ia=ianpn(ip),ianpn(ip)+napt(ipt)-1
      ialoc=ialoc+1
      r1=sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      if (r1 >= (rcap-radat(ipt)) .and. r1 <= (rcap+dcap+radat(ipt))) then
         lhepoverlap = .true.
         exit
      end if
   end do
end subroutine DUExternalCapsidShell

!........................................................................

subroutine DUExternalUniformShell
   real(8), save :: rcapchage = 2.0d0   ! radial location of capside charge from inner surface
   do ia=ianpn(ip),ianpn(ip)+napt(ipt)-1
      ialoc=ialoc+1
      r1=sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      if (r1 <= (rcap-radat(ipt))) then
         du%external = du%external + zat(ipt)*zat(1)*EpsiFourPi/(rcap+rcapchage)
      else if (r1 >= (rcap-radat(ipt)) .and. r1 <= (rcap+dcap+radat(ipt))) then
         lhepoverlap = .true.
         exit
      else if (r1 > (rcap+dcap+radat(ipt))) then
         du%external = du%external + zat(ipt)*zat(1)*EpsiFourPi/r1
      end if
      r1=sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if (r1 <= (rcap-radat(ipt))) then
         du%external = du%external - zat(ipt)*zat(1)*EpsiFourPi/(rcap+rcapchage)
      else if (r1 >= (rcap-radat(ipt)) .and. r1 <= (rcap+dcap+radat(ipt))) then
         lhepoverlap = .true.
         exit
      else if (r1> (rcap+dcap+radat(ipt))) then
         du%external = du%external - zat(ipt)*zat(1)*EpsiFourPi/r1
      end if
   end do
end subroutine DUExternalUniformShell

!........................................................................

subroutine DUExternalSphDielBoundary_q
   real(8) :: rfacnew, rrnew, thetanew, phinew, rrationew
   real(8) :: rfacold, rrold, thetaold, phiold, rratioold
   real(8) :: sum, term
   complex(8) xCCLM
   integer(4) :: l, m
   real(8) :: sumnew, sumold
   real(8) :: termnew, termold

   if (slave) return  ! generation and use of QQ are made on master

! ... calculate trial multipole moments

   QQtm(0:lmaxdiel,0:lmaxdiel) = QQ(0:lmaxdiel,0:lmaxdiel)
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1

      rrationew = boundaryrad/sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      rratioold = boundaryrad/sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      call CarToSph('rad',rtm(1,ialoc),rtm(2,ialoc),rtm(3,ialoc),rrnew,thetanew,phinew)
      call CarToSph('rad',r(1,ia),r(2,ia),r(3,ia),rrold,thetaold,phiold)
      rfacnew = One
      rfacold = One
      do l = 0, lmaxdiel
         rfacnew = rfacnew*rrationew
         rfacold = rfacold*rratioold
         do m = 0, l
            QQtm(m,l) = QQtm(m,l) + az(ia)*(rfacnew*xCCLM(l,m,thetanew,phinew,0)-rfacold*xCCLM(l,m,thetaold,phiold,0))
         end do
      end do
   end do

! ... calculate external energy change

   sumold = zero
   sumnew = zero
   do l = 0, lmaxdiel
      sum = zero
      do m = 0,l
         term = mfac(m)*((real(QQtm(m,l))**2+aimag(QQtm(m,l))**2) - (real(QQ(m,l))**2+aimag(QQ(m,l))**2))
         termold = mfac(m)*((real(QQ(m,l))**2+aimag(QQ(m,l))**2))
         termnew = mfac(m)*((real(QQtm(m,l))**2+aimag(QQtm(m,l))**2))
         sum = sum + term
         sumold = sumold + lfac(l)*termold
         sumnew = sumnew + lfac(l)*termnew
      if (itest == 98 .and. l == 1) write(*,'(a,2f14.8)') 'DUExternaĺSphDielr Boundary_q: QQ(m,1)', QQ(m,l)
      end do
      du%external = du%external + lfac(l)*sum
   end do
   if (itest == 98) write(*,'(a,f14.8)') 'DUExternaĺSphDielr Boundary_q: du%external', du%external
   if (itest == 98) write(*,'(a,f14.8)') 'DUExternaĺSphDielr Boundary_q: u_old_ext', sumold
   if (itest == 98) write(*,'(a,f14.8)') 'DUExternaĺSphDielr Boundary_q: u_new_ext', sumnew

end subroutine DUExternalSphDielBoundary_q

!........................................................................

subroutine DUExternalSphDielBoundary_p(str)
   character(40), parameter :: txroutine ='DUExternalSphDielBoundary_p'

   character(1), intent(in) :: str           ! select equations to be used
                                             ! '0': short expansion
                                             ! '1': long expansion
   real(8) :: PL, sum, r2(1:3), ufac
   real(8) :: r1new, r2new, cosanew, tnew, facnew, rfacnew, term1new, term2new, term3new
   real(8) :: r1old, r2old, cosaold, told, facold, rfacold, term1old, term2old, term3old
   integer(4) :: ia, ialoc, ja, l

! ... calculate external energy

   sum = zero
   ialoc = 0
   do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
      ialoc = ialoc+1
      r1new = sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      if (r1new < boundaryrad) call Stop(txroutine, 'r1new < boundaryrad', uout)
      r1old = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      do ja = 1, np
         r2old = sqrt(r(1,ja)**2+r(2,ja)**2+r(3,ja)**2)
         if (r2old**2 < boundaryrad**2) cycle ! exclue any atoms inside the dielectric sphere
         if (ia == ja) then                   ! tentative
            r2new = r1new
            r2(1:3) = rtm(1:3,ialoc)
            ufac = one
         else
            if (lptm(ja)) then
              ! r2new = r1new  ?           ! not yet finalized
              ! ufac = ?
               call stop(txroutine,' ltpm', uout)
            else
               r2new = r2old
               r2(1:3) = r(1:3,ja)
               ufac = two
            end if
         end if

         facnew = one/(r1new*r2new)
         facold = one/(r1old*r2old)
         cosanew = (rtm(1,ialoc)*r2(1)+rtm(2,ialoc)*r2(2)+rtm(3,ialoc)*r2(3))*facnew
         cosanew = max(-one, min(cosanew, one))
         cosaold = (r(1,ia)*r(1,ja)+r(2,ia)*r(2,ja)+r(3,ia)*r(3,ja))*facold
         cosaold = max(-one, min(cosaold, one))
         tnew = boundaryrad**2*facnew
         told = boundaryrad**2*facold
         facnew = sqrt(one-two*cosanew*tnew+tnew**2)
         facold = sqrt(one-two*cosaold*told+told**2)

         if (str == '1') then
            term1new = lfac1*tnew*(one/facnew - one)
            term1old = lfac1*told*(one/facold - one)
            term2new = lfac2*tnew*log(half*(facnew+one-cosanew*tnew))
            term2old = lfac2*told*log(half*(facold+one-cosaold*told))
            term3new = zero
            term3old = zero
         else if (str == '2') then
            term1new = lfac1*tnew*(one/facnew - one)
            term1old = lfac1*told*(one/facold - one)
            if (ja /= ia) then
               term2new = lfac2*log((one-cosanew)/(facnew+tnew-cosanew))
               term2old = lfac2*log((one-cosaold)/(facold+told-cosaold))
            else
               term2new = lfac2*log(one-tnew)
               term2old = lfac2*log(one-told)
            end if
            term3new = lfac2*tnew
            term3old = lfac2*told
         end if

         rfacnew = tnew
         rfacold = told
         do l = 1, lmaxdiel
            rfacnew = rfacnew*tnew
            rfacold = rfacold*told
            term3new = term3new - lfac3(l)*rfacnew*PL(l,cosanew)
            term3old = term3old - lfac3(l)*rfacold*PL(l,cosaold)
         end do
         if (itest == 98) then
            write(*,*)
            write(*,'(a)') 'DUExternalSphDielr Boundary_p'
            write(*,'(a,3f14.6)') 'r1old,    r1new   ', r1old, r1new
            write(*,'(a,3f14.6)') 'r2old,    r2new   ', r2old, r2new
            write(*,'(a,3f14.6)') 'cosaold,  cosanew ', cosaold, cosanew
            write(*,'(a,3f14.6)') 'facold,   facnew  ', facold, facnew
            write(*,'(a,3f14.6)') 'term1old, term1new', term1old, term1new
            write(*,'(a,3f14.6)') 'term2old, term2new', term2old, term2new
            write(*,'(a,3f14.6)') 'term3old, term3new', term3old, term3new
            write(*,'(a,3f14.6)') 'termnew', term1new+ term2new + term3new
            write(*,'(a,3f14.6)') 'termold', term1old+ term2old + term3old
            write(*,'(a,3f14.6)') 'd(term)', (term1new + term2new + term3new - term1old - term2old - term3old)
         end if
         sum = sum + ufac*(az(ia)*az(ja))*(term1new + term2new + term3new - term1old - term2old - term3old)
      end do
   end do
   du%external = du%external + sum
   if (itest == 98) write(*,'(a,g15.5)') 'du%external',du%external

end subroutine DUExternalSphDielBoundary_p

!........................................................................

subroutine  DUExternalCoreShell
   do ia = ipnan(ip), ipnan(ip)+napt(ipt)-1
      ialoc = ialoc+1
      r1 = sqrt(rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      if ((r1 < rChargeIn) .or. (r1 > rChargeOut)) then
         lhepoverlap = .true.
         exit
      end if
   end do
end subroutine DUExternalCoreShell

!........................................................................

subroutine DUExternalInsulatingSphere
   real(8):: sum
   sum = Zero
   do ia = ipnan(ip), ipnan(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)
      r2=rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2
      r1=sqrt(r2)
      if (r1 < rInsSphere) then
         sum = sum + zat(iat)/(two*rInsSphere)*(three-r2/rInsSphere**2)
      else
         sum = sum + zat(iat)/r1
      end if
      r2=r(1,ia)**2+r(2,ia)**2+r(3,ia)**2
      r1=sqrt(r2)
      if (r1 < rInsSphere) then
         sum = sum - zat(iat)/(two*rInsSphere)*(three-r2/rInsSphere**2)
      else
         sum = sum - zat(iat)/r1
      end if
   end do
   du%external = du%external + (zInsSphere*EpsiFourPi)*sum
end subroutine DUExternalInsulatingSphere

!........................................................................

subroutine DUExternalHollowSphere
   real(8):: sum
   sum = Zero
   do ia = ipnan(ip), ipnan(ip)+napt(ipt)-1
      ialoc = ialoc+1
      iat = iatan(ia)

      r2=rtm(1,ialoc)**2+rtm(2,ialoc)**2+rtm(3,ialoc)**2
      r1=sqrt(r2)
      if (r1 < rInSphere) then
         sum = sum - zat(iat)*threehalf*(rInSphere**2-rOutSphere**2)/(rOutSphere**3-rInSphere**3)
      else if (r1 <= rOutSphere) then
         sum = sum - zat(iat)/(rOutSphere**3-rInSphere**3) * &
            (half*r2+rInSphere**3/r1-threehalf*rInSphere**2+threehalf*(rInSphere**2-rOutSphere**2))
      else if (r1 > rOutSphere) then
         sum = sum + zat(iat)/r1
      end if

      r2=r(1,ia)**2+r(2,ia)**2+r(3,ia)**2
      r1=sqrt(r2)
      if (r1 < rInSphere) then
         sum = sum + zat(iat)*threehalf*(rInSphere**2-rOutSphere**2)/(rOutSphere**3-rInSphere**3)
      else if (r1 <= rOutSphere) then
         sum = sum + zat(iat)/(rOutSphere**3-rInSphere**3) * &
            (half*r2+rInSphere**3/r1-threehalf*rInSphere**2+threehalf*(rInSphere**2-rOutSphere**2))
      else if (r1 > rOutSphere) then
         sum = sum - zat(iat)/r1
      end if
   end do
   du%external = du%external + (zInsSphere*EpsiFourPi)*sum
end subroutine DUExternalHollowSphere

!........................................................................

subroutine DUExternalOrientingField
      du%external = du%external - ofstrength* &
      (ofaxis(1)*(oritm(3,1,iploc)-ori(3,1,ip)) &
     + ofaxis(2)*(oritm(3,2,iploc)-ori(3,2,ip)) &
     + ofaxis(3)*(oritm(3,3,iploc)-ori(3,3,ip)))
end subroutine DUExternalOrientingField

!........................................................................

subroutine  DUExternalHardCylinder
   real(8) :: r2new, r2old, rlnew, rlold
   do ia = ipnan(ip), ipnan(ip)+napt(ipt)-1
      ialoc = ialoc+1
      r2new = rtm(1,ialoc)**2+rtm(2,ialoc)**2
      if (r2new < rCylinder**2) then
         lhepoverlap = .true.
         exit
      end if
      if (zCylinder /= zero) then
         r2old = r(1,ia)**2+r(2,ia)**2
         rlnew = sqrt(r2new + cyllen2**2)
         rlold = sqrt(r2old + cyllen2**2)
         du%external = du%external + EpsiFourPi*zCylinder*zat(iatan(ia)) * &
            log((rlnew+cyllen2)/(rlnew-cyllen2)*(rlold-cyllen2)/(rlold+cyllen2))
      end if
   end do
end subroutine DUExternalHardCylinder

!........................................................................

subroutine DUExternalLTDepletion
   character(40), parameter :: txroutine ='DUExternalLTDepletion'
   real(8) :: h, LinInter
   integer(4) :: ih
   do ia = ipnan(ip), ipnan(ip)+napt(ipt)-1
      ialoc = ialoc+1
      h = boxlen2(3) - abs(rtm(3,ialoc))        ! select nearest surface
      if (h < Zero) then
         lhepoverlap = .true.
         exit
      end if
      ih = int(h*dhexti)                        ! get bin
      if (ih < nhext) du%external = du%external + LinInter(hext(ih), hext(ih+1), h, uext(ih), uext(ih+1))
      h = boxlen2(3) - abs(r(3,ia))             ! select nearest surface
      ih = int(h*dhexti)                        ! get bin
      if (ih < nhext) du%external = du%external - LinInter(hext(ih), hext(ih+1), h, uext(ih), uext(ih+1))
   end do
end subroutine DUExternalLTDepletion

!........................................................................

end subroutine DUExternal

! ........................... Under Development ...............................

!************************************************************************
!> \page denergy denergy.F90
!! **UWeakChargeANew**
!! *calculate two-body potential energy for new configuration*
!************************************************************************

!     only monoatomic particles
!     weak charges, limited to charged hard spheres

subroutine UWeakChargeANew(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='UWeakChargeANew'

   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

if (itest == 90) then
   call writehead(3,txroutine, uout)                           !cc
   write(uout,*) 'nptm, ipnptm(1:nptm)', nptm, ipnptm(1:nptm)  !cc
end if

   utwobnew(0:nptpt) = Zero
   lhsoverlap =.true.

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
         end if
         if (r2 > rcut2) cycle
         if (r2 < r2atat(iptjpt)) goto 400

         if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end

         if (.not.laztm(iploc)) cycle  ! iploc uncharged
         if (.not.laz(jp)) cycle       ! jp uncharged
         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                           d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         if (itest == 90) write(uout,'(a,2i5,2l5,f10.5)') ' ip, jp, laztm, laz, usum', ip,jp,laztm(iploc),laz(jp), usum  !cc

         utwobnew(iptjpt) = utwobnew(iptjpt) + usum
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then
      do iploc = myid+1, nptm, nproc ! adapted for _PAR_
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-rotm(1,jploc)
            dy = rotm(2,iploc)-rotm(2,jploc)
            dz = rotm(3,iploc)-rotm(3,jploc)
            call PBCr2(dx,dy,dz,r2)
            if (lellipsoid) Then
                if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc))) goto 400
            end if
            if (r2 > rcut2) cycle
            if (r2 < r2atat(iptjpt)) goto 400


            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end

            if (.not.laztm(iploc)) cycle  ! iploc uncharged
            if (.not.laztm(jploc)) cycle  ! jploc uncharged
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

          if (itest == 90) write(uout,'(a,2i5,2l5,f10.5)') ' ip ,jp, laztm, laztm, usum', ip,jp,laztm(iploc),laztm(jploc), usum  !cc

            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         end do
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)
   lhsoverlap =.false.

  400 continue

      if (itest == 90) write(*,'(a,10f10.5)') '   utwobnew(0:nptpt)', utwobnew(0:nptpt) !cc

end subroutine UWeakChargeANew

!************************************************************************
!> \page denergy denergy.F90
!! **UWeakChargeAOld**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     only monoatomic particles
!     weak charges, limited to charged hard spheres

subroutine UWeakChargeAOld

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UWeakChargeAOld'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz, r2, d, usum
   logical    :: EllipsoidOverlap, SuperballOverlap

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

if (itest == 90) then
   call writehead(3,txroutine, uout)                           !cc
   write(uout,*) 'nptm, ipnptm(1:nptm)', nptm, ipnptm(1:nptm)  !cc
end if

   utwobold(0:nptpt) = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 > rcut2) cycle
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
         end if

         if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end

         if (.not.laz(ip)) cycle  ! ip uncharged
         if (.not.laz(jp)) cycle  ! jp uncharged
         ibuf = iubuflow(iptjpt)
         do
            if (r2 >= ubuf(ibuf)) exit
            ibuf = ibuf+12
         end do
         d = r2-ubuf(ibuf)
         usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                           d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         if (itest == 90) write(uout,'(a,2i5,2l5,f10.5)') ' ip, jp, laz, laz, usum', ip,jp,laz(ip),laz(jp), usum  !cc

         utwobold(iptjpt) = utwobold(iptjpt) + usum
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then
      do iploc = myid+1, nptm, nproc ! adapted for _PAR_
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
              if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) goto 400
            end if
            if (r2 < r2atat(iptjpt)) goto 400
            if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end

            ibuf = iubuflow(iptjpt)
            if (.not.laz(ip)) cycle  ! ip uncharged
            if (.not.laz(jp)) cycle  ! jp uncharged
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

          if (itest == 90)  write(uout,'(a,2i5,2l5,f10.5)') ' both moved: ip, jp, usum', ip,jp,laz(ip),laz(jp), usum  !cc

            utwobold(iptjpt) = utwobold(iptjpt) + usum
         end do
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

     if (itest == 90)  write(*,'(a,10f10.5)') '   utwobold(0:nptpt)', utwobold(0:nptpt) !cc

end subroutine UWeakChargeAOld

!************************************************************************
!> \page denergy denergy.F90
!! **UWeakChargePNew**
!! *calculate two-body potential energy for new configuration*
!************************************************************************

!     weak charges, limited to charged hard spheres

subroutine UWeakChargePNew(lhsoverlap,jp)

   use EnergyModule
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='UWeakChargePNew'
   logical, parameter :: lintrapartint = .true.                   ! enable intraparticle interaction
   integer(4) :: ip, iploc, ipt, jploc, jpt, iptjpt, ibuf
   integer(4) :: ia, ialoc, ialow, iaupp, kialow, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum
   integer(4) :: jaloc, kjalow
   logical    :: EllipsoidOverlap, SuperballOverlap

if (itest == 90) then
   call writehead(3,txroutine, uout)                           !cc
   write(uout,'(a,6i5)') 'nptm, ipnapm(iploc)', nptm, (ipnptm(1:nptm))
   write(uout,'(a,6i5)') 'natm, ianatm(ialoc)', natm, (ianatm(1:natm))
   write(uout,'(a,2l5)') 'laz(ialoc)', (laz(ianatm(1:natm)))
   write(uout,'(a,2l5)') 'laztm(ialoc)', (laztm(1:natm))
end if

   lhsoverlap =.true.
   utwobnew(0:nptpt) = Zero

   kialow = 0
   do iploc = 1, nptm
      ip = ipnptm(iploc)
         !   write(uout,*) 'ip', ip
      ipt = iptpn(ip)
      if (nptm == natm) then           ! adaptation for charge change move involving one atom per particle
         ialow = ianatm(iploc)
         iaupp = ialow
      else
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
      end if
      do jploc = 0, nneighpn(ip)
         if (jploc == 0) then
            if (.not.lintrapartint) cycle
            !do not calculate the intrapartenergy for a non charge change move
            if (nptm /= natm) cycle

            !if nptm == natm and the move was not a ChargeChangeMove then the
            !number of atoms for the moved particles is one. Therefore the
            !interaction for ip == jp is not calculated. Explanation:
            !ialow = ianatm(iploc) = ianpn(ip) (as 1 Atom per particle, see also SetTrialAtomProp
            !as iaupp = ialow: ia = ianpn(ip) (see loop)
            !jalow = ianpn(jp) = ianpn(ip) (see below)
            !as jaupp = jalow: ja = ianpn(ip) (see loop)
            !therefore ja = ia and no intraparticle energies are calculated
            !as this case is always skipped

            jp = ip
         else
            jp = jpnlist(jploc,ip)
         end if
         if (lptm(jp) .and. jp /= ip) cycle
       !  write(uout,*) ' jp', jp
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = rotm(1,iploc)-ro(1,jp)
         dy = rotm(2,iploc)-ro(2,jp)
         dz = rotm(3,iploc)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle
         if (lellipsoid) Then
            if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         end if
         if (lsuperball) Then
            if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
         end if

         usum = Zero
         ialoc = kialow
         do ia = ialow, iaupp
            iat = iatan(ia)
            ialoc = ialoc+1
        !    write(uout,*) '  ia, ialoc, laztm(ialoc)', ia, laztm(ialoc), ialoc
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
            if (jp == ip .and. (nptm /= natm)) jalow = ia + 1
            do ja = jalow, jaupp
               if (ja == ia) cycle
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = rtm(1,ialoc)-r(1,ja)-dxopbc
               dy = rtm(2,ialoc)-r(2,ja)-dyopbc
               dz = rtm(3,ialoc)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               if (jp /= ip .and. r2 < r2atat(iatjat)) goto 400
               if (.not.laztm(ialoc)) cycle  ! ia uncharged
               if (.not.laz(ja)) cycle       ! ja uncharged
               ! do not check for overlap within one particle

               if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                       d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

           if (itest == 90) write(uout,'(a,2i5,2l5,2f10.5)') '   ia, ja, laztm, laz, usum', ia,ja,laztm(ialoc),laz(ja), usum  !cc

            end do
         end do
         utwobnew(iptjpt) = utwobnew(iptjpt) + usum
      end do
      if (nptm == natm) then
         kialow = kialow+1
      else
         kialow = kialow+napt(ipt)
      end if
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
  !   write(uout,*) 'HERE lptmdutwob'
      kialow = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
      if (nptm == natm) then           ! adaptation for charge change move involving one atom per particle
         ialow = ianatm(iploc)
         iaupp = ialow
      else
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
      end if
   !         write(uout,*) 'ip, ialow, iaupp', ip, ialow, iaupp
      if (nptm == natm) then
         kialow = kialow+1
      else
         kialow = kialow+napt(ipt)
      end if
         kjalow = kialow
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = rotm(1,iploc)-rotm(1,jploc)
            dy = rotm(2,iploc)-rotm(2,jploc)
            dz = rotm(3,iploc)-rotm(3,jploc)
            call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
            dx = dx-dxopbc
            dy = dy-dyopbc
            dz = dz-dzopbc
            r2 = dx**2+dy**2+dz**2
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
            end if

            usum = Zero
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
      if (nptm == natm) then           ! adaptation for charge change move involving one atom per particle
         jalow = ianatm(jploc)
         jaupp = jalow
      else
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
      end if
    !        write(uout,*) 'jp, jalow, jaupp', jp, jalow, jaupp
            ialoc = 0
      if (nptm == natm) then
         kjalow = kialow+1
      else
         kjalow = kialow+napt(ipt)
      end if
            do ia = ialow, iaupp
               iat = iatan(ia)
               ialoc = ialoc+1
           !       write(uout,*) 'ialoc, laztm(ialoc)', ialoc, laztm(ialoc)
               if (.not.laztm(ialoc)) cycle  ! iploc uncharged   !xnewx
               jaloc = 1  ! temorary fix
               do ja = jalow, jaupp
                  jaloc = jaloc+1
           !       write(uout,*) 'jaloc, laztm(jaloc)', jaloc, laztm(jaloc)
                  if (.not.laztm(jaloc)) cycle  ! jploc uncharged  !xnewx
                  jat = iatan(ja)
                  iatjat = iatat(iat,jat)
                  dx = rtm(1,ialoc)-rtm(1,jaloc)-dxopbc
                  dy = rtm(2,ialoc)-rtm(2,jaloc)-dyopbc
                  dz = rtm(3,ialoc)-rtm(3,jaloc)-dzopbc
                  r2 = dx**2+dy**2+dz**2
                  if (r2 < r2atat(iatjat)) goto 400

                  if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
                  ibuf = iubuflow(iatjat)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                          d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

           if (itest == 90) write(uout,'(a,2i5,2l5,2f10.5)') 'xx   ia, ja, laztm, laztm, usum', &
              ia, ja, laztm(ialoc), laztm(jaloc), usum
               end do
            end do
      if (nptm == natm) then
         kjalow = kjalow+1
      else
         kjalow = kjalow+napt(jpt)
      end if
            utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         end do
      if (nptm == natm) then
         kialow = kialow+1
      else
         kialow = kialow+napt(ipt)
      end if
      end do
   end if

   utwobnew(0) = sum(utwobnew(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)

   lhsoverlap =.false.

  400 continue

   if (itest == 90)   write(uout,'(a,10f10.5)') '   utwobnew(0:nptpt)', utwobnew(0:nptpt) !cc

end subroutine UWeakChargePNew

!************************************************************************
!> \page denergy denergy.F90
!! **UWeakChargePOld**
!! *calculate two-body potential energy for old configuration*
!************************************************************************

!     weak charges, limited to charged hard spheres

subroutine UWeakChargePOld

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UWeakChargePOld'
   logical, parameter :: lintrapartint = .true.                   ! enable intraparticle interaction
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf
   integer(4) :: ia, ialoc, ialow, iaupp, kialow, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum
   integer(4) :: jaloc, kjalow
   logical    :: EllipsoidOverlap, SuperballOverlap

if (itest == 90) then
   call writehead(3,txroutine, uout)                           !cc
   write(uout,*) 'nptm, ipnptm(1:nptm)', nptm, ipnptm(1:nptm)  !cc
end if

   utwobold(0:nptpt) = Zero

   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      if (nptm == natm) then           ! adaptation for charge change move involving one atom per particle
         ialow = ianatm(iploc)
         iaupp = ialow
      end if
      do jploc = 0, nneighpn(ip)
         if (jploc == 0) then
            if (.not.lintrapartint) cycle
            if (nptm /= natm) cycle
            jp = ip
         else
            jp = jpnlist(jploc,ip)
         end if
         if (lptm(jp) .and. jp /= ip) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         usum = Zero
         do ia = ialow, iaupp
            if (.not.laz(ia)) cycle  ! ia uncharged
            iat = iatan(ia)
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
            if (jp == ip .and. (nptm /= natm)) jalow = ia + 1
            do ja = jalow, jaupp
               if (ja == ia) cycle
               if (.not.laz(ja)) cycle  ! ja uncharged
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2

               if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                       d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

           if (itest == 90) write(uout,'(a,2i5,2l5,2f10.5)') '   ia, ja, laz, laz, usum', ia,ja,laz(ia),laz(ja), usum !cc

            end do
         end do
         utwobold(iptjpt) = utwobold(iptjpt) + usum
      end do
   end do

! ... contribution from pairs where both particles are displaced

   if (lptmdutwob) then      ! not adapted for _PAR_ !!
      kialow = 0
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
         kjalow = kialow+napt(ipt)
         do jploc = iploc+1, nptm
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
            dx = dx-dxopbc
            dy = dy-dyopbc
            dz = dz-dzopbc
            r2 = dx**2+dy**2+dz**2
            if (r2 > rcut2) cycle
            if (lellipsoid) Then
               if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
            end if
            if (lsuperball) Then
               if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
            end if

            usum = Zero
            jalow = ianpn(jp)
            jaupp = jalow+napt(jpt)-1
            ialoc = kialow
            do ia = ialow, iaupp
               if (.not.laz(ia)) cycle  ! ip uncharged              !xnewx
               iat = iatan(ia)
               ialoc = ialoc+1
               jaloc = kjalow
               do ja = jalow, jaupp
                  if (.not.laz(ja)) cycle  ! jp uncharged           !xnewx
                  jaloc = jaloc+1
                  jat = iatan(ja)
                  iatjat = iatat(iat,jat)
                  dx = r(1,ia)-r(1,ja)-dxopbc
                  dy = r(2,ia)-r(2,ja)-dyopbc
                  dz = r(3,ia)-r(3,ja)-dzopbc
                  r2 = dx**2+dy**2+dz**2
                  if (r2 < r2atat(iatjat)) goto 400

                  if (r2 < r2umin(iatjat)) goto 400       ! outside lower end
                  ibuf = iubuflow(iatjat)
                  do
                     if (r2 >= ubuf(ibuf)) exit
                     ibuf = ibuf+12
                  end do
                  d = r2-ubuf(ibuf)
                  usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                          d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

               end do
               kjalow = kjalow+napt(jpt)
            end do
            utwobold(iptjpt) = utwobold(iptjpt) + usum
         end do
         kialow = kialow+napt(ipt)
      end do
   end if

   utwobold(0) = sum(utwobold(1:nptpt))
   du%twob(0:nptpt) = du%twob(0:nptpt) - utwobold(0:nptpt)

  400 continue

      if (itest == 90) write(uout,'(a,10f10.5)') '   utwobold(0:nptpt)', utwobold(0:nptpt)  !cc

end subroutine UWeakChargePOld
