module DirectDEnergyModule

implicit none
private
public UTwoBodyNewDirect

subroutine ULJNew(lhsoverlap,jp)

   use EnergyModule, only:
   use MolModule, only: uout
   use MolModule, only: lmonoatom
   implicit none

   logical,    intent(out) :: lhsoverlap
   integer(4), intent(out) :: jp

   character(40), parameter :: txroutine ='UTwoBodyANew'


   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)
   if (any(zat(1:npt)).ne. 0.0d0) txpot=='(1,6,12)'

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

end subroutine ULJNew
subroutine UTwoBodyAOld

end module DirectDEnergyModule
