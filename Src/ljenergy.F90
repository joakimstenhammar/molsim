module FlexLJEnergyModule

!calculate Potential of Form 4*eps*[(sigma/r)**12 - (sigma/r)**6]

implicit none
private
public IOPotFlexLJ

   real(8)  , allocatable  :: 4eps(:,:)
   real(8)  , allocatable  :: sig2(:,:)

subroutine IOPotFlexLJ(iStage)

   use MolModule, only: nptpt, npt, iptpt
   use MolModule, only: uin, uout
   use MolModule, only: iWriteInput, iReadInput
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='IOPotFlexLJ'
   character(80), parameter :: txheading ='flexible Lennard-Jones potential data'
   real(8), allocatable, save :: epsilonLJ(:), sigmaLJ(:)
   integer  :: ipt, jpt, ptpt
   integer  :: io

   namelist /nmlLennardJones/ epsilonLJ, sigmaLJ

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      !allocate memory for raeding the input variables
      if (.not.allocated(epsilonLJ)) then
         allocate(epsilonLJ(nptpt))
         allocate(sigmaLJ(nptpt))
      end if
      epsilonLJ = 0.0d0
      sigmaLJ = 0.0d0

      !reading the input variables
      rewind(uin)
      read(uin,nmlPotential, iostat=io)
      if(io .ne. 0) then
            !print Warning when the namelist is not given
            call Warn(txroutine ,'nmlLennardJones is not given in input file. Using default values.', uout)
      end if

   case(iWriteInput)

      !allocate and calculate the variables of the LJ Potential
      if (.not.allocated(eps)) then
         allocate(eps(npt,npt))
         allocate(sig6(npt,npt))
      end if

      do ipt = 1, npt
         do jpt = 1, npt
            4eps(ipt,jpt) = 4*epsilonLJ(iptpt(ipt,jpt)) ! only 4*epsilon is needed
            sig2(ipt,jpt) = sigmaLJ(iptpt(ipt,jpt))**2 ! only sigma**2 is needed
         end do
      end do

      if (master) then

         call WriteHead(2, txheading, uout)
         write(uout,'()')
         write(uout,'(a)') "Using flexible LJ Potential"
         write(uout,'(a)') "Note that no other interactions can be used."
         write(uout,'()')
         write(uout,'(a)') 'LJ parameters'
         write(uout,'(a54)') repeat('-',54)
         write(uout,'(a20,2(tr2,a15))') 'Type-Type','Epsilon', 'Sigma'
         write(uout,'(a20,2(tr2,a15))') repeat('-',20), repeat('-',15), repeat('-',15), repeat('-',15)
         do ptpt = 1, nptpt
            write(uout,'(a20,2*(tr2,ES15.8))') txptpt(ptpt), epsilonLJ(ptpt), sigmaLJ(ptpt)
         end do
         write(uout,'(a54)') repeat('-',54)
      end if
   end select

end subroutine IOPotFlexLJ


!subroutine ULJNewA(lhsoverlap,jp)

   !!use EnergyModule, only:
   !!use MolModule, only: uout
   !!use MolModule, only: lmonoatom
   !implicit none

   !logical,    intent(out) :: lhsoverlap
   !integer(4), intent(out) :: jp

   !real(8)  :: dr(3), r2

!!  todo (pascal): have it checked in potential routine
   !!if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)
   !!if (any(zat(1:npt)).ne. 0.0d0) txpot=='(1,6,12)'

!!   write(uout,*) txroutine

   !utwobnew(0:nptpt) = Zero
   !lhsoverlap =.true.

   !do iploc = 1, nptm
      !ip = ipnptm(iploc)
      !ipt = iptpn(ip)
!!      write(uout,'(a,i5,3f10.5)') 'ip,rotm(1:3,iploc)',ip, rotm(1:3,iploc)

      !do jploc = 1, nneighpn(ip) !using standard neighbor list
         !jp = jpnlist(jploc,ip)
         !if (lptm(jp)) cycle ! do not calculate interactions where both particles moves
         !jpt = iptpn(jp)
         !iptjpt = iptpt(ipt,jpt)
         !dr(1:3) = rotm(1:3,iploc) - ro(1:3,jp)
         !call PBCr2(dr(1), dr(2), dr(3) ,r2)
         !if (lellipsoid) Then
            !if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp),radellipsoid2,aellipsoid)) goto 400
         !end if
         !if (lsuperball) Then
            !if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),ori(1,1,jp))) goto 400
         !end if
         !if (r2 > rcut2) cycle
         !if (r2 < r2atat(iptjpt)) goto 400

         !if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
         !ibuf = iubuflow(iptjpt)
         !do
            !if (r2 >= ubuf(ibuf)) exit
            !ibuf = ibuf+12
         !end do
         !d = r2-ubuf(ibuf)
         !usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                           !d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

         !utwobnew(iptjpt) = utwobnew(iptjpt) + usum
!!         write(uout,'(a,i5,6f10.5)') 'jp,ro(1:3,jp), r2, usum, utwobnew(iptjpt)',jp, ro(1:3,jp), r2, usum, utwobnew(iptjpt)
      !end do
   !end do

!! ... contribution from pairs where both particles are displaced

   !if (lptmdutwob) then      ! not adapted for _PAR_ !!
      !do iploc = 1, nptm
         !ip = ipnptm(iploc)
         !ipt = iptpn(ip)
         !do jploc = iploc+1, nptm
            !jp = ipnptm(jploc)
            !jpt = iptpn(jp)
            !iptjpt = iptpt(ipt,jpt)
            !dx = rotm(1,iploc)-rotm(1,jploc)
            !dy = rotm(2,iploc)-rotm(2,jploc)
            !dz = rotm(3,iploc)-rotm(3,jploc)
            !call PBCr2(dx,dy,dz,r2)
            !if (lellipsoid) Then
                !if (EllipsoidOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc),radellipsoid2,aellipsoid)) goto 400
            !end if
            !if (lsuperball) Then
               !if (SuperballOverlap(r2,[dx,dy,dz],oritm(1,1,iploc),oritm(1,1,jploc))) goto 400
            !end if
            !if (r2 > rcut2) cycle
            !if (r2 < r2atat(iptjpt)) goto 400

            !if (r2 < r2umin(iptjpt)) goto 400       ! outside lower end
            !ibuf = iubuflow(iptjpt)
            !do
               !if (r2 >= ubuf(ibuf)) exit
               !ibuf = ibuf+12
            !end do
            !d = r2-ubuf(ibuf)
            !usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                              !d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

            !utwobnew(iptjpt) = utwobnew(iptjpt) + usum
         !end do
      !end do
   !end if

   !utwobnew(0) = sum(utwobnew(1:nptpt))
   !du%twob(0:nptpt) = du%twob(0:nptpt) + utwobnew(0:nptpt)
   !lhsoverlap =.false.

  !400 continue

!end subroutine ULJNew
!subroutine UTwoBodyAOld

end module flexLJEnergyModule
