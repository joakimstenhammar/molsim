module FlexLJEnergyModule

!calculate Potential of Form 4*eps*[(sigma/r)**12 - (sigma/r)**6]

implicit none
private
public IOPotFlexLJ
public UFlexLJ

   real(8)  , allocatable  :: 4eps(:,:)
   real(8)  , allocatable  :: sig2(:,:)

contains

pure function uLJ(r2, ipt, jpt) result(energy)
   real(8), intent(in)  :: r2
   integer(4), intent(in)  :: ipt
   integer(4), intent(in)  :: jpt
   real(8)  :: energy
   real(8)  :: invr6

   invr6 = (sig2(ipt,jpt)/r2)**3
   energy = 4eps(ipt,jpt)*(invr6**2 - invr6)

end function uLJ

pure function fLJ(r2, ipt, jpt) result(force)
   real(8), intent(in)  :: r2
   integer(4), intent(in)  :: ipt
   integer(4), intent(in)  :: jpt
   real(8)  :: energy
   real(8)  :: invr6, invr1

   invr6 = (sig2(ipt,jpt)/r2)**3
   invr1 = -6.0d0/sqrt(r2)
   force = 4eps(ipt,jpt)*(Two*invr6**2-invr6)*invr1

end function fLJ

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

   namelist /nmlFlexLJ/ epsilonLJ, sigmaLJ

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
      read(uin,nmlFlexLJ, iostat=io)
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

subroutine UFlexLJ(utwob, utot, force, virial)
   use MolModule, only lclist
   use MolModule, only lmonoatom

   real(8), allocatable, intent(inout)  :: utwob(:)
   real(8), intent(inout)               :: utot
   real(8), allocatable, intent(inout)  :: force(:,:)
   real(8), intent(inout)               :: virial

   !if(lclist) then
      !if(lmonoatom) then
         !call UFlexLJCellMono(utwob, utot, force, virial)
      !else
         !call UFlexLJCellPoly(utwob, utot, force, virial)
      !end if
   !else
      !if(lmonoatom) then
         call UFlexLJMono(utwob, utot, force, virial)
      !else
         !call UFlexLJPoly(utwob, utot, force, virial)
      !end if
   !end if
end subroutine

subroutine UFlexLJMono(utwob, utot, force, virial)

   use MolModule, only: np, iptpn, iptpt, nptpt
   use MolModule, only: myid
   use MolModule, only: ro
   use MolModule, only: jpnlist
   real(8), allocatable, intent(inout)  :: utwob(:)
   real(8), intent(inout)               :: utot
   real(8), allocatable, intent(inout)  :: force(:,:)
   real(8), intent(inout)               :: virial

   integer(4)  :: ip, jp, iptjpt, iploc, jploc
   real(8)  :: dr(3), r2
   real(8)  :: uloc, floc, virtwob

   virtwob         = 0.0d0
   do iploc = 1, Getnpmyid() !get particles from nlist
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (lmc) then
            if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)

         dr = ro(1:3,ip)-ro(1:3,jp)
         call PBCr2(dr(1), dr(2), dr(3),r2)
         if (r2 > rcut2) cycle
         if (r2 < r2umin(iptjpt)) call StopUTwoBodyA
         if (r2 < r2atat(iptjpt)) then
            uloc = 1d10                ! emulate hs overlap
         else
            uloc = uLJ(r2, ipt, jpt)
            floc = fLJ(r2, ipt, jpt)
         end if

         utwob(iptjpt) = utwob(iptjpt) + uloc
         force(1:3,ip) = force(1:3,ip) + floc*dr(1:3)
         force(1:3,jp) = force(1:3,jp) + floc*dr(1:3)
         virial     = virial     - (floc * r2)
      end do
   end do

   utwob(0) = sum(utwob(1:nptpt))
   utot     = utot + utwob(0)

end subroutine UFlexLJMono



end module flexLJEnergyModule
