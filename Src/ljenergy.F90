module FlexLJEnergyModule

!calculate Potential of Form 4*eps*[(sigma/r)**12 - (sigma/r)**6]

implicit none
private
public IOPotFlexLJ
public UFlexLJ, DUFlexLJ

   real(8)  , allocatable  :: foureps(:)
   real(8)  , allocatable  :: sig2(:)

contains

subroutine IOPotFlexLJ(iStage)

   use MolModule, only: natat, txatat
   use MolModule, only: uin, uout, ltrace, master
   use MolModule, only: iWriteInput, iReadInput
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='IOPotFlexLJ'
   character(80), parameter :: txheading ='flexible Lennard-Jones potential data'
   real(8), allocatable, save :: epsilonLJ(:), sigmaLJ(:)
   integer  :: iatjat
   integer  :: io

   namelist /nmlFlexLJ/ epsilonLJ, sigmaLJ

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      !allocate memory for raeding the input variables
      if (.not.allocated(epsilonLJ)) then
         allocate(epsilonLJ(natat))
         allocate(sigmaLJ(natat))
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
      if (.not.allocated(foureps)) then
         allocate(foureps(natat))
         allocate(sig2(natat))
      end if

      foureps = 4.0d0*epsilonLJ
      sig2 = sigmaLJ**2

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
         do iatjat = 1, natat
            write(uout,'(a20,2(tr2,ES15.8))') txatat(iatjat), epsilonLJ(iatjat), sigmaLJ(iatjat)
         end do
         write(uout,'(a54)') repeat('-',54)
      end if
   end select

end subroutine IOPotFlexLJ

pure function uLJ(r2, iatjat) result(energy)
   implicit none
   real(8), intent(in)  :: r2
   integer(4), intent(in)  :: iatjat
   real(8)  :: energy
   real(8)  :: invr6

   invr6 = (sig2(iatjat)/r2)**3
   energy = foureps(iatjat)*(invr6**2 - invr6)

end function uLJ

pure function fLJ(r2, iatjat) result(force) !note that this subroutine is not tuned for efficiency
   implicit none
   real(8), intent(in)  :: r2
   integer(4), intent(in)  :: iatjat
   real(8)  :: force
   real(8)  :: invr6, invr1

   invr6 = (sig2(iatjat)/r2)**3
   invr1 = -6.0d0/sqrt(r2)
   force = foureps(iatjat)*(2.0d0*invr6**2-invr6)*invr1

end function fLJ

subroutine ufvLJdr(dr, iatjat, energy, forcea, forceb, virial, loverlap) !calculated the force, virial and energy from a distance and atom-type combination
   use MolModule, only: rcut2, r2umin, r2atat
   implicit none

   real(8), intent(in)  :: dr(3) !note that the PBC are neede to be already applied!
   integer(4), intent(in)  :: iatjat
   real(8), intent(inout)  :: energy
   real(8), intent(inout)  :: forcea(3)
   real(8), intent(inout)  :: forceb(3)
   real(8), intent(inout)  :: virial
   logical, intent(out)  :: loverlap
   real(8)  :: totforce, r2

   loverlap = .false.
   r2 = sum(dr**2)
   if (r2 > rcut2) return
   if ((r2 < r2umin(iatjat)) .or. (r2 < r2atat(iatjat))) then
      loverlap = .true.
      return
   end if
   energy = energy + uLJ(r2, iatjat)
   totforce = fLJ(r2,iatjat)
   forcea(1:3) = forcea(1:3) + totforce*dr(1:3)
   forceb(1:3) = forceb(1:3) - totforce*dr(1:3)
   virial     = virial - (totforce * r2)
end subroutine ufvLJdr

subroutine uLJdr(dr, iatjat, energy, loverlap) !calculated the energy from a distance and atom-type combination
   use MolModule, only: rcut2, r2umin, r2atat
   implicit none

   real(8), intent(in)  :: dr(3) !note that the PBC are neede to be already applied!
   integer(4), intent(in)  :: iatjat
   real(8), intent(inout)  :: energy
   logical, intent(out)  :: loverlap
   real(8)  :: r2

   loverlap = .false.
   r2 = sum(dr**2)
   if ((r2 < r2umin(iatjat)) .or. (r2 < r2atat(iatjat))) then
      loverlap = .true.
      return
   end if
   if (r2 .le. rcut2) then
      energy = energy + uLJ(r2, iatjat)
   end if
end subroutine uLJdr


subroutine UFlexLJ(utwob, utot, force, virial)
   !todo(pascal)  : add lclist and lmonoatom
   !use MolModule, only: lclist
   !use MolModule, only: lmonoatom
   use MolModule, only: nptpt
   implicit none

   real(8), allocatable, intent(inout)  :: utwob(:)
   real(8), intent(inout)               :: utot
   real(8), allocatable, intent(inout)  :: force(:,:)
   real(8), intent(inout)               :: virial
   real(8)  :: dutotold

   dutotold = utwob(0)
   !if(lclist) then
      !if(lmonoatom) then
         !call UFlexLJCellMono(utwob, utot, force, virial)
      !else
         !call UFlexLJCellPoly(utwob, utot, force, virial)
      !end if
   !else
      !if(lmonoatom) then
         call UFlexLJMono(utwob, force, virial)
      !else
         !call UFlexLJPoly(utwob, utot, force, virial)
      !end if
   !end if
   utwob(0) = sum(utwob(1:nptpt))
   utot = utot + utwob(0) - dutotold
end subroutine

subroutine StopUFlexLJ(i, j, ri, rj, itjt, dr)
   use MolModule, only: boxlen2, dpbc, uout, myid, r2umin, r2atat
   implicit none
   integer(4), intent(in)  :: i
   integer(4), intent(in)  :: j
   real(8), intent(in)     :: ri(3)
   real(8), intent(in)     :: rj(3)
   integer(4), intent(in)  :: itjt
   real(8), intent(in)     :: dr(3)

   character(40), parameter :: txroutine ='StopUFlexLJ'
   real(8)  :: r2

   r2 = sum(dr**2)

   write(uout,'(a,i5)')  'i', i
   write(uout,'(a,i5)')  'j', j
   write(uout,'(a,3e15.5)') 'r(i)         = ', ri(1:3)
   write(uout,'(a,3e15.5)') 'r(j)         = ', rj(1:3)
   write(uout,'(a,3e15.5)') 'boxlen2      = ', boxlen2
   write(uout,'(a,3e15.5)') 'dpbc         = ', dpbc
   write(uout,'(a,i15)')    'uout         = ', uout
   write(uout,'(a,i15)')    'myid         = ', myid
   write(uout,'(a,i15)')    'itjt         = ', itjt
   write(uout,'(a,e15.5)')  'r2           = ', r2
   write(uout,'(a,e15.5)')  'r2umin(itjt) = ', r2umin(itjt)
   write(uout,'(a,e15.5)')  'r2atat(itjt) = ', r2atat(itjt)
   call Stop(txroutine, 'r2 < r2umin or r2 < r2atat', uout)

end subroutine StopUFlexLJ

subroutine UFlexLJMono(utwob, force, virial)

   use MolModule, only: ro, iptpn, iptpt
   use MolModule, only: ipnploc, nneighpn
   use MolModule, only: lmc
   use MolModule, only: jpnlist, npmyid
   implicit none
   real(8), allocatable, intent(inout)  :: utwob(:)
   real(8), allocatable, intent(inout)  :: force(:,:)
   real(8), intent(inout)               :: virial

   integer(4)  :: ip, jp, ipt, jpt, iptjpt, iploc, jploc
   real(8)  :: dr(3)
   logical  :: loverlap

   do iploc = 1, npmyid !get particles from nlist
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
         call PBC(dr(1), dr(2), dr(3))
         call ufvLJdr(dr, iptjpt, utwob(iptjpt), force(1:3,ip), force(1:3,jp), virial, loverlap)
         if(loverlap)   call StopUFlexLJ(ip, jp, ro(1:3,ip), ro(1:3,jp), iptjpt, dr(1:3))
      end do
   end do

end subroutine UFlexLJMono

subroutine DUFlexLJ(dutwob, dutot, lhsoverlap)
   !todo(pascal)  : add lclist and lmonoatom
   !use MolModule, only: lclist
   !use MolModule, only: lmonoatom
   use MolModule, only: nptpt

   implicit none

   real(8), allocatable, intent(inout) :: dutwob(:)
   real(8), intent(inout)              :: dutot
   logical, intent(inout)              :: lhsoverlap
   real(8)  :: dutwbold

   lhsoverlap =.true.
   dutwbold = dutwob(0)
   !if(lclist) then
      !if(lmonoatom) then
         !call UFlexLJCellMono(utwob, utot, force, virial)
      !else
         !call UFlexLJCellPoly(utwob, utot, force, virial)
      !end if
   !else
      !if(lmonoatom) then
         call DUFlexLJMono(dutwob, lhsoverlap)
      !else
         !call UFlexLJPoly(utwob, utot, force, virial)
      !end if
   !end if

   dutwob(0) = sum(dutwob(1:nptpt))
   dutot = dutot + dutwob(0) - dutwbold
end subroutine DUFlexLJ

subroutine DUFlexLJMono(dutwob, lhsoverlap)

   use MolModule, only: ro, rotm, iptpn, iptpt, rcut2
   use MolModule, only: jpnlist, nneighpn
   use MolModule, only: ipnptm, nptm, lptm, lptmdutwob
   use MolModule, only: nproc
   implicit none
   real(8), allocatable, intent(inout)  :: dutwob(:)
   logical, intent(inout)               :: lhsoverlap

   integer(4)  :: ip, jp, ipt, jpt, iptjpt, iploc, jploc
   real(8)  :: dr(3), r2

! do first new then old energy as, overlaps will ocure in the new energy, and then the old energy is not needed
! ----------------- NEW ENERGY ------------

! ... contribution from pairs where only one particle is moved
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dr(1:3) = rotm(1:3,iploc)-ro(1:3,jp)
         call PBC(dr(1), dr(2), dr(3))
         call uLJdr(dr, iptjpt, dutwob(iptjpt), lhsoverlap)
         if(lhsoverlap) return
      end do
   end do

! ... contribution from pairs where both particle is moved

   if (lptmdutwob) then
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm, nproc ! increment with nproc for parallel simulations
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dr(1:3) = rotm(1:3,iploc)-rotm(1:3,jp)
            call PBC(dr(1), dr(2), dr(3))
            call uLJdr(dr, iptjpt, dutwob(iptjpt), lhsoverlap)
            if(lhsoverlap) return
         end do
      end do
   end if

   lhsoverlap =.false.

! ----------------- OLD ENERGY ------------

! ... contribution from pairs where only one particle is moved
   do iploc = 1, nptm
      ip = ipnptm(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(ip)
         jp = jpnlist(jploc,ip)
         if (lptm(jp)) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dr(1:3) = ro(1:3,iploc)-ro(1:3,jp)
         call PBCr2(dr(1), dr(2), dr(3), r2)
         !as the old configuration should be free of overlaps one can skip the checks
         if(r2 < rcut2) then
            dutwob(iptjpt) = dutwob(iptjpt) - uLJ(r2,iptjpt)
         end if
      end do
   end do

! ... contribution from pairs where both particle is moved

   if (lptmdutwob) then
      do iploc = 1, nptm
         ip = ipnptm(iploc)
         ipt = iptpn(ip)
         do jploc = iploc+1, nptm, nproc ! increment with nproc for parallel simulations
            jp = ipnptm(jploc)
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dr(1:3) = ro(1:3,iploc)-ro(1:3,jp)
            call PBCr2(dr(1), dr(2), dr(3), r2)
            !as the old configuration should be free of overlaps one can skip the checks
            if(r2 < rcut2) then
               dutwob(iptjpt) = dutwob(iptjpt) - uLJ(r2,iptjpt)
            end if
         end do
      end do
   end if

end subroutine DUFlexLJMono

end module flexLJEnergyModule
