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

!************************************************************************
!*                                                                      *
!*     CoordinateModule                                                 *
!*                                                                      *
!************************************************************************

! ... module for coordinate

module CoordinateModule

   use MolModule

   integer(4), parameter :: ntrydef = 100     ! default number of trials per particle

   character(20), allocatable :: txsetconf(:) ! select way of generating start configuration
   integer(4), allocatable    :: nucell(:,:)  ! number of unit cells in x-, y-, and z-direction
   real(8), allocatable       :: rclow(:,:)   ! box in which particles are set (one corner)
   real(8), allocatable       :: rcupp(:,:)   ! box in which particles are set (diagonal corner)

   real(8), allocatable       :: roshift(:,:) ! shift of lattice points in a unit cell in fraction of lattice length
   real(8), allocatable       :: radatset(:)  ! hard-core radius of atom type used when setting particles
   logical, allocatable       :: lranori(:)   ! logical flag for random particle orientation (some options of txsetconf)
   real(8), allocatable       :: bondscl(:)   ! bond length scaling factor for chains
   real(8), allocatable       :: anglemin(:)  ! minimum angle between consecutive particles in a chain

   integer(4)    :: iptnode                   ! particle type of nodes (for SetPeriodicNetwork)
   integer(4)    :: ictstrand                 ! chain type of strands

   real(8),      allocatable :: rnwt(:)       ! radius of network type inwt [Allocate with nnwt]
   character(8), allocatable :: txoriginnwt(:)! network center of network type inwt ("origin" and "random") [Allocate with nnwt]

   real(8)       :: radlimit(2)               ! lower and upper radial limit of particles (for SetCoreCell)

   integer(4)    :: itestcoordinate           ! =1, call of TestMakeCrossLink

end module CoordinateModule

!************************************************************************
!*                                                                      *
!*     Coordinate                                                       *
!*                                                                      *
!************************************************************************

! ... handle coordinates and velocities

subroutine Coordinate(iStage)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Coordinate'
   character(80), parameter :: txheading ='initial configuration data'
   logical    :: GetlSetVel, GetlZeroMom
   integer(4) :: ip

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

   case (iWriteInput)

      if (master) call WriteHead(2, txheading, uout)

! ... get initial coordinates

      if (txstart == 'setconf') then                    ! configuration from SetConfiguration

         if (lweakcharge) call SetChargeWeakChargeCase
         call SetConfiguration

         if (lmd) then
            call OriToQua(np, 1, np, ori, qua)
            if (master) write(uout,'()')
            if (GetlSetVel()) then
               call SetVel
               if (GetlZeroMom()) call SetZeroMom
               if (master) write(uout,'(a)') 'velocities: Maxwell distribution'
            else
               rod = Zero
               quad = Zero
               angvelo = Zero
               if (master) write(uout,'(a)') 'velocities: zero'
            end if
         end if

      else if (txstart == 'readfin') then               ! coordinates from fin

         if (lweakcharge) call SetChargeWeakChargeCase

         rewind(uin)
         read(uin,*) (ro(1:3,ip),ori(1:3,1:3,ip),ip = 1,np)

         if (lclink) call CrossLink_Steffi        ! fix for reading cross-link bonds (used for periodic networks)

         call OrthoOri(np, 1, np, ori, 1.0d-4, uout)
         call SetAtomProp(1, np, .false.)
         if (master) write(uout,'(a)') 'coordinates: read from fin'
         if (lmd) then
            call OriToQua(np, 1, np, ori, qua)
            if (master) write(uout,'()')
            if (GetlSetVel()) then
               call SetVel
               if (GetlZeroMom()) call SetZeroMom
               if (master) write(uout,'(a)') 'velocities: maxwell distribution'
           else
              read(uin,*) (rod(1:3,ip),angvelo(1:3,ip),ip = 1,np)
              call AngVelToQuaVel(np, 1, np, qua, angvelo, quad)
              if (master) write(uout,'(a)') 'velocities: read from fin'
            end if
         end if

      else if (txstart == 'zero') then                  ! coordinates from fcnf

         call IOCnf('read')
         call SetAtomProp(1, np, .false.)
         if (master) then
            write(uout,'(a)') 'box lengths and coordinates: read from fcnf'
            write(uout,'(a,t35,i15)') 'current seed                   = ', iseed
            if (lbcbox) then
               write(uout,'(a,t35,3(f10.3,2x))') 'current box lengths (x,y,z)    = ', boxlen
            else if (lbcrd .or. lbcto) then
               write(uout,'(a,t35,3(f10.3,2x))') 'current side length    = ', cellside
            end if
         end if
         if (lmd) then
            if (master) write(uout,'()')
            if (GetlSetVel()) then
               call SetVel
               if (GetlZeroMom()) call SetZeroMom
               if (master) write(uout,'(a)') 'velocities: maxwell distribution'
            else
               if (master) write(uout,'(a)') 'velocities: read from fcnf'
            end if
         end if

!        if (txuser == 'niklas') call Coordinate_Niklas   ! paper 1

      else if (txstart == 'continue') then           ! coordinates from fcnf

         call IOCnf('read')
         call SetAtomProp(1, np, .false.)
         if (master) then
            write(uout,'(a)') 'seed, box lengths, and coordinates: read from fcnf'
            write(uout,'(a,t35,i15)') 'current seed                   = ', iseed
            if (lbcbox) then
               write(uout,'(a,t35,3(f10.3,2x))') 'current box lengths (x,y,z)    = ', boxlen
            else if (lbcrd .or. lbcto) then
               write(uout,'(a,t35,3(f10.3,2x))') 'current side length    = ', cellside
            end if
         end if
         if (lmd.and.master) then
            write(uout,'()')
            write(uout,'(a)') 'velocities: read from fcnf'
         end if

      else

         call Stop(txroutine, 'unsupported value of start', uout)

      end if

! ... calculate initial linear and angular moments

      if (lmd.and.master) then
         call GetLinMom
         call GetAngMom
         write(uout,'()')
         call WriteLinMom
         call WriteAngMom
      end if

! ... scale box length and coordinates ?

      if (txstart == 'readfin' .or. txstart == 'zero') then
         if (abs(lenscl-One) > 1.0d-4) then
            if (lbcbox) then
               boxlen = boxlen*lenscl
            else if (lbcrd .or. lbcto) then
               cellside = cellside*lenscl
            else
               call Stop(txroutine,'invalid geometry for length rescaling',uout)
            end if
            call SetBoxParam
            if (lewald) call EwaldSetup
            ro = ro*lenscl
            call SetAtomPos(1, np, .false.)
            if (master) then
                write(uout,'()')
                write(uout,'(a,t35,f5.2)') 'box length and coord. scaled with = ', lenscl
            end if
         end if
      end if

! ... get nstep1beg

      if (txstart == 'continue') then
         nstep1beg = nstep1done + 1
      else
         nstep1beg = 1
      endif

      if (master) then
         write(uout,'()')
         write(uout,'(a,t35,i5)') 'first macrostep                = ', nstep1beg
      end if
      call WritePartAtomProp

      if (master) call FileFlush(uout)

   case (iAfterMacrostep)

      call IOCnf('write')

   case (iAfterSimulation)

      call WritePartAtomProp

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

! chain translate such that the particle with smallest z-coordinate becomes located at zaim

subroutine Coordinate_Niklas   ! Niklas 2006-12-20
   real :: zaim, zlow
   if (lbd) then
     write(uout, *) 'chain translation towards surface enabled'
     zaim = -boxlen2(3) + 12.0
     zlow = minval(ro(3,1:np))
     ro(3,1:np) = (zaim-zlow) + ro(3,1:np)
   end if
end subroutine Coordinate_Niklas

subroutine CrossLink_Steffi    ! Steffi
    integer(4) :: jp,i
    bondcl = 0
    iptnode = 1
    nbondcl = 0
    read(uin,*) (nbondcl(ip),bondcl(1:4,ip),ip=1,nppt(iptnode))
    do ip=1,nppt(iptnode)
       do i=1,nbondcl(ip)
          jp=bondcl(i,ip)
          if (jp==0) cycle
          ncl = ncl + 1
          nbondcl(jp)=nbondcl(jp)+1
          bondcl(nbondcl(jp),jp)=ip
       end do
    end do
end subroutine CrossLink_Steffi

end subroutine Coordinate

!************************************************************************
!*                                                                      *
!*     SetChargeWeakChargeCase                                          *
!*                                                                      *
!************************************************************************

! ... initiate charge of for the weak-charge case

subroutine SetChargeWeakChargeCase
   use CoordinateModule
   implicit none
   where (abs(az(1:na)) > 1d-10)
      laz(1:na) = .true.  ! initiate laz for weak and strong charges
   elsewhere
      laz(1:na) = .false. !
   endwhere
end subroutine SetChargeWeakChargeCase

!************************************************************************
!*                                                                      *
!*     SetConfiguration                                                 *
!*                                                                      *
!************************************************************************

! ... generate a start configuration

subroutine SetConfiguration

   use CoordinateModule
   implicit none

   character(40), parameter :: txroutine ='SetConfiguration'
   character     :: charfmt*1
   logical       :: lsetconf
   integer(4)    :: ipt, ict, iptgen0, igen, isign, imagn, inwt

   external SetOrigin, SetSq2D, SetHex2D, SetPC, SetBCC, SetFCC, SetSM2, SetDiamond, SetH2O, SetN2, SetBenzene

   namelist /nmlSetConfiguration/ txsetconf, nucell, rclow, rcupp, roshift,                          &
                                  radatset, lranori,  bondscl, anglemin,                             &
                                  ngen, ictgen, nbranch, ibranchpbeg, ibranchpinc,                   &
                                  iptnode, ictstrand,                                                &
                                  rnwt, txoriginnwt,                                                 &
                                  radlimit,                                                          &
                                  itestcoordinate

   if (.not.allocated(txsetconf)) then 
      allocate(txsetconf(npt), nucell(3,npt), rclow(3,npt), rcupp(3, npt),       &
      roshift(3,npt), radatset(nat), lranori(npt), bondscl(nct), anglemin(nct))
      txsetconf   = ""
      nucell      = 0
      rclow       = 0.0E+00
      rcupp       = 0.0E+00
      roshift     = 0.0E+00
      radatset    = 0.0E+00
      lranori     = .false.
      bondscl     = 0.0E+00
      anglemin    = 0.0E+00
   end if
   if (.not.allocated(radatset)) then 
      allocate(radatset(nat))
      radatset = 0.0E+00
   end if

   if (lnetwork) then
      if(.not.allocated(rnwt)) then
         allocate(rnwt(nnwt))
      end if
      if(.not.allocated(txoriginnwt)) then
         allocate(txoriginnwt(nnwt))
      end if
   end if

! ... read input data

   nucell(1:3,1:npt) = 0
   if (lbcbox) then
      rclow(1,1:npt) = -boxlen2(1)
      rclow(2,1:npt) = -boxlen2(2)
      rclow(3,1:npt) = -boxlen2(3)
      rcupp(1,1:npt) = +boxlen2(1)
      rcupp(2,1:npt) = +boxlen2(2)
      rcupp(3,1:npt) = +boxlen2(3)
   else if (lbcrd) then
      rclow(1,1:npt) = -boxlen2(1)
      rclow(2,1:npt) = -boxlen2(2)
      rclow(3,1:npt) = -boxlen2(3)
      rcupp(1,1:npt) = +boxlen2(1)
      rcupp(2,1:npt) = +boxlen2(2)
      rcupp(3,1:npt) = +boxlen2(3)
   else if (lbcto) then
      rclow(1,1:npt) = -boxlen2(1)
      rclow(2,1:npt) = -boxlen2(2)
      rclow(3,1:npt) = -boxlen2(3)
      rcupp(1,1:npt) = +boxlen2(1)
      rcupp(2,1:npt) = +boxlen2(2)
      rcupp(3,1:npt) = +boxlen2(3)
   else if (lbcsph) then
      rclow(1:3,1:npt) = -sphrad
      rcupp(1:3,1:npt) = +sphrad
   else if (lbccyl) then
      rclow(1,1:npt) = -cylrad
      rcupp(1,1:npt) = +cylrad
      rclow(2,1:npt) = -cylrad
      rcupp(2,1:npt) = +cylrad
      rclow(3,1:npt)   = -(Half*cyllen)
      rcupp(3,1:npt)   = +(Half*cyllen)
   else if (lbcell) then
      rclow(1,1:npt) = -ellrad(1)
      rclow(2,1:npt) = -ellrad(2)
      rclow(3,1:npt) = -ellrad(3)
      rcupp(1,1:npt) = +ellrad(1)
      rcupp(2,1:npt) = +ellrad(2)
      rcupp(3,1:npt) = +ellrad(3)
   end if
   roshift           = Zero
   radatset          = radat
   lranori           =.false.
   lfixedori         =.false.
   bondscl           = One
   anglemin          = Zero
   iptnode           = 0
   ictstrand         = 0
   itestcoordinate   = 0

   if (lnetwork) then
      rnwt(1:nnwt)        = 10.0
      txoriginnwt(1:nnwt) = 'origin'
   end if

   rewind(uin)
   read(uin,nmlSetConfiguration)

   do ipt = 1, npt
     call LowerCase(txsetconf(ipt))
   end do
   do inwt = 1, nnwt
      call LowerCase(txoriginnwt(inwt))
   end do

   anglemin = anglemin*sclang

   if (.not.allocated(lpset)) then 
      allocate(lpset(np_alloc))
      lpset = .false.
   end if
   lpset =.false.

! ... back compatibility

   where ((ictpt(1:npt) > 0) .and. (txsetconf(1:npt) == 'random')) txsetconf(1:npt) = 'chainrandom'
   where (txsetconf(1:npt) == 'chainrandompos')  txsetconf(1:npt) = 'chainrandomintori'

! ... if lclink: check consistencies and crosslinks variables and if diamongel set ncl

  if (lclink) then
     if (count(txsetconf == 'periodicnetwork') + count(txsetconf(:)(1:12) == 'hierarchical') + count(txsetconf == 'network') == 0) &
        call Stop(txroutine, 'lclink: no call of SetPeriodicNetwork, SetNetwork or SetHierarchical', uout)
     if (count(txsetconf == 'periodicnetwork') > 0) ncl = 4*nppt(iptnode)
  end if

! ... if hierarchical structures, copy particle type data from generation zero

   if (count(txsetconf(:)(1:12) == 'hierarchical') > 0) then
      igen = 0
      iptgen0 = iptpn(ipnsegcn(1,icnct(ictgen(igen))))
      do igen = 1, ngen
         ipt = iptpn(ipnsegcn(1,icnct(ictgen(igen))))
!         nucell(1:3,ipt) = nucell(1:3,iptgen0)
!         rclow(1:3,ipt) = rclow(1:3,iptgen0)
!         rcupp(1:3,ipt) = rcupp(1:3,iptgen0)
!         roshift(1:3,ipt) = roshift(1:3,iptgen0)
      end do
   end if

! ... write input data

   if (master) then
      call SignMagn(maxval(abs(rclow)), isign, imagn)
      write(charfmt,'(i1)') 3 - imagn
      write(uout,'(a)') 'start configuration from SetConfiguration'
      write(uout,'(a)') '-----------------------------------------'
      write(uout,'(a,t15,a,t32,a,t60,a)')                        &
        'part. type', 'routine', 'no of unit cells (x,y,z)', ' lower (x,y,z)      upper (x,y,z)      shift (x,y,z)'
      write(uout,'(a,t15,a,t32,a,t60,a)')                        &
        '----------', '-------', '------------------------', ' -------------      -------------      -------------'
      do ipt = 1, npt
         write(uout,'(i2,2x,a,t15,a,t35,3(i5),t58,3(3f6.'//charfmt//',x))') &
         ipt, txpt(ipt), txsetconf(ipt), nucell(1:3,ipt), rclow(1:3,ipt), rcupp(1:3,ipt), roshift(1:3,ipt)
      end do
      if (count(radatset(1:nat) /= radat(1:nat)) > 0) then
         write(uout,'()')
         write(uout,'(a,t45,10f8.3)') 'radius of atom types when setting pos.  =', radatset(1:nat)
      end if
      if (count(lranori(1:npt)) > 0) then
         write(uout,'()')
         write(uout,'(a,t45,8g15.5)') 'random particle oritentation            = ',lranori(1:npt)
      end if
      if (lchain) then
         write(uout,'()')
         write(uout,'(a,t45,6f8.3)') 'bond length scaling factor              = ', bondscl(1:nct)
         write(uout,'(a,t45,6f8.3)') 'minimum angle between consecutive beads = ', anglemin(1:nct)/sclang
      end if
      if (lclink) then
         if (count(txsetconf == 'periodicnetwork') > 0) then
            write(uout,'()')
            write(uout,'(a,t45,i8)') 'chain type of strands                    = ', ictstrand
            write(uout,'(a,t45,i8)') 'particle type of nodes                   = ', iptnode
            write(uout,'(a,t45,i8)') 'number of nodes                          = ', nppt(iptnode)
            write(uout,'(a,t45,i8)') 'maximum number of crosslinks of the nodes= ', maxnbondcl(iptnode)
         end if
         if (count(txsetconf == 'network') > 0) then
            write(uout,'()')
            write(uout,'(a,t60,2i8)')    'number of network types                                 = ', nnwt
            write(uout,'(a,t60,2i8)')    'number of networks of the network types                 = ', nnwnwt(1:nnwt)
            write(uout,'(a,t60,2f8.2)')  'radius of networks of the network types                 = ', rnwt(1:nnwt)
            write(uout,'(a,t60,2i8)')    'particle type of nodes of the network types             = ', iptclnwt(1:nnwt)
            write(uout,'(a,t60,2x,2a8)') 'position of gel center of the different network types   = ', txoriginnwt(1:nnwt)
            do ipt = 1, npt   ! one particle type can only be used in one network type
               if (count(ipt == iptclnwt(1:nnwt)) > 1) call stop(txroutine, 'error in iptclnwt', uout)
            end do
            do ict = 1, nct   ! one chain type can only be used in one network type
               if (count(ncctnwt(ict,1:nnwt) > 0) > 1) call stop(txroutine, 'chain type used for more than one network type', uout)   ! Cornelius Hofzumahaus 
            end do
         end if
      end if
   end if

! ... initialize for crosslinked structures made by MakeCrossLink

    if (lclink .and. .not.lhierarchical) then
       ncl = 0
       nbondcl = 0
       bondcl = 0
    end if

! ... set particles

   do ipt = 1, npt
      if (txsetconf(ipt) == 'origin') then
         nucell(1:3,ipt) = 1
         call SetLattice(ipt, SetOrigin)
      else if (txsetconf(ipt) == 'sq2dlattice') then
         call SetLattice(ipt, SetSq2D)
      else if (txsetconf(ipt) == 'hex2dlattice') then
         call SetLattice(ipt, SetHex2D)
      else if (txsetconf(ipt) == 'pclattice') then
         call SetLattice(ipt, SetPC)
      else if (txsetconf(ipt) == 'bcclattice') then
         call SetLattice(ipt, SetBCC)
      else if (txsetconf(ipt) == 'fcclattice') then
         call SetLattice(ipt, SetFCC)
      else if (txsetconf(ipt) == 'sm2lattice') then
         call SetLattice(ipt, SetSM2)
      else if (txsetconf(ipt) == 'diamondlattice') then
         call SetLattice(ipt, SetDiamond)
      else if (txsetconf(ipt) == 'h2olattice') then
         call SetLattice(ipt, SetH2O)
      else if (txsetconf(ipt) == 'n2lattice') then
         call SetLattice(ipt, SetN2)
      else if (txsetconf(ipt) == 'benzenelattice') then
         call SetLattice(ipt, SetBenzene)
      else if (txsetconf(ipt) == 'random') then
         call SetRandom(ipt)
      else if (txsetconf(ipt) == 'randomfixori') then
         call SetRandomFixOri(ipt)

      else if (txsetconf(ipt) =='chainline') then
         call SetChainLine(ipt)
      else if (txsetconf(ipt) =='chaincircle') then
         call SetChainCircle(ipt)
      else if (txsetconf(ipt) == 'chainrandom') then
         call SetChainRandom(ipt)
      else if (txsetconf(ipt) == 'chainrandompos' .or. txsetconf(ipt) == 'chainrandomintori') then
         call SetChainRandomIntOri(ipt)
      else if (txsetconf(ipt) == 'sphbrushlattice') then
         call SetSphBrush('lattice')
      else if (txsetconf(ipt) == 'sphbrushrandom') then
         call SetSphBrush('random')
      else if (txsetconf(ipt) == 'planbrushrandom') then
         call SetPlanarBrush('random')
      else if (txsetconf(ipt) == 'hierarchicallattice') then
         call SetHierarchical('lattice')
      else if (txsetconf(ipt) == 'hierarchicalrandom') then
         call SetHierarchical('random')
      else if (txsetconf(ipt) == 'periodicnetwork') then
         call SetPeriodicNetwork(ipt)
      else if(txsetconf(ipt) == 'network') then
         call SetNetwork(ipt)
      else if (txsetconf(ipt) == 'coreshell') then
         call SetCoreShell(ipt)
      else
         call SetConfigurationUser(ipt, lsetconf)
         if (.not.lsetconf) call Stop(txroutine, 'unsupported value of txsetconf(ipt)', uout)
      end if
   end do

   if (master) then
      if (lclink) then
         if (count(txsetconf == 'periodicnetwork') > 0) then
            write(uout,'(a,t45,i8)') 'number of crosslinks                     = ', ncl
         end if
         if (count(txsetconf == 'network') > 0) then
            write(uout,'(a,t60,i8)') 'number of crosslinks                     = ', ncl
         end if
      end if
   end if

   if (lpolarization .or. ldipole .or. ldipolesph) call SetAtomDipMom(1, np)  ! set dipole moments
   if (lpolarization) call SetAtomPolTens(1, np)                          ! set polarization tensor

end subroutine SetConfiguration

!************************************************************************
!*                                                                      *
!*     SetLattice                                                       *
!*                                                                      *
!************************************************************************

! ... generate a lattice configuration

subroutine SetLattice(ipt, latticesub)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt       ! particle type
   external latticesub                 ! name of subroutine providing nlp, rol, and oril

   character(40), parameter :: txroutine ='SetLattice'
   integer(4) :: ip                    ! number of particle to be set
   integer(4) :: nset                  ! number of setted particles of a given type
   integer(4) :: nlp                   ! number of particles in a unit cell
   real(8)    :: rol(3,8)              ! coordinates of the lattice points in a unit cell ranging from (0,0,0) to (1,1,1)
   real(8)    :: oril(3,3,8)           ! orientation of particle frame at different lattice points in a unit cell
   integer(4) :: ix, iy, iz, ilp
   real(8)    :: rlen(3), rorigin(3), dxpbc, dypbc, dzpbc, ddelta

   real(8), parameter :: theta = 0.40d0
   real(8), parameter :: sfac = 1.348d0
   real(8), parameter ::xys = (sfac-1.0d0)/1.414213562d0
   real(8) :: tmp, cosa, sina, x, y, gamma

   if (txsetconf(1) == 'sm2lattice') then          ! special care for sm2 lattice
      if (npt > 1) call stop(txroutine, 'sm2 lattice and npt > 1',uout)
      tmp = xys*(boxlen2(1)+boxlen2(2))
      sina = tmp/sqrt(tmp**2+(tmp+boxlen(1))**2)
      cosa = sqrt(one-sina**2)
      boxlen(1) = sqrt(tmp**2+(tmp+boxlen(1))**2)
      boxlen(2) = boxlen(1)*(cosa**2-sina**2)
      call SetBoxParam
   end if

   ddelta = 1.0d-10
   if (lbcrd .or. lbcto) ddelta = 1.0d-10    ! to ensure correct round off (important for RD and TO b.c.)

   rol = Zero                                                  ! initiate parameter
   oril = Zero                                                 ! initiate parameter
   call latticesub(nlp, rol, oril, cosa, sina, theta)          ! get parameters describing a unit cell
   rlen(1:3) = (rcupp(1:3,ipt)-rclow(1:3,ipt))/nucell(1:3,ipt) ! prepare rlen
   rorigin(1:3) = Half*(rclow(1:3,ipt)+rcupp(1:3,ipt))         ! prepare rorigin

   gamma = 30*DegToRad
   gamma = Zero

! ... loop over all lattice points and try to set particles at those

   nset = 1
   do ilp = 1, nlp
      do iz = 0, nucell(3,ipt)-1
         do iy = 0, nucell(2,ipt)-1
            do ix = 0, nucell(1,ipt)-1
               ip = nset-1+ipnpt(ipt)

! ... set rotation invariant point

               ro(1,ip) = ix-Half*nucell(1,ipt) + (rol(1,ilp) + roshift(1,ipt))
               ro(2,ip) = iy-Half*nucell(2,ipt) + (rol(2,ilp) + roshift(2,ipt))
               ro(3,ip) = iz-Half*nucell(3,ipt) + (rol(3,ilp) + roshift(3,ipt))

   !           write(*,*)  iy, gamma, tan(gamma), rlen(2), tan(gamma)*rlen(2)*iy
               ro(1,ip) = rorigin(1) + ddelta + ro(1,ip)*rlen(1) + tan(gamma)*rlen(2)*iy
               ro(2,ip) = rorigin(2) + ddelta + ro(2,ip)*rlen(2)
               ro(3,ip) = rorigin(3) + ddelta + ro(3,ip)*rlen(3)

               if (txsetconf(1) == 'sm2lattice') then          ! special care for sm2 lattice
                  tmp = (ro(1,ip)+ro(2,ip))
                  x = ro(1,ip)+xys*tmp
                  y = ro(2,ip)+xys*tmp
                  ro(1,ip) = cosa*x+sina*y
                  ro(2,ip) =-sina*x+cosa*y
                  call PBC(ro(1,ip), ro(2,ip), ro(3,ip))
               end if

! ... test if outside the simulation cell

              if (lbcbox) then
                  call PBC(ro(1,ip), ro(2,ip), ro(3,ip))
              !   if (abs(ro(1,ip)) > boxlen2(1)) cycle
              !   if (abs(ro(2,ip)) > boxlen2(2)) cycle
              !   if (abs(ro(3,ip)) > boxlen2(3)) cycle
              else if (lbcsph) then
                 if (sum(ro(1:3,ip)**2) > sphrad2) cycle
              else if (lbccyl) then
                 if (sum(ro(1:2,ip)**2) >= cylrad2-0.1d0) cycle
                 if (abs(ro(3,ip)) > Half*cyllen-0.1d0) cycle
              else if (lbcell) then
                 if (sum((ro(1:3,ip)*ellradi(1:3))**2) >= 0.9999d0) cycle
              else
                 call PBC2(ro(1,ip),ro(2,ip),ro(3,ip),dxpbc,dypbc,dzpbc)
                 if (abs(dxpbc) > Zero) cycle
                 if (abs(dypbc) > Zero) cycle
                 if (abs(dzpbc) > Zero) cycle
              end if

! ... set particle frame and atoms

               ori(1:3,1:3,ip) = oril(1:3,1:3,ilp)
               if (lranori(ipt)) call SetPartOriRandom(iseed,ori(1,1,ip))   ! set random particle orientation
               call OrthoOri(np, ip, ip, ori, 1.0d-4, uout)
               call SetAtomPos(ip, ip, .false.)

! ... configuration accepted

               lpset(ip) =.true.
               if (nset == nppt(ipt)) return
               nset = nset+1

            end do
         end do
      end do
   end do

   if (master) then
      write(uout,'(a,i5)') 'ipt       = ', ipt
      write(uout,'(a,i5)') 'nppt(ipt) = ', nppt(ipt)
      write(uout,'(a,i5)') 'nset      = ', nset
   end if
   call Stop(txroutine, 'lattice set failed', uout)

end subroutine SetLattice

!************************************************************************
!*                                                                      *
!*     SetOrigin                                                        *
!*                                                                      *
!************************************************************************

! ... generate a point in the center of the unit cell

subroutine SetOrigin(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   nlp = 1
   rol(1:3,1) = 0.5d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
end subroutine SetOrigin

!************************************************************************
!*                                                                      *
!*     SetSq2D                                                          *
!*                                                                      *
!************************************************************************

! ... generate a 2d square lattice

subroutine SetSq2D(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   nlp = 1
   rol(1:3,1) = 0.0d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
end subroutine SetSq2D

!************************************************************************
!*                                                                      *
!*     SetHex2D                                                         *
!*                                                                      *
!************************************************************************

! ... generate a 2d hexagonal lattice

subroutine SetHex2D(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)

   nlp = 2
   rol(1:3,1) = 0.0d0
   rol(1,2) = 0.5d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.0d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
   oril(1,1,2) = 1.0d0
   oril(2,2,2) = 1.0d0
   oril(3,3,2) = 1.0d0
end subroutine SetHex2D

!************************************************************************
!*                                                                      *
!*     SetPC                                                            *
!*                                                                      *
!************************************************************************

! ... generate a primitive cubic lattice

subroutine SetPC(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   nlp = 1
   rol(1:3,1) = 0.0d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
end subroutine SetPC

!************************************************************************
!*                                                                      *
!*     SetBCC                                                           *
!*                                                                      *
!************************************************************************

! ... generate a body-centered cubic lattice

subroutine SetBCC(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   nlp = 2
   rol(1:3,1) = 0.0d0
   rol(1:3,2) = 0.5d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
   oril(1,1,2) = 1.0d0
   oril(2,2,2) = 1.0d0
   oril(3,3,2) = 1.0d0
end subroutine SetBCC

!************************************************************************
!*                                                                      *
!*     SetFCC                                                           *
!*                                                                      *
!************************************************************************

! ... generate a face-centered cubic lattice

subroutine SetFCC(nlp, rol, oril)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)

   nlp = 4
   rol(1:3,1) = 0.0d0
   rol(1,2) = 0.0d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.5d0
   rol(1,3) = 0.5d0
   rol(2,3) = 0.0d0
   rol(3,3) = 0.5d0
   rol(1,4) = 0.5d0
   rol(2,4) = 0.5d0
   rol(3,4) = 0.0d0
   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
   oril(1,1,2) = 1.0d0
   oril(2,2,2) = 1.0d0
   oril(3,3,2) = 1.0d0
   oril(1,1,3) = 1.0d0
   oril(2,2,3) = 1.0d0
   oril(3,3,3) = 1.0d0
   oril(1,1,4) = 1.0d0
   oril(2,2,4) = 1.0d0
   oril(3,3,4) = 1.0d0
end subroutine SetFCC

!************************************************************************
!*                                                                      *
!*     SetSM2                                                           *
!*                                                                      *
!************************************************************************

! ... generate a sm2 lattice

subroutine SetSM2(nlp, rol, oril, cosa, sina, theta)
   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   real(8),    intent(in)  :: cosa
   real(8),    intent(in)  :: sina
   real(8),    intent(in)  :: theta
   real(8) :: x, y
   real(8), parameter :: fac = 0.785398164d0
   nlp = 4
   rol(1:3,1) = 0.0d0
   rol(1,2) = 0.0d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.5d0
   rol(1,3) = 0.5d0
   rol(2,3) = 0.0d0
   rol(3,3) = 0.5d0
   rol(1,4) = 0.5d0
   rol(2,4) = 0.5d0
   rol(3,4) = 0.0d0

   x = cosa*cos(0.785398164-theta)-sina*sin(0.785398164-theta)
   y = sina*cos(0.785398164-theta)+cosa*sin(0.785398164-theta)
   oril(3,1,1) = 1.0d0
   oril(1,2,1) = x
   oril(2,2,1) = -y
   oril(1,3,1) = x
   oril(2,3,1) = y
   oril(3,1,4) = 1.0d0
   oril(1,2,4) = x
   oril(2,2,4) = -y
   oril(1,3,4) = x
   oril(2,3,4) = y

   x = cosa*cos(0.785398164+theta)-sina*sin(0.785398164+theta)
   y = sina*cos(0.785398164+theta)+cosa*sin(0.785398164+theta)
   oril(3,1,2) = 1.0d0
   oril(1,2,2) = x
   oril(2,2,2) = -y
   oril(1,3,2) = x
   oril(2,3,2) = y
   oril(3,1,3) = 1.0d0
   oril(1,2,3) = x
   oril(2,2,3) = -y
   oril(1,3,3) = x
   oril(2,3,3) = y

end subroutine SetSM2

!************************************************************************
!*                                                                      *
!*     SetDiamond                                                       *
!*                                                                      *
!************************************************************************

! ... generate a cubic diamond lattice

subroutine SetDiamond(nlp, rol, oril)

   implicit none
   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,8)
   real(8),    intent(out) :: oril(3,3,8)

   nlp = 8
   rol(1:3,1) = 0.0d0
   rol(1,2) = 0.0d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.5d0
   rol(1,3) = 0.5d0
   rol(2,3) = 0.0d0
   rol(3,3) = 0.5d0
   rol(1,4) = 0.5d0
   rol(2,4) = 0.5d0
   rol(3,4) = 0.0d0
   rol(1,5) = 0.25d0
   rol(2,5) = 0.25d0
   rol(3,5) = 0.25d0
   rol(1,6) = 0.25d0
   rol(2,6) = 0.75d0
   rol(3,6) = 0.75d0
   rol(1,7) = 0.75d0
   rol(2,7) = 0.25d0
   rol(3,7) = 0.75d0
   rol(1,8) = 0.75d0
   rol(2,8) = 0.75d0
   rol(3,8) = 0.25d0

   oril(1,1,1) = 1.0d0
   oril(2,2,1) = 1.0d0
   oril(3,3,1) = 1.0d0
   oril(1,1,2) = 1.0d0
   oril(2,2,2) = 1.0d0
   oril(3,3,2) = 1.0d0
   oril(1,1,3) = 1.0d0
   oril(2,2,3) = 1.0d0
   oril(3,3,3) = 1.0d0
   oril(1,1,4) = 1.0d0
   oril(2,2,4) = 1.0d0
   oril(3,3,4) = 1.0d0
   oril(1,1,5) = 1.0d0
   oril(2,2,5) = 1.0d0
   oril(3,3,5) = 1.0d0
   oril(1,1,6) = 1.0d0
   oril(2,2,6) = 1.0d0
   oril(3,3,6) = 1.0d0
   oril(1,1,7) = 1.0d0
   oril(2,2,7) = 1.0d0
   oril(3,3,7) = 1.0d0
   oril(1,1,8) = 1.0d0
   oril(2,2,8) = 1.0d0
   oril(3,3,8) = 1.0d0

end subroutine SetDiamond

!************************************************************************
!*                                                                      *
!*     SetH2O                                                           *
!*                                                                      *
!************************************************************************

! ... generate a cubic lattice, im3m (ice viii)

subroutine SetH2O(nlp, rol, oril)

   implicit none

   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   real(8) :: s2

   nlp = 2
   rol(1:3,1) = 0.0d0
   rol(1:3,2) = 0.5d0
   s2 = 1.0d0/sqrt(2.0d0)
   oril(1,1,1) = s2
   oril(2,1,1) = s2
   oril(1,2,1) =-s2
   oril(2,2,1) = s2
   oril(3,3,1) = 1.0d0
   oril(1,1,2) = s2
   oril(2,1,2) =-s2
   oril(1,2,2) = s2
   oril(2,2,2) = s2
   oril(3,3,2) = 1.0d0

end subroutine SetH2O

!************************************************************************
!*                                                                      *
!*     SetN2                                                            *
!*                                                                      *
!************************************************************************

! ... generate a cubic lattice, pa3 (solid n2)

subroutine SetN2(nlp, rol, oril)

   implicit none

   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   real(8)    :: s2, s3, s6

   nlp = 4
   rol(1:3,1) = 0.0d0
   rol(1,2) = 0.0d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.5d0
   rol(1,3) = 0.5d0
   rol(2,3) = 0.0d0
   rol(3,3) = 0.5d0
   rol(1,4) = 0.5d0
   rol(2,4) = 0.5d0
   rol(3,4) = 0.0d0

   s2 = 1.0d0/sqrt(2.0d0)
   s3 = 1.0d0/sqrt(3.0d0)
   s6 = 1.0d0/sqrt(6.0d0)
   oril(1,1,1) = s2
   oril(2,1,1) =-s2
   oril(3,1,1) = 0.0d0
   oril(1,2,1) = s6
   oril(2,2,1) = s6
   oril(3,2,1) =-2.0d0*s6
   oril(1,3,1) = s3
   oril(2,3,1) = s3
   oril(3,3,1) = s3
   oril(1,1,2) = s2
   oril(2,1,2) = 0.0d0
   oril(3,1,2) =-s2
   oril(1,2,2) =-s6
   oril(2,2,2) =-2.0d0*s6
   oril(3,2,2) =-s6
   oril(1,3,2) =-s3
   oril(2,3,2) = s3
   oril(3,3,2) =-s3
   oril(1,1,3) = s2
   oril(2,1,3) =-s2
   oril(3,1,3) = 0.0d0
   oril(1,2,3) = s6
   oril(2,2,3) = s6
   oril(3,2,3) = 2.0d0*s6
   oril(1,3,3) =-s3
   oril(2,3,3) =-s3
   oril(3,3,3) = s3
   oril(1,1,4) = 0.0d0
   oril(2,1,4) = s2
   oril(3,3,4) =-s2
   oril(1,2,4) = 2.0d0*s6
   oril(2,2,4) = s6
   oril(3,2,4) = s6
   oril(1,3,4) = s3
   oril(2,3,4) =-s3
   oril(3,3,4) =-s3

end subroutine SetN2

!************************************************************************
!*                                                                      *
!*     SetBenzene                                                       *
!*                                                                      *
!************************************************************************

! ... generate a cubic lattice, pbca (solid benzene)

!     a molecular orthogonal frame is attached to a benzene molecule.
!     the benzene molecule is orientated in this frame occording to
!         x'-axis  along the vector origin-carbon no 1
!         y'-axis  along the vector origin-the middle of carbon 2 and 3
!         z'-axis  perpendicular to the plane of the benzene molecule
!
!     this routine
!     1) uses coordinates of carbon no 1,2, and 3 given in the fixed
!        coordinate system. thses carbons belongs to the benzene
!        molecule at the corner of the unit cell. then evaluates f(m),
!        g(m), and h(m)
!     2) generates coordinates of the three axes; oril which
!        describe the orientation of the four molecular frames in a
!        unit cell. the coordinates are given in the fixed frame.

!     oril    orientation of particle frame in the lab frame
!             at different lattcie points in a unit cell
!             oril(m,1,2) reflexion of oril(m,1,2) in the plane (1,0,0)
!             oril(m,1,3) reflexion of oril(m,1,3) in the plane (0,1,0)
!             oril(m,1,4) reflexion of oril(m,1,4) in the plane (0,0,1)
!             oril(m,2,1:4) and oril(m,3,1:4) are obtained in the same way.
!             ref: cox et al., proc. r. soc. a 247, 1(1958)

subroutine SetBenzene(nlp, rol, oril)

   implicit none

   integer(4), intent(out) :: nlp
   real(8),    intent(out) :: rol(3,4)
   real(8),    intent(out) :: oril(3,3,4)
   real(8)    :: f(3), g(3), h(3)                            ! orientation of the benezene mol in the corner
   real(8)    :: x1 = -0.38340, y1 = 1.33980, z1 = -0.04010  ! coordinate of carbon 1
   real(8)    :: x2 = -0.97130, y2 = 0.46500, z2 =  0.87320  ! coordinate of carbon 2
   real(8)    :: x3 = -0.59140, y3 =-0.87110, z3 =  0.92160  ! coordinate of carbon 3
   integer(4) :: ip, mm

   nlp = 4
   rol(1,1) = 0.0d0
   rol(2,1) = 0.0d0
   rol(3,1) = 0.0d0
   rol(1,2) = 0.5d0
   rol(2,2) = 0.5d0
   rol(3,2) = 0.0d0
   rol(1,3) = 0.0d0
   rol(2,3) = 0.5d0
   rol(3,3) = 0.5d0
   rol(1,4) = 0.5d0
   rol(2,4) = 0.0d0
   rol(3,4) = 0.5d0

   f(1) = x1
   f(2) = y1
   f(3) = z1
   f = f/sqrt(sum(f(1:3)**2))

   g(1) = x2+x3
   g(2) = y2+y3
   g(3) = z2+z3
   g = g/sqrt(sum(g(1:3)**2))

   h(1) = f(2)*g(3)-f(3)*g(2)
   h(2) = f(3)*g(1)-f(1)*g(3)
   h(3) = f(1)*g(2)-f(2)*g(1)
   h = h/sqrt(sum(h(1:3)**2))

   do ip = 1, 4
      oril(1:3,1,ip) = f(1:3)
      oril(1:3,2,ip) = g(1:3)
      oril(1:3,3,ip) = h(1:3)
   end do
   do mm = 1, 3
      oril(mm,1,mm+1) = -f(mm)
      oril(mm,2,mm+1) = -g(mm)
      oril(mm,3,mm+1) = -h(mm)
   end do
!  write(*,*) ' -- from SetBenzene --'
!  do ip = 1, 4
!     write(*,*) ' molecule', ip, '       x         y       z'
!     write(*,*) oril(1:3,1,ip)
!     write(*,*) oril(1:3,2,ip)
!     write(*,*) oril(1:3,3,ip)
!  end do

end subroutine SetBenzene

!************************************************************************
!*                                                                      *
!*     SetRandom                                                        *
!*                                                                      *
!************************************************************************

! ... generate random positions and orientations

subroutine SetRandom(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt                   ! particle type

   character(40), parameter :: txroutine ='SetRandom'
   integer(4) :: ntry, itry, ip, iploc
   logical    :: lWarnHCOverlap

   if (nppt(ipt) == 0) return
   ntry = ntrydef

   do iploc = 1, nppt(ipt)                         ! loop over particles of type ipt
      ip = iploc-1+ipnpt(ipt)                      ! particle to be set

      do itry = 1, ntry                            ! loop over attempts to set the particle
         call SetPartPosRandom(ip)                 ! set particle position
         call SetPartOriRandom(iseed,ori(1,1,ip))  ! set particle orientation
         call SetAtomPos(ip,ip,.false.)            ! set atom positions
         if (luext) then
            if (txuext(ipt) == 'sphdielboundary') then!
              if (CheckHSDielBoundaryOverlap(ip)) cycle! under development
            end if
         end if
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
         lpset(ip) = .true.                        ! configuration accepted
         exit
      end do

      if (itry > ntry) then                         ! number of  attempts exceeds the maximal one ?
         if (master) write(uout,'(4(a,i5,2x))') 'particle type =', ipt,'iploc =', iploc, 'nppt(ipt) =', nppt(ipt)
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if

   end do

contains

!........................................................................

function CheckHSDielBoundaryOverlap(ip)  ! exclude HS--dielectric-boundary overlap
   use CoordinateModule
   implicit none
   character(40), parameter :: txroutine ='CheckHSDielBoundaryOverlap'
   integer(4), intent(in) :: ip
   real(8) :: r1
   logical :: CheckHSDielBoundaryOverlap
   if(na /= np) call Stop(txroutine, ' na /= np', uout)
   r1 = sqrt(ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2)
   if (abs(r1-boundaryrad) < radatset(iptpn(ip))) then
      CheckHSDielBoundaryOverlap =.true.
   else
      CheckHSDielBoundaryOverlap =.false.
   end if
end function CheckHSDielBoundaryOverlap

end subroutine SetRandom

!************************************************************************
!*                                                                      *
!*     SetRandomFixOri                                                  *
!*                                                                      *
!************************************************************************

! ... generate random positions and fixed orientations

subroutine SetRandomFixOri(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt                   ! particle type
   integer(4), parameter :: ntry = 100
   integer(4) :: itry, ip, iploc
   logical    :: lWarnHCOverlap

   if (nppt(ipt) == 0) return

   do iploc = 1, nppt(ipt)                         ! loop over particles of type ipt
      ip = iploc-1+ipnpt(ipt)                      ! particle to be set

      do itry = 1, ntry                            ! loop over attempts to set the particle
         call SetPartPosRandom(ip)                 ! set particle position
         call SetPartOriLab(ori(1,1,ip))           ! set particle orientation
         call SetAtomPos(ip,ip,.false.)            ! set atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
         lpset(ip) = .true.                        ! configuration accepted
         exit
      end do

      if (itry > ntry) then                         ! number of  attempts exceeds the maximal one ?
         if (master) write(uout,'(4(a,i5,2x))') 'particle type =', ipt,'iploc =', iploc, 'nppt(ipt) =', nppt(ipt)
         call Stop('SetRandom', 'random configuration failed, itry > ntry', uout)
      end if

   end do

end subroutine SetRandomFixOri

!************************************************************************
!*                                                                      *
!*     SetChainLine                                                     *
!*                                                                      *
!************************************************************************

! ... generate a linear configuration for chain particles (only for one chain)
!     along the x-axis, lab orientation

subroutine SetChainLine(iptset)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainLine'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, jp, jseg
   real(8)    :: bondloc
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   if (lclink) call Stop(txroutine, 'lclink is true', uout)
   ntry = ntrydef

   if (.not.first) return                                   ! should be called only once
   first =.false.

   do ic = 1, nc                                            ! loop over chains
      ict = ictcn(ic)                                       ! chain type
      bondloc = bondscl(ict)*bond(ict)%eq                   ! bond length to be used
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)                             ! particle to be set

         do itry = 1, ntry                                  ! loop over attempts to set the particle
            if (iseg == 1) then                             ! a first segment
               ro(1,ip) = -half*(npct(ict)-1)*bondloc
               ro(2:3,ip) = zero
            else                                            ! a remaining segment
               jp = ipnsegcn(iseg-1,ic)
               ro(1,ip) = ro(1,jp) + bondloc
               ro(2:3,ip) = ro(2:3,jp)
            end if
            if (CheckPartOutsideBox(ip)) cycle              ! check if particle is outside the box
            if (CheckTooFoldedChain(ip, bondloc)) cycle     ! check if a too folded chain
            call SetPartOriLoc(ori(1,1,ip))                 ! set particle orientation
            call SetAtomPos(ip,ip,.false.)                  ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                               ! configuration accepted
            exit
         end do

         if (itry > ntry) then                              ! number of attempts exceeds the maximal one ?
            if (master) write(uout,*)
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            if (master) write(uout,'(i4,3g12.5)') (ipnsegcn(jseg,ic), ro(:,ipnsegcn(jseg,ic)), jseg = 1, iseg)
            call Stop(txroutine, 'linear chain configuration failed, itry > ntry', uout)
         end if

      end do

   end do

contains

!........................................................................

subroutine SetPartOriLoc(ori)  ! x' in -z direction, y' in y direction, z' in x direction
   real(8),    intent(out)    :: ori(3,3)  ! oritentation matrix
   ori(1,1) = 0.0d0
   ori(2,1) = 0.0d0
   ori(3,1) =-1.0d0
   ori(1,2) = 0.0d0
   ori(2,2) = 1.0d0
   ori(3,2) = 0.0d0
   ori(1,3) = 1.0d0
   ori(2,3) = 0.0d0
   ori(3,3) = 0.0d0
end subroutine SetPartOriLoc

end subroutine SetChainLine

!************************************************************************
!*                                                                      *
!*     SetChainCircle                                                   *
!*                                                                      *
!************************************************************************

! ... generate a circular configuration for chain particles (only for one chain)
!     in the xy-plane

subroutine SetChainCircle(iptset)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainCircle'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, jseg
   real(8)    :: bondloc
   real(8)    :: radcir, angle0
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, lWarnHCOverlap

   if (lclink) call Stop(txroutine, 'lclink is true', uout)
   ntry = 1
   write(*,*) txroutine

   if (.not.first) return                                   ! should be called only once
   first =.false.

   if (nc /= 1) call stop(txroutine, 'nc =/1', uout)

   do ic = 1, nc                                            ! loop over chains
      ict = ictcn(ic)                                       ! chain type
      bondloc = bondscl(ict)*bond(ict)%eq                   ! bond length to be used
      angle0 = TwoPi/npct(ict)
      radcir = half*bondloc/sin(half*angle0)
      radcir = half*bondloc*(one+one/tan(half*angle0))
       if (lsuperball) call SuperballSub(radcir)             ! adjust radcir for superballs
!      write(*,'(a,f10.3)') 'bondscl(ict) ',bondscl(ict)
!      write(*,'(a,f10.3)') 'bond(ict)%eq ',bond(ict)%eq
!      write(*,'(a,f10.3)') 'bondloc ',bondloc
!      write(*,'(a,i10)') 'npart',npct(ict)
!      write(*,'(a,f10.3)') 'radcir ',radcir
!      write(*,'(a,f10.3)') 'angle0 ',angle0
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)                             ! particle to be set
         do itry = 1, ntry                                  ! loop over attempts to set the particle
            ro(1,ip) = radcir*cos((iseg-1)*angle0)
            ro(2,ip) = radcir*sin((iseg-1)*angle0)
            ro(3,ip) = zero
            if (CheckPartOutsideBox(ip)) cycle              ! check if particle is outside the box
            call SetOriCircle(ro(1,ip),ori(1,1,ip))         ! set orientation, z'-axis in the tangient
            call SetAtomPos(ip,ip,.false.)                  ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                               ! configuration accepted
            exit
         end do

         if (itry > ntry) then                              ! number of attempts exceeds the maximal one ?
            if (master) write(uout,*)
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            if (master) write(uout,'(i4,3g12.5)') (ipnsegcn(jseg,ic), ro(:,ipnsegcn(jseg,ic)), jseg = 1, iseg)
            call Stop(txroutine, 'circular chain configuration failed, itry > ntry', uout)
         end if

      end do
   end do

contains

!........................................................................

subroutine SuperballSub(radcir)            ! calculate radius of the circle for superballs
   real(8), intent(inout) :: radcir
   integer(4) :: iter
   real(8) :: r21(3), r2, of, SuperballOverlapOF

   write(*,*) 'SuperballSub: qsuperball',qsuperball
   write(*,*) 'SuperballSub: radcir (orig)',radcir
   if (qsuperball < qsuperball_max_nr) then
      do iter = 1, 2
         if (iter == 2) radcir = radcir*(1d0+1d-10)/sqrt(of)     ! scale with 1/sqrt(of) to get hc-contact start configuration
         do iseg = 1, 2
            ip = ipnsegcn(iseg,ic)
            ro(1,ip) = radcir*cos((iseg-1)*angle0)
            ro(2,ip) = radcir*sin((iseg-1)*angle0)
            ro(3,ip) = zero
            call SetOriCircle(ro(1,ip),ori(1,1,ip))
      !       write(*,'(a,3(3f10.3,5x))') 'setoricircle',  ro(1:3,ip), ori(1:3,2,ip), ori(1:3,3,ip)
         end do
         r21(1:3) = ro(1:3,ipnsegcn(1,ic))-ro(1:3,ipnsegcn(2,ic))
         r2 = r21(1)**2 + r21(2)**2 + r21(3)**2
         of = SuperballOverlapOF(r2,r21,ori(1,1,ipnsegcn(1,ic)),ori(1,1,ipnsegcn(2,ic)))
         write(*,'(a,i5,2f15.10)') 'SuperballSub (qsuperball < qsuperball_max_nr): iter, radcir, of', iter, radcir, of
      end do
   else
      radcir = radcir*(1d0+1d-6)
   write(*,*) 'SuperballSub (qsuperball > qsuperball_max_nr): radcir',radcir
   end if
end subroutine SuperballSub

!........................................................................

subroutine SetOriCircle(ro,ori)
   real(8), intent(inout) :: ro(3)
   real(8), intent(inout) :: ori(3,3)
   ori(:,:) = zero
   ori(3,1) = -one                                          ! x'-axis
   ori(1:3,2) = ro(1:3)                                     ! y'-axis
   ori(1,3) = ori(2,1)*ori(3,2) - ori(3,1)*ori(2,2)         ! z'-axis
   ori(2,3) = ori(3,1)*ori(1,2) - ori(1,1)*ori(3,2)
   ori(3,3) = ori(1,1)*ori(2,2) - ori(2,1)*ori(1,2)
   ori(1:3,1) = ori(1:3,1)/sqrt(ori(1,1)**2+ori(2,1)**2+ori(3,1)**2)   ! normalize
   ori(1:3,2) = ori(1:3,2)/sqrt(ori(1,2)**2+ori(2,2)**2+ori(3,2)**2)
   ori(1:3,3) = ori(1:3,3)/sqrt(ori(1,3)**2+ori(2,3)**2+ori(3,3)**2)
end subroutine SetOriCircle

end subroutine SetChainCircle

!************************************************************************
!*                                                                      *
!*     SetChainRandom                                                   *
!*                                                                      *
!************************************************************************

! ... generate random positions and orientations for chain particles

subroutine SetChainRandom(iptset)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainRandom'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, ipt, jp
   real(8)    :: bondloc
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   ntry = ntrydef

   if (.not.first) return                                    ! should be called only once
   first =.false.

   do ic = 1, nc                                            ! loop over chains
      ict = ictcn(ic)                                       ! chain type
      bondloc = bondscl(ict)*bond(ict)%eq                   ! bond length to be used
      if (lclink) then
         ip = ipnsegcn(1,ic)
         ipt = iptpn(ip)
         if ((ipt == iptnode) .or. (ict == ictstrand)) cycle ! exclude chains which are nodes or strands
         if(ihnpn(ip) /= 0) cycle                           ! exclude chains which are in hierarchical structure
      ! MORE
      end if
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)                             ! particle to be set

         do itry = 1, ntry                                  ! loop over attempts to set the particle
            if (iseg == 1) then                             ! a first segment
               call SetPartPosRandom(ip)
            else                                            ! a remaining segment
               jp = ipnsegcn(iseg-1,ic)
               call SetPartPosRandomN(ip, jp, bondloc)
            end if
            if (CheckPartOutsideBox(ip)) cycle              ! check if particle is outside the box
            if (CheckTooFoldedChain(ip, bondloc)) cycle     ! check if a too folded chain
            call SetPartOriRandom(iseed,ori(1,1,ip))        ! set random particle orientation
            call SetAtomPos(ip,ip,.false.)                  ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                               ! configuration accepted
            exit
         end do

         if (itry > ntry) then                              ! number of attempts exceeds the maximal one ?
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
         end if

      end do

   end do

end subroutine SetChainRandom

!************************************************************************
!*                                                                      *
!*     SetChainRandomIntOri                                             *
!*                                                                      *
!************************************************************************

! ... generate random positions and internally fixed orientations for chain particles

subroutine SetChainRandomIntOri(iptset)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: iptset                ! type of particle type in chain

   character(40), parameter :: txroutine ='SetChainRandomIntOri'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, ipt, jp, ipprev
   real(8)    :: bondloc, r21(3), r23(3)
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   ntry = ntrydef

   lfixedori = .true.

   if (.not.first) return                                   ! should be called only once
   first =.false.

   do ic = 1, nc                                            ! loop over chains
      ict = ictcn(ic)                                       ! chain type
      bondloc = bondscl(ict)*bond(ict)%eq                   ! bond length to be used
      if (lclink) then
         ip = ipnsegcn(1,ic)
         ipt = iptpn(ip)
         if ((ipt == iptnode) .or. (ict == ictstrand)) cycle ! exclude chains which are nodes or strands
         if(ihnpn(ip) /= 0) cycle                           ! exclude chains which are in hierarchical structure
      ! MORE
      end if
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)                             ! particle to be set

         do itry = 1, ntry                                  ! loop over attempts to set the particle
            if (iseg == 1) then                             ! a first segment
               call SetPartPosRandom(ip)
            else                                            ! a remaining segment
               jp = ipnsegcn(iseg-1,ic)
               call SetPartPosRandomN(ip, jp, bondloc)
            end if

            if (CheckPartOutsideBox(ip)) cycle              ! check if particle is outside the box

            if (CheckTooFoldedChain(ip, bondloc)) cycle     ! check if a too folded chain

!jvr Set atom orientations

            if (iseg == 1) then
               call SetPartOriLab(ori(1,1,ip))                        ! set first particle in frame orientation
               call SetAtomPos(ip,ip,.false.)
            else if (iseg == 2) then
               ipprev = ipnsegcn(iseg-1,ic)

               r21(1) = ro(1,ipprev)-ro(1,ip)
               r21(2) = ro(2,ipprev)-ro(2,ip)
               r21(3) = ro(3,ipprev)-ro(3,ip)
               call PBC(r21(1), r21(2), r21(3))
               call SetPartOriBond('end',-r21,-r21, ori(1,1,ipprev))  ! set first particle to align with r12 vector			
               call SetPartOriBond('end', r21, r21, ori(1,1,ip))      ! set second paricle to aligne with r12 vector

               call SetAtomPos(ipprev,ip,.false.)
            else
               ipprev = ipnsegcn(iseg-1,ic)
               r21(1) = ro(1,ipnsegcn(iseg-2,ic))-ro(1,ipprev)
               r21(2) = ro(2,ipnsegcn(iseg-2,ic))-ro(2,ipprev)
               r21(3) = ro(3,ipnsegcn(iseg-2,ic))-ro(3,ipprev)
               r23(1) = ro(1,ip)-ro(1,ipprev)
               r23(2) = ro(2,ip)-ro(2,ipprev)
               r23(3) = ro(3,ip)-ro(3,ipprev)
               call PBC(r21(1), r21(2), r21(3))
               call PBC(r23(1), r23(2), r23(3))
               call SetPartOriBond('mid', r21, r23, ori(1,1,ipprev))  ! set one but last particle to orient with chain bonds
               call SetPartOriBond('end',-r23,-r23, ori(1,1,ip))      ! set last paricle to aligne with r23 vector

               call SetAtomPos(ipprev,ip,.false.)
            end if

            if (lWarnHCOverlap(ip, radatset, .true.)) cycle           ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                                         ! configuration accepted
            exit
         end do

         if (itry > ntry) then                                        ! number of attempts exceeds the maximal one ?
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
         end if

      end do

   end do

end subroutine SetChainRandomIntOri

!************************************************************************
!*                                                                      *
!*     SetSphBrush                                                      *
!*                                                                      *
!************************************************************************

! ... generate a configuration for a spherical brush

!     special considerations:
!        lattice arrangement of grafted particles
!        the particle at which the grafted chains are attached should be set first

subroutine SetSphBrush(txfirstseg)

   use CoordinateModule
   implicit none

   character(*), intent(in) :: txfirstseg

   character(40), parameter :: txroutine ='SetSphBrush'
   integer(4), parameter :: nlatticepoint = 20
   integer(4), save :: jptsph = 1        ! type of particle at which the grafted chains are attached
   integer(4) :: jpsph, iploc
   integer(4) :: ntry, itry, iseg, ic, ict, ip, ipt, jp, ictemp
   real(8)    :: bondloc
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap

   if (ipnsegcn(1,1) < 2) call Stop(txroutine, 'central particle not set', uout)
   ntry = ntrydef

   if (.not.first) return                                    ! should be called only once
   first =.false.

   iseg = 0
   do
      iseg = iseg + 1

      do ic = 1, nc                                          ! loop over chains
         ict = ictcn(ic)                                     ! chain type
         ictemp = ic-1
         jpsph = 1+mod(ictemp,nppt(jptsph))                  ! particle at which the chain should be grafted
         bondloc = bondscl(ict)*bond(ict)%eq                 ! bond length to be used

         ip = ipnsegcn(iseg,ic)                              ! particle to be set
         ipt = iptpn(ip)

         do itry = 1, ntry                                   ! loop over attempts to set the particle
            if (iseg == 1) then                              ! a first segment
               call SetFirstSegment
            else                                             ! a remaining segment
               jp = ipnsegcn(iseg-1,ic)
               call SetPartPosRandomN(ip, jp, bondloc)
            end if
            if (CheckPartOutsideBox(ip)) cycle               ! check if particle outside box
            if (CheckTooFoldedChain(ip,bondloc)) cycle       ! check if a too folded chain
            call SetPartOriRandom(iseed,ori(1,1,ip))         ! set random particle orientation
            call SetAtomPos(ip,ip,.false.)                   ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle  ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                                ! configuration accepted
            exit
         end do

         if (itry > ntry) then                               ! number of attempts exceeds the maximal one ?
            if (master) write(uout,'(4(a,i5,2x))') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
         end if

      end do

      if ((ic >= nc) .and. (iseg == npct(ict))) exit
   end do

contains

!........................................................................

subroutine SetFirstSegment
   character(40), parameter :: txroutine ='SetFirstSegment'
   real(8)    :: dx, dy, dz, norm, Random, distance
   real(8)    :: rosph(1:3,1:nlatticepoint) = RESHAPE( [ -0.1753,   0.8443,   0.5052,    &
                                                         -0.7715,  -0.4709,  -0.4275,    &
                                                         -0.4656,  -0.8756,   0.1268,    &
                                                          0.2180,  -0.0537,  -0.9743,    &
                                                          0.4813,   0.4430,   0.7563,    &
                                                          0.2930,  -0.4969,   0.8166,    &
                                                         -0.6226,  -0.4433,   0.6443,    &
                                                         -0.5361,   0.8195,  -0.2015,    &
                                                          0.2623,  -0.9568,   0.1239,    &
                                                         -0.9814,   0.1484,  -0.1210,    &
                                                          0.9439,   0.2543,   0.2103,    &
                                                          0.4272,   0.9040,  -0.0087,    &
                                                          0.0694,   0.7101,  -0.7006,    &
                                                         -0.7800,   0.3004,   0.5490,    &
                                                         -0.1457,   0.0874,   0.9852,    &
                                                          0.7021,  -0.4912,  -0.5150,    &
                                                         -0.0645,  -0.7476,  -0.6610,    &
                                                         -0.5142,   0.1592,  -0.8426,    &
                                                          0.8492,  -0.4574,   0.2633,    &
                                                          0.7860,   0.3167,  -0.5307 ],  &
                                                        [ 3, 20 ])
!                                                [ 3, nlatticepoint ])
   distance = (One+1d-10)*(radatset(ipt)+radatset(iptpn(jpsph)))
   if (txfirstseg == 'lattice') then
      iploc = 1 + (ic-1)/nppt(jptsph)
      if (iploc > nlatticepoint) call Stop(txroutine,'iploc > nlatticepoint',uout)
      norm = distance/sqrt(rosph(1,iploc)**2+rosph(2,iploc)**2+rosph(3,iploc)**2)
      ro(1,ip) = ro(1,jpsph) + rosph(1,iploc)*norm
      ro(2,ip) = ro(2,jpsph) + rosph(2,iploc)*norm
      ro(3,ip) = ro(3,jpsph) + rosph(3,iploc)*norm
   else if (txfirstseg == 'random') then
      dx = Random(iseed)-Half
      dy = Random(iseed)-Half
      dz = Random(iseed)-Half
      norm = distance/sqrt(dx**2+dy**2+dz**2)
      ro(1,ip) = ro(1,jpsph) + dx*norm
      ro(2,ip) = ro(2,jpsph) + dy*norm
      ro(3,ip) = ro(3,jpsph) + dz*norm
   else
      call stop(txroutine, 'error in txfirstseg',uout)
   end if
end subroutine SetFirstSegment

end subroutine SetSphBrush

!************************************************************************
!*                                                                      *
!*     SetPlanarBrush                                                   *
!*                                                                      *
!************************************************************************

! ... generate a configuration for a planar brush

!     special considerations:
!        particles belonging to grafted chains should be set first
!        txbc == 'xy'

subroutine SetPlanarBrush(txfirstseg)

   use CoordinateModule
   implicit none

   character(*), intent(in) :: txfirstseg            ! 'random' random distribution of first segment

   character(40), parameter :: txroutine ='SetPlanarBrush'
   integer(4) :: ntry, itry, iseg, ic, ict, ip, ipt, jp
   real(8)    :: bondloc
   logical    :: first =.true.
   logical    :: CheckTooFoldedChain, lWarnHCOverlap

   if (txbc /= 'xy') call Stop(txroutine, 'txbc /= ''xy''', uout)
   if (nct > 1) call Stop(txroutine, 'nct > 1', uout)
   ntry = ntrydef

   if (.not.first) return                              ! should be called only once
   first =.false.

   iseg = 0
   do
      iseg = iseg + 1

      do ic = 1, nc                                    ! loop over chains
         ict = ictcn(ic)                               ! chain type
         bondloc = bondscl(ict)*bond(ict)%eq           ! bond length to be used

         ip = ipnsegcn(iseg,ic)                        ! particle to be set
         ipt = iptpn(ip)

         do itry = 1, ntry                             ! loop over attempts to set the particle
            if (iseg == 1) then                        ! a first segment (random position on the surface)
               call SetFirstSegment
            else                                       ! a remaining segment (straight chain i z-dir)
               jp = ipnsegcn(iseg-1,ic)
               ro(1,ip) = ro(1,jp)
               ro(2,ip) = ro(2,jp)
               ro(3,ip) = ro(3,jp) + bondloc
            end if
            if (CheckTooFoldedChain(ip,bondloc)) cycle ! check if a too folded chain
            call SetPartOriRandom(iseed,ori(1,1,ip))   ! set random particle orientation
            call SetAtomPos(ip,ip,.false.)             ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
            lpset(ip) =.true.                          ! configuration accepted
            exit
         end do

         if (itry > ntry) then                         ! number of attempts exceeds the maximal one ?
            if (master) write(uout,'(4(a,i5),1x)') 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
            call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
         end if

      end do

      if ((ic >= nc) .and. (iseg == npct(ict))) exit
   end do

contains

!........................................................................

subroutine SetFirstSegment
   character(40), parameter :: txroutine ='SetFirstSegment'
   real(8) :: Random
   if (txfirstseg == 'random') then
      ro(1,ip) = 0.99d0*rclow(1,ipt)+(rcupp(1,ipt)-rclow(1,ipt))*Random(iseed)*0.99d0
      ro(2,ip) = 0.99d0*rclow(2,ipt)+(rcupp(2,ipt)-rclow(2,ipt))*Random(iseed)*0.99d0
      ro(3,ip) = rclow(3,ipt)
   else
      call stop(txroutine ,'error in txfirstseg',uout)
   end if
end subroutine SetFirstSegment

end subroutine SetPlanarBrush

!************************************************************************
!*                                                                      *
!*     SetHierarchical                                                  *
!*                                                                      *
!************************************************************************

! ... generate a configuration for hierarchical polymers

subroutine SetHierarchical(txfirstseg)

   use CoordinateModule
   implicit none

   character(*), intent(in) :: txfirstseg         ! 'random', 'lattice'

   character(40), parameter :: txroutine ='SetHierarchical'
   real(8)    :: bondloc
   integer(4) :: ntry, itry, iseg, ic, ict, ip, ipt, jp, igen
   logical    :: first =.true.
   logical    :: CheckPartOutsideBox, CheckTooFoldedChain, lWarnHCOverlap
   integer(4) :: iattempt

   if (.not.first) return                                           ! should be called only once
   first =.false.
   ntry = ntrydef

   attempt: do iattempt = 1, ntry                                      ! attempt multiple times to set the structure
      !initialize all relevant lpset to .false.
      do igen = 0, ngen                                                ! loop over generations
         ict = ictgen(igen)                                            ! chain type
         do ic = icnct(ict), icnct(ict)+ncct(ict)-1                    ! loop over chains
            do iseg = 1, npct(ict)                                     ! loop over segments
               ip = ipnsegcn(iseg,ic)                                  ! particle to be set
               lpset(ip) = .false.
            end do
         end do
      end do

      do igen = 0, ngen                                                ! loop over generations
         ict = ictgen(igen)                                            ! chain type
         bondloc = bondscl(ict)*bond(ict)%eq                           ! bond length to be used
         do ic = icnct(ict), icnct(ict)+ncct(ict)-1                    ! loop over chains
            do iseg = 1, npct(ict)                                     ! loop over segments
               ip = ipnsegcn(iseg,ic)                                  ! particle to be set
               ipt = iptpn(ip)                                         ! particle type (needed for SetFirstSegment)

               do itry = 1, ntry                                       ! loop over attempts to set the particle
                  if (iseg == 1) then                                  ! a first segment?
                     if (igen == 0) then                               ! zeroth generation?
                        call SetFirstSegment                           ! set first particle of a component
                     else                                              ! generation 1 or higher
                        jp = bondcl(1,ip)                              ! particle to which ip is crosslinked to
                        call SetPartPosRandomN(ip, jp, bondloc)
                     end if
                  else                                                 ! not a first segment
                     jp = ipnsegcn(iseg-1,ic)                          ! id of particle to which jp will be a neigbor
                     call SetPartPosRandomN(ip, jp, bondloc)
                  end if
                  if (CheckPartOutsideBox(ip)) cycle                   ! check if particle is outside the box
                  if (CheckTooFoldedChain(ip,bondloc)) cycle           ! check if a too folded chain
                  call SetPartOriRandom(iseed,ori(1,1,ip))             ! set random particle orientation
                  call SetAtomPos(ip,ip,.false.)                       ! set atom positions
                  if (lWarnHCOverlap(ip, radatset, .true.)) cycle      ! check if atom-atom hard-core overlap
                  lpset(ip) =.true.                                    ! configuration accepted
                  exit
               end do

               if (itry > ntry) then                                   ! number of attempts exceeds the maximal one?
                  cycle attempt                                        ! have another attempt at setting the structure
               end if
            end do
         end do
      end do
      exit                                                             ! if the attempt was successfull: accept configuration
   end do attempt

   if (iattempt > ntry) then                                   ! number of attempts exceeds the maximal one?
      if (master) then
         write(uout,'(5(a,i5,2x))') 'igen =', igen, 'ic =', ic, 'segment =', iseg, 'npct(ict) =', npct(ict)
         write(*,'(i5,3f10.3)') (jp, ro(1:3,jp), jp = 1, ip-1)
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if
   end if

contains

!........................................................................

subroutine SetFirstSegment                                           ! set first segment
   character(40), parameter :: txroutine ='SetFirstSegment'
   integer(4) :: ic = 0
   integer(4) :: ix, iy, iz
   real(8) :: rlen(3), rorigin(3)

   if (txfirstseg == 'lattice') then                                 ! pc lattice
      ic = ic + 1
      if (ic == 1) then
         if (product(nucell(1:3,ipt)) < ncct(ictgen(0))) call stop(txroutine, 'nucell too small', uout)
         rlen(1:3) = (rcupp(1:3,ipt)-rclow(1:3,ipt))/nucell(1:3,ipt) ! prepare rlen
         rorigin(1:3) = Half*(rclow(1:3,ipt)+rcupp(1:3,ipt))         ! prepare rorigin
      end if
      iz =  mod(ic-1,nucell(3,ipt))                                  ! inner loop: z-direction
      iy =  mod((ic-1)/nucell(3,ipt), nucell(2,ipt))
      ix = (ic-1)/(nucell(2,ipt)*nucell(3,ipt))                      ! outer loop: x-direction
      ro(1,ip) = ((ix-Half*nucell(1,ipt)+roshift(1,ipt)) * rlen(1)) + rorigin(1)
      ro(2,ip) = ((iy-Half*nucell(2,ipt)+roshift(2,ipt)) * rlen(2)) + rorigin(2)
      ro(3,ip) = ((iz-Half*nucell(3,ipt)+roshift(3,ipt)) * rlen(3)) + rorigin(3)
   else if (txfirstseg == 'random') then                             ! random
      call SetPartPosRandom(ip)
   else
      call stop(txroutine, 'error in txfirstseg',uout)
   end if
end subroutine SetFirstSegment

end subroutine SetHierarchical

!************************************************************************
!*                                                                      *
!*     SetPeriodicNetwork                                               *
!*                                                                      *
!************************************************************************

! ... generate a periodic network

!  unit cell: diamond-like containing 8 nodes and 16 strands
!  input variables: nucell, rcupp, rclow, iptnode, ictstrand
!  boundary condition: txbc == 'xyz'

subroutine SetPeriodicNetwork(ipt)

   use CoordinateModule
   implicit none

   character(40), parameter :: txroutine ='SetPeriodicNetwork'

   integer(4), intent(in)  :: ipt        ! particle type

   real(8), parameter :: rol(1:3,1:4) = RESHAPE([ 0.25, 0.25, 0.25,    &
                                                  0.25, 0.75, 0.75,    &
                                                  0.75, 0.25, 0.75,    &
                                                  0.75, 0.75, 0.25 ],  &
                                                [ 3, 4 ])
   real(8), parameter :: signgel(1:3,1:4) = RESHAPE([ -One, -One, -One,   &
                                                      -One, +One, +One,   &
                                                      +One, -One, +One,   &
                                                      +One, +One, -One ], &
                                                [ 3, 4 ])

   integer(4) :: ip, ic, idir, inode, iseg, ix, iy, iz, nnode, nstrand, nclmade
   real(8)    :: rlen(3), rorigin(3), nnsep(0:3), nnsepreq
   real(8), allocatable :: rbond(:,:), r2bond(:)
   integer(4) :: npclend                  ! number of particles that ends a crosslink
   integer(4), allocatable :: ipclbeg(:)  ! id of particles that begin crosslinks
   integer(4), allocatable :: ipclend(:)  ! id of particles that ends crosslinks

   external SetDiamond

   if (ipt /= iptnode) return
   if (txbc /= 'xyz') call Stop(txroutine, 'txbc /= ''xyz''', uout)

! ... number of nodes and strands

   nnode = nppt(iptnode)
   nstrand = ncct(ictstrand)

! ... check requirements and consistencies

   if (nnode /= 8*product(nucell(1:3,iptnode))) &
       call Stop (txroutine, 'number of nodes is not consistent with number of unit cells', uout)
   if (nstrand /= 16*product(nucell(1:3,iptnode))) &
       call Stop (txroutine, 'number of strands is not consistent with number of unit cells', uout)

! ... allocate memory

    allocate(ipclbeg(nnode), ipclend(2*nstrand), rbond(3,nct), r2bond(nct))
    ipclbeg = 0
    ipclend = 0
    rbond = 0.0E+00
    r2bond = 0.0E+00

! ... store all particles of type iptnode in ipclbeg

   ipclbeg(1:nnode) = [ (ip, ip = ipnpt(iptnode) , ipnpt(iptnode)+(nppt(iptnode)-1)) ]

! ... set node particles

   call SetLattice(iptnode, SetDiamond)

! ... set strand particles

   rlen(1:3) = (rcupp(1:3,iptnode)-rclow(1:3,iptnode))/nucell(1:3,iptnode)      ! prepare rlen
   rorigin(1:3) = Half*(rclow(1:3,iptnode)+rcupp(1:3,iptnode))                  ! prepare rorigin

   nnsep(1:3) = boxlen(1:3)/(Four*nucell(1:3,iptnode))                          ! node-node separation along the axes
   nnsep(0) = sqrt(sum(nnsep(1:3)**2))                                          ! node-node separation
   nnsepreq = sum(radatset(1:npt)*npptct(1:npt,ictstrand)) + radatset(iptnode)  ! required node-node separation
   rbond(1:3,ictstrand) = nnsep(1:3)/(npct(ictstrand)+1)                        ! initial bond length along the axes
   if (nnsepreq > nnsep(0)) call Stop(txroutine, 'too small node-node separation', uout) ! (no further hard-core overlap check)
   r2bond(ictstrand) = (nnsep(0)/(npct(ictstrand) + 1))**2

   npclend = 0                                   ! initialize npclend
   ic=icnct(ictstrand)-1                         ! set chain offset
   do idir = 1, 4                                ! loop over the four differnt directions of strands from the crosslinks
      do inode = 1, 4                            ! loop over the four nodes in the unit cell
         do iz = 0, nucell(3,iptnode)-1          ! loop over number of unit cells in z direction
            do iy = 0, nucell(2,iptnode)-1       ! loop over number of unit cells in y direction
               do ix = 0, nucell(1,iptnode)-1    ! loop over number of unit cells in x direction
                  ic=ic+1
                  do iseg = 1, npct(ictstrand)   ! loop over number of particles of chain ic
                     ip = ipnsegcn(iseg,ic)      ! particle to be set
                     ro(1,ip) = ((ix-Half*nucell(1,iptnode)+(rol(1,inode)+roshift(1,iptnode)))*rlen(1))+rorigin(1)
                     ro(2,ip) = ((iy-Half*nucell(2,iptnode)+(rol(2,inode)+roshift(2,iptnode)))*rlen(2))+rorigin(2)
                     ro(3,ip) = ((iz-Half*nucell(3,iptnode)+(rol(3,inode)+roshift(3,iptnode)))*rlen(3))+rorigin(3)
                     ro(1,ip) = ro(1,ip)+iseg*signgel(1,idir)*rbond(1,ictstrand)
                     ro(2,ip) = ro(2,ip)+iseg*signgel(2,idir)*rbond(2,ictstrand)
                     ro(3,ip) = ro(3,ip)+iseg*signgel(3,idir)*rbond(3,ictstrand)
                     call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
                     call SetPartOriRandom(iseed,ori(1,1,ip))             ! set random particle orientation
                     call SetAtomPos(ip,ip,.false.)                       ! set atom positions
                     lpset(ip) =.true.                                    ! position accepted
                     if ((iseg == 1) .or. (iseg == npct(ictstrand))) then ! store particle ip as clend
                         npclend = npclend + 1
                         if (npclend > 2*nstrand) call stop(txroutine, 'npclend > 2*nstrand', uout)
                         ipclend(npclend) = ip
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do

   call MakeCrossLink(nnode, npclend, ipclbeg, ipclend, r2bond, nclmade) ! make crosslinks between node and strand ends
!!!   write(*,*) txroutine, 'nclmade',nclmade
   ncl = ncl + nclmade
!!!   write(*,*) txroutine, 'ncl (updated)',ncl

   deallocate(ipclbeg, ipclend, rbond, r2bond)

end subroutine SetPeriodicNetwork

!************************************************************************
!*                                                                      *
!*     SetNetwork                                                       *
!*                                                                      *
!************************************************************************

! ... generate a nonperiodic network

!  unit cell: diamond-like containing 8 nodes and 16 strands
!  input variables: nnwt (nmlParticle), nnwnwt (nmlParticle), ncctnwt (nmlParticle), rnwt, iptclnwt (nmlParticle), txoriginnwt
!  boundary condition: txbc == 'xyz' or txbc == 'sph'

subroutine SetNetwork(ipt)

   use CoordinateModule
   implicit none

   character(40), parameter :: txroutine ='SetNetwork'

   integer(4), intent(in)  :: ipt        ! particle type

   integer(4)     :: iploc, jploc
   integer(4)     :: ic                  ! chain counter
   integer(4)     :: iseg                ! segment counter
   integer(4)     :: ip                  ! particle counter
   integer(4)     :: npclend             ! counter crosslinkable strand particles
   integer(4)     :: ict                 ! chain type of strand
   integer(4)     :: inwt                ! network type counter
   integer(4)     :: inw                 ! network counter
   integer(4)     :: ntry, itry          ! counter attempts to set gel
   integer(4)     :: nclmade             ! number of crossliks formed
   integer(4), allocatable :: ipclbeg(:) ! id of particles that begin crosslinks
   integer(4), allocatable :: ipclend(:) ! id of particles that ends crosslinks

   real(8)        :: rorigin(3)          ! coordinates of origin of network

! ... calculated by SetNetworkPos

   integer(4)     :: nnode, nstrand, nclnode
   real(8), allocatable :: ronode(:,:)
   real(8), allocatable :: rostrand(:,:)

   real(8)        :: Random
   logical        :: CheckPartOutsideBox, lWarnHCOverlap

   ntry = ntrydef

! ... determine network type

   if (count(ipt == iptclnwt(1:nnwt)) == 0) return   ! ipt is not a node
   do inwt = 1, nnwt                                 ! determine network type
      if (ipt == iptclnwt(inwt)) exit                ! ipt is a node of network type inwt
   end do

! ... determine chain type of strand

   do ict = 1, nct   ! Currently only one chain type per network is allowed 
      if (ncctnwt(ict,inwt) > 0) exit
   end do

! ... check condition

   if ((txoriginnwt(inwt) == 'origin') .and. (nnwnwt(inwt) > One)) call Stop(txroutine, 'txoriginnwt == origin .and. nnwnwt(inwt) > One', uout)

! ... allocate memory

   allocate(ronode(3,np_alloc), rostrand(3,np_alloc))
   ronode = 0.0E+00
   rostrand = 0.0E+00

! ... determine nnode, ronode, nstrand, rostrand, nclnode

   call SetNetworkPos(rnwt(inwt), bond(ict)%eq, npct(ict), nnode, ronode, nstrand, rostrand, nclnode)

   if(nnode*nnwnwt(inwt) /= nppt(ipt) .or. nstrand*nnwnwt(inwt) /= ncct(ict)) then ! when particle and chain number don't accord to the neccesary ones: stop excecution and write numbers
      call FileOpen(922, 'topo.dat', 'form/noread')      ! Cornelius Hofzumahaus 
      write(922,'(a10,i)') 'nchain'  , nstrand           !
      write(922,'(a10,i)') 'nnode'   , nnode             !
      write(922,'(a10,i)') 'npchain' , nstrand*npct(ict) !
      close(922)                                         ! 

      write(*,'(a50,i)') "failed to set network of network type:"     , inwt
      write(*,'(a50,i)') "required number of chains:"                   , nstrand*nnwnwt(inwt)
      write(*,'(a50,i)') "required number of node particles:"          , nnode*nnwnwt(inwt)
      write(*,'(a50,i)') "required (total) number of strand particles:" , nstrand*nnwnwt(inwt)*npct(ict)

      call Stop (txroutine, 'Please adjust nppt(node) and ncct(strand)', uout)
   end if

! ... allocate memory

   allocate(ipclbeg(nnode), ipclend(2*nstrand))
   ipclbeg = 0
   ipclend = 0

! ... set gels

   maxnbondcl(ipt) = nclnode

   do inw = 1, nnwnwt(inwt)! loop over all networks of network type inwt
try:  do itry = 1, ntry    ! loop over attempts to set the gel

!  ... set network origin

         if (txoriginnwt(inwt) == 'origin') then
            rorigin = Zero
         else
            if (lbcsph) then
               do
                  rorigin(1) = Two*sphrad*(Random(iseed) - Half)
                  rorigin(2) = Two*sphrad*(Random(iseed) - Half)
                  rorigin(3) = Two*sphrad*(Random(iseed) - Half)
                  if (sum(rorigin(1:3)**2) <= sphrad2) exit
               end do
            else if(lbcbox) then
               rorigin(1) = boxlen(1)*(Random(iseed) - Half)
               rorigin(2) = boxlen(2)*(Random(iseed) - Half)
               rorigin(3) = boxlen(3)*(Random(iseed) - Half)
               call PBC(rorigin(1), rorigin(2), rorigin(3))
            else
               call Stop (txroutine, 'txbc is not compatible with origin position of network', uout)
            end if
         end if

! ... set node particles

         do iploc = 1 , nnode
            ip = ipnpt(ipt) + iploc - 1 + (inw - 1)*nnode
            ro(1:3,ip) = rorigin(1:3) + ronode(1:3,iploc)
            call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
            call SetPartOriRandom(iseed,ori(1,1,ip))              ! set random particle orientation
            call SetAtomPos(ip, ip, .false.)                      ! set atom positions
            if (lWarnHCOverlap(ip, radatset, .true.)) cycle try   ! check for hard-core overlap
            if (CheckPartOutsideBox(ip)) cycle try                ! check that particle is inside box
            lpset(ip) =.true.                                     ! position accepted
            ipclbeg(iploc) = ip
         end do

! ...  set strand particles

         jploc = 0
         npclend = 0
         do ic = icnct(ict) + (inw - 1)*nstrand, icnct(ict) + (inw)*nstrand - 1
            do iseg = 1, npct(ict)
               jploc = jploc + 1
               ip = ipnsegcn(iseg,ic)
               ro(1:3,ip) = rorigin(1:3) + rostrand(1:3,jploc)
               call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
               call SetPartOriRandom(iseed,ori(1,1,ip))              ! set random particle orientation
               call SetAtomPos(ip, ip, .false.)                      ! set atom positions
               if (lWarnHCOverlap(ip, radatset, .true.)) cycle try   ! check for hard-core overlap
               if (CheckPartOutsideBox(ip)) cycle try                ! check that particle is inside bo
               lpset(ip) =.true.                                     ! position accepted
               if((iseg == 1) .or. (iseg == npct(ict))) then
                  npclend = npclend + 1
                  if (npclend > 2*nstrand) call stop(txroutine, 'npclend > 2*nstrand', uout)
                  ipclend(npclend) = ip
                  maxnbondcl(iptpn(ip)) = 1
               end if
            end do
         end do

! ... make crosslinks

!!!         write(*,*) txroutine
!!!         write(*,*) 'nnode, npclend',nnode,npclend
!!!         write(*,*) 'ipclbeg',ipclbeg(nnode)
!!!         write(*,*) 'ipclend',ipclend(npclend)

         call MakeCrossLink(nnode, npclend, ipclbeg, ipclend, (bond(1:nct)%eq)**2, nclmade)
!!!         write(*,*) txroutine, 'nclmade',nclmade
         ncl = ncl + nclmade
!!!         write(*,*) txroutine, 'ncl (updated)',ncl

         exit

      end do try

      if (itry > ntry) then                         ! number of  attempts exceeds the maximal one ?
         if (master) write(uout,'(2(a,i5,2x))') 'network type =', inwt,'network number =', inw
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if

   end do

   deallocate(ronode, rostrand)
   deallocate(ipclbeg, ipclend)

end subroutine SetNetwork

!........................................................................

subroutine SetNetworkPos(radgel, bondlen, npstrand, nnode, ronodeout, nstrand, rostrandout, nclnode)

   use CoordinateModule
   implicit none

   real(8),    intent(in)  :: radgel                  ! radius of gel
   real(8),    intent(in)  :: bondlen                 ! length of strand bond-length
   integer(4), intent(in)  :: npstrand                ! number of particles in strand
   integer(4), intent(out) :: nnode                   ! number of nodes
   real(8),    allocatable :: ronode(:,:)             ! coordinates of nodes
   real(8),    intent(out) :: ronodeout(3,np_alloc)   ! coordinates of nodes for output
   integer(4), intent(out) :: nstrand                 ! number of stand chains
   real(8),    allocatable :: rostrand(:,:)           ! coordinates of strand particles
   real(8),    intent(out) :: rostrandout(3,np_alloc) ! coordinates of strand particles for output
   integer(4), intent(out) :: nclnode                 ! number of crosslinks of a node (diamond: nclnode = 4)

   character(40), parameter :: txroutine ='SetNetworkPos'

!  ... network parameters

   real(8)     :: xbondlen       ! length of bond in strand along axis of unit cell
   real(8)     :: celllen        ! length of one unit cell
   integer(4)  :: ncell          ! number of unit cells in x,y,z direction
   real(8)     :: radgel2        ! square of network radius

! ... from SetDiamond

   integer(4)     :: nlp         ! number of particles in a unit cell
   real(8)        :: rol(3,8)    ! coordinates of the lattice point in a unit cell ranging from (0,0,0) to (1,1,1)
   real(8)        :: oril(3,3,8) ! orientations of particle frame at differen lattice points in a unit cell

   integer(4)     :: ilp, ix, iy, iz, idir, iseg, jp, iploc, jploc, npart_strand
   real(8)        :: InvFlt

   real(8), allocatable :: vhelp(:,:)

   real(8), parameter   :: ddelta = 1.0d-5

   real(8), parameter   :: signgel(1:3,1:4) = RESHAPE([ -One, -One, -One,    &   ! corresponds to symmetry of node particles on octahedral positions
                                                        -One, +One, +One,    &
                                                        +One, -One, +One,    &
                                                        +One, +One, -One ],  &
                                                                     [ 3, 4 ])
   real(8), parameter   :: signgel2(1:3,1:4) = RESHAPE([ +One, -One, -One,   &   ! corresponds to symmetry of node particles on corner and plane positions
                                                         +One, +One, +One,   &
                                                         -One, -One, +One,   &
                                                         -One, +One, -One ], &
                                                                     [ 3, 4 ])

   nclnode = 4                                     ! each node has 4 crosslinks
   call SetDiamond(nlp,rol,oril)                   ! get diamond unit cell informations
   radgel2 = radgel**2                             ! network radius squared
   celllen = Four*sqrt(Third)*bondlen*(npstrand+1) ! length of one cubic unit cell
   xbondlen = sqrt(Third) * bondlen                ! bond length projected on an external axis
   ncell = int(Two*radgel*InvFlt(celllen)) + 1     ! number of required unit cells
   if(modulo(ncell,2) == 1) ncell = ncell + 1      ! odd ncell -> even ncell in order to guarantee particle at ( 0 0 0 )

   nnode = 0
   nstrand = 0
   npart_strand = 0
   if(.not.allocated(rostrand)) then 
      allocate(rostrand(3,np_alloc), ronode(3,np_alloc))
      rostrand = 0.0E+00
      ronode = 0.0E+00
   end if
   rostrand = 0.00E+00
   ronode   = 0.00E+00
   ronodeout   = 0.00E+00
   rostrandout   = 0.00E+00

! ... loop over all lattice points and try to set particles at those

   do ilp = 1, nlp               ! loop over all lattice positions of one unit cell
      do iz = 0, ncell - 1       ! loop over number of unit cells in z-direction
         do iy = 0, ncell - 1    ! loop over number of unit cells in y-direction
            do ix = 0, ncell - 1 ! loop over number of unit cells in x-direction

! ... set node particles

               iploc = nnode+1
               if(iploc > size(ronode,2)) then
                  if (allocated(vhelp)) deallocate(vhelp)
                  allocate(vhelp(3,size(ronode,2)))
                  vhelp = 0.0E+00
                  vhelp = ronode
                  deallocate(ronode)
                  allocate(ronode(3,2*size(vhelp,2)))
                  ronode = 0.0E+00
                  ronode = vhelp
               end if
               ronode(1,iploc) = (ix - Half*ncell + rol(1,ilp))*celllen
               ronode(2,iploc) = (iy - Half*ncell + rol(2,ilp))*celllen
               ronode(3,iploc) = (iz - Half*ncell + rol(3,ilp))*celllen
               if(sum(ronode(1:3,iploc)**2) > radgel2) cycle    ! restrict to inside radgel
               nnode = nnode + 1                                ! update nnode

! ... set strand particles

               do idir = 1, 4                                   ! loop over all four directions from each node particle
     segment:     do iseg = 1, npstrand                         ! loop over all particles of each chain
                  jploc = npart_strand+1
                  if(jploc > size(rostrand,2)) then
                     if (allocated(vhelp)) deallocate(vhelp)
                     allocate(vhelp(3,size(rostrand,2)))
                     vhelp = 0.0E+00
                     vhelp = rostrand
                     deallocate(rostrand)
                     allocate(rostrand(3,2*size(vhelp,2)))
                     rostrand = 0.0E+00
                     rostrand = vhelp
                  end if 
                  if(ilp >= 5) then                             ! distinguish between node symmetry
                        if((idir == 1) .or. (idir == 3)) then   ! chains starting at node
                           rostrand(1:3,jploc) = ronode(1:3,iploc) + iseg*signgel(1:3,idir)*xbondlen
                        else                                    ! chains arriving at node
                           rostrand(1:3,jploc) = ronode(1:3,iploc) + (npstrand-iseg+1)*signgel(1:3,idir)*xbondlen
                        end if
                     else
                        if((idir == 2) .or. (idir == 4)) then   ! chains starting at node
                           rostrand(1:3,jploc) = ronode(1:3,iploc) + iseg*signgel2(1:3,idir)*xbondlen
                        else                                    ! chains arriving at node
                           rostrand(1:3,jploc) = ronode(1:3,iploc) + (npstrand-iseg+1)*signgel2(1:3,idir)*xbondlen
                        end if
                     end if
                     do jp = 1, jploc -1                        ! check if strand particle already exist
                        if(sum((rostrand(1:3,jp) - rostrand(1:3,jploc))**2) < ddelta ) exit segment
                     end do
                     npart_strand = npart_strand + 1            ! update npart_strand
                     if(iseg == npstrand) nstrand = nstrand + 1 ! update nstrand
                  end do segment
               end do
            end do
         end do
      end do
   end do

   if (allocated(vhelp)) deallocate(vhelp)

   rostrandout = rostrand(1:3,1:np_alloc)
   ronodeout   = ronode(1:3,1:np_alloc)
   
end subroutine SetNetworkPos

!***********************************************************************
!*                                                                     *
!*    SetCoreShell                                                     *
!*                                                                     *
!***********************************************************************

! ... generate a configuration with particles located in a spherical shell
!     bounded radially by radlimit(1:2) (Steffi Schneider)

subroutine SetCoreShell(ipt)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ipt            ! particle type

   character(40), parameter :: txroutine ='SetCoreShell'
   integer(4), parameter :: ntry = 1000000
   integer(4) :: itry, ip, iploc
   real(8)    :: r1
   logical    :: lWarnHCOverlap

! ... check consistencies

   if (txbc /= 'sph') call Stop(txroutine, 'txbc /= ''sph''', uout)
   if (rChargeIn > rChargeOut)  call stop(txroutine,'rChargeIn > rChargeOut',uout)
   if (rchargeOut > sphrad) call stop (txroutine,'rChargeOut > sphrad',uout)
   if (nppt(ipt) == 0) return

   do iploc = 1, nppt(ipt)                         ! loop over particles of type ipt
      ip = iploc-1+iptpn(ipt)                      ! particle to be set

      do itry = 1, ntry                            ! loop over attempts to set the particle
         call SetPartPosRandom(ip)                 ! trial particle position
         call SetPartOriRandom(iseed,ori(1,1,ip))  ! set particle orientation
         call SetAtomPos(ip,ip,.false.)            ! set atom positions
         if (lWarnHCOverlap(ip, radatset, .true.)) cycle ! check if atom-atom hard-core overlap
         r1=sqrt(ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2)
         if ((r1 < radlimit(1)) .or. (r1 > radlimit(2))) cycle ! test if atom-atom hard-core overlap
         lpset(ip) = .true.                        ! trial configuration accepted
         exit
      end do

      if (itry > ntry) then                        ! number of trial attempts exceeds the maximal one ?
         if (master) write(uout,'(4(1x,a,i5))') 'particle type =', ipt,'iploc =', iploc, 'nppt(ipt) =', nppt(ipt)
         call Stop(txroutine, 'random configuration failed, itry > ntry', uout)
      end if

   end do

 end subroutine SetCoreShell

!************************************************************************
!*                                                                      *
!*     SetPartPosRandom                                                 *
!*                                                                      *
!************************************************************************

! ... generate a random particle position

subroutine SetPartPosRandom(ip)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ip           ! particle

   real(8), parameter :: FAC = 0.99d0
   integer(4) :: ipt
   real(8) :: Random

   ipt = iptpn(ip)
   if (lbcsph) then
      do
         ro(1,ip) = Two*sphrad*(Random(iseed)-Half)
         ro(2,ip) = Two*sphrad*(Random(iseed)-Half)
         ro(3,ip) = Two*sphrad*(Random(iseed)-Half)
         if ((ro(1,ip) < rclow(1,ipt)) .or. (ro(1,ip) > rcupp(1,ipt))) cycle
         if ((ro(2,ip) < rclow(2,ipt)) .or. (ro(2,ip) > rcupp(2,ipt))) cycle
         if ((ro(3,ip) < rclow(3,ipt)) .or. (ro(3,ip) > rcupp(3,ipt))) cycle
         if (sum(ro(1:3,ip)**2) <= FAC*sphrad2) exit
      end do
   else if (lbccyl) then
      do
         ro(1,ip) = Two*cylrad*(Random(iseed)-Half)
         ro(2,ip) = Two*cylrad*(Random(iseed)-Half)
         if ((ro(1,ip) < rclow(1,ipt)) .or. (ro(1,ip) > rcupp(1,ipt))) cycle
         if ((ro(2,ip) < rclow(2,ipt)) .or. (ro(2,ip) > rcupp(2,ipt))) cycle
         if (sum(ro(1:2,ip)**2) <= FAC*cylrad2) exit
      end do
      ro(3,ip) = FAC*(rclow(3,ipt) + (rcupp(3,ipt)-rclow(3,ipt))*Random(iseed))
   else if (lbcell) then
      do
         ro(1,ip) = Two*ellrad(1)*(Random(iseed)-Half)
         ro(2,ip) = Two*ellrad(2)*(Random(iseed)-Half)
         ro(3,ip) = Two*ellrad(2)*(Random(iseed)-Half)
         if ((ro(1,ip) < rclow(1,ipt)) .or. (ro(1,ip) > rcupp(1,ipt))) cycle
         if ((ro(2,ip) < rclow(2,ipt)) .or. (ro(2,ip) > rcupp(2,ipt))) cycle
         if ((ro(3,ip) < rclow(3,ipt)) .or. (ro(3,ip) > rcupp(3,ipt))) cycle
         if (sum((ro(1:3,ip)*ellradi(1:3))**2) <= FAC) exit
      end do
   else
      ro(1,ip) = FAC*(rclow(1,ipt) + (rcupp(1,ipt)-rclow(1,ipt))*Random(iseed))
      ro(2,ip) = FAC*(rclow(2,ipt) + (rcupp(2,ipt)-rclow(2,ipt))*Random(iseed))
      ro(3,ip) = FAC*(rclow(3,ipt) + (rcupp(3,ipt)-rclow(3,ipt))*Random(iseed))
      call PBC(ro(1,ip),ro(2,ip),ro(3,ip))
   end if

end subroutine SetPartPosRandom

!************************************************************************
!*                                                                      *
!*     SetPartPosRandomN                                                *
!*                                                                      *
!************************************************************************

! ... generate a random position at a given distance to a previously set particle

subroutine SetPartPosRandomN(ip,jp,distance)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ip           ! particle to be set
   integer(4), intent(in) :: jp           ! neighbouring particle
   real(8),    intent(in) :: distance     ! distance to the particle

   real(8) :: dx, dy, dz, anorm, Random

   dx = Random(iseed) - half
   dy = Random(iseed) - half
   dz = Random(iseed) - half
   anorm = distance/sqrt(dx**2 + dy**2 + dz**2)
   ro(1,ip) = ro(1,jp) + dx*anorm
   ro(2,ip) = ro(2,jp) + dy*anorm
   ro(3,ip) = ro(3,jp) + dz*anorm

end subroutine SetPartPosRandomN

!************************************************************************
!*                                                                      *
!*     CheckPartOutsideBox                                              *
!*                                                                      *
!************************************************************************

! ... check if a particle is outside the box

logical function CheckPartOutsideBox(ip)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ip           ! particle

   integer(4) :: ipt
   logical :: nocrossing(3)

   CheckPartOutsideBox = .false.
   ipt = iptpn(ip)

! ... check if outside the box defined by rclow and rcupp

   if (lbcbox) then
      where (rclow(1:3,ipt) > -boxlen2(1:3) .or. rcupp(1:3,ipt) < boxlen2(1:3))
         nocrossing(1:3) =.true.
      elsewhere
         nocrossing(1:3) =.false.
      end where
   else if (lbcsph .or. lbccyl .or. lbcell) then
      nocrossing(1:3) =.true.
   else
      nocrossing(1:3) = .false.
   end if

   if (nocrossing(1) .and. (ro(1,ip) < rclow(1,ipt) .or. ro(1,ip) > rcupp(1,ipt))) CheckPartOutsideBox = .true.
   if (nocrossing(2) .and. (ro(2,ip) < rclow(2,ipt) .or. ro(2,ip) > rcupp(2,ipt))) CheckPartOutsideBox = .true.
   if (nocrossing(3) .and. (ro(3,ip) < rclow(3,ipt) .or. ro(3,ip) > rcupp(3,ipt))) CheckPartOutsideBox = .true.

! ... checks related to the standard simulation cell
                                                             ! ??? lbcrd and lbcto
   if (lbcbox) then
      if (txbc == 'xyz' .or. lbcrd .or. lbcto) then
         call PBC(ro(1,ip), ro(2,ip), ro(3,ip))
      else if (txbc == 'xy') then
         if (abs(ro(3,ip)) > boxlen2(3)) CheckPartOutsideBox = .true.
         call PBC(ro(1,ip), ro(2,ip), ro(3,ip))
      else if (txbc == 'z') then
         if (abs(ro(1,ip)) > boxlen2(1)) CheckPartOutsideBox = .true.
         if (abs(ro(2,ip)) > boxlen2(2)) CheckPartOutsideBox = .true.
         call PBC(ro(1,ip), ro(2,ip), ro(3,ip))
      end if
   else if (lbcsph) then
      if (sum(ro(1:3,ip)**2) >= sphrad2-0.1d0) CheckPartOutsideBox = .true.
   else if (lbccyl) then
      if (sum(ro(1:2,ip)**2) >= cylrad2-0.1d0) CheckPartOutsideBox = .true.
      if (abs(ro(3,ip)) > Half*cyllen-0.1d0) CheckPartOutsideBox = .true.
   else if (lbcell) then
      if (sum((ro(1:3,ip)*ellradi(1:3))**2) >= 0.9999d0) CheckPartOutsideBox = .true.
   end if

end function CheckPartOutsideBox

!************************************************************************
!*                                                                      *
!*     CheckTooFoldedChain                                              *
!*                                                                      *
!************************************************************************

! ... check if a too folded chain

logical function CheckTooFoldedChain(ip,bondeq)

   use CoordinateModule
   implicit none

   integer(4), intent(in) :: ip
   real(8),    intent(in) :: bondeq

   integer(4) :: ic, ict, iseg, jp
   real(8) :: dx, dy, dz, r2

   CheckTooFoldedChain = .false.
   iseg = isegpn(ip)                            ! segment
   if (iseg > 2) then
      ic = icnpn(ip)                            ! chain
      ict = ictcn(ic)                           ! chain type
      jp = ipnsegcn(iseg-2,ic)                  ! particle two bonds back
      dx = ro(1,ip)-ro(1,jp)
      dy = ro(2,ip)-ro(2,jp)
      dz = ro(3,ip)-ro(3,jp)
      call PBC(dx,dy,dz)
      r2 = dx**2+dy**2+dz**2                    ! separation squared
      if (r2 < Two*bondeq**2*(One-cos(anglemin(ict)))) CheckTooFoldedChain = .true.
   end if

end function CheckTooFoldedChain

!************************************************************************
!*                                                                      *
!*     MakeCrossLink                                                    *
!*                                                                      *
!************************************************************************

! ... make crosslinks between particles labeled clbeg and particles labeled clend
!     on the basis of their separation

!                     npclbeg, ipclbeg
!                     npclend, lpclend                nbondcl, bondcl
!     SetPeriodicNetwork ----------> MakeCrossLink -------------------> UCrossLink, etc
!

subroutine MakeCrossLink(npclbeg, npclend, ipclbeg, ipclend, r2bond, nclmade)

   use CoordinateModule
   implicit none

   character(40), parameter :: txroutine ='MakeCrossLink'

   integer(4), intent(in) :: npclbeg               ! current number of potential particles that begins a crosslink
   integer(4), intent(in) :: npclend               ! current number of potential particles that ends a crosslink
   integer(4), intent(in) :: ipclbeg(*)
   integer(4), intent(in) :: ipclend(*)
   real(8)   , intent(in) :: r2bond(*)
   integer(4), intent(out) :: nclmade              ! number of formed crosslinks

   real(8), parameter :: delta = 1.0d-5, faclow = One - delta, facupp = One + delta
   integer(4) :: ip, ipt, jp, jpt, jct, iclbeg, iclend, ic, jc
   real(8)    :: dx, dy, dz, r2

   nclmade = 0

   do iclbeg = 1, npclbeg
      ip = ipclbeg(iclbeg)            ! potential crosslink begining
      ipt = iptpn(ip)
      ic = icnpn(ip)
      do iclend = 1, npclend
         jp = ipclend(iclend)         ! potential crosslink end
         jpt = iptpn(jp)
         jc = icnpn(jp)
         if (ic == jc) cycle          ! exclude particles residing on the same chain (for chains with two particles)
         jct = ictpn(jp)
         dx = ro(1,jp)-ro(1,ip)
         dy = ro(2,jp)-ro(2,ip)
         dz = ro(3,jp)-ro(3,ip)
         call PBC(dx,dy,dz)
         r2 = dx**2+dy**2+dz**2
         if ((r2 > faclow*r2bond(jct)) .and. (r2 < facupp*r2bond(jct))) then                 ! bond separated (within delta)?
            if ((nbondcl(ip) < maxnbondcl(ipt)) .and. (nbondcl(jp) < maxnbondcl(jpt))) then  ! below the upper limit of crosslinks?
               nclmade = nclmade + 1
               nbondcl(ip) = nbondcl(ip) + 1
               nbondcl(jp) = nbondcl(jp) + 1
               bondcl(nbondcl(ip),ip) = jp
               bondcl(nbondcl(jp),jp) = ip
            end if
         end if
      end do
   end do

   if (master .and. itestcoordinate == 1) call TestMakeCrossLink(uout)

contains

!........................................................................

subroutine TestMakeCrossLink(unit)
   integer(4),   intent(in) :: unit
   integer(4) :: ibond
   call WriteHead(3,'Test'//trim(txroutine), unit)
   write(unit,'(a)') 'iclbeg       ic        ip'
   write(unit,'(i5,2i10)') (iclbeg, icnpn(ipclbeg(iclbeg)) ,ipclbeg(iclbeg), iclbeg = 1, npclbeg)
   write(unit,'(a)')
   write(unit,'(a)') 'iclend       ic        ip'
   write(unit,'(i5,2i10)') (iclend, icnpn(ipclend(iclend)) ,ipclend(iclend), iclend = 1, npclend)
   write(unit,'(a)')
   write(unit,'(a)') 'ip    ic          maxnbondcl(ipt)         nbondcl(ip)          bondcl(1,ip)  ...  '
   do ip = 1, np
      write(unit,'(2i5,2i18,10x,4i10)') ip, icnpn(ip), maxnbondcl(iptpn(ip)), nbondcl(ip), (bondcl(ibond,ip), ibond = 1, nbondcl(ip))
   end do
end subroutine TestMakeCrossLink

end subroutine MakeCrossLink

!************************************************************************
!*                                                                      *
!*     TestCoordinate                                                   *
!*                                                                      *
!************************************************************************

! ... write test coordinates

subroutine TestCoordinate(unit)
   use CoordinateModule
   implicit none
   character(80), parameter :: txheading ='TestCoordinate'
   character(22), parameter :: fmt3   = '(i4,t10,a,t23,3es12.4)'
   character(44), parameter :: fmt333 = '(i4,t10,a,t23,3es12.4,2x,3es12.4,2x,3es12.4)'
   integer(4),   intent(in) :: unit
   integer(4) :: ip
   call WriteHead(3, txheading, unit)
   write(unit,'(a,i5,a,i5)') 'macrostep = ', istep1, '   microstep = ', istep2
   write(unit,'(a)') 'particle position'
   write(unit,fmt3) (ip,txpt(iptpn(ip)),ro(1:3,ip),ip = 1,np)
   if (lpolyatom) then
      write(unit,'(a)') 'particle orientation'
      write(unit,fmt333) (ip,txpt(iptpn(ip)),ori(1:3,1:3,ip),ip = 1,np)
   end if
end subroutine TestCoordinate

!************************************************************************
!*                                                                      *
!*     StoreInteger                                                     *
!*                                                                      *
!************************************************************************

! ... store integer

subroutine StoreInteger(mnxxx, ip, nxxx, ipxxx)

   use CoordinateModule
   implicit none
   character(40), parameter :: txroutine ='StoreInteger'

   integer(4), intent(in) :: mnxxx           ! dimension
   integer(4), intent(in) :: ip              ! = 0 : initialize accmulator
                                             ! /=0 : store ip
   integer(4), intent(out) :: nxxx           ! accumulator (number of stored integers)
   integer(4), intent(out) :: ipxxx(mnxxx)   ! storage of ip

   if (ip == 0) then
      nxxx = 0
   else
      nxxx = nxxx + 1
      if (nxxx > mnxxx) then
           write(*,'(a,i6)') 'maximum number of integers that can be stored  =', mnxxx
           call Stop(txroutine, 'nxxx > mnxxx', uout)
      end if
      ipxxx(nxxx) = ip
   end if

end subroutine StoreInteger

