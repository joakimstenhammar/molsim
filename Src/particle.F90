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


!> \page nmlParticle
!! The namelist  \ref nmlParticle contains variables that describe the number of particles, their geometries, masses, etc. Examples of
!! namelist  \ref nmlParticle describing different types of particles are given in Appendix B.
!! * Variables:
!!  * \subpage txelec
!!  * \subpage lclink
!!  * \subpage lmultigraft
!!  * \subpage lhierarchical
!!  * \subpage maxnbondcl
!!  * \subpage ngen
!!  * \subpage ictgen
!!  * \subpage nbranch
!!  * \subpage ibranchpbeg
!!  * \subpage ibranchpinc
!!  * \subpage nnwt
!!  * \subpage nct
!!  * \subpage txct
!!  * \subpage ncct
!!  * \subpage npptct
!!  * \subpage txcopolymer
!!  * \subpage nblockict
!!  * \subpage npt
!!  * \subpage txpt
!!  * \subpage nppt
!!  * \subpage natpt
!!  * \subpage txat
!!  * \subpage massat
!!  * \subpage radat
!!  * \subpage zat
!!  * \subpage zatalpha
!!  * \subpage sigat
!!  * \subpage epsat
!!  * \subpage latweakcharge
!!  * \subpage jatweakcharge
!!  * \subpage pK
!!  * \subpage pH
!!  * \subpage naatpt
!!  * \subpage txaat
!!  * \subpage rain
!!  * \subpage dipain
!!  * \subpage polain
!!  * \subpage lintsite
!!  * \subpage raintin
!!  * \subpage lradatbox
!!  * \subpage itestpart

!************************************************************************
!> \page particle particle.F90
!! **ParticleModule**
!! *module for particle*
!************************************************************************

module ParticleModule

   use MolModule
!> \page txaat
!! `character(10)`(1:napt,1:\ref npt )
!! * \ref txaat (ialoc,ipt) is a text label for atom of no ialoc (local list) on particle of type ipt. ialoc runs from 1 to napt(ipt), the no of atoms of particle type ipt.
   character(10) :: txaat(mnapt,mnpt)
!> \page rain
!! `real`(1:3,1:napt,1:\ref npt )
!! * \ref rain (1:3,ialoc,ipt) is the x:y:z-coordinate of atom no ialoc (local list) in a particle of type ipt. ialoc runs from 1 to
!!   napt(ipt), the no of atoms of particle type ipt. \ref rain need not necessarily be given in the principle frame axes.
   real(8)       :: rain(3,mnapt,mnpt)
!> \page dipain
!! `real`(1:3,1:napt,1:\ref npt )
!! * \ref dipain (1:3,ialoc,ipt) is the x:y:z-component of the dipole moment of atom no ialoc (local list) in a particle of type ipt.
!!   ialoc runs from 1 to napt(ipt), the no of atoms of particle type ipt. \ref dipain has to be given in the same frame as \ref rain.
   real(8)       :: dipain(0:3,mnapt,mnpt)
!> \page polain
!! `real`(1:napt,1:\ref npt )
!! * \ref polain (1:6,ialoc,ipt) is the xx:yy:zz:xy:xz:yz-component of the symmetric polarizability of atom no ialoc (local list) in a
!!   particle of type ipt. ialoc runs from 1 to napt(ipt), the number of atoms of particle type ipt. \ref polain has to be given in the same
!!   frame as \ref rain.
   real(8)       :: polain(6,mnapt,mnpt)
!> \page raintin
!! `real`(1:3,1:napt,1:\ref npt )
!! * \ref raintin (1:3,ialoc,ipt) is the x:y:z-coordinate of interaction size no ialoc (local list) in a particle of type ipt. ialoc runs from 1 to napt(ipt), the no of atoms of particle type ipt.
   real(8)       :: raintin(3,mnapt,mnpt)  !*interaction site coordinates in input frame, cf rain
!> \page itestpart
!! `integer`
!! **default:** `0`
!! * Flag for test output. This possibility is for maintenance purposes.
!! * `0`: Nothing. The normal option.
!! * `10`: Chain pointers
   integer(4)    :: itestpart
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type :: block_type
      integer(4)  :: pt  !particle type
      integer(4)  :: np  !number of particles
   end type block_type

!> \page rep_block_ict
!! `block_type`(1:pt,1:np)
!! **default:** pt*np*`0`
!! * Particle type and number of particles in in each block and chain type.
   type(block_type), allocatable :: rep_iblock_ict(:,:)
!> \page nblockict
!! `integer`(1:\ref nct)
!! **default:** \ref nct*`0`
!! * Number of blocks in each repeating of a chain of a chain type.
   integer(4)  :: nblockict(mnct)
!> \page iptsegct
!! `integer`(maxval(npct(1:\ref nct)),1:\ref nct)
!! **default:** npct*\ref nct*`0`
!! * Particle type ipt of segments iseg of chain type ict
   integer(4), allocatable :: iptsegct(:,:)

end module ParticleModule



!> \page nmlNetworkConfiguration
!! The namelist \ref nmlNetworkConfiguration contains variables that describe the number of networks the particles and chains, which form the
!! networks and the the network topology.
!! * Variables:
!!  * \subpage nnwnwt
!!  * \subpage ncctnwt
!!  * \subpage txnwt
!!  * \subpage txtoponwt
!!  * \subpage iptclnwt

!> \page nmlCopolymerSequence
!! The namelist \ref nmlCopolymerSequence contains variables that describe the sequence of copolymers.
!! * Variables:
!!  * \subpage iptsegct


!> \page nmlRepeating
!! The namelist \ref nmlRepeating contains variables that define the repeating block structure of copolymers. The copolymers consist of
!! repeating units, each consisting of blocks of one particle type.
!! * Variables:
!!  * \subpage rep_block_ict

!************************************************************************
!> \page particle particle.F90
!!  **Particle(iStage)**
!! *particle variables*
!************************************************************************

subroutine Particle(iStage)

   use ParticleModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Particle'
   character(80), parameter :: txheading ='particle data'
   integer(4) :: igen, ialoc, inwt, ict, ic, jc, iseg, ipt, iat, iatloc, m
   logical                   :: luniformsequence
   character(10)             :: txhelp  ! auxiliary

   namelist /nmlParticle/ txelec,                                                 &
                          lclink, lmultigraft, maxnbondcl,                        &
                          nnwt,                                                   &
                          ngen, ictgen, nbranch, ibranchpbeg, ibranchpinc,        &
                          nct, txct, ncct, npptct, txcopolymer, lspma,            &
                          nblockict,                                              &
                          npt, txpt, nppt, natpt,                                 &
                          txat, massat, radat, zat, zatalpha, sigat, epsat,       &
                          naatpt, txaat, rain, dipain, polain, lintsite, raintin, &
                          latweakcharge, pK, pH, jatweakcharge,                   &
                          lradatbox, itestpart

   namelist /nmlRepeating/  rep_iblock_ict

   namelist /nmlCopolymerSequence/ iptsegct

   namelist /nmlNetworkConfiguration/ nnwnwt, ncctnwt, txnwt, txtoponwt, iptclnwt

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

! ... set initial values

      txelec                ='charge'
      lclink                =.false.
      lmultigraft           =.false.
      nnwt                  = 0
      maxnbondcl            = 1
      ngen                  =-1
      ictgen(0:mngen)       = 1
      nbranch(0:mngen-1)    = 0
      ibranchpbeg(0:mngen-1)= 1
      ibranchpinc(0:mngen-1)= 1
      nct                   = 0
      npptct          = 0
      txcopolymer     = 'block'
      lspma           = .false.
      nblockict       = 0
      massat          = One
      zat             = Zero
      zatalpha        = Zero
      sigat           = Zero
      epsat           = Zero
      rain            = Zero
      dipain          = Zero
      polain          = Zero
      lintsite        =.false.
      raintin         = Zero
      lradatbox       =.false.
      itestpart       = 0
      latweakcharge   = .false.
      pK              = Zero
      pH              = Zero
      jatweakcharge   = 0

! ... read input data (nmlParticle)

      rewind(uin)
      read(uin,nmlParticle)

      call LowerCase(txelec)
      do ict = 1, nct
         call LowerCase(txcopolymer(ict))
      end do

! ... read input data (nmlRepeating)

      if(any(txcopolymer(1:nct) == 'repeating')) then
         if(.not. allocated(rep_iblock_ict)) then
            allocate(rep_iblock_ict(maxval(nblockict(1:nct)),nct))
         end if

         rep_iblock_ict = block_type(0,0)
         rewind(uin)
         read(uin,nmlRepeating)
      end if

! ... read input data (nmlCopolymerSequence)

      if(any(txcopolymer(1:nct) == 'sequence')) then
         do ict = 1, nct
            npct(ict) = sum(npptct(1:npt,ict))  ! npct needed for following allocation
         end do
         if(.not. allocated(iptsegct)) then
            allocate(iptsegct(maxval(npct(1:nct)),nct))
         end if
         iptsegct = 0
         rewind(uin)
         read(uin,nmlCopolymerSequence)
      end if

! ... read input data (nmlNetworkConfiguration)

      if (nnwt > 0) then
         if (.not.allocated(nnwnwt)) then
            allocate(nnwnwt(nnwt),ncctnwt(nct,nnwt),iptclnwt(nnwt),txtoponwt(nnwt), &
                     txnwt(nnwt),nctnwt(nnwt),ncnwt(nnwt),npnwt(nnwt),nclnwt(nnwt), &
                     npweakchargenwt(nnwt),lptnwt(npt,nnwt))
            nctnwt          = 0
            ncnwt           = 0
            npnwt           = 0
            nclnwt          = 0
            npweakchargenwt = 0
            lptnwt          = .false.
         end if

         ! ... default values of input variables
         nnwnwt(1:nnwt)        = 0
         ncctnwt(1:nct,1:nnwt) = 0
         iptclnwt(1:nnwt)      = 0
         txtoponwt(1:nnwt)     = 'default'
         txnwt(1:nnwt)         = 'network'

         ! ... read input
         rewind(uin)
         read(uin,nmlNetworkConfiguration)
      end if

! ... determine types of atoms

      lcharge = .false.
      lweakcharge = .false.
      ldipole = .false.
      lpolarization = .false.
      laimage = .false.
      ldipolesph = .false.
      ldieldis = .false.
      if (txelec == 'charge') then
         lcharge =.true.
      else if (txelec == 'weakcharge') then
         lweakcharge =.true.
      else if (txelec == 'dip') then
         ldipole =.true.
      else if (txelec == 'pol') then
         lpolarization =.true.
      else if (txelec(1:6) == 'dipsph') then
         ldipolesph =.true.
         if (txelec(7:11) == 'image') laimage = .true.
      else if (txelec == 'dieldis') then
         ldieldis =.true.
      else
         call Stop(txroutine, 'error in txelec', uout)
      end if

! ... check input

      if (lmc .and. lintsite) call Stop(txroutine, 'nonatomic interaction sites not implemented for lmc', uout)

      if (any(sum(npptct(:,1:nct),DIM=1) == 0)) then
         call Stop(txroutine, 'a chain without any particles was defined', uout)
      end if
      if (any(nppt(1:npt) <= 0)) then
         call Stop(txroutine, 'a particle type without any particles was defined', uout)
      end if

! ... take case of double meaning of dipain

      do ipt = 1, mnpt
         do ialoc = 1, mnapt
            if (dipain(0,ialoc,ipt)/=Zero) &
            dipain(1:3,ialoc,ipt) = dipain(0,ialoc,ipt)*dipain(1:3,ialoc,ipt)/sqrt(sum(dipain(1:3,ialoc,ipt)**2))
         end do
      end do

! ... change from one to several sizes in a unit cube

      if (txuser(1:5) == "cube_") call cube

      if (txuser == 'jos') call JosUser(1)

! ... set object parameters

      call SetObjectParam1

! ... check condition

      if (lweakcharge .and. lewald .and. lpolyatom) then
         call Stop(txroutine,'(lweakcharge .and. lnetwork) not adapted for polyatomic systems',uout)
      end if
      if (lweakcharge) then
         if (.not. any(latweakcharge(1:nat))) then
            call Stop(txroutine, 'no titratable atom type',uout)
         end if
      end if

! ... allocate memory

      if (.not.allocated(ro)) then
         allocate(ro(3,np_alloc))
         ro = 0.0E+00
      end if
      if (.not.allocated(rotm)) then
         allocate(rotm(3,np_alloc))
         rotm = 0.0E+00
      end if
      if (.not.allocated(drotm)) then
         allocate(drotm(3,np_alloc))
         drotm = 0.0E+00
      end if
      if (.not.allocated(ori)) then
         allocate(ori(3,3,np_alloc))
         ori = 0.0E+00
      end if
      if (.not.allocated(oritm)) then
         allocate(oritm(3,3,np_alloc))
         oritm = 0.0E+00
      end if
      if (.not.allocated(rod)) then
         allocate(rod(3,np_alloc))
         rod = 0.0E+00
      end if
      if (.not.allocated(qua)) then
         allocate(qua(0:3,np_alloc))
         qua = 0.0E+00
      end if
      if (.not.allocated(quad)) then
         allocate(quad(0:3,np_alloc))
         quad = 0.0E+00
      end if
      if (.not.allocated(r)) then
         allocate(r(3,na_alloc))
         r = 0.0E+00
      end if
      if (.not.allocated(rtm)) then
         allocate(rtm(3,na_alloc))
         rtm = 0.0E+00
      end if

      if (.not.allocated(dip)) then
         allocate(dip(3,na_alloc))
         dip = 0.0E+00
      end if
      if (.not.allocated(diptm)) then
         allocate(diptm(3,na_alloc))
         diptm = 0.0E+00
      end if
      if (.not.allocated(idm)) then
         allocate(idm(3,na_alloc))
         idm = 0.0E+00
      end if
      if (.not.allocated(idm1)) then
         allocate(idm1(3,na_alloc))
         idm1 = 0.0E+00
      end if
      if (.not.allocated(idm2)) then
         allocate(idm2(3,na_alloc))
         idm2 = 0.0E+00
      end if
      if (.not.allocated(idmo)) then
         allocate(idmo(3,np_alloc))
         idmo = 0.0E+00
      end if
      if (.not.allocated(diptot)) then
         allocate(diptot(3,na_alloc))
         diptot = 0.0E+00
      end if
      if (.not.allocated(vaux)) then
         allocate(vaux(3,max(1000,2*na_alloc)))
         vaux = 0.0E+00
      end if
      if (.not.allocated(ivaux)) then
         allocate(ivaux(3,max(1000,2*na_alloc)))
         ivaux = 0
      end if

      if (.not.allocated(angvelo)) then
         allocate(angvelo(3,np_alloc))
         angvelo = 0.0E+00
      end if
      if (.not.allocated(forceo)) then
         allocate(forceo(3,np_alloc))
         forceo = 0.0E+00
      end if
      if (.not.allocated(torqueo)) then
         allocate(torqueo(3,np_alloc))
         torqueo = 0.0E+00
      end if
      if (.not.allocated(force)) then
         allocate(force(3,na_alloc), torque(3,na_alloc))
         force = 0.0E+00
         torque = 0.0E+00
      end if
      if (.not.allocated(drostep)) then
         allocate(drostep(3,np_alloc))
         drostep = 0.0E+00
      end if

      if (.not.allocated(laztm)) then
         allocate(laztm(na_alloc))
         laztm = .false.
      end if
      if (.not.allocated(aztm)) then
         allocate(aztm(na_alloc))
         aztm = 0.0E+00
      end if

! ... set object parameters

      call SetObjectParam2

      if (lpolarization) idm = zero                ! initialize induced dipole moment

   case (iWriteInput)

      if (master) then
         do igen = 0, ngen-1
            ipt = ipnsegcn(1,icnct(ictgen(igen)))
            if (ibranchpbeg(igen) + ibranchpinc(igen)*(nbranch(igen)-1) > npct(ictgen(igen))) call Stop(txroutine,'chain too short', uout)
            if (nbranch(igen) * ncct(ictgen(igen)) /= ncct(ictgen(igen+1))) call Stop(txroutine, 'no matching of number of chains',uout)
         end do

   ! ... particle data

         call WriteHead(2, txheading, uout)
         if (txelec == 'charge') then
            write(uout,'(a)') 'charged atoms enabled'
         else if (txelec == 'dip') then
            write(uout,'(a)') 'charged and dipolar atoms enabled'
         else if (txelec == 'pol') then
            write(uout,'(a)') 'charged, dipolar, and polarizable atoms enabled'
         end if
         write(uout,'()')

         if (lnetwork) then
         write(uout,'(a,i7)')             'number of networks             = ', nnw
         end if
         if (lhierarchical) then
         write(uout,'(a,i7)')             'number of hierarchical struc.  = ', nh
         end if
         if (lchain) then
         write(uout,'(a,i7)')             'number of chains               = ', nc
         end if
         write(uout,'(a,i7)')             'number of particles            = ', np
         write(uout,'(a,i7)')             'number of atoms                = ', na
         if (lnetwork) then
         write(uout,'(a,i7)')             'number of network types        = ', nnwt
         end if
         if (lchain) then
         write(uout,'(a,i7,5x,a,i7,x,a)') 'number of chains types         = ', nct, '(',mnct,')'
         end if
         write(uout,'(a,i7,5x,a,i7,x,a)') 'number of particle types       = ', npt, '(',mnpt,')'
         write(uout,'(a,i7,5x,a,i7,x,a)') 'number of atom types           = ', nat, '(',mnat,')'

         if (lmultigraft) then
            write(uout,'()')
            write(uout,'(a,t55,l8)')  'multigrafted chain                                  = ', lmultigraft
         end if

         if (lnetwork) then
            write(txhelp,'(i3)') nct
            write(uout,'()')
            write(uout,'(a,t14,a,t24,a,t32,a,t46,'//trim(adjustl(txhelp))//'(a8,i0,a6,2x))') &
            ' network  ',' type ', ' no ', ' topology ', ('ncctnwt(',ict,',inwt)',ict = 1, nct)
            write(uout,'(a,t14,a,t24,a,t32,a,t46,'//trim(adjustl(txhelp))//'(a15,2x))')       &
            '----------','------', '----', '----------', ('---------------',ict = 1, nct)
            do inwt = 1, nnwt
               write(uout,'(1x,a,t14,1x,i2,t24,i2,t32,a,t46,'//trim(adjustl(txhelp))//'(i5,11x))') &
               txnwt(inwt), inwt, nnwnwt(inwt), txtoponwt(inwt), ncctnwt(1:nct,inwt)
            end do
            write(uout,'()')
         end if

         if (lhierarchical) then
            write(uout,'()')
            write(uout,'(a)')         'hierarchical data'
            write(uout,'(a)')         '-----------------'
            write(uout,'(a,t55,i8)')  'maximum number of crosslinks of a particle          = ', maxval(maxnbondcl(1:npt))
            write(uout,'(a,t55,i8)')  'number of generations                               = ', ngen
            write(uout,'(a,t55,8i8)') 'chain type of the generations                       = ', ictgen(0:ngen)
            write(uout,'(a,t55,8i8)') 'number of branches                                  = ', nbranch(0:ngen-1)
            write(uout,'(a,t55,8i8)') 'particle number of chain for first branch point     = ', ibranchpbeg(0:ngen-1)
            write(uout,'(a,t55,8i8)') 'particle increment of chain between branch points   = ', ibranchpinc(0:ngen-1)
            write(uout,'(a,t55,5i8)') 'number of crosslinks                                = ', ncl
         end if

         if (lchain) then
            write(txhelp,'(i3)') npt
            write(uout,'()')
            write(uout,'(a,t20,a,t45,a,t55,'//trim(adjustl(txhelp))//'(a7,i0,a5,2x))') &
            'chain type', 'no of chains', 'topology', ('npptct(',ipt,',ict) ',ipt = 1, npt)
            write(uout,'(a,t20,a,t45,a,t55,'//trim(adjustl(txhelp))//'(a13,2x))')       &
            '----------', '------------', '--------', ('-------------', ipt = 1, npt)
            do ict = 1, nct
               write(uout,'(a,t20,i5,t45,a,t55,'//trim(adjustl(txhelp))//'(i5,9x))') &
               txct(ict), ncct(ict), txcopolymer(ict), npptct(1:npt,ict)
            end do
            write(uout,'()')
            do ict = 1, nct
               if (txcopolymer(ict) == 'block') cycle ! No sequence output for block-like chain types
               write(txhelp,'(i3)') ict
               write(uout,'(a)') repeat('- ',13)//'sequence of chain type '//trim(adjustl(txhelp))&
                                 &//' '''//trim(adjustl(txct(ict)))//''''//repeat(' -',13)
               write(uout,'()')
               luniformsequence = .true.
  testuniform: do ic = icnct(ict), icnct(ict) + ncct(ict) - 1
                  do jc = ic + 1, icnct(ict) + ncct(ict) - 1
                     if (any(iptpn(ipnsegcn(1:npct(ict),ic)) /= iptpn(ipnsegcn(1:npct(ict),jc)))) then
                        luniformsequence = .false.
                        exit testuniform
                     end if
                  end do
               end do testuniform
               do ic = icnct(ict), icnct(ict)+ncct(ict)-1
                  if (.not. luniformsequence) write(uout,'(a5,i0)') 'ic = ', ic
                  do iseg = 1, npct(ict)
                     write(txhelp,'(i10)') iptpn(ipnsegcn(iseg,ic))
                     txhelp = merge(trim(adjustl(txhelp))//' ',trim(adjustl(txhelp))//'-',iseg == npct(ict))
                     if (modulo(iseg,50) == 0 .or. iseg == npct(ict)) then
                        write(uout,fmt='(a)',advance='yes') trim(adjustl(txhelp))
                     else
                        write(uout,fmt='(a)',advance='no')  trim(adjustl(txhelp))
                     end if
                  end do
                  if (luniformsequence) exit
                  write (uout,'()')
               end do
               write (uout,'()')
               if (ict == nct) write(uout,'(a)') repeat('- ',45)
            end do
            if (lspma) write(uout,'(a)')
            if (lspma) write(uout,'(a,l5)') 'lspma                                               = ', lspma
         end if

         write(uout,'()')
         write(uout,'(a,t16,a,t26,a,t35,a,t42,a,t50,a,t63,a,t76,a,t86,a,t96,a,t105,a,t114,a)')       &
         'particle', 'no', 'atom', 'no', 'mass', 'hard core', 'no of unit', 'weak  ', 'pK', 'sigma', 'epsilon', 'zatalpha'
         write(uout,'(a,t16,a,t26,a,t35,a,t42,a,t50,a,t63,a,t76,a,t86,a,t96,a,t105,a,t114,a)')       &
         '        ', '  ', '    ', '  ', '    ', 'radius   ', 'charges   ', 'charge', '      ','     ', '       '
         write(uout,'(a,t16,a,t26,a,t35,a,t42,a,t50,a,t63,a,t76,a,t86,a,t96,a,t105,a,t114,a)')       &
         '--------', '--', '----', '--', '----', '---------', '----------', '------', '--', '-----', '-------', '--------'

         iat = 0
         do ipt = 1, npt
            write(uout,'(a,t12,i6)') txpt(ipt), nppt(ipt)
            do iatloc = 1, natpt(ipt)
               iat = iat+1
               write(uout,'(t26,a,t32,i4,t39,f8.3,t48,f8.3,t60,f10.4,t72,l7,t81,f10.4,t92,f10.4,t102,f10.4,t112,f10.4)') &
                  txat(iat), naat(iat), massat(iat), radat(iat), zat(iat), latweakcharge(iat), pK(iat), sigat(iat), epsat(iat), zatalpha(iat)
            end do
         end do

         if (lweakcharge) then
            write(uout,'()')
            write(uout,'(a,t6,f8.3)') 'pH = ', pH
            write(uout,'(a,t25,a,t34,a)') 'titratable atom type', 'pK' ,'atom type of counterion'
            write(uout,'(a,t25,a,t34,a)') '--------------------', '-----', '-----------------------'
            do iat = 1, nat
               if (.not. latweakcharge(iat)) cycle
               write(uout,'(t10,i0,t22,f8.3,t40,i0)') iat, pK(iat) , jatweakcharge(iat)
            end do
         end if

         if (lradatbox) then
            write(uout,'()')
            write(uout,'(a)') 'atoms and atom radii are used to examine if particles are inside the box'
            write(uout,'()')
         end if

         do ipt = 1, npt

            write(uout,'()')
            write(uout,'(a,a)') 'molecular input frame: ', txpt(ipt)
            write(uout,'(a,40a)') '---------------------- ', ('-',m = 1,len(trim(txpt(ipt))))
            write(uout,'()')
            write(uout,'(a,t20,a)') 'atom', 'coordinate (x,y,z)'
            write(uout,'(a,t20,a)') '----', '------------------'
            write(uout,'(a,3f10.5)') (txaat(ialoc,ipt), rain(1:3,ialoc,ipt), ialoc = 1,napt(ipt))

            if (ldipole .or. ldipolesph) then
               write(uout,'()')
               write(uout,'(a,t20,a,t60,a)') 'atom', 'dipole moment (x,y,z)'
               write(uout,'(a,t20,a,t60,a)') '----', '---------------------'
               write(uout,'(a,3f10.5)') (txaat(ialoc,ipt), dipain(1:3,ialoc,ipt), ialoc = 1, napt(ipt))
            end if

            if (lpolarization) then
               write(uout,'()')
               write(uout,'(a,t20,a,t60,a)') 'atom', 'dipole moment (x,y,z)', 'polarizability (xx,yy,zz,xy,xz,yz)'
               write(uout,'(a,t20,a,t60,a)') '----', '---------------------', '----------------------------------'
                write(uout,'(a,3f10.5,5x,6f10.5)')                                            &
                (txaat(ialoc,ipt), dipain(1:3,ialoc,ipt), polain(1:6,ialoc,ipt), ialoc = 1, napt(ipt))
            end if

            write(uout,'()')
            write(uout,'(a,a)') 'principal axis frame: ', txpt(ipt)
            write(uout,'(a,40a)') '--------------------- ', ('-',m = 1,len(trim(txpt(ipt))))
            write(uout,'()')
            write(uout,'(a,t20,a)') 'atom', 'coordinate (x,y,z)'
            write(uout,'(a,t20,a)') '----', '------------------'
            write(uout,'(a,3f10.5)') (txaat(ialoc,ipt), ra(1:3,ialoc,ipt), ialoc = 1,napt(ipt))

            if (ldipole .or. ldipolesph) then
               write(uout,'()')
               write(uout,'(a,t20,a,t60,a)') 'atom', 'dipole moment (x,y,z)'
               write(uout,'(a,t20,a,t60,a)') '----', '---------------------'
               write(uout,'(a,3f10.5)') (txaat(ialoc,ipt), dipa(1:3,ialoc,ipt), ialoc = 1,napt(ipt))
            end if

            if (lpolarization) then
               write(uout,'()')
               write(uout,'(a,t20,a,t60,a)') 'atom', 'dipole moment (x,y,z)', 'polarizability (xx,yy,zz,xy,xz,yz)'
               write(uout,'(a,t20,a,t60,a)') '----', '---------------------', '----------------------------------'
               write(uout,'(a,3f10.5,5x,6f10.5)') (txaat(ialoc,ipt), dipa(1:3,ialoc,ipt), poltensa(1:6,ialoc,ipt), ialoc = 1,napt(ipt))
            end if

         end do

         write(uout,'()')
         write(uout,'(a,t22,a,t37,a,t70,a)') 'particle', 'mass', 'moment of inertia (x,y,z)', 'maximal atom-com distance'
         write(uout,'(a,t22,a,t37,a,t70,a)') '--------', '----', '-------------------------', '-------------------------'
         write(uout,'(a,5x,f10.3,5x,3f10.3,10x,f8.3)') (txpt(ipt), masspt(ipt), mompt(1:3,ipt), racom(ipt), ipt = 1,npt)
      end if

   end select

contains

!........................................................................

subroutine cube                     ! change from one to nside**3 homongeneously distributed atoms in a cube
   integer(4) :: nside, ix, iy, iz, nn
   real(8)    :: dip(3)
   read(txuser(6:6),'(i1)') nside
   dip(1:3) = dipain(1:3,1,1)/nside**3
   nn = 0
   do ix = 1, nside
      do iy = 1, nside
         do iz = 1, nside
            nn = nn + 1
            rain(1,nn,1) = -half + (ix-half)/nside
            rain(2,nn,1) = -half + (iy-half)/nside
            rain(3,nn,1) = -half + (iz-half)/nside
         end do
      end do
   end do
   txaat(1:nn,1) = txaat(1,1)
   naatpt(1,1) = nn
   dipain(1,1:nn,1) = dip(1)
   dipain(2,1:nn,1) = dip(2)
   dipain(3,1:nn,1) = dip(3)
end subroutine cube

!........................................................................

end subroutine Particle

!************************************************************************
!> \page particle particle.F90
!!  **SetObjectParam1**
!! *setobject parameters*
!************************************************************************

subroutine SetObjectParam1

   use ParticleModule
   implicit none

   character(40), parameter :: txroutine ='SetObjectParam1'
   integer(4) :: ia, ialoc, iat, iatloc         ! atoms
   integer(4) :: ip, iploc, ipt, nptemp         ! particles
   integer(4) :: iseg, ic, icloc, ict, nctemp   ! chains
   integer(4) :: igen                           ! hierarchical
   integer(4) :: inw, inwt, inwloc, iclloc      ! network

! ... check some input variables

   if (ngen > mngen)                   call Stop(txroutine, 'ngen > mngen', uout)
   if (count(natpt(1:npt) > mnat) > 0) call Stop(txroutine, 'natpt > mnat', uout)
   if (nct > mnct)                     call Stop(txroutine, 'nct > mnct', uout)
   if (npt < 1)                        call Stop(txroutine, 'npt < 1', uout)
   if (npt > mnpt)                     call Stop(txroutine, 'npt > mnpt', uout)

! ... check consistence among ncct, npptct, and nppt

   do ipt = 1, npt
      nptemp = sum(ncct(1:nct)*npptct(ipt,1:nct))
      if ( .not.((nptemp == 0).or.(nptemp == nppt(ipt)))) then
         write(uout,'(a,i5)') 'ipt = ', ipt
         write(uout,'(a,i5)') 'sum(ncct(1:nct)*npptct(ipt,1:nct)) = ', sum(ncct(1:nct)*npptct(ipt,1:nct))
         write(uout,'(a,i5)') 'nppt(ipt) = ', nppt(ipt)
         call Stop(txroutine, 'inconsistency among ncct, npptct, and nppt', uout)
      end if
   end do

! ... check consistence of npptct and iptsegct (only if txcopolymer = 'sequence')
   do ict = 1, nct
      if (txcopolymer(ict) == 'sequence') then
         do ipt = 1, npt
            if (npptct(ipt,ict) /= count(iptsegct(1:npct(ict),ict) == ipt)) then
               write(uout,'(a,i5)') 'ict = ', ict
               write(uout,'(a,i5)') 'ipt = ', ipt
               write(uout,'(a,i5)') 'count(iptsegct(1:npct(ict),ict) == ipt) = ', &
                                     count(iptsegct(1:npct(ict),ict) == ipt)
               write(uout,'(a,i5)') 'npptct(ipt,ict) = ', npptct(ipt,ict)
               call Stop(txroutine, 'copolymer sequence: inconsistency among npptct and iptsegct', uout)
            end if
         end do
      end if
   end do

! ... check consistence among nnwnwt and ncct

   if(nnwt > 0) then
      do ict = 1, nct
         nctemp = sum(nnwnwt(1:nnwt)*ncctnwt(ict,1:nnwt))
         if (.not.((nctemp == 0).or.(nctemp == ncct(ict)))) then
            write(uout,'(a,i5)') 'ict = ', ict
            write(uout,'(a,i5)') 'sum(nnwnwt(1:nnwt)*ncctnwt(ict,1:nnwt)) = ', sum(nnwnwt(1:nnwt)*ncctnwt(ict,1:nnwt))
            write(uout,'(a,i5)') 'ncct(ict) = ', ncct(ict)
            call Stop(txroutine, 'inconsistency among nnwnwt, ncctnwt, and ncct', uout)
         end if
      end do
   end if

! ... check whether nppt(iptclnwt(inwt)) is a multiple of nnwnwt(inwt)

   do inwt = 1, nnwt
      if (.not.(modulo(nppt(iptclnwt(inwt)),nnwnwt(inwt)) == 0)) then
         write(uout,'(a,i5)') 'inwt = ', inwt
         write(uout,'(a,i5)') 'nppt(iptclnwt(inwt) = ', nppt(iptclnwt(inwt))
         write(uout,'(a,i5)') 'nnwnwt(inwt) = ', nnwnwt(inwt)
         call Stop(txroutine, 'number of node particles no multiple of number of corresponding network', uout)
      end if
   end do

   call Set_lnetwork      ! flag for network structures
   call Set_nh            ! number of hierarchical structures
   call Set_lhierarchical ! flag for hierarchical structures
   call Set_lchain        ! flag for chains
   call Set_ncl           ! number of cross-links
   call Set_nc            ! number of chains
   call Set_nptct         ! number of particle types of a chain type
   call Set_npct          ! number of particles a chain type
   call Set_nphn          ! number of particles of the hierarchical structure
   call Set_napt          ! number of atoms of a particle type
   call Set_np            ! number of particles
   call Set_np_alloc      ! number of particles for memory allocation
   call Set_na            ! number of atoms
   call Set_na_alloc      ! number of atoms for memory allocation
   call Set_nat           ! number of atom types
   call Set_naat          ! number of atoms of a given atom type
   call Set_lmonoatom     ! flag for monoatomic particles
   call Set_nctct         ! number of chain type pairs
   call Set_nptpt         ! number of particle type pairs
   call Set_natat         ! number of atom type pairs

   call Set_ipnsegcn      ! chain and segment  -> particle
   call Set_ictcn         ! chain              -> its chain type
   call Set_ictpt         ! particle type      -> its chain type
   call Set_ictpn         ! particle           -> its chain type
   call Set_icnpn         ! particle           -> its chain
   call Set_iptpn         ! particle           -> its particle type
   call Set_iptat         ! atom type          -> its particle type
   call Set_iptan         ! atom               -> its particle type
   call Set_ipnan         ! atom               -> its particle
   call Set_iatan         ! atom               -> its atom type
   call Set_icnct         ! chain type         -> its first chain
   call Set_ipnpt         ! particle type      -> its first particle
   call Set_iatpt         ! particle type      -> its first atom type
   call Set_ianpn         ! particle           -> its first atom
   call Set_ianat         ! atom type          -> its first atom
   call Set_ictct         ! two chain types    -> chain type pair
   call Set_iptpt         ! two particle types -> particle type pair
   call Set_iatat         ! two atom types     -> atom type pair
   call Set_isegpn        ! particle number    -> segment number
   call Set_bondnn        ! bond and particle  -> bonded particle

   if (lnetwork) then
      call Set_lptnwt     ! particle type ipt used for network type inwt?
      call Set_nnw        ! number of networks
      call Set_nclnwt     ! number of cross-links of network type
      call Set_nctnwt     ! number of chain types of network type inwt
      call Set_ncnwt      ! number of chains of network type inwt
      call Set_npnwt      ! number of particles of network type inwt
      call Set_nnwtnwt    ! number of network type pairs
      call Set_npweakchargenwt ! number of titratable particles in network type inwt

      call Set_inwtnwn    ! network                                      -> its network type
      call Set_inwtct     ! chain type (1:nct)                           -> its network type (1:nnwt)
      call Set_icnclocnwn ! local chain (1:ncnwt) and network (1:nnw)    -> its chain (1:nc)
      call Set_inwtcn     ! chain (1:nc)                                 -> its network type (1:nnwt)
      call Set_inwncn     ! chain (1:nc)                                 -> its network (1:nnw)
      call Set_inwnnwt    ! network type (1:nnwt)                        -> its first network (1:nnw)
      call Set_inwtnwt    ! two network types (1:nnwt)                   -> network type pair (1:nnwt)
      call Set_ipncllocnwn! cross-link (1:nclnwt) and network (1:nnw)    -> its particle (1:np)
      call Set_ipnplocnwn ! local particle (1:npnwt) and network (1:nnw) -> its particle (1:np)
      call Set_lpnnwn     ! particle (1:np) part of network (1:nnw)?
   end if

   if (lhierarchical) then
       call Set_ipnhn  ! hierarchical strcture -> its first particle
       call Set_genic  ! chain number -> generation number
       call Set_bondcl ! crosslink and particle -> crosslinked particle
       call Set_ihnpn  ! particle -> its hierarchical structure
   end if

   if (lclink .and. .not.lhierarchical) then
      maxvalnbondcl = maxval(maxnbondcl(1:npt))
      if (.not.allocated(nbondcl)) then
         allocate(nbondcl(np_alloc))
         nbondcl = 0
      end if
      if (.not.allocated(bondcl)) then
         allocate(bondcl(maxvalnbondcl,np_alloc))
         bondcl = 0
      end if
   end if

   if (master .and. itestpart == 10) then
      call TestChainPointer(uout)
      if (lnetwork) call TestNetworkPointer(uout)
      if (lhierarchical) call TestCrosslinkPointer(uout)
   end if

contains

!........................................................................

subroutine Set_lnetwork ! flag for network structures
   lnetwork =.false.
   if (nnwt > 0) lnetwork =.true.
end subroutine Set_lnetwork

!........................................................................

subroutine Set_nh  ! number of hierarchical structures
   nh = 0
   if (ngen > 0) nh = ncct(ictgen(0))
end subroutine Set_nh

!........................................................................

subroutine Set_lhierarchical  ! flag for hierarchical structures
   lhierarchical =.false.
   if (nh > 0) lhierarchical =.true.
end subroutine Set_lhierarchical

!........................................................................

subroutine Set_lchain  ! flag for chains
   lchain =.false.
   if (nct > 0) lchain =.true.
end subroutine Set_lchain

!........................................................................

subroutine Set_ncl ! number of cross-links
   ncl = 0
   if(ngen > -1) then
      ncl = sum(nbranch(0:ngen-1)*ncct(ictgen(0:ngen-1)))
   end if
end subroutine Set_ncl

!........................................................................

subroutine Set_nc ! number of chains
   nc = sum(ncct(1:nct))
end subroutine Set_nc

!........................................................................

subroutine Set_nptct  ! number of particle types of a chain type
   do ict = 1, nct
      nptct(ict) = count(npptct(1:npt,ict) > 0)
   end do
end subroutine Set_nptct

!........................................................................

subroutine Set_npct  ! number of particle types of a chain type
   do ict = 1, nct
      npct(ict) = sum(npptct(1:npt,ict))
   end do
end subroutine Set_npct

!........................................................................

subroutine Set_nphn  ! number of particles of the hierarchical structure
   nphn = 0
   do igen = 0, ngen
      nphn = nphn + npct(ictgen(igen))*(ncct(ictgen(igen))/nh) ! normalize by nh to get the numbers in one hierarchical structure
   end do
end subroutine Set_nphn

!........................................................................

subroutine Set_napt  ! number of atoms of a particle type
   if (.not.allocated(napt)) then
      allocate(napt(npt))
      napt = 0
   end if
   do ipt = 1, npt
      napt(ipt) = sum(naatpt(1:natpt(ipt),ipt))
   end do
   if (count(napt(1:npt) < 1    ) > 0) call Stop(txroutine, 'napt(ipt) < 1    ', uout)
   if (count(napt(1:npt) > mnapt) > 0) call Stop(txroutine, 'napt(ipt) > mnapt', uout)
end subroutine Set_napt

!........................................................................

subroutine Set_np  ! number of particles
   np = sum(nppt(1:npt))
   if (np < 1) call Stop(txroutine, 'np < 1', uout)
end subroutine Set_np

!........................................................................

subroutine Set_np_alloc ! number of particles for memory allocation
   logical :: first = .true.
   if (lmvt) then
      if (first) then
        np_alloc = facmvt*np
        first = .false.
      end if
   else
      np_alloc = np
   end if
   if (np > np_alloc) call Stop(txroutine, 'np > np_alloc', uout)
end subroutine Set_np_alloc

!........................................................................

subroutine Set_na  ! number of atoms
   na = sum(napt(1:npt)*nppt(1:npt))
   if (na < 1)  call Stop(txroutine, 'na < 1', uout)
end subroutine Set_na

!........................................................................

subroutine Set_na_alloc ! number of atoms for memory allocation
   logical :: first = .true.
   if (lmvt) then
      if (first) then
         na_alloc = facmvt*na
        first = .false.
      end if
   else
      na_alloc = na
   end if
   if (na > na_alloc) call Stop(txroutine, 'na > na_alloc', uout)
end subroutine Set_na_alloc

!........................................................................

subroutine Set_nat  ! number of atom types
   nat = sum(natpt(1:npt))
   if (nat < 1)  call Stop(txroutine, 'nat < 1', uout)
end subroutine Set_nat

!........................................................................

subroutine Set_naat  ! number of atoms of a given atom type
   if (.not.allocated(naat)) then
      allocate(naat(nat))
      naat = 0
   end if
   iat = 0
   do ipt = 1, npt
      do iatloc = 1, natpt(ipt)
         iat = iat+1
         naat(iat) = naatpt(iatloc,ipt)
      end do
   end do
end subroutine Set_naat

!........................................................................

subroutine Set_lmonoatom  ! flag for monoatomic particles
   lmonoatom =.true.
   if (count(napt > 1) > 0)   lmonoatom =.false.
   if (maxval(dipain) > Zero) lmonoatom =.false.
   lpolyatom =.not.lmonoatom
end subroutine Set_lmonoatom

!........................................................................

subroutine Set_nctct  ! number of chain type pairs
   nctct = (nct*(nct+1))/2
end subroutine Set_nctct

!........................................................................

subroutine Set_nptpt  ! number of particle type pairs
   nptpt = (npt*(npt+1))/2
end subroutine Set_nptpt

!........................................................................

subroutine Set_natat  ! number of atom type pairs
   natat = (nat*(nat+1))/2
end subroutine Set_natat

!........................................................................

subroutine Set_ipnsegcn  ! chain and segment -> particle

   character(40), parameter :: txroutine ='Set_ipnsegcn'
   integer(4) :: nrep, irep, nreplen
   integer(4) :: iblock
   integer(4) :: iplow
   integer(4) :: ipt
   integer(4), allocatable :: npset(:)
   integer(4), allocatable :: ipstart(:)
   integer(4), allocatable :: iptiseg(:)

   if (.not.allocated(ipnsegcn)) then
      allocate(ipnsegcn(maxval(npct(1:nct)),nc))   ! defined in MolModule
      ipnsegcn = 0
   end if
   ic = 0
   do ict = 1, nct
      if (txcopolymer(ict) == 'block') then
         do icloc = 1, ncct(ict)                               ! loop over chains of type ict
            ic = ic+1                                          ! global chain number
            iseg = 0                                           ! initiate segment counter
            do ipt = 1, npt                                    ! loop over particle types
               iplow = sum(nppt(1:ipt-1)) + sum(ncct(1:ict-1)*npptct(ipt,1:ict-1)) + (icloc-1)*npptct(ipt,ict)
               do iploc = 1, npptct(ipt,ict)
                  iseg = iseg+1
                  ipnsegcn(iseg,ic) = iplow + iploc
               end do
            end do
         end do
      else if (txcopolymer(ict) == 'regular') then
         nrep = maxval(npptct(1:npt,ict))                      ! number of repetitions
         do ipt = 1, npt
            if (npptct(ipt,ict) > 0 .and. npptct(ipt,ict) < nrep) nrep = npptct(ipt,ict)
         end do
         nreplen = sum(npptct(1:npt,ict))/nrep                 ! length of a repetition (number of particles)
         do ipt = 1, npt
            if (mod(npptct(ipt,ict),nrep)/=0) call stop(txroutine,'error in npptct(ipt,ipt) for making regular copolymer', uout)
         end do
         do icloc = 1, ncct(ict)                               ! loop over chains of type ict
            ic = ic+1                                          ! global chain number
            iseg = 0                                           ! initiate segment counter
            do irep = 1, nrep                                  ! loop over number of repititions
               do ipt = 1, npt                                 ! loop over particle types
                  iplow = (irep-1)*npptct(ipt,ict)/nrep + sum(nppt(1:ipt-1)) + sum(ncct(1:ict-1)*npptct(ipt,1:ict-1)) + (icloc-1)*npptct(ipt,ict)
                  do iploc = 1, npptct(ipt,ict)/nrep
                     iseg = iseg+1
                     ipnsegcn(iseg,ic) = iplow + iploc
                  end do
               end do
            end do
         end do
      else if (txcopolymer(ict) == 'repeating') then
         if(.not. allocated(npset)) allocate(npset(npt))
         if(.not. allocated(ipstart)) allocate(ipstart(npt))
         npset = 0
         ipstart = 0
         if(any( rep_iblock_ict(1:nblockict(ict),ict)%np .le. 0 ) ) call stop(txroutine,'block of 0 length in repetition', uout)
         if(any( rep_iblock_ict(1:nblockict(ict),ict)%pt .le. 0 ) ) call stop(txroutine,'block without pt in repetition', uout)

         do icloc = 1, ncct(ict)                               ! loop over chains of type ict
            ic = ic + 1
            !repeating structure
            do ipt = 1, npt
               ipstart(ipt) = sum(nppt(1:ipt-1)) + sum(ncct(1:ict-1)*npptct(ipt,1:ict-1)) + (icloc - 1)*npptct(ipt,ict)
            end do
            iseg = 0
            npset = 0
            do while (iseg < sum(npptct(1:npt,ict)))
               do iblock = 1, nblockict(ict)
                  ipt = rep_iblock_ict(iblock,ict)%pt
                  do iploc = 1, min(rep_iblock_ict(iblock,ict)%np , npptct(ipt,ict) - npset(ipt))
                     iseg = iseg + 1
                     npset(ipt) = npset(ipt) + 1
                     ipnsegcn(iseg,ic) = npset(ipt) + ipstart(ipt)
                  end do
               end do
            end do

         end do

         deallocate(npset, ipstart)

      else if (txcopolymer(ict) == 'random') then


         !prepare allocatable variables
         if(.not. allocated(ipstart)) allocate(ipstart(npt))
         allocate(iptiseg(npct(ict)))

         !loob over chains
         do icloc = 1, ncct(ict)
            ic = ic+1                                          ! global chain number

            !create fresh list of iptiseg and iplowipt
            iseg = 0
            do ipt = 1, npt
               iptiseg((iseg+1):(iseg+npptct(ipt,ict))) = ipt
               ipstart(ipt) = sum(nppt(1:ipt-1)) + sum(ncct(1:ict-1)*npptct(ipt,1:ict-1)) + (icloc-1)*npptct(ipt,ict)
               iseg = iseg+npptct(ipt,ict)
            end do

            !shuffle iptiseg
            call KnuthShuffle(iptiseg, size(iptiseg), iseed)

            !assign particles
            do iseg = 1, npct(ict)
               ipt = iptiseg(iseg)
               ipstart(ipt) = ipstart(ipt) + 1
               ipnsegcn(iseg,ic) = ipstart(ipt)
            end do

         end do

         deallocate(iptiseg, ipstart)

      else if (txcopolymer(ict) == 'sequence') then

         ! ... prepare allocatable var
         if(.not. allocated(ipstart))  then
            allocate(ipstart(npt))
            ipstart = 0
         end if

         ! ... loop over chains of type ict
         do icloc = 1, ncct(ict)
            ic = ic+1                            ! global chain number
            do ipt = 1, npt
               ipstart(ipt) = sum(nppt(1:ipt-1)) + sum(ncct(1:ict-1)*npptct(ipt,1:ict-1)) + (icloc-1)*npptct(ipt,ict)
            end do
            do iseg = 1, npct(ict)               ! loop over chain segments
               ipt = iptsegct(iseg,ict)          ! particle type of segment iseg
               ipstart(ipt) = ipstart(ipt) + 1   ! increment ipstart
               ipnsegcn(iseg,ic) = ipstart(ipt)  ! set ipnsegcn
            end do
         end do

         deallocate(ipstart)

      end if
    end do

end subroutine Set_ipnsegcn

!........................................................................

subroutine Set_ictcn  ! chain -> its chain type
   if (.not.allocated(ictcn)) then
      allocate(ictcn(nc))
      ictcn = 0
   end if
   ic = 0
   do ict = 1, nct
      do icloc = 1, ncct(ict)
         ic = ic+1
         ictcn(ic) = ict
      end do
   end do
end subroutine Set_ictcn

!........................................................................

subroutine Set_ictpt  ! partcle type -> its chain type
   if (.not.allocated(ictpt)) then
      allocate(ictpt(npt))
      ictpt = 0
   end if
   ictpt(1:npt) = 0
   do ict = 1, nct
      do ipt = 1, npt
         if (npptct(ipt,ict) > 0) ictpt(ipt) = ict
      end do
   end do
end subroutine Set_ictpt

!........................................................................

subroutine Set_ihnpn  ! particle -> its hierarchical structure
   integer(4) :: ih, igen, ict, ic, iseg, ip
   if (.not.allocated(ihnpn)) then
      allocate(ihnpn(np_alloc))
      ihnpn = 0
   end if
   if (.not.allocated(icihigen)) then
      allocate(icihigen(nh,0:ngen))
      icihigen = 0
   end if
   do ih = 1, nh
      do igen = 0, ngen
         ict = ictgen(igen)
         nch(igen) = ncct(ict)/nh                  ! number of chains of generation igen in a h structure
         icihigen(ih,igen) = icnct(ict) + nch(igen)*(ih-1)
         do ic = icihigen(ih,igen), icihigen(ih,igen) + nch(igen) -1
            do iseg = 1, npct(ict)
               ip = ipnsegcn(iseg,ic)
               ihnpn(ip) = ih                      ! particle -> its hierarchical strcutre
            end do
         end do
      end do
   end do
end subroutine Set_ihnpn

!........................................................................

subroutine Set_ictpn  ! particle -> its chain type
   if (.not.allocated(ictpn)) then
      allocate(ictpn(np_alloc))
      ictpn = 0
   end if
   ictpn(1:np) = 0
   ic = 0
   do ict = 1, nct
      do icloc = 1, ncct(ict)
         ic = ic+1
         do iseg = 1, npct(ict)
            ictpn(ipnsegcn(iseg,ic)) = ict
         end do
      end do
   end do
end subroutine Set_ictpn

!........................................................................

subroutine Set_icnpn  ! particle -> its chain
   if (.not.allocated(icnpn)) then
      allocate(icnpn(np_alloc))
      icnpn = 0
   end if
   icnpn(1:np) = 0
   ic = 0
   do ict = 1, nct
      do icloc = 1, ncct(ict)
         ic = ic+1
         do iseg = 1, npct(ict)
            icnpn(ipnsegcn(iseg,ic)) = ic
         end do
      end do
   end do
end subroutine Set_icnpn

!........................................................................

subroutine Set_iptpn  ! particle -> its particle type
   if (.not.allocated(iptpn)) then
      allocate(iptpn(np_alloc))
      iptpn = 0
   end if
   ip = 0
   do ipt = 1, npt
      do iploc = 1, nppt(ipt)
         ip = ip+1
         iptpn(ip) = ipt
      end do
   end do
end subroutine Set_iptpn

!........................................................................

subroutine Set_iptat  ! atom type -> its particle type
   if (.not.allocated(iptat)) then
      allocate(iptat(nat))
      iptat = 0
   end if
   iat = 0
   do ipt = 1, npt
      do iatloc = 1, natpt(ipt)
         iat = iat+1
         iptat(iat) = ipt
      end do
   end do
end subroutine Set_iptat

!........................................................................

subroutine Set_iptan  ! atom -> its particle type
   if (.not.allocated(iptan)) then
      allocate(iptan(na_alloc))
      iptan = 0
   end if
   ia = 0
   do ipt = 1, npt
      do ialoc = 1, napt(ipt)*nppt(ipt)
         ia = ia+1
         iptan(ia) = ipt
      end do
   end do
end subroutine Set_iptan

!........................................................................

subroutine Set_ipnan  ! atom -> its particle
   if (.not.allocated(ipnan)) then
      allocate(ipnan(na_alloc))
      ipnan = 0
   end if
   ia = 0
   ip = 0
   ia = 0
   do ipt = 1, npt
      do iploc = 1, nppt(ipt)
         ip = ip+1
         do ialoc = 1, napt(ipt)
            ia = ia+1
            ipnan(ia) = ip
         end do
      end do
   end do
end subroutine Set_ipnan

!........................................................................

subroutine Set_iatan  ! atom -> its atom type
   integer(4)  :: ia0
   integer(4)  :: ia1
   integer(4)  :: ia2
   integer(4)  :: iaat
   if (.not.allocated(iatan)) then
      allocate(iatan(na_alloc))
      iatan = 0
   end if
   ia0 = 0
   iat = 0
   do ipt = 1, npt
      ia1 = ia0
      iaat = 0
      do iatloc = 1, natpt(ipt)
         ia2 = ia1
         iaat = iaat+1
         iat = iat+1
         do iploc = 1, nppt(ipt)
            ia = ia2
            do ialoc = 1, naatpt(iaat,ipt)
               ia = ia+1
               iatan(ia) = iat
            end do
            ia2 = ia2+napt(ipt)
         end do
         ia1 = ia1+naatpt(iaat,ipt)
      end do
      ia0 = ia0+napt(ipt)*nppt(ipt)
   end do
end subroutine Set_iatan

!........................................................................

subroutine Set_icnct  ! chain type -> its first chain
   if (.not.allocated(icnct)) then
      allocate(icnct(nct))
      icnct = 0
   end if
   ic = 1
   do ict = 1, nct
      icnct(ict) = ic
      ic = ic+ncct(ict)
   end do
end subroutine Set_icnct

!........................................................................

subroutine Set_ipnpt  ! particle type -> its first particle
   if (.not.allocated(ipnpt)) then
      allocate(ipnpt(npt))
      ipnpt = 0
   end if
   ip = 1
   do ipt = 1, npt
      ipnpt(ipt) = ip
      ip = ip+nppt(ipt)
   end do
end subroutine Set_ipnpt

!........................................................................

subroutine Set_iatpt  ! particle type -> its first atom type
   if (.not.allocated(iatpt)) then
      allocate(iatpt(npt))
      iatpt = 0
   end if
   iat = 1
   do ipt = 1, npt
      iatpt(ipt) = iat
      iat = iat+natpt(ipt)
   end do
end subroutine Set_iatpt

!........................................................................

subroutine Set_ianpn  ! particle -> its first atom
   if (.not.allocated(ianpn)) then
      allocate(ianpn(np_alloc))
      ianpn = 0
   end if
   ip = 0
   ia = 1
   do ipt = 1, npt
      do iploc = 1, nppt(ipt)
         ip = ip+1
         ianpn(ip) = ia
         ia = ia+napt(ipt)
      end do
   end do
end subroutine Set_ianpn

!........................................................................

subroutine Set_ianat  ! atom type -> its first atom
   if (.not.allocated(ianat)) then
      allocate(ianat(nat))
      ianat = 0
   end if
!   ia = 1
!   do iat = 1, nat
!      ianat(iat) = ia
!      ipt = iptat(iat)
!      ia = ia + nppt(ipt)*naatpt(iat,ipt)
!         write(*,*) iat, ipt, nppt(ipt), naatpt(iat,ipt)
!   end do
!   write(*,*) ianat(1:nat)

   iat = 0
   ia = 1
   do ipt = 1, npt
      do iatloc = 1, nat
         if (naatpt(iatloc,ipt) > 0) then
            iat = iat + 1
            ianat(iat) = ia
            ia = ia + nppt(ipt)*naatpt(iatloc,ipt)
         end if
      end do
   end do
!   write(*,*) ianat(1:nat); stop 55

end subroutine Set_ianat

!........................................................................

subroutine Set_ictct  ! two chain types -> chain type pair
   integer(4)  :: jct
   integer(4)  :: ictjct
   if (.not.allocated(ictct)) then
      allocate(ictct(nct,nct))
      ictct = 0
   end if
   ictjct = 0
   do ict = 1, nct
      do jct = ict, nct
         ictjct = ictjct+1
         ictct(ict,jct) = ictjct
         ictct(jct,ict) = ictjct
      end do
   end do
end subroutine Set_ictct

!........................................................................

subroutine Set_iptpt  ! two particle types -> particle type pair
   integer(4)  :: jpt
   integer(4)  :: iptjpt
   if (.not.allocated(iptpt)) then
      allocate(iptpt(npt,npt))
      iptpt = 0
   end if
   iptjpt = 0
   do ipt = 1, npt
      do jpt = ipt, npt
         iptjpt = iptjpt+1
         iptpt(ipt,jpt) = iptjpt
         iptpt(jpt,ipt) = iptjpt
      end do
   end do
end subroutine Set_iptpt

!........................................................................

subroutine Set_iatat  ! two atom types -> atom type pair
   integer(4)  :: jat
   integer(4)  :: iatjat
   if (.not.allocated(iatat)) then
      allocate(iatat(nat,nat))
      iatat = 0
   end if
   iatjat = 0
   do iat = 1, nat
      do jat = iat, nat
         iatjat = iatjat+1
         iatat(iat,jat) = iatjat
         iatat(jat,iat) = iatjat
      end do
   end do
end subroutine Set_iatat

!........................................................................

subroutine Set_isegpn   ! particle number -> segment number
   if (.not.allocated(isegpn)) then
      allocate(isegpn(np_alloc))
      isegpn = 0
   end if
   isegpn(1:np) = 0
   do ic = 1, nc
      ict = ictcn(ic)
      do iseg = 1, npct(ict)
         ip = ipnsegcn(iseg,ic)
         isegpn(ip) = iseg
      end do
   end do
end subroutine Set_isegpn

!........................................................................

subroutine Set_bondnn   ! bond and particle -> bonded particle
   if (.not.allocated(bondnn)) then
      allocate(bondnn(2,np))
      bondnn = 0
   end if
   do ic = 1, nc
      ict = ictcn(ic)
      iseg = 1                                     ! first segment
      ip = ipnsegcn(iseg,ic)
      if (npct(ict) == 1) then
         bondnn(1,ip) = 0                          ! first segment no bond below
         bondnn(2,ip) = 0                          ! end segment no bond uppwards
      else
         bondnn(1,ip) = 0                          ! first segment no bond below
         bondnn(2,ip) = ipnsegcn(iseg+1,ic)        ! first segment bonds upwards
         do iseg = 2, npct(ict)-1                  ! middle segment
            ip = ipnsegcn(iseg,ic)
            bondnn(1,ip) = ipnsegcn(iseg-1,ic)     ! middle segment bonds downwards
            bondnn(2,ip) = ipnsegcn(iseg+1,ic)     ! middle segment bonds upwards
         end do
         iseg = npct(ict)                          ! end segment
         ip = ipnsegcn(iseg,ic)
         bondnn(1,ip) = ipnsegcn(iseg-1,ic)        ! end segment bonds downwards
         bondnn(2,ip) = 0                          ! end segment no bond uppwards
      end if
   end do
end subroutine Set_bondnn

!........................................................................

subroutine Set_nnw ! number of networks
   nnw = sum(nnwnwt(1:nnwt))
end subroutine Set_nnw

!........................................................................

subroutine Set_nclnwt ! number of cross-links per network of network type inwt
   do inwt = 1, nnwt
      nclnwt(inwt) = nppt(iptclnwt(inwt))/nnwnwt(inwt)
   end do
end subroutine Set_nclnwt

!........................................................................

subroutine Set_lptnwt ! particle type ipt used for network type inwt?
   lptnwt = .false.
   do ipt = 1, npt
      do inwt = 1, nnwt
         if (ipt == iptclnwt(inwt)) then ! ... ipt forms the nodes of networks of type inwt
            lptnwt(ipt,inwt) = .true.
            exit
         else if (ictpt(ipt) > 0) then ! ... ipt forms chains ...
            if (ncctnwt(ictpt(ipt),inwt) > 0) lptnwt(ipt,inwt) = .true. ! ... and is part of networks of type inwt
         end if
      end do
   end do
end subroutine Set_lptnwt

!........................................................................

subroutine Set_nctnwt ! number of chain types of network type inwt
   do inwt = 1, nnwt
      nctnwt(inwt) = count(ncctnwt(1:nct,inwt) > 0)
   end do
end subroutine Set_nctnwt

!........................................................................

subroutine Set_ncnwt ! number of chains of network type inwt
   do inwt = 1, nnwt
      ncnwt(inwt) = sum(ncctnwt(1:nct,inwt))
   end do
end subroutine Set_ncnwt

!........................................................................

subroutine Set_npnwt ! number of particles of network type inwt
   do inwt = 1, nnwt
      npnwt(inwt) = nclnwt(inwt)
      do ict = 1, nct
         npnwt(inwt) = npnwt(inwt) + npct(ict)*ncctnwt(ict,inwt)
      end do
   end do
end subroutine Set_npnwt

!........................................................................

subroutine Set_nnwtnwt ! number of network type pairs
   nnwtnwt = (nnwt*(nnwt+1))/2
end subroutine Set_nnwtnwt

!........................................................................

subroutine Set_npweakchargenwt ! number of titratable particles in networks of type inwt
   if (lweakcharge) then
      do inwt = 1, nnwt
         npweakchargenwt(inwt) = Zero
         do ipt = 1, npt
            if ((.not. lptnwt(ipt,inwt)) .or. (.not. latweakcharge(iatpt(ipt)))) cycle
            npweakchargenwt(inwt) = npweakchargenwt(inwt) + nppt(ipt)/nnwnwt(inwt)
         end do
      end do
   end if
end subroutine Set_npweakchargenwt

!........................................................................

subroutine Set_inwtnwn ! network -> its network type
   if (.not.allocated(inwtnwn)) then
      allocate(inwtnwn(nnw))
      inwtnwn = 0
   end if
   inw = 0
   do inwt = 1, nnwt
      do inwloc = 1, nnwnwt(inwt)
         inw = inw+1
         inwtnwn(inw) = inwt
      end do
   end do
end subroutine Set_inwtnwn

!........................................................................

subroutine Set_inwtct  ! chain type -> its network type
   if (.not.allocated(inwtct)) then
      allocate(inwtct(nct))
      inwtct = 0
   end if
   inwtct(1:nct) = 0
   do inwt = 1, nnwt
      where (ncctnwt(1:nct,inwt) > 0) inwtct(1:nct) = inwt
   end do
end subroutine Set_inwtct

!........................................................................

subroutine Set_icnclocnwn ! local chain and network -> its chain
   character(40), parameter :: txroutine = 'Set_icnclocnwn'
   integer(4)               :: iclow
   integer(4)               :: icctloc
   if (.not.allocated(icnclocnwn)) then
      allocate(icnclocnwn(maxval(ncnwt(1:nnwt)),nnw))
      icnclocnwn = 0
   end if
   inw = 0
   do inwt = 1, nnwt
      if (txtoponwt(inwt) == 'default') then   ! currently only txtoponwt = 'default' is supported
         do inwloc = 1, nnwnwt(inwt)
            inw = inw+1
            icloc = 0
            do ict = 1, nct
               iclow = sum(ncct(1:ict-1)) + sum(nnwnwt(1:inwt-1)*ncctnwt(ict,1:inwt-1)) + (inwloc-1)*ncctnwt(ict,inwt)
               do icctloc = 1, ncctnwt(ict,inwt)
                  icloc = icloc+1
                  icnclocnwn(icloc,inw) = iclow + icctloc
               end do
            end do
         end do
      else
         call Stop(txroutine, 'txtoponwt /= "default"', uout)
      end if
   end do
end subroutine Set_icnclocnwn

!........................................................................

subroutine Set_inwtcn  ! chain -> its network type
   if (.not.allocated(inwtcn)) then
      allocate(inwtcn(nc))
      inwtcn = 0
   end if
   inwtcn(1:nc) = 0
   inw = 0
   do inwt = 1, nnwt
      do inwloc = 1, nnwnwt(inwt)
         inw = inw+1
         do icloc = 1, ncnwt(inwt)
            inwtcn(icnclocnwn(icloc,inw)) = inwt
         end do
      end do
   end do
end subroutine Set_inwtcn

!........................................................................

subroutine Set_inwncn ! chain -> its network
   if (.not.allocated(inwncn)) then
      allocate(inwncn(nc))
      inwncn = 0
   end if
   inwncn(1:nc) = 0
   inw = 0
   do inwt = 1, nnwt
      do inwloc = 1, nnwnwt(inwt)
         inw = inw+1
         do icloc = 1, ncnwt(inwt)
            inwncn(icnclocnwn(icloc,inw)) = inw
         end do
      end do
   end do
end subroutine Set_inwncn

!........................................................................

subroutine Set_inwnnwt ! network type -> its first network
   if (.not.allocated(inwnnwt)) then
      allocate(inwnnwt(nnwt))
      inwnnwt = 0
   end if
   inw = 1
   do inwt = 1, nnwt
      inwnnwt(inwt) = inw
      inw = inw + nnwnwt(inwt)
   end do
end subroutine Set_inwnnwt

!........................................................................

subroutine Set_inwtnwt  ! two network types -> network type pair
   integer(4)  :: inwtjnwt
   integer(4)  :: jnwt
   if (.not.allocated(inwtnwt)) then
      allocate(inwtnwt(nnwt,nnwt))
      inwtnwt = 0
   end if
   inwtjnwt = 0
   do inwt = 1, nnwt
      do jnwt = inwt, nnwt
         inwtjnwt = inwtjnwt+1
         inwtnwt(inwt,jnwt) = inwtjnwt
         inwtnwt(jnwt,inwt) = inwtjnwt
      end do
   end do
end subroutine Set_inwtnwt

!........................................................................

subroutine Set_ipncllocnwn ! cross-link and network -> its particle
   if (.not.allocated(ipncllocnwn)) then
      allocate(ipncllocnwn(maxval(nclnwt(1:nnwt)),nnw))
      ipncllocnwn = 0
   end if
   inw = 0
   do inwt = 1, nnwt
      ip = ipnpt(iptclnwt(inwt))-1
      do inwloc = 1, nnwnwt(inwt)
         inw = inw+1
         do iclloc = 1, nclnwt(inwt)
            ip = ip+1
            ipncllocnwn(iclloc,inw) = ip
         end do
      end do
   end do
end subroutine Set_ipncllocnwn

!........................................................................

subroutine Set_ipnplocnwn ! local particle and network -> its particle
   if (.not.allocated(ipnplocnwn)) then
      allocate(ipnplocnwn(maxval(npnwt(1:nnwt)),nnw))
      ipnplocnwn = 0
   end if
   inw = 0
   do inwt = 1, nnwt
      do inwloc = 1, nnwnwt(inwt)
         inw = inw+1
         iploc = 0
         do iclloc = 1, nclnwt(inwt)
            iploc = iploc+1
            ip = ipncllocnwn(iclloc,inw)
            ipnplocnwn(iploc,inw) = ip
         end do
         do icloc = 1, ncnwt(inwt)
            ic = icnclocnwn(icloc,inw)
            do iseg = 1, npct(ictcn(ic))
               iploc = iploc+1
               ip = ipnsegcn(iseg,ic)
               ipnplocnwn(iploc,inw) = ip
            end do
         end do
      end do
   end do
end subroutine Set_ipnplocnwn

!........................................................................

subroutine Set_lpnnwn ! particle part of network?
   if (.not.allocated(lpnnwn)) then
      allocate(lpnnwn(np,nnw))
      lpnnwn = .false.
   end if
   do inw = 1, nnw
      do ip = 1, np
         if (icnpn(ip) > 0) then
            if (inwncn(icnpn(ip)) == inw) lpnnwn(ip,inw) =.true.
         end if
      end do
      inwt = inwtnwn(inw)
      do iclloc = 1, nclnwt(inwt)
         ip = ipncllocnwn(iclloc,inw)
         lpnnwn(ip,inw) =.true.
      end do
   end do
end subroutine Set_lpnnwn

!........................................................................

subroutine Set_ipnhn  ! hierarchical strcture -> its first particle
   ! get particle with lowest number in hierarchical strucutre
   ipnhn = minval(ipnsegcn(1,icnct(ictgen(0:ngen))))

   !old code gets number of first particle in strucutre:
      !igen = 0
      !ict = ictgen(igen)
      !ic = icnct(ict)
      !iseg = 1
      !ipnhn = ipnsegcn(iseg,ic)
end subroutine Set_ipnhn

!........................................................................

subroutine Set_genic  ! chain number -> generation number
   if (.not.allocated(genic)) then
      allocate(genic(nc))
      genic = 0
   end if
   do igen = 0, ngen                                                    ! loop over generations
      ict = ictgen(igen)                                                ! chain type
      do ic = icnct(ict), icnct(ict)+ncct(ict)-1                        ! loop over chains
         genic(ic) = igen
      end do
   end do
end subroutine Set_genic

!........................................................................

subroutine Set_bondcl   ! crosslink and particle -> crosslinked particle
   character(40), parameter :: txroutine ='Set_bondcl'
   integer(4) :: iseg, ic, ict, ip, ipt, jp, jpt, isegmc, icmc, icmcloc, igen, icscloc, isidechain
   integer(4) :: jploc

   maxvalnbondcl = maxval(maxnbondcl(1:npt))
   if (.not.allocated(nbondcl)) then
      allocate(nbondcl(np_alloc))
      nbondcl = 0
   end if
   if (.not.allocated(bondcl)) then
      allocate(bondcl(maxvalnbondcl,np_alloc))
      bondcl = 0
   end if

   nbondcl = 0
   do igen = 0, ngen                                                    ! loop over generations
      ict = ictgen(igen)                                                ! chain type
      do ic = icnct(ict), icnct(ict)+ncct(ict)-1                        ! loop over chains
         do iseg = 1, npct(ict)                                         ! loop over segments
            ip = ipnsegcn(iseg,ic)                                      ! particle to be set
            ipt = iptpn(ip)                                             ! particle type
            if (iseg == 1) then                                          ! a first segment?
               if (igen > 0) then                                        ! not zeroth generation?
                  icscloc = ic-(icnct(ictgen(igen)) - 1)                ! local side-chain id
                  icmcloc = (icscloc-1)/nbranch(igen-1) + 1             ! local main-chain id
                  icmc = icnct(ictgen(igen-1)) + (icmcloc-1)            ! main-chain id
                  isidechain = mod(icscloc-1,nbranch(igen-1)) + 1       ! order number of the side chain of main chain
                                                                        ! main-chain segment to form a crosslink
                  if (ibranchpinc(igen-1) >= 0) then                    ! fixed increment (for placing the side chains)
                     isegmc = ibranchpbeg(igen-1) + ibranchpinc(igen-1)*(isidechain-1)
                  else    ! fixed length of block of main-chain segments with side chains, separated by one main-chain segment without side chain
                     isegmc = ibranchpbeg(igen-1)+isidechain-1-(isidechain-1)/ibranchpinc(igen-1)
                  end if

                  jp = ipnsegcn(isegmc,icmc)                            ! id of particle to which ip will form a crosslink
                  jpt = iptpn(jp)                                             ! particle type
                  if ((nbondcl(ip) < maxnbondcl(ipt)) .and. (nbondcl(jp) < maxnbondcl(jpt))) then    ! below the upper limit of crosslinks?
                     nbondcl(ip) = nbondcl(ip) + 1
                     nbondcl(jp) = nbondcl(jp) + 1
                     if (nbondcl(ip) > maxvalnbondcl) call Stop (txroutine, 'nbondcl > maxvalnbondcl', uout)
                     if (nbondcl(jp) > maxvalnbondcl) call Stop (txroutine, 'nbondcl > maxvalnbondcl', uout)
                     bondcl(nbondcl(ip),ip) = jp
                     bondcl(nbondcl(jp),jp) = ip
                     if(lmultigraft) then    !   shift all grafting point to the backbone
!             write(*,*) ip, jp, nbondcl(ip), nbondcl(jp), bondcl(nbondcl(ip),ip), bondcl(nbondcl(jp),jp)
!             write(*,*) bondcl(nbondcl(ip)-1,ip), bondcl(nbondcl(jp)-1,jp)
                         if(nbondcl(jp)==2) then
!                write(*,*) ip, jp, bondcl(nbondcl(ip),ip), bondcl(nbondcl(jp),jp), bondcl(nbondcl(jp)-1,jp)
                            jploc = bondcl(nbondcl(jp)-1,jp)               ! jploc
                            bondcl(nbondcl(ip),ip) = jploc                 ! bind ip to jploc
                            nbondcl(jploc) = nbondcl(jploc) + 1            ! update jploc
                            bondcl(nbondcl(jploc),jploc) = ip
                            nbondcl(jp) = nbondcl(jp) - 1                  ! update jp
                         end if
                      end if
                  else
                     call Warn(txroutine, 'maxnbondcl hitted', 6)
                     write(*,*) 'ip, nbondcl(ip)', ip, nbondcl(ip)
                     write(*,*) 'jp, nbondcl(jp)', jp, nbondcl(jp)
                  end if
               end if
            end if
         end do
      end do
   end do

end subroutine Set_bondcl

!........................................................................

subroutine TestChainPointer(unit)
   integer(4) :: ic, iseg, ip, ict, jct
   integer(4),   intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//' chain', unit)
   write(unit,'(a)') 'ict, jct, ictct'
   write(unit,'(3i5)') ((ict, jct, ictct(ict, jct), jct = ict, nct), ict = 1, nct)
   write(unit,'()')
   write(unit,'(a)') '        ic, ictcn(ic),      iseg,   ipnsegcn(iseg,ic) = i'
   write(unit,'(4i11)') ((ic,ictcn(ic),iseg,ipnsegcn(iseg,ic),iseg = 1,npct(ictcn(ic))),ic = 1,nc)
   write(unit,'()')
   write(unit,'(a)') '        ip, iptpn(ip), ictpn(ip), icnpn(ip), bondnn(1,ip), bondnn(2,ip), isegpn(ip)'
   write(unit,'(7i11)') (ip, iptpn(ip), ictpn(ip), icnpn(ip), bondnn(1,ip), bondnn(2,ip), isegpn(ip) ,ip = 1,np)

end subroutine TestChainPointer

!........................................................................

subroutine TestNetworkPointer(unit)
   integer(4),   intent(in) :: unit
   integer(4)               :: inw, inwt, jnwt, icloc, ic, iclloc
   character(len=3)         :: txfmt
   call WriteHead(3, 'Test'//trim(txroutine)//' network', unit)
   write(unit,'(a)') 'inwt, jnwt, inwtnwt'
   write(unit,'(3i5)') ((inwt, jnwt, inwtnwt(inwt, jnwt), jnwt = inwt, nnwt), inwt = 1, nnwt)
   write(txfmt,'(i3)') nnwt
   write(unit,'(a,'//trim(adjustl(txfmt))//'i5)') 'inwnnwt(1:nnwt)', inwnnwt(1:nnwt)
   write(txfmt,'(i3)') nct
   write(unit,'(a,'//trim(adjustl(txfmt))//'i5)') 'inwtct(1:nct)', inwtct(1:nct)
   write(unit,'()')
   write(unit,'(a)') '        inw, inwtnwn(inw),      icloc,   icnclocnwn(icloc,inw) = i'
   write(unit,'(4i11)') ((inw,inwtnwn(inw),icloc,icnclocnwn(icloc,inw),icloc = 1,ncnwt(inwtnwn(inw))),inw = 1,nnw)
   write(unit,'()')
   write(unit,'(a)') '        ic, ictcn(ic), inwtcn(ic), inwncn(ic)'
   write(unit,'(4i11)') (ic, ictcn(ic), inwtcn(ic), inwncn(ic),ic = 1,nc)
   write(unit,'()')
   write(unit,'(a)') '       iclloc,          inw,    iptclnwt, ipncllocnwn(iclloc,inw)'
   write(unit,'(4i13)') ((iclloc, inw, iptclnwt(inwtnwn(inw)), ipncllocnwn(iclloc,inw), iclloc = 1,nclnwt(inwtnwn(inw))), inw = 1, nnw)
   write(unit,'()')
   write(unit,'(a)') '        ipt,       inwt,       lptnwt(ipt,inwt)'
   write(unit,'(2i11,12x,l)') ((ipt, inwt, lptnwt(ipt,inwt), ipt = 1,npt), inwt = 1, nnwt)
   write(unit,'()')
   write(unit,'(a)') '       ip,      inw,       lpnnwn(ip,inw)'
   write(unit,'(2i11,12x,l)') ((ip, inw, lpnnwn(ip,inw), ip = 1,np), inw = 1, nnw)
   write(unit,'()')
   write(unit,'(a)') '       iploc,         inw,        inwt,  ipnplocnwn'
   write(unit,'(4i12)') ((iploc, inw, inwtnwn(inw), ipnplocnwn(iploc,inw), iploc = 1, npnwt(inwtnwn(inw))), inw = 1, nnw)
end subroutine TestNetworkPointer

!........................................................................

subroutine TestCrosslinkPointer(unit)
   integer(4),   intent(in) :: unit
   integer(4) :: ibond, ip
   call WriteHead(3, 'Test'//trim(txroutine)//' crosslink', unit)
   write(unit,'(a,40i5)') 'nch(0:ngen)', nch(0:ngen)
   write(unit,'()')
   write(unit,'(a)') 'ip, ihnpn(ip)'
   do ip = 1, np
      write(unit,'(2i5)') ip, ihnpn(ip)
   end do
   write(unit,'()')
   write(unit,'(a)') 'ip    ic           maxnbondcl(ipt)         nbondcl(ip)          bondcl(1,ip)  ...  '
   do ip = 1, np
      write(unit,'(2i5,2i18,10x,4i10)') ip, icnpn(ip), maxnbondcl(iptpn(ip)), nbondcl(ip), (bondcl(ibond,ip), ibond = 1, nbondcl(ip))
   end do
   if (sum(nbondcl(1:np))/2 /= ncl) call Stop('Test'//trim(txroutine), 'error in the number of crosslinks made', unit)
end subroutine TestCrosslinkPointer

!........................................................................

end subroutine SetObjectParam1

!************************************************************************
!> \page particle particle.F90
!!  **SetObjectParam2**
!! *set object parameters*
!************************************************************************

subroutine SetObjectParam2

   use ParticleModule
   implicit none

   character(40), parameter :: txroutine ='SetObjectParam2'
   integer(4) :: inwt, jnwt, ict, jct, ip, ipt, jpt, ia, iat, jat, iatjat, iatloc, ntemp, ja
   real(8)    :: r2

! ... set txnwtnwt, txctct, txptpt, and txatat

   if (.not.allocated(txnwtnwt)) then
      allocate(txnwtnwt(nnwtnwt))
      txnwtnwt = ''
   end if
   if (.not.allocated(txctct)) then
      allocate(txctct(nctct))
      txctct = ''
   end if
   if (.not.allocated(txptpt)) then
      allocate(txptpt(nptpt))
      txptpt = ''
   end if
   if (.not.allocated(txatat)) then
      allocate(txatat(natat))
      txatat = ''
   end if

   do inwt = 1, nnwt
      do jnwt = inwt, nnwt
         txnwtnwt(inwtnwt(inwt,jnwt)) = trim(txnwt(inwt))//'-'//trim(txnwt(jnwt))
      end do
   end do
   do ict = 1, nct
       do jct = ict, nct
          txctct(ictct(ict,jct)) = trim(txct(ict))//'-'//trim(txct(jct))
       end do
    end do
    do ipt = 1, npt
       do jpt = ipt, npt
          txptpt(iptpt(ipt,jpt)) = trim(txpt(ipt))//'-'//trim(txpt(jpt))
       end do
    end do
    do iat = 1, nat
       do jat = iat, nat
          txatat(iatat(iat,jat)) = trim(txat(iat))//'-'//trim(txat(jat))
       end do
    end do

! ... set az

    if (.not.allocated(az)) then
       allocate(az(na_alloc))
       az = 0.0E+00
    end if
    do ia = 1, na
       iat     = iatan(ia)
       az(ia)  = zat(iat)
    end do

! ... weak charge: allocate memory and set pointer iananweakcharge

   if (lweakcharge) then

! ... allocate and initialize weak charge related parameters

      if(.not.allocated(laz)) then
         allocate(laz(na_alloc))
         laz = .false.
      end if
      if (.not. allocated(iananweakcharge)) then
         allocate(iananweakcharge(na_alloc))
         iananweakcharge = 0
      end if
      pHmpK(1:nat) = pH - pK(1:nat)

! ... test whether specified number of counterions matches the corresponding number of ionizable units

      if (any(jatweakcharge(1:nat) /= 0)) then ! any counterions specified?
         do jat = 1, nat
            if (any(jatweakcharge(1:nat) == jat)) then ! jat is a counterion

! ... check some conditions

               if (naat(jat) /= 1) call Stop(txroutine,'naat(jat) /= 1', uout) ! particle with counter charge has to be an ion

! ... sum up number of ionizable units for which jat shall be used as counterion type

               ntemp = 0
               do iat = 1, nat
                  if (iat == jat) cycle
                  if (latweakcharge(iat) .and. (jatweakcharge(iat) == jat)) then ! iat is the corresponding ionizable atom type
                     if ((zat(iat)+zat(jat)) /= Zero) call Stop(txroutine,'zat(iat) /= -zat(jat)',uout)
                     ntemp = ntemp + (naat(iat)*nppt(iptat(iat)))
                  end if
               end do

! ... test whether numbers match

               if (ntemp /= (naat(jat)*nppt(iptat(jat)))) then
                  call Stop(txroutine,'inconsistency among number of titrating units and counterions', uout)
               end if
            end if
         end do

! ... assign every ionizable unit a counterion (store in iananweakcharge)

         do jat = 1, nat
            if (any(jatweakcharge(1:nat) == jat)) then
               ja = ianat(jat) - 1
               do ia = 1, na
                  if (latweakcharge(iatan(ia)) .and. jatweakcharge(iatan(ia)) == jat) then
                     ja = ja + 1
                     iananweakcharge(ia) = ja
                  end if
               end do
            end if
         end do
      end if
   end if

! ... calculate particle masses, masspt, and massipt

    if (.not.allocated(masspt)) then
       allocate(masspt(npt), massipt(npt))
       masspt = 0.0E+00
       massipt = 0.0E+00
    end if
    iat = 0
    do ipt = 1, npt
       masspt(ipt) = Zero
       do iatloc = 1, natpt(ipt)
          iat = iat+1
          masspt(ipt) = masspt(ipt)+naatpt(iatloc,ipt)*massat(iat)
       end do
       if (masspt(ipt) > Zero) massipt(ipt) = One/masspt(ipt)
    end do

! ... calculate network mass massnwt and inverse mass massinwt

   if (lnetwork) then
      if (.not.allocated(massnwt)) then
         allocate(massnwt(nnwt), massinwt(nnwt))
         massnwt = 0.0E+00
         massinwt = 0.0E+00
      end if
      do inwt = 1, nnwt
         massnwt(inwt)  = sum(masspt(iptpn(1:np)), MASK=lpnnwn(1:np,inwnnwt(inwt)))
         if (massnwt(inwt) == Zero) then
            call Warn(txroutine,'NetworkProperties flawed for massnwt(inwt) == Zero',uout)
         end if
         if (massnwt(inwt) > Zero) massinwt(inwt) = One/massnwt(inwt)
      end do
   end if

! ... calculate atom coordinates in the principal axes system and mompt and momipt

    if (.not.allocated(poltens)) then
       allocate(poltens(6,na_alloc))
       poltens = 0.0E+00
    end if
    ntemp = maxval(napt(1:npt))
    if (.not.allocated(ra)) then
       allocate(ra(3,ntemp,npt), dipa(3,ntemp,npt), poltensa(6,ntemp,npt), rasite(3,ntemp,npt))
       ra = 0.0E+00
       dipa = 0.0E+00
       poltensa = 0.0E+00
       rasite = 0.0E+00
    end if
    if (.not.allocated(mompt)) then
       allocate(mompt(3,npt), momipt(3,npt))
       mompt = 0.0E+00
       momipt = 0.0E+00
    end if
    do ipt = 1, npt
       call PAxesSystem(ipt)
    end do

! ... for project with Gunnar

   if (txuser(1:9) == 'md_dipole') call md_dipole

! ... set lapolarization

    if (lpolarization) then
       if (.not.allocated(lapolarization)) then
          allocate(lapolarization(na_alloc))
          lapolarization = .false.
       end if
       lapolarization(1:na) =.false.
       do ia = 1, na
          ip = ipnan(ia)
          ipt = iptan(ia)
          iatloc = ia-(ianpn(ip)-1)
          if (sum(poltensa(1:3,iatloc,ipt)) > Zero) lapolarization(ia) =.true.
       end do
    end if

! ... calculate massp, massip, momp, and momip

    if (.not.allocated(massp)) then
       allocate(massp(np_alloc), massip(np_alloc), momp(3,np_alloc), momip(3,np_alloc))
       massp = 0.0E+00
       massip = 0.0E+00
       momp = 0.0E+00
       momip = 0.0E+00
    end if
    do ip = 1, np
       ipt = iptpn(ip)
       massp(ip) = masspt(ipt)
       massip(ip) = massipt(ipt)
       momp(1:3,ip) = mompt(1:3,ipt)
       momip(1:3,ip) = momipt(1:3,ipt)
    end do

! ... calculate the maximal atom-com distance, racom

    if (.not.allocated(racom)) then
       allocate(racom(npt))
       racom = 0.0E+00
    end if
    do ipt = 1, npt
       racom(ipt) = Zero
       do iatloc = 1, natpt(ipt)
          r2 = sqrt(ra(1,iatloc,ipt)**2+ra(2,iatloc,ipt)**2+ra(3,iatloc,ipt)**2)
          racom(ipt) = max(racom(ipt),r2)
       end do
    end do

! ... set r1atat and r2atat

    if (.not.allocated(r1atat)) then
       allocate(r1atat(natat), r2atat(natat))
       r1atat = 0.0E+00
       r2atat = 0.0E+00
    end if
    do iat = 1, nat
       do jat = iat, nat
          iatjat = iatat(iat,jat)
          r1atat(iatjat) = (radat(iat)+radat(jat))
          r2atat(iatjat) = (radat(iat)+radat(jat))**2
       end do
    end do

! ... calculate the number of translational and rotational degrees of freedom

    if (.not.allocated(itradegfree)) then
       allocate(itradegfree(1:npt), irotdegfree(1:npt))
       itradegfree = 0
       irotdegfree = 0
    end if
    do ipt = 1, npt
       itradegfree(ipt) = 3*nppt(ipt)
       irotdegfree(ipt) = 3*nppt(ipt)
       if (real(min(mompt(1,ipt), mompt(2,ipt),mompt(3,ipt))) < 1.0d-6) irotdegfree(ipt) = 2*nppt(ipt)
       if (real(max(mompt(1,ipt), mompt(2,ipt),mompt(3,ipt))) < 1.0d-6) irotdegfree(ipt) = 0
    end do

contains

!........................................................................

! set moments of inertia and lpolyatom = .true. for angular motion, project with Gunnar 2006-12

subroutine md_dipole
   mompt(1:3,1) = 1.0
   momipt(1:3,1) = One/mompt(1:3,1)
   lpolyatom = .true.
   lmonoatom = .false.
end subroutine md_dipole

!........................................................................

end subroutine SetObjectParam2

!************************************************************************
!> \page particle particle.F90
!!  **PAxesSystem**
!! *transform from input to principal axes system*
!************************************************************************

!     modified for use of expansion centers not situated at atoms/peo

subroutine PAxesSystem(ipt)

   use ParticleModule
   implicit none

   integer(4) :: ipt, ia, ialoc, ialow, iat, m, nrot
   real(8) :: rot(3), tpi(3,3), eivr(3,3), atemp(3,3), btemp(3,3), diagonal(3)
   real(8) :: InvFlt

   ialow = ianpn(ipnpt(ipt))

! ... calculate center of mass of particle type ipt

   rot = Zero
   ia = ialow-1
   do ialoc = 1, napt(ipt)
      ia = ia+1
      iat = iatan(ia)
      rot(1) = rot(1)+massat(iat)*rain(1,ialoc,ipt)
      rot(2) = rot(2)+massat(iat)*rain(2,ialoc,ipt)
      rot(3) = rot(3)+massat(iat)*rain(3,ialoc,ipt)
   end do
   rot(1:3) = rot(1:3)*massipt(ipt)

! ... set up the moment of inertia tensor

   tpi= Zero
   ia = ialow-1
   do ialoc = 1, napt(ipt)
      ia = ia+1
      iat = iatan(ia)
      rtm(1,ialoc) = rain(1,ialoc,ipt)-rot(1)
      rtm(2,ialoc) = rain(2,ialoc,ipt)-rot(2)
      rtm(3,ialoc) = rain(3,ialoc,ipt)-rot(3)
      tpi(1,1) = tpi(1,1)+massat(iat)*(rtm(2,ialoc)**2+rtm(3,ialoc)**2)
      tpi(1,2) = tpi(1,2)-massat(iat)*rtm(1,ialoc)*rtm(2,ialoc)
      tpi(1,3) = tpi(1,3)-massat(iat)*rtm(1,ialoc)*rtm(3,ialoc)
      tpi(2,2) = tpi(2,2)+massat(iat)*(rtm(1,ialoc)**2+rtm(3,ialoc)**2)
      tpi(2,3) = tpi(2,3)-massat(iat)*rtm(2,ialoc)*rtm(3,ialoc)
      tpi(3,3) = tpi(3,3)+massat(iat)*(rtm(1,ialoc)**2+rtm(2,ialoc)**2)
   end do
   tpi(2,1) = tpi(1,2)
   tpi(3,1) = tpi(1,3)
   tpi(3,2) = tpi(2,3)

!  write(*,'(a,3(3f22.18))') 'tpi', tpi(1:3,1:3)
   where(abs(tpi) < 1d-10) tpi = zero    ! security measure: treat accidentially nondiagonal as diagonal

! ... calculate moment of intertia tensor, mompt and  momipt

! ... eivr(m,i) m:th component of the i:th eigenvalue

   call Diag(3, tpi, diagonal, eivr, nrot)
   do m = 1, 3
      mompt(m,ipt) = diagonal(m)
      momipt(m,ipt) = InvFlt(mompt(m,ipt))
   end do

!   write(*,'(a,3(3f10.4))') 'diag', tpi(1,1), tpi(2,2), tpi(3,3)
!   write(*,'(a,3(3f10.4))') 'eivr', eivr(1:3,1:3)

! ... calculate the atom coordinates in the principal axes system
! ... (as well as site coordinates)

   do ialoc = 1, napt(ipt)
      ra(1,ialoc,ipt) = eivr(1,1)*rtm(1,ialoc)+eivr(2,1)*rtm(2,ialoc)+eivr(3,1)*rtm(3,ialoc)
      ra(2,ialoc,ipt) = eivr(1,2)*rtm(1,ialoc)+eivr(2,2)*rtm(2,ialoc)+eivr(3,2)*rtm(3,ialoc)
      ra(3,ialoc,ipt) = eivr(1,3)*rtm(1,ialoc)+eivr(2,3)*rtm(2,ialoc)+eivr(3,3)*rtm(3,ialoc)
      if (lintsite) then
         rasite(1,ialoc,ipt) = eivr(1,1)*(raintin(1,ialoc,ipt)-rot(1))          &
                              +eivr(2,1)*(raintin(2,ialoc,ipt)-rot(2))          &
                              +eivr(3,1)*(raintin(3,ialoc,ipt)-rot(3))
         rasite(2,ialoc,ipt) = eivr(1,2)*(raintin(1,ialoc,ipt)-rot(1))          &
                              +eivr(2,2)*(raintin(2,ialoc,ipt)-rot(2))          &
                              +eivr(3,2)*(raintin(3,ialoc,ipt)-rot(3))
         rasite(3,ialoc,ipt) = eivr(1,3)*(raintin(1,ialoc,ipt)-rot(1))          &
                              +eivr(2,3)*(raintin(2,ialoc,ipt)-rot(2))          &
                              +eivr(3,3)*(raintin(3,ialoc,ipt)-rot(3))
      end if
   end do

   if (ldipole .or. ldipolesph) then

! ... calculate the dipole moments in the principal axes system

      do ialoc = 1, napt(ipt)
         dipa(1,ialoc,ipt) = eivr(1,1)*dipain(1,ialoc,ipt)+eivr(2,1)*dipain(2,ialoc,ipt)+eivr(3,1)*dipain(3,ialoc,ipt)
         dipa(2,ialoc,ipt) = eivr(1,2)*dipain(1,ialoc,ipt)+eivr(2,2)*dipain(2,ialoc,ipt)+eivr(3,2)*dipain(3,ialoc,ipt)
         dipa(3,ialoc,ipt) = eivr(1,3)*dipain(1,ialoc,ipt)+eivr(2,3)*dipain(2,ialoc,ipt)+eivr(3,3)*dipain(3,ialoc,ipt)
      end do

   end if

   if (lpolarization) then

! ... calculate the dipole moments in the principal axes system

      do ialoc = 1, napt(ipt)
         dipa(1,ialoc,ipt) = eivr(1,1)*dipain(1,ialoc,ipt)+eivr(2,1)*dipain(2,ialoc,ipt)+eivr(3,1)*dipain(3,ialoc,ipt)
         dipa(2,ialoc,ipt) = eivr(1,2)*dipain(1,ialoc,ipt)+eivr(2,2)*dipain(2,ialoc,ipt)+eivr(3,2)*dipain(3,ialoc,ipt)
         dipa(3,ialoc,ipt) = eivr(1,3)*dipain(1,ialoc,ipt)+eivr(2,3)*dipain(2,ialoc,ipt)+eivr(3,3)*dipain(3,ialoc,ipt)
      end do

! ... calculate the polarizablity tensor in the principal axes system

      do ialoc = 1, napt(ipt)
         atemp(1,1) = polain(1,ialoc,ipt)
         atemp(1,2) = polain(4,ialoc,ipt)
         atemp(1,3) = polain(5,ialoc,ipt)
         atemp(2,1) = polain(4,ialoc,ipt)
         atemp(2,2) = polain(2,ialoc,ipt)
         atemp(2,3) = polain(6,ialoc,ipt)
         atemp(3,1) = polain(5,ialoc,ipt)
         atemp(3,2) = polain(6,ialoc,ipt)
         atemp(3,3) = polain(3,ialoc,ipt)

! ... calculate btemp = eivr*atemp

         do m = 1, 3
            btemp(1,m) = eivr(1,1)*atemp(1,m)+eivr(2,1)*atemp(2,m)+eivr(3,1)*atemp(3,m)
            btemp(2,m) = eivr(1,2)*atemp(1,m)+eivr(2,2)*atemp(2,m)+eivr(3,2)*atemp(3,m)
            btemp(3,m) = eivr(1,3)*atemp(1,m)+eivr(2,3)*atemp(2,m)+eivr(3,3)*atemp(3,m)
         end do

! ... calculate atemp = btemp*eivr(-1)

         do m = 1, 3
            atemp(1,m) = btemp(1,1)*eivr(1,m)+btemp(1,2)*eivr(2,m)+btemp(1,3)*eivr(3,m)
            atemp(2,m) = btemp(2,1)*eivr(1,m)+btemp(2,2)*eivr(2,m)+btemp(2,3)*eivr(3,m)
            atemp(3,m) = btemp(3,1)*eivr(1,m)+btemp(3,2)*eivr(2,m)+btemp(3,3)*eivr(3,m)
         end do

! ... restore

         poltensa(1,ialoc,ipt) = atemp(1,1)
         poltensa(2,ialoc,ipt) = atemp(2,2)
         poltensa(3,ialoc,ipt) = atemp(3,3)
         poltensa(4,ialoc,ipt) = atemp(1,2)
         poltensa(5,ialoc,ipt) = atemp(1,3)
         poltensa(6,ialoc,ipt) = atemp(2,3)
      end do

   end if

end subroutine PAxesSystem

!************************************************************************
!> \page particle particle.F90
!  **SetAtomProp**
! *set atom properties*
!************************************************************************


subroutine SetAtomProp(iplow, ipupp, lint)

   use MolModule
   implicit none

   integer(4), intent(in) :: iplow
   integer(4), intent(in) :: ipupp
   logical,    intent(in) :: lint

   call SetAtomPos(iplow, ipupp, lint)
   if (ldipole .or. lpolarization .or. ldipolesph) call SetAtomDipMom(iplow, ipupp)
   if (lpolarization) call SetAtomPolTens(iplow, ipupp)
   if (laimage) call SetImageSph(iplow,ipupp,1)

end subroutine SetAtomProp

