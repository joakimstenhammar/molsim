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


!> \page nmlStatic
!! The namelist  \ref nmlStatic contains variables that control the interval of the analysis and static analysis
!! routines used.
!! * Variables:
!!  * \subpage istatic
!!  * \subpage lspdf
!!  * \subpage lrdf
!!  * \subpage lrdfchain
!!  * \subpage lrdfsph
!!  * \subpage lg3
!!  * \subpage lrdfcond
!!  * \subpage lsf
!!  * \subpage langdf
!!  * \subpage langextdf
!!  * \subpage loridipdf
!!  * \subpage llsphharaver
!!  * \subpage lradangdf
!!  * \subpage lkirkwoodgk
!!  * \subpage loripoldf
!!  * \subpage lnnhb
!!  * \subpage lnndf
!!  * \subpage lchaindf
!!  * \subpage lchaintypedf
!!  * \subpage lchaintypeextdf
!!  * \subpage lcbpc
!!  * \subpage lltt
!!  * \subpage lcluster
!!  * \subpage lzerosecondmoment
!!  * \subpage lmultipoledf
!!  * \subpage lenergydf
!!  * \subpage lwidom1
!!  * \subpage lwidom2
!!  * \subpage lmeanforce1
!!  * \subpage lmeanforce2
!!  * \subpage lpotmeanforce
!!  * \subpage lsurfacearea
!!  * \subpage lcrystalformat
!!  * \subpage ltrajectory
!!  * \subpage lsubstructuredf
!!  * \subpage lnetworkdf
!!  * \subpage lnetworkradialdf
!!  * \subpage lstaticuser

!> \page lspdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Single particle distribution functions are calculated. Further specification is given in namelist  \ref nmlSPDF.
!! * `.false`: No calculation.

!> \page lrdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Radial distribution functions or running coordination numbers or are calculated.  Further specification is given in
!! namelist  \ref nmlRDF.
!! * `.false`: No calculation.

!> \page lrdfchain
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Radial distribution functions or running coordination numbers are calculated between center of masses of chains. Further specification is given in namelist \ref nmlRDFChain.
!! * `.false`: No calculation.

!> \page lrdfsph
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Radial distribution functions or running coordination numbers are calculated for particles which positions are projected on a sphere. Further specification is given in namelist  \ref nmlRDFSph.
!! * `.false`: No calculation.

!> \page lg3
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Normalized triplet correlation functions are calculated. Further specification is given in namelist  \ref nmlG3Dist.
!! * `.false`: No calculation.

!> \page lrdfcond
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Conditional radial distribution functions are calculated. Futher specification is given in namelist  \ref nmlRDFCond.
!! * `.false.`: No calculation.

!> \page lsf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Partial structure factors are calculated. Further specification is given in namelist \ref nmlSF.
!! * `.false.`: No calculation.

!> \page langdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Angular distribution functions are calculated. Further specification is given in namelist  \ref nmlAngDF.
!! * `.false.`: No calculation.

!> \page langextdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Angular 2d distribution functions with respect to external frame are calculated. Further specification is given in namelist  \ref nmlAngExtDF.
!! * `.false.`: No calculation.

!> \page loridipdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Orientation/dipole distribution functions are calculated. Further specification is given in namelist  \ref nmlOriDipDF.
!! * `.false.`: No calculation.

!> \page llsphharaver
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Averages of unnormalized spherical harmonics are calculated.
!! * `.false.`: No calculation.

!> \page lradangdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Radial-angular 2d distribution functions are calculated. Further specification is given in namelist  \ref nmlRadAngDF.
!! * `.false.`: No calculation.

!> \page lkirkwoodgk
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Kirkwood gk-factors are calculated. Further specification is given in namelist \ref nmlKirkwoodgk.
!! * `.false.`: No calculation.

!> \page loripoldf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Orientation polarization distribution functions are calculated. Further specification is given in namelist  \ref nmlOriPolDF.
!! * `.false.`: No calculation.

!> \page lnnhb
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Number of nearest neighbours and no of hydrogen bonds are calculated. Further specification is given in namelist  \ref nmlNNHB.
!! * `.false.`: No calculation.

!> \page lnndf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Nearest neighbour distribution functions are calculated. Further specification is given in namelist  \ref nmlNNDF.
!! * `.false.`: No calculation.

!> \page lchaindf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Chain distribution functions are calculated. Further specification is given in namelist  \ref nmlChainDF.
!! * `.false.`: No calculation.

!> \page lchaintypedf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Chain distribution functions are calculated. Further specification is given in namelist  \ref nmlChainTypeDF.
!! * `.false.`: No calculation.

!> \page lchaintypeextdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Chain type distribution functions with respect to the lab frame are calculated. Further specification is given in namelist  \ref nmlChainTypeExtDF.
!! * `.false.`: No calculation.

!> \page lcbpc
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Probabilities of chain particles to be near particles of another type are calculated.  Further specification is given in namelist  \ref nmlCBPC.
!! * `.false.`: No calculation.

!> \page lltt
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Loop, tail, and train statistics are calculated. Further specification is given in namelist  \ref nmlLoopTailTrain.
!! * `.false.`: No calculation.

!> \page lcluster
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Cluster size distribution functions are calculated. Further specification is given in namelist  \ref nmlCluster.
!! * `.false.`: No calculation.

!> \page lzerosecondmoment
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Zero and second moment for an ionic, neutral system is calculated.
!! * `.false.`: No calculation.

!> \page lmultipoledf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Electrostatic multipole moment distribution functions are calculated. Further specification is given in namelists  \ref nmlMultipoleDF.
!! * `.false.`: No preparation.

!> \page lenergydf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Energy distribution functions are calculated. Further specification is given in namelist  \ref nmlEnergyDF.
!! * `.false.`: No calculation.

!> \page lwidom1
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Excess chemical potentials are prepared. Further specification is given in namelists \ref nmlWidom1.
!! * `.false.`: No preparation.

!> \page lwidom2
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Excess chemical potentials are prepared. Further specification is given in namelists \ref nmlWidom2.
!! * `.false.`: No preparation.

!> \page lmeanforce1
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Mean force between two particles is calculated. Further specification is given in namelist  \ref nmlMeanForce1.
!! * `.false.`: No calculation.

!> \page lmeanforce2
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Mean force between two particles is calculated. Further specification is given in namelist  \ref nmlMeanForce2.
!! * `.false.`: No calculation.

!> \page lpotmeanforce
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Potential of mean force between two particles is calculated. Further specification is given in namelist  \ref nmlPotMeanForce.
!! * `.false.`: No calculation.

!> \page lsurfacearea
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Surface area available around atoms residing in particles of a specified type is calculated. Further specification is given in namelist  \ref nmlSurfaceArea.
!! * `.false.`: No calculation.

!> \page lcrystalformat
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Write atoms in the crystallographic format stating with atoms belong to molecules closest to the origin. Hard-code limitation of spatial distance and number of particles are in action.
!! * `.false.`: No calculation.

!> \page ltrajectory
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Write trajectory on FLIST. Further specification is given in namelist \ref nmlTrajectory.
!! * `.false.`: No calculation.

!> \page lsubstructuredf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Substructures of chains, hierarchical structures and networks are analyses. Further specification is given in namelist  nmlSustructureDF.
!! * `.false.`: No calculation.

!> \page lnetworkdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Network distribution functions are calculated. Further specification is given in \ref nmlNetworkDF.
!! * `.false.`: No calculation.

!> \page lnetworkradialdf
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Radial network distribution functions are calculated. Further specification is given in \ref nmlNetworkRadialDF.
!! * `.false.`: No calculation.

!> \page lstaticuser
!! `logical`
!! **default:** `.false.`
!! * `.true.`: StaticUser is called and from where user-provided static analysis routines are called in file moluser.F90.
!! * `.false.`: No call.

!************************************************************************
!> \page static static.F90
!! **StaticDriver**
!! *driver of static analysis routines*
!************************************************************************

subroutine StaticDriver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='StaticDriver'
   character(80), parameter :: txheading ='static analysis: general'

   logical,       save :: lspdf, lrdf, lrdfchain, lrdfsph, lg3, lrdfcond, lsf,                  &
                          langdf, langextdf, loridipdf, lsphharaver, lradangdf,                 &
                          lkirkwoodgk, loripoldf, lnnhb, lnndf,                                 &
                          lchaindf, lchaintypedf, lchaintypeextdf, lcbpc, lltt,                 &
                          lcluster, lzerosecondmoment, lmultipoledf,                            &
                          lenergydf, lwidom1, lwidom2, lmeanforce1, lmeanforce2, lpotmeanforce, &
                          lsurfacearea, lcrystalformat, ltrajectory, lsubstructuredf,           &
                          lnetworkdf, lnetworkradialdf,                                         &
                          lstaticuser
   type(scalar_var), allocatable, save :: var(:)
   integer(4), save :: nvar = 1

   namelist /nmlStatic/ istatic,                                                               &
                        lspdf, lrdf, lrdfchain, lrdfsph, lg3, lrdfcond, lsf,                   &
                        langdf, langextdf, loridipdf, lsphharaver, lradangdf,                  &
                        lkirkwoodgk, loripoldf, lnnhb, lnndf,                                  &
                        lchaindf, lchaintypedf, lchaintypeextdf, lcbpc, lltt,                  &
                        lcluster, lzerosecondmoment, lmultipoledf,                             &
                        lenergydf, lwidom1, lwidom2, lmeanforce1, lmeanforce2, lpotmeanforce,  &
                        lsurfacearea, lcrystalformat, ltrajectory, lsubstructuredf,            &
                        lnetworkdf, lnetworkradialdf,                                          &
                        lstaticuser

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

      istatic          = 1
      lspdf            = .false.
      lrdf             = .false.
      lrdfchain        = .false.
      lrdfsph          = .false.
      lrdfcond         = .false.
      lg3              = .false.
      lsf              = .false.
      langdf           = .false.
      langextdf        = .false.
      loridipdf        = .false.
      lsphharaver      = .false.
      lradangdf        = .false.
      lkirkwoodgk      = .false.
      loripoldf        = .false.
      lnnhb            = .false.
      lnndf            = .false.
      lchaindf         = .false.
      lchaintypedf     = .false.
      lchaintypeextdf  = .false.
      lcbpc            = .false.
      lltt             = .false.
      lcluster         = .false.
      lzerosecondmoment= .false.
      lmultipoledf     = .false.
      lenergydf        = .false.
      lwidom1          = .false.
      lwidom2          = .false.
      lmeanforce1      = .false.
      lmeanforce2      = .false.
      lpotmeanforce    = .false.
      lsurfacearea     = .false.
      lcrystalformat   = .false.
      ltrajectory      = .false.
      lsubstructuredf  = .false.
      lnetworkdf       = .false.
      lnetworkradialdf = .false.
      lstaticuser      = .false.

      rewind(uin)
      read(uin,nmlStatic)

      call StaticDriverSub

   case (iWriteInput)

      nvar = 1
      allocate(var(nvar))
      var(1)%label = 'volume'

      if(istatic > nstep2) then
         call Stop(txroutine, 'istatic > nstep2', uout)
      end if


      if (.not.lgroup) then
         if (lspdf           ) call Stop(txroutine, 'spdf is selected, but no group division', uout)
         if (lrdf            ) call Stop(txroutine, 'rdf is selected, but no group division', uout)
         if (langdf          ) call Stop(txroutine, 'angdf is selected, but no group division', uout)
         if (lnnhb           ) call Stop(txroutine, 'nnhb is selected, but no group division', uout)
         if (lnndf           ) call Stop(txroutine, 'nndf is selected, but no group division', uout)
         if (lenergydf       ) call Stop(txroutine, 'energydf is selected, but no group division', uout)
         if (lsubstructuredf ) call Stop(txroutine, 'substructuredf is selected, but no group division', uout)
         if (lnetworkdf      ) call Stop(txroutine, 'networkdf is selected, but no group division', uout)
         if (lnetworkradialdf) call Stop(txroutine, 'networkradialdf is selected, but no group division', uout)
      end if

      call StaticDriverSub

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var
      call StaticDriverSub

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      call StaticDriverSub

   case (iSimulationStep)

      if (mod(istep2,istatic) == 0) then
         var(1)%value = vol
         call ScalarSample(iStage, 1, nvar, var)
         call StaticDriverSub
      end if

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call StaticDriverSub

   case (iAfterSimulation)

      if (master) then
         call ScalarSample(iStage, 1, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,i10)')     'sampling interval              = ', istatic
         write(uout,'(a,t35,i10)')     'no of time steps/passes used   = ', var(1)%nsamp2*var(1)%nsamp1
         write(uout,'(a,t35,2es12.3)') 'average volume                 = ', var(1)%avs1, var(1)%avsd
         write(uout,'()')
         write(uout,'(a)') 'static analysis routines used'
         write(uout,'(a)') '-----------------------------'
         if (lspdf)             write(uout,'(a)') '   spdf          '
         if (lrdf)              write(uout,'(a)') '   rdf           '
         if (lrdfchain)         write(uout,'(a)') '   rdfchain      '
         if (lrdfsph)           write(uout,'(a)') '   rdfsph        '
         if (lrdfcond)          write(uout,'(a)') '   rdfcond       '
         if (lg3)               write(uout,'(a)') '   g3            '
         if (lsf)               write(uout,'(a)') '   sf            '
         if (langdf)            write(uout,'(a)') '   angdf         '
         if (langextdf)         write(uout,'(a)') '   angextdf      '
         if (loridipdf)         write(uout,'(a)') '   oridipdf      '
         if (lsphharaver)       write(uout,'(a)') '   sphharaver    '
         if (lradangdf)         write(uout,'(a)') '   radangdf      '
         if (lkirkwoodgk)       write(uout,'(a)') '   kirkwoodgk    '
         if (loripoldf)         write(uout,'(a)') '   oripoldf      '
         if (lnnhb)             write(uout,'(a)') '   nnhb          '
         if (lnndf)             write(uout,'(a)') '   nndf          '
         if (lchaindf)          write(uout,'(a)') '   chaindf       '
         if (lchaintypedf)      write(uout,'(a)') '   chaintypedf   '
         if (lchaintypeextdf)   write(uout,'(a)') '   chaintypedflab'
         if (lcbpc)             write(uout,'(a)') '   chainbead_part_contact'
         if (lltt)              write(uout,'(a)') '   ltt           '
         if (lcluster)          write(uout,'(a)') '   cluster       '
         if (lzerosecondmoment) write(uout,'(a)') '   zerosecondmoment'
         if (lmultipoledf)      write(uout,'(a)') '   multipoledf   '
         if (lenergydf)         write(uout,'(a)') '   energydf      '
         if (lwidom1)           write(uout,'(a)') '   widom1        '
         if (lwidom2)           write(uout,'(a)') '   widom2        '
         if (lmeanforce1)       write(uout,'(a)') '   meanforce1    '
         if (lmeanforce2)       write(uout,'(a)') '   meanforce2    '
         if (lpotmeanforce)     write(uout,'(a)') '   potmeanforce  '
         if (lsurfacearea)      write(uout,'(a)') '   surfacearea   '
         if (lcrystalformat)    write(uout,'(a)') '   crystalformat '
         if (ltrajectory)       write(uout,'(a)') '   trajectory    '
         if (lstaticuser)       write(uout,'(a)') '   staticuser    '
         if (lnetworkdf)        write(uout,'(a)') '   networkdf     '
         if (lnetworkradialdf)  write(uout,'(a)') '   networkradialdf'
         if (lsubstructuredf)   write(uout,'(a)') '   substructuredf'
      end if

      if(master) call fileflush(uout)
      call StaticDriverSub

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine StaticDriverSub
   if (lspdf)             call SPDF(iStage)
   if (lrdf)              call RDF(iStage)
   if (lrdfchain)         call RDFChain(iStage)
   if (lrdfsph)           call RDFSph(iStage)
   if (lrdfcond)          call RDFCond(iStage)
   if (lg3)               call G3Dist(iStage)
   if (lsf)               call SFDriver(iStage)
   if (langdf)            call AngDF(iStage)
   if (langextdf)         call AngExtDF(iStage)
   if (loridipdf)         call OriDipDF(iStage)
   if (lsphharaver)       call SphHarAver(iStage)
   if (lradangdf)         call RadAngDF(iStage)
   if (lkirkwoodgk)       call Kirkwoodgk(iStage)
   if (loripoldf)         call OriPolDF(iStage)
   if (lnnhb)             call NNHB(iStage)
   if (lnndf)             call NNDF(iStage)
   if (lchaindf)          call ChainDF(iStage)
   if (lchaintypedf)      call ChainTypeDF(iStage)
   if (lchaintypeextdf)   call ChainTypeExtDF(iStage)
   if (lcbpc)             call ChainBeadPartContact(iStage)
   if (lltt)              call LoopTailTrain(iStage)
   if (lcluster)          call ClusterSD(iStage)
   if (lzerosecondmoment) call ZeroSecondMoment(iStage)
   if (lmultipoledf)      call MultipoleDF(iStage)
   if (lenergydf)         call EnergyDF(iStage)
   if (lwidom1)           call Widom1(iStage)
   if (lwidom2)           call Widom2(iStage)
   if (lmeanforce1)       call MeanForce1(iStage)
   if (lmeanforce2)       call MeanForce2(iStage)
   if (lpotmeanforce)     call PotMeanForce(iStage)
   if (lsurfacearea)      call SurfaceArea(iStage)
   if (lcrystalformat)    call Crystalformat(iStage)
   if (ltrajectory)       call Trajectory(iStage)
   if (lsubstructuredf)   call SubStructureDF(iStage)
   if (lnetworkdf)        call NetworkDF(iStage)
   if (lnetworkradialdf)  call NetworkRadialDF(iStage)
   if (lstaticuser)       call StaticUser(iStage)
end subroutine StaticDriverSub

!........................................................................

end subroutine StaticDriver

!************************************************************************
!> \page static static.F90
!! **SPDF**
!! *calculate single particle distribution functions*
!************************************************************************



!     group division of particles is required

!> \page nmlSPDF
!! The namelist  \ref nmlSPDF contains variables that control the calculation of single particle distribution functions. Any combination
!! of the types of distribution functions listed below may be selected through vtype\%l.
!!    | type | label           | Quantity                                                                    |
!!    | :--: | :-------------: | :-------------------------------------------------------------------------: |
!!    | 1    | rrden           | Reduced radial number density                                               |
!!    | 2    | rren            | Reduced running coordination number                                         |
!!    | 3    | zden            | Number density in the z-direction                                           |
!!    | 4    | zden2           | Number density in the z-direction (special)                                 |
!!    | 5    | z'*z            | Projection of molecular z'-axis on the box z-axis                           |
!!    | 6    | z'*x            | Projection of molecular z'-axis on the box x-axis                           |
!!    | 7    | z'*y            | Projection of molecular z'-axis on the box y-axis                           |
!!    | 8    | \format{opt_d}  | Radial density projected on the z = 0 plane                                 |
!!    | 9    | elpot           | radial electrostaic potential evaluated at particle positions               |
!!    | 10   | m*E             | cos(m(ia)*E); m(ia) is the dipole moment of atom ia and E an external field |
!! * Variables:
!!  * \subpage nmlSPDF_vtype

!> \page nmlSPDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real) (1:10)`
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * min: /0.0,0.0,-X,-X,-1.0,Y/
!! * X = lcyl/2 (only \ref txbc = 'cyl'), box(3)/2 (else)
!! * Y = rsph (only \ref txbc = 'sph), rcyl (only \ref txbc = 'cyl')
!! * max /10.0,10.0,X,X,1.0,0.0/


subroutine SPDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SPDF'
   character(80), parameter :: txheading ='single particle distribution functions'
   integer(4)   , parameter :: ntype = 10
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4) :: itype, ivar, ibin, ip, igr, jp
   real(8)    :: r1, r2, ac, norm, vsum, dvol, darea, ui, uuu, fdum(3), pot, InvFlt, angcos

   namelist /nmlSPDF/ vtype

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      if (lbcrd .or. lbcto) call Stop(txroutine,'Erroneous boundary conditions',uout)

      vtype%l   =.false.
      vtype(1:2)%min = Zero
      vtype(1:2)%max = 10.0d0
      if (lbccyl) then
         vtype(3:4)%min =-Half*cyllen
         vtype(3:4)%max =+Half*cyllen
      else
         vtype(3:4)%min =-boxlen2(3)
         vtype(3:4)%max =+boxlen2(3)
      end if
      vtype(5:7)%min =-One
      vtype(5:7)%max =+One
      if (lbcsph) then
         vtype(8)%max = sphrad
      else if (lbccyl) then
         vtype(8)%max = cylrad
      end if
      vtype(9)%min = Zero
      vtype(9)%max = 10.0d0
      vtype%nbin   = 100
      vtype(10)%min =-One
      vtype(10)%max =+One

      rewind(uin)
      read(uin,nmlSPDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

      if (txuser ==  'jos' .and. lbccyl) vtype(8)%max = cylrad  ! pick up change

! ... set remaining elements of vtype

      vtype%label = ['rrden','rrdf ','zden ','zden2','z''*z ','     ','     ','opt_d','elpot','m*E  ']
      vtype%nvar = ngr(1)

! ... check some conditions

      if (vtype(10)%l) then
         if (count(abs(efield_ext(1:3)) /= Zero) == 0) call stop(txroutine, 'vtype(10)%l and efield_ext = 0' , uout)
      end if

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(ngr(1),ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do igr = 1, ngr(1)
               ivar = ivar+1
               ipnt(igr,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = 1, np
         igr = igrpn(ip,1)
         if (igr <= 0) cycle
         r2 = ro(1,ip)**2+ro(2,ip)**2+ro(3,ip)**2

! ... sample type 1

         itype = 1
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            r1 = sqrt(r2)
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 2

         itype = 2
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            r1 = sqrt(r2)
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 3

         itype = 3
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
!           ac = ro(1,ip)    ! x-direction
            ac = ro(3,ip)
            ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 4

         itype = 4
         if (vtype(itype)%l) then
            if (icnpn(ip) /= 0) then                          ! temporary
            if (isegpn(ip) == 20 .or. isegpn(ip) == 40) then  ! temporary
            ivar = ipnt(igr,itype)
            ac = ro(3,ip)
            ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if                                           ! temporary
            end if                                           ! temporary
         end if

! ... sample type 5

         itype = 5
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            ac = ori(3,3,ip)
            ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 6

         itype = 6
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            ac = ori(1,3,ip)
            ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 7

         itype = 7
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            ac = ori(2,3,ip)
            ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 8

         itype = 8
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            r1 = sqrt(ro(1,ip)**2 + ro(2,ip)**2)
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end if

! ... sample type 9

         itype = 9
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            ui = Zero
            do jp = 1, np
               if (ip == jp) cycle
               call UTwoBodyPair(ip, jp, uuu, fdum)
               ui = ui+uuu
            end do
            pot = Zero
            if (abs(az(ip)) > 1d-10) pot=ui*sclene/(az(ip)*Ech*AvNo)
            r1 = sqrt(r2)
            if (luext) then                                                ! Steffi: not checked from here 2013-07-01
               if (txuext(iptpn(ip)) == 'insulating_sphere') then
                  if (r1 >= rInsSphere) then
                     pot = pot+zInsSphere*EpsiFourPi/r1/ech
                  else
                     pot = pot+zInssphere*EpsiFourPi/ech/(Two*rInsSphere)*(Three-r2/rInsSphere**2)
                  end if
               end if
            end if                                                         ! to here  /Per
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + pot
            var(ivar)%nsampbin(ibin)=var(ivar)%nsampbin(ibin) + One
         end if

! ... sample type 10

            itype = 10
            if (vtype(itype)%l) then
               ivar = ipnt(igr,itype)
               ac = angcos(dip(1,ip),dip(2,ip),dip(3,ip),efield_ext(1),efield_ext(2),efield_ext(3))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
      end do

   case (iAfterMacrostep)

      do igr = 1, ngr(1)
         itype = 1
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            vsum = sum(var(ivar)%avs2(0:var(ivar)%nbin-1))
            norm = var(ivar)%nsamp2*(var(ivar)%max**3-var(ivar)%min**3)*InvFlt(vsum)
            do ibin = 0, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
         itype = 2
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            norm = InvFlt(grvar(igrpnt(1,igr))%avs2)
            var(ivar)%avs2(0) = var(ivar)%avs2(0)*norm
            do ibin = 1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1) + var(ivar)%avs2(ibin)*norm
            end do
         end if
         itype = 3
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            if (lbccyl) then
               norm = One/(Pi*cylrad**2*var(ivar)%bin)
            else
               norm = One/(boxlen(1)*boxlen(2)*var(ivar)%bin)
            end if
            do ibin = 0, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm
            end do
         end if
         itype = 4
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            norm = One/(boxlen(1)*boxlen(2)*var(ivar)%bin)
            do ibin = 0, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm
            end do
         end if
         itype = 5
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            call DistFuncNorm(ivar, ivar, var)
         end if
         itype = 8
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            do ibin = 0, var(ivar)%nbin
               if (txuser == 'jos') then ! Jos' project
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)/(pi*darea(ibin,var(ivar)%min,var(ivar)%bin)*cyllen)
               else ! orig
                   var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)/(pi*darea(ibin,var(ivar)%min,var(ivar)%bin))
               end if
            end do
         end if
         itype = 9
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            norm = var(ivar)%nsamp2              ! factor too counteract the normalization in DistFuncSample
            do ibin = -1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = norm*var(ivar)%avs2(ibin)*InvFlt(var(ivar)%nsampbin(ibin))
            end do
         end if
         itype = 10
         if (vtype(itype)%l) then
            ivar = ipnt(igr,itype)
            call DistFuncNorm(ivar, ivar, var)
         end if
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

! ... calculate excess ammount

      itype = 3
      if (vtype(itype)%l) call ExcessAmount(vtype(itype)%nbin, ngr(1), txgr, ipnt(1:ngr(1),itype), var, boxlen, uout)

      if (txuser == 'jos' .and. vtype(8)%l) call SPDF_jos

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine SPDF_jos
   real(8) :: rho1_0, rho2_0, rho2_av, n_mon
   itype = 8
   ivar = ipnt(1,itype)
   rho1_0 = var(ivar)%avs1(var(ivar)%nbin-1)
   ivar = ipnt(2,itype)
   rho2_0 = var(ivar)%avs1(var(ivar)%nbin-1)
   rho2_av = nppt(2)/(cyllen*pi*(cylrad**2-rCylinder**2))
   n_mon = float(nppt(1)-nppt(2))/(10*abs(zCylinder))  ! 10 denotes the length of one momomer
   call FileOpen(uuser, fuser, 'form/noread')
   write(*,*) 'n_mon, nppt(2)',n_mon, nppt(2)
   write(uuser,'(2f10.4)') n_mon/nppt(2), rho1_0*rho2_0/rho2_av**2 - one ! output
   close(uuser)
end subroutine SPDF_jos

!........................................................................

end subroutine SPDF

!************************************************************************
!> \page static static.F90
!! **RDF**
!! *calculate running coordination number or radial distribution function*
!************************************************************************


!> \page nmlRDF
!! The namelist  \ref nmlRDF contains variables that control the calculation of running coordination number or radial distribution
!! functions. Any combination of the types of distribution functions listed below may be selected through vtype\%l.
!!    | type | label   |   quantity        |
!!    | :--: | :-----: | :---------------: |
!!    | 1    | com-com | particle-particle |
!!    | 2    | com-xxx | particle-atom     |
!!    | 3    | xxx-xxx | atom-atom         |
!! * Variables:
!!  * \subpage nmlRDF_vtype
!!  * \subpage nmlRDF_rmax
!!  * \subpage nmlRDF_ndim
!!  * \subpage nmlRDF_nbin
!!  * \subpage nmlRDF_func
!!  * \subpage l2dtwo

!     xxx is the first three characters of txat
!     group division of particles is required

!> \page nmlRDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)` (1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /0.0,0.0,0.0/ Max: /10.0,10.0,10.0/

!> \page nmlRDF_rmax rmax
!! `real`
!! **default:** `10.0`
!! * Upper distance for particle separation considered.

!> \page nmlRDF_ndim ndim
!! `integer`
!! **default:** `3`
!! * `2`: Sample distribution function in the xy-plane.
!! * `3`: Sample distribution function in the xyz-space.

!> \page nmlRDF_nbin nbin
!! `integer`
!! **default:** `100`
!! * Number of bins used to sample the distribution functions.

!> \page nmlRDF_func func
!! `character(3)`
!! **default:** `rdf`
!! * `rdf`: Radial distribution functions are calculated.
!! * `rcn`: Running coordination numbers are calculated.

!> \page l2dtwo
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Sampling of 2d radial distribution functions separately for z>0 and z<0. Useful for a system with \ref txbc='xy' and two equivalent surfaces at z=±box2(3) (only ndim = 2)
!! * `.false.`: Nothing.

subroutine RDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='RDF'
   character(80), parameter :: txheading ='running coordinate number/radial distribution function'
   integer(4)   , parameter :: ntype = 3
   type(static1D_var),         save :: vtype(ntype)
   logical,                    save :: l2dtwo
   real(8),                    save :: rmax
   integer(4),                 save :: ndim
   character(3),               save :: func
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   type(df_var),  allocatable, save :: tauvar(:)    !jvr storage for <rab.rij/|rab|>shell
   type(df_var),  allocatable, save :: tauvarvar(:) !jvr storage for var*<rab.rij/|rab|>shell
   integer(4),    allocatable, save :: ipnt(:,:,:)
   integer(4), parameter :: nvar_s = 2
   type(scalar_var), save :: var_s(nvar_s)

   real(8),       save :: rdfvols2, rdfvolfac
   integer(4) :: itype, ivar, ibin, ip, ipt, jp, jpt, ia, iat, ja, jat, iatjat, igr, jgr, igrjgr
   real(8)    :: dx, dy, dz, dropbc(3), r1, r2, dvol, darea, InvFlt, dxo, dyo, dzo

   namelist /nmlRDF/ vtype, rmax, ndim, func, l2dtwo

   if (ltrace) call WriteTrace(2, txheading, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = Zero
      vtype%max = 10.0d0
      vtype%nbin = 100
      rmax   = 10.0d0
      ndim   = 3
      func   ='rdf'
      l2dtwo =.false.

      rewind(uin)
      read(uin,nmlRDF)

      call LowerCase(func)

      if (txbc == 'xyz') rmax = min(rmax, minval(boxlen2(1:3)))    ! Maybe switch off Jos
      if (txbc == 'xy') rmax = min(rmax, minval(boxlen2(1:2)))
      if (txbc == 'z') rmax = min(rmax, boxlen2(3))
      if (lbcrd) rmax = sqrt(Two/Three)*cellside
      if (lbcto) rmax = sqrt(Three/Two)*cellside
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, ' vtype%nbin > mnbin_df', uout)
      if (ndim == 3) l2dtwo =.false.
      if ((func /= 'rdf') .and. (func /= 'rcn')) call Stop(txheading, 'check value of func', uout)

      if (ndim == 3) rdfvolfac = FourPiThird
      if (ndim == 2) rdfvolfac = Pi

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%nvar = [ngr(1)*ngr(2), ngr(1)*sum(natgr(1:ngr(2),2)), sum(natgr(1:ngr(1),1))*sum(natgr(1:ngr(2),2))]

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nat**2,ngrgr,ntype))
      ipnt = 0

!jvr ... if sample tau for pressure calculations, allocate memory
!     possible further improvement: move the complete hard sphere pressure calculation to this routine, as for single particles type 1 distribution functions can be used

      allocate(tauvar(nvar), tauvarvar(nvar))

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      itype = 1
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            do jgr = 1, ngr(2)
               igrjgr = igrgr(igr,jgr)
               ivar = ivar+1
               ipnt(1,igrjgr,itype) = ivar
               var(ivar)%label = trim(txgrgr(igrjgr))//' '//'com'//'-'//'com'
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
!jvr dummy
               tauvar(ivar)%label = trim(txgrgr(igrjgr))//' '//'com'//'-'//'com'
               tauvar(ivar)%min = vtype(itype)%min
               tauvar(ivar)%max = vtype(itype)%max
               tauvar(ivar)%nbin = vtype(itype)%nbin
            end do
         end do
      end if
      itype = 2
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            do jgr = 1, ngr(2)
               igrjgr = igrgr(igr,jgr)
               do jat = iatgr(jgr,2), iatgr(jgr,2)+natgr(jgr,2)-1
                  ivar = ivar+1
                  ipnt(jat,igrjgr,itype) = ivar
                  var(ivar)%label = trim(txgrgr(igrjgr))//' '//'com'//'-'//txat(jat)(1:3)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
!jvr dummy
                  tauvar(ivar)%label = trim(txgrgr(igrjgr))//' '//'com'//'-'//txat(jat)(1:3)
                  tauvar(ivar)%min = vtype(itype)%min
                  tauvar(ivar)%max = vtype(itype)%max
                  tauvar(ivar)%nbin = vtype(itype)%nbin
               end do
            end do
         end do
      end if
      itype = 3
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            do jgr = 1, ngr(2)
               igrjgr = igrgr(igr,jgr)
               do iat = iatgr(igr,1), iatgr(igr,1)+natgr(igr,1)-1
                  do jat = iatgr(jgr,2), iatgr(jgr,2)+natgr(jgr,2)-1
                     iatjat = iat+nat*(jat-1)
                     ivar = ivar+1
                     ipnt(iatjat,igrjgr,itype) = ivar
                     var(ivar)%label = trim(txgrgr(igrjgr))//' '//txat(iat)(1:3)//'-'//txat(jat)(1:3)
                     var(ivar)%min = vtype(itype)%min
                     var(ivar)%max = vtype(itype)%max
                     var(ivar)%nbin = vtype(itype)%nbin
!jvr setup tau
                     tauvar(ivar)%label = trim(txgrgr(igrjgr))//' '//txat(iat)(1:3)//'-'//txat(jat)(1:3)
                     tauvar(ivar)%min = vtype(itype)%min
                     tauvar(ivar)%max = vtype(itype)%max
                     tauvar(ivar)%nbin = vtype(itype)%nbin
                  end do
               end do
            end do
         end do
      end if

!jvr setup scalar samples
      var_s(1)%label = 'volume'
      var_s(2)%label = 'pressure (hs)'
      call DistFuncSample(iStage, nvar, var)
      call DistFuncSample(iStage, nvar, tauvar) !jvr

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      call DistFuncSample(iStage, nvar, tauvar)
      call ScalarSample(iStage, 1, nvar_s, var_s)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var, tauvar, var_s

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      call DistFuncSample(iStage, nvar, tauvar)    !jvr
      call ScalarSample(iStage, 1, nvar_s, var_s)
      rdfvols2 = Zero

   case (iSimulationStep)

      if (ndim == 3) rdfvols2 = rdfvols2+vol
      if (ndim == 2) rdfvols2 = rdfvols2+boxlen(1)*boxlen(2)

      var_s(1)%value = vol
      call ScalarSample(iStage,1,1,var_s)              ! volume

      var%nsamp2 = var%nsamp2 + 1
      tauvar%nsamp2 = 1

       do ip = ipmyid(1), ipmyid(2)
         igr = igrpn(ip,1)
         if (igr <= 0) cycle
         ipt = iptpn(ip)
         do jp = 1, np
            if (jp == ip) cycle
            jgr = igrpn(jp,2)
            if (jgr <= 0) cycle
            jpt = iptpn(jp)
            igrjgr = igrgr(igr,jgr)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC2(dx,dy,dz,dropbc(1),dropbc(2),dropbc(3))
            dx = dx-dropbc(1)
            dy = dy-dropbc(2)
            dz = dz-dropbc(3)
            if (ndim == 2) dz = Zero
            if (l2dtwo .and. (ro(3,ip)*ro(3,jp) < Zero)) cycle  ! two separate slabs
            r2 = dx**2+dy**2+dz**2
            if (r2 < rmax**2) then

! ... sample type 1

               itype = 1
               if (vtype(itype)%l) then
                  r1 = sqrt(r2)
                  ivar = ipnt(1,igrjgr,itype)
                  ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if

! ... sample type 2

               itype = 2
               if (vtype(itype)%l) then
                  do ja = ianpn(jp), ianpn(jp)+napt(jpt)-1
                     jat = iatan(ja)
                     dx = ro(1,ip)-r(1,ja)-dropbc(1)
                     dy = ro(2,ip)-r(2,ja)-dropbc(2)
                     dz = ro(3,ip)-r(3,ja)-dropbc(3)
                     if (ndim == 2) dz = Zero
                     r1 = sqrt(dx**2+dy**2+dz**2)
                     ivar = ipnt(jat,igrjgr,itype)
                     ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
                  end do
               end if

! ... sample type 3

               itype = 3
               if (vtype(itype)%l) then
                  dxo = dx
                  dyo = dy
                  dzo = dz
                  do ia = ianpn(ip), ianpn(ip)+napt(ipt)-1
                     iat = iatan(ia)
                     do ja = ianpn(jp), ianpn(jp)+napt(jpt)-1
                        jat = iatan(ja)
                        iatjat = iat+nat*(jat-1)
                        dx = r(1,ia)-r(1,ja)-dropbc(1)
                        dy = r(2,ia)-r(2,ja)-dropbc(2)
                        dz = r(3,ia)-r(3,ja)-dropbc(3)
                        if (ndim == 2) dz = Zero
                        r1 = sqrt(dx**2+dy**2+dz**2)
                        ivar = ipnt(iatjat,igrjgr,itype)
                        ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                        var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
                        tauvar(ivar)%avs2(ibin) = tauvar(ivar)%avs2(ibin) + (dxo*dx+dyo*dy+dzo*dz)/r1  ! tau = dro.dr/r1
                        tauvar(ivar)%nsampbin(ibin) = tauvar(ivar)%nsampbin(ibin) + One
                     end do
                  end do
               end if

            end if
         end do
      end do

   case (iAfterMacrostep)

#if defined (_PAR_)
      do ivar = 1, nvar
          call par_allreduce_reals(var(ivar)%avs2(-1), vaux, var(ivar)%nbin+2)
      end do
#endif

      call ScalarSample(iStage,1,1,var_s)           ! volume

      rdfvols2 = rdfvols2/var(1)%nsamp2
      do igr = 1, ngr(1)
         do jgr = 1, ngr(2)
            igrjgr = igrgr(igr,jgr)
            if (vtype(1)%l) then
               ivar = ipnt(1,igrjgr,1)
               var(ivar)%norm = InvFlt(grvar(igrpnt(1,igr))%avs2)
               if (func == 'rdf') then
                  var(ivar)%norm = var(ivar)%norm*rdfvols2*InvFlt(rdfvolfac*grvar(igrpnt(2,jgr))%avs2)
                  if (igr == jgr) var(ivar)%norm = var(ivar)%norm*(grvar(igrpnt(2,jgr))%avs2/(grvar(igrpnt(2,jgr))%avs2-One))
               end if
               if (l2dtwo) var(ivar)%norm = Two*var(ivar)%norm                ! to get correct normalization
            end if
            if (vtype(2)%l) then
               do jat = iatgr(jgr,2), iatgr(jgr,2)+natgr(jgr,2)-1
                  ivar = ipnt(jat,igrjgr,2)
                  var(ivar)%norm = InvFlt(grvar(igrpnt(1,igr))%avs2)
                  if (func == 'rdf') var(ivar)%norm = var(ivar)%norm*rdfvols2*InvFlt(rdfvolfac*naat(jat)*grvar(igrpnt(2,jgr))%avs2)
               end do
            end if
            if (vtype(3)%l) then
               do iat = iatgr(igr,1), iatgr(igr,1)+natgr(igr,1)-1
                  do jat = iatgr(jgr,2), iatgr(jgr,2)+natgr(jgr,2)-1
                     iatjat = iat+nat*(jat-1)
                     ivar = ipnt(iatjat,igrjgr,3)
                     var(ivar)%norm = InvFlt(grvar(igrpnt(1,igr))%avs2*naat(iat))
                     if (func == 'rdf') var(ivar)%norm = var(ivar)%norm*rdfvols2*InvFlt(rdfvolfac*naat(jat)*grvar(igrpnt(2,jgr))%avs2)
                  end do
               end do
            end if
         end do
      end do

      do ivar = 1, nvar
         if (func == 'rcn') then
            var(ivar)%avs2(0) = var(ivar)%avs2(0)*var(ivar)%norm
            do ibin = 1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1)+var(ivar)%avs2(ibin)*var(ivar)%norm
            end do
         else if (func == 'rdf') then
            do ibin = 0, var(ivar)%nbin
               if (ndim == 3) then
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*var(ivar)%norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
                  tauvar(ivar)%avs2(ibin) = tauvar(ivar)%avs2(ibin)*InvFlt(tauvar(ivar)%nsampbin(ibin))
               end if
               if (ndim == 2) var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*var(ivar)%norm/darea(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
      end do

      call DistFuncSample(iStage, nvar, var)
      call DistFuncSample(iStage, nvar, tauvar)

      if (vtype(3)%l) then
         call PressureContact(var_s(1)%avs2, var_s(2)%value)
         call ScalarSample(iSimulationStep, 2, 2, var_s)
         call ScalarSample(iStage, 2, 2, var_s)
      end if

      if (lsim .and. master) write(ucnf) var, tauvar, var_s

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call DistFuncSample(iStage, nvar, tauvar)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f8.2)')            'upper separation distance      = ', rmax
         write(uout,'(a,t35,i8  )')            'number of dimensions           = ', ndim
         write(uout,'(a,t35,l   )')            'l2dtwo (two separated slabs)   = ', l2dtwo

        if (vtype(3)%l) then
            call ScalarSample(iStage,2,2,var_s)
            write(uout, '(a, t35, 10f8.4)')      'pressure (NkT/V)(cont. cont.) S= ', var_s(2)%avs1, var_s(2)%avsd
            write(uout, '(a, t35, 10f8.4)')      'pressure (NkT/V)(total      ) S= ', prsrreds3+var_s(2)%avs1, sqrt(prsrredsd**2+var_s(2)%avsd**2)
         endif
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
!         call DistFuncWrite('tau', nvar, tauvar, uout, ulist, ishow, iplot, ilist)
!         tauvarvar = tauvar
!         do ivar = 1, nvar
!            do ibin = 1, tauvarvar(ivar)%nbin
!               tauvarvar(ivar)%avs1(ibin) = tauvarvar(ivar)%avs1(ibin) * var(ivar)%avs1(ibin) ! Simplified ? Jos
!            end do
!         end do
!         call DistFuncWrite('tau.gr', nvar, tauvarvar, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   contains

!........................................................................

      subroutine PressureContact(vols, prsrhs)
         use MolModule
         implicit none

         real(8), intent(in)  :: vols
         real(8), intent(out) :: prsrhs

         real(8)     :: denst, densiat, densjat, rcont, gcont, term
         integer(4)  :: i, igr, jgr, iat, jat

         prsrhs = Zero
         do igr = 1, ngr(1)
            do jgr = 1, ngr(2)
               igrjgr = igrgr(igr,jgr)
               do iat = iatgr(igr,1), iatgr(igr,1)+natgr(igr,1)-1
                  densiat=naat(iat)*nppt(iptat(iat))/vols
                  do jat = iatgr(jgr,2), iatgr(jgr,2)+natgr(jgr,2)-1
                     densjat=naat(jat)*nppt(iptat(jat))/vols
                     iatjat = iat+nat*(jat-1)
                     ivar = ipnt(iatjat,igrjgr,3)
                     rcont = radat(iat)+radat(jat)
                     call gcontact(var(ivar)%min, var(ivar)%bin, var(ivar)%avs2(-1), tauvar(ivar)%avs2(-1), rcont, gcont) ! tauvar(ivar)%min and tauvar(ivar)%bin are not needed and therefore removed
                     term = rcont**2*densiat*densjat*gcont
                     prsrhs = prsrhs + term
                  end do
               end do
            end do
         end do

         denst = sum([ (nppt(i), i=1, npt) ])/vols
         prsrhs = TwoPi/(Three*denst)*prsrhs
      end subroutine PressureContact

      !subroutine GContact(vlow, vbin, vs2, tlow, tbin, tauvar, rcont, gcont) ! tlow and tbin is not used
      subroutine GContact(vlow, vbin, vs2, tauvar, rcont, gcont)
         implicit none
     real(8), intent(in) :: vlow
     real(8), intent(in) :: vbin
     real(8), intent(in) :: vs2(-1:*)
     !real(8), intent(in) :: tlow
     !real(8), intent(in) :: tbin
     real(8), intent(in) :: tauvar(-1:*)
         real(8), intent(in) :: rcont
         real(8), intent(out):: gcont

         integer(4), parameter :: npol = 2, ndp = npol+1
          integer(4) :: i, ihs

     real(8) :: xx(ndp), yy(ndp), tt(ndp), ww(ndp), a(0:npol), dum1, dum2
     real(8), external :: PolVal

         ihs = int((rcont-vlow-1.0d-10)/vbin)
         xx(1:ndp) = vlow + vbin * ( 0.5d0+ [ ( ihs+i, i=1, ndp )] )
         yy(1:ndp) = [ (vs2(ihs+i),i=1,ndp ) ]
         tt(1:ndp) = [ (tauvar(ihs+i),i=1,ndp ) ]
         yy = yy*tt
         ww(1:ndp) = 1.0d0
         if (count(yy(1:ndp) < 1.0d-10) > 0) then
            gcont = 0.0d0
         else if (yy(1)>0.5d0) then
            yy(1:ndp) = log(yy(1:ndp))
            call PolFit(npol, ndp, xx, yy, ww, 0, 6, a, dum1, dum2)
            gcont = exp(PolVal(npol, a, rcont))
         else
            call PolFit(npol, ndp, xx, yy, ww, 0, 6, a, dum1, dum2)
            gcont = PolVal(npol, a, rcont)
         end if
         gcont = max(0.0d0, gcont)
      end subroutine GContact

!........................................................................

end subroutine RDF

!************************************************************************
!> \page static static.F90
!! **RDFChain**
!! *calculate running coordination number or radial distribution function*
!************************************************************************


!> \page nmlRDFChain
!! The namelist  \ref nmlRDFChain contains variables that control the calculation of running coordination number or radial distribution
!! functions. Any combination of the types of distribution functions listed below may be selected through vtype\%l.
!!    | type | label   | quantity                      |
!!    | :--: | :-----: | :---------------------------: |
!!    | 1    | com-com | center of mass-center of mass |
!! Only chains within a separation of rmax are considered.
!! * Variables:
!!  * \subpage nmlRDFChain_vtype
!!  * \subpage nmlRDFChain_rmax
!!  * \subpage nmlRDFChain_func

!> \page nmlRDFChain_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)` (1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /0.0/ Max: /10.0/

!> \page nmlRDFChain_rmax rmax
!! `real`
!! **default:** `10.0`
!! * Upper distance for particle separation considered.

!> \page nmlRDFChain_func func
!! `character(3)`
!! **default:** `rdf`
!! * `rdf`: Radial distribution functions are calculated.
!! * `rcn`: Running coordination numbers are calculated.
subroutine RDFChain(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='RDFChain'
   character(80), parameter :: txheading ='running coordination number/radial distribution function for chains'
   integer(4)   , parameter :: ntype = 1
   type(static1D_var),         save :: vtype(ntype)
   real(8),                    save :: rmax
   character(3),               save :: func
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:)
   real(8),       allocatable, save :: rcom(:,:)
   real(8), save :: rdfvols2
   type(chainprop_var) :: ChainProperty
   integer(4) :: itype, ivar, ibin, ic, jc, ict, jct, ictjct
   real(8)    :: r1, r2, dvol, dcx, dcy, dcz, dropbc(3), InvFlt

   namelist /nmlRDFChain/ vtype, rmax, func

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = Zero
      vtype%max = 10.0d0
      vtype%nbin= 100
      rmax   = 10.0d0
      func   = 'rdf'

      rewind(uin)
      read(uin,nmlRDFChain)

      call LowerCase(func)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, ' vtype%nbin > mnbin_df', uout)
      if ((func /= 'rdf') .and. (func /= 'rcn')) call Stop(txheading, 'check value of func', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%nvar = [(nct*(nct+1))/2]

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nctct), rcom(3,nc))
      ipnt = 0
      rcom = 0.0E+00

      ivar = 0
      itype = 1
      do ict = 1, nct
         do jct = ict, nct
            ictjct = ictct(ict,jct)
            ivar = ivar+1
            ipnt(ictjct) = ivar
            var(ivar)%label = trim(txctct(ictjct))//' '//'com'//'-'//'com'
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end do
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      rdfvols2 = Zero

   case (iSimulationStep)

      rdfvols2 = rdfvols2+vol

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc                              ! get com for all chains
         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)
         rcom(1:3,ic) = ChainProperty%ro(1:3)
      end do

      do ic = 1, nc                              ! loop over all pairs of chains
         ict = ictcn(ic)
         do jc = ic+1, nc
            jct = ictcn(jc)
            dcx = rcom(1,ic) - rcom(1,jc)
            dcy = rcom(2,ic) - rcom(2,jc)
            dcz = rcom(3,ic) - rcom(3,jc)
            call PBC2(dcx,dcy,dcz,dropbc(1),dropbc(2),dropbc(3))
            dcx = dcx-dropbc(1)
            dcy = dcy-dropbc(2)
            dcz = dcz-dropbc(3)
            r2 = dcx**2 + dcy**2 + dcz**2
            if (r2 < rmax**2) then

! ... sample type 1

               itype = 1
               if (vtype(itype)%l) then
                  r1 = sqrt(r2)
                  ictjct = ictct(ict,jct)
                  ivar = ipnt(ictjct)
                  ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if

            end if
         end do
      end do

   case (iAfterMacrostep)

      rdfvols2 = rdfvols2/var(1)%nsamp2
      do ict = 1, nct
         do jct = 1, nct
            ictjct = ictct(ict,jct)
            if (vtype(1)%l) then
               ivar = ipnt(ictjct)
               var(ivar)%norm = InvInt(ncct(ict))
               if (func == 'rdf') var(ivar)%norm = var(ivar)%norm*rdfvols2*InvFlt(FourPiThird*ncct(jct))
               if (ict == jct) var(ivar)%norm = var(ivar)%norm*Two
            end if
         end do
      end do

      do ivar = 1, nvar
         if (func == 'rcn') then
            var(ivar)%avs2(0) = var(ivar)%avs2(0)*var(ivar)%norm
            do ibin = 1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1)+var(ivar)%avs2(ibin)*var(ivar)%norm
            end do
         else if (func == 'rdf') then
            do ibin = 0, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*var(ivar)%norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end if
      end do

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,f8.2)')            'upper separation distance      = ', rmax
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, ipnt, rcom)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine RDFChain

!************************************************************************
!> \page static static.F90
!! **RDFSph**
!! *calculate rcn or rdf of particles which positions are projected on a sphere*
!************************************************************************


! rcn(r*t') = 2*pi*r**2*sigma*int(0 to t')(g(t)*sin(t)dt), sigma = n/(4*pi*r**2)

!> \page nmlRDFSph
!! The namelist  \ref nmlRDFSph contains variables that control the calculation of running coordination number or radial distribution functions of particles which positions are projected on a sphere.
!! * Variables:
!!  * \subpage ipsph
!!  * \subpage iptrcnsph
!!  * \subpage jptrcnsph
!!  * \subpage nmlRDFSph_nbin
!!  * \subpage nmlRDFSph_func

!> \page ipsph
!> `integer`
!! * The identity of the particle on which surface the projection is made.

!> \page iptrcnsph
!! `integer`
!! * Type of particle for which distribution function should be calculated.

!> \page jptrcnsph
!! `integer`
!! * Type of particle for which distribution function should be calculated.

!> \page nmlRDFSph_nbin nbin
!! `integer`
!! **default:** `100`
!! * Number of bins used to sample the distribution functions.

!> \page nmlRDFSph_func func
!! `character(3)`
!! **default:** `rdf`
!! * `rdf`: Radial distribution functions are calculated.
!! * `rcn`: Running coordination numbers are calculated.

subroutine RDFSph(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(80), parameter :: txheading ='run/rdf of particles which positions are projected on a sphere'
   character(40), parameter :: txroutine ='RDFSph'
   character(3),  save :: func
   integer(4),    save :: nbin, nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    save :: ipsph                   ! particle on which the projection should be made
   integer(4),    save :: iptrdfsph, jptrdfsph    ! particle types for rdf
   real(8),       save :: radsph, radsphi
   integer(4) :: ivar, ibin, ip, jp
   real(8) :: dx1, dy1, dz1, dx2, dy2, dz2, theta, thetalow, thetaupp, arclen, fac, angcos

   namelist /nmlRDFSph/ ipsph, iptrdfsph, jptrdfsph, nbin, func

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (slave) return                   ! only master

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      nbin    = 100
      func    = 'rdf'

      rewind(uin)
      read(uin,nmlRDFSph)

      call LowerCase(func)

      if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)
      if (iptrdfsph == iptpn(ipsph)) call Stop(txroutine, 'error in input', uout)
      if (jptrdfsph == iptpn(ipsph)) call Stop(txroutine, 'error in input', uout)
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)
      if ((func /= 'rdf') .and. (func /= 'rcn')) call Stop(txheading, 'check value of func', uout)
      radsph = radat(iptpn(ipsph)) + Half*(radat(iptrdfsph)+radat(jptrdfsph))
      radsphi = One/radsph

   case (iWriteInput)

! ... set nvar as well as allocate memory

      nvar = 1
      allocate(var(nvar))

      var(nvar)%label = trim(txpt(iptrdfsph))//"-"//trim(txpt(jptrdfsph))
      var(nvar)%min = Zero
      var(nvar)%max = Pi * radsph
      var(nvar)%nbin = nbin
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = ipnpt(iptrdfsph), ipnpt(iptrdfsph) + nppt(iptrdfsph) - 1
         dx1 = ro(1,ip)-ro(1,ipsph)
         dy1 = ro(2,ip)-ro(2,ipsph)
         dz1 = ro(3,ip)-ro(3,ipsph)
         call PBC(dx1,dy1,dz1)
         do jp = ipnpt(jptrdfsph), ipnpt(jptrdfsph) + nppt(jptrdfsph) - 1
            if (jp == ip) cycle
            dx2 = ro(1,jp)-ro(1,ipsph)
            dy2 = ro(2,jp)-ro(2,ipsph)
            dz2 = ro(3,jp)-ro(3,ipsph)
            call PBC(dx2,dy2,dz2)
            theta = acos(angcos(dx1,dy1,dz1,dx2,dy2,dz2))
            arclen = radsph*theta              !arc length between projected particle ip and jp
            ivar = 1
            ibin = max(-1,min(floor(var(ivar)%bini*(arclen-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end do

   case (iAfterMacrostep)

      ivar = 1
      var(ivar)%norm = One/nppt(iptrdfsph)
      if (func == 'rdf') var(ivar)%norm = var(ivar)%norm*(FourPi*radsph**2/nppt(jptrdfsph))

      if (func == 'rcn') then
         var(ivar)%avs2(0) = var(ivar)%avs2(0)*var(ivar)%norm
         do ibin = 1, var(ivar)%nbin
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1)+var(ivar)%avs2(ibin)*var(ivar)%norm
         end do
      else if (func == 'rdf') then
         do ibin = 0, var(ivar)%nbin
            thetalow = (var(ivar)%min+ibin*var(ivar)%bin)*radsphi
            thetaupp = (var(ivar)%min+(ibin+1)*var(ivar)%bin)*radsphi
            fac = TwoPi*radsph**2*(cos(thetalow) - cos(thetaupp))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*var(ivar)%norm/fac
         end do
      end if

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,i8  )')            'ipsph                          = ', ipsph
      write(uout,'(a,t35,i8  )')            'iptrdfsph                      = ', iptrdfsph
      write(uout,'(a,t35,i8  )')            'jptrdfsph                      = ', jptrdfsph
      write(uout,'(a,t35,i5,5x,a,i5,2x,a)') 'number of grid points          = ', nbin, '(',mnbin_df,')'
      write(uout,'(a,t35,f8.3)')            'projected radius               = ', radsph
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine RDFSph

!************************************************************************
!> \page static static.F90
!! **RDFCond**
!! *calculate conditional radial distribution function*
!************************************************************************


!> \page nmlRDFCond
!! The namelist  \ref nmlRDFCond contains variables that control the calculation of conditional radial distribution functions. The
!! distribution functions are made for five different bins of arccos(theta) making up -1 to 1.
!! * Variables:
!!  * \subpage nmlRDFCond_vtype
!!  * \subpage nmlRDFCond_rmax

!> \page nmlRDFCond_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`
!! Min: /0.0/ Max: /10.0/
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.

!>\page nmlRDFCond_rmax rmax
!! `real`
!! **default:** `10.0`
!! * Upper distance for particle separation considered.
subroutine RDFCond(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='RDFCond'
   character(80), parameter :: txheading ='Conditional rdf [g(r;abs(cos(z''*r)))]'
   integer(4)   , parameter :: ntype = 5
   type(static1D_var),         save :: vtype(ntype)
   real(8),                    save :: rmax
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),       save :: rdfvols2, rdfvolfac
   integer(4) :: itype, ivar, ibin, ip, jp
   real(8)    :: dx, dy, dz, dropbc(3), r1, r2, dvol, norm, angcos, ac, InvFlt

   namelist /nmlRDFCond/ vtype, rmax

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.true.
      vtype%min = Zero
      vtype%max = 10.0d0
      vtype%nbin = 100
      rmax = 10.0d0

      rewind(uin)
      read(uin,nmlRDFCond)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

      rdfvolfac = FourPiThird

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['0.0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']
      vtype%nvar = 1

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(1,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            ivar = ivar + 1
            ipnt(1,itype) = ivar
            var(ivar)%label = vtype(itype)%label
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      rdfvols2 = Zero

   case (iSimulationStep)

      rdfvols2 = rdfvols2+vol
      var%nsamp2 = var%nsamp2+1

      do ip = ipmyid(1), ipmyid(2)
         do jp = 1, np
            if (jp == ip) cycle
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC2(dx,dy,dz,dropbc(1),dropbc(2),dropbc(3))
            dx = dx-dropbc(1)
            dy = dy-dropbc(2)
            dz = dz-dropbc(3)
            r2 = dx**2+dy**2+dz**2
            ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),-dx,-dy,-dz)
            ivar = max(0,min(int(1 + nvar*abs(ac)), ntype))
            if (r2 < rmax**2) then
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
         end do
      end do

   case (iAfterMacrostep)

#if defined (_PAR_)
      do ivar = 1, nvar
          call par_allreduce_reals(var(ivar)%avs2(-1), vaux, var(ivar)%nbin+2)
      end do
#endif

      rdfvols2 = rdfvols2/var(1)%nsamp2
      do ivar = 1, nvar
         norm = rdfvols2*InvFlt(rdfvolfac*(real(np))**2/nvar)
         do ibin = 0, var(ivar)%nbin
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
         end do
      end do

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f8.2)') 'upper separation distance      = ', rmax
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine RDFCond

!************************************************************************
!> \page static static.F90
!! **G3Dist**
!! *calculate g(r12,r13,r23)/g(r12)*
!************************************************************************


!     instead of r23, use cos(theta) where theta is the angle between r12 and r13
!     no n/(n-1) etc correction

!     Note, data are not saved after each macrostep

!> \page nmlG3Dist
!! The namelist  \ref nmlG3Dist contains variables that control the calculation of normalized triplet correlation functions g(r12, r13,
!! r23)/g(r12). r23 is represented by the angle theta formed by the vectors r12 and r13.
!! * Variables:
!!  * \subpage ipt1
!!  * \subpage ipt2
!!  * \subpage ipt3
!!  * \subpage snbin
!!  * \subpage slow
!!  * \subpage supp
!!  * \subpage tnbin
!!  * \subpage tlow
!!  * \subpage tupp
!!  * \subpage anbin
!!  * \subpage alow
!!  * \subpage aupp
!!  * \subpage nskip

!> \page ipt1
!! `integer`
!! * Particle type of the first particle.

!> \page ipt2
!! `integer`
!! * Particle type of the second particle.

!> \page ipt3
!! `integer`
!! * Particle type of the third particle.

!> \page snbin
!! `integer`
!! * Number of bins to sample the distance r12.

!> \page slow
!! `real`
!! * Lower end of r12 to be sampled.

!> \page supp
!! `real`
!! * Upper end of r12 to be sampled.

!> \page tnbin
!! `integer`
!! * Number of bins to sample the distance r13.

!> \page tlow
!! `real`
!! * Lower end of r13 to be sampled.

!> \page tupp
!! `real`
!! * Upper end of r13 to be sampled.

!> \page anbin
!! `integer`
!! * Number of bins to sample the angle theta.

!> \page alow
!! `real`
!! * Lower end of theta to be sampled.

!> \page aupp
!! `real`
!! * Upper end of theta to be sampled.

!> \page nskip
!! `integer`
!! **default:** `1`
!! * Interval of listing the distance r12.

subroutine G3Dist(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='G3Dist'
   character(80), parameter :: txheading ='three-particle distribution functions'
   real(8), allocatable :: g2(:), g3(:,:,:)
   integer(4) :: snbin, tnbin, anbin, nsamp1
   integer(4) :: nskip
   integer(4) :: ip1, ipt1, ip2, ipt2, ip3, ipt3
   integer(4) :: ip1low, ip1upp, ip2low, ip2upp, ip3low, ip3upp
   integer(4) :: is, it, ia
   real(8)    :: norm, anorm, t, a, costh, dvol
   real(8)    :: xo1, yo1, zo1, xo2, yo2, zo2, xo3, yo3, zo3
   real(8)    :: dx12, dy12, dz12, dx13, dy13, dz13, r12, r13, r212, r213, r223
   real(8)    :: dxpbc, dypbc, dzpbc
   real(8)    :: slow, supp, sbin, sbini, supp2
   real(8)    :: tlow, tupp, tbin, tbini, tupp2
   real(8)    :: alow, aupp, abin, abini
   save

   namelist /nmlG3Dist/ ipt1, ipt2, ipt3, snbin, slow, supp, tnbin, tlow, tupp, anbin, alow, aupp, nskip

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      nskip = 1

      rewind(uin)
      read(uin,nmlG3Dist)

   case (iWriteInput)

      ip1low = ipnpt(ipt1)
      ip1upp = ipnpt(ipt1)+nppt(ipt1)-1
      ip2low = ipnpt(ipt2)
      ip2upp = ipnpt(ipt2)+nppt(ipt2)-1
      ip3low = ipnpt(ipt3)
      ip3upp = ipnpt(ipt3)+nppt(ipt3)-1

      sbin  = (supp-slow)/snbin
      sbini = One/sbin
      supp2 = supp**2

      tbin  = (tupp-tlow)/tnbin
      tbini = One/tbin
      tupp2 = tupp**2

      abin  = (aupp-alow)/anbin
      abini = One/abin

      allocate(g2(-1:snbin), g3(-1:anbin,-1:tnbin,-1:snbin))
      g2 = 0.0E+00
      g3 = 0.0E+00

   case (iBeforeSimulation)

      nsamp1 = Zero
      g2(-1:snbin) = Zero
      g3(-1:anbin,-1:tnbin,-1:snbin) = Zero

   case (iSimulationStep)

      nsamp1 = nsamp1+1
      do ip1 = ip1low, ip1upp
         xo1 = ro(1,ip1)
         yo1 = ro(2,ip1)
         zo1 = ro(3,ip1)
         do ip2 = ip2low, ip2upp
            if (ip2 == ip1) cycle
            xo2 = ro(1,ip2)
            yo2 = ro(2,ip2)
            zo2 = ro(3,ip2)

! ... position particle 2 in the box of particle 1

            dx12 = xo1-xo2
            dy12 = yo1-yo2
            dz12 = zo1-zo2
            call PBC2(dx12,dy12,dz12,dxpbc,dypbc,dzpbc)
            xo2 = xo2 + dxpbc
            yo2 = yo2 + dypbc
            zo2 = zo2 + dzpbc
            r212 = (xo1-xo2)**2+(yo1-yo2)**2+(zo1-zo2)**2
            if (r212 < supp2) then
               r12 = sqrt(r212)

               is = max(-1,min(floor(sbini*(r12-slow)),int(snbin)))
               g2(is) = g2(is) + One

               do ip3 = ip3low, ip3upp
                  if (ip3 == ip1 .or. ip3 == ip2) cycle
                  xo3 = ro(1,ip3)
                  yo3 = ro(2,ip3)
                  zo3 = ro(3,ip3)

! ... position particle 3 in the box of particle 1

                  dx13 = xo1-xo3
                  dy13 = yo1-yo3
                  dz13 = zo1-zo3
                  call PBC2(dx13,dy13,dz13,dxpbc,dypbc,dzpbc)
                  xo3 = xo3 + dxpbc
                  yo3 = yo3 + dypbc
                  zo3 = zo3 + dzpbc
                  r213 = (xo1-xo3)**2+(yo1-yo3)**2+(zo1-zo3)**2
                  if (r213 < tupp2) then
                     r13 = sqrt(r213)
                     r223 = (xo2-xo3)**2+(yo2-yo3)**2+(zo2-zo3)**2
                     costh = (r213+r212-r223)/(2.0*r13*r12)
                     it = max(-1,min(floor(tbini*(r13-tlow)),int(tnbin)))
                     ia = max(-1,min(floor(abini*(costh-alow)),int(anbin)))
                     g3(ia,it,is) = g3(ia,it,is) + One
                  end if
               end do
            end if
         end do
      end do

   case (iAfterSimulation)

      do is = 0, snbin
         norm = Zero
         if (g2(is) /= 0) norm = volst/(nppt(ipt3)*g2(is))
         do it = 0, tnbin
            anorm = norm/(TwoPi*(tlow+(it+0.5)*tbin)**2*tbin*abin)
            g3(0:anbin,it,is) = g3(0:anbin,it,is)*anorm
         end do
         norm = volst/(FourPiThird*nppt(ipt1)*nppt(ipt2)*nsamp1)
         g2(is) = g2(is)*norm/dvol(is,slow,sbin)
      end do

! ............... write output  ..............

      call WriteHead(2, txheading, uout)
      write(uout,'(a,a)') 'particle type 1: ', txpt(ipt1)
      write(uout,'(a,a)') 'particle type 2: ', txpt(ipt2)
      write(uout,'(a,a)') 'particle type 3: ', txpt(ipt3)
      write(uout,'()')
      write(uout,'(a)') 'variable   nbin      lower       upper         bin        bini'
      write(uout,'(a)') '--------   ----      -----       -----         ---        ----'
      write(uout,'(4x,a,i10,4f12.3)')'s', snbin, slow, supp, sbin, sbini
      write(uout,'(4x,a,i10,4f12.3)')'r', tnbin, tlow, tupp, tbin, tbini
      write(uout,'(4x,a,i10,4f12.3)')'a', anbin, alow, aupp, abin, abini

      if (ishow/= 0) then
         write(uout,'()')
         write(uout,'(a)') 'g2(r12)'
         write(uout,'(a)') '-------'
         do is = 0, snbin-1
            write(uout,'(2f10.4)') slow+(is+0.5)*sbin, g2(is)
         end do
         write(uout,'()')
         write(uout,'(a)') 'g3(r12,t,a)'
         write(uout,'(a)') '-----------'
         do is = 0, snbin-1, nskip
            write(uout,'()')
            write(uout,'(a,f10.4)') 'r12 = ', slow+(is+0.5)*sbin
            write(uout,'()')
            write(uout,'(a)')  't   a = cos(theta)   g3(r12,t,a)'
            write(uout,'(a)')  '------------------------------'
            do it = 0, tnbin-1
               t = tlow+(it+0.5)*tbin
               do ia = 0, anbin-1
                  a = alow+(ia+0.5)*abin
                  write(uout,'(3f10.4)') t, a, g3(ia,it,is)
               end do
            end do
         end do
      end if

      if (ilist/= 0) then
         write(ulist,*) 1+(snbin-1)/nskip, tnbin, anbin
         do is = 0, snbin-1, nskip
            write(ulist,'(f10.4)') slow+(is+0.5)*sbin
            do it = 0, tnbin-1
               t = tlow+(it+0.5)*tbin
               do ia = 0, anbin-1
                  a = alow+(ia+0.5)*abin
                  write(ulist,'(3f10.4)') t, a, g3(ia,it,is)
               end do
            end do
         end do
      end if

      deallocate(g2, g3)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine G3Dist

!************************************************************************
!> \page static static.F90
!! **SFDriver**
!! *calculate partial structure factors*
!************************************************************************


subroutine SFDriver(iStage)
   use MolModule
   implicit none
   character(40), parameter :: txroutine ='SFDriver'

   integer(4), intent(in) :: iStage

   if (txbc == 'xyz' .or. txbc == 'xy') then
      call SFPBC(iStage)
   else if (lbcsph .or. lbccyl .or. lbcell) then
      call SFNoPBC(iStage)
   else
      call Stop(txroutine,'Unsupported boundary condition',uout)
   end if

end subroutine SFDriver

!************************************************************************
!> \page static static.F90
!! **SFPBC**
!! *calculate partial structure factors*
!************************************************************************


!     calculate scattering intensities

!     type  label    quantity
!     ----  -----    --------
!     1     100      sf in the 100 direction
!     2     110      sf in the 110 direction
!     3     111      sf in the 111 direction

!> \page nmlSF
!! The namelist  \ref nmlSF contains variables that control the calculation of partial structure factors.
!! * (i.) \ref txbc='sph' or \ref txbc='cyl' Variables:
!!  * \subpage txkscale
!!  * \subpage klow
!!  * \subpage logklow
!!  * \subpage logkupp
!!  * \subpage nmlSFnoPBC_nbin
!! * (ii.) Cubic box and \ref txbc='xyz'. The largest k-vector is 2Pi/box. Variables:
!!  * \subpage nmlSF_ndim
!!  * \subpage nmlSF_nbin
!!  * \subpage lqsorted
!!  * \subpage lsi

!> \page nmlSF_ndim ndim
!! `integer`
!! **default:** `3`
!! * `2`:  Structure factor in the xy-plane. It is averaged over 2 100 and 2 110 directions.
!! * `3`:  Structure factor in the xyz-space. It is averaged over 3 100, 6 110 directions, and 4 110 directions.

!> \page nmlSF_nbin nbin
!! `integer`
!! **default:** `100`
!! * Number of bins used to sample the partial structure factors.

!> \page lqsorted
!! `logical`
!! **default:** `.true.`
!! * `.false.`:  List separately the structure factors of the different types directions.
!! * `.true.`:  Sort the structure factors of the different types of directions and list the sorted structure factor.

!> \page lsi
!! `logical`
!! **default:** `.false.`
!! * `.false.`:  Nothing.
!! * `.true.`:  Scattering intensities are calculated. Further specification is given in \ref nmlScatIntens.

subroutine SFPBC(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SFPBC'
   character(80), parameter :: txheading ='structure factor'
   integer(4)   , parameter :: ntype = 3
   character(14), parameter :: txtype(ntype) = ['100 directions', '110 directions', '111 directions']
   real(8)      , parameter :: sffac(ntype) = [ One, sqrt(Two), sqrt(Three) ]
   integer(4)   , save :: dirlow(ntype), dirupp(ntype)
   integer(4)   , save :: ndim, nbin, nvar, nvartype(ntype)
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   logical      , save :: lqsorted, lsi
   integer(4)   , save :: ndiv = 1   ! factor by which the number of bins is reduced (in FLIST)
   integer(4) :: ip, ipt, jpt, ivar, iptjpt, ibin, jbin(ntype), itype, nq, m, dupp, dlow
   real(8)    :: norm, kx, ky, kz
   real(8)   , allocatable, save :: q(:), sfpar(:,:), sfparsd(:,:)
   complex(8), allocatable, save :: eikr(:,:,:), eikr1(:)
   complex(8) :: eikrtemp

   namelist /nmlSF/ ndim, nbin, lqsorted, lsi

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ndim = 3
      nbin = 100
      lqsorted =.true.
      lsi =.false.

      rewind(uin)
      read(uin,nmlSF)

      if (ndim == 3) then
         dirlow(1:ndim) = [ 1, 4, 10 ]
         dirupp(1:ndim) = [ 3, 9, 13 ]
         if (boxlen(1) /= boxlen(2)) call Stop(txroutine, 'boxlen(1) /= boxlen(2)', uout)
         if (boxlen(1) /= boxlen(3)) call Stop(txroutine, 'boxlen(1) /= boxlen(3)', uout)
      else if (ndim == 2) then
         dirlow(1:ndim) = [ 1, 3 ]
         dirupp(1:ndim) = [ 2, 4 ]
         if (boxlen(1) /= boxlen(2)) call Stop(txroutine, 'boxlen(1) /= boxlen(2)', uout)
      else
         call Stop(txroutine, 'ndim out of range', uout)
      end if
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

! ... allocate memory

      allocate(eikr(0:nbin,13,npt), eikr1(13), q(3*nbin), sfpar(3*nbin,nptpt), sfparsd(3*nbin,nptpt))
      eikr = cmplx(Zero,Zero)
      eikr1 = cmplx(Zero,Zero)
      q = 0.0E+00
      sfpar = 0.0E+00
      sfparsd = 0.0E+00

   case (iWriteInput)

! ... set nvartype and nvar as well as allocate memory

      nvartype(1:ntype) = [(npt*(npt+1))/2, npt*(npt+1)/2, npt*(npt+1)/2]
      nvar = sum(nvartype(1:ndim))
      allocate(var(nvar), ipnt(nptpt,ndim))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1,ndim
         do ipt = 1, npt
            do jpt = ipt, npt
               iptjpt = iptpt(ipt,jpt)
               ivar = ivar + 1
               ipnt(iptjpt,itype) = ivar
               var(ivar)%label = trim(txptpt(iptjpt))//' '//txtype(itype)
               var(ivar)%min = Half*TwoPiBoxi(1)*sffac(itype)
               var(ivar)%max = (nbin+Half)*TwoPiBoxi(1)*sffac(itype)
               var(ivar)%nbin = nbin
            end do
         end do
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      eikr(0:nbin,1:dirupp(ndim),1:npt) = cmplx(Zero,Zero)

! ... calculate sum( exp(i (k * r ))
!                i               i

      do ip = 1, np
         ipt = iptpn(ip)
         kx = TwoPiBoxi(1)*ro(1,ip)
         ky = TwoPiBoxi(2)*ro(2,ip)
         kz = TwoPiBoxi(3)*ro(3,ip)
         if (ndim == 3) then
            eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
            eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
            eikr1(3) = cmplx(cos(kz),sin(kz))              ! ( 0 0 1) direction
            eikr1(4) = cmplx(cos(kx+ky),sin(kx+ky))        ! ( 1 1 0) direction
            eikr1(5) = cmplx(cos(ky+kz),sin(ky+kz))        ! ( 0 1 1) direction
            eikr1(6) = cmplx(cos(kz+kx),sin(kz+kx))        ! ( 1 0 1) direction
            eikr1(7) = cmplx(cos(kx-ky),sin(kx-ky))        ! ( 1-1 0) direction
            eikr1(8) = cmplx(cos(ky-kz),sin(ky-kz))        ! ( 0 1-1) direction
            eikr1(9) = cmplx(cos(kz-kx),sin(kz-kx))        ! (-1 0 1) direction
            eikr1(10) = cmplx(cos(kx+ky+kz),sin(kx+ky+kz)) ! ( 1 1 1) direction
            eikr1(11) = cmplx(cos(kx+ky-kz),sin(kx+ky-kz)) ! ( 1 1-1) direction
            eikr1(12) = cmplx(cos(kx-ky+kz),sin(kx-ky+kz)) ! ( 1-1 1) direction
            eikr1(13) = cmplx(cos(kz+ky-kx),sin(kz+ky-kx)) ! (-1 1 1) direction
         else if (ndim == 2) then
            eikr1(1) = cmplx(cos(kx),sin(kx))              ! ( 1 0 0) direction
            eikr1(2) = cmplx(cos(ky),sin(ky))              ! ( 0 1 0) direction
            eikr1(3) = cmplx(cos(kx+ky),sin(kx+ky))        ! ( 1 1 0) direction
            eikr1(4) = cmplx(cos(kx-ky),sin(kx-ky))        ! ( 1-1 0) direction
         end if
         do m = 1, dirupp(ndim)
            eikrtemp = cmplx(One,Zero)
            do ibin = 0, nbin-1
               eikrtemp = eikrtemp*eikr1(m)
               eikr(ibin,m,ipt) = eikr(ibin,m,ipt)+eikrtemp
            end do
         end do
      end do

! ... sum over different equivalent k-vectors

      do itype = 1,ndim
         dlow=dirlow(itype)
         dupp=dirupp(itype)
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(iptpt(ipt,jpt),itype)
               do ibin = 0, nbin-1
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ &
                  sum(real(eikr(ibin,dlow:dupp,ipt))*real(eikr(ibin,dlow:dupp,jpt))  &
                      +aimag(eikr(ibin,dlow:dupp,ipt))*aimag(eikr(ibin,dlow:dupp,jpt)))
               end do
            end do
         end do
      end do

   case (iAfterMacrostep)

! ... normalizing

      do itype = 1,ndim
         do ipt = 1, npt
            do jpt = ipt, npt
               ivar = ipnt(iptpt(ipt,jpt),itype)
               norm = One/((dirupp(itype)-dirlow(itype)+1)*sqrt(real(nppt(ipt))*real(nppt(jpt))))
               var(ivar)%avs2(-1:nbin) = var(ivar)%avs2(-1:nbin)*norm
            end do
         end do
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      if (ndim == 3) write(uout,'(a,t35,i5,5x,a)') 'number of dimensions           = ', ndim, &
       'averaged over 3x(100), 6x(110), and 4x(111) directions'
      if (ndim == 2) write(uout,'(a,t35,i5,5x,a)') 'number of dimensions           = ', ndim, &
       'averaged over 2x(100) and 2x(110) directions'
      write(uout,'(a,t35,i5,5x,a,i5,2x,a)') 'number of grid points          = ', nbin, '(',mnbin_df,')'
      write(uout,'(a,t35,g10.3)')           'q-sorted                       = ', lqsorted
      write(uout,'(a,t35,f8.2)')            'lower wavevector limit         = ', var(1)%min+0.5*TwoPiBoxi(1)
      write(uout,'(a,t35,f8.2)')            'upper wavevector limit         = ', var(1)%max-0.5*TwoPiBoxi(1)

      if (.not.lqsorted) then

         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      else

! ... store sf data in q-sorted order

         do ipt = 1, npt
            do jpt = ipt, npt
               iptjpt = iptpt(ipt,jpt)
               jbin(1:ndim) = 1
               nq = 0
               if (ndim == 3) then
                  do ibin = 1, ndim*nbin
                     if (jbin(1)*sffac(1) < jbin(2)*sffac(2) .and. jbin(1)*sffac(1) < jbin(3)*sffac(3) .and. jbin(1) <= nbin) then
                        call sortsub(1)
                     else if (jbin(2)*sffac(2) < jbin(3)*sffac(3) .and. jbin(2) <= nbin) then
                        call sortsub(2)
                     else if (jbin(3) <= nbin) then
                        call sortsub(3)
                     end if
                  end do
               else if (ndim == 2) then
                  do ibin = 1, ndim*nbin
                     if (jbin(1)*sffac(1) < jbin(2)*sffac(2) .and. jbin(1) <= nbin) then
                        call sortsub(1)
                     else if (jbin(2) <= nbin) then
                        call sortsub(2)
                     end if
                  end do
               end if
            end do
         end do

! ... write sf data in q-sorted order

         if (ishow > 0) then
            write(uout,'()')
            write(uout,'(a)') 'q-sorted structure factor'
            do ipt = 1, npt
               do jpt = ipt, npt
                  iptjpt = iptpt(ipt,jpt)
                  write(uout,'(a)') txptpt(iptjpt)
                  write(uout,'(i5,3g15.5)') (ibin,q(ibin),sfpar(ibin,iptjpt),sfparsd(ibin,iptjpt),ibin = 1,nq)
               end do
            end do
         end if
         if (ilist > 0) then
            write(ulist,*) 'structure factor'
            write(ulist,*) nptpt
            do ipt = 1, npt
               do jpt = ipt, npt
                  iptjpt = iptpt(ipt,jpt)
                  write(ulist,*) 'sf '//txptpt(iptjpt)
                  write(ulist,*) nq/ndiv
                  write(ulist,'(g15.5,a,g15.5,a,g15.5)') (q(ibin),tab,sfpar(ibin,iptjpt),tab,sfparsd(ibin,iptjpt),ibin = 1,nq/ndiv)
               end do
            end do
         end if

         if (lsi) call ScatIntens(nq, q, sfpar)

      end if

      deallocate(eikr, eikr1, q, sfpar, sfparsd)
      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine sortsub(itype)
   integer(4), intent(in) :: itype                           ! type of direction
   ivar = ipnt(iptjpt,itype)                                 ! get variable id
   nq = nq+1                                                 ! number of items in q-sorted list
   q(nq) = jbin(itype)*TwoPiBoxi(itype)*sffac(itype)         ! get wave vector
   sfpar(ibin,iptjpt) = var(ivar)%avs1(jbin(itype)-1)        ! store value of sf
   sfparsd(ibin,iptjpt) = var(ivar)%avsd(jbin(itype)-1)      ! store uncertainty of sf
   jbin(itype) = jbin(itype) + 1                             ! update jbin
end subroutine sortsub

!........................................................................

end subroutine SFPBC

!************************************************************************
!> \page static static.F90
!! **ScatIntens**
!! *calculate scattering intensity, form factor, and structure factor*
!************************************************************************


!> \page nmlScatIntens
!! The namelist  \ref nmlScatIntens contains variables that control the calculation of the scattering intensities using a multi shell profile of constant scattering properties of each particle.
!! * Variables:
!!  * \subpage nshell
!!  * \subpage rshell
!!  * \subpage cshell

!> \page nshell
!! `integer`(1:\ref npt)
!! **default:** `nshell(1:`\ref npt `) = 1`
!! * Number of shells

!> \page rshell
!! `real`(1:\ref nshell,1:\ref npt)
!! ** default:** `rshell(1,1:`\ref npt `) = [ (`\ref radat (ipt) `, ipt = 1,`\ref npt `) ]`
!! * Radius of the shells.

!> \page cshell
!! `real`(1:\ref nshell ,1:\ref npt )
!! **default:** `cshell(1,1:`\ref npt `) = One`


subroutine ScatIntens(nbin, q, sfpar)

   use MolModule
#if !defined (_NOIEEE_)
   use, intrinsic :: IEEE_ARITHMETIC
#endif
   implicit none

   character(40), parameter :: txroutine ='ScatIntens'
   integer(4)   , parameter :: mnshell = 4

   integer(4), intent(in) :: nbin
   real(8),    intent(in) :: q(nbin)
   real(8),    intent(in) :: sfpar(nbin,nptpt)

   integer(4), allocatable :: nshell(:)
   real(8), allocatable :: rshell(:,:), cshell(:,:)

   real(8) :: fnorm(0:nbin,npt)              ! normalized field amplitude
   real(8) :: fazero(npt)                    ! field amplitude at k = 0
   real(8) :: si(nbin)                       ! scattering intensity
   real(8) :: ff(nbin)                       ! effective form factor
   real(8) :: sftot(nbin)                    ! total structure factor
   real(8) :: fac, arg, fnorm2aver
   integer(4) :: ipt, jpt, iptjpt, ishell, ibin

   namelist /nmlScatIntens/ nshell, rshell, cshell

    if (.not.allocated(nshell)) then
       allocate(nshell(npt), rshell(mnshell,npt), cshell(mnshell,npt))
       nshell = 0
       rshell = 0.0E+00
       cshell = 0.0E+00
    end if

    nshell(1:npt) = 1
    rshell(1,1:npt) = [ (radat(ipt), ipt = 1,npt) ]
    cshell(1,1:npt) = One

   if (master) call FileOpen(uin, fin , 'form/noread')

   rewind(uin)
   read(uin,nmlScatIntens)
   close(uin)

! ... check conditions

   do ipt = 1, npt
   if (nshell(ipt) > mnshell) call Stop(txroutine, 'nshell(ipt) > mnshell', uout)
      do ishell = 2, nshell(ipt)
         if (rshell(ishell,ipt) < rshell(ishell-1,ipt)) call Stop(txroutine, 'rshell is not in right order', uout)
      end do
   end do

! ... calculate normalized field amplitudes and field amplitudes at k = 0

   fnorm = Zero
   do ipt = 1, npt
      do ishell = 1, nshell(ipt)
         fac = cshell(ishell,ipt)
         fac = cshell(ishell,ipt)-sum(cshell(ishell+1:nshell(ipt),ipt))
         fac = fac*FourPiThird*rshell(ishell,ipt)**3
         do ibin = 1, nbin
            arg = q(ibin)*rshell(ishell,ipt)
            if (arg == Zero) then
               fnorm(ibin,ipt) = fnorm(ibin,ipt) + fac
            else
               fnorm(ibin,ipt) = fnorm(ibin,ipt) + fac*(Three*(sin(arg)-arg*cos(arg))/arg**3)
            end if
         end do
         fazero(ipt) = fnorm(0,ipt)
      end do
      if(fazero(ipt) .ne. 0.0d0) then
         fac = One/fazero(ipt)
      else
#if !defined (_NOIEEE_)
         fac = IEEE_VALUE(fac,IEEE_QUIET_NAN)
#else
         fac = huge(fac)
#endif
      end if


      fnorm(0:nbin,ipt) = fnorm(0:nbin,ipt)*fac
   end do

! ... calculate scattering intenstiy

   si = Zero
   iptjpt = 0
   do ipt = 1, npt
      do jpt = ipt, npt
         iptjpt = iptjpt+1
         if (ipt == jpt) then
            fac = real(nppt(ipt))/np
         else
            fac = Two*sqrt(real(nppt(ipt)*nppt(jpt)))/np
         end if
         si(1:nbin) = si(1:nbin) +                                     &
         fac*(fazero(ipt)*fnorm(1:nbin,ipt))*(fazero(jpt)*fnorm(1:nbin,jpt))*sfpar(1:nbin,iptjpt)
      end do
   end do

! ... calculate effective form factor

   ff = Zero
   do ipt = 1, npt
      fac = real(nppt(ipt))/np
      ff(1:nbin) = ff(1:nbin)+fac*(fazero(ipt)*fnorm(1:nbin,ipt))**2
   end do
   fnorm2aver = sum(nppt(1:npt)*fazero(1:npt)**2)/np
   ff = ff/fnorm2aver

! ... calculate effective total structure factor

   sftot = si/(fnorm2aver*ff)

   write(ulist,'(4(a,a))')  'wave vector', char(9), 'scatt. iten.', char(9), 'form factor', char(9), 'total structure factor'
   write(ulist,'(i5)') nbin
   write(ulist,'(4(g15.5,a))') (q(ibin),char(9),si(ibin),char(9),ff(ibin),char(9),sftot(ibin),char(9),ibin = 1,nbin)

end subroutine ScatIntens

!************************************************************************
!> \page static static.F90
!! **SFNoPBC**
!! *calculate partial structure factors*
!************************************************************************


!     no periodical boundary conditions; averaged over 3x(100) directions

!> \page txkscale
!! `character(3)`
!! **default:** `lin`
!! * ``'lin'``:  Linear k-scale
!! * ``'log'``:  Logarithmic k-scale

!> \page klow
!! `real`
!! **default:** if \ref sphrad \f$ \neq \f$ 0: `2Pi/(10*`\ref sphrad `)`, else `0.01`
!! * Lower k-vector (linear scale)

!> \page logklow
!! `real`
!! **default:** `-4`
!! * Lower k-vector (logarithmic scale)

!> \page logkupp
!! `real`
!! **default:** `1`
!! * Upper k-vector (logarithmic scale)

!> \page nmlSFnoPBC_nbin nbin
!! `integer`
!! **default:** `100`
!! * Number of bins used to sample the partial structure factors.

subroutine SFNoPBC(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SFNoPBC'
   character(80), parameter :: txheading ='structure factor'
   character(3),  save :: txkscale
   integer(4) ,   save :: nvar, nbin
   real(8),       save :: klow, logklow, logkupp, dlogk
   type(df_var), allocatable, save :: var(:)
   integer(4) :: ip, ipt, jpt, ivar, ibin, m
   real(8)    :: norm, kx, ky, kz, inc
   complex(8) :: eikr(0:mnbin_df,1:3,1:npt), eikr1(1:3), eikrtemp

   namelist /nmlSF/ txkscale, klow, logklow, logkupp, nbin

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      txkscale = 'lin'                              ! linear k scale
      if(sphrad .ne. 0.0d0) then          !prevent division by zero
         klow    = TwoPi/(10.0d0*sphrad)               ! lower k-vector (linear scale)
      else
         klow = 0.01
      end if

      logklow =-Four                                ! lower k-vector (logarithmic scale)
      logkupp =+One                                 ! upper k-vector (logarithmic scale)
      nbin = 100

      rewind(uin)
      read(uin,nmlSF)

      call lowercase(txkscale)

      if (.not.lbcsph .and. .not.lbccyl .and. .not.lbcell) &
         call Stop(txroutine, '.not.lbcsph .and. not.lbccyl .and. .not.lbcell', uout)
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

      dlogk = (logkupp-logklow)/nbin

   case (iWriteInput)

! ... set nvar as well as allocate memory

      nvar = nptpt
      allocate(var(nvar))

! ... set label, min, max, and nbin

      do ipt = 1, npt
         do jpt = ipt, npt
            ivar = iptpt(ipt,jpt)
            var(ivar)%label = trim(txptpt(ivar))//' 100 directions'
            if (txkscale == 'lin') then
               var(ivar)%min = Half*klow
               var(ivar)%max = (nbin+Half)*klow
            else if (txkscale == 'log') then
               var(ivar)%min = logklow-Half*dlogk
               var(ivar)%max = logkupp+Half*dlogk
            else
               call Stop(txroutine, 'error in txkscale', uout)
            end if
            var(ivar)%nbin = nbin
         end do
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      eikr(0:nbin,1:3,1:npt) = cmplx(Zero,Zero)

! ... sample

      if (txkscale == 'lin') then
         do ip = 1, np
            ipt = iptpn(ip)
            kx = klow*ro(1,ip)
            ky = klow*ro(2,ip)
            kz = klow*ro(3,ip)
            eikr1(1) = cmplx(cos(kx),sin(kx))
            eikr1(2) = cmplx(cos(ky),sin(ky))
            eikr1(3) = cmplx(cos(kz),sin(kz))
            do m = 1, 3
               eikrtemp = cmplx(One,Zero)
               do ibin = 0, nbin-1
                  eikrtemp = eikrtemp*eikr1(m)
                  eikr(ibin,m,ipt) = eikr(ibin,m,ipt)+eikrtemp
               end do
            end do
         end do
      else if (txkscale == 'log') then
         do ip = 1, np
            ipt = iptpn(ip)
            do ibin = 0, nbin-1
               inc = 10.d0**(logklow+ibin*dlogk)
               kx = inc*ro(1,ip)
               ky = inc*ro(2,ip)
               kz = inc*ro(3,ip)
               eikr(ibin,1,ipt) = eikr(ibin,1,ipt)+cmplx(cos(kx),sin(kx))
               eikr(ibin,2,ipt) = eikr(ibin,2,ipt)+cmplx(cos(ky),sin(ky))
               eikr(ibin,3,ipt) = eikr(ibin,3,ipt)+cmplx(cos(kz),sin(kz))
            end do
         end do
      end if

      do ipt = 1, npt
         do jpt = ipt, npt
            ivar = iptpt(ipt,jpt)
            do ibin = 0, nbin-1
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ &
               sum(real(eikr(ibin,1:3,ipt))*real(eikr(ibin,1:3,jpt))+aimag(eikr(ibin,1:3,ipt))*aimag(eikr(ibin,1:3,jpt)))
            end do
         end do
      end do

   case (iAfterMacrostep)

! ... normalizing

      do ipt = 1, npt
         do jpt = ipt, npt
            ivar = iptpt(ipt,jpt)
            norm = One/(Three*sqrt(real(nppt(ipt)*nppt(jpt))))
            var(ivar)%avs2(-1:var(ivar)%nbin) = var(ivar)%avs2(-1:var(ivar)%nbin)*norm
         end do
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,i5,5x,a,i5,2x,a)') 'number of grid points          = ', nbin, '(',mnbin_df,')'
      if (txkscale == 'lin') then
         write(uout,'(a,t35,f8.2)')         'lower wavevector limit         = ', var(1)%min+0.5*klow
         write(uout,'(a,t35,f8.2)')         'upper wavevector limit         = ', var(1)%max-0.5*klow
      else if (txkscale == 'log') then
         write(uout,'(a,t35,f8.2)')         'log lower wavevector limit     = ', var(1)%min+0.5*dlogk
         write(uout,'(a,t35,f8.2)')         'log upper wavevector limit     = ', var(1)%max-0.5*dlogk
      end if

      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine SFNoPBC

!************************************************************************
!> \page static static.F90
!! **AngDF**
!! *calculate angular distribution functions*
!************************************************************************


!     p(a) =< delta(a-a') > where a is the angle cosine given (for particles within [rmin:rmax]) by
!
!     type  label  quantity
!     ----  -----  --------
!     1     z'*z'  cos( z'(ip)*z'(jp)     )
!     2     z'*r   cos( z'(ip)*r(ijp)     )
!     3     r*z'   cos( r(ijp)*z'(jp)     )
!     4     oh*r   cos( oh(ip)*r(ijp)     )   ip has to be water
!     5     r*oh   cos( r(ijp)*oh(jp)     )   jp has to be water
!     6     oh.o   cos( oh(ip).o(jp)      )   ip and jp have to be water
!     7     o.o.o  cos( o(jp).o(ip).o(kp) )
!     8     m*m    cos( m(ia)*m(ja)       )   m(ia) is the dipole vector of atom ia = ianpn(ip) + iashift
!     9     m*r    cos( m(ia)*r(ijp)      )                         (the value of iashift is currently hard coded)
!     10    r*m    cos( r(ijp)*m(ja)      )

!  and where particle ip belongs to group igr and jp as well as kp to jgr

!> \page nmlAngDF
!! The namelist  \ref nmlAngDF contains variables that control the calculation of angular distribution functions. Any combination of the
!! types of distribution functions listed below may be selected through vtype\%l.
!!
!!    | type |  label          |  quantity                                                                                      |
!!    | ---- | --------------  |------------------------------------------------------------------------------------------------|
!!    | 1    | z'*z'           |cos( z'(ip)*z'(jp)     )                                                                        |
!!    | 2    | z'*r            |cos( z'(ip)*r(ijp)     )                                                                        |
!!    | 3    | r*z'            |cos( r(ijp)*z'(jp)     )                                                                        |
!!    | 4    | oh*r            |cos( oh(ip)*r(ijp)     )   ip has to be water                                                   |
!!    | 5    | r*oh            |cos( r(ijp)*oh(jp)     )   jp has to be water                                                   |
!!    | 6    | oh.o            |cos( oh(ip).o(jp)      )   ip and jp have to be water                                           |
!!    | 7    | \format{o.o.o}  |cos( o(jp).o(ip).o(kp) )                                                                        |
!!    | 8    | m*m             |cos( m(ia)*m(ja)       )   m(ia) is the dipole vector of atom ia = ianpn(ip) + iashift          |
!!    | 9    | m*r             |cos( m(ia)*r(ijp)      )                         (the value of iashift is currently hard coded) |
!!    | 10   | r*m             |cos( r(ijp)*m(ja)      )                                                                        |
!!
!! h(ip) denotes the z'-axis of particle ip, r(ijp) the normalized vector between particle ip and jp defined as r(jp)-r(ip), oh(ip)
!! the direction of a oh-bond in water, oh(ip).o(jp) the largest angle of the four possible formed by the oh(ip) direction and the
!! h(ip)-o(jp) vector (hydrogen bond angle), and o(jp).o(ip).o(jp) the angle formed by the location of three particles. The sampling
!! of types 4-6 assumes that the first three sites of water are oxygen, hydrogen, and hydrogen. Generally ip refers to reference
!! particles and jp to field particles. Only field particles within a separation of rmax from a reference particle are considered.
!! * Variables:
!!  * \subpage nmlAngDF_vtype
!!  * \subpage nmlAngDF_rmin
!!  * \subpage nmlAngDF_rmax

!> \page nmlAngDF_vtype vtype
!! `static2D_var(logical, 2*real, 2*real, 2*integer, logical, character, real)`(1:10)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /(-1.0, 0.0)/ Max: /(1.0, 2Pi)/

!> \page nmlAngDF_rmin rmin
!! `real`
!! **default:** `0.0`
!! * Lower distance for particle separation considered.

!> \page nmlAngDF_rmax rmax
!! `real`
!! **default:** `10.0`
!! * Upper distance for particle separation considered.
subroutine AngDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AngDF'
   character(80), parameter :: txheading ='angular distribution functions'
   integer(4)   , parameter :: ntype = 10
   type(static1D_var),         save :: vtype(ntype)
   real(8),                    save :: rmin
   real(8),                    save :: rmax
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4) :: ivar, ibin, itype, m
   integer(4) :: ip, jp, kp, ia, ialow, ja, jalow, kl, igr, jgr, igrjgr, jjgr, id, jd
   real(8)    :: dx, dy, dz, r2, ac, actemp, dxoh, dyoh, dzoh, dkx, dky, dkz, dlx, dly, dlz
   real(8)    :: dxx, dyy, dzz, angcos, dmx, dmy, dmz
   integer(4), save :: klstep(4) = [ 1,1,0,0 ]
   integer(4), save :: kstep(4)  = [ 1,2,0,0 ]
   integer(4), save :: lstep(4)  = [ 0,0,1,2 ]
   integer(4), save :: iashift = 0

   namelist /nmlAngDF/ vtype, rmin, rmax

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = -One
      vtype%max = +One
      vtype%nbin = 100
      rmin  = Zero
      rmax  = 3.5d0

      rewind(uin)
      read(uin,nmlAngDF)
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, ' vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['z''*z''','z''*r ','r*z'' ','oh*r ','r*oh ','oh.o ','o.o.o','m*m  ','m*r  ','r*m  ']
      vtype%nvar = ngrgr

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(ngrgr,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do igrjgr = 1, ngrgr
               ivar = ivar+1
               ipnt(igrjgr,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txgrgr(igrjgr)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = ipmyid(1), ipmyid(2)
         igr = igrpn(ip,1)
         if (igr <= 0) cycle
         ia = ianpn(ip)
         id = ia+iashift
         do jp = 1, np
            if (ip == jp) cycle
            jgr = igrpn(jp,2)
            if (jgr <= 0) cycle
            igrjgr = igrgr(igr,jgr)
            ja = ianpn(jp)
            jd = ja+iashift
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBC(dx,dy,dz)
            if ((id <= na) .and. (jd <= na)) then
               dmx = r(1,id)-r(1,jd)
               dmy = r(2,id)-r(2,jd)
               dmz = r(3,id)-r(3,jd)
               call PBC(dmx,dmy,dmz)
            end if
            r2 = dx**2+dy**2+dz**2
            if (r2 > rmax**2) cycle
            if (rmin**2 > r2) cycle
! wall      if (r2 < 306) cycle

! ... sample type 1

            itype = 1
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 2

            itype = 2
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),-dx,-dy,-dz)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 3

            itype = 3
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(+dx,+dy,+dz,ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 4

            itype = 4
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ialow = ianpn(ip)
               dxoh = r(1,ialow+1)-r(1,ialow)
               dyoh = r(2,ialow+1)-r(2,ialow)
               dzoh = r(3,ialow+1)-r(3,ialow)
               ac = angcos(dxoh,dyoh,dzoh,-dx,-dy,-dz)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               dxoh = r(1,ialow+2)-r(1,ialow)
               dyoh = r(2,ialow+2)-r(2,ialow)
               dzoh = r(3,ialow+2)-r(3,ialow)
               ac = angcos(dxoh,dyoh,dzoh,-dx,-dy,-dz)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 5

            itype = 5
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               jalow = ianpn(jp)
               dxoh = r(1,jalow+1)-r(1,jalow)
               dyoh = r(2,jalow+1)-r(2,jalow)
               dzoh = r(3,jalow+1)-r(3,jalow)
               ac = angcos(+dx,+dy,+dz,dxoh,dyoh,dzoh)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               dxoh = r(1,jalow+2)-r(1,jalow)
               dyoh = r(2,jalow+2)-r(2,jalow)
               dzoh = r(3,jalow+2)-r(3,jalow)
               ac = angcos(+dx,+dy,+dz,dxoh,dyoh,dzoh)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 6

            itype = 6
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac =+One
               do m = 1, 4
                  kl = klstep(m)*(ia+kstep(m))+(1-klstep(m))*(ja+lstep(m))
                  dkx = r(1,ia)-r(1,kl)
                  dky = r(2,ia)-r(2,kl)
                  dkz = r(3,ia)-r(3,kl)
                  call PBC(dkx,dky,dkz)
                  dlx = r(1,ja)-r(1,kl)
                  dly = r(2,ja)-r(2,kl)
                  dlz = r(3,ja)-r(3,kl)
                  call PBC(dlx,dly,dlz)
                  actemp = angcos(dkx,dky,dkz,dlx,dly,dlz)
                  ac = min(ac,actemp)
               end do
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if

! ... sample type 7

            itype = 7
            if (vtype(itype)%l) then
              ivar = ipnt(igrjgr,itype)
               do kp = jp+1, np
                  if (kp == ip) cycle
                  jjgr = igrpn(kp,2)
                  if (jjgr/= jgr) cycle
                  dxx = ro(1,ip)-ro(1,kp)
                  dyy = ro(2,ip)-ro(2,kp)
                  dzz = ro(3,ip)-ro(3,kp)
                  call PBC(dxx,dyy,dzz)
                  r2 = dxx**2+dyy**2+dzz**2
                  if (r2 > rmax**2) cycle
                  ac = angcos(-dx,-dy,-dz,-dxx,-dyy,-dzz)
                  ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
                 var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
              end do
           end if

! ... sample type 8

            if ((id <= na) .and. (jd <= na)) then
            itype = 8
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(dip(1,id),dip(2,id),dip(3,id),dip(1,jd),dip(2,jd),dip(3,jd))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
            end if

! ... sample type 9

            if ((id <= na) .and. (jd <= na)) then
            itype = 9
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(dip(1,id),dip(2,id),dip(3,id),-dmx,-dmy,-dmz)
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
            end if

! ... sample type 10

            if ((id <= na) .and. (jd <= na)) then
            itype = 10
            if (vtype(itype)%l) then
               ivar = ipnt(igrjgr,itype)
               ac = angcos(dip(1,id),dip(2,id),dip(3,id),dmx,dmy,dmz)                  ! id -> jd ? /Per
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
            end if

         end do
      end do

   case (iAfterMacrostep)

#if defined (_PAR_)
      do ivar = 1, nvar
          call par_allreduce_reals(var(ivar)%avs2(-1), vaux, var(ivar)%nbin+2)
      end do
#endif

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f8.2)') 'lower separation distance      = ', rmin
         write(uout,'(a,t35,f8.2)') 'upper separation distance      = ', rmax
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AngDF

!************************************************************************
!> \page static static.F90
!! **AngExtDF**
!! *calculate angular 2d distribution functions with respect to external frame*
!************************************************************************


!> \page nmlAngExtDF
!! The namelist  \ref nmlAngDF contains variables that control the calculation of 2D angular distribution functions with respect to an external frame.
!!
!!    | type | label  |  quantity                              |
!!    | ---- | -------|  --------------------------------------|
!!    | 1    | z'-axis|  cos(theta)-phi 2d distributin function|
!!
!! * Variables:
!! * \subpage nmlAngExtDF_vtype

!> \page nmlAngExtDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /7*(-1.0)/ Max: /7*0.0/

subroutine AngExtDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='AngExtDF'
   character(80), parameter :: txheading ='angular-external frame 2d distribution functions'
   integer(4)   , parameter :: ntype = 1
   type(static2D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df2D_var),allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4)    :: itype, ivar, ibin1, ibin2, ip, ipt
   real(8)       :: r1, theta, phi, ac

   namelist /nmlAngExtDF/ vtype

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype(1)%min(1:2) = [ -One, Zero ]
      vtype(1)%max(1:2) = [ One, Two*Pi ]
      vtype(1)%nbin(1:2) = [ 20, 20 ]

      rewind(uin)
      read(uin,nmlAngExtDF)

      if (maxval(vtype%nbin(1)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(1) > mnbin_df2D', uout)
      if (maxval(vtype%nbin(2)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(2) > mnbin_df2D', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['z''-axis']
      vtype%nvar = npt

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(npt,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         do ipt = 1, npt
            ivar = ivar+1
            ipnt(ipt,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)//' '//txpt(ipt)
            var(ivar)%min(1:2) = vtype(itype)%min(1:2)
            var(ivar)%max(1:2) = vtype(itype)%max(1:2)
            var(ivar)%nbin(1:2)= vtype(itype)%nbin(1:2)
         end do
      end do
      call DistFunc2DSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do itype = 1, ntype
         do ipt = 1, npt
            ivar = ipnt(ipt,itype)
            do ip = ipnpt(ipt), ipnpt(ipt)+(nppt(ipt)-1)
               call CarToSph('rad',ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),r1,theta,phi)
               ac = cos(theta)
               ibin1 = max(-1,min(floor(var(ivar)%bini(1)*(ac-var(ivar)%min(1))),int(var(ivar)%nbin(1))))
               ibin2 = max(-1,min(floor(var(ivar)%bini(2)*(phi-var(ivar)%min(2))),int(var(ivar)%nbin(2))))
               var(ivar)%avs2(ibin1,ibin2) = var(ivar)%avs2(ibin1,ibin2) + One
            end do
         end do
      end do

   case (iAfterMacrostep)

! ... normalization

      call DistFunc2DNorm(1, nvar, var)
      call DistFunc2DSample(iStage, nvar, var)

      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      do ivar = 1, nvar
         var(ivar)%avs1 = FourPi*var(ivar)%avs1   ! change normalization
         var(ivar)%avsd = FourPi*var(ivar)%avsd   ! change normalization
      end do
      call WriteHead(2, txheading, uout)
      write(uout,'(a                     )') 'first dimension (column) is cosine theta'
      write(uout,'(a                     )') 'second dimension (row) is phi'
      call DistFunc2DHead(nvar, var, uout)
      call DistFunc2DShow(1, txheading, nvar, var, uout)
      call DistFunc2DList(1, txheading, nvar, var, ulist)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AngExtDF

!************************************************************************
!> \page static static.F90
!! **OriDipDF**
!! *calculate orientation/dipole distribution function*
!************************************************************************


!> \page nmlOriDipDF
!! The namelist  \ref nmlOriDipDF contains variables that control the calculation of orientation/dipole distribution functions.
!!
!!    | type | label    | distribution functions               |
!!    | ---- | -------- | ------------------------------------ |
!!    | 1    | dir      | orientation df based on ori(1:3,1:3) |
!!    | 2    | dir_aver | symmetry-averaged orientation df     |
!!    | 3    | dip      | dipole df                            |
!!    | 4    | dip_tot  | total dipole df                      |
!!
!! * Variables:
!!  * \subpage nmlOriDipDF_vtype

!> \page nmlOriDipDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /4*(-1)/ Max: /4*1.0/

subroutine OriDipDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='OriDipDF'
   character(80), parameter :: txheading ='orientation/dipole distribution function'
   integer(4)   , parameter :: ntype = 4
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4)   , parameter :: ndir(1:ntype) = [13, 3, 3, 3]
   character(5),  save      :: txdir(ntype,max(ndir(1),ndir(2)))
   real(8)     ,  parameter :: dir(3,ndir(1)) = reshape( &
           [One,Zero,Zero,  Zero, One,Zero,  Zero,Zero, One, &
            One, One,Zero,   One,Zero, One,  Zero, One, One, -One, One,Zero,  -One,Zero, One,  Zero,-One, One, &
            One, One, One,  -One, One, One,   One,-One, One,  One, One,-One] , [3,ndir(1)] )
   integer(4)   , parameter :: ipntdir(1:ndir(1)) = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3]
   integer(4) :: itype, ivar, idir, ibin, ip
   real(8)    :: ac, angcos

   namelist /nmlOriDipDF/ vtype

   if (slave) return   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.true.
      vtype(1:2)%min = -1.0
      vtype(1:2)%max =  1.0
      vtype(1:2)%nbin = 100
      vtype(3:4)%min = -1.0
      vtype(3:4)%max =  1.0
      vtype(3:4)%nbin = 100

      rewind(uin)
      read(uin,nmlOriDipDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

      txdir(1,1:ndir(1)) = ['100 ','010 ','001 ','110 ','101 ','011 ','-110','-101','0-11', '111 ','-111','1-11','11-1']
      txdir(2,1:ndir(2)) = [ '100','110','111' ]
      txdir(3,1:3)       = [ 'x','y','z' ]
      txdir(4,1:3)       = [ 'x','y','z' ]

! ... set remaining elements of vtype

      vtype%label = ['dir:     ', 'dir_aver:', 'dip:     ', 'dip_tot: ']
      vtype%nvar = [ndir(1), ndir(2), 3, 3]

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(maxval(ndir(1:ntype)),ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do idir = 1, ndir(itype)
               ivar = ivar + 1
               ipnt(idir,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//trim(txdir(itype,idir))
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+np

! ... sample type 1

      itype = 1
      if (vtype(itype)%l) then
         do ip = ipmyid(1), ipmyid(2)
            do idir = 1, ndir(1)
               ivar = ipnt(idir,itype)
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip), dir(1,idir),dir(2,idir),dir(3,idir))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end do
         end do
      end if

! ... sample type 2

      itype = 2
      if (vtype(itype)%l) then
         do ip = ipmyid(1), ipmyid(2)
            do idir = 1, ndir(1)
               ivar = ipnt(ipntdir(idir),itype)
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip), dir(1,idir),dir(2,idir),dir(3,idir))
               ibin = max(-1,min(floor(var(ivar)%bini*(ac-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end do
         end do
      end if

! ... sample type 3

      itype = 3
      if (vtype(itype)%l) then
         do ip = ipmyid(1), ipmyid(2)
            do idir = 1, 3
               ivar = ipnt(idir,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(dip(idir,ip)-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end do
         end do
      end if

! ... sample type 4

      itype = 4
      if (vtype(itype)%l) then
         do idir = 1, 3
            ivar = ipnt(idir,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(sum(dip(idir,1:np))-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
         end do
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OriDipDF

!************************************************************************
!> \page static static.F90
!! **SphHarAver**
!! *calculate averages of unnormalized spherical harmonics*
!************************************************************************


subroutine SphHarAver(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SphHarAver'
   character(80), parameter :: txheading ='averages of unnormalized spherical harmonics'
   integer(4),                    save :: lmax               ! maximal value of l
   integer(4),                    save :: nvar
   type(scalar_var), allocatable, save :: var(:)

   integer(4) :: ivar, l, m, ip
   real(8) :: r1, theta, phi
   character(3) :: txl, txm
   complex(8) :: clm, cclm

   if (slave) return   ! master only

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

! ... set lmax (can be changed)

      lmax = 4

! ... set nvar as well as allocate memory

      nvar = (lmax+1)**2
      allocate(var(nvar))

! ... set label and norm

      ivar = 0
      do l = 0, lmax
         write(txl,'(i3)') l
         do m = 0, l
            write(txm,'(i3)') m
            ivar = ivar + 1
            var(ivar)%label = txl//txm//' r'
            if (m == 0) cycle
            ivar = ivar + 1
            var(ivar)%label = txl//txm//' i'
         end do
      end do

      var%norm = one/np
      var%norm = one

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var%value = Zero
      do ip = 1, np               ! loop over all particles
         call CarToSph('rad',ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),r1,theta,phi)
         ivar = 0
         do l = 0, lmax
            do m = 0, l
               clm = CCLM(l,m,theta,phi,0)
               ivar = ivar + 1
               var(ivar)%value = real(clm)
               if (m == 0) cycle
               ivar = ivar + 1
               var(ivar)%value = imag(clm)
            end do
         end do
         call ScalarSample(iStage, 1, nvar, var)
      end do

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var
      call ScalarNorm(iStage, 1, nvar, var, 1)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)

      call WriteHead(2, txheading, uout)
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine SphHarAver

!************************************************************************
!> \page static static.F90
!! **RadAngDF**
!! *calculate radial-angular 2d distribution functions*
!************************************************************************


!> \page nmlRadAngDF
!! The namelist  \ref nmlRadAngDF contains variables that control the calculation of radial-angular 2d distribution functions.
!!
!!    | type |  label        |  quantity                                              |
!!    | ---- | ------------- | ------------------------------------------------------ |
!!    | 1    | x'y'-z'-plane | two particle radial - angular r2d distributin function |
!!
!! * Variables:
!!  * \subpage nmlRadAngDF_vtype
!!  * \subpage txCoordSys
!!  * \subpage txAngle

!> \page nmlRadAngDF_vtype vtype
!! `static2D_var(logical, 2*real, 2*real, 2*integer, logical, character, real)`
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Default values depend on \ref txCoordSys

!> \page txCoordSys
!! `character(5)`
!! **default:** `'polar'`
!! * Coordinate system of distribution functions.
!! * ``'polar'``:  r and cos(theta)
!! * ``'cart'``:  sqrt(x**2+y**2) and z

!> \page txAngle
!! `character(6)`
!! **default:** `'theta1'`
!! * Choice of angle.
!! * ``'theta1'``:  theta 1
!! * ``'theta2'``:  tehta 2
!! * ``'psi'``:  psi


subroutine RadAngDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='RadAngDF'
   character(80), parameter :: txheading ='radial-angular 2d distribution functions'
   integer(4)   , parameter :: ntype = 1
   type(static2D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df2D_var),  allocatable, save :: var(:)          ! radial-angular 2d df
   type(df2D_var),  allocatable, save :: var2(:)         ! radial-angular 2d polarization df
   integer(4),      allocatable, save :: ipnt(:,:)
   character(5),                 save :: txCoordSys      ! 'polar' or 'cart'
   character(6),                 save :: txAngle         ! 'theta1', 'theta2', or 'psi'
   integer(4)    :: itype, ivar, ibin1, ibin2, ip, jp, ipt, jpt
   real(8)       :: dx, dy, dz, r1, r2, zz, xy, ac, an, angcos, dvol, darea, InvFlt

   namelist /nmlRadAngDF/ vtype, txCoordSys, txAngle

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype(1)%min(1:2) = [ Pi, Pi ]
      vtype(1)%max(1:2) = [ Pi, Pi ]
      txCoordSys = 'polar'
      txAngle = 'theta1'

      rewind(uin)
      read(uin,nmlRadAngDF)

      call LowerCase(txCoordSys)

! ... check condition

      if (maxval(vtype%nbin(1)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(1) > mnbin_df2D', uout)
      if (maxval(vtype%nbin(2)) > mnbin_df) call Stop(txroutine, 'vtype%nbin(2) > mnbin_df2D', uout)
      if (txCoordSys /= 'polar' .and. txCoordSys /= 'cart') call stop(txroutine, 'error in txCoordSys', uout)
      if (txAngle /= 'theta1' .and. txAngle /= 'angle2' .and. txAngle /= 'psi') call stop(txroutine, 'error in txAngle', uout)
      if (txCoordSys == 'cart' .and. txAngle == 'psi') call stop(txroutine, 'erroneous combindation of txCoordSys and txAngle', uout)

! ... set defult values of vmin and vmax

      if (vtype(1)%min(1) == Pi .and. vtype(1)%min(2) == Pi .and. vtype(1)%max(1) == Pi .and. vtype(1)%max(2) == Pi) then
         if (txCoordSys == 'polar') then
            vtype(1)%min(1:2) = [ 0.0d0, -1.d0 ]
            vtype(1)%max(1:2) = [ 10.0d0, 1.0d0 ]
            vtype(1)%nbin(1:2) = [ 20, 10 ]
         else if (txCoordSys == 'cart') then
            vtype(1)%min(1:2) = [ 0.0d0, -10.0d0 ]
            vtype(1)%max(1:2) = [ 10.0d0, 10.0d0 ]
            vtype(1)%nbin(1:2) = [ 40, 40 ]
         end if
      end if

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['x''y''-z''--plane']
      vtype%nvar = 1

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), var2(nvar), ipnt(1,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min(1:2) = vtype(itype)%min(1:2)
         var(ivar)%max(1:2) = vtype(itype)%max(1:2)
         var(ivar)%nbin(1:2)= vtype(itype)%nbin(1:2)
         var2(ivar)%label = var(ivar)%label
         var2(ivar)%min(1:2)  = var(ivar)%min(1:2)
         var2(ivar)%max(1:2)  = var(ivar)%max(1:2)
         var2(ivar)%nbin(1:2) = var(ivar)%nbin(1:2)
      end do
      call DistFunc2DSample(iStage, nvar, var)
      call DistFunc2DSample(iStage, nvar, var2)

   case (iBeforeSimulation)

      call DistFunc2DSample(iStage, nvar, var)
      call DistFunc2DSample(iStage, nvar, var2)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var2

   case (iBeforeMacrostep)

      call DistFunc2DSample(iStage, nvar, var)
      call DistFunc2DSample(iStage, nvar, var2)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1
      var2%nsamp2 = var2%nsamp2 + 1

      do itype = 1, ntype
      ivar = ipnt(1,itype)
      do ip = ipmyid(1), ipmyid(2)              ! loop over all particles
         ipt = iptpn(ip)
         do jp = 1, np
            if (jp == ip) cycle
            jpt = iptpn(jp)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > (2*var(1)%max(1))**2) cycle
            r1 = sqrt(r2)
            if (txAngle == 'theta1') then
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),-dx,-dy,-dz)
            else if (txAngle == 'theta2') then
               ac = angcos(dx, dy, dz, ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
            else if (txAngle == 'psi') then
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
            end if
            if (txCoordSys == 'polar') then
               ibin1 = max(-1,min(floor(var(ivar)%bini(1)*(r1-var(ivar)%min(1))),int(var(ivar)%nbin(1))))
               if (txAngle == 'theta1' .or. txAngle == 'theta2') then
                  ibin2 = max(-1,min(floor(var(ivar)%bini(2)*(ac-var(ivar)%min(2))),int(var(ivar)%nbin(2))))
               else if (txAngle == 'psi') then
                  ibin2 = max(-1,min(floor(var(ivar)%bini(2)*(ac-var(ivar)%min(2))),int(var(ivar)%nbin(2))))
               end if
            else if (txCoordSys == 'cart') then
               an = acos(ac)
               xy = r1*sin(an)       ! distance separation perpendcular to the molecular axis
               zz = r1*cos(an)       ! distance separation in the direction of the molecular z-axis
               ibin1 = max(-1,min(floor(var(ivar)%bini(1)*(xy-var(ivar)%min(1))),int(var(ivar)%nbin(1))))
               ibin2 = max(-1,min(floor(var(ivar)%bini(2)*(zz-var(ivar)%min(2))),int(var(ivar)%nbin(2))))
            end if
            var(ivar)%avs2(ibin1,ibin2) = var(ivar)%avs2(ibin1,ibin2) + One
            var2(ivar)%avs2(ibin1,ibin2) = var2(ivar)%avs2(ibin1,ibin2)+sum(ori(1:3,3,ip)*ori(1:3,3,jp))
         end do
      end do
      end do

   case (iAfterMacrostep)

! ... reduce var%avs2 to master

#if defined (_PAR_)
      if (ltime) call CpuAdd('start', 'comm', 1, uout)
      do ivar = 1, nvar
         call par_reduce_reals(var(ivar)%avs2(-1,-1), vaux, (mnbin_df2D+2)**2)
         call par_reduce_reals(var2(ivar)%avs2(-1,-1), vaux, (mnbin_df2D+2)**2)
      end do
      if (ltime) call CpuAdd('stop', 'comm', 1, uout)
#endif

! ... normalization

      if (master) then
         if (txCoordSys == 'polar') then
            do ivar = 1, nvar
               var(ivar)%norm = var(ivar)%nbin(2)*vol/(FourPiThird*np*(np-1))
               do ibin2 = -1, var(ivar)%nbin(2)
                  do ibin1 = 0, var(ivar)%nbin(1)
                     var2(ivar)%avs2(ibin1,ibin2) = var2(ivar)%avs2(ibin1,ibin2)*InvFlt(var(ivar)%avs2(ibin1,ibin2))*var(ivar)%nsamp2
                     var(ivar)%avs2(ibin1,ibin2) = var(ivar)%avs2(ibin1,ibin2)*var(ivar)%norm/dvol(ibin1,var(ivar)%min(1),var(ivar)%bin(1))
                  end do
               end do
            end do
         else if (txCoordSys == 'cart') then
            do ivar = 1, nvar
               var(ivar)%norm = var(ivar)%bini(2)*vol/(Pi*np*(np-1))
               do ibin2 = -1, var(ivar)%nbin(2)
                  do ibin1 = 0, var(ivar)%nbin(1)
                     var2(ivar)%avs2(ibin1,ibin2) = var2(ivar)%avs2(ibin1,ibin2)/var(ivar)%avs2(ibin1,ibin2)*var2(ivar)%nsamp2
                     var(ivar)%avs2(ibin1,ibin2) = var(ivar)%avs2(ibin1,ibin2)*var(ivar)%norm/darea(ibin1,var(ivar)%min(1),var(ivar)%bin(1))
                  end do
               end do
            end do
         end if
         call DistFunc2DSample(iStage, nvar, var)
         call DistFunc2DSample(iStage, nvar, var2)
      end if

      if (lsim .and. master) write(ucnf) var
      if (lsim .and. master) write(ucnf) var2

   case (iAfterSimulation)

      if (master) then
         call DistFunc2DSample(iStage, nvar, var)
         call DistFunc2DSample(iStage, nvar, var2)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,a               )') 'coord system of indep vars     = ', txCoordSys
         write(uout,'(a,t35,a               )') 'angle                          = ', txAngle
         if (txCoordSys == 'polar') then
            write(uout,'(a                     )') 'first dimension (row) is the radial distrance'
            write(uout,'(a                     )') 'second dimension (row) is the cosine of the angle'
         else if (txCoordSys == 'cart') then
            write(uout,'(a                     )') 'first dimension (column) is the distance in the x''-y''--plane'
            write(uout,'(a                     )') 'second dimension (row) is the distance along the z''--axis'
         end if
         call DistFunc2DHead(nvar, var, uout)
         call DistFunc2DShow(1, txheading, nvar, var, uout)
         call DistFunc2DList(1, txheading, nvar, var, ulist)
         call DistFunc2DHead(nvar, var2, uout)
         call DistFunc2DShow(1, trim(txheading)//' polarization density', nvar, var2, uout)
         call DistFunc2DList(1, trim(txheading)//' polarization density', nvar, var2, ulist)
      end if

      deallocate(var, var2)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine RadAngDF

!************************************************************************
!> \page static static.F90
!! **Kirkwoodgk**
!! *calculate Kirkwood gk-factor*
!************************************************************************


!> \page nmlKirkwoodgk
!! The namelist  \ref nmlKirkwoodgk contains variables that control the calculation of running orientation averages. Any combination of
!! the types of distribution functions listed below may be selected through vtype\%l.
!!
!!    | type | label | quantity           |
!!    | ---- | ----- | ------------------ |
!!    | 1    | h*h   | sum_jp h(ip)*h(jp) |
!!
!! * Variables:
!!  * \subpage nmlKirkwoodgk_vtype
!!  * \subpage nmlKirkwoodgk_rmax

!> \page nmlKirkwoodgk_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization (by radius**3/2), title, and number of variable of vtype.
!! * Min: /0.0/ Max: /10.0/

!> \page nmlKirkwoodgk_rmax rmax
!! `real`
!! **default:** `10.0`
!! * Upper distance for particle separation considered.

subroutine Kirkwoodgk(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Kirkwoodgk'
   character(80), parameter :: txheading ='Kirkwood gk-factor'
   integer(4)   , parameter :: ntype = 1
   type(static1D_var),         save :: vtype(ntype)
   real(8),                    save :: rmax
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4) :: itype, ivar, ibin, ip, ipt, jp, jpt, iptjpt
   real(8)    :: dx, dy, dz, r2, ac, angcos

   namelist /nmlKirkwoodgk/ vtype, rmax

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (npt /= 1) call Stop(txroutine, 'npt /=1', uout) ! not yet fully adapted for npt /= 1

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.false.
      vtype(1)%min = Zero
      vtype(1)%max = 10.0d0
      vtype%nbin   = 100
      rmax  = 10.0

      rewind(uin)
      read(uin,nmlKirkwoodgk)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['z''*z''']
      vtype%nvar = nptpt

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nptpt,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do iptjpt = 1, nptpt
               ivar = ivar+1
               ipnt(iptjpt,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txptpt(iptjpt)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = ipmyid(1), ipmyid(2)
         ipt = iptpn(ip)
         do jp = 1, np
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rmax**2) cycle

! ... sample type 1

            itype = 1
            if (vtype(itype)%l) then
               ivar = ipnt(iptjpt,itype)
               ac = angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),ori(1,3,jp),ori(2,3,jp),ori(3,3,jp))
               ibin = max(-1,min(floor(var(ivar)%bini*(sqrt(r2)-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+ac
            end if

         end do
      end do

   case (iAfterMacrostep)

#if defined (_PAR_)
      ivar = 1
      call par_allreduce_reals(var(ivar)%avs2(-1), vaux, var(ivar)%nbin+2)
#endif

      do ipt = 1, npt
         do jpt = ipt, npt
            iptjpt = iptpt(ipt,jpt)
            if (vtype(1)%l) then
               ivar = ipnt(iptjpt,1)
               var(ivar)%norm = InvInt(nppt(ipt))
            end if
         end do
      end do
      do ivar = 1, nvar
         ibin = -1
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*var(ivar)%norm
         do ibin = 0, var(ivar)%nbin
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1) + var(ivar)%avs2(ibin)*var(ivar)%norm
         end do
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,f8.2)') 'upper separation distance      = ', rmax
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine Kirkwoodgk

!************************************************************************
!> \page static static.F90
!! **OriPolDF**
!! *calculate orientation polarization distribution functions*
!************************************************************************



!> \page nmlOriPolDF
!! The namelist  \ref nmlOriPolDF contains variables that control the calculation of orientation polarization distribution functions of
!! molecular dipole moments. Any combination of the types of distribution functions listed below may be selected through vtype\%l.
!!
!!    | type | label | distribution functions                                                                 |
!!    | ---- | ----- | -------------------------------------------------------------------------------------- |
!!    | 1    | tot   | orientation polarization, total                                                        |
!!    | 2    | par   | orientation polarization, parallel to central dipole (related to Kirkwood's Gk factor) |
!!    | 3    | per   | orientation polarization, perpendicular to central dipole                              |
!!    | 4    | pevpa | orientation polarization, perpendicular versus parallel                                |
!!
!! * The analysis is made for a set of different radii
!! * The orientation polarization is normalized by (dipole moment)**2 or radius**(3/2)*(dipole moment)**2
!!
!! * Variables:
!!  * \subpage nmlOriPolDF_vtype
!!  * \subpage lnorm
!!  * \subpage nmlOriPolDF_nrad
!!  * \subpage nmlOriPolDF_radius

!> \page nmlOriPolDF_vtype vtype
!! `static2D_var(logical, real, real, integer, logical, character, real)`(1:4)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization (by radius**3/2), title, and number of variable of vtype.
!! * Min: /0.0, -10.0, 0.0, 0.0/ Max: /4*(10.0)/

!> \page lnorm
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Invoking normalization.
!! * `.false.`:  No normalization invoked.

!> \page nmlOriPolDF_nrad nrad
!! `integer`
!! **default:** `2`
!! * Number of radii to consider

!> \page nmlOriPolDF_radius radius
!> `real`(1:nrad)
!! **default:** `20.0`
!! * Radii to consider

subroutine OriPolDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='OriPolDF'
   character(80), parameter :: txheading ='orientation polarization distribution functions'
   integer(4)   , parameter :: mnrad = 8
   integer(4)   , parameter :: ntype = 4
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   logical,                    save :: lnorm
   integer(4),                 save :: nrad
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),                    save :: radius(mnrad), radius2(mnrad)
   integer(4)    :: itype, ivar, irad, ibin, ivar2, ip, jp
   real(8)       :: dx, dy, dz, r2, angcos, pol(mnrad,3), InvFlt
   real(8)       :: value(3)                ! total, parallel, and perpendicular polarization
   character(5)  :: txrad

   namelist /nmlOriPolDF/ vtype, lnorm, nrad, radius

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = [ 0.0d0, -10.0d0,  0.0d0,  0.0d0]
      vtype%max = [10.0d0,  10.0d0, 10.0d0, 10.0d0]
      vtype%nbin = 100
      lnorm =.false.
      nrad = 2
      radius = 20.0

      rewind(uin)
      read(uin,nmlOriPolDF)

      radius2(1:nrad) = radius(1:nrad)**2

! ... check condition

      if (vtype(4)%l) then                           ! normalization dependence
         if (.not.vtype(2)%l) call Stop(txroutine, 'vtype(4)%l .and. .not.vtype(2)%l', uout)
         vtype(4)%min = vtype(2)%min
         vtype(4)%max = vtype(2)%max
      end if
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)
      if (nrad > mnrad) call Stop(txroutine, 'nrad > mnrad', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['total','par  ','per  ','pevpa']
      vtype%nvar = nrad

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nrad,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
         do irad  = 1, nrad
            ivar = ivar+1
            ipnt(irad,itype) = ivar
            write(txrad,'(f5.1)') radius(irad)
            var(ivar)%label = trim(vtype(itype)%label)//' '//txrad
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = ipmyid(1), ipmyid(2)              ! loop over all particles

         pol = Zero               ! calculate the polairization within distance radius(irad) from ip
         do jp = 1, np
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            do irad = 1, nrad
               if (r2 < radius2(irad)) pol(irad,1:3) = pol(irad,1:3) + ori(1:3,3,jp)
            end do
         end do

         do irad = 1, nrad
            value(1) = sqrt(pol(irad,1)**2 + pol(irad,2)**2 + pol(irad,3)**2)
            value(2) = value(1)*angcos(ori(1,3,ip),ori(2,3,ip),ori(3,3,ip),pol(irad,1),pol(irad,2),pol(irad,3))
            value(3) = sqrt(max(Zero,value(1)**2-value(2)**2))
            if (lnorm) value(1:3) = value(1:3)/radius(irad)**1.5d0        ! normalization with respect to r**3/2

! ... sample type 1, 2, and 3

            do itype = 1, ntype-1
               if (vtype(itype)%l) then
               ivar = ipnt(irad,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(value(itype)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if
            end do

! ... sample type 4

            itype = 4
            if (vtype(itype)%l) then
               ivar = ipnt(irad,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(value(2)-var(ivar)%min)),var(ivar)%nbin))
!                write(*,'(a,i3,3f12.2)') 'ibin, value(2), value(3)',ibin, value(2), value(3)
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+value(3)
            end if
         end do

      end do

   case (iAfterMacrostep)

! ... reduce var%avs2 to master

#if defined (_PAR_)
      if (ltime) call CpuAdd('start', 'comm', 1, uout)
      do ivar = 1, nvar
         call par_reduce_reals(var(ivar)%avs2(-1), vaux, mnbin_df+2)
      end do
      if (ltime) call CpuAdd('stop', 'comm', 1, uout)
#endif

! ... normalization

      if (master) then
         if (vtype(4)%l) then
            do irad = 1, nrad
               ivar = ipnt(irad,4)
               ivar2 = ipnt(irad,2)
               do ibin = 1, var(ivar)%nbin
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*InvFlt(var(ivar2)%avs2(ibin))*var(ivar)%nsamp2
               end do
            end do
         end if
         do itype = 1, 3                            ! area normalized to unity
            if (vtype(itype)%l) call DistFuncNorm(ipnt(1,itype), ipnt(nrad,itype), var)
         end do
         call DistFuncSample(iStage, nvar, var)
      end if

   case (iAfterSimulation)

      if (master) then
         call DistFuncSample(iStage, nvar, var)
         call WriteHead(2, txheading, uout)
         write(uout,'(a,t35,l8)')              'normalization by radius**(3/2) = ', lnorm
         write(uout,'(a,t35,i8  )')            'number of radiii               = ', nrad
         write(uout,'(a,t35,9f8.2)')           'radii                          = ', radius(1:nrad)
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
         call DistFuncAverValue(3*nrad, var, uout)
      end if

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OriPolDF

!************************************************************************
!> \page static static.F90
!! **NNHB**
!! *calculate number of nearest neighbours and number of hydrogen bonds*
!************************************************************************


!     local variables

!     thnn             distance thresholds
!     thhb             energy thresholds
!     inns2(n,ith,igr) number of particles belonging to group igr having n nearest neighbours of type jpt
!                      within the distance given by thnn(ith)
!     ihbs2(n,ith,igr) number of particles belonging to group igr having n nearest neighbours of type jpt
!                      within rmax interacting with at least thhb(ith)

!> \page nmlNNHB
!! The namelist  \ref nmlNNHB contains variables that control the calculation of no of nearest neighbours and no of hydrogen bonds. Field particles within \ref thnn are considered as nearest neighbours to a reference particle. Field particles within rmax and with an interaction energy below a threshold value are considered as hydrogen bonded to a reference particle.
!! * Variables:
!!  * \subpage nthnn
!!  * \subpage thnn
!!  * \subpage nthhb
!!  * \subpage thhb
!!  * \subpage nnnhb
!!  * \subpage nmlNNHB_rmax

!> \page nthnn
!! `integer`
!! **default:** `2`
!! * Number of separations considered.

!> \page thnn
!! `real`(1:\ref nthnn)
!! **default:** [  3.3,  3.5]
!! * Distances considered. The largest \ref thnn should not be larger than rmax given below.

!> \page nthhb
!! `integer`
!! **default:** `2`
!! * Number of energy thresholds.

!> \page thhb
!! `real`(1:\ref nthhb)
!! **default:** [-10.0,-16.0]
!! * Energy thresholds.

!> \page nnnhb
!! `integer`
!! **default:** `10`
!! * Maximal number of nearest neighbours and hydrogen bonds considered in the analysis.

!> \page nmlNNHB_rmax rmax
!! `real`
!! **default:** `3.5`
!! * Upper separation distance for hydrogen bonded particles.
subroutine NNHB(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NNHB'
   character(80), parameter :: txheading ='number of nearest neighbours and hydrogen bonds'
   integer(4)   , parameter :: mnth = 2
   integer(4), save :: nnn(mnth), nthnn, nhb(mnth), nthhb, nnnhb
   real(8),    save :: thnn(mnth), thhb(mnth)
   integer(4), allocatable, save :: inns2(:,:,:), inns1(:,:,:), ihbs2(:,:,:), ihbs1(:,:,:)
   real(8),    allocatable, save :: nnav(:,:), nnsd(:,:), hbav(:,:), hbsd(:,:), nnnorm(:,:), hbnorm(:,:)
   integer(4), save :: nsamp1
   real(8),    save :: rmax
   integer(4) :: ip, jp, igr, jgr, ith, inn, ihb
   real(8)    :: dx, dy, dz, r2, uuu, fdum(3), aver, term, norm, norm1

   namelist /nmlNNHB/ nthnn, thnn, nthhb, thhb, nnnhb, rmax

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      nthnn= 2
      nthhb= 2
      thnn =[  3.3,  3.5]
      thhb =[-10.0,-16.0]
      nnnhb = 10
      rmax = 3.5

      rewind(uin)
      read(uin,nmlNNHB)

      if (nthnn > mnth) call Stop(txroutine, 'nthnn > mnth', uout)
      if (nthhb > mnth) call Stop(txroutine, 'nthhb > mnth', uout)

   case (iWriteInput)

! ... allocate memory

      allocate(inns2(0:nnnhb,nthnn,ngr(1)), inns1(0:nnnhb,nthnn,ngr(1)), ihbs2(0:nnnhb,nthnn,ngr(1)), ihbs1(0:nnnhb,nthnn,ngr(1)))
      inns2 = 0
      inns1 = 0
      ihbs2 = 0
      ihbs1 = 0
      allocate(nnav(nthnn,ngr(1)), nnsd(nthnn,ngr(1)), hbav(nthnn,ngr(1)), hbsd(nthnn,ngr(1)), nnnorm(nthnn,ngr(1)), hbnorm(nthnn,ngr(1)))
      nnav = 0.0E+00
      nnsd = 0.0E+00
      hbav = 0.0E+00
      hbsd = 0.0E+00
      nnnorm = 0.0E+00
      hbnorm = 0.0E+00

! ... check that only one field group is selected

      if (ngr(2) /= 1) call Stop(txroutine, 'ngr(2) /= 1', uout)

      thnn(1:nthnn) = thnn(1:nthnn)**2

   case (iBeforeSimulation)

      nsamp1 = 0
      inns1 = 0
      ihbs1 = 0
      nnav = Zero
      nnsd = Zero
      hbav = Zero
      hbsd = Zero
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) nsamp1, inns1, ihbs1, nnav, nnsd, hbav, hbsd

   case (iBeforeMacrostep)

      inns2 = 0
      ihbs2 = 0

   case (iSimulationStep)

      do ip = 1, np
         igr = igrpn(ip,1)
         if (igr <= 0) cycle
         nnn = 0
         nhb = 0
         do jp = 1, np
            if (ip == jp) cycle
            jgr = igrpn(jp,2)
            if (jgr <= 0) cycle
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < rmax**2) then
               call UTwoBodyPair(ip, jp, uuu, fdum)
               do ith = 1, nthnn
                  if (r2 < thnn(ith)) nnn(ith) = nnn(ith) + 1
                  if (uuu < thhb(ith)) nhb(ith) = nhb(ith) + 1
               end do
            end if
         end do
         do ith = 1, nthnn
            inn = min(nnnhb,nnn(ith))
            inns2(inn,ith,igr) = inns2(inn,ith,igr) + 1
            ihb = min(nnnhb,nhb(ith))
            ihbs2(ihb,ith,igr) = ihbs2(ihb,ith,igr) + 1
         end do
      end do

   case (iAfterMacrostep)

      nsamp1 = nsamp1 + 1
      do igr = 1, ngr(1)
         do ith = 1, nthnn
            inns1(0:nnnhb,ith,igr) = inns1(0:nnnhb,ith,igr)+inns2(0:nnnhb,ith,igr)
            aver = sum( [ (ihb, ihb = 0,nnnhb) ]*inns2(0:nnnhb,ith,igr))
            term = aver*InvInt(sum(inns2(0:nnnhb,ith,igr)))
            nnav(ith,igr) = nnav(ith,igr) + term
            nnsd(ith,igr) = nnsd(ith,igr) + term**2
            ihbs1(0:nnnhb,ith,igr) = ihbs1(0:nnnhb,ith,igr)+ihbs2(0:nnnhb,ith,igr)
            aver = sum( [ (ihb, ihb = 0,nnnhb) ]*ihbs2(0:nnnhb,ith,igr))
            term = aver*InvInt(sum(ihbs2(0:nnnhb,ith,igr)))
            hbav(ith,igr) = hbav(ith,igr) + term
            hbsd(ith,igr) = hbsd(ith,igr) + term**2
         end do
      end do

      if (lsim .and. master) write(ucnf) nsamp1, inns1, ihbs1, nnav, nnsd, hbav, hbsd

   case (iAfterSimulation)

      norm = InvInt(nsamp1)
      norm1 = InvInt(nsamp1-1)
      do igr = 1, ngr(1)
         do ith = 1, nthnn
            nnnorm(ith,igr) = InvInt(sum(inns1(0:nnnhb,ith,igr)))
            nnav(ith,igr) = nnav(ith,igr)*norm
            nnsd(ith,igr) = sqrt(max(Zero,nnsd(ith,igr)*norm-nnav(ith,igr)**2)*norm1)
            hbnorm(ith,igr) = One*InvInt(sum(ihbs1(0:nnnhb,ith,igr)))
            hbav(ith,igr) = hbav(ith,igr)*norm
            hbsd(ith,igr) = sqrt(max(Zero,hbsd(ith,igr)*norm-hbav(ith,igr)**2)*norm1)
         end do
      end do

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,i8  )') 'max number of hydrogen bonds   = ', nnnhb
      write(uout,'(a,t35,f8.2)') 'upper separation distance      = ', rmax
      do igr = 1, ngr(1)
         write(uout,'()')
         write(uout,'()')
         write(uout,'(a)') txgr(igr)
         write(uout,'()')
         write(uout,'(a,t30,a)') 'nearest neighbour', 'hydrogen bonds'
         write(uout,'(a,t30,a)') '-----------------', '--------------'
         write(uout,'()')
         write(uout,'(a,2f8.1,t30,a,2f8.1)') 'no', sqrt(thnn(1:nthnn)), 'no', thhb(1:nthhb)
         write(uout,'()')
         do ihb = 0, nnnhb
            write(uout,'(i2,2f8.2,t30,i2,2f8.2)')               &
               ihb, inns1(ihb,1:nthnn,igr)*nnnorm(1:nthnn,igr), ihb, ihbs1(ihb,1:nthnn,igr)*hbnorm(1:nthnn,igr)
         end do
         write(uout,'()')
         write(uout,'(a,2f8.2,t30,a,2f8.2)') 'av', nnav(1:nthnn,igr), 'av', hbav(1:nthhb,igr)
         write(uout,'(a,1x,2f8.3,t30,a,1x,2f8.3)') 'sd', nnsd(1:nthnn,igr), 'sd', hbsd(1:nthhb,igr)
      end do

      deallocate(inns2, inns1, ihbs2, ihbs1)
      deallocate(nnav, nnsd, hbav, hbsd, nnnorm, hbnorm)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine NNHB

!************************************************************************
!> \page static static.F90
!! **NNDF**
!! *calculate the nearest neighbour distribution*
!************************************************************************


! ... p(r) =< delta(r-r') > for particles belonging to group igr
!     where r is the distance between the particle ip and
!     the nearest particle belonging to group jgr

!> \page nmlNNDF
!! The namelist  \ref nmlNNDF contains variables that control the calculation of nearest neighbour distribution functions. The R-F
!! nearest neighbour distribution function gives the distribution of distances between a reference particle R and the nearest field
!! particle F.
!! * Variables:
!!  * \subpage nmlNNDF_vtype

!> \page nmlNNDF_vtype vtype
!! `static2D_var(logical, real, real, integer, logical, character, real)`
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization title, and number of variable of vtype.
!! * Min: /0.0/ Max: /10.0/

   subroutine NNDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='NNDF'
   character(80), parameter :: txheading ='nearest-neighbour distribution function'
   integer(4)   , parameter :: ntype = 1
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   real(8)     , allocatable,  save :: r2low(:)
   integer(4) :: ivar, ibin, ip, jp, igr, jgr
   real(8)    :: dx, dy, dz, r1, r2

   namelist /nmlNNDF/ vtype

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = Zero
      vtype%max = 10.0d0
      vtype%nbin   = 100

      rewind(uin)
      read(uin,nmlNNDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%nvar = ngrgr

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), r2low(maxngr))
      r2low = 0.0E+00

! ... set ipnt, label, min, max, and nbin

      var%label = txgrgr
      var%min = vtype(1)%min
      var%max = vtype(1)%max
      var%nbin = vtype(1)%nbin
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ip = 1, np
         igr = igrpn(ip,1)
         if (igr == 0) cycle
         r2low(1:ngr(2)) = huge(One)
         do jp = 1, np
            if (jp == ip) cycle
            jgr = igrpn(jp,2)
            if (jgr == 0) cycle
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            r2low(jgr) = min(r2low(jgr),r2)
         end do
         do jgr = 1, ngr(2)
            ivar = igrgr(igr,jgr)
            r1 = sqrt(r2low(jgr))
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, r2low)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine NNDF

!************************************************************************
!> \page static static.F90
!! **ChainDF**
!! *calculate chain distribution functions*
!************************************************************************


!     also final average over chains of same type and the spread of their df

!> \page nmlChainDF
!! The namelist  \ref nmlChainDF contains variables that control the calculation of chain distribution functions. Distribution functions
!! are calculated for each chain. Any combination of the types of distribution functions listed below may be selected through vtype\%l.
!!
!!    | type | label | quantity                         |
!!    | ---- | ------| -------------------------------- |
!!    | 1    | rbb   | bead-bead separation             |
!!    | 2    | ree   | end-to-end separation            |
!!    | 3    | rg    | radius of gyration               |
!!    | 4    | angle | angle between consequtive beads  |
!!    | 5    | cos   | cos(180-angle)                   |
!!    | 6    | shape | ree**2/rg**2 ratio               |
!!    | 7    | asph  | asphericity                      |
!!    | 8    | torp  | toroidicity                      |
!!
!! * Variables:
!!  * \subpage nmlChainDF_vtype

!> \page nmlChainDF_vtype vtype
!! `static2D_var(logical, real, real, integer, logical, character, real)`(1:8)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization title, and number of variable of vtype.
!! * Min: /0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0/ Max: /12.0, 100.0, 100.0, 180.0, 1.0, 12.0, 1.0, 1.0/

subroutine ChainDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainDF'
   character(80), parameter :: txheading ='chain distribution functions'
   integer(4)   , parameter :: ntype = 8
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:)
   integer(4),                 save :: nvar2
   type(df_var), allocatable,  save :: var2(:)
   integer(4),   allocatable        :: ilow(:), iupp(:)
   real(8),      allocatable        :: vspread(:)
   type(chainprop_var) :: ChainProperty

   integer(4)     :: itype, ivar, ibin, ic, ict, ivar2
   real(8)        :: value
   character(3)   :: txcn

   namelist /nmlChainDF/ vtype

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (slave) return                   ! only master

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = [ Zero   , Zero  , Zero  ,  Zero , -One  , Zero  , Zero   , Zero  ]
      vtype%max = [ 12.d0  , 100.d0, 100.d0, 180.d0,  One  , 12.d0 ,  One   ,  One  ]
      vtype%nbin = 100

      rewind(uin)
      read(uin,nmlChainDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['rbb  ','ree  ','rg   ','angle','cos  ','shape','asph ' ,'torp ']
      vtype%nvar = nc

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nc,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ic = 1, nc

               if (lhierarchical) then          ! if hierarchical structure only chain type 1
                  if((ihnpn(ipnsegcn(1,ic)) /= 0) .and. (ictcn(ic) /= 1)) cycle
               end if

               ivar = ivar+1
               ipnt(ic,itype) = ivar
               write(txcn,'(i3)') ic
               var(ivar)%label = trim(vtype(itype)%label)//' '//txcn
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      ! set nvar accoding to ivar, to skipped chain if hierachical structure
      nvar = ivar
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc

         if (lhierarchical) then          ! if hierarchical structure only chain type 1
            if((ihnpn(ipnsegcn(1,ic)) /= 0) .and. (ictcn(ic) /= 1)) cycle
         end if

         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)

! ... sample type 1 to type 8

         do itype = 1, 8
            if (vtype(itype)%l) then
               ivar = ipnt(ic,itype)
               if (itype == 1) then
                  value = sqrt(ChainProperty%rbb2)
               else if (itype == 2) then
                  value = sqrt(ChainProperty%ree2)
               else if (itype == 3) then
                  value = sqrt(ChainProperty%rg2)
               else if (itype == 4) then
                  value = ChainProperty%angle
               else if (itype == 5) then
                  value = ChainProperty%cos
               else if (itype == 6) then
                  value = ChainProperty%ree2/ChainProperty%rg2
               else if (itype == 7) then
                  value = ChainProperty%asph
               else if (itype == 8) then
                  value = ChainProperty%torp
               end if
               ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do

      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

! ..............    make average of distribution functions over chains of same type   .............

! ... set nvartype2 and nvar2 as well as allocate memory

      vtype%nvar = nct
      nvar2 = sum(vtype%nvar, 1, vtype%l)
      allocate(ilow(nvar2), iupp(nvar2), var2(nvar2), vspread(nvar2))
      ilow = 0
      iupp = 0
      vspread = 0.0E+00

! ... set label, min, max, nbin, and bin

      ivar2 = 0
      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar2 = ivar2+1
               ivar = ivar+1
               ilow(ivar2) = ivar                    ! lower index of chains of type ict
               ivar = ivar+ncct(ict)-1
               iupp(ivar2) = ivar                    ! upper index of chains of type ict
               var2(ivar2)%label = trim(vtype(itype)%label)//' '//txct(ict)
               var2(ivar2)%min = vtype(itype)%min
               var2(ivar2)%max = vtype(itype)%max
               var2(ivar2)%nbin = vtype(itype)%nbin
               var2(ivar2)%bin = (var2(ivar2)%max-var2(ivar2)%min)/var2(ivar2)%nbin
            end do
         end if
      end do

      call DistFuncAverDist(nvar2, ilow, iupp, var, var2, vspread)

      call DistFuncWrite(trim(txheading)//': aver', nvar2, var2, uout, ulist, ishow, iplot, ilist)
      write(uout,'()')
      write(uout,'(a,(t30,f10.4))') 'rms spread = ', vspread

      call DistFuncAverValue(nvar2, var2, uout)

      deallocate(var, ipnt, ilow, iupp, var2, vspread)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainDF

!************************************************************************
!> \page static static.F90
!! **ChainTypeDF**
!! *calculate chain type distribution functions*
!************************************************************************


!     type  label  quantity
!     ----  -----  --------
!     1     rbb    bead-bead separation
!     2     ree    end-to-end separation
!     3     rg     radius of gyration
!     4     angle  angle between consequtive beads
!     5     cos    cos(180-angle)
!     6     shape  ree**2/rg**2 ratio
!     7     asph   asphericity
!     8     torp   toroidicity

!> \page nmlChainTypeDF
!! The namelist  \ref nmlChainTypeDF contains variables that control the calculation of chain distribution functions. Distribution
!! functions are calculated for each type of chains. Any combination of the types of distribution functions listed below may be
!! selected through vtype\%l. Same input variables as for \ref nmlChainDF.

subroutine ChainTypeDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainTypeDF'
   character(80), parameter :: txheading ='chain type distribution functions'
   integer(4)   , parameter :: ntype = 8
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:)
   type(chainprop_var) :: ChainProperty
   integer(4)     :: itype, ivar, ibin, ic, ict
   real(8)        :: value

   namelist /nmlChainTypeDF/ vtype

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (slave) return                   ! only master

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = [ Zero   , Zero  , Zero  ,  Zero , -One  , Zero  , Zero   , Zero  ]
      vtype%max = [ 12.d0  , 100.d0, 100.d0, 180.d0,  One  , 12.d0 ,  One   ,  One  ]
      vtype%nbin = 100

      rewind(uin)
      read(uin,nmlChainTypeDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['rbb  ','ree  ','rg   ','angle','cos  ','shape','asph ' ,'torp ']
      vtype%nvar = nct

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nct,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ivar+1
               ipnt(ict,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc
         ict = ictcn(ic)
         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)

! ... sample type 1 to type 8

         do itype = 1, 8
            if (vtype(itype)%l) then
               ivar = ipnt(ict,itype)
               if (itype == 1) then
                  value = sqrt(ChainProperty%rbb2)
               else if (itype == 2) then
                  value = sqrt(ChainProperty%ree2)
               else if (itype == 3) then
                  value = sqrt(ChainProperty%rg2)
               else if (itype == 4) then
                  value = ChainProperty%angle
               else if (itype == 5) then
                  value = ChainProperty%cos
               else if (itype == 6) then
                  value = ChainProperty%ree2/ChainProperty%rg2
               else if (itype == 7) then
                  value = ChainProperty%asph
               else if (itype == 8) then
                  value = ChainProperty%torp
               end if
               ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do

      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainTypeDF

!************************************************************************
!> \page static static.F90
!! **ChainTypeExtDF**
!! *calculate chain type distribution function relative external frame*
!************************************************************************


!> \page nmlChainTypeExtDF
!! The namelist  \ref nmlChainTypeExtDF contains variables that control the calculation of chain distribution functions where the
!! distribution functions are with respect to the lab frame. Distribution functions are calculated for each type of chains. Any
!! combination of the types of distribution functions listed below may be selected through vtype\%l.
!!
!!    | type | label          | quantity                                                |
!!    | ---- | -------------- | ------------------------------------------------------- |
!!    | 1    | \format{com_x} | com, x-direction                                        |
!!    | 2    | com_y          | com, y-direction                                        |
!!    | 3    | com_z          | com, z-direction                                        |
!!    | 4    | rg_z           | square of radius of gyration projected on the z-axis    |
!!    | 5    | rg_xy          | square of radius of gyration projected on the xy-plane  |
!!
!! * Variables:
!!  * \subpage nmlChainTypeExtDF_vtype

!> \page nmlChainTypeExtDF_vtype vtype
!! `static2D_var(logical, real, real, integer, logical, character, real)`(1:5)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization title, and number of variable of vtype.
!! * Min: /-box2(1),-box2(2),-box2(3),0.0,0.0/
!! * Max: /box2(1),box2(2),box2(3),100.0,100.0/

subroutine ChainTypeExtDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainTypeExtDF'
   character(80), parameter :: txheading ='chain type distribution functions'
   integer(4)   , parameter :: ntype = 5
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:)
   type(chainprop_var) :: ChainProperty
   integer(4)     :: itype, ivar, ibin, ic, ict
   real(8)        :: value

   namelist /nmlChainTypeExtDF/ vtype

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l =.false.
      vtype%min = [-boxlen2(1), -boxlen2(2), -boxlen2(3),  Zero ,  Zero  ]
      vtype%max = [ boxlen2(1),  boxlen2(2),  boxlen2(3),100.d0 ,100.d0  ]
      vtype%nbin = 100

      rewind(uin)
      read(uin,nmlChainTypeExtDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['com_x','com_y','com_z','rg_z ','rg_xy']
      vtype%nvar = nct

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(nct,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ivar+1
               ipnt(ict,itype) = ivar
               var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
               var(ivar)%min = vtype(itype)%min
               var(ivar)%max = vtype(itype)%max
               var(ivar)%nbin = vtype(itype)%nbin
            end do
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc
         ict = ictcn(ic)

         call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
         call CalcChainProperty(ic, vaux, ChainProperty)

! ... sample type 1 to type 5

         do itype = 1, 5
            if (vtype(itype)%l) then
               ivar = ipnt(ict,itype)
               if (itype == 1) then
                  value = ChainProperty%ro(1)
               else if (itype == 2) then
                  value = ChainProperty%ro(2)
               else if (itype == 3) then
                  value = ChainProperty%ro(3)
               else if (itype == 4) then
                  value = sqrt(ChainProperty%rg2z)
               else if (itype == 5) then
                  value = sqrt(ChainProperty%rg2xy)
               end if
               ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do
      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      call DistFuncAverValue(nvar, var, uout)

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainTypeExtDF

!************************************************************************
!> \page static static.F90
!! **ChainBeadPartContact**
!! *calculate probability of chain bead-particle contact*
!************************************************************************


!> \page nmlCBPC
!! The namelist  \ref nmlCBPC contains variables that control the calculation of probabilities that chain particles are near particles of
!! other types. Distribution functions are calculated for each pair of chain molecules and particles of type \ref iptpart.
!! * Variables:
!!  * \subpage iptpart
!!  * \subpage rcontact

!> \page iptpart
!! `integer`
!! * Type of particle to be considered.

!> \page rcontact
!! `real`
!! * Largest separation between chain particles and particles of the other type for considering the chain particle near the other particle.

   subroutine ChainBeadPartContact(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ChainBeadPartContact'
   character(80), parameter :: txheading ='probability of chain bead-particle contacts'
   integer(4)   , parameter :: ntype = 10
   integer(4),                 save :: nbin
   integer(4),                 save :: nvar
   integer(4),                 save :: iptpart
   real(8),                    save :: rcontact
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:)
   integer(4) :: iseg, ic, ict, ip, iploc, jp, jploc, ivar, ibin
   real(8) ::  dx, dy, dz, r2
   character(6) :: str1, str2

   namelist /nmlCBPC/ iptpart, rcontact

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      rewind(uin)
      read(uin,nmlCBPC)

      if (nct > 1) call Stop(txroutine, 'nct > 1', uout)
      if (nppt(iptpart) > ntype) call Stop(txroutine, 'nppt(iptpart) > ntype', uout)

   case (iWriteInput)

      nbin = npct(1)
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

! ... set nvar as well as allocate memory

      nvar = nc*nppt(iptpart)
      allocate(var(nvar), ipnt(nc, ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do ic = 1, nc                         ! loop over all chains
         ict = ictcn(ic)
         do iploc = 1, nppt(iptpart)        ! loop over all particles of type iptpart
            ivar = ivar+1
            ipnt(ic,iploc) = ivar
            write(str1,'(i6)') ic
            write(str2,'(i6)') ipnpt(iptpart)+(iploc-1)
            var(ivar)%label = 'c.'//trim(adjustl(str1))//' p.'//trim(adjustl(str2))
            var(ivar)%min = 1-0.5
            var(ivar)%max = npct(ict)+0.5
            var(ivar)%nbin = nbin
         end do
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do ic = 1, nc
         ict = ictcn(ic)
         do iseg = 1, npct(ict)
            ip = ipnsegcn(iseg,ic)
            do jploc = 1, nppt(iptpart)
               jp = ipnpt(iptpart)+(jploc-1)
               if (ip == jp) call Stop(txroutine, 'ip == jp', uout)
               dx = ro(1,ip)-ro(1,jp)
               dy = ro(2,ip)-ro(2,jp)
               dz = ro(3,ip)-ro(3,jp)
               call PBCr2(dx,dy,dz,r2)
               if (r2 < rcontact**2) then
                  ivar = ipnt(ic,jploc)
                  ibin = iseg-1
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if
            end do
         end do
      end do

   case (iAfterMacrostep)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,i8  )') 'particle type used             = ', iptpart
      write(uout,'(a,t35,f8.2)') 'upper contact distance         = ', rcontact
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

     deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ChainBeadPartContact

!************************************************************************
!> \page static static.F90
!! **LoopTailTrain**
!! *calculate loop, tail, and tail characteristics for chains*
!************************************************************************


!     type  label
!     ----  -----
!     1     number of loops
!     2     number of tails
!     3     number of trains
!     4     length of loops
!     5     length of tails
!     6     length of trains
!     7     number of segments in loops
!     8     number of segments in tails
!     9     number of segments in trains
!     10    number of adsorbed chains
!     11    number of adsorbed segments

!> \page nmlLoopTailTrain
!! The namelist  \ref nmlLoopTailTrain contains variables that control the calculation loop, tail, and train characteristics. Properties are calculated for each type of chains separately.
!! * Variables:
!!  * \subpage adscond

!> \page adscond
!! `adscond_var(character, character, real)`
!! * Adsorption condition('xy-plane'), selection of box ends, threshold distance of adsorption.

subroutine LoopTailTrain(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='LoopTailTrain'
   character(80), parameter :: txheading ='loop, tail, and train characteristics'
   integer(4)   , parameter :: ntype = 11
   integer(4)   , parameter :: mnobj = 100
   type(adscond_var),             save :: adscond
   integer(4)      ,              save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   integer(4),       allocatable, save :: nadschain(:), nadsseg(:)

   logical    :: ladschain                           ! .true. if chain is adsorbed
   logical,          allocatable, save :: ladsseg(:) ! .true. if segment is adsorbed
   integer(4) :: nobj(3)                             ! number of objects
   integer(4) :: lenobj(3,mnobj)                     ! length of objects
   integer(4) :: nsegobj(3)                          ! number of segments involved
   integer(4) :: ic, ict, i, ioffset, iobj

   namelist /nmlLoopTailTrain/ adscond

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      rewind(uin)
      read(uin,nmlLoopTailTrain)

      call LowerCase(adscond%txobject)
      call LowerCase(adscond%txsurface)

      if (.not.allocated(ladsseg)) then
         allocate(ladsseg(maxval(npct(1:nct))))
         ladsseg = .false.
      end if

! ... set nvar as well as allocate memory

      nvar = nct*ntype
      allocate(var(nvar), nadschain(nct), nadsseg(nct))
      nadschain = 0
      nadsseg = 0

! ... set label

      do ict = 1, nct
         ioffset = ntype*(ict-1)                                     ! offset for different chain types
         var(1+ioffset )%label = 'number of loops                = '
         var(2+ioffset )%label = 'number of tails                = '
         var(3+ioffset )%label = 'number of trains               = '
         var(4+ioffset )%label = 'length of loops                = '
         var(5+ioffset )%label = 'length of tails                = '
         var(6+ioffset )%label = 'length of trains               = '
         var(7+ioffset )%label = 'number of segments in loops    = '
         var(8+ioffset )%label = 'number of segments in tails    = '
         var(9+ioffset )%label = 'number of segments in trains   = '
         var(10+ioffset)%label = 'number of adsorbed chains      = '
         var(11+ioffset)%label = 'number of adsorbed segments    = '
      end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      nadsseg = 0
      nadschain = 0
      do ic = 1, nc                                                      ! loop over chains
         ict = ictcn(ic)
         ioffset = ntype*(ict-1)
         call CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)
         if (ladschain) then
            call CalcLTT(npct(ict),ladsseg,mnobj,nobj,lenobj,nsegobj)    ! calculate nobj, lenobj, and nsegobj
            do i = 1, 3                                                  ! loop over the different characteristics
               var(i+ioffset)%value = nobj(i)                            ! assign number of objects
               call ScalarSample(iStage, i+ioffset, i+ioffset, var)
               do iobj = 1, nobj(i)
                  var(i+3+ioffset)%value = lenobj(i,iobj)                ! assign the lengths of the objects
                  call ScalarSample(iStage, i+3+ioffset, i+3+ioffset, var)
               end do
               var(i+6+ioffset)%value = nsegobj(i)                       ! assign the total number of segments in these objects
               call ScalarSample(iStage, i+6+ioffset, i+6+ioffset, var)
            end do
            nadschain(ict) = nadschain(ict) + 1
            nadsseg(ict) = nadsseg(ict) + nsegobj(3)
         end if
      end do
      do ict = 1, nct
         ioffset = ntype*(ict-1)
         var(10+ioffset)%value = nadschain(ict)
         call ScalarSample(iStage, 10+ioffset, 10+ioffset, var)
         var(11+ioffset)%value = nadsseg(ict)
         if (nadsseg(ict) > 0) call ScalarSample(iStage, 11+ioffset, 11+ioffset, var)
      end do

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t40,a    )') 'adsorbing object               = ', adscond%txobject
      if (adscond%txobject == 'plane_xy') &
      write(uout,'(a,t40,a    )') 'adsorbing surface              = ', adscond%txsurface
      write(uout,'(a,t35,f10.3)') 'adsorbing distance             = ', adscond%dist
      do ict = 1, nct
         ioffset = ntype*(ict-1)
         write(uout,'()')
         if (nct > 1) write(uout,'(2a)') 'chain: ', txct(ict)
         if (nct > 1) write(uout,'(a)')  '------------------'
         call ScalarWrite(iStage, 1+ioffset, ntype+ioffset, var, 1, '(a,t35,4f15.3,f15.0)', uout)
      end do

      deallocate(var)
      deallocate(ladsseg)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine LoopTailTrain

!************************************************************************
!> \page static static.F90
!! **CheckAdsChainSeg**
!! *check if a chain and its segments are adsorbed*
!************************************************************************


subroutine CheckAdsChainSeg(ic, adscond, ladschain, ladsseg)

   use MolModule
   implicit none

   integer(4)       , intent(in)  :: ic           ! chain to be considered
   type(adscond_var), intent(in)  :: adscond      ! adsorbing conditions
   logical          , intent(out) :: ladschain    ! flag for adsorbed chain
   logical          , intent(out) :: ladsseg(*)   ! flag for adsorbed segments

   character(3) :: txno
   real(8)      :: zloc, dist2, dx, dy, dz, r2
   integer(4)   :: npchain, iseg, ip, jp

   npchain = npct(ictcn(ic))                  ! number of particles in the chain
   ladsseg(1:npchain) = .false.               ! initialize ladsseg
   if (adscond%txobject == 'plane_xy') then   ! segment adsorbed according to its z-coordinate
      zloc = boxlen2(3) - adscond%dist
      do iseg = 1, npchain                    ! loop over segments, check if segment is adsorbed
         ip = ipnsegcn(iseg,ic)
         if (adscond%txsurface == '-') then                   ! check lower surface
            if (ro(3,ip) < -zloc) ladsseg(iseg) = .true.
         else if (adscond%txsurface == '+') then              ! check upper surface
            if (ro(3,ip) >  zloc) ladsseg(iseg) = .true.
         else if (adscond%txsurface == '-|+') then            ! check both surfaces
            if (abs(ro(3,ip)) > zloc) ladsseg(iseg) = .true.
         else
            call Stop('CheckAdsChainSeg', 'error in adscond%txsurface', uout)
         end if
      end do
   else if (adscond%txobject(1:4) == 'part') then
      txno = adscond%txobject(6:8)
      read(txno,'(i3)') jp                    ! jp is the id of the "adsorbing" particle
      dist2 = adscond%dist**2
      do iseg = 1, npchain                    ! loop over segments, check if segment is adsorbed at jp
         ip = ipnsegcn(iseg,ic)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 < dist2) ladsseg(iseg) = .true.
      end do
   end if

   ladschain = .false.
   if (count(ladsseg(1:npchain)) > 0) ladschain = .true. ! check if chain is adsorbed

end subroutine CheckAdsChainSeg

!************************************************************************
!> \page static static.F90
!! **ClusterSD**
!! *calculate cluster size distribution*
!************************************************************************


!> \page nmlCluster
!! The namelist  \ref nmlCluster contains variables that control the calculation of cluster size distribution functions. Objects within
!! the distance \ref rcluster from a given object forms a virtual bond and all objects directly or indirectly bonded form a cluster.
!! * Variables:
!!  * \subpage txobj
!!  * \subpage l1d
!!  * \subpage l2d
!!  * \subpage lpercolation
!!  * \subpage nobjt
!!  * \subpage iobjt
!!  * \subpage rcluster
!!  * \subpage txweight
!!  * \subpage itestcluster

!> \page txobj
!! `character(8)`
!! **default:** `'particle'`
!! * Object for which cluster analysis should be performed.
!! * ``'particle'``:  Particles, \ref rcluster applies between particles.
!! * ``'chain'``:  Chains, \ref rcluster applies between chain particles.
!! * ``'chaincom'``: Chains, \ref rcluster applies between chain com's.

!> \page l1d
!! `logical`
!! **default:** `.true.`
!! * `.true.`:  Calculation of 1d cluster size distribution functions.
!! * `.false.`:  No such calculation.

!> \page l2d
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Calculation of 2d cluster size distribution functions.nptcluster=2 is required.
!! * `.false.`:  No such calculation.


!> \page lpercolation
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Percolation analysis is made.
!! * `.false.`:  No such analysis.

!> \page nobjt
!! `integer`
!! **default:** `1`
!! * Number of object types to be considered.

!> \page iobjt
!! `integer`(1:\ref npt)
!! **default:** \ref npt*`0`
!! * Object types to be considered.

!> \page rcluster
!! `real`
!! * Upper separation for bonded objects (particles if \ref txobj='particle' or 'chain' and chain com if \ref txobj='chaincom').

!> \page txweight
!! `character(6)`
!! **default:** `'mass'`
!! * ``'mass'``: Mass-weighted distribution.
!! * ``'number'``: Number-weighted distribution.

!> \page itestcluster
!! `integer`
!! **default:** `0`
!! * Flag for test output. This possibility is for maintenance purposes.
!! * `0`:  Nothing. The normal option.
!! * `1`:  Intermediate cluster data (several routines in static.F90). Type of particles to be inserted in the different sets.

   subroutine ClusterSD(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ClusterSD'
   character(80), parameter :: txheading ='cluster size distribution'
   integer(4)   , parameter :: ntype = 1                 ! number of types of properties
   integer(4)   , parameter :: mnpvar = 1                ! number of variables
   integer(4)   , parameter :: mnsizemax = mnbin_df      ! maximal size of largest cluster for 1D smapling
   integer(4)   , parameter :: mnsizemax2d = mnbin_df2d  ! maximal size of largest cluster for 2D sampling
   integer(4)   , parameter :: mnobj = 1000              ! maximal number of objects to be considered
   integer(4)   , parameter :: mnpair = 100000           ! maximal number of pairs
   character(8),     save :: txobj                       !*='particle' particle cluster size distribution
                                                         !*='chain' chain cluster size distribtuion based on particle-particle separation
                                                         !*='chaincom' chain cluster size distribtuion based on com separation
   logical,          save :: l1d                         !*flag for 1d cluster size distribution function
   logical,          save :: l2d                         !*flag for 2d cluster size distribution function
   logical,          save :: lpercolation                !*flag for perculation calculation
   integer(4),       save :: nvartype(ntype)
   integer(4),       save :: nobjt                       !*number of object types to be considered
   integer(4), allocatable, save :: iobjt(:)             !*object types to be considered (global number)
   integer(4),       save :: nobj                        ! number of objects to be considered
   integer(4),       save :: iobj(mnobj)                 ! objects to be considered (global number)
   integer(4),       save :: iobjtloc(mnobj)             ! object types to be considered (local number)
   real(8),          save :: rcluster, rcluster2         !*separation for "connected" objects
   character(6),     save :: txweight                    !*='number' number-weighted distribution
                                                         !*='mass' mass-weighted distribution
   integer(4),       save :: itestcluster                !*=1, test output

   integer(4),       save :: npvar
   type(scalar_var), allocatable, save :: pvar(:)
   character(10),    allocatable, save :: txobjt(:)      ! label of the object types
   integer(4),       save :: nvar
   type(df_var), allocatable,   save :: var(:)
   type(df2d_var), allocatable, save :: vvar(:)

   logical, save :: ltest
   integer(4) :: npair, n1(mnpair), n2(mnpair)
   integer(4) :: iclusteriobj(mnobj)                     ! pointer: object -> its cluster
   integer(4) :: ncluster, nsizemax, nsizedist(mnsizemax)
   integer(4) :: nsizemax2d(2), nsizedist2d(0:mnsizemax2d,0:mnsizemax2d)
   integer(4) :: iptloc, ipt, ip, ictloc, ict, ic, isize, isize1, isize2
   integer(4), save :: nsizemax_save = 0, nsizemax2d_save(2) = 0
   logical :: Perculation, lreturn

   namelist /nmlCluster/ txobj, l1d, l2d, lpercolation, nobjt, iobjt, rcluster, txweight, itestcluster

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(iobjt)) then
         allocate(iobjt(npt))
         iobjt = 0
      end if

      txobj = 'particle'
      l1d = .true.
      l2d = .false.
      lpercolation = .false.
      nobjt = 1
      txweight = 'mass'
      itestcluster = 0

      rewind(uin)
      read(uin,nmlCluster)

      call LowerCase(txobj)
      call LowerCase(txweight)

      ltest =.false.
      if (itestcluster == 1) ltest =.true.

   case (iWriteInput)

     if (.not.allocated(txobjt)) then
        allocate(txobjt(nobjt))
     end if

     if (nobjt > npt) call Stop(txroutine, 'nobjt > npt', uout)
     if (l2d .and. nobjt/= 2)  call Stop(txroutine, 'l2d .and. nobjt/= 2', uout)

! ... calculate nobj and set iobj and iobjtloc

      if (txobj == 'particle') then
         nobj = 0
         do iptloc = 1, nobjt
            ipt = iobjt(iptloc)
            do ip = ipnpt(ipt), ipnpt(ipt)-1+nppt(ipt)
               nobj = nobj + 1
               if (nobj > mnobj) call Stop(txroutine, 'nobj > mnobj', uout)
               iobj(nobj) = ip                 ! global object number
               iobjtloc(nobj) = iptloc         ! local object type
            end do
         end do
         txobjt(1:nobjt) = txpt(iobjt(1:nobjt))
      else if (txobj(1:5) == 'chain') then
         nobj = 0
         do ictloc = 1, nobjt
            ict = iobjt(ictloc)
            do ic = icnct(ict), icnct(ict)-1+ncct(ict)
               nobj = nobj + 1
               if (nobj > mnobj) call Stop(txroutine, 'nobj > mnobj', uout)
               iobj(nobj) = ic                 ! global object number
               iobjtloc(nobj) = ictloc         ! local object type
            end do
         end do
         txobjt(1:nobjt) = txct(iobjt(1:nobjt))
      else
         call stop(txroutine, 'error in txobj', uout)
      end if

! ... set nvartype and nvar as well as allocate memory

      nvartype = 1
      nvar = sum(nvartype)
      allocate(var(nvar), vvar(nvar))

! ... set label

      var(1)%label = '1d: '
      do iptloc = 1, nobjt
         var(1)%label = trim(var(1)%label)//' '//txobjt(iptloc)
      end do
      var(1)%min =+Half
      var(1)%max = mnbin_df+var(1)%min
      var(1)%nbin = mnbin_df
      call DistFuncSample(iStage, nvar, var)

      vvar(1)%label = '2d: '
      do iptloc = 1, nobjt
         vvar(1)%label = trim(vvar(1)%label)//' '//txobjt(iptloc)
      end do
      vvar(1)%min =+Half
      vvar(1)%max = mnbin_df2d+vvar(1)%min
      vvar(1)%nbin = mnbin_df2d
      call DistFunc2DSample(iStage, nvar, vvar)

      npvar = nvar
      allocate(pvar(nvar))
      pvar(1)%label = 'perculation'
      pvar(1)%norm  = One

      rcluster2 = rcluster**2

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      call DistFunc2DSample(iStage, nvar, vvar)
      call ScalarSample(iStage, 1, npvar, pvar)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) vvar
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) pvar

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      call DistFunc2DSample(iStage, nvar, vvar)
      call ScalarSample(iStage, 1, npvar, pvar)

   case (iSimulationStep)

      if (itestcluster == 1) call WriteHead(3, 'TestCluster', uout)

      var%nsamp2 = var%nsamp2 + 1
      vvar%nsamp2 = vvar%nsamp2 + 1

! ... get a list of object pairs

      if (txobj == 'particle') then
         call CalcPartPairListGeneral(nobj, iobj, rcluster2, mnpair, npair, n1, n2)
      else if (txobj == 'chain') then
         call CalcChainPairListGeneral(nobj, iobj, rcluster2, mnpair, npair, n1, n2)
      else if (txobj == 'chaincom') then
         call CalcChainComPairListGeneral(nobj, iobj, rcluster2, mnpair, npair, n1, n2)
      end if

! ... make clusters

      call MakeCluster(nobj, npair, n1, n2, iclusteriobj)

! ... calculate cluster size distribution

      call ClusterSD1D(mnobj, nobj, iclusteriobj, mnsizemax, ltest, uout, ncluster, nsizemax, nsizedist)

! ... calculate 2d cluster size distribution

      if (l2d) &
      call ClusterSD2D(mnobj, nobj, iclusteriobj, iobjtloc, mnsizemax2d, ltest, uout, nsizemax2d, nsizedist2d)

! ... number or mass average

      if (txweight == 'mass') then
         do isize = 1, nsizemax
            nsizedist(isize) = isize*nsizedist(isize)
         end do
         do isize1 = 0, nsizemax2d(1)
            do isize2 = 0, nsizemax2d(2)
               nsizedist2d(isize1,isize2) = (isize1+isize2)*nsizedist2d(isize1,isize2)
            end do
         end do
      else if (txweight == 'number') then
      else
         call Stop(txroutine, 'error in txweight', uout)
      end if

      if (itestcluster == 1) then
         write(uout,'(a)') '2d cluster size distribution'
         do isize1 = 0, nsizemax2d(1)
            write(uout,'((21i4))') nsizedist2d(isize1,0:nsizemax2d(2))
         end do
      end if

! ... perculation

      lreturn = .false.
      if (lpercolation) then

         if (ltest) then
            write(uout,*)
            write(uout,*) 'In TestPerculation'
            write(uout,*) 'nobj = ',nobj
            write(uout,*) 'iclusteriobj(1:nobj)', iclusteriobj(1:nobj)
            write(uout,*) 'npair =', npair
            write(uout,*) 'n1(1:npair)', n1(1:npair)
            write(uout,*) 'n2(1:npair)', n2(1:npair)
         end if
         lreturn = Perculation(nobj, iclusteriobj, npair, n1, n2, r, boxlen, ltest, uout)
         if (ltest) then
            write(uout,*)
            write(uout,*) 'In TestPerculation'
            write(uout,*) 'lpercolation', lreturn
         end if
      end if

! ... accumulate

      var(1)%avs2(0:nsizemax-1) = var(1)%avs2(0:nsizemax-1)+nsizedist(1:nsizemax)
      vvar(1)%avs2(0:nsizemax2d(1),0:nsizemax2d(2)) = &
         vvar(1)%avs2(0:nsizemax2d(1),0:nsizemax2d(2))+nsizedist2d(0:nsizemax2d(1),0:nsizemax2d(2))
      pvar(1)%value = 0
      if (lreturn) pvar(1)%value = 1
      call ScalarSample(iStage, 1, nvar, pvar)

! ... save largest cluster sizes

      nsizemax_save = max(nsizemax_save,nsizemax)
      nsizemax2d_save(1) = max(nsizemax2d_save(1),nsizemax2d(1))
      nsizemax2d_save(2) = max(nsizemax2d_save(2),nsizemax2d(2))

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      call DistFunc2DNorm(1, nvar, vvar)
      call DistFunc2DSample(iStage, nvar, vvar)
      call ScalarSample(iStage, 1, npvar, pvar)
      call ScalarNorm(iStage, 1, npvar, pvar, 0)
      if (lsim .and. master) write(ucnf) var
      if (lsim .and. master) write(ucnf) vvar
      if (lsim .and. master) write(ucnf) pvar

      call WriteAverage    !!! 2005-05-01

   case (iAfterSimulation)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,a   )') 'object                         = ', txobj
      write(uout,'(a,t35,i8  )') 'number of object types         = ', nobjt
      write(uout,'(a,t35,5i8 )') 'object types (number)          = ', iobjt(1:nobjt)
      write(uout,'(a,t35,5a  )') 'object types (label)           = ', txobjt(1:nobjt)
      write(uout,'(a,t35,f8.3)') 'rcluster                       = ', rcluster
      write(uout,'(a,t35,a   )') 'weight                         = ', txweight

      if (l1d) then
         call DistFuncSample(iStage, nvar, var)
         call DistFuncHead(nvar, var, uout)
         var%nbin = nsizemax_save                ! truncate to nsizemax_save
         call DistFuncWrite('1d '//txheading, nvar, var, uout, ulist, ishow, iplot, ilist)
      end if

      if (l2d) then
         call DistFunc2DSample(iStage, nvar, vvar)
         call DistFunc2DHead(nvar, vvar, uout)
         vvar%nbin(1) = nsizemax2d_save(1)       !  truncate to nsizemax2d_save(1)
         vvar%nbin(2) = nsizemax2d_save(2)       !  truncate to nsizemax2d_save(2)
         if (ishow /= 0) call DistFunc2DShow(ishow, '2d '//txheading, nvar, vvar, uout)
         if (ilist /= 0) call DistFunc2DList(ilist, '2d '//txheading, nvar, vvar, ulist)
      end if

      if (lpercolation) then
         call ScalarSample(iStage, 1, npvar, pvar)
         call ScalarNorm(iStage, 1, npvar, pvar, 0)
         call ScalarWrite(iStage, 1, npvar, pvar, 1, '(a,t35,4f15.5,f15.0)', uout)
      end if

      deallocate(var, vvar)
      deallocate(txobjt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine WriteAverage
   integer(4) :: i, ibin
   real(8)    :: aver
   do i = 1, nvar
      aver = 0.0d0
      do ibin =-1, var(i)%nbin
         aver = aver + var(i)%avs2(ibin)*(var(i)%min+(ibin+0.5)*var(i)%bin)*var(i)%bin
      end do
      write(uout,'(a,i2,t27,a,t50,f15.5)') 'Average cluster size', i, var(i)%label, aver
   end do
end subroutine WriteAverage

!........................................................................

end subroutine ClusterSD

!************************************************************************
!> \page static static.F90
!! **ClusterSD1D**
!! *calculate cluster size distribution*
!************************************************************************


subroutine ClusterSD1D(mnobj, nobj, iclusteriobj, mnsizemax, ltest, unit, ncluster, nsizemax, nsizedist)

   implicit none
   character(40), parameter :: txroutine ='ClusterSD1D'

   integer(4), intent(in)  :: mnobj                ! maximal number of objects
   integer(4), intent(in)  :: nobj                 ! number of objects
   integer(4), intent(in)  :: iclusteriobj(mnobj)  ! pointer: object -> its cluster
   integer(4), intent(in)  :: mnsizemax            ! max size of largest cluster
   logical,    intent(in)  :: ltest                ! logical flag for test output
   integer(4), intent(in)  :: unit                 ! output unit
   integer(4), intent(out) :: ncluster             ! number of clusters
   integer(4), intent(out) :: nsizemax             ! largest cluster size
   integer(4), intent(out) :: nsizedist(mnsizemax) ! number of clusters of a given size (given by the argument)

   integer(4) :: nsize(mnobj)                      ! sizes of the clusters
   integer(4) :: iobj, icluster

! ... get sizes of the nobj clusters (including the empty clusters)

   nsize(1:nobj) = 0
   do iobj = 1, nobj
      nsize(iclusteriobj(iobj)) = nsize(iclusteriobj(iobj)) + 1
   end do

   if (ltest) write(unit,'(a,(20i4))') 'nsize(icluster) = ', nsize(1:nobj)

! ... get number of non-empty clusters

   ncluster = count(nsize(1:nobj) > 0)

   if (ltest) write(unit,'(a,20i4)') 'ncluster = ', ncluster

! ... get largest cluster size

   nsizemax = maxval(nsize(1:nobj))

   if (ltest) write(unit,'(a,20i4)') 'nsizemax = ', nsizemax
   if (nsizemax > mnsizemax) call Stop(txroutine, 'nsizemax > mnsizemax', unit)

! ... get cluster size distribution of the ncluster non-empty clusters

   nsizedist = 0
   do icluster = 1, nobj
      if (nsize(icluster) > 0) nsizedist(nsize(icluster)) = nsizedist(nsize(icluster)) + 1
   end do

   if (ltest) write(unit,'(a,(20i4))') 'nsizedist = ', nsizedist(1:nsizemax)

end subroutine ClusterSD1D

!************************************************************************
!> \page static static.F90
!! **ClusterSD2D**
!! *calculate 2d cluster size distribution*
!************************************************************************


subroutine ClusterSD2D(mnobj, nobj, iclusteriobj, iobjprop, mnsizemax, ltest, unit, nsizemax2d, nsizedist2d)

   implicit none
   character(40), parameter :: txroutine ='ClusterSD2D'
   integer(4), intent(in)  :: mnobj                ! maximal number of objects
   integer(4), intent(in)  :: nobj                 ! number of objects
   integer(4), intent(in)  :: iclusteriobj(mnobj)  ! pointer: object -> its cluster
   integer(4), intent(in)  :: iobjprop(mnobj)      ! property of the objects (should be 1 or 2) for 2d distribution
   integer(4), intent(in)  :: mnsizemax            ! max size of largest cluster
   logical,    intent(in)  :: ltest                ! logical flag for test output
   integer(4), intent(in)  :: unit                 ! output unit
   integer(4), intent(out) :: nsizemax2d(2)        ! largest cluster size in the two dimensions
   integer(4), intent(out) :: nsizedist2d(0:mnsizemax,0:mnsizemax)  ! number of clusters with different composition

   integer(4) :: nobjpropcluster(2,mnobj)          ! number of objects of differernt property of clusters
   integer(4) :: iobj, iprop, icluster, nobj1, nobj2, isize

! ... check that iobjprop is valid

   if (ltest) write(unit,'(a,(20i4))')'iobjprop = ', iobjprop(1:nobj)
   if (minval(iobjprop(1:nobj)) < 0) call Stop(txroutine, 'minval(iobjprop(1:nobj)) < 0)', unit)
   if (maxval(iobjprop(1:nobj)) > 2) call Stop(txroutine, 'maxval(iobjprop(1:nobj)) > 2)', unit)

! ... get the composition of each cluster (including empty clusters)

   nobjpropcluster = 0
   do iobj = 1, nobj
      iprop = iobjprop(iobj)
      icluster = iclusteriobj(iobj)
      nobjpropcluster(iprop,icluster) = nobjpropcluster(iprop,icluster) + 1
   end do
   if (ltest) write(unit,'(a,(20i4))')'nobjpropcluster(1,:) = ', nobjpropcluster(1,1:nobj)
   if (ltest) write(unit,'(a,(20i4))')'nobjpropcluster(2,:) = ', nobjpropcluster(2,1:nobj)

! ... determine largest numbers of the different object types in a cluster

   nsizemax2d(1) = maxval(nobjpropcluster(1,1:mnobj))
   nsizemax2d(2) = maxval(nobjpropcluster(2,1:mnobj))
   if (nsizemax2d(1) > mnsizemax) call Stop(txroutine, 'nsizemax2d(1) > mnsizemax', unit)
   if (nsizemax2d(2) > mnsizemax) call Stop(txroutine, 'nsizemax2d(2) > mnsizemax', unit)
   if (ltest) write(unit,'(a,(20i4))')'nsizemax2d(1:2) = ', nsizemax2d(1:2)

! ... determine number of clusters with different composition (including empty clusters)

   nsizedist2d(0:mnsizemax,0:mnsizemax) = 0
   do icluster = 1, nobj
      nobj1 = nobjpropcluster(1,icluster)
      nobj2 = nobjpropcluster(2,icluster)
      nsizedist2d(nobj1,nobj2) = nsizedist2d(nobj1,nobj2) + 1
   end do
   nsizedist2d(0,0) = 0   ! remove empty clusters

   if (ltest) then
      write(unit,*) 'nsizedist2d(0:nsizemax2d(1),0:nsizemax2d(2))'
      do isize = 0, nsizemax2d(1)
         write(unit,'((21i4))') nsizedist2d(isize,0:nsizemax2d(2))
      end do
   end if

end subroutine ClusterSD2D

!************************************************************************
!> \page static static.F90
!! **ZeroSecondMoment**
!! *calculate zero and second moment for an ionic, neutral system*
!************************************************************************


subroutine ZeroSecondMoment(iStage)

   use PotentialModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='ZeroSecondMoment'
   character(80), parameter :: txheading ='radial distribution function and zero and second moment (Coulomb system)'
   integer(4)   , parameter :: ntype = 1
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4),   allocatable,  save :: ipnt(:,:,:)
   real(8),                    save :: rcutdist, rcutdist2
   real(8),       allocatable, save :: vsum(:)
   real(8),       allocatable, save :: amcond(:,:)
   real(8),       allocatable, save :: ugr(:)
   real(8),       allocatable, save :: dens(:)
   integer(4) :: itype, ivar, ibin, ip, jp, ipt, jpt, iptjpt, ism, ibuf, m
   real(8)    :: dx, dy, dz, r2, usum, d, norm, InvFlt
   real(8)    :: summ, sums, sumz, dum, fac, facs, facz, dvol, volrdf

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   select case (iStage)
   case (iReadInput)

      vtype%l =.true.
      vtype%min = Zero
      vtype%max = 50.0d0
      vtype%nbin= 100
      rcutdist  = 50

      rcutdist2 = rcutdist**2
      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, ' vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['rdf']
      vtype%nvar = nptpt

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(npt,npt,ntype), vsum(nvar), amcond(-1:vtype(1)%nbin,npt+1), ugr(0:nptpt), dens(npt))
      ipnt = 0
      vsum = 0.0E+00
      amcond = 0.0E+00
      ugr = 0.0E+00
      dens = 0.0E+00

      nvar = 0
      do itype = 1, ntype
         do ipt  = 1, npt
            do jpt  = ipt, npt
               nvar = nvar + 1
               ipnt(ipt,jpt,itype) = nvar
               var(nvar)%label = trim(vtype(itype)%label)//' '//txpt(ipt)//' '//txpt(jpt)
               var(nvar)%min = vtype(itype)%min
               var(nvar)%max = vtype(itype)%max
               var(nvar)%nbin = vtype(itype)%nbin
            end do
         end do
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2+1
      itype = 1
      do ip = 1, np
         ipt = iptpn(ip)
         do jp = ip + 1, np
            jpt = iptpn(jp)
            ivar = ipnt(ipt,jpt,itype)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 < rcutdist2) then
               ibin = max(-1,min(floor(var(ivar)%bini*(sqrt(r2)-var(ivar)%min)),var(ivar)%nbin))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+One
            end if
         end do
      end do

   case (iAfterMacrostep)

! ... normalizing using the rdf volume

      itype = 1
      do ipt = 1, npt
         do jpt = ipt, npt
            ivar = ipnt(ipt,jpt,itype)
            vsum(ivar) = sum(var(ivar)%avs2(-1:var(ivar)%nbin))
            volrdf = min(FourPiThird*rcutdist**3,volst)
            norm = var(ivar)%nsamp2*volrdf/FourPiThird*InvFlt(vsum(ivar))
            if (ipt == jpt) norm = norm*(nppt(ipt)-1)/nppt(ipt)
            do ibin = 0, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end do
      end do
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)

! ... calculate energy and pressure from rdf

      ugr(0:nptpt) = Zero
      do ipt = 1, npt
         do jpt = ipt, npt
            iptjpt = iptpt(ipt,jpt)
            ivar = iptjpt
            fac = nppt(ipt)*nppt(jpt)*FourPiThird/volst
            if (ipt == jpt) fac = fac*0.5
            summ = Zero
            do ibin = 0, var(ivar)%nbin
               r2 = (var(ivar)%min+(0.5+ibin)*var(ivar)%bin)**2
               if (r2 < r2umin(iptjpt)) cycle
               ibuf = iubuflow(iptjpt)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
                  if (ibuf > nbuf) call Stop(txroutine, 'ibuf > nbuf', uout)
               end do
               d = r2-ubuf(ibuf)
               usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                                 d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
               summ = summ + usum*var(ivar)%avs1(ibin)*dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
            ugr(iptjpt) = ugr(iptjpt) + summ * fac
            ugr(0)      = ugr(0)      + summ * fac
         end do
      end do

! ............... write output  ..............

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,i5,5x,a,i5,2x,a)') 'number of grid points          = ', var(1)%nbin, '(',mnbin_DF,')'
      write(uout,'(a,t35,10f8.2)')          'cutoff distance                = ', rcutdist
      write(uout,'()')
      write(uout,'(a)') 'quantity (coulombic contribution from rdf''s)'
      write(uout,'(a)') '--------------------------------------------'
      write(uout,'(a,t20,f20.5)') 'total config. energy', ugr(0)/np
      write(uout,'(a,t20,f20.5)') (txptpt(iptjpt),ugr(iptjpt)/np,iptjpt = 1,nptpt)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

! ............. zero and second moment calcuations based on g(r) .............

!     linse and andersen, JCP, 85, 3027 (1986).

      ism = npt+1
      fac = FourPi*var(1)%bin
      dens(1:npt) = nppt(1:npt)/volst
      amcond(0:var(1)%nbin,1:ism) =-One
      do ipt = 1, npt
         do jpt = 1, npt
            iptjpt = iptpt(ipt,jpt)
            ivar = iptjpt
            facz =-fac*ucoffx(1,jpt)/ucoffx(1,ipt)*dens(jpt)
            facs =-fac*4.0*Pi/(6.0*GasConstant*temp*scltem)*ucoffx(1,iptjpt)*sclene*dens(ipt)*dens(jpt)
            sumz = Zero
            sums = Zero
            do ibin = 0, var(1)%nbin
               r2 = (var(ivar)%min+(0.5+ibin)*var(ivar)%bin)**2
               sumz = sumz+var(ivar)%avs1(ibin)*r2
               sums = sums+var(ivar)%avs1(ibin)*r2**2
               amcond(ibin,ipt) = amcond(ibin,ipt)+facz*sumz
               amcond(ibin,ism) = amcond(ibin,ism)+facs*sums
            end do
         end do
      end do

      if (ishow/= 0) then
         write(uout,'()')
         write(uout,'(a)') 'zero and second moment conditions'
         write(uout,'(a)') '---------------------------------'
         write(uout,'(5x,a,t20,a,t43,a)') 'distance', 'zero ...', 'second'
         write(uout,'(5x,a,t20,a,t43,a)') '--------', '--------', '------'
         do ibin = 0, var(1)%nbin
            r2 = (var(1)%min+(0.5+ibin)*var(1)%bin)
            write(uout,'(6f12.3)') r2, (amcond(ibin,ipt),ipt = 1,npt+1)
         end do
      end if
      if (iplot/= 0) then
         do m = 1, npt
            call Plot(txpt(m), var(1)%nbin, amcond(0,m), 'none', var(1)%min, var(1)%max, dum, dum, uout)
         end do
         call Plot('sm', var(1)%nbin, amcond(0,ism), 'none', var(1)%min, var(1)%max, dum, dum, uout)
      end if

      deallocate(var, ipnt, vsum, amcond, ugr, dens)

   end select

end subroutine ZeroSecondMoment

!************************************************************************
!> \page static static.F90
!! **MultipoleDF**
!! *calculate multipole moment distribution function and mean square average*
!************************************************************************


!> \page nmlMultipoleDF
!! The namelist  \ref nmlMultipoleDF contains variables that control the calculation of electrostatic multipole moment distribution functions.
!! * Variables:
!!  * \subpage lmax
!!  * \subpage vmin
!!  * \subpage vmax
!!  * \subpage nmlMultipoleDF_nbin
!!  * \subpage nmlMultipoleDF_nrad
!!  * \subpage nmlMultipoleDF_radius

!> \page lmax
!! `integer`
!! **default:** `2`
!! * The largest multipole moment considered is 2**\ref lmax Real and imaginary components are analyzed separately and only for nonnegative m.

!> \page vmin
!! `real`(1:\ref lmax)
!! **default:** `0.0`
!! * Lower end of the sampled distribution function of each type.

!> \page vmax
!! `real`(1:\ref lmax)
!! **default:** `10.0`
!! * Upper end of the sampled distribution function of each type.

!> \page nmlMultipoleDF_nbin nbin
!! `integer`
!! **default:** `100`
!! * Number of bins used to sample the distribution functions.

!> \page nmlMultipoleDF_nrad nrad
!! `integer`
!! **default:** `1`
!! * Number of radii for which multipole moment df are calculated.  If \ref txbc='sph', then nrad=1 and radius=rsph are forced; otherwise multipole moment are calculated over spheres centered on all particles.

!> \page nmlMultipoleDF_radius radius
!! `real`(1:nrad)
!! **default:** nrad*`20.0`
!! * Radii.

subroutine MultipoleDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MultipolDF'
   character(80), parameter :: txheading1 ='multipole moment distribution function'
   character(80), parameter :: txheading2 ='mean square multipole moment'
   character(80), parameter :: txheading3 ='mean square multipole moment (m averaged)'
   integer(4)   , parameter :: mlmax = 4               ! maximal number of lmax
   integer(4)   , parameter :: mnrad = 3               ! maximal number of nrad

   integer(4), save      :: lmax                       ! largest multipole moment index of l
   real(8),    save      :: vmin(mlmax)                ! lower distance limit of histogram
   real(8),    save      :: vmax(mlmax)                ! upper distance limit of histogram
   integer(4), save      :: nbin                       ! number of bins
   integer(4), save      :: nrad                       ! number of radii to be considered
   real(8),    save      :: radius(mnrad)              ! radius of the volume
   integer(4), save      :: iplow                      ! lower particle number
   integer(4), save      :: ipupp                      ! upper particle number
   integer(4), save      :: nvar, nlvar                ! number of variables
   integer(4), save      :: mfac(0:mlmax)              ! 1 for m = 0, 2 for other m
   integer(4), allocatable, save :: index(:,:,:)       ! pointer
   integer(4), allocatable, save :: lindex(:,:)        ! pointer
   complex(8), allocatable, save :: mpm(:,:,:)         ! multipole moment of the volume
   character(2)          :: txSph(0:1) = [' R', ' I']
   character(2)          :: str1, str2
   character(6)          :: str3
   integer(4)            :: l, m, s, ibin, ip, ivar, ilvar, irad
   real(8)               :: origin(1:3), value(0:1)

   type(df_var),      allocatable, save :: var(:)      ! distribution function
   type(scalar_var),  allocatable, save :: svar(:)     ! mean square average over configurations
   type(scalar_var),  allocatable, save :: lvar(:)     ! mean square average over configurations and m

   namelist /nmlMultipoleDF/ lmax, vmin, vmax, nbin, nrad, radius

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      lmax          = 2
      vmin(1:mlmax) = Zero
      vmax(1:mlmax) = 10.0d0
      nbin          = 100
      nrad          = 1
      radius(mnrad) = 20.0d0

      rewind(uin)
      read(uin,nmlMultipoleDF)

      if (lmax > mlmax) call Stop(txroutine, 'lmax > mlmax', uout)
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)
      if (nrad > mnrad) call Stop(txroutine, 'nrad > mnrad', uout)

   case (iWriteInput)

      if (txbc == 'xyz' .or. lbcrd .or. lbcto) then
         iplow = ipmyid(1)
         ipupp = ipmyid(2)
         iplow = 1
         ipupp = np
      else if (lbcsph) then
         iplow = 1
         ipupp = 1
         nrad = 1
         radius(1) = sphrad
      else
         call Stop(txroutine, 'error in txbc', uout)
      end if

! ... set nvar and nlvar

      nvar = (lmax+1)**2*mnrad
      nlvar = (lmax+1)*mnrad

! ... allocate memory

      allocate(var(nvar),svar(nvar),lvar(nvar))
      allocate(index(0:2*lmax+1,0:lmax,nrad), lindex(0:lmax,nrad))
      index = 0
      lindex = 0
      allocate(mpm(0:lmax,0:lmax,nrad))
      mpm = cmplx(Zero,Zero)

      ivar = 0
      ilvar = 0
      do irad = 1, nrad
         write(str3,'(f6.1)') radius(irad)
         do l = 1, lmax
            write(str1,'(i2)') l
            ilvar = ilvar + 1
            lindex(l, irad) = ilvar
            lvar(ilvar)%label = 'rad = '//trim(adjustl(str3))//' l = '//trim(adjustl(str1))
            lvar(ilvar)%norm  = One
            do m = 0, l
               write(str2,'(i2)') m
               do s = 0, 1
                  if ((m == 0) .and. (s == 1)) cycle
                  ivar = ivar + 1
                  index(2*m+s, l, irad) = ivar
                  var(ivar)%label  = 'rad = '//trim(adjustl(str3))//' l = '//trim(adjustl(str1))//' m = '//trim(adjustl(str2))//txSph(s)
                  var(ivar)%min    = vmin(l)
                  var(ivar)%max    = vmax(l)
                  var(ivar)%nbin   = nbin
                  svar(ivar)%label = var(ivar)%label
                  svar(ivar)%norm  = One
               end do
            end do
         end do
      end do
      nvar = ivar
      nlvar = ilvar

      mfac(0) = One           ! factor of two in averaging over nonaxial components
      mfac(1:lmax) = Two

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, svar)
      call ScalarSample(iStage, 1, nlvar, lvar)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, svar)
      call ScalarSample(iStage, 1, nlvar, lvar)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) svar
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) lvar

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, svar)
      call ScalarSample(iStage, 1, nlvar, lvar)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1
      do ip = iplow, ipupp                                            ! loop over particles
         if (lbcsph) then                                             ! spherical boundary contion
            origin(1:3) = [Zero, Zero, Zero]
         else                                                         ! periodic boundary conditions
            origin(1:3) = ro(1:3,ip)
         end if
         call CalcMultipole(lmax, mnrad, lmax, nrad, radius, origin, mpm)
         do irad = 1, nrad                                            ! loop over radii
            do l = 1, lmax                                            ! loop over multipole moment index l
               ilvar = lindex(l, irad)
               lvar(ilvar)%value = Zero
               do m = 0, l                                            ! loop over multipole moment index m
                  value(0) = real(mpm(m,l,irad))
                  value(1) = aimag(mpm(m,l,irad))
                  do s = 0, 1                                         ! assign real and imaginary values
                     if ((m == 0) .and. (s == 1)) cycle               ! no imaginary value for m==0
                     ivar = index(2*m+s, l, irad)
                     ibin = max(-1,min(floor(var(ivar)%bini*(value(s)-var(ivar)%min)),int(var(ivar)%nbin)))
                     var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
                     svar(ivar)%value = value(s)**2
                     call ScalarSample(iStage, ivar, ivar, svar)
                     lvar(ilvar)%value = lvar(ilvar)%value + mfac(m)*value(s)**2
                  end do
               end do
               call ScalarSample(iStage, ilvar, ilvar, lvar)
            end do
         end do
      end do

   case (iAfterMacrostep)

#if defined (_PAR_)
!       if (txbc == 'xyz' .or. lbcrd .or. lbcto) call DistFuncAllreduce(iStage, 1, nvar, var, vaux) ! OK
!       if (txbc == 'xyz' .or. lbcrd .or. lbcto) call ScalarAllreduce(iStage, 1, nvar, svar, vaux)  ! OK
!       if (txbc == 'xyz' .or. lbcrd .or. lbcto) call ScalarAllreduce(iStage, 1, nlvar, lvar, vaux) ! OK
#endif
      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      call ScalarSample(iStage, 1, nvar, svar)
      call ScalarSample(iStage, 1, nlvar, lvar)
      if (lsim .and. master) write(ucnf) var
      if (lsim .and. master) write(ucnf) svar
      if (lsim .and. master) write(ucnf) lvar
      call ScalarNorm(iStage, 1, nvar, svar, 0)
      call ScalarNorm(iStage, 1, nlvar, lvar, 0)

   case (iAfterSimulation)

      if (master) then
         call WriteHead(2, txheading1, uout)
         write(uout,'(a,i8    )') 'lmax                           = ', lmax
         write(uout,'(a,4es8.1)') 'vmin                           = ', vmin(1:mlmax)
         write(uout,'(a,4es8.1)') 'vmax                           = ', vmax(1:mlmax)
         write(uout,'(a,i8    )') 'nbin                           = ', nbin
         write(uout,'(a,i8    )') 'nrad                           = ', nrad
         write(uout,'(a,5f8.1 )') 'radius                         = ', radius(1:nrad)

#if defined (_PAR_)
!        if (txbc == 'xyz' .or. lbcrd .or. lbcto) call ScalarAllreduce(iStage, 1, nvar, svar, vaux)   ! problem with statistics.F90 routine
!        if (txbc == 'xyz' .or. lbcrd .or. lbcto) call ScalarAllreduce(iStage, 1, nlvar, lvar, vaux)  ! problem with statistics.F90 routine
#endif
         call DistFuncSample(iStage, nvar, var)
         call DistFuncHead(nvar, var, uout)
         call DistFuncWrite(txheading1, nvar, var, uout, ulist, ishow, iplot, ilist)

         call WriteHead(2, txheading2, uout)
         call ScalarSample(iStage, 1, nvar, svar)
         call ScalarNorm(iStage, 1, nvar, svar, 0)
         call ScalarWrite(iStage, 1, nvar, svar, 1, '(a,t35,4es15.4,f15.0)', uout)

         call WriteHead(2, txheading3, uout)
         call ScalarSample(iStage, 1, nlvar, lvar)
         call ScalarNorm(iStage, 1, nlvar, lvar, 0)
         call ScalarWrite(iStage, 1, nlvar, lvar, 1, '(a,t35,4es15.4,f15.0)', uout)

         deallocate(var, svar, lvar, index, lindex, mpm)
      end if

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine MultipoleDF

!************************************************************************
!> \page static static.F90
!! **CalcMultipole**
!! *calculate multipole moments of spherical volumes arising from dipoles*
!************************************************************************


subroutine CalcMultipole(mlmax, mnrad, lmax, nrad, radius, origin, mpm)

   use MolModule
   implicit none

   real(8), parameter :: qfac = sqrt(30d0), ofac = sqrt(105d0), hfac = sqrt(252d0)

   integer(4), intent(in)   :: mlmax
   integer(4), intent(in)   :: mnrad
   integer(4), intent(in)   :: lmax                       ! largest multipole moment index l
   integer(4), intent(in)   :: nrad                       ! number of radii to be considered
   real(8), intent(in)      :: radius(mnrad)              ! radius of the volume
   real(8), intent(in)      :: origin(1:3)                ! origin of the mulitpole moment expansion
   complex(8), intent(out)  :: mpm(0:mlmax,0:mlmax,mnrad) ! multipole moment of the volume

   real(8)    :: w3j(0:mlmax-1,0:mlmax,0:1,-(mlmax-1):mlmax-1,-mlmax:mlmax,-1:1)
   logical    :: first = .true.
   integer(4) :: jp, irad, m
   real(8)    :: fac, radius2(mnrad), dx, dy, dz, r2, r1, theta, phi
   complex(8) :: q(-1:1), CCLM

   if (first) then
      call SetW3J(mlmax-1,mlmax,1,w3j)
      first = .false.
   end if
   radius2(1:nrad) = radius(1:nrad)**2

   mpm(0:mlmax,0:mlmax,1:nrad) = cmplx(Zero,Zero)

   do jp = 1, np                                                 ! loop over particles
      dx = ro(1,jp)-origin(1)
      dy = ro(2,jp)-origin(2)
      dz = ro(3,jp)-origin(3)
      call PBCr2(dx,dy,dz,r2)
      do irad = 1, nrad                                           ! loop over radii
         if (r2 > radius2(irad)) cycle                            ! distance to particle too large
         call CarToStd1(dip(1,jp),dip(2,jp),dip(3,jp),q(-1:1))    ! dipole moment of particle jp in standard form
         if (lmax>=1) then                                        ! add dipole moment contribution
            mpm(0,1,irad) = mpm(0,1,irad) + q(0)
            mpm(1,1,irad) = mpm(1,1,irad) + q(1)
         endif
         if (r2 < 1.0d-40) cycle                                  ! remove r2 == zero
         call CarToSph('rad',dx,dy,dz,r1,theta,phi)
         if (lmax>=2) then                                        ! add quadrupole moment contribution
            fac = qfac*r1
            do m = -1, 1
                         mpm(0,2,irad) = mpm(0,2,irad) + cmplx(fac*w3j(1,2,1,m,   0,-m))*q(-m)*CCLM(1,m  ,theta,phi,0)
               if (m/=1)  mpm(1,2,irad) = mpm(1,2,irad) - cmplx(fac*w3j(1,2,1,m+1,-1,-m))*q(-m)*CCLM(1,m+1,theta,phi,0)
               if (m==-1) mpm(2,2,irad) = mpm(2,2,irad) + cmplx(fac*w3j(1,2,1,m+2,-2,-m))*q(-m)*CCLM(1,m+2,theta,phi,0)
            end do
         end if
         if (lmax>=3) then                                        ! add octupole moment contribution
            fac = ofac*r2
            do m = -1, 1
                         mpm(0,3,irad) = mpm(0,3,irad) - cmplx(fac*w3j(2,3,1,m   ,0,-m))*q(-m)*CCLM(2,m  ,theta,phi,0)
                         mpm(1,3,irad) = mpm(1,3,irad) + cmplx(fac*w3j(2,3,1,m+1,-1,-m))*q(-m)*CCLM(2,m+1,theta,phi,0)
               if (m/=1)  mpm(2,3,irad) = mpm(2,3,irad) - cmplx(fac*w3j(2,3,1,m+2,-2,-m))*q(-m)*CCLM(2,m+2,theta,phi,0)
               if (m==-1) mpm(3,3,irad) = mpm(3,3,irad) + cmplx(fac*w3j(2,3,1,m+3,-3,-m))*q(-m)*CCLM(2,m+3,theta,phi,0)
            end do
         end if
         if (lmax>=4) then                                        ! add hexadecapole moment contribution
            fac = hfac*r2*r1
            do m = -1, 1
                         mpm(0,4,irad) = mpm(0,4,irad) + cmplx(fac*w3j(3,4,1,m  , 0,-m))*q(-m)*CCLM(3,m  ,theta,phi,0)
                         mpm(1,4,irad) = mpm(1,4,irad) - cmplx(fac*w3j(3,4,1,m+1,-1,-m))*q(-m)*CCLM(3,m+1,theta,phi,0)
                         mpm(2,4,irad) = mpm(2,4,irad) + cmplx(fac*w3j(3,4,1,m+2,-2,-m))*q(-m)*CCLM(3,m+2,theta,phi,0)
               if (m/=1)  mpm(3,4,irad) = mpm(3,4,irad) - cmplx(fac*w3j(3,4,1,m+3,-3,-m))*q(-m)*CCLM(3,m+3,theta,phi,0)
               if (m==-1) mpm(4,4,irad) = mpm(4,4,irad) + cmplx(fac*w3j(3,4,1,m+4,-4,-m))*q(-m)*CCLM(3,m+4,theta,phi,0)
            end do
         end if
      end do
   end do

end subroutine CalcMultipole

!************************************************************************
!> \page static static.F90
!! **EnergyDF**
!! *calculate energy distribution functions*
!************************************************************************


!     group division of particles is required

!> \page nmlEnergyDF
!! The namelist  \ref nmlEnergyDF contains variables that control the calculation of energy distribution functions. Any combination of
!! the types of distribution functions listed below may be selected through vtype. Only contributions from two-body potential terms
!! are included.
!!
!!    | type | label | quantity                                            |
!!    | ---- | ----- | --------------------------------------------------- |
!!    | 1    | bindu | binding energy, total                               |
!!    | 2    | bindu | binding energy with particles in group jgr          |
!!    | 3    | pairu | pair energy with particles in group jgr within rmax |
!! * Variables:
!!  * \subpage nmlEnergyDF_vtype
!!  * \subpage nmlEnergyDF_rmax

!> \page nmlEnergyDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /-150.0,-150.0,-25.0/ Max: /0.0,0.0,0.0/

!> \page nmlEnergyDF_rmax rmax
!! `real`
!! **default:** `1.0d10`
!! * Upper distance for particle separation for pair energy distribution function.


subroutine EnergyDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='EnergyDF'
   character(80), parameter :: txheading ='energy distribution function'
   integer(4)   , parameter :: ntype = 3
   type(static1D_var),         save :: vtype(ntype)
   real(8),                    save :: rmax
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   real(8),       allocatable, save :: uigr(:,:), ui(:)
   integer(4) :: itype, ivar, ibin, ip, jp, igr, jgr, igrjgr
   real(8)    :: dx, dy, dz, r2, uuu, fdum(3)

   namelist /nmlEnergyDF/ vtype, rmax

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l   =.false.
      vtype%min = [ -100.d0, -100d0, -25.d0 ]
      vtype%max = [  Zero  ,  Zero ,  25.d0 ]
      vtype%nbin = 100
      rmax = 1.0d10

      rewind(uin)
      read(uin,nmlEnergyDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['bindu','bindu','pairu']
      vtype%nvar = [ ngr(1), ngrgr, ngrgr]

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(ngrgr,ntype), uigr(np_alloc,ngr(2)), ui(np_alloc))
      ipnt = 0
      uigr = 0.0E+00
      ui = 0.0E+00

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            if (itype == 1) then
               do igr = 1, ngr(1)
                  ivar = ivar+1
                  ipnt(igr,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            else if (itype == 2) then
               do igrjgr = 1, ngrgr
                  ivar = ivar+1
                  ipnt(igrjgr,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txgrgr(igrjgr)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            else if (itype == 3) then
               do igrjgr = 1, ngrgr
                  ivar = ivar+1
                  ipnt(igrjgr,itype) = ivar
                  var(ivar)%label = trim(vtype(itype)%label)//' '//txgrgr(igrjgr)
                  var(ivar)%min = vtype(itype)%min
                  var(ivar)%max = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
               end do
            end if
         end if
      end do
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1
      ui(1:np) = Zero
      uigr(1:np,1:ngr(2)) = Zero

      do ip = 1, np
         igr = igrpn(ip,1)
         if (igr <= 0) cycle
         do jp = 1, np
            if (ip == jp) cycle
            jgr = igrpn(jp,2)
            if (vtype(1)%l .or. jgr >= 1) call UTwoBodyPair(ip, jp, uuu, fdum)
            if (vtype(1)%l) ui(ip) = ui(ip)+uuu
            if (jgr <= 0) cycle
            if (vtype(2)%l) uigr(ip,jgr) = uigr(ip,jgr)+uuu

! ... sample type 3

            itype = 3
            if (vtype(itype)%l) then
               dx = ro(1,ip)-ro(1,jp)
               dy = ro(2,ip)-ro(2,jp)
               dz = ro(3,ip)-ro(3,jp)
               call PBCr2(dx,dy,dz,r2)
               if (r2 < rmax**2) then
                 ivar = ipnt(igrgr(igr,jgr),itype)
                 ibin = max(-1,min(floor(var(ivar)%bini*(uuu-var(ivar)%min)),int(var(ivar)%nbin)))
                 var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               end if
            end if

         end do
      end do

! ... sample type 1

      itype = 1
      if (vtype(itype)%l) then
         do ip = 1, np
            igr = igrpn(ip,1)
            if (igr <= 0) cycle
            ivar = ipnt(igr,1)
            ibin = max(-1,min(floor(var(ivar)%bini*(ui(ip)-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         end do
      end if

! ... sample type 2

      itype = 2
      if (vtype(itype)%l) then
         do ip = 1, np
            igr = igrpn(ip,1)
            if (igr <= 0) cycle
            do jgr = 1, ngr(2)
               ivar = ipnt(igrgr(igr,jgr),itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(uigr(ip,jgr)-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end do
         end do
      end if

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,f8.2)') 'upper separation for pair      = ', rmax
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, ipnt, uigr, ui)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine EnergyDF

!************************************************************************
!> \page static static.F90
!! **Widom1**
!! *calculate excess chemical potential using Widom's method (neutral set)*
!************************************************************************

!     see JCP 39, 2808 (1963) and Mol. Phys. 64, 247 (1988)

!     System may contain polyatom particles, but routine is only checked for insertion of monoatom particles

!> \page nmlWidom1
!! The namelist  \ref nmlWidom1 contains variables that control the calculation of excess chemical potentials using Widom's insertion method. Only charge neutral combinations of particles should be inserted.
!! * Variables:
!!  * \subpage ntimes
!!  * \subpage nset
!!  * \subpage nptset
!!  * \subpage iptset

!> \page ntimes
!! `integer`
!! **default:** `1`
!! * Number of samplings per occasion.

!> \page nset
!! `integer`
!! **default:** `1`
!! * Number of sets to evaluate.

!> \page nptset
!! `integer`(1:\ref nset)
!! * Number of particles to be inserted in a set.

!> \page iptset
!! `integer`(1:\ref nptset,1:\ref nset)
!! * Types of particles to be inserted in a sets
   subroutine Widom1(iStage)

   use MolModule
   use MollibModule, only: InvInt
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Widom1'
   character(80), parameter :: txheading ='excess chemical potential: Widom''s method, neutral set'
   integer(4)   , parameter :: mninset = 2            ! max number of sets
   integer(4)   , parameter :: mnptset = 3            ! max number of particles in each set
   integer(4),       save :: ntimes                   !*number of insertions per occasion
   integer(4),       save :: nset                     !*number of sets to evaluate
   integer(4),       save :: nptset(mninset)          !*number of particles to insert in each set
   integer(4),       save :: iptset(mnptset,mninset)  !*particle types to be inserted in the sets
   type(scalar_var), allocatable, save :: var(:)
   integer(4),       save :: nvar, ivol, iacc, iide, iexe
   integer(4) :: iset, itimes, iseedWidom
   real(8)    :: duwidom(0:npt)
   logical    :: lhsoverlap

   namelist /nmlWidom1/ ntimes, nset, nptset, iptset

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      ntimes = 1
      nset   = 1

      rewind(uin)
      read(uin,nmlWidom1)

!      if (lnvt .or. lntp) then
!      else
!         call Stop(txroutine, 'unsupported ensemble', uout)
!      end if
!     if (lpolyatom) call Stop(txroutine, 'lpolyatom not allowed', uout)
      if (nset > mninset) call Stop(txroutine, 'nset > mninset', uout)
      if (count(nptset > mnptset) > 0) call Stop(txroutine, 'nptset > mnptset', uout)

   case (iWriteInput)

!! ... set nvar as well as allocate memory

      ivol = 1
      iacc = ivol + 1
      iide = iacc + nset
      iexe = iide + nset
      nvar = iexe + nset - 1
      allocate(var(nvar))

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      iseedWidom = iseed        ! brach off seed (to preserve seed among master and slaves)

      var(ivol)%value = vol
      call ScalarSample(iStage, ivol, ivol, var)
      do itimes = 1, ntimes
         do iset = 1, nset
            var(iide+iset-1)%value = log(product(nppt(iptset(1:nptset(iset),iset))/vol))  ! ideal contribution
            call DUWidomx(nptset(iset), iptset(1,iset), iseedWidom, duwidom, lhsoverlap)  ! excess contribution
            if (.not.lhsoverlap) then
               var(iacc+iset-1)%value = One
               var(iexe+iset-1)%value = exp(-beta*duwidom(0))
            else
               var(iacc+iset-1)%value = Zero
               var(iexe+iset-1)%value = Zero
            end if
         end do
         var(iexe:iexe+nset-1)%value = vol*var(iexe:iexe+nset-1)%value
         call ScalarSample(iStage, iacc, nvar, var)
      end do

   case (iAfterMacrostep)

      var%nsamp1 = var%nsamp1 + 1
      var(ivol)%avs2 =  var(ivol)%avs2/var(ivol)%nsamp2
      var(iacc:iacc+nset-1)%avs2 = var(iacc:iacc+nset-1)%avs2/(var(iacc:iacc+nset-1)%nsamp2)
      var(iide:iide+nset-1)%avs2 = var(iide:iide+nset-1)%avs2/(var(iacc:iacc+nset-1)%nsamp2)
      var(iexe:iexe+nset-1)%avs2 = var(iexe:iexe+nset-1)%avs2/var(ivol)%avs2
      var(iexe:iexe+nset-1)%avs2 =-log(var(iexe:iexe+nset-1)%avs2*InvInt(var(iexe:iexe+nset-1)%nsamp2))
      var%avs1 = var%avs1 + var%avs2
      var%avsd = var%avsd + var%avs2**2
      if (lsim .and. master) write(ucnf) var

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,es10.3)')      'average volume                 = ', var(ivol)%avs2
      write(uout,'(a,t35,i10)')         'number of insertions per conf. = ', ntimes
      write(uout,'(a,t35,i10)')         'total number of insertions     = ', var(iacc)%nsamp2
      write(uout,'()')
      write(uout,'(a,t35,8i10)')        'particle type(s)               = ', (iptset(1:nptset(iset),iset), iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'fraction accepted              = ', (var(iacc+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., ideal (kT)      = ', (var(iide+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., excess (kT)     = ', (var(iexe+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., total (kT)      = ', (var(iide+iset-1)%avs2+var(iexe+iset-1)%avs2, iset = 1, nset)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,es10.3)')       'average volume                 = ', var(ivol)%avs2
      write(uout,'(a,t35,i10)')          'number of insertions per conf. = ', ntimes
      write(uout,'(a,t35,i10)')          'total number of insertions     = ', var(iacc)%nsamp1*var(iacc)%nsamp2
      write(uout,'()')
      write(uout,'(a,t35,8i10)')         'particle type(s)               = ', (iptset(1:nptset(iset),iset), iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'fraction accepted              = ', (var(iacc+iset-1)%avs1, var(iacc+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., ideal (kT)      = ', (var(iide+iset-1)%avs1, var(iide+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., excess (kT)     = ', (var(iexe+iset-1)%avs1, var(iexe+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., total (kT)      = ', (var(iide+iset-1)%avs1+var(iexe+iset-1)%avs1, &
                                    sqrt(var(iide+iset-1)%avsd**2+var(iexe+iset-1)%avsd**2), iset = 1, nset)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine Widom1

!************************************************************************
!> \page static static.F90
!! **Widom2**
!! *calculate excess chemical potential using Widom's method (charge integration)*
!************************************************************************

!     see JCP 39, 2808 (1963) and Mol. Phys. 64, 247 (1988)

! ... not generalized to npt ensemble

!> \page nmlWidom2
!!The namelist  \ref nmlWidom2 contains variables that control the calculation of excess chemical potentials using Widom's insertion.
!! Single particles are inserted and charge integration according to Svensson and Jönsson.
!! * Variables:
!!  * \subpage nmlWidom2_ntimes
!!  * \subpage nmlWidom2_nset
!!  * \subpage nmlWidom2_iptset

!> \page nmlWidom2_ntimes ntimes
!! `integer`
!! **default:** `1`
!! * Number of samplings per occasion.

!> \page nmlWidom2_nset nset
!! `integer`
!! **default:** `1`
!! * Number of sets to evaluate.

!> \page nmlWidom2_iptset iptset
!! `integer`(1:\ref nptset,1:\ref nset)
!! * Types of particles to be inserted in a sets


   subroutine Widom2(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='Widom2'
   character(80), parameter :: txheading ='excess chemical potential: Widom''s method, charge integration'
   integer(4)   , parameter :: mninset     = 2          ! max number of sets
   integer(4),       save :: ntimes
   integer(4),       save :: nset
   integer(4),       save :: iptset(mninset)
   integer(4),       save :: nvar, ivol, iacc, iide, iexe
   integer(4),       save :: nlam
   type(scalar_var), allocatable, save :: var(:)
   real(8),          allocatable, save :: utop(:,:), ubot(:,:), uratio(:)
   real(8),          allocatable, save :: zati(:)
   real(8),          save :: dlam
   integer(4) :: itimes, iset, ilam, iseedWidom
   real(8)    :: lam, duel, duwidom(0:npt), npi, Trap, fac
   logical    :: lhsoverlap

   namelist /nmlWidom2/ ntimes, nset, iptset

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if (lnvt) then

   else if (lntp) then
      return           ! temporary exit
   else
      call Stop(txroutine, 'unsupported ensemble', uout)
   end if

   select case (iStage)
   case (iReadInput)

      ntimes = 1
      nset   = 1

      rewind(uin)
      read(uin,nmlWidom2)

      if (lpolyatom) call Stop(txroutine, 'lpolyatom not allowed', uout)
      if (nset > mninset) call Stop(txroutine, 'nset > mninset', uout)
      do iset = 1, nset
         if (iptset(iset) < 0 .or. iptset(iset) > npt) call Stop(txroutine, 'iptset out of range', uout)
      end do

   case (iWriteInput)

! ... set nlam, the number of integration points - 1, and allocate memory

      nlam = 5
      dlam = One/nlam
      allocate(utop(0:nlam,nset), ubot(0:nlam,nset), uratio(0:nlam), zati(nat))
      utop = 0.0E+00
      ubot = 0.0E+00
      uratio = 0.0E+00
      zati = 0.0E+00

! ... set nvar as well as allocate memory

      ivol = 1
      iacc = ivol + 1
      iide = iacc + nset
      iexe = iide + nset
      nvar = iexe + nset - 1                                       ! is correct
      allocate(var(nvar))

! ... check

      zati = Zero
      where(zat(1:npt) /= Zero) zati(1:npt) = One/zat(1:npt)
      if (count(zat(1:npt) /= Zero) == 0) call Stop(txroutine, 'all zat are Zero', uout)

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      utop(0:nlam,1:nset) = Zero
      ubot(0:nlam,1:nset) = Zero

   case (iSimulationStep)

      iseedWidom = iseed        ! brach off seed (to preserve seed among master and slaves)

      npi = One/np
      var(ivol)%value = vol
      call ScalarSample(iStage, ivol, ivol, var)
      do itimes = 1, ntimes
         do iset = 1, nset
            var(iide+iset-1)%value = log(nppt(iptset(iset))/vol)  ! ideal contribution
            call DUWidomx(1, iptset(iset), iseedWidom, duwidom, lhsoverlap)
            if (.not.lhsoverlap) then
               var(iacc+iset-1)%value = One
               do ilam = 0, nlam
                  lam = ilam*dlam
                  duel = beta*sum((One-lam*zat(iptset(iset))*zati(1:npt)*npi)*duwidom(1:npt))
                  fac = exp(-lam*duel)
                  utop(ilam,iset) = utop(ilam,iset) + duel*fac
                  ubot(ilam,iset) = ubot(ilam,iset) +      fac
!                  if (iset == 2 .and. ilam == 1) then
!                  write(*,*) 'itimes, iset, ilam, exp(), duel*exp(), utop, ubot', &
!                     itimes, iset, ilam, fac, duel*fac, utop(ilam,iset), ubot(ilam,iset)
!                  end if
               end do
            else
               var(iacc+iset-1)%value = Zero
            end if
         end do
         var(iexe:iexe+nset-1)%value = Zero     ! added below
 !     write(*,*) 'itimes 1, utop(1,2)', itimes, utop(1,2)
         call ScalarSample(iStage, 1, nvar, var)
 !     write(*,*) 'itimes 2, utop(1,2)', itimes, utop(1,2)
      end do
 !     write(*,*) '1 , utop(1,2)', utop(1,2)
   case (iAfterMacrostep)

      var%nsamp1 = var%nsamp1 + 1
      var(ivol)%avs2 = var(ivol)%avs2/var(ivol)%nsamp2
      var(iacc:iacc+nset-1)%avs2 = var(iacc:iacc+nset-1)%avs2/(var(iacc:iacc+nset-1)%nsamp2)
      var(iide:iide+nset-1)%avs2 = var(iide:iide+nset-1)%avs2/(var(iacc:iacc+nset-1)%nsamp2)
      do iset = 1, nset
         uratio(0:nlam) = utop(0:nlam,iset)/ubot(0:nlam,iset)
         var(iexe+iset-1)%avs2 = -log(var(iacc+iset-1)%avs2) + Trap(nlam+1,uratio(0),dlam)   ! mu(excess)
 !        write(*,*) 'iset, utop', iset, utop(0:nlam,iset)
 !        write(*,*) 'iset, ubot', iset, ubot(0:nlam,iset)
      end do
      var%avs1 = var%avs1 + var%avs2
      var%avsd = var%avsd + var%avs2**2
      if (lsim .and. master) write(ucnf) var

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,es10.3)')      'average volume                 = ', var(ivol)%avs2
      write(uout,'(a,t35,i10)')         'number of insertions per conf. = ', ntimes
      write(uout,'(a,t35,i10)')         'total number of insertions     = ', var(iacc)%nsamp2
      write(uout,'()')
      write(uout,'(a,t35,2(i10,5x))')   'particle type(s)               = ', (iptset(iset), iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'fraction accepted              = ', (var(iacc+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., ideal (kT)      = ', (var(iide+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., excess (kT)     = ', (var(iexe+iset-1)%avs2, iset = 1, nset)
      write(uout,'(a,t35,2(f10.3,5x))') 'chemical pot., total (kT)      = ', (var(iide+iset-1)%avs2+var(iexe+iset-1)%avs2, iset = 1, nset)

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,es10.3)')       'average volume                 = ', var(ivol)%avs1
      write(uout,'(a,t35,i10)')          'number of insertions per conf. = ', ntimes
      write(uout,'(a,t35,i10)')          'total number of insertions     = ', var(iacc)%nsamp1*var(iacc)%nsamp2
      write(uout,'()')
      write(uout,'(a,t35,2(i10,5x))')    'particle type(s)               = ', (iptset(iset), iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'fraction accepted              = ', (var(iacc+iset-1)%avs1, var(iacc+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., ideal (kT)      = ', (var(iide+iset-1)%avs1, var(iide+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., excess (kT)     = ', (var(iexe+iset-1)%avs1, var(iexe+iset-1)%avsd, iset = 1, nset)
      write(uout,'(a,t35,2(2f10.3,5x))') 'chemical pot., total (kT)      = ', (var(iide+iset-1)%avs1+var(iexe+iset-1)%avs1, &
                                    sqrt(var(iide+iset-1)%avsd**2+var(iexe+iset-1)%avsd**2), iset = 1, nset)

      deallocate(utop, ubot, uratio, var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine Widom2

!************************************************************************
!> \page static static.F90
!! **DUWidomx**
!! *calculate the change of the potential energy of the insertion*
!************************************************************************


subroutine DUWidomx(nptset, iptset, iseedWidom, duwidom, lhsoverlap)

   use MolModule
   implicit none

   integer(4), intent(in) :: nptset
   integer(4), intent(in) :: iptset(*)
   integer(4), intent(inout) :: iseedWidom
   real(8), intent(out)   :: duwidom(0:npt)
   logical, intent(out)   :: lhsoverlap

   integer(4) :: ip, jp, ipt, jpt
   integer(4) :: ja, iat, jat
   real(8) :: roset(3,3), uuu, dxpbc, dypbc, dzpbc
   real(8) :: Random

! ... get random positions

   if (lbcbox) then
      do ip = 1, nptset
         roset(1,ip) = boxlen(1)*(Random(iseedWidom)-Half)
         roset(2,ip) = boxlen(2)*(Random(iseedWidom)-Half)
         roset(3,ip) = boxlen(3)*(Random(iseedWidom)-Half)
      end do
   else if (lbcrd) then
      do ip = 1, nptset
         do
            roset(1,ip) = boxlen(1)*(Random(iseedWidom)-Half)
            roset(2,ip) = boxlen(2)*(Random(iseedWidom)-Half)
            roset(3,ip) = boxlen(3)*(Random(iseedWidom)-Half)
            call PBC2(roset(1,ip),roset(2,ip),roset(3,ip),dxpbc,dypbc,dzpbc)
            if ((dxpbc == Zero) .and. (dypbc == Zero) .and. (dzpbc == Zero)) exit
         end do
      end do
   else if (lbcto) then
      do ip = 1, nptset
         do
            roset(1,ip) = boxlen(1)*(Random(iseedWidom)-Half)
            roset(2,ip) = boxlen(2)*(Random(iseedWidom)-Half)
            roset(3,ip) = boxlen(3)*(Random(iseedWidom)-Half)
            call PBC2(roset(1,ip),roset(2,ip),roset(3,ip),dxpbc,dypbc,dzpbc)
            if ((dxpbc == Zero) .and. (dypbc == Zero) .and. (dzpbc == Zero)) exit
         end do
      end do
   else if (lbcsph) then
      do ip = 1, nptset
         do
            roset(1,ip) = Two*sphrad*(Random(iseedWidom)-Half)
            roset(2,ip) = Two*sphrad*(Random(iseedWidom)-Half)
            roset(3,ip) = Two*sphrad*(Random(iseedWidom)-Half)
            if (roset(1,ip)**2+roset(2,ip)**2+roset(3,ip)**2 < sphrad2) exit
         end do
      end do
   else if (lbccyl) then
      do ip = 1, nptset
         do
            roset(1,ip) = Two*cylrad*(Random(iseedWidom)-Half)
            roset(2,ip) = Two*cylrad*(Random(iseedWidom)-Half)
            if (roset(1,ip)**2+roset(2,ip)**2 < cylrad2) exit
         end do
         roset(3,ip) = cyllen*(Random(iseedWidom)-Half)
      end do
   end if

! ... calculate energy difference

   if (lpolyatom) then             ! polyatomic system BUT yet monoatomic inserted particle
      duwidom(0:npt) = Zero
      do ip = 1, nptset
         ipt = iptset(ip)
         iat = iatpt(ipt)
         do jp = 1, np                ! inserted particle - other particle pairs
            jpt = iptpn(jp)
            do ja = ianpn(jp), ianpn(jp)+napt(jpt)-1
               jat = iatan(ja)
               call UTwoBodyWidomP(iat, roset(1,ip), jat, r(1,ja), uuu, lhsoverlap)
               if (lhsoverlap) goto 400
               duwidom(jpt) = duwidom(jpt)+uuu
            end do
         end do

         do jp = ip+1, nptset          ! inserted particle - inserted particle pairs
            jpt = iptset(jp)
            jat = iatpt(jpt)
            call UTwoBodyWidomP(iat, roset(1,ip), jat, roset(1,jp), uuu, lhsoverlap)
            if (lhsoverlap) goto 400
            duwidom(jpt) = duwidom(jpt)+uuu
         end do
      end do
      duwidom(0) = sum(duwidom(1:npt))
   else
      duwidom(0:npt) = Zero
      do ip = 1, nptset
         ipt = iptset(ip)
         do jp = 1, np                 ! inserted particle - other particle pairs
            jpt = iptpn(jp)
            call UTwoBodyWidom(ipt, roset(1:3,ip), jpt, ro(1:3,jp), uuu, lhsoverlap)
            if (lhsoverlap) goto 400
            duwidom(jpt) = duwidom(jpt)+uuu
         end do
         do jp = ip+1, nptset          ! inserted particle - inserted particle pairs
            jpt = iptset(jp)
            call UTwoBodyWidom(ipt, roset(1:3,ip), jpt, roset(1:3,jp), uuu, lhsoverlap)
            if (lhsoverlap) goto 400
            duwidom(jpt) = duwidom(jpt)+uuu
         end do
      end do
      duwidom(0) = sum(duwidom(1:npt))
   end if

400 continue

end subroutine DUWidomx

!************************************************************************
!> \page static static.F90
!! **UTwoBodyWidom**
!! *calculate two-body potential between two particles; monoatomic particles*
!************************************************************************


subroutine UTwoBodyWidom(ipt, roi, jpt, roj, usum, lhsoverlap)

   use MolModule
   implicit none
   character(8), parameter :: txroutine ='UtwoBodyWidom'

   integer(4), intent(in) :: ipt
   real(8), intent(in)    :: roi(3)
   integer(4), intent(in) :: jpt
   real(8), intent(in)    :: roj(3)
   real(8), intent(out)   :: usum
   logical, intent(out)   :: lhsoverlap
   integer(4) :: iptjpt, ibuf
   real(8) :: dx, dy, dz, r2, d

   lhsoverlap =.true.
   usum = Zero

   iptjpt = iptpt(ipt,jpt)
   dx = roi(1)-roj(1)
   dy = roi(2)-roj(2)
   dz = roi(3)-roj(3)
   call PBCr2(dx,dy,dz,r2)
   if (r2 < r2atat(iptjpt)) goto 400
   if (r2 < r2umin(iptjpt)) call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)

   if (lewald) then  ! temporary fix for lewald (then mi convention for the chem potential)

      usum = EpsiFourPi*zat(ipt)*zat(jpt)/sqrt(r2)

   else             ! standard part

      if (r2 > rcut2) goto 300
      ibuf = iubuflow(iptjpt)
      do
         if (r2 >= ubuf(ibuf)) exit
         ibuf = ibuf+12
         if (ibuf > nbuf) call Stop(txroutine, 'ibuf > nbuf', uout)
      end do
      d = r2-ubuf(ibuf)
      usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                        d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

   end if

300 continue

   lhsoverlap =.false.

400 continue

end subroutine UTwoBodyWidom

!************************************************************************
!> \page static static.F90
!! **UTwoBodyWidomP**
!! *calculate two-body potential between two particles, polyatom particles*
!************************************************************************


subroutine UTwoBodyWidomP(iat, roi, jat, roj, usum, lhsoverlap)

   use MolModule
   implicit none
   character(8), parameter :: txroutine ='UTwoBodyWidomP'

   integer(4), intent(in) :: iat
   real(8), intent(in)    :: roi(3)
   integer(4), intent(in) :: jat
   real(8), intent(in)    :: roj(3)
   real(8), intent(out)   :: usum
   logical, intent(out)   :: lhsoverlap
   integer(4) :: iatjat, ibuf
   real(8) :: dx, dy, dz, r2, d

   lhsoverlap =.true.
   usum = Zero

   iatjat = iatat(iat,jat)
   dx = roi(1)-roj(1)
   dy = roi(2)-roj(2)
   dz = roi(3)-roj(3)
   call PBCr2(dx,dy,dz,r2)
   if (r2 < r2atat(iatjat)) goto 400
   if (r2 < r2umin(iatjat)) call Stop(txroutine, 'r2 < r2umin(iatjat)', uout)

   if (lewald) then  ! temporary fix for lewald (then mi convention for the chem potential)

      call Stop(txroutine, 'Not available for Ewald summation', uout)

   else             ! standard part

      if (r2 > rcut2) goto 300
      ibuf = iubuflow(iatjat)
      do
         if (r2 >= ubuf(ibuf)) exit
         ibuf = ibuf+12
         if (ibuf > nbuf) call Stop(txroutine, 'ibuf > nbuf', uout)
      end do
      d = r2-ubuf(ibuf)
      usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                        d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))

   end if

300 continue

   lhsoverlap =.false.

400 continue

end subroutine UTwoBodyWidomP

!************************************************************************
!> \page static static.F90
!! **MeanForce1**
!! *calculate mean energies and forces between two particles*
!************************************************************************


!     f = fmean + fhs
!     system: cylinder with particles 1 and 2 fixed, z(1) = -z(2) is required

!> \page nmlMeanForce1
!! The namelist  \ref nmlMeanForce1 contains variables that control the calculation of mean force between two particles (\ref txbc='cyl' is
!! required) according to F = Fmean + Fhs (the surface approach).
!! * Variables:
!!  * \subpage dr

!> \page dr
!! `real`
!! * Displacement for evaluation of Fhs

subroutine MeanForce1(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MeanForce1'
   character(80), parameter :: txheading ='mean energies and forces (f = fmean + fhs)'
   integer(4), save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   real(8),    save :: n12(3)                      ! normal vector from part 2 and 1
   real(8),    save :: dr                          ! displacement for evaluation of fhs (need to be optimized)
   real(8),    save :: dri, dx, dy, dz
   integer(4) :: ip, jp
   real(8) :: uuu, fff(3), usum,  fsum(3), fhs_1, fhs_2, result1, result2, result3, result4, data(8)

   namelist /nmlMeanForce1/ dr

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom)', uout)
      if (.not.lbccyl) call Stop(txroutine, '.not.lbccyl', uout)

      rewind(uin)
      read(uin,nmlMeanForce1)

      dr = abs(dr)          ! ensure positive dr

   case (iWriteInput)

! ... set nvar as well as allocate memory

      nvar = 11
      allocate(var(nvar))

! ... set label

      var(1)%label = 'umean/kT                       = '
      var(2)%label = 'fmean_x/kT                     = '
      var(3)%label = 'fmean_y/kT                     = '
      var(4)%label = 'fmean_z/kT                     = '
      var(5)%label = 'fmean_n/kT                     = '
      var(6)%label = 'fhs_n/kT                       = '
      var(7)%label = 'f/kT                           = '
      var(8)%label = 'aver number of overlap (fhs)   = '
      var(9)%label = 'aver number of overlap (fhs)   = '
      var(10)%label = 'aver number of overlap (fhs)   = '
      var(11)%label = 'aver number of overlap (fhs)   = '

! ... determine the normal between particle 1 and 2 in positive direction

      n12(1:3) = r(1:3,1)-r(1:3,2)
      uuu = One/sqrt( n12(1)**2 + n12(2)**2 + n12(3)**2)
      n12(1:3) = n12(1:3)*uuu

      dri = One/dr
      dx = n12(1)*dr
      dy = n12(2)*dr
      dz = n12(3)*dr

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

! ... sample the energy and force on particle 1 and 2

      usum = Zero
      fsum(1:3) = Zero
      do ip = 1, 2
         do jp = 1, np
            if (jp == ip) cycle
            call UTwoBodyPair(ip,jp,uuu,fff)      ! force on ip from jp
            usum = usum + uuu
            if (ip == 1) fsum(1:3) = fsum(1:3) + fff(1:3)
            if (ip == 2) fsum(1:3) = fsum(1:3) - fff(1:3)
         end do
      end do
      var(1)%value   = usum * (Half*beta)
      var(2:4)%value = fsum(1:3) * (Half*beta)
      var(5)%value = n12(1)*var(2)%value + n12(2)*var(3)%value + n12(3)*var(4)%value

! ... sample the force on particle 1 and 2 arising from hs overlap by displacing the particle
!     ref Wu et al. JCP 111, 7084 (1999)

      call hsoverlap_singlemove(1,+dx,+dy,+dz,result1)    ! move part 1 in positive n-dir
      call hsoverlap_singlemove(1,-dx,-dy,-dz,result2)    ! move part 1 in negative n-dir
      call hsoverlap_singlemove(2,+dx,+dy,+dz,result3)    ! move part 2 in positive n-dir
      call hsoverlap_singlemove(2,-dx,-dy,-dz,result4)    ! move part 2 in negative n-dir

      fhs_1 = -(result1-result2)                          ! fhs on part 1 in positive n-dir
      fhs_2 = -(result3-result4)                          ! fhs on part 2 in positive n-dir
      var(6)%value = Half*(fhs_1-fhs_2)*dri

      var(8)%value = result1
      var(9)%value = result2
      var(10)%value = result3
      var(11)%value = result4

! ... total mean force

      var(7)%value = var(5)%value + var(6)%value

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,3f10.4)') 'displacement parameter for fhs = ', dr
      write(uout,'(a,t35,3f10.4)') 'position of part 1             = ', r(1:3,1)
      write(uout,'(a,t35,3f10.4)') 'position of part 2             = ', r(1:3,2)
      write(uout,'(a,t35,3f10.4)') 'normal vector (from 2 to 1)    = ', n12(1:3)
      write(uout,'()')
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

      data(1)   = Two*r(3,2)                        ! separation
      data(2)   = dr                                ! displacement parameter for fhs
      data(3:4) = [ var(7)%avs1, var(7)%avsd ]      ! f/kT
      data(5:6) = [ var(5)%avs1, var(5)%avsd ]      ! fmean/kT
      data(7:8) = [ var(6)%avs1, var(6)%avsd ]      ! fhs/kT
      call WriteVecAppend(8, data(1:8), 'meanforce1.data')

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

subroutine hsoverlap_singlemove(ip, dx, dy, dz, result) ! check hs-overlap: ip is displaced (dx,dy,dz)
   use MolModule
   implicit none
   integer(4), intent(in)  :: ip
   real(8)   , intent(in)  :: dx, dy, dz
   real(8)   , intent(out) :: result
   integer(4) :: iploc
   integer(4) :: jpdum
   logical :: lhsoverlap

   nptm  = 1
   iploc = 1
   ipnptm(iploc) = ip
   rotm(1,iploc) = ro(1,ip) + dx
   rotm(2,iploc) = ro(2,ip) + dy
   rotm(3,iploc) = ro(3,ip) + dz
   call SetTrialAtomProp
   call UTwoBodyANew(lhsoverlap,jpdum)

   result = Zero
   if (lhsoverlap) result = One  ! hard-core overlap for particle ip

end subroutine hsoverlap_singlemove

!........................................................................

end subroutine MeanForce1

!************************************************************************
!> \page static static.F90
!! **MeanForce2**
!! *calculate mean energies and forces between two particles*
!************************************************************************


!     f = fcorr + fbond + fhs(midplane) + fideal(midplane and an end)
!     system: cylinder with particles 1 and 2 fixed, z(1) = -z(2) is required

!> \page nmlMeanForce2
!! The namelist  \ref nmlMeanForce2 contains variables that control the calculation of mean force between two particles (only \ref txbc='cyl')
!! according to F = Fcorr + Fideal(midplane and an end) + Fhs(midplane) (midplane approach).
!! * Variables:
!!  * \subpage thickness
!!  * \subpage dz

!> \page thickness
!! `real`
!! * Thickness of volume for sampling Fideal

!> \page dz
!! `real`
!! * Displacement for evaluation of fhs

subroutine MeanForce2(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='MeanForce2'
   character(80), parameter :: txheading ='mean energies and forces (f = fcorr + fideal(midplane and an end) + fhs(midplane))'
   integer(4), save :: nvar
   type(scalar_var), allocatable, save :: var(:)
   real(8),    save :: thickness                    ! thickness of volume for sampling fideal (need to be optimized)
   real(8),    save :: dz                           ! displacement for evaluation of fhs (need to be optimized)
   real(8),    save :: thickness2, thicknessi, dzi
   real(8),    save :: rmin, rmax, zlow, zupp

   integer(4) :: ip, jp, ipt, ia, iat
   real(8) :: uuu, fff(3), usum,  fsum(3), ubondacross, fbondacross(3), partsum, shs1, shs2
   real(8) :: data(12), huge

   namelist /nmlMeanForce2/ thickness, dz

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      if (.not.lmonoatom) call Warn(txroutine, '.not.lmonoatom', uout)
      if (.not.lbccyl) call Stop(txroutine, '.not.lbccyl', uout)

      rewind(uin)
      read(uin,nmlMeanForce2)

      thickness2 = Half*thickness
      thicknessi = One/thickness

      dzi = One/dz

   case (iWriteInput)

! ... set nvar as well as allocate memory

      nvar = 14
      allocate(var(nvar))

! ... set label

      var(1)%label = 'ucorr/kT                       = '
      var(2)%label = 'fcorr_x/kT                     = '
      var(3)%label = 'fcorr_y/kT                     = '
      var(4)%label = 'fcorr_z/kT                     = '
      var(5)%label = 'ubond/kT                       = '
      var(6)%label = 'fbond_x/kT                     = '
      var(7)%label = 'fbond_y/kT                     = '
      var(8)%label = 'fbond_z/kT                     = '
      var(9)%label = 'fhs/kT                         = '
      var(10)%label = 'area*density(z = 0)            = '
      var(11)%label = 'area*density(z = cyllen/2)       = '
      var(12)%label = 'area*diff(density)             = '
      var(13)%label = 'f/kT                           = '
      var(14)%label = 'aver number of overlap (fhs)   = '

! ... set zlow and zupp used to determin fhs (not completely generalized)

      rmin = huge(One)
      rmax = Zero
      do ip = 3, np
         ipt = iptpn(ip)
         if (napt(ipt) > 1) call Stop(txroutine, 'napt(ipt) > 1', uout)
         ia = ianpn(ip)
         iat = iatan(ia)
         rmin = min(rmin, radat(iat))
         rmax = max(rmax, radat(iat))
      end do
      if ((rmax > Zero) .and. (rmin > Zero)) then
         if (rmax/rmin - One > 0.001d0) call Stop(txroutine, 'rmin /= rmax', uout)
      end if
      zlow = -Two*rmax
      zupp = +Two*rmax

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

! ... sample the correlation energy and force across z = 0

      usum = Zero
      fsum(1:3) = Zero
      do ip = 1, np
         do jp = ip+1, np
            if (ro(3,ip)*ro(3,jp) < Zero) then ! ip and jp are residing in different Half
               call UTwoBodyPair(ip, jp, uuu, fff)
               usum = usum + uuu
               if (ro(3,ip) < Zero) fsum(1:3) = fsum(1:3) - fff(1:3)
               if (ro(3,ip) > Zero) fsum(1:3) = fsum(1:3) + fff(1:3)
            end if
         end do
      end do
      var(1)%value   = usum * beta
      var(2:4)%value = fsum(1:3) * beta

! ... sample the bond energy and force from bonds across z = 0

      if (lchain) then
         call UBondAcrossZ(ubondacross, fbondacross)
         var(5)%value   = ubondacross * beta
         var(6:8)%value = fbondacross(1:3) * beta
      else
         var(5)%value   = Zero
         var(6:8)%value = Zero
      end if

! ... sample the hard-core contribution across the midplane

      if (rmax > Zero) then
         shs1 = Zero
         shs2 = Zero
         do ip = 3, np
            if (ro(3,ip) < zlow) cycle
            if (ro(3,ip) > zupp) cycle
            if (ro(3,ip) < Zero) then
                if (MFHCOverlapSingleMove(ip, Zero, Zero, +dz, jp)) then  ! overlap with jp
                   if (ro(3,ip)*ro(3,jp) < Zero) then                     ! overlap across z = 0
                      var(14)%value = var(14)%value + One
                      shs1 = shs1 + One
!               write(*,'(a,2(i5,3f10.3,3x))') 'hs overlap: ip, jp:', ip, ro(1:3,ip), jp, ro(1:3,jp)
!               do jp = 3, np
!                  if (jp == ip) cycle
!                  r1 = sqrt(sum((ro(1:3,ip)-ro(1:3,jp))**2))
!                  if (r1 < 10*dz+2*radat(iatan(jp))) then
!                     write(*,'(5x,a,i5,a,3f10.3)') 'jp = ', jp, ' ro = ', ro(1:3,jp)
!                  end if
!               end do
                   end if
                end if
            end if
            if (ro(3,ip) > Zero) then
                if (MFHCOverlapSingleMove(ip, Zero, Zero, -dz, jp)) then ! overlap with jp
                   if (ro(3,ip)*ro(3,jp) < Zero) then                    ! overlap across z = 0
                      shs2 = shs2 + One
!               write(*,'(a,2(i5,3f10.3,3x))') 'hs overlap: ip, jp:', ip, ro(1:3,ip), jp, ro(1:3,jp)
!               do jp = 3, np
!                  if (jp == ip) cycle
!                  r1 = sqrt(sum((ro(1:3,ip)-ro(1:3,jp))**2))
!                  if (r1 < 10*dz+2*radat(iatan(jp))) &
!                  write(*,'(5x,a,i5,a,3f10.3)') 'jp = ', jp, ' ro = ', ro(1:3,jp)
!               end do
                   end if
                end if
            end if
         end do
         if (shs1 /= shs2) write(*,*) 'shsx', shs1, shs2   ! test of symmetry
         var(9)%value = shs1*dzi
      else
         var(9)%value = Zero
      end if

! ... sample the area * density of the remaing particles abs(z) = 0

      partsum = count(abs(ro(3,3:np)) < thickness2)
      var(10)%value = partsum*thicknessi

! ... sample the area * density of the remaining particles abs(z) = cyllen/2

      partsum = count(abs(ro(3,3:np)) > Half*cyllen - thickness2)
      var(11)%value = partsum*thicknessi

      var(12)%value = var(10)%value - var(11)%value

! ... total mean force

      var(13)%value = var(4)%value + var(8)%value + var(9)%value + var(12)%value

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,3f10.4)') 'tickness for fideal            = ', thickness
      write(uout,'(a,t35,3f10.4)') 'displacement parameter for fhs = ', dz
      write(uout,'(a,t35,3f10.4)') 'position of part 1             = ', ro(1:3,1)
      write(uout,'(a,t35,3f10.4)') 'position of part 2             = ', ro(1:3,2)

      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)

      data(1)     = Two*ro(3,2)                    ! separation
      data(2)     = thickness                      ! thickness for evaluation of fideal
      data(3:4)   = [ var(13)%avs1, var(13)%avsd ] ! f/kT
      data(5:6)   = [ var(4)%avs1, var(4)%avsd ]   ! fcorr/kT
      data(7:8)   = [ var(8)%avs1, var(8)%avsd ]   ! fbond/kT
      data(9:10)  = [ var(9)%avs1, var(9)%avsd ]   ! fhs
      data(11:12) = [ var(12)%avs1, var(12)%avsd ] ! fideal/kT (net)
      call WriteVecAppend(12, data(1:12), 'meanforce2.data')

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

!........................................................................

function MFHCOverlapSingleMove(ip, dx, dy, dz, jp)
   use MolModule
   implicit none
   integer(4), intent(in)  :: ip
   real(8)   , intent(in)  :: dx, dy, dz
   integer(4), intent(out) :: jp
   integer(4) :: iploc
   logical :: lhsoverlap
   logical :: MFHCOverlapSingleMove

   nptm  = 1
   iploc = 1
   ipnptm(iploc) = ip
   rotm(1,iploc) = ro(1,ip) + dx
   rotm(2,iploc) = ro(2,ip) + dy
   rotm(3,iploc) = ro(3,ip) + dz
   call SetTrialAtomProp
   if (lmonoatom) then
      call UTwoBodyANew(lhsoverlap, jp)
   else
      call UTwoBodyPNew(lhsoverlap, jp)
   endif
   MFHCOverlapSingleMove = lhsoverlap

end function MFHCOverlapSingleMove

!........................................................................

end subroutine MeanForce2

!************************************************************************
!> \page static static.F90
!! **UBondAcrossZ**
!! *calculate potential energy from bonds across z = 0*
!************************************************************************


subroutine UBondAcrossZ(ubondacross, fbondacross)

   use MolModule
   implicit none

   real(8), intent(out) :: ubondacross
   real(8), intent(out) :: fbondacross(1:3)

   integer(4) :: ic, ict, iseg, ip, jp
   real(8)    :: dx, dy, dz, r1, r2, term, fac

   ubondacross = Zero
   fbondacross(1:3) = Zero
   do ic = 1, nc
      ict = ictcn(ic)
      iseg = 1
      ip = ipnsegcn(iseg,ic)
      do iseg = 2, npct(ict)
         jp = ipnsegcn(iseg,ic)
         if (r(3,ip)*r(3,jp) < Zero) then   ! ip and jp are residing in different box halves
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            r1 = sqrt(r2)

            term = r1-bond(ict)%eq
            fac  = bond(ict)%k*term**(bond(ict)%p-1)
            ubondacross = ubondacross + fac*term
            term =-bond(ict)%p*fac/r1
            fbondacross(1) = fbondacross(1) + (term * dx)*sign(One,ro(3,ip))
            fbondacross(2) = fbondacross(2) + (term * dy)*sign(One,ro(3,ip))
            fbondacross(3) = fbondacross(3) + (term * dz)*sign(One,ro(3,ip))
         end if
         ip = jp
      end do
   end do

end subroutine UBondAcrossZ

!************************************************************************
!> \page static static.F90
!! **PotMeanForce**
!! *calculate potential of mean force*
!************************************************************************


!     system: cylindrical geometry
!     pmf calculated between two particles of type iptz. they are assumed to be on the z-axis
!     pmf is set to Zero at cyllen/2

!> \page nmlPotMeanForce
!! The namelist  \ref nmlPotMeanForce contains variables that control the calculation of the potential mean force between two particles
!! (only \ref txbc='cyl'). The two particles of type \ref iptpmf are assumed to be located on the z-axis. The pmf is set to zero at the
!! separation \ref cyllen /2.
!! * Variables:
!!  * \subpage iptpmf
!!  * \subpage rpmfzero

!> \page iptpmf
!! `integer`
!! **default:** `1`
!! * Types of particles for which the potential of mean force is calculated. \ref nppt(\ref iptpmf) = 2 is required.

!> \page rpmfzero
!! `real`
!! **default:** \ref cyllen `/2`
!! * Separation at which pmf is set to zero.

subroutine PotMeanForce(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='PotMeanForce'
   character(80), parameter :: txheading ='potential of mean force between two particles'
   integer(4),    save :: iptpmf
   real(8),       save :: rpmfzero
   integer(4),    save :: nbin, nvar
   type(df_var), allocatable,  save :: var(:)
   integer(4) ::  ivar, ibin
   real(8) :: value, weight, MCWeightInverse

   namelist /nmlPotMeanForce/ iptpmf, rpmfzero

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      if (.not.lbccyl) call Stop(txroutine, '.not.lbccyl', uout)

      iptpmf = 1
      rpmfzero = half*cyllen

      rewind(uin)
      read(uin,nmlPotMeanForce)

      if (nppt(iptpmf) /= 2) call Stop (txroutine,'nppt(iptpmf) /= 2',uout) ! 2 particles of type iptpmf ?

   case (iWriteInput)

      if (radat(iptpmf) > 0) then
         nbin = 20*(cyllen/radat(iptpmf))+0.5
      else
         nbin = 200
      end if
      if (nbin > mnbin_df) call Stop(txroutine, 'nbin > mnbin_df', uout)

! ... set nvar as well as allocate memory

      nvar = 2
      allocate(var(nvar))

! ... set label, min, max, and nbin

      var(1)%label = 'pmf'
      var(1)%min = Zero
      var(1)%max = cyllen
      var(1)%nbin = nbin
      var(2)%label = 'z-mean'
      var(2)%min =-Half*cyllen
      var(2)%max =+Half*cyllen
      var(2)%nbin = nbin
      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      weight = One
      if (lmcweight) weight = MCWeightInverse()                            ! unbrella sampling

      ivar = 1
      value = abs(ro(3,ipnpt(iptpmf))-ro(3,ipnpt(iptpmf)+1))              ! distance between the particles
      ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
      var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+weight

      ivar = 2
      value = sum(ro(3,ipnpt(iptpmf):ipnpt(iptpmf)+nppt(iptpmf)-1))/nppt(iptpmf)  ! com of the particles
      ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
      var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)+weight

   case (iAfterMacrostep)

      call DistFuncNorm(1, nvar, var)
      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var)

! ... transform from probabilty to potential of mean force

      ivar = 1
      nbin = var(ivar)%nbin
      where(var(ivar)%avs1(1:nbin-1) > Zero)
         var(ivar)%avsd(1:nbin-1) = abs(log(var(ivar)%avs1(1:nbin-1)+var(ivar)%avsd(1:nbin-1))-log(var(ivar)%avs1(1:nbin-1)))
         var(ivar)%avs1(1:nbin-1) =-log(var(ivar)%avs1(1:nbin-1))                  ! pmf = -ln p
      endwhere
      ibin = max(-1,min(floor(var(ivar)%bini*(rpmfzero-var(ivar)%min)),int(nbin))) ! ibin of pmfzero
      var(ivar)%avs1(1:nbin-1) = var(ivar)%avs1(1:nbin-1)-var(ivar)%avs1(ibin)     ! vertical shift

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,3i10)') 'type of particles for pmf calc = ', iptpmf
      write(uout,'(a,t35,f10.3)')'distance where pmf is zero     = ', rpmfzero
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine PotMeanForce

!************************************************************************
!> \page static static.F90
!! **SurfaceArea**
!! *calculate surface area exposed by all particles of one type*
!************************************************************************


!> \page nmlSurfaceArea
!! The namelist  \ref nmlSurfaceArea contains variables that control the calculation of the surface area available for a spherical prob.
!! * Variables:
!!  * \subpage nmlSurfaceArea_ipt
!!  * \subpage rprobe
!!  * \subpage wradat
!!  * \subpage nrandom

!> \page nmlSurfaceArea_ipt ipt
!! `integer`
!! * Type of particles of interest for area determination.

!> \page rprobe
!! `real`
!! **default:** `1.7`
!! * Radius of probe particle

!> \page wradat
!! `real`(1:nat)
!! **default:** nat*`0.0`
!! * van der Waal radius of atoms.

!> \page nrandom
!! `integer`
!! **default:** `10000`
!! * Number of random numbers.

subroutine SurfaceArea(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SurfaceArea'
   character(80), parameter :: txheading ='surface area available around atoms residing in particles of specified type'
   type(scalar_var), allocatable, save :: var(:)
   real(8),    allocatable,       save :: xran(:), yran(:), zran(:)
   integer(4), allocatable,       save :: nsurf(:)         ! list of atoms beloning to particles of type ipt
   integer(4),                    save :: ipt              ! particle type of interest
   real(8)   ,                    save :: rprobe           ! radius of probe particle
   real(8)   , allocatable,       save :: wradat(:)        ! van der waal radius of atom types
   integer(4),                    save :: nrandom          ! number of random numbers
   integer(4),                    save :: nvar
   integer(4) :: ip, ia, ja, ialoc, ialow, iat, jat, ivar, jvar, m
   real(8) :: xtry, ytry, ztry, radiat, dx, dy, dz, r2
   character(4) :: charip, charia

   namelist /nmlSurfaceArea/ ipt, rprobe, wradat, nrandom

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(wradat)) then
         allocate(wradat(nat))
         wradat = 0.0E+00
      end if

      rprobe = 1.7
      wradat = zero
      nrandom = 10000

      read(uin,nmlSurfaceArea)

   case (iWriteInput)

! ... allocate memory

      allocate(var(na_alloc), xran(nrandom), yran(nrandom), zran(nrandom), nsurf(na_alloc))
      xran = 0.0E+00
      yran = 0.0E+00
      zran = 0.0E+00
      nsurf = 0

! ... set the number and the list of surface atoms as well as label and norm

      ivar = 0
      do ip = ipnpt(ipt), ipnpt(ipt)+nppt(ipt)-1
         write(charip,'(i4)') ip
         ialow = ianpn(ip)-1
         do ialoc = 1, napt(ipt)
            ia = ialow+ialoc
            write(charia,'(i4)') ia
            ivar = ivar + 1
            nsurf(ivar) = ia
            var(ivar)%label = charip//' '//txpt(iptpn(ip))//'     '//charia//' '//txat(iatan(ia))
            var(ivar)%norm = 4.0d0*Pi*(wradat(iatan(ia))+rprobe)**2/real(nrandom)
         end do
      end do
      nvar = ivar

! ... generate the random coordinates

      do m = 1, nrandom
         call SphRandom(iseed, xran(m), yran(m), zran(m))
      end do

   case (iBeforeSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call ScalarSample(iStage, 1, nvar, var)

   case (iSimulationStep)

      var%value = zero
      do ivar = 1, nvar                                    ! loop over all atoms making up the surface
         ia = nsurf(ivar)
         iat = iatan(ia)
         radiat = wradat(iat)+rprobe
         do m = 1, nrandom                                 ! loop over all points of the surface of atom ia
            xtry = r(1,ia)+radiat*xran(m)
            ytry = r(2,ia)+radiat*yran(m)
            ztry = r(3,ia)+radiat*zran(m)
            do jvar = 1, nvar                              ! check whether the point is too close another atom
               if (jvar == ivar) cycle
               ja = nsurf(jvar)
               jat = iatan(ja)
               dx = r(1,ja)-xtry
               dy = r(2,ja)-ytry
               dz = r(3,ja)-ztry
               call PBCr2(dx,dy,dz,r2)
               if (r2 < (wradat(jat)+rprobe)**2) goto 100 ! check to see if this point does not contribute to the area
            end do
            var(ivar)%value = var(ivar)%value + One
  100       continue
         end do
      end do

      call ScalarSample(iStage, 1, nvar, var)

   case (iAfterMacrostep)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 1)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call ScalarSample(iStage, 1, nvar, var)
      call ScalarNorm(iStage, 1, nvar, var, 0)
      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,a)')    'particle type                  = ', txpt(ipt)
      write(uout,'(a,t35,f8.2)') 'radius of probe                = ', rprobe
      write(uout,'(a,t35,i8  )') 'number of random points        = ', nrandom
      write(uout,'()')
      write(uout,'(a,t15,a,t30,a)') 'atom no', 'atom type', 'van der waal radius'
      write(uout,'(a,t15,a,t30,a)') '-------', '---------', '-------------------'
      write(uout,'(i3,t15,a,t30,f8.3)') (iat, txat(iat), wradat(iat), iat = 1, nat)
      write(uout,'()')
      call ScalarWrite(iStage, 1, nvar, var, 1, '(a,t35,4f15.5,f15.0)', uout)
      write(uout,'()')
      write(uout,'(a,f8.2)') 'total surface area exposed     = ', sum(var%avs1)

      deallocate(var, xran, yran, zran, nsurf)

   end select

end subroutine SurfaceArea

!************************************************************************
!> \page static static.F90
!! **Crystalformat**
!! *write atoms in the crystalographic format starting with atoms*
!************************************************************************

!     belonging to molecules closest to (xx,yy,zz) first and the other in increasing distance

subroutine Crystalformat(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='crystalformat'
   character(90), parameter :: txheading ='atoms in particles written in crystalographic format at increasing particle separation'
   integer(4),  allocatable :: index(:)
   real(8)   ,  allocatable :: r2(:)
   integer(4) :: npmax, m, ip, ipt, ia, ialow, ialoc, iat
   real(8) :: rr(3), rmax, dx, dy, dz

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iAfterSimulation)

! ... set parameters

      rr = [Zero, Zero, Zero]
      rmax = 12.0
      npmax = 50

      allocate(r2(np_alloc), index(np_alloc))
      r2 = 0.0E+00
      index = 0

! ... set up the order

      do ip = 1, np
         dx = ro(1,ip)-rr(1)
         dy = ro(2,ip)-rr(2)
         dz = ro(3,ip)-rr(3)
         call PBCr2(dx,dy,dz,r2(ip))
      end do
      call HeapSortIndex(np, r2, index)

! ... write results

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,3f10.3)') 'origin                         = ', rr(1:3)
      write(uout,'(a,t35,f8.2)')   'cutoff distance                = ', rmax
      write(uout,'(a,t35,i8  )')   'maximum number of particles    = ', npmax
      write(uout,*)
      if (ilist > 0)  write(ulist,*)
      write(uout,'(a,a)') 'c ', txtitle(1:70)
      if (ilist > 0)  write(ulist,'(a,a)') 'c ', txtitle(1:70)
      do m = 1, min(npmax,np)
         ip = index(m)
         if (r2(ip) > rmax**2) cycle
         ipt = iptpn(ip)
         ialow = ianpn(ip)-1
         do ialoc = 1, napt(ipt)
            ia = ialow+ialoc
            iat = iatan(ia)
            write(uout,'(a,3f12.4)') txat(iat)(1:1), r(1:3,ia)
            if (ilist > 0) write(ulist,'(a,3f12.4)') txat(iat)(1:1), r(1:3,ia)
         end do
      end do
      write(uout,'(a)') 'end'
      if (ilist > 0) write(ulist,'(a)') 'end'

      deallocate(r2, index)

   end select

end subroutine Crystalformat

!************************************************************************
!> \page static static.F90
!! **Trajectory**
!! *write the trajectory on flist using every iskip timestep*
!************************************************************************


!> \page nmlTrajectory
!! The namelist  \ref nmlTrajectory contains variables that control the output of a trajectory.
!! * Variables:
!!  * \subpage nmlTrajectory_iskip
!!  * \subpage nmlTrajectory_rmin
!!  * \subpage nmlTrajectory_rmax

!> \page nmlTrajectory_iskip iskip
!! `integer`
!! **default:** `10`
!! * Write every trajectory data for every iskip timestep.

!> \page nmlTrajectory_rmin rmin
!! `real`(1:3)
!! **default:** `-boxlen2(1:3)`
!! * Lower, left, and front corner of box bounding the particles.

!> \page nmlTrajectory_rmax rmax
!! `real`(1:3)
!! **default:** `boxlen2(1:3)`
!! * Upper, right, and back corner of box bounding the particles.

   subroutine Trajectory(iStage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='trajectory'
   character(90), parameter :: txheading ='write trajectory on flist using every iskip timesteps, spatial restriction'
   integer(4),              save :: iskip, nsamp, nhit
   real(8),                 save :: rmin(3), rmax(3), dx, dy, dz
   integer(4), allocatable, save :: indx(:)
   real(8),    allocatable, save :: roold(:,:), ronopbc(:,:)
   integer(4) :: ip
   real(8)    :: GetTStep

   namelist /nmlTrajectory/ iskip, rmin, rmax

   select case (iStage)
   case (iReadInput)

      iskip = 10
      rmin =-boxlen2
      rmax =+boxlen2

      rewind(uin)
      read(uin,nmlTrajectory)

   case (iBeforeSimulation)

      allocate(indx(np_alloc), roold(3,np_alloc), ronopbc(3,np_alloc))
      indx = 0
      roold = 0.0E+00
      ronopbc = 0.0E+00

      nsamp = 0
      nhit = 0
      do ip = 1, np
         if ((ro(3,ip) > rmin(3)) .and. (ro(3,ip) < rmax(3))) then
            if ((ro(2,ip) > rmin(2)) .and. (ro(2,ip) < rmax(2))) then
               if ((ro(1,ip) > rmin(1)) .and. (ro(1,ip) < rmax(1))) then
                  nhit = nhit+1
                  indx(nhit) = ip
               end if
            end if
         end if
      end do

! ... save initial coordiantes

      do ip = 1, nhit
         ronopbc(1:3,indx(ip)) = ro(1:3,indx(ip))
         roold(1:3,indx(ip)) = ro(1:3,indx(ip))
      end do
      if (ilist > 0) write(ulist,'(a)') 'trajectory output'

   case (iBeforeMacrostep)

   case (iSimulationStep)

      nsamp = nsamp+1
      if (mod(nsamp,iskip) == 0) then
         do ip = 1, nhit
            dx = ro(1,indx(ip))-roold(1,indx(ip))            ! displacement between adjacent time internals (not PBC corrected)
            dy = ro(2,indx(ip))-roold(2,indx(ip))
            dz = ro(3,indx(ip))-roold(3,indx(ip))
            call PBC(dx,dy,dz)                               ! displacement two adjacent time internals (PBC corrected)
            roold(1:3,indx(ip)) = ro(1:3,indx(ip))           ! save actual position
            ronopbc(1,indx(ip)) = ronopbc(1,indx(ip))+dx     ! update the position (no pbc)
            ronopbc(2,indx(ip)) = ronopbc(2,indx(ip))+dy     ! update the position (no pbc)
            ronopbc(3,indx(ip)) = ronopbc(3,indx(ip))+dz     ! update the position (no pbc)
         end do
         if (ilist > 0) write(ulist,'(i5,2(3f8.3,2x))') (ip, ronopbc(1:3,indx(ip)),ro(1:3,indx(ip)),ip = 1,nhit)
      end if

   case (iAfterMacrostep)

   case (iAfterSimulation)

      call WriteHead(2, txheading, uout)
      write(uout,'(a,t35,3f12.3)') 'rmin (x, y, z)                 ', rmin(1:3)
      write(uout,'(a,t35,3f12.3)') 'rmax (x, y, z)                 ', rmax(1:3)
      write(uout,'(a,t35,i12)')    'number of particles selected   ', nhit
      write(uout,'(a,t35,i12)')    'number of time intervals       ', nsamp/iskip
      write(uout,'(a,t35,g12.3)')  'time step                      ', GetTStep()
      write(uout,'(a,t35,g12.3)')  'time interval of trajectory    ', GetTStep()*iskip
      write(uout,'(a,t35,g12.3)')  '(tstep*iskip)                 '
      write(uout,'()')
      write(uout,'(a)') 'particles selected:'
      write(uout,'(40i5)') indx(1:nhit)

      deallocate(indx, roold, ronopbc)

   end select

end subroutine Trajectory

!************************************************************************
!> \page static static.F90
!! **ExcessAmount**
!! *calculate excess amount*
!************************************************************************


subroutine ExcessAmount(nbin, not, txot, ipntot, var, boxlen, uout)

   use StatisticsModule
   implicit none
   character(80), parameter :: txheading ='excess adsorbed amount calculated from z-number density distribution'

   integer(4), intent(in)   :: nbin          ! number of bins
   integer(4), intent(in)   :: not           ! number of object types
   character(*), intent(in) :: txot(*)       ! name of object types
   integer(4), intent(in)   :: ipntot(*)     ! pointer from object type to variable number
   type(df_var),intent(in)  :: var(*)        ! variable contaning number density function in the z-direction
   real(8), intent(in)      :: boxlen(3)     ! box length
   integer(4), intent(in)   :: uout          ! output unit

   integer(4) :: ibinlow, ibinupp, iot, ivar
   real(8)    :: fac, nBulk, nTot, nExe

   if (minval(boxlen(1:3)) <= 0.0d0) return   ! not a box

   ibinlow = nbin/4                          ! lower division between "surface" and "bulk" layer
   ibinupp = nbin-1-ibinlow                  ! upper division between "surface" and "bulk" layer
   fac = real(nbin)/real(ibinupp-ibinlow+1)

   call WriteHead(2, txheading, uout)
   write(uout,'(a)') 'surfaces at z = -Lz/2 and Lz/2 is assumed'
   ivar = ipntot(1)
   write(uout,'(a,t35,f8.2,a,f8.2)') 'bulk density calculated between = ', &
     var(ivar)%min+var(ivar)%bin*ibinlow, ' and', var(ivar)%min+var(ivar)%bin*(ibinupp+1)
   write(uout,*)
   write(uout,'(a)') 'object type      nBulk       nTot        nExe       nExe/2      rhoBulk     nExe/(2*Lx*Ly)'
   write(uout,'(a)') '-----------      -----       ----        ----       ------      -------     --------------'
   do iot = 1, not
      ivar = ipntot(iot)
      nBulk = var(ivar)%bin*sum(var(ivar)%avs1(ibinlow:ibinupp))*(boxlen(1)*boxlen(2))
      nTot  = var(ivar)%bin*sum(var(ivar)%avs1(0:nbin-1))*(boxlen(1)*boxlen(2))
      nExe  = nTot-nBulk*fac
      write(uout,'(a,t15,6g12.4)') txot(iot), nBulk, nTot, nExe, nExe/2, &
      nBulk/(boxlen(1)*boxlen(2)*boxlen(3)/fac), nExe/(2*boxlen(1)*boxlen(2))
   end do
end subroutine ExcessAmount

!************************************************************************
!> \page static static.F90
!! **SubStructureDF**
!! *calculate distribution functions of properties of a SubStructure*
!************************************************************************


!> \page nmlSubStructureDF
!! The namelist  nmlSustructureDF contains variables that control the calculation of substructures.
!!
!!    | type | label    | quantity                                                                                |
!!    | ---- | -----    | ----------------------------------------------------------------------------------------|
!!    | 1    | rgsub    | radius of gyration of substructure                                                      |
!!    | 2    | rden     | radial number density                                                                   |
!!    | 3    | rrden    | reduced radial number density                                                           |
!!    | 4    | com-sub  | distance of com of chain to center of mass of substructure                              |
!!    | 5    | rg-sub   | average rg of a chain a a function of the distance to the center of mass of substructure|
!!    | 6    | z(r)     | sum of all charges as a function of the distance to the center of mass of substructure  |
!!    | 7    | zcum(r)  | cumulative charge as a function of the distance to the center of mass of substructure   |
!!    | 8    | alpha    | global degree of ionization distribution                                                |
!!    | 9    | a(r)     | reduced radial degree of ionization                                                     |
!!
!! * Variables:
!!  * \subpage nmlSubStructureDF_vtype
!!  * \subpage lptinsub

!> \page nmlSubStructureDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:3)
!! * Flag for engagement, lower end, upper end, number of bins, flag for normalization, title, and number of variable of vtype.
!! * Min: /0.0/ Max: /100.0/

!> \page lptinsub
!! `logical`(1:\ref npt)
!! * `.true.`:  mass of particles of type ipt are used to evaluate mass center of substructure.


subroutine SubStructureDF(iStage)

   use MolModule
   use MollibModule, only: InvInt

   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SubStructureDF'
   character(80), parameter :: txheading ='distribution functions of a substructure'
   integer(4)   , parameter :: ntype = 9
   type(static1D_var),         save :: vtype(ntype)
   integer(4),                 save :: nvar
   type(df_var),  allocatable, save :: var(:)
   type(scalar_var), allocatable, save :: sclvar(:)
   integer(4),    allocatable, save :: ipnt(:,:)
   integer(4), save                 :: lp
   logical, save  :: lptinsub(mnpt)                  ! particle type ipt part of substructure?

   type(chainprop_var) :: ChainProperty

   real(8)           :: msub
   real(8)   , save  :: msubi, mpt(mnpt)
   integer(4)        :: ntitr, ncharged   ! CHZA: number of titratable groups in substructure (type 8)
   real(8)   , save  :: ntitri            ! CHZA: inverse number of titratable groups in substructure (type 8)

   real(8)    :: rcom(3), rgsub, alpha
   integer(4) :: ipt, ip, ivar, ibin, itype, igr, ict, ic
   real(8)    :: InvFlt
   real(8)    :: r2, r1, dr(3), vsum, norm, dvol

   namelist /nmlSubStructureDF/ vtype, lptinsub

   if (slave) return                   ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

! ... default parameters

      vtype%l    = .false.
      vtype%min  = Zero
      vtype%max  = 100.0d0
      vtype%nbin = 100
      lptinsub = .false.

      rewind(uin)
      read(uin,nmlSubStructureDF)

! ... condition: vtype(7) needs vtype(6) to be evaluated   ! CHZA

      if (vtype(7)%l .and. (.not. vtype(6)%l)) then
         call Warn(txroutine,'vtype(6)%l .and. (.not. vtype(5)%l): vtype(5) = .true.',uout)
         vtype(6) = vtype(7)
      else if (vtype(6)%l .and. vtype(7)%l) then
         vtype(7) = vtype(6)
      end if   ! CHZA

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

      vtype%label = ['<Rg>            ','<rdens>         ','<rrden>         ','<com-sub>       ','<Rg_chain> - com',&
                    &'<Z(r)>          ','<Z_cum(r)>      ','<alpha>         ','<a(r)>          ']   ! CHZA

      vtype(1)%nvar = 1
      vtype(2)%nvar = ngr(1)
      vtype(3)%nvar = ngr(1)
      vtype(4)%nvar = nct
      vtype(5)%nvar = nct
      vtype(6)%nvar = 1    ! CHZA
      vtype(7)%nvar = 1    ! CHZA
      vtype(8)%nvar = 1    ! CHZA
      vtype(9)%nvar = 1    ! CHZA

! ... set nvar and allocate memory

      nvar = sum(vtype%nvar, 1, vtype%l)
      allocate(var(nvar), ipnt(maxval(vtype(1:ntype)%nvar),ntype)) !for distribution function
      ipnt = 0
      allocate(sclvar(nvar)) ! for scalar measurement

! ... set ipnt, label, min, max, and nbin

      ivar = 0

      itype = 1
      if (vtype(itype)%l) then
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min = vtype(itype)%min
         var(ivar)%max = vtype(itype)%max
         var(ivar)%nbin = vtype(itype)%nbin
         sclvar(ivar)%label = trim(vtype(itype)%label)
         sclvar(ivar)%norm = One
      end if

      itype = 2
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            ivar = ivar+1
            ipnt(igr,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
            sclvar(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
            sclvar(ivar)%norm = One
         end do
      end if

      itype = 3    ! CHZA
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            ivar = ivar+1
            ipnt(igr,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
            sclvar(ivar)%label = trim(vtype(itype)%label)//' '//txgr(igr)
            sclvar(ivar)%norm = One
         end do
      end if

      itype = 4
      if (vtype(itype)%l) then
         do ict = 1, nct
            ivar = ivar+1
            ipnt(ict,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
            sclvar(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
            sclvar(ivar)%norm = InvInt(ncct(ict))
         end do
      end if

      itype = 5
      if (vtype(itype)%l) then
         do ict = 1, nct
            ivar = ivar+1
            ipnt(ict,itype) = ivar
            var(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
            var(ivar)%min = vtype(itype)%min
            var(ivar)%max = vtype(itype)%max
            var(ivar)%nbin = vtype(itype)%nbin
            sclvar(ivar)%label = trim(vtype(itype)%label)//' '//txct(ict)
            sclvar(ivar)%norm = InvInt(ncct(ict))
         end do
      end if

      itype = 6   ! CHZA
      if (vtype(itype)%l) then
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min = vtype(itype)%min
         var(ivar)%max = vtype(itype)%max
         var(ivar)%nbin = vtype(itype)%nbin
         sclvar(ivar)%label = trim(vtype(itype)%label)
         sclvar(ivar)%norm = One
      end if

      itype = 7   ! CHZA
      if (vtype(itype)%l) then
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min = vtype(itype)%min
         var(ivar)%max = vtype(itype)%max
         var(ivar)%nbin = vtype(itype)%nbin
         sclvar(ivar)%label = trim(vtype(itype)%label)
         sclvar(ivar)%norm = One
      end if

      itype = 8
      if (vtype(itype)%l) then
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min = vtype(itype)%min
         var(ivar)%max = vtype(itype)%max
         var(ivar)%nbin = vtype(itype)%nbin
         sclvar(ivar)%label = trim(vtype(itype)%label)
         sclvar(ivar)%norm = One
         ntitr  = sum(nppt(1:npt),1,MASK=(lptinsub(1:npt) .and. latweakcharge(1:npt)))
         ntitri = One/real(ntitr)
      end if

      itype = 9   ! CHZA
      if (vtype(itype)%l) then
         ivar = ivar+1
         ipnt(1,itype) = ivar
         var(ivar)%label = trim(vtype(itype)%label)
         var(ivar)%min = vtype(itype)%min
         var(ivar)%max = vtype(itype)%max
         var(ivar)%nbin = vtype(itype)%nbin
         sclvar(ivar)%label = trim(vtype(itype)%label)
         sclvar(ivar)%norm = One
      end if

      call DistFuncSample(iStage, nvar, var) ! iStage: iWriteInput
      ! -> Initiate bin and bini

   case (iBeforeSimulation)

      msub = Zero   ! CHZA: Mass calculation of substructure
      mpt(1:npt) = masspt(1:npt)
      if (sum(mpt(1:npt),MASK=lptinsub(1:npt)) == Zero) mpt(1:npt) = One
      msub = sum(mpt(1:npt)*nppt(1:npt),MASK=lptinsub(1:npt))
      msubi = InvFlt(msub)   ! CHZA: Inverse mass of substructure

      do ipt = 1, npt   ! locate one particle in substructure to be used as origin of the substructure
         if(lptinsub(ipt)) then
            lp = ipnpt(ipt)
            exit
         end if
      end do

      call DistFuncSample(iStage, nvar, var)    ! iStage: iBeforeSimulation
      ! -> Initiate nsamp1, avs1, avsd
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

      call ScalarSample(iStage, 1, nvar, sclvar)! iStage: iBeforeSimulation
      ! -> Initiate nsamp1, avs1, avsd, fls1, flsd
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) sclvar

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)    ! iStage: iBeforeMacrostep
      ! -> Initiate nsamp2, avs2, nsampbin
      call ScalarSample(iStage, 1, nvar, sclvar)! iStage: iBeforeMacrostep
      ! -> Inititate nsamp2, avs2, fls2

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1
      sclvar%value = Zero

! ... calculate coordinates of center of mass

      rcom = Zero
      do ipt = 1, npt
         if (.not. lptinsub(ipt)) cycle
         do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
            dr(1:3) = ro(1:3,ip) - ro(1:3,lp)
            call PBC(dr(1),dr(2),dr(3))
            rcom(1:3) = rcom(1:3) + mpt(ipt)*dr(1:3)
         end do
      end do
      rcom(1:3) = rcom(1:3)*msubi + ro(1:3,lp)

      call PBC(rcom(1),rcom(2),rcom(3))   ! CHZA: apply PBC on rcom after shift by coordinates of lp

! ... sample type 1

      itype = 1
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         rgsub = Zero
         do ipt = 1, npt
            if (.not. lptinsub(ipt)) cycle
            do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
               dr(1:3) = ro(1:3,ip) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               rgsub = rgsub + mpt(ipt)*r2
            end do
         end do
         rgsub = sqrt(rgsub * msubi)
         ibin = max(-1,min(floor(var(ivar)%bini*(rgsub-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         sclvar(ivar)%value = sclvar(ivar)%value + rgsub
      end if

! ... sample type 2

      itype = 2
      if (vtype(itype)%l) then
         do ip = 1, np
            igr = igrpn(ip,1)
            if (igr <= 0) cycle
            dr(1:3) = ro(1:3,ip) - rcom(1:3)
            call PBCr2(dr(1), dr(2), dr(3), r2)
            r1 = sqrt(r2)
            ivar = ipnt(igr,itype)
            ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
            var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            sclvar(ivar)%value = sclvar(ivar)%value + r1
         end do

         do igr = 1, ngr(1)
            ivar = ipnt(igr,itype)
            sclvar(ivar)%value = sclvar(ivar)%value*InvFlt(grvar(igrpnt(1,igr))%value)   ! CHZA: to check
         end do
      end if

! ... sample type 3

      itype = 3
      if (vtype(itype)%l) then
         do ip = 1, np
            igr = igrpn(ip,1)
            if (igr <= 0) cycle
            ivar = ipnt(igr,itype)
            do ibin = -1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar - ngr(1))%avs2(ibin)
            end do
            sclvar(ivar)%value = sclvar(ivar - ngr(1))%value
         end do
      end if

! ... sample type 4 / 5

      if (vtype(4)%l .or. vtype(5)%l) then

         do ic = 1, nc
            ict = ictcn(ic)
            call UndoPBCChain(ro(1,ipnsegcn(1,ic)), ic, 1, vaux)
            call CalcChainProperty(ic, vaux, ChainProperty)
            dr(1:3) = ChainProperty%ro(1:3) - rcom(1:3)
            call PBCr2(dr(1), dr(2), dr(3), r2)
            r1 = sqrt(r2)

            if(vtype(4)%l) then
               itype = 4
               ivar = ipnt(ict,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               sclvar(ivar)%value = sclvar(ivar)%value + r1
            end if

            if(vtype(5)%l) then
               itype = 5
               ivar = ipnt(ict,itype)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + sqrt(ChainProperty%rg2)
               sclvar(ivar)%value = sclvar(ivar)%value + sqrt(ChainProperty%rg2)
               var(ivar)%nsampbin(ibin)=var(ivar)%nsampbin(ibin) + One
            end if
         end do
      end if

! ... sample type 6

      itype = 6
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         do ipt = 1, npt
            if (zat(iatpt(ipt)) == Zero) cycle
            do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
               if (lweakcharge .and. (.not. laz(ip))) cycle
               dr(1:3) = ro(1:3,ip) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + zat(iatpt(iptpn(ip)))
               var(ivar)%nsampbin(ibin) = var(ivar)%nsampbin(ibin) + One
               sclvar(ivar)%value = sclvar(ivar)%value + zat(iatpt(iptpn(ip)))
            end do
         end do
      end if

! ... sample type 7

      itype = 7
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         do ibin = -1, var(ivar)%nbin
            var(ivar)%avs2(ibin) = sum(var(ivar-1)%avs2(-1:ibin))
            var(ivar)%nsampbin(ibin) = sum(var(ivar-1)%nsampbin(-1:ibin))
            sclvar(ivar)%value = sclvar(ivar-1)%value
         end do
      end if

! ... sample type 8

      itype = 8
      if (vtype(itype)%l) then
         ivar     = ipnt(1,itype)
         alpha    = Zero
         ncharged = Zero
         do ipt = 1, npt
            if ((.not. lptinsub(ipt)) .or. (.not. latweakcharge(iatpt(ipt)))) cycle
            ncharged = ncharged + count(laz(ipnpt(ipt):ipnpt(ipt)+nppt(ipt)-1))
         end do
         alpha = real(ncharged) * ntitri
         ibin = max(-1,min(floor(var(ivar)%bini*(alpha-var(ivar)%min)),int(var(ivar)%nbin)))
         var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
         sclvar(ivar)%value = sclvar(ivar)%value + alpha
      end if

! ... sample type 9

      itype = 9
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         do ipt = 1, npt
            if ((.not. lptinsub(ipt)) .or. (.not. latweakcharge(iatpt(ipt)))) cycle
            do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
               dr(1:3) = ro(1:3,ip) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               if (laz(ip)) var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               var(ivar)%nsampbin(ibin) = var(ivar)%nsampbin(ibin) + One
            end do
         end do
         sclvar(ivar)%value = sclvar(ivar)%value + sum(var(ivar)%avs2(-1:var(ivar)%nbin)/var(ivar)%nsampbin(-1:var(ivar)%nbin))
      end if

      call ScalarSample(iStage, 1, nvar, sclvar)! iStage: iSimulationStep
      ! -> sample nsamp2, avs2, fls2

   case (iAfterMacrostep)

! ... normalisation

      call DistFuncNorm(1,sum(vtype(1:2)%nvar,1,vtype(1:2)%l), var)   ! type 1-2
      call DistFuncNorm(sum(vtype(1:3)%nvar,1,vtype(1:3)%l)+1, sum(vtype(1:4)%nvar,1,vtype(1:4)%l), var)   ! type 4
      call DistFuncNorm(sum(vtype(1:7)%nvar,1,vtype(1:7)%l)+1, sum(vtype(1:8)%nvar,1,vtype(1:8)%l), var)   ! type 8

      itype = 3
      if (vtype(itype)%l) then
         do igr = 1, ngr(1)
            ivar = ipnt(igr,itype)
            vsum = sum(var(ivar)%avs2(-1:var(ivar)%nbin))
            norm = var(ivar)%nsamp2 * (var(ivar)%max**3-var(ivar)%min**3) * InvFlt(vsum)
            do ibin = -1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) * norm / dvol(ibin,var(ivar)%min,var(ivar)%bin)
            end do
         end do
      end if

      itype = 5
      if (vtype(itype)%l) then
         do ict = 1, nct
            ivar = ipnt(ict,itype)
            norm = var(ivar)%nsamp2             ! factor to counteract the normalization in DistFuncSample
            do ibin = -1, var(ivar)%nbin
               var(ivar)%avs2(ibin) = norm*var(ivar)%avs2(ibin)*InvFlt(var(ivar)%nsampbin(ibin))
            end do
         end do
      end if

      itype = 9
      if (vtype(itype)%l) then
         ivar = ipnt(1,itype)
         vsum = sum(var(ivar)%avs2(-1:var(ivar)%nbin))
         norm = var(ivar)%nsamp2 * sum(var(ivar)%nsampbin(-1:var(ivar)%nbin)) * InvFlt(vsum) ! *nsamp2 in order to counteract wrong normalization in distfuncsample
         do ibin = -1, var(ivar)%nbin
            if (var(ivar)%nsampbin(ibin) /= Zero) var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) * norm / var(ivar)%nsampbin(ibin)
         end do
      end if

      call DistFuncSample(iStage, nvar, var)  ! iStage: iAfterMacrostep
      ! -> update nsamp1, divide avs2 by nsamp2, sum up avs1 and avsd
      if (lsim .and. master) write(ucnf) var

      call ScalarSample(iStage, 1, nvar, sclvar) ! iStage: iAfterMacrostep
      call ScalarNorm(iStage, 1, nvar, sclvar, 0) ! iStage: iAfterMacrostep
      if (lsim .and. master) write(ucnf) sclvar

      call WriteHead(2, txheading, uout)
      if (nct > 1) write(uout,'(a)') 'substructure '
      if (nct > 1) write(uout,'(a)') '-------------'
      call ScalarWrite(iStage, 1, nvar, sclvar, 1, '(a,t35,4f15.5,f15.0)', uout) ! iStage: iAfterMacrostep

   case (iAfterSimulation)

      call DistFuncSample(iStage, nvar, var) ! iStage: iAfterSimulation

      call ScalarSample(iStage, 1, nvar, sclvar)
      call ScalarNorm(iStage, 1, nvar, sclvar, 0)

      call WriteHead(2, txheading, uout)
      if (nct > 1) write(uout,'(a)') 'substructure '
      if (nct > 1) write(uout,'(a)') '-------------'
      call ScalarWrite(iStage, 1, nvar, sclvar, 1, '(a,t35,4f15.5,f15.0)', uout)

      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, ipnt, sclvar)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine SubStructureDF

!************************************************************************
!> \page static static.F90
!! **NetworkDF**
!! *calculate network distribution functions*
!************************************************************************


!     also final average over networks of same type and the spread of their df

!> \page nmlNetworkDF
!! The namelist  \ref nmlNetworkDF contains variables that control the calculation of network distribution functions. Distribution
!! functions are calculated for each network. Any combination of the types of distribution functions listed below may be selected
!! through vtype\%l.
!!
!!    | type | label |   quantity               |
!!    | ---- | ----- |   ---------------------- |
!!    | 1    | rg    |   radius of gyration     |
!!    | 2    | asph  |   asphericity            |
!!    | 3    | alpha |   degree of ionization   |
!!
!! * Variables:
!!  * \subpage nmlNetworkDF_vtype

!> \page nmlNetworkDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:3)
!! * Flag for engagement, lower end, upper end, number of bins. Other flags are not used.
!! * Min: /0.0,0.0,0.0/ Max: /100.0,1.0,1.0/


subroutine NetworkDF(iStage)

   use MolModule
   implicit none

   integer(4), intent(in)           :: iStage

   character(40), parameter         :: txroutine ='NetworkDF'
   character(80), parameter         :: txheading ='network distribution functions'
   integer(4)   , parameter         :: ntype = 3
   type(static1D_var),        save  :: vtype(ntype)
   integer(4),                save  :: nvar
   type(df_var), allocatable, save  :: var(:)
   integer(4),   allocatable, save  :: ipnt(:,:,:)
   integer(4),                save  :: nvar2
   type(df_var), allocatable, save  :: var2(:)
   integer(4),   allocatable        :: ilow(:), iupp(:)
   real(8),      allocatable        :: vspread(:)
   type(networkprop_var)            :: NetworkProperty
   integer(4), save                 :: ngrloc(ntype)

   integer(4)     :: itype, ivar, ivar2, ibin, inw, inwt, igrloc
   real(8)        :: value
   character(3)   :: txnwn

   namelist /nmlNetworkDF/ vtype

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (slave) return ! only master

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

      vtype%l    = .false.
      vtype%min  = [ Zero   , Zero   , Zero   ]
      vtype%max  = [ 100.d0 , One    , One    ]
      vtype%nbin = 100

      rewind(uin)
      read(uin,nmlNetworkDF)

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine,'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype

      vtype%label = ['rg   ','asph ','alpha' ]
      vtype%nvar  = [ nnw   , nnw   , nnw    ]
      ngrloc(1:ntype) = vtype(1:ntype)%nvar / nnw

! ... set nvartype and nvar as well as allocate memory

      nvar = sum(vtype(1:ntype)%nvar,dim=1,mask=vtype%l)
      allocate(var(nvar),ipnt(maxval(ngrloc(1:ntype)),nnw,ntype))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do inw = 1, nnw
               write (txnwn,"(i3)") inw
               do igrloc = 1, ngrloc(itype)
                  ivar = ivar+1
                  ipnt(igrloc,inw,itype) = ivar
                  var(ivar)%min  = vtype(itype)%min
                  var(ivar)%max  = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
                  var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txnwn))    ! for those with ngrloc = 1
               end do
            end do
         end if
      end do

      call DistFuncSample(iStage, nvar, var)

   case (iBeforeSimulation)

      call DistFuncSample(iStage, nvar, var)
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage, nvar, var)

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do inw = 1, nnw
         call CalcNetworkProperty(inw,NetworkProperty)

! ... sample type 1 to 3

         do itype = 1, 3
            if (vtype(itype)%l) then
               igrloc = 1
               ivar = ipnt(igrloc,inw,itype)
               if (itype == 1) then
                  value = sqrt(NetworkProperty%rg2)
               else if (itype == 2) then
                  value = NetworkProperty%asph
               else if (itype == 3) then
                  value = NetworkProperty%alpha
               end if
               ibin = max(-1,min(floor(var(ivar)%bini*(value-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end if
         end do

      end do

   case (iAfterMacrostep)

      call DistFuncNorm(1,nvar,var)
      call DistFuncSample(iStage,nvar,var)
      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      call DistFuncSample(iStage,nvar,var)
      call WriteHead(2,txheading,uout)
      call DistFuncHead(nvar,var,uout)
      call DistFuncWrite(txheading,nvar,var,uout,ulist,ishow,iplot,ilist)
      call DistFuncAverValue(nvar,var,uout)

! ..............    make average of distribution functions over networks of same type   .............

! ... set nvartype2 and nvar2 as well as allocate memory

      vtype%nvar = nnwt
      nvar2 = sum(vtype%nvar,1,vtype%l)
      allocate(ilow(nvar2),iupp(nvar2),var2(nvar2),vspread(nvar2))
      ilow = 0
      iupp = 0
      vspread = 0.0E+00

! ... set label, min, max, nbin, and bin

      ivar2 = 0
      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do inwt = 1, nnwt
               ivar2 = ivar2+1
               ivar = ivar+1
               ilow(ivar2) = ivar                    ! lower index of networks of type inwt
               ivar = ivar+nnwnwt(inwt)-1
               iupp(ivar2) = ivar                    ! upper index of networks of type inwt
               var2(ivar2)%label = trim(vtype(itype)%label)//' '//txnwt(inwt)
               var2(ivar2)%min = vtype(itype)%min
               var2(ivar2)%max = vtype(itype)%max
               var2(ivar2)%nbin = vtype(itype)%nbin
               var2(ivar2)%bin = (var2(ivar2)%max-var2(ivar2)%min)/var2(ivar2)%nbin
            end do
         end if
      end do

      call DistFuncAverDist(nvar2,ilow,iupp,var,var2,vspread)

      call DistFuncWrite(trim(txheading)//': aver',nvar2,var2,uout,ulist,ishow,iplot,ilist)
      write(uout,'()')
      write(uout,'(a,(t30,f10.4))') 'rms spread = ', vspread

      call DistFuncAverValue(nvar2,var2,uout)

      deallocate(var,ipnt,ilow,iupp,var2,vspread)

   end select

   if (ltime) call CpuAdd('stop',txroutine,1,uout)

end subroutine NetworkDF

!************************************************************************
!> \page static static.F90
!! **NetworkRadialDF**
!! *calculate radial distribution functions of properties of networks*
!************************************************************************


!> \page nmlNetworkRadialDF
!! The namelist  \ref nmlNetworkRadialDF contains variables that control the calculation of radial network distribution functions. Distribution
!! functions are calculated for each network. Any combination of the types of radial distribution functions listed below may be selected
!! through vtype\%l.
!!
!!    | type | label     |  quantity                              |
!!    | ---- | -----     |  ------------------------------------- |
!!    |  1   | rpart(r)  |  radial particle number distribution   |
!!    |  2   | rdens(r)  |  radial particle density distribution  |
!!    |  3   | rgchain(r)|  radial chain radius of gyration       |
!!    |  4   | q(r)      |  radial sum of all charges             |
!!    |  5   | qcum(r)   |  cumulative radial sum of all charges  |
!!    |  6   | alpha(r)  |  reduced degree of ionization          |
!!    |  7   | rchain(r) |  radial chain number distribution      |
!!
!! * Variables:
!!  * \subpage nmlNetworkRadialDF_vtype

!> \page nmlNetworkRadialDF_vtype vtype
!! `static1D_var(logical, real, real, integer, logical, character, real)`(1:7)
!! * Flag for engagement, lower end, upper end, number of bins. Other flags are not used.
!! * Min: /0.0,0.0,0.0/ Max: /100.0,100.0,100.0/

subroutine NetworkRadialDF(iStage)

   use MolModule

   implicit none

   integer(4),          intent(in) :: iStage

   character(40),        parameter :: txroutine ='NetworkRadialDF'
   character(80),        parameter :: txheading ='radial distribution functions of networks'
   integer(4)   ,        parameter :: ntype = 7
   type(static1D_var),        save :: vtype(ntype)
   integer(4),                save :: nvar
   integer(4),                save :: ngrloc(ntype)
   real(8),      allocatable, save :: nsampbin1(:,:)  ! nsampbin1 is raised by One if a property was assigned to a bin within one
                                                      ! macrostep. After the simulation is done a property was sampled nsampbin1
                                                      ! times. By dividing the sample by nsampbin1, the actual average of the
                                                      ! property in that bin is obtained.
   type(df_var), allocatable, save :: var(:)
   integer(4),   allocatable, save :: ipnt(:,:,:)
   type(networkprop_var)           :: NetworkProperty
   type(chainprop_var)             :: ChainProperty

   character(3)                    :: txinw
   integer(4)                      :: itype, ivar, ibin
   integer(4)                      :: inwt, inw, ict, ic, icloc, ipt, ip, iploc, igr, igrloc
   real(8)                         :: InvFlt
   real(8)                         :: rcom(1:3), r2, r1, vsum, norm, dvol, dr(3)

   namelist /nmlNetworkRadialDF/ vtype

   if (slave) return ! only master

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

! ... default parameters

      vtype%l    = .false.
      vtype%min  = Zero
      vtype%max  = 100.0d0
      vtype%nbin = 100

      rewind(uin)
      read(uin,nmlNetworkRadialDF)

! ... condition: vtype(2) needs vtype(1) to be evaluated

      if (vtype(2)%l .and. (.not. vtype(1)%l)) then
         call Warn(txroutine,'vtype(2)%l .and. (.not. vtype(1)%l): vtype(1) = .true.',uout)
         vtype(1) = vtype(2)
      else if (vtype(1)%l .and. vtype(2)%l) then
         vtype(2) = vtype(1)
      end if

! ... condition: vtype(5) needs vtype(4) to be evaluated

      if (vtype(5)%l .and. (.not. vtype(4)%l)) then
         call Warn(txroutine,'vtype(5)%l .and. (.not. vtype(4)%l): vtype(4) = .true.',uout)
         vtype(4) = vtype(5)
      else if (vtype(4)%l .and. vtype(5)%l) then
         vtype(5) = vtype(4)
      end if

      if (maxval(vtype%nbin) > mnbin_df) call Stop(txroutine, 'vtype%nbin > mnbin_df', uout)

   case (iWriteInput)

! ... set remaining elements of vtype and ngrloc

      vtype%label = ['<rpart>   ','<rdens>   ','<rgchain> ', &
                    &'<sum(q)>  ','<q_cum>   ','<alpha>   ', &
                    &'<rchain>  ' ]

      vtype%nvar  = [ ngr(1)*nnw,   ngr(1)*nnw,     nct*nnw, &
                    &        nnw,          nnw,         nnw, &
                    & ngr(1)*nnw  ]

      ngrloc(1:ntype) = vtype(1:ntype)%nvar / nnw   ! = ngr(1), ngr(1), nct, 1, 1, 1, ngr(1)

! ... set nvar and allocate memory

      nvar = sum(vtype(1:ntype)%nvar,1,vtype%l)
      allocate(var(nvar),ipnt(maxval(ngrloc(1:ntype)),nnw,ntype),nsampbin1(nvar,-1:maxval(vtype(1:ntype)%nbin)))
      ipnt = 0

! ... set ipnt, label, min, max, and nbin

      ivar = 0
      do itype = 1, ntype
         if (vtype(itype)%l) then
            do inw = 1, nnw
               write (txinw,"(i3)") inw
               do igrloc = 1, ngrloc(itype)
                  ivar = ivar+1
                  ipnt(igrloc,inw,itype) = ivar
                  var(ivar)%min  = vtype(itype)%min
                  var(ivar)%max  = vtype(itype)%max
                  var(ivar)%nbin = vtype(itype)%nbin
                  if      (itype == 1) then
                     var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txinw))//' '//txgr(igrloc)
                  else if (itype == 2) then
                     var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txinw))//' '//txgr(igrloc)
                  else if (itype == 3) then
                     var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txinw))//' '//txct(igrloc)
                  else if (itype == 7) then
                     var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txinw))//' '//txgr(igrloc)
                  else
                     var(ivar)%label = trim(vtype(itype)%label)//' inw:'//trim(adjustl(txinw))    ! for those with ngrloc = 1
                  end if
               end do
            end do
         end if
      end do

      call DistFuncSample(iStage,nvar,var) ! -> Initiate bin and bini

   case (iBeforeSimulation)

      call DistFuncSample(iStage,nvar,var) ! -> Initiate nsamp1, avs1, avsd
      nsampbin1 = Zero
      if (lsim .and. master .and. (txstart == 'continue')) read(ucnf) var

   case (iBeforeMacrostep)

      call DistFuncSample(iStage,nvar,var) ! -> Initiate nsamp2, avs2, nsampbin

   case (iSimulationStep)

      var%nsamp2 = var%nsamp2 + 1

      do inw = 1, nnw
         inwt = inwtnwn(inw)

! ... get center of mass of network inw and network type inwt
         call CalcNetworkProperty(inw,NetworkProperty)

         rcom(1:3) = NetworkProperty%ro(1:3)

! ... sample type 1

         itype = 1
         if (vtype(itype)%l) then
            do ip = 1, np
               igr = igrpn(ip,1)
               if (igr <= 0) cycle
               ivar = ipnt(igr,inw,itype)
               dr(1:3) = ro(1:3,ip) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end do
         end if

! ... sample type 2 in iStage == 'iAfterMacrostep'

! ... sample type 3

         itype = 3
         if (vtype(itype)%l) then
            do icloc = 1, ncnwt(inwt)
               ic  = icnclocnwn(icloc,inw)
               ict = ictcn(ic)
               ivar = ipnt(ict,inw,itype)
               call UndoPBCChain(ro(1:3,ipnsegcn(1,ic)),ic,1,vaux)
               call CalcChainProperty(ic,vaux,ChainProperty)
               dr(1:3) = ChainProperty%ro(1:3) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + sqrt(ChainProperty%rg2)
               var(ivar)%nsampbin(ibin) = var(ivar)%nsampbin(ibin) + One
            end do
         end if

! ... sample type 4

         itype = 4
         if (vtype(itype)%l) then
            ivar = ipnt(1,inw,itype)
            do ipt = 1, npt
               if (zat(iatpt(ipt)) == Zero) cycle
               do ip = ipnpt(ipt), ipnpt(ipt) + nppt(ipt) - 1
                  if (lweakcharge) then
                     if (.not.laz(ip)) cycle
                  end if
                  dr(1:3) = ro(1:3,ip) - rcom(1:3)
                  call PBCr2(dr(1), dr(2), dr(3), r2)
                  r1 = sqrt(r2)
                  ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + zat(iatpt(iptpn(ip)))
               end do
            end do
         end if

! ... sample type 5 in iStage == 'iAfterMacrostep'

! ... sample type 6

         itype = 6
         if (vtype(itype)%l) then
            ivar = ipnt(1,inw,itype)
            do iploc = 1, npnwt(inwt)
               ip = ipnplocnwn(iploc,inw)
               ipt = iptpn(ip)
               if (.not. latweakcharge(iatpt(ipt))) cycle
               dr(1:3) = ro(1:3,ip) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               if (laz(ip)) var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
               var(ivar)%nsampbin(ibin) = var(ivar)%nsampbin(ibin) + One
            end do
         end if

! ... sample type 7

         itype = 7
         if (vtype(itype)%l) then
            do ic = 1, nc
               igr = igrpn(ipnsegcn(1,ic),1)
               if (igr <= 0) cycle
               ivar = ipnt(igr,inw,itype)
               call UndoPBCChain(ro(1,ipnsegcn(1,ic)),ic,1,vaux)
               call CalcChainProperty(ic,vaux,ChainProperty)
               dr(1:3) = ChainProperty%ro(1:3) - rcom(1:3)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               r1 = sqrt(r2)
               ibin = max(-1,min(floor(var(ivar)%bini*(r1-var(ivar)%min)),int(var(ivar)%nbin)))
               var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) + One
            end do
         end if

      end do

   case (iAfterMacrostep)

! ... sample dependent distribution functions

         do inw = 1, nnw

! ... sample type 2 by copying data from type 1
            itype = 2
            if (vtype(itype)%l) then
               do igr = 1, ngr(1)
                  ivar = ipnt(igr,inw,itype)
                  var(ivar)%avs2 = var(ipnt(igr,inw,1))%avs2  ! reference to itype 1 in ipnt
               end do
            end if

! ... sample type 5 by copying and summing up data from type 4

            itype = 5
            if (vtype(itype)%l) then
               ivar = ipnt(1,inw,itype)
               var(ivar)%avs2(-1) = var(ipnt(1,inw,4))%avs2(-1)
               do ibin = 0, var(ivar)%nbin
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin-1)+var(ipnt(1,inw,4))%avs2(ibin) ! reference to itype 4 in ipnt
               end do
            end if

         end do

! ... normalisation

      do inw = 1, nnw

         itype = 2
         if (vtype(itype)%l) then
            do igr = 1, ngr(1)
               ivar = ipnt(igr,inw,itype)
               vsum = sum(var(ivar)%avs2(-1:var(ivar)%nbin))
               norm = var(ivar)%nsamp2 * (var(ivar)%max**3-var(ivar)%min**3) * InvFlt(vsum)
               do ibin = -1, var(ivar)%nbin
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin) * norm / dvol(ibin,var(ivar)%min,var(ivar)%bin)
               end do
            end do
         end if

         itype = 3
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ipnt(ict,inw,itype)
               norm = var(ivar)%nsamp2             ! factor to counteract the normalization in DistFuncSample
               do ibin = -1, var(ivar)%nbin
                  if (var(ivar)%nsampbin(ibin) > Zero) then
                     var(ivar)%avs2(ibin) = norm*var(ivar)%avs2(ibin)/var(ivar)%nsampbin(ibin)
                     nsampbin1(ivar,ibin) = nsampbin1(ivar,ibin)+One
                  end if
               end do
            end do
         end if

         itype = 6
         if (vtype(itype)%l) then
            ivar = ipnt(1,inw,itype)
            vsum = sum(var(ivar)%avs2(-1:var(ivar)%nbin))
            norm = var(ivar)%nsamp2 * sum(var(ivar)%nsampbin(-1:var(ivar)%nbin)) * InvFlt(vsum) ! *nsamp2 in order to counteract wrong normalization in distfuncsample
            do ibin = -1, var(ivar)%nbin
               if (var(ivar)%nsampbin(ibin) > Zero) then
                  var(ivar)%avs2(ibin) = var(ivar)%avs2(ibin)*norm/var(ivar)%nsampbin(ibin)
                  nsampbin1(ivar,ibin) = nsampbin1(ivar,ibin)+One
               end if
            end do
         end if

      end do

      call DistFuncSample(iStage, nvar, var) ! -> update nsamp1, divide avs2 by nsamp2, sum up avs1 and avsd

      if (lsim .and. master) write(ucnf) var

   case (iAfterSimulation)

      do inw = 1, nnw

         itype = 3
         if (vtype(itype)%l) then
            do ict = 1, nct
               ivar = ipnt(ict,inw,itype)
               norm = var(ivar)%nsamp1             ! factor to counteract the normalization in DistFuncSample
               do ibin = -1, var(ivar)%nbin
                  if (nsampbin1(ivar,ibin) > Zero) var(ivar)%avs1(ibin) = norm*var(ivar)%avs1(ibin)/nsampbin1(ivar,ibin)
               end do
            end do
         end if

         itype = 6
         if (vtype(itype)%l) then
            ivar = ipnt(1,inw,itype)
            norm = var(ivar)%nsamp1 ! *nsamp1 in order to counteract wrong normalization in distfuncsample
            do ibin = -1, var(ivar)%nbin
               if (nsampbin1(ivar,ibin) > Zero) var(ivar)%avs1(ibin) = norm*var(ivar)%avs1(ibin)/nsampbin1(ivar,ibin)
            end do
         end if

      end do

      call DistFuncSample(iStage, nvar, var)
      call DistFuncHead(nvar, var, uout)
      call DistFuncWrite(txheading, nvar, var, uout, ulist, ishow, iplot, ilist)

      deallocate(var, ipnt)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine NetworkRadialDF
