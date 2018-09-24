! ... 'version 6.4.5, Sep 18, 2015'

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
!> \page dynamic dynamic.F90
!! **DynamicModule**
!! *module for dynamic*
!************************************************************************


!> \page nmlDynamic
!! * The namelist  \ref nmlDynamic contains variables that control the dynamic analyses. The detailed description of the use of
!!   stochastic data and evaluation of the correlation functions are given by the the data structures sf_var and cf_input_var, see
!!   chapter 8.
!! * Variables:
!!  * \subpage sfnplow
!!  * \subpage sfnpupp
!!  * \subpage lmsd
!!  * \subpage lorix
!!  * \subpage loriy
!!  * \subpage loriz
!!  * \subpage lliv
!!  * \subpage lanv
!!  * \subpage lfor
!!  * \subpage ltor
!!  * \subpage lidm
!!  * \subpage lutot
!!  * \subpage itestdyn

module DynamicModule

   use MolModule

! ... data structure for stochastic data
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type sf_var
      integer(4)     :: nplow                            ! lower particle number
      integer(4)     :: npupp                            ! upper particle number
      integer(4)     :: ndim                             ! number of dimensions of the stochastic variable
      integer(4)     :: fac                              ! factor describing the separation of stochasitic data
      integer(4)     :: ngr                              ! number of groups
   end type sf_var

!     multiple-tau correation function; ref.  Ramirez, J. Chem. Phys. 133, 154103 (2010)

!                     idata
!     level     1        2       ....      nmean   ......     nlevel-1   nlevel
!
!       1       x        x                   x                   x          x
!               ------------------------------
!                 /
!                /
!       2       x        x                   x                   x          x
!               ------------------------------
!                 /
!                /
!                .
!                 /
!                /
!     nolevel   x        x                   x                   x          x

! ... data structure for correlation data, input
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type cf_input_var
      integer(4)     :: nmean                            ! length (number of values) forming a mean
      integer(4)     :: nlevel                           ! length (number of values) forming a level
      integer(4)     :: nolevel                          ! number of levels
      integer(4)     :: legendre                         ! order of Legendre polynomial (1, 2, or 3)
      logical        :: lsvalue                          ! logical flag if for single sf value
      logical        :: lsubmean                         ! logical flag if subtraction of mean of cf
      logical        :: lnorm                            ! logical flag if normalisation of cf
   end type cf_input_var

! ... data structure for correlation data
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type cf_var
      integer(4)     :: nmean                            ! length (number of values) forming a mean
      integer(4)     :: nlevel                           ! length (number of values) forming a level
      integer(4)     :: ratio                            ! nlevel/nmean
      integer(4)     :: nolevel                          ! number of levels
      integer(4)     :: legendre                         ! order of Legendre polynomial (1, 2, or 3)
      logical        :: lsvalue                          ! logical flag if for single sf value
      logical        :: lsubmean                         ! logical flag if subtraction of mean of cf
      logical        :: lnorm                            ! logical flag if normalisation of cf

      real(8), allocatable    :: sf(:,:,:,:)             ! stochastic function
      real(8), allocatable    :: sf_aver(:,:,:)          ! stochastic function, local average (of particle)
      real(8), allocatable    :: sf_mean(:,:)            ! stochastic function, gobal mean
      real(8), allocatable    :: cf(:,:,:)               ! time correlation function
      real(8), allocatable    :: cf2(:,:,:)              ! time correlation function squared
      integer(4), allocatable :: Np(:,:,:)               ! counter of correlation sampling (of particle)
      integer(4), allocatable :: Ngr(:,:,:)              ! counter of correlation sampling (of group)
      integer(8), allocatable :: Nlev(:,:)               ! counter of level (of particle)
   end type cf_var
!> \page sfnplow
!! `integer`
!! **default:** `1`
!! * Lower id of particle to be used in the dynamic analysis.
   integer(4) :: sfnplow
!> \page sfnpupp
!! `integer`
!! **default:** `np`
!! * Upper id of particle to be used in the dynamic analysis.
   integer(4) :: sfnpupp
!> \page itestdyn
!! `integer`
!! **default:** `0`
!! * Flag for test output.
   integer(4), save :: itestdyn                          ! for test output

end module DynamicModule

!************************************************************************
!> \page dynamic dynamic.F90
!! **DynamicDriver**
!! *driver of dynamic analysis routines*
!************************************************************************


!> \page lmsd
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Mean square displacements are calculated. Further specification is given in namelist  \ref nmlMSD .
!! * `.false.`:  Nothing.

!> \page lorix
!! `logical`
!! **default:** `.false.`
!! *

!> \page lorix
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Orientational tcf for the molecular x'-axis is calculated. Further specification is given in namelist  nmlOriXTCF .
!! * `.false.`:  Nothing.

!> \page loriy
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Orientational tcf for the molecular y'-axis is calculated. Further specification is given in namelist  nmlOriYTCF .
!! * `.false.`:  Nothing.

!> \page loriz
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Orientational tcf for the molecular z'-axis is calculated. Further specification is given in namelist  nmlOriZTCF .
!! * `.false.`:  Nothing.

!> \page lliv
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Velocity tcfs for each molecular axis (molecular frame) and total velocity tcf are calculated. Further specification is given in namelist  nmlLinVelCF.
!! * `.false.`: Nothing.

!> \page lanv
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Angular velocity tcfs for each molecular axis are calculated (molecular frame). Further specification is given in namelist  nmlAnvVelCF.
!! * `.false.`: Nothing.

!> \page lfor
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Force tcfs for each molecular axis and total force tcf are calculated (box frame). Further specification is given in namelist nmlForCF.
!! * `.false.`: Nothing.

!> \page ltor
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Torque tcfs for each molecular axis and total torque tcf are calculated (box frame). Further specification is given in namelist nmlTorCF.
!! * `.false.`: Nothing.

!> \page lidm
!! `logical`
!! **default:** `.false.`
!! * `.true.`:  Induced dipole moment tcfs for each molecular axis and total induced dipole moment tcf are calculated (molecular frame). Further specification is given in namelist nmlIDPMCF.
!! * `.false.`: Nothing.

!> \page lutot
!! `logical`
!! **default:** `.false.`
!! * `.true.`: Total energy tcf are calculated. Further specification is given in namelist nmlUtotCF.
!! * `.false.`: Nothing.

subroutine DynamicDriver(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='DynamicDriver'
   character(80), parameter :: txheading ='dynamic analysis: general'
   logical,       save :: lmsd, lorix, loriy, loriz, lliv, lanv, lfor, ltor, lidm, lutot

   namelist /nmlDynamic/ sfnplow, sfnpupp, &
                         lmsd, lorix, loriy, loriz, lliv, lanv, lfor, ltor, lidm, lutot, itestdyn

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

      sfnplow  = 1
      sfnpupp  = np
      lmsd     = .false.
      lorix    = .false.
      loriy    = .false.
      loriz    = .false.
      lliv     = .false.
      lanv     = .false.
      lfor     = .false.
      ltor     = .false.
      lidm     = .false.
      lutot    = .false.
      itestdyn = 0

      rewind(uin)
      read(uin,nmlDynamic)

      call DynamicDriverSub

   case (iWriteInput)

      call DynamicDriverSub

   case (iBeforeSimulation)

      call DynamicDriverSub

   case (iBeforeMacrostep)

      call DynamicDriverSub

   case (iSimulationStep)

      call DynamicDriverSub

   case (iAfterMacrostep)

      call DynamicDriverSub

   case (iAfterSimulation)

      if (master) then
         call WriteHead(2, txheading, uout)
         write(uout,'(a)') 'dynamic analysis routines used'
         write(uout,'(a)') '------------------------------'
         if (lmsd)         write(uout,'(a)') '   msd        '
         if (lorix)        write(uout,'(a)') '   orix       '
         if (loriy)        write(uout,'(a)') '   oriy       '
         if (loriz)        write(uout,'(a)') '   oriz       '
         if (lliv)         write(uout,'(a)') '   lin vel    '
         if (lanv)         write(uout,'(a)') '   ang vel    '
         if (lfor)         write(uout,'(a)') '   force      '
         if (ltor)         write(uout,'(a)') '   torque     '
         if (lidm)         write(uout,'(a)') '   idm        '
         if (lutot)        write(uout,'(a)') '   utot       '
      end if
      write(uout,*)
      write(uout,'(a,t35,i5)') 'itestdyn                       =', itestdyn

      if(master) call fileflush(uout)

      call DynamicDriverSub

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine DynamicDriverSub
   if (lmsd)  call MSD(iStage)
   if (lorix) call OrixTCF(iStage)
   if (loriy) call OriyTCF(iStage)
   if (loriz) call OrizTCF(iStage)
   if (lliv)  call LinVelTCF(iStage)
   if (lanv)  call AngVelTCF(iStage)
   if (lfor)  call ForTCF(iStage)
   if (ltor)  call TorTCF(iStage)
   if (lidm)  call IDMTCF(iStage)
   if (lutot) call UtotTCF(iStage)
end subroutine DynamicDriverSub

!........................................................................

end subroutine DynamicDriver

!************************************************************************
!> \page dynamic dynamic.F90
!! **MSD**
!! *mean square displacement*
!************************************************************************


!> \page nmlMSD
!! The namelist  \ref nmlMSD contains variables that control the calculation of the mean square displacement.
!! * Variables:
!!  * \subpage nmlMSD_sf
!!  * \subpage nmlMSD_cfin

!> \page nmlMSD_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** \ref sfnplow, \ref sfnpupp, 3, 1, -

!> \page nmlMSD_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, .false., .false.`

subroutine MSD(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='MSD'
   character(80), parameter :: txheading ='mean square displacement'
   character(3), parameter :: txfunc ='msd'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlMSD/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      cfin%lsubmean = .false.
      cfin%lnorm = .false.
      rewind(uin)
      read(uin,nmlMSD)
      rewind(uin)
      if ((cfin%nmean /= 1) .and. .not.cfin%lsvalue) then
          call Warn(txroutine, '(cfin%nmean /= 1) .and. .not.cfin%lsvalue', uout)
          return
      end if
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      if ((cfin%nmean /= 1) .and. .not.cfin%lsvalue) return
      call CFCalc(iStage, txfunc, ro(1,1), sf, cf)
   case (iAfterSimulation)
      if ((cfin%nmean /= 1) .and. .not.cfin%lsvalue) return
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine MSD

!************************************************************************
!> \page dynamic dynamic.F90
!! **OrixTCF**
!! *orientation time correlation function of paricle x'-axis*
!************************************************************************

!> \page nmlOriXTCF nmlOriXTCF, nmlOriYTCF and nmlOriZTCF
!! * The namelists  nmlOriXTCF ,  nmlOriYTCF , and  nmlOriZTCF contain variables which control the calculation of the three orientational tcf:s of the molecular axes.
!! * Variables:
!!  * \subpage nmlOriXTCF_sf
!!  * \subpage nmlOriXTCF_cfin

!> \page nmlOriXTCF_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** \ref sfnplow, \ref sfnpupp, 3, 3, -

!> \page nmlOriXTCF_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, -, -`
subroutine OrixTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='OrixTCF'
   character(80), parameter :: txheading ='orientation, x, time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlOrixTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 3
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlOrixTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, ori(1,1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OrixTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **OriyTCF**
!! *orientation time correlation function of paricle y'-axis*
!************************************************************************


subroutine OriyTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='OriyTCF'
   character(80), parameter :: txheading ='orientation, y, time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlOriyTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 3
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlOriyTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, ori(1,2,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OriyTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **OrizTCF**
!! *orientation time correlation function of paricle z'-axis*
!************************************************************************


subroutine OrizTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='OrizTCF'
   character(80), parameter :: txheading ='orientation, z, time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlOrizTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 3
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlOrizTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, ori(1,3,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine OrizTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **LinVelTCF**
!! *linear velocity time correlation function*
!************************************************************************


!> \page nmlLinVelTCF nmlLinVelTCF and nmlAngVelTCF
!! The namelists  nmlLinVelTCF and  nmlAngVelTCF contain variables that control the calculation of the velocity and angular velocity, respectively, tcf:s for the molecular axes.
!! * Variables:
!!  * \subpage nmlLinVelTCF_sf
!!  * \subpage nmlLinVelTCF_cfin

!> \page nmlLinVelTCF_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** \ref sfnplow, \ref sfnpupp, 3, 1, -

!> \page nmlLinVelTCF_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, -, -`

subroutine LinVelTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='LinVelTCF'
   character(80), parameter :: txheading ='linear velocity time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlLinVelTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlLinVelTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
 !!! later: change to particle coordinate system
      call CFCalc(iStage, txfunc, rod(1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine LinVelTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **AngVelTCF**
!! *angular velocity time correlation function*
!************************************************************************


subroutine AngVelTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='AngVelTCF'
   character(80), parameter :: txheading ='angular velocity time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlAngVelTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlAngVelTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, angvelo(1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine AngVelTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **ForTCF**
!! *force time correlation function*
!************************************************************************



!> \page nmlForTCF nmlForTCF and nmlTorTCF
!! The namelists  nmlForTCF and  nmlTorTCF contain variables that control the calculation of the force and torque, respectively, tcf:s about the box axes.
!! * Variables:
!!  * \subpage nmlForTCF_sf
!!  * \subpage nmlForTCF_cfin

!> \page nmlForTCF_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** \ref sfnplow, \ref sfnpupp, 3, 1, -

!> \page nmlForTCF_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, -, -`
subroutine ForTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='ForTCF'
   character(80), parameter :: txheading ='force time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlForTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlForTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, forceo(1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine ForTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **TorTCF**
!! *torque time correlation function*
!************************************************************************


subroutine TorTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='TorTCF'
   character(80), parameter :: txheading ='torque time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlTorTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlTorTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, torqueo(1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine TorTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **IDMTCF**
!! *induced dipole moment time correlation function*
!************************************************************************

!> \page nmlIDMTCF
!! The namelist  \ref nmlIDMTCF contains variables that control the calculation of the induced dipole moment tcf for the molecular axes.
!! * Variables:
!!  * \subpage nmlIDMTCF_sf
!!  * \subpage nmlIDMTCF_cfin

!> \page nmlIDMTCF_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** \ref sfnplow, \ref sfnpupp, 3, 1, -

!> \page nmlIDMTCF_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, -, -`

subroutine IDMTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='IDMTCF'
   character(80), parameter :: txheading ='induced dipole moment time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlIDMTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = sfnplow
      sf%npupp = sfnpupp
      sf%ndim = 3
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlIDMTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate',sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, idm(1,1), sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine IDMTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **UtotTCF**
!! *total potential time correlation function*
!************************************************************************


!> \page nmlUtotTCF
!! The namelists  \ref nmlUtotTCF contain variables that control the calculation of tcf of the total potential energy.
!! * Variables:
!!  * \subpage nmlUtotTCF_sf
!!  * \subpage nmlUtotTCF_cfin

!> \page nmlUtotTCF_sf sf
!! `sf_var(int, int, int, int, int)`
!! **default:** `1, 1, 1, 1, -`

!> \page nmlUtotTCF_cfin cfin
!! `cf_input_var(int, int, int, int, logical, logical, logical)`
!! **default:** `-, -, -, 1, -, -, -`

subroutine UtotTCF(iStage)

   use DynamicModule
   implicit none

   integer(4), intent(in) :: iStage
   character(40), parameter :: txroutine ='UtotTCF'
   character(80), parameter :: txheading ='total potential energy time correlation function'
   character(3), parameter :: txfunc ='tcf'
   type(sf_var), save :: sf
   type(cf_input_var) :: cfin
   type(cf_var), save :: cf

   namelist /nmlUtotTCF/ sf, cfin

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   select case (iStage)
   case (iWriteInput)
      sf%nplow = 1
      sf%npupp = 1
      sf%ndim = 1
      sf%fac = 1
      cfin%legendre = 1
      rewind(uin)
      read(uin,nmlUtotTCF)
      call PrepareDynamic(sf, cfin, cf)
      call MemoryDynamic('allocate', sf, cf)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
   case (iSimulationStep)
      call CFCalc(iStage, txfunc, [u%tot*sclene/(np*GasConstant*temp*scltem)], sf, cf)
   case (iAfterSimulation)
      call CFCalc(iStage, txfunc, [zero], sf, cf)
      call CFWrite(txheading, sf, cf)
      call MemoryDynamic('deallocate', sf, cf)
   end select
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine UtotTCF

!************************************************************************
!> \page dynamic dynamic.F90
!! **PrepareDynamic**
!! *prepare for dynamic analyse*
!************************************************************************


subroutine PrepareDynamic(sf, cfin, cf)

   use DynamicModule
   implicit none

   character(40), parameter :: txroutine ='PrepareDynamic'

   type(sf_var), intent(inout) :: sf
   type(cf_input_var), intent(in) :: cfin
   type(cf_var), intent(out) :: cf

! ... copy from input variable

   sf%Ngr      = ngr(1)
   cf%nmean    = cfin%nmean
   cf%nlevel   = cfin%nlevel
   cf%nolevel  = cfin%nolevel
   cf%legendre = cfin%legendre
   cf%lsvalue  = cfin%lsvalue
   cf%lsubmean = cfin%lsubmean
   cf%lnorm    = cfin%lnorm

! ... check variables

   if (sf%nplow < 1) call Stop(txroutine, 'sf%nplow < 1', uout)
   if (sf%nplow > sf%npupp) call Stop(txroutine, 'sf%nplow > sf%npupp', uout)
   if (sf%ndim < 1) call Stop(txroutine, 'sf%ndim < 1', uout)
   if (sf%ndim > 3) call Stop(txroutine, 'sf%ndim > 3', uout)
   if (sf%ngr < 1) call Stop(txroutine, 'sf%ngr < 1', uout)
   if (sf%fac < 1) call Stop(txroutine, 'sf%fac < 1', uout)

   if (cf%nmean < 1) call Stop(txroutine, 'cf%nmean < 1', uout)
   if (cf%nlevel < 1) call Stop(txroutine, 'cf%nlevel < 1', uout)
   if (cf%nolevel < 1) call Stop(txroutine, 'cf%nolevel < 1', uout)
   if (cf%nmean > cf%nlevel) call Stop(txroutine, 'cf%nmean > cf%nlevel', uout)
   if (cf%legendre < 1) call Stop(txroutine, 'cf%legendre < 1', uout)
   if (cf%legendre > 3) call Stop(txroutine, 'cf%legendre > 3', uout)

end subroutine PrepareDynamic

!************************************************************************
!> \page dynamic dynamic.F90
!! **MemoryDynamic**
!! *allocate memory*
!************************************************************************


subroutine MemoryDynamic(str, sf, cf)

   use DynamicModule
   implicit none

   type(sf_var), intent(inout) :: sf
   type(cf_var), intent(inout) :: cf

   character(40), parameter :: txroutine ='MemoryDynamic'
   character(*) :: str
   integer(4) :: ierr

   if (str(1:1) == 'a') then
      if(.not. allocated(cf%sf)) then
         allocate(cf%sf(cf%nolevel,cf%nlevel,sf%ndim,np), stat = ierr)
      end if
      if(.not. allocated(cf%sf_aver)) then
         allocate(cf%sf_aver(sf%ndim,cf%nolevel,np), stat = ierr)
      end if
      if(.not. allocated(cf%sf_mean)) then
         allocate(cf%sf_mean(sf%ndim,sf%ngr), stat = ierr)
      end if
      if(.not. allocated(cf%cf)) then
         allocate(cf%cf(cf%nlevel,cf%nolevel,sf%ngr), stat = ierr)
      end if
      if(.not. allocated(cf%cf2)) then
         allocate(cf%cf2(cf%nlevel,cf%nolevel,sf%ngr), stat = ierr)
      end if
      if(.not. allocated(cf%Np)) then
         allocate(cf%Np(cf%nlevel,cf%nolevel,np), stat = ierr)
      end if
      if(.not. allocated(cf%Ngr)) then
         allocate(cf%Ngr(cf%nlevel,cf%nolevel,sf%ngr), stat = ierr)
      end if
      if(.not. allocated(cf%Nlev)) then
         allocate(cf%Nlev(cf%nolevel,np), stat = ierr)
      end if
      if (ierr /= 0) call WriteIOStat(txroutine, 'memory allocation failed', ierr, 2, 6)
   else if (str(1:1) == 'd') then
      if(allocated(cf%sf)) deallocate(cf%sf)
      if(allocated(cf%sf_aver)) deallocate(cf%sf_aver)
      if(allocated(cf%sf_mean)) deallocate(cf%sf_mean)
      if(allocated(cf%cf)) deallocate(cf%cf)
      if(allocated(cf%cf2)) deallocate(cf%cf2)
      if(allocated(cf%Np)) deallocate(cf%Np)
      if(allocated(cf%Ngr)) deallocate(cf%Ngr)
      if(allocated(cf%Nlev)) deallocate(cf%Nlev)
   end if

end subroutine MemoryDynamic

!************************************************************************
!> \page dynamic dynamic.F90
!! **CFCalc**
!! *calculate correlation function*
!************************************************************************


subroutine CFCalc(iStage, txfunc, sfdata, sf, cf)

   use DynamicModule
   use MollibModule, only: InvInt
   implicit none

   integer(4),   intent(in) :: iStage
   character(3), intent(in) :: txfunc    ! 'msd', or 'tcf'
   real(8),      intent(in) :: sfdata(*) ! stochastic data
   type(sf_var), intent(in) :: sf        ! stochastic data
   type(cf_var), intent(inout) :: cf     ! correlation data

   character(40), parameter :: txroutine ='CFCalc'
   integer(4) :: ndim, ilev, idata, igr, ip, ishift

   select case (iStage)
   case (iWriteInput)

      cf%ratio = cf%nlevel/cf%nmean
      cf%sf = Zero
      cf%sf_aver = Zero
      cf%sf_mean = Zero
      cf%cf = Zero
      cf%cf2 = Zero
      cf%Np = 0
      cf%Ngr = 0
      cf%Nlev = 0
      if (mod(cf%nlevel,cf%nmean) /= 0) call stop(txroutine, 'mod(cf%nlevel,cf%nmean)=/0', 6)
      if (cf%legendre < 1 .or. cf%legendre > 3) call stop(txroutine, 'cf%legendre < 1 .or. cf%legendre > 3', 6)
      if (cf%nmean == 1 .and. cf%nolevel > 1) call warn(txroutine, 'cf%nmean == 1: cf%nolevel = 1 is sufficient', 6)

      if (itestdyn == 1) call FileOpen(uuser, fuser, 'form/noread')

   case (iSimulationStep)                    ! sample time correlation function one time

      ndim = sf%ndim
      do ip = sf%nplow, sf%npupp
         igr = igrpn(ip,1)                        ! determine group
         ishift = (ip-1)*sf%fac*ndim              ! separation between stochasitc data for successive ip
         if (itestdyn == 1) call TestCFCalc(1, uuser)
         if(cf%lnorm) cf%sf_mean(1:ndim,igr) = cf%sf_mean(1:ndim,igr) + sfdata(ishift+1:ishift+ndim)
         call CFCalcSub(txfunc, 1, ip, igr, sfdata(ishift+1:ishift+ndim), sf, cf)
      end do

   case (iAfterSimulation)                   ! normalize correlation function by the number of samplings

      if (itestdyn == 1) call TestCFCalc(2, uuser)
      do igr = 1, sf%ngr
         do ilev = 1, cf%nolevel
            do idata = 1, cf%nlevel
               if (itestdyn == 1) call TestCFCalc(3, uuser)
               cf%cf(idata,ilev,igr) = cf%cf(idata,ilev,igr)*InvInt(cf%Ngr(idata,ilev,igr))
               cf%cf2(idata,ilev,igr) = InvInt(cf%Ngr(idata,ilev,igr)) * &
                     ((cf%cf2(idata,ilev,igr)*InvInt(cf%Ngr(idata,ilev,igr))) - cf%cf(idata,ilev,igr)**2)
               cf%cf2(idata,ilev,igr) = (max(zero,cf%cf2(idata,ilev,igr)))
            end do
         end do
      end do

   end select

contains

!........................................................................

subroutine TestCFCalc(iopt, unit)
   integer(4),   intent(in) :: iopt
   integer(4),   intent(in) :: unit
   if (iopt == 1) then
      call WriteHead(3, 'Test'//trim(txroutine), uuser)
      write(unit,'(a,t20,2i5,3(3f10.5,2x))') 'ipt, ip, sfdata', iptpn(ip), ip, sfdata(ishift+1:ishift+3*ndim)
   else if (iopt == 2) then
      call WriteHead(3, 'Test'//trim(txroutine), uuser)
      write(uuser,'(a)')'igr, ilevel, idata, cf%cf(), cf%Ngr(idata,ilevel,igr)'
   else if (iopt == 3) then
      write(uuser,'(3i5,f12.6,i6)') igr, ilev, idata, cf%cf(idata,ilev,igr), cf%Ngr(idata,ilev,igr)
   end if
end subroutine TestCFCalc

!........................................................................

end subroutine CFCalc

!************************************************************************
!> \page dynamic dynamic.F90
!! **CFCalcSub**
!! *calculate time correlation function, subroutine*
!************************************************************************


recursive subroutine CFCalcSub(txfunc, ilev, ip, igr, data, sf, cf)

   use DynamicModule
   implicit none

   character(3), intent(in)    :: txfunc   ! 'msd', or 'tcf'
   integer(4), intent(in)      :: ilev     ! level of time correlation function
   integer(4), intent(in)      :: ip       ! particle
   integer(4), intent(in)      :: igr      ! group to which particle belongs to
   real(8), intent(in)         :: data(*)  ! stoichastic data
   type(sf_var), intent(in)    :: sf       ! stoichastic function
   type(cf_var), intent(inout) :: cf       ! time correlation function

   real(8), save :: dro(3), value0(3)
   real(8)     :: dx, dy, dz, corr, value(3)
   integer(4)  :: ndim, idata

   ndim = sf%ndim

! ... add stochastic data to start of S(ilev,:,1:ndim)

   cf%sf(ilev,:,:,ip) = EOSHIFT(cf%sf(ilev,:,:,ip), SHIFT = -1, DIM=1, BOUNDARY = data(1:ndim))

   if (itestdyn == 1) call TestCFCalcSub(1, uuser)

   if (txfunc == 'tcf') then

! ... sample time correlation function invoking data

      do idata = 1, min(cf%Np(1,ilev,ip) + 1, cf%nlevel)             ! loop over data in level idata
         if ((ilev == 1) .or. (idata > cf%ratio)) then               ! check if sample correlation
            corr = sum(data(1:ndim)*cf%sf(ilev,idata,1:ndim,ip))      ! evaluate correlation
            if(cf%legendre == 1) then
            else if(cf%legendre == 2) then
               corr = half*(three*corr**2 - one)
            else if(cf%legendre == 3) then
               corr = half*corr*(three*corr**2 - three)
            end if
            if (itestdyn == 1) call TestCFCalcSub(2, uuser)
            cf%cf(idata,ilev,igr) = cf%cf(idata,ilev,igr) + corr       ! sum correlation
            cf%cf2(idata,ilev,igr) = cf%cf2(idata,ilev,igr) + corr**2  ! sum correlation squared
         end if
         cf%Np(idata,ilev,ip) =  cf%Np(idata,ilev,ip) + 1            ! update number of terms per particle
         cf%Ngr(idata,ilev,igr) = cf%Ngr(idata,ilev,igr) + 1         ! update number of terms per group
      end do

   else if (txfunc == 'msd') then

! ... sample mean square displacement invoking data

!      write(*,*) "ilev, idata, cf%sf(ilev,idata,1,ip), value, dro(1)"
      dro(1:3) = zero                                                 ! set displacement holder to zero
      do idata = 1, min(cf%Np(1,ilev,ip) + 1, cf%nlevel)              ! loop over data in level idata
         if (idata == 1) then
            if(ilev == 1) value0 = data(1:ndim)
            value = value0
         else
            value = cf%sf(ilev,idata-1,:,ip)
         end if
         dx = value(1) - cf%sf(ilev,idata,1,ip)
         dy = value(2) - cf%sf(ilev,idata,2,ip)
         dz = value(3) - cf%sf(ilev,idata,3,ip)
         call PBC(dx,dy,dz)
         dro(1) = dro(1) + dx
         dro(2) = dro(2) + dy
         dro(3) = dro(3) + dz
         corr = dro(1)**2+dro(2)**2+dro(3)**2
!     write(*,'(2i10,3f12.4)') ilev, idata, cf%sf(ilev,idata,1,ip), value(1), dro(1)
         if (itestdyn == 1) call TestCFCalcSub(2, uuser)
         if ((ilev == 1) .or. (idata > cf%ratio)) then               ! check if sample correlation
            cf%cf(idata,ilev,igr) = cf%cf(idata,ilev,igr) + corr       ! sum correlation
            cf%cf2(idata,ilev,igr) = cf%cf2(idata,ilev,igr) + corr**2  ! sum correlation squared
         end if
         cf%Np(idata,ilev,ip) =  cf%Np(idata,ilev,ip) + 1            ! update number of terms per particle
         cf%Ngr(idata,ilev,igr) = cf%Ngr(idata,ilev,igr) + 1         ! update number of terms per group
      end do

   end if

! ... sample of mean

   if (.not.cf%lsvalue) then
      cf%sf_aver(1:ndim,ilev,ip) = cf%sf_aver(1:ndim,ilev,ip) + data(1:ndim) ! sum data (nmean times); for the next level
   else if(cf%Nlev(ilev,ip) == 0) then
      cf%sf_aver(1:ndim,ilev,ip) = cf%sf_aver(1:ndim,ilev,ip) + data(1:ndim)
   end if
   cf%Nlev(ilev,ip) = cf%Nlev(ilev,ip) + 1                             ! update number of terms
   if (itestdyn == 1) call TestCFCalcSub(3, uuser)

! ... pass data to next level

   if (cf%Nlev(ilev,ip) == cf%nmean) then                              ! nmean data sampled
      if (ilev < cf%nolevel) then                                      ! check that not last level
         if (.not.cf%lsvalue) cf%sf_aver(:,ilev,ip) = cf%sf_aver(:,ilev,ip)/cf%nmean  ! make local average
         if (itestdyn == 1) call TestCFCalcSub(4, uuser)
         call CFCalcSub(txfunc, ilev+1, ip, igr, cf%sf_aver(:,ilev,ip), sf, cf)   ! pass mean to next level
      end if
      cf%Nlev(ilev,ip) = 0                                             ! initiate for the next averaging
      cf%sf_aver(:,ilev,ip) = Zero                                     ! initiate for the next averaging
   end if

contains

!........................................................................

subroutine TestCFCalcSub(iopt, unit)
   integer(4),   intent(in) :: iopt
   integer(4),   intent(in) :: unit
   if (iopt == 1) then
      write(unit,*)
      write(unit,'(a,t45,i5)')     'CFCalcSub: ip', ip
      write(unit,'(a,t45,i5)')     'CFCalcSub: ilevel', ilev
      write(unit,'(a,t45,5f12.5)') 'CFCalcSub: cf%sf(ilevel,1:nlevel,1,ip)', cf%sf(ilev,:,1,ip)
   else if (iopt == 2) then
      write(unit,'(a,t45,i5,f12.6)') 'CFCalcSub: idata, corr', idata, corr
   else if (iopt == 3) then
      write(unit,'(a,t45,5f12.5)') 'CFCalcSub: local sum of sf', cf%sf_aver(:,ilev,ip)
   else if (iopt == 4) then
      write(unit,'(a,t85,5f12.5)') 'CFCalcSub: pass local average sf to next level',cf%sf_aver(:,ilev,ip)
   end if
end subroutine TestCFCalcSub

!........................................................................

end subroutine CFCalcSub

!************************************************************************
!> \page dynamic dynamic.F90
!! **CFWrite**
!! *write time correlation function*
!************************************************************************


subroutine CFWrite(txheading, sf, cf)

   use DynamicModule
   use MollibModule, only: InvInt
   implicit none

   character(*), intent(in) :: txheading ! heading
   type(sf_var), intent(in) :: sf        ! stoichastic function
   type(cf_var), intent(in) :: cf        ! time correlation function

   real(8), allocatable :: CF0(:), CF0sd(:), CFinf(:),  CFinfsd(:), CFnorm(:), CFnormsd(:)
   real(8) :: CFsub, CFsubsd, value, valuesd
   real(8) :: InvFlt
   real(8) :: t_low, t_mid, t_upp
   integer(4), allocatable :: nbin(:)
   integer(4) :: ilev, idata, igr, nzero, unit, i
   character(5) :: str

   allocate(CF0(1:sf%ngr), CF0sd(1:sf%ngr))
   CF0 = 0.0E+00
   CF0sd = 0.0E+00
   allocate(CFinf(1:sf%ngr), CFinfsd(1:sf%ngr))
   CFinf = 0.0E+00
   CFinfsd = 0.0E+00
   allocate(CFnorm(1:sf%ngr), CFnormsd(1:sf%ngr))
   CFnorm = 0.0E+00
   CFnormsd = 0.0E+00
   allocate(nbin(1:sf%ngr))
   nbin = 0

! ... initiate t_low, t_mid, and t_upp

   t_low = zero
   t_mid = half
   t_upp = one
   if (cf%lsvalue) then  ! the single stochastic data for the interval are taken at the upper end
      t_low = one
      t_mid = one
   end if


! ... calculate normalization variables

   do igr = 1, sf%ngr
      CF0(igr) = cf%cf(1,1,igr)
      CF0sd(igr) = cf%cf2(1,1,igr)
      CFinf(igr) = (InvInt(cf%Ngr(1,1,igr)))**2 * sum(cf%sf_mean(1:sf%ndim,igr)**2)
      CFinfsd(igr) =(max(zero,InvInt(cf%Ngr(1,1,igr)) * (CF0(igr) - CFinf(igr))))*4.0d0*abs(CFinf(igr))
      CFnorm(igr) = InvFlt(CF0(igr) - CFinf(igr))
      CFnormsd(igr) = CF0sd(igr) + CFinfsd(igr)
   end do

! ... calculate number of bins

      nbin = 0
      do igr = 1, sf%ngr
         do ilev = 1, cf%nolevel
            do idata = 1, cf%nlevel
               if ((ilev == 1) .or. (idata > cf%ratio)) then
                  if (cf%Ngr(idata,ilev,igr) > 0) nbin(igr) = nbin(igr) + 1
               end if
            end do
         end do
      end do

! ... write time correlation function of FOUT

   if (ishow /= 0) then
      unit = uout
!      if (itestdyn == 1) unit = uuser
      call WriteHead(2, txheading, unit)
      write(unit,'(a,t35,i5)') 'lower particle number          =', sf%nplow
      write(unit,'(a,t35,i5)') 'upper particle number          =', sf%npupp
      write(unit,'(a,t35,i5)') 'maximal number of groups       =', sf%ngr
      write(unit,'(a,t35,i5)') 'dimension of stochastic func   =', sf%ndim
      write(unit,'()')
      write(unit,'(a,t35,i5)') 'length of aver (nmean, m)      =', cf%nmean
      write(unit,'(a,t35,i5)') 'length of level (nlevel, p)    =', cf%nlevel
      write(unit,'(a,t35,i5)') 'number of levels (nolevel, l)  =', cf%nolevel
      write(unit,'(a,t35,i5)') 'order of Legendre polynomial   =', cf%legendre
      write(unit,'()')
      write(unit,'(a,t35,l )') 'single cf value                =', cf%lsvalue
      write(unit,'(a,t35,l )') 'subraction of average          =', cf%lsubmean
      write(unit,'(a,t35,l )') 'normalization                  =', cf%lnorm
      write(unit,'()')
      write(unit,'(a,t35,i5)') 'number of time steps           =', cf%nlevel*cf%nmean**(cf%nolevel-1)
      write(unit,'(a,t35,9i5)')'number of bins                 =', nbin(1:sf%ngr)
      do igr = 1, sf%ngr
         if (sum(cf%Ngr(:,:,igr)) == 0) cycle    ! emty correlation function
         write(unit,'()')
         write(str,'(i5)') igr
         write(unit,'(a)') 'group:'//trim(adjustl(str))
         write(unit,'(a)') '-------'
         write(unit,'()')
         write(unit,'(2(a,f10.4))') 'CF(0)   = ',CF0(igr),   '   +/-', sqrt(CF0sd(igr))
         write(unit,'(2(a,f10.4))') 'CF(inf) = ',CFinf(igr), '   +/-', CFinfsd(igr)
         write(unit,'()')
         write(unit,900) 'time (middle, lower, upper)', 'n samp', ('value', 'precision', i = 1,3)
         write(unit,900) '---------------------------', '------', ('-----', '---------', i = 1,3)
 900     format(a,t32,a,t48,a,t54,a,t73,a,t79,a,t98,a,t104,a)
         do ilev = 1, cf%nolevel
            do idata = 1, cf%nlevel
               if ((ilev == 1) .or. (idata > cf%ratio)) then
                  if (cf%Ngr(idata,ilev,igr) > 0) then
                     write(unit,'(f10.2,2f8.0,t30,i8)',advance='no') &
                     tau_val(t_mid), tau_val(t_low), tau_val(t_upp), cf%Ngr(idata,ilev,igr)
                     call cf_evaluate(.false., .false., value, valuesd)         ! lsubmean = .false., lnorm = .false.
                     write(unit,'(f15.4,f10.4)', advance='no') value, valuesd
                     call cf_evaluate(.true., .false., value, valuesd)          ! lsubmean = .true., lnorm = .false.
                     write(unit,'(f15.4,f10.4)', advance='no') value, valuesd
                     call cf_evaluate(.true., .true., value, valuesd)           ! lsubmean = .true., lnorm = .true.
                     write(unit,'(f15.4,f10.4)') value, valuesd
                  end if
               end if
            end do
         end do
      end do
   end if

! ... write time correlation function of FLIST

   if (ilist /= 0) then
      nzero = 0
      do igr = 1, sf%ngr
         if (sum(cf%Ngr(:,:,igr)) == 0) nzero = nzero + 1
      end do
      write(ulist,'(a)') txheading
      write(ulist,'(i5)') sf%ngr - nzero
      do igr = 1, sf%ngr
         if (sum(cf%Ngr(:,:,igr)) == 0) cycle   ! emty correlation function
         write(str,'(i5)') igr
         write(ulist,'(a)') 'gr:'//trim(adjustl(str))
         write(ulist,'(i6)') nbin(igr)
         do ilev = 1, cf%nolevel
            do idata = 1, cf%nlevel
               if ((ilev == 1) .or. (idata > cf%ratio)) then
                  if (cf%Ngr(idata,ilev,igr) > 0) then
                     call cf_evaluate(cf%lsubmean, cf%lnorm, value, valuesd)
                     write(ulist,'(f10.2,a,2(e20.4,a))') tau_val(t_mid), tab, value, tab, valuesd
                  end if
               end if
            end do
         end do
      end do
   end if

  deallocate(CF0, CF0sd, CFinf, CFinfsd, CFnorm, CFnormsd, nbin)

contains

!........................................................................

subroutine cf_evaluate(lsubmean,lnorm, value, valuesd)
   logical, intent(in) :: lsubmean
   logical, intent(in) :: lnorm
   real(8), intent(out) :: value
   real(8), intent(out) :: valuesd
         value = cf%cf(idata,ilev,igr)
         valuesd = sqrt(cf%cf2(idata,ilev,igr))
         CFsub = cf%cf(idata,ilev,igr)
         CFsubsd = cf%cf2(idata,ilev,igr)
         if(lsubmean) then
            value = value - CFinf(igr)
            valuesd = sqrt(valuesd**2+CFinfsd(igr))
            CFsub = CFsub - CFinf(igr)
            CFsubsd = CFsubsd + CFinfsd(igr)
         end if
         if(lnorm) then
            value = value*CFnorm(igr)
            valuesd = abs(CFsub*CFnorm(igr))*sqrt(CFsubsd*InvFlt(CFsub**2) + CFnormsd(igr)*CFnorm(igr)**2)
         end if
end subroutine cf_evaluate

real(8) function tau_val(fac)
      real(8), intent(in) :: fac
      real(8) :: taulow, deltau
      taulow = 1.0d0*(idata-1)*cf%nmean**(ilev-1)  ! tau value, lower end of interval
      deltau = cf%nmean**(ilev-1)-1                ! tau value, interval width
      tau_val = taulow + fac*deltau                ! tau value in the interval determined by fac
end function tau_val

!........................................................................

end subroutine CFWrite


