! ... 'version 6.4.7, Sep 18, 2015'

!************************************************************************
!> \page statistics statistics.F90
!! **StatisticsModule**
!! *module for statistics*
!************************************************************************


module StatisticsModule

! ... data structure for scalar quantities

   integer(4), parameter :: mnblocklen = 20
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type scalar_var
      character(40) :: label                             ! label
      real(8)       :: norm                              ! normalization factor
      integer(4)    :: nsamp1                            ! number of macrosteps sampled
      integer(4)    :: nsamp2                            ! number of values sampled per macrostep
      real(8)       :: avs1                              ! average of the run
      real(8)       :: avsd                              ! precision of average of the run
      real(8)       :: avs2                              ! average of a macrostep
      real(8)       :: fls1                              ! fluctuation of the run
      real(8)       :: flsd                              ! precision of fluctuation of the run
      real(8)       :: fls2                              ! fluctuation of a macrostep
      real(8)       :: value                             ! value of a configuration
      integer(4)    :: nsamp                             ! number of samplings
      integer(4)    :: nblocklen                         ! for sampling with variable blocklen
      integer(4)    :: nblock(mnblocklen)                ! for sampling with variable blocklen
      real(8)       :: av_s1(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: av_sd(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: av_s2(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: fl_s1(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: fl_sd(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: fl_s2(mnblocklen)                 ! for sampling with variable blocklen
      real(8)       :: av_sd_extrap                      ! for sampling with variable blocklen
      real(8)       :: av_sd_stateff                     ! for sampling with variable blocklen
      real(8)       :: fl_sd_extrap                      ! for sampling with variable blocklen
      real(8)       :: fl_sd_stateff                     ! for sampling with variable blocklen
   end type scalar_var

! ... data structure for one-dimensional distribution functions

   integer(4), parameter :: mnbin_df =1000
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type df_var
      character(27) :: label                             ! label
      real(8)       :: norm                              ! normalization factor
      integer(4)    :: nsamp1                            ! number of macrosteps sampled
      integer(4)    :: nsamp2                            ! number of values sampled per macrostep
      real(8)       :: min                               ! minimum value of df
      real(8)       :: max                               ! maximum value of df
      integer(4)    :: nbin                              ! number of grid points
      real(8)       :: bin                               ! grid length of df
      real(8)       :: bini                              ! inverse of bin
      real(8)       :: nsampbin(-1:mnbin_df)             ! number of values sampled in each bin during macrostep
      real(8)       :: avs1(-1:mnbin_df)                 ! average of the run
      real(8)       :: avsd(-1:mnbin_df)                 ! precision of average of the run
      real(8)       :: avs2(-1:mnbin_df)                 ! average of a macrostep
   end type df_var

! ... data structure for two-dimensional distribution functions

   integer(4), parameter :: mnbin_df2d = 40
!   integer(4), parameter :: mnbin_df2d = 100            ! Elliposids (Erik W)
! These are documented in the manual in Chapter 7 (file datastructures.md)
   type df2d_var
      character(27) :: label                             ! label of df
      real(8)       :: norm                              ! normalization factor
      integer(4)    :: nsamp1                            ! number of macrosteps sampled
      integer(4)    :: nsamp2                            ! number of values sampled per macrostep
      real(8)       :: min(2)                            ! minimum value of df
      real(8)       :: max(2)                            ! maximum value of df
      integer(4)    :: nbin(2)                           ! number of grid points
      real(8)       :: bin(2)                            ! grid length of df
      real(8)       :: bini(2)                           ! 1/bini
      real(8)       :: avs1(-1:mnbin_df2d,-1:mnbin_df2d) ! average of the run
      real(8)       :: avsd(-1:mnbin_df2d,-1:mnbin_df2d) ! precision of average of the run
      real(8)       :: avs2(-1:mnbin_df2d,-1:mnbin_df2d) ! average of a macrostep
   end type df2d_var

   logical    :: lblockaverwrite =.false.                ! control writing on blockaver.data

end module StatisticsModule

!************************************************************************
!> \page statistics statistics.F90
!! **SetlBlockAver**
!! *set variable \ref lblockaver*
!************************************************************************


subroutine SetlBlockAver(lwrite)

   use StatisticsModule
   implicit none
   logical, intent(in) :: lwrite

   lblockaverwrite = lwrite

end subroutine SetlBlockAver

!************************************************************************
!> \page statistics statistics.F90
!! **ScalarSample**
!! *sample scalar quantities*
!************************************************************************


subroutine ScalarSample(iStage, ilow, iupp, var)

   use StatisticsModule
   implicit none

   integer(4),       intent(in)  :: iStage
   integer(4),       intent(in)  :: ilow      ! lower variable
   integer(4),       intent(in)  :: iupp      ! upper variable
   type(scalar_var), intent(out) :: var(*)    ! scalar variable

   call ScalarSampleBlock(iStage, ilow, iupp, var)
   call ScalarSampleExtrap(iStage, ilow, iupp, var)

end subroutine ScalarSample

!************************************************************************
!> \page statistics statistics.F90
!! **ScalarSampleBlock**
!! *sample scalar quantities; fixed block length*
!************************************************************************


subroutine ScalarSampleBlock(iStage, ilow, iupp, var)

   use StatisticsModule
   use MollibModule, only: InvInt
   implicit none

   real(8),    parameter :: Zero = 0.0d0, One = 1.0d0
   integer(4),       intent(in)  :: iStage
   integer(4),       intent(in)  :: ilow      ! lower variable
   integer(4),       intent(in)  :: iupp      ! upper variable
   type(scalar_var), intent(out) :: var(*)    ! scalar variable
   integer(4) :: i
   real(8)    :: norm, norm1

   select case (iStage)
   case (3)  ! before simulation

      do i = ilow, iupp
         var(i)%nsamp1 = 0                                                      ! initiate nsamp1
         var(i)%avs1 = Zero                                                     ! initiate avs1
         var(i)%avsd = Zero                                                     ! initiate avsd
         var(i)%fls1 = Zero                                                     ! initiate fls1
         var(i)%flsd = Zero                                                     ! initiate flsd
      end do

   case (4)  ! before a macrostep

      do i = ilow, iupp
         var(i)%nsamp2 = 0                                                      ! initiate nsamp2
         var(i)%avs2 = Zero                                                     ! initiate avs2
         var(i)%fls2 = Zero                                                     ! initiate fls2
      end do

   case (5)  ! one simulation time step/pass

      do i = ilow, iupp
         var(i)%nsamp2 = var(i)%nsamp2 + 1                                      ! sample nsamp2
! if(i == 1 .or. i == 7) write(*,'(a,2i5,2f10.5)') 'iStage,i,var(i)%value',iStage,i,sngl(var(i)%value)
         var(i)%avs2   = var(i)%avs2   + var(i)%value                           ! sample avs2
         var(i)%fls2   = var(i)%fls2   + var(i)%value**2                        ! sample fls2
      end do

   case (6)  ! after a macrostep

      do i = ilow, iupp
         norm = InvInt(var(i)%nsamp2)
         var(i)%nsamp1 = var(i)%nsamp1 + 1                                      ! sample nsamp1
         var(i)%avs2   = var(i)%avs2*norm                                       ! form average avs2
         var(i)%avs1   = var(i)%avs1 + var(i)%avs2                              ! sample avs1
         var(i)%avsd   = var(i)%avsd + var(i)%avs2**2                           ! sample avsd
         var(i)%fls2   = sqrt(max(Zero, var(i)%fls2*norm - var(i)%avs2**2))     ! sample fls2
         var(i)%fls1   = var(i)%fls1 + var(i)%fls2                              ! sample fls1
         var(i)%flsd   = var(i)%flsd + var(i)%fls2**2                           ! sample flsd
      end do

   case (7)  ! after simulation

      do i = ilow, iupp
         norm  = InvInt(var(i)%nsamp1)
         norm1 = InvInt(var(i)%nsamp1-1)
         var(i)%avs1 = var(i)%avs1*norm                                         ! form average avs1
         var(i)%avsd = sqrt(max(Zero,(var(i)%avsd*norm-var(i)%avs1**2)*norm1))  ! form average avsd
         var(i)%fls1 = var(i)%fls1*norm                                         ! form average fls1
         var(i)%flsd = sqrt(max(Zero,(var(i)%flsd*norm-var(i)%fls1**2)*norm1))  ! form average flsd
      end do

   end select

end subroutine ScalarSampleBlock

!************************************************************************
!> \page statistics statistics.F90
!! **ScalarSampleExtrap**
!! *sample scalar quantities; extrapolate to large blocklen*
!************************************************************************


subroutine ScalarSampleExtrap(iStage, ilow, iupp, var)

   use StatisticsModule
   use MollibModule, only: InvInt
   implicit none

   real(8),    parameter :: Zero = 0.0d0, One = 1.0d0
   integer(4),       intent(in)  :: iStage
   integer(4),       intent(in)  :: ilow      ! lower variable
   integer(4),       intent(in)  :: iupp      ! upper variable
   type(scalar_var), intent(out) :: var(*)    ! scalar variable
   integer(4), save :: blocklen(mnblocklen)
   real(8),    save :: blockleni(mnblocklen)
   integer(4) :: i, ibl, nblocklen
   real(8)    :: InvFlt, norm, norm1
   real(8)    :: xfit(50), yfit_av(50), yfit_fl(50), wfit(50), afit_av(0:2), afit_fl(0:2), PolVal, dum1, dum2

   select case (iStage)
   case (3)  ! before simulation

      blocklen = [ (2**(ibl-1), ibl = 1,mnblocklen) ]                           ! initiate blocklen
      blockleni = One/blocklen                                                  ! initiate blockleni

      do i = ilow, iupp
         var(i)%nsamp = 0
         var(i)%nblock = 0                                                      ! initiate nblock
         var(i)%av_s1 = Zero                                                    ! initiate av_s1
         var(i)%av_sd = Zero                                                    ! initiate av_sd
         var(i)%av_s2 = Zero                                                    ! initiate av_s2
         var(i)%fl_s1 = Zero                                                    ! initiate fl_s1
         var(i)%fl_sd = Zero                                                    ! initiate fl_sd
         var(i)%fl_s2 = Zero                                                    ! initiate fl_s2
      end do

   case (5)  ! one simulation time step/pass

      do i = ilow, iupp
         var(i)%nsamp = var(i)%nsamp + 1
         do ibl = 1, mnblocklen                                                 ! loop over block lenghts
            var(i)%av_s2(ibl) = var(i)%av_s2(ibl) + var(i)%value                ! add data to av_s2
            var(i)%fl_s2(ibl) = var(i)%fl_s2(ibl) + var(i)%value**2             ! add data to fl_s2
            if (mod(var(i)%nsamp,blocklen(ibl)) == 0) then                      ! check for end of block length
               var(i)%nblock(ibl) = var(i)%nblock(ibl) + 1                      ! update number of blocks
               var(i)%av_s2(ibl) = var(i)%av_s2(ibl) * blockleni(ibl)           ! get average of the block
               var(i)%av_s1(ibl) = var(i)%av_s1(ibl) + var(i)%av_s2(ibl)        ! sum up av_sd for run average
               var(i)%av_sd(ibl) = var(i)%av_sd(ibl) + var(i)%av_s2(ibl)**2     ! sum up av_s2 precision of run average
               var(i)%fl_s2(ibl) = sqrt(max(Zero, var(i)%fl_s2(ibl) * blockleni(ibl) - var(i)%av_s2(ibl)**2)) ! get average of the block
               var(i)%fl_s1(ibl) = var(i)%fl_s1(ibl) + var(i)%fl_s2(ibl)        ! sum up fl_sd for run average
               var(i)%fl_sd(ibl) = var(i)%fl_sd(ibl) + var(i)%fl_s2(ibl)**2     ! sum up fl_s2 precision of run average
               var(i)%av_s2(ibl) = Zero                                         ! initialize av_s2 for next block
               var(i)%fl_s2(ibl) = Zero                                         ! initialize fl_s2 for next block
            end if
         end do
      end do

   case (7)  ! after simulation

      if (lblockaverwrite) call ScalarSampleWrite(1)  ! open unit and write head
      do i = ilow, iupp
         nblocklen = max(1,min(int(1+log(float(max(var(i)%nsamp,1)))/log(2.0d0)),mnblocklen)) ! number of block lenghts
         var(i)%nblocklen = nblocklen
         if ((nblocklen > 1) .and. minval(var(i)%nblock(1:nblocklen)) > 0) then ! extrapolate by linear fit of av_sd vs log(nblock) to nblock = 1
            do ibl = 1, nblocklen                              ! make block averages for varying block length
               norm = InvInt(var(i)%nblock(ibl))
               norm1 = InvInt(var(i)%nblock(ibl)-1)
               xfit(ibl) = log10(real(var(i)%nblock(ibl)))                                            ! save x for linear fit
               wfit(ibl) = sqrt(real(var(i)%nblock(ibl)))                                             ! save w for linear fit
               var(i)%av_s1(ibl) = var(i)%av_s1(ibl)*norm                                             ! average
               var(i)%av_sd(ibl) = sqrt(max(Zero,var(i)%av_sd(ibl)*norm-var(i)%av_s1(ibl)**2)*norm1)  ! its sd
               yfit_av(ibl) = var(i)%av_sd(ibl)                                                       ! save average for linear fit
               var(i)%fl_s1(ibl) = var(i)%fl_s1(ibl)*norm                                             ! fluctuation
               var(i)%fl_sd(ibl) = sqrt(max(Zero,var(i)%fl_sd(ibl)*norm-var(i)%fl_s1(ibl)**2)*norm1)  ! its sd
               yfit_fl(ibl) = var(i)%fl_sd(ibl)                                                       ! save fluctuation for linear fit
            end do
            call PolFit(1, nblocklen, xfit, yfit_av, wfit, 0, 6, afit_av, dum1, dum2)                 ! fit
            call PolFit(1, nblocklen, xfit, yfit_fl, wfit, 0, 6, afit_fl, dum1, dum2)                 ! fit
            var(i)%av_sd_extrap = max(Zero,PolVal(1,afit_av,0.0d0))                                   ! extrapolated sd of average
            var(i)%fl_sd_extrap = max(Zero,PolVal(1,afit_fl,0.0d0))                                   ! extrapolated sd of fluctuation
            var(i)%av_sd_stateff = (var(i)%nblock(1)*InvInt(var(i)%nblock(1)-1))*(var(i)%av_sd_extrap*InvFlt(var(i)%av_sd(1)))**2
            var(i)%fl_sd_stateff = (var(i)%nblock(1)*InvInt(var(i)%nblock(1)-1))*(var(i)%fl_sd_extrap*InvFlt(var(i)%fl_sd(1)))**2
         else
            var(i)%av_sd_extrap = Zero
            var(i)%fl_sd_extrap = Zero
            var(i)%av_sd_stateff = Zero
            var(i)%fl_sd_stateff = Zero
         end if
         if (lblockaverwrite) call ScalarSampleWrite(2)
      end do
      if (lblockaverwrite) call ScalarSampleWrite(3)

   end select

contains

!........................................................................

subroutine ScalarSampleWrite(iStage)

   integer(4), intent(in) :: iStage
   integer(4), parameter :: unit = 50
   logical, save :: first = .true.

   select case (iStage)
   case (1)                          ! open file for blockaverage data
      if (first) then
         open(unit,file = 'blockaver.data', status = 'unknown')
         first = .false.
         write(unit,'(a)') 'results from block averaging with variable block length'
         write(unit,'(a)') '-------------------------------------------------------'
      else
         open(unit,file = 'blockaver.data', status = 'unknown', position = 'append')
      end if
   case (2)
      norm=var(i)%norm               ! output blockaverage data
      write(unit,'()')
      write(unit,'(a,a)') 'Variable: ', var(i)%label(1:len_trim(var(i)%label)-1)
      write(unit,'()')
      write(unit,'(t5,a,t19,a,t36,a,t52,a,t71,a,t86,a,t101,a,t115,a,t125,a)') 'nblock','blocklen','product','average','sd','fit','average','sd','fit'
      write(unit,'(t5,a,t19,a,t36,a,t52,a,t71,a,t86,a,t101,a,t115,a,t125,a)') '------','--------','-------','-------','--','---','-------','--','---'
      write(unit,'(3(i10,a),g14.7,a,g12.5,a,g12.5,a,g14.7,a,g12.5,a,g12.5)') &
        (var(i)%nblock(ibl), char(9), blocklen(ibl), char(9), var(i)%nblock(ibl)*blocklen(ibl), char(9), &
        var(i)%av_s1(ibl)*norm, char(9), var(i)%av_sd(ibl)*norm, char(9), PolVal(1,afit_av,xfit(ibl))*norm, char(9), &
        var(i)%fl_s1(ibl)*norm, char(9), var(i)%fl_sd(ibl)*norm, char(9), PolVal(1,afit_fl,xfit(ibl))*norm, ibl = 1,var(i)%nblocklen)
      write(unit,'()')
      write(unit,'(a,t35,g12.5,t90,g12.5)') 'average value            = ', var(i)%av_s1(1)*norm,                                var(i)%fl_s1(1)*norm
      write(unit,'(a,t35,g12.5,t90,g12.5)') 'fluctuation              = ', var(i)%av_sd(1)*sqrt(real(var(i)%nblock(1)-1))*norm, var(i)%fl_sd(1)*sqrt(real(var(i)%nblock(1)-1))*norm
      write(unit,'(a,t35,g12.5,t90,g12.5)') 'extrapolated uncertainty = ', var(i)%av_sd_extrap*norm,                            var(i)%fl_sd_extrap*norm
      write(unit,'(a,t35,g12.5,t90,g12.5)') 'statistical inefficiency = ', var(i)%av_sd_stateff,                                var(i)%fl_sd_stateff
   case (3)
      close(unit)
   end select

end subroutine ScalarSampleWrite

!........................................................................

end subroutine ScalarSampleExtrap

!************************************************************************
!> \page statistics statistics.F90
!! **ScalarNorm**
!! *normalize scalar quantities*
!************************************************************************


subroutine ScalarNorm(iStage, ilow, iupp, var, iopt)

   use StatisticsModule
   implicit none

   integer(4),       intent(in)    :: iStage
   integer(4),       intent(in)    :: ilow      ! lower variable
   integer(4),       intent(in)    :: iupp      ! upper variable
   type(scalar_var), intent(inout) :: var(*)    ! scalar variable
   integer(4),       intent(in)    :: iopt      ! = 0, standard (multiplication by norm)
                                                ! = 1, fluctuation quantities are normalized with sqrt(norm)
                                                ! = 2, for rms quantities

   integer(4) :: i, nblocklen

   select case (iStage)
   case (6)  ! after a macrostep

      if (iopt == 0) then
         do i = ilow, iupp
            var(i)%avs2 = var(i)%avs2*var(i)%norm
            var(i)%fls2 = var(i)%fls2*var(i)%norm
  !          if(i == 1 .or. i == 7) write(*,'(a,3i5,2f10.5)') 'iStage,iopt,i,%avs2, fls2',iStage,iopt,i,sngl(var(i)%avs2), sngl(var(i)%fls2)
         end do
      else if (iopt == 1) then
         do i = ilow, iupp
            var(i)%avs2 = var(i)%avs2*var(i)%norm
            var(i)%fls2 = sqrt(var(i)%fls2*var(i)%norm)   ! sqrt(...)
         end do
      else if (iopt == 2) then
         do i = ilow, iupp
   !         if(i == 1 .or. i == 7) write(*,'(a,3i5,2f10.5)') '1 iStage,iopt,i,%avs2, fls2',iStage,iopt,i,sngl(var(i)%avs2), sngl(var(i)%fls2)
            var(i)%avs2 = sqrt(var(i)%avs2*var(i)%norm)   ! sqrt(...)
            var(i)%fls2 = sqrt(var(i)%fls2*var(i)%norm)   ! sqrt(...)
    !        if(i == 1 .or. i == 7) write(*,'(a,3i5,2f10.5)') '2 iStage,iopt,i,%avs2, fls2',iStage,iopt,i,sngl(var(i)%avs2), sngl(var(i)%fls2)
         end do
      end if

   case (7)  ! after simulation

! ... fixed blocklen

       if (iopt == 0) then
          do i = ilow, iupp
             var(i)%avs1 = var(i)%avs1*var(i)%norm
             var(i)%avsd = var(i)%avsd*var(i)%norm
             var(i)%fls1 = var(i)%fls1*var(i)%norm
             var(i)%flsd = var(i)%flsd*var(i)%norm
!             if(i == 1 .or. i == 7) then
!                write(*,'(a,3i5,2f10.5)') 'iStage,iopt,i,%avs1,%avsd',iStage,iopt,i,sngl(var(i)%avs1),sngl(var(i)%avsd)
!                write(*,'(a,3i5,2f10.5)') 'iStage,iopt,i,%fls1,%flsd',iStage,iopt,i,sngl(var(i)%fls1),sngl(var(i)%flsd)
!             end if
          end do
       else if (iopt == 1) then
          do i = ilow, iupp
             var(i)%avs1 = var(i)%avs1*var(i)%norm
             var(i)%avsd = var(i)%avsd*var(i)%norm
             var(i)%fls1 = var(i)%fls1*sqrt(var(i)%norm)  ! sqrt(norm)
             var(i)%flsd = var(i)%flsd*sqrt(var(i)%norm)  ! sqrt(norm)
          end do
       else if (iopt == 2) then
          do i = ilow, iupp
!             if(i == 1 .or. i == 7) then
!                write(*,'(a,3i5,2f10.5)') '1 iStage,iopt,i,%avs1,%avsd',iStage,iopt,i,sngl(var(i)%avs1),sngl(var(i)%avsd)
!                write(*,'(a,3i5,2f10.5)') '1 iStage,iopt,i,%fls1,%flsd',iStage,iopt,i,sngl(var(i)%fls1),sngl(var(i)%flsd)
!             end if
             var(i)%avs1 = sqrt(var(i)%avs1*var(i)%norm)  ! sqrt(...)
             var(i)%avsd = sqrt(var(i)%avs1**2+var(i)%avsd*var(i)%norm) - var(i)%avs1  ! sd of sqrt of average
             var(i)%fls1 = sqrt(var(i)%fls1*var(i)%norm)  ! sqrt(...)
             var(i)%flsd = sqrt(var(i)%fls1**2+var(i)%flsd*var(i)%norm) - var(i)%fls1  ! sd of sqrt of fluktuation
!             if(i == 1 .or. i == 7) then
!                write(*,'(a,3i5,2f10.5)') '2 iStage,iopt,i,%avs1,%avsd',iStage,iopt,i,sngl(var(i)%avs1),sngl(var(i)%avsd)
!                write(*,'(a,3i5,2f10.5)') '2 iStage,iopt,i,%fls1,%flsd',iStage,iopt,i,sngl(var(i)%fls1),sngl(var(i)%flsd)
!             end if
          end do
       end if

! ... variable blocklen

       if (iopt == 0) then
          do i = ilow, iupp
             var(i)%av_s1 = var(i)%av_s1*var(i)%norm
             var(i)%av_sd = var(i)%av_sd*var(i)%norm
             var(i)%av_sd_extrap = var(i)%av_sd_extrap*var(i)%norm
             var(i)%fl_s1 = var(i)%fl_s1*var(i)%norm
             var(i)%fl_sd = var(i)%fl_sd*var(i)%norm
          end do
       else if (iopt == 1) then
          do i = ilow, iupp
             var(i)%av_s1 = var(i)%av_s1*var(i)%norm
             var(i)%av_sd = var(i)%av_sd*var(i)%norm
             var(i)%av_sd_extrap = var(i)%av_sd_extrap*var(i)%norm
             var(i)%fl_s1 = var(i)%fl_s1*sqrt(var(i)%norm)
             var(i)%fl_sd = var(i)%fl_sd*sqrt(var(i)%norm)
          end do
       else if (iopt == 2) then
          do i = ilow, iupp
             nblocklen = var(i)%nblocklen
             var(i)%av_s1 = sqrt(var(i)%av_s1*var(i)%norm)  ! sqrt(...)
             var(i)%av_sd = sqrt(var(i)%av_s1**2+var(i)%av_sd*var(i)%norm) - var(i)%av_s1  ! sd of sqrt of average
             var(i)%av_sd_extrap = sqrt(var(i)%av_s1(1)**2+var(i)%av_sd_extrap*var(i)%norm) - var(i)%av_s1(1)  ! sd of sqrt of average (extrapolated value)
             var(i)%fl_s1 = sqrt(var(i)%fl_s1*var(i)%norm)  ! sqrt(...)
             var(i)%fl_sd = sqrt(var(i)%fl_s1**2+var(i)%fl_sd*var(i)%norm) - var(i)%fl_s1  ! sd of sqrt of flucutation
             var(i)%fl_sd_extrap = sqrt(var(i)%fl_s1(nblocklen)**2+var(i)%fl_sd_extrap*var(i)%norm) - var(i)%fl_s1(nblocklen)  ! sd of sqrt of flucutation (extrapolated value)
          end do
       end if

   end select

end subroutine ScalarNorm

!************************************************************************
!> \page statistics statistics.F90
!! **ScalarWrite**
!! *write scalar quantities*
!************************************************************************


subroutine ScalarWrite(iStage, ilow, iupp, var, iopt, fmt, unit)

   use StatisticsModule
   implicit none

   integer(4),       intent(in) :: iStage
   integer(4),       intent(in) :: ilow      ! lower variable
   integer(4),       intent(in) :: iupp      ! upper variable
   type(scalar_var), intent(in) :: var(*)    ! scalar variable
   integer(4),       intent(in) :: iopt      ! =0, averages are written
                                             ! =1, averages and fluctuations are written
   character(*),     intent(in) :: fmt       ! format
   integer(4),       intent(in) :: unit      ! output unit number
   integer(4) :: i, nblocklen
   character(40), parameter :: format = '(a,t43,a,t55,a,t69,a,t86,a,t102,a)'

   select case (iStage)
   case (6)  ! after a macrostep

      if (iopt == 0) then
         write(unit,format) 'quantity', 'average'
         write(unit,format) '--------', '-------'
         do i = ilow, iupp
            write(unit,fmt) trim(var(i)%label), var(i)%avs2
         end do
      else if (iopt == 1) then
         write(unit,format) 'quantity', 'average', 'fluctuation'
         write(unit,format) '--------', '-------', '-----------'
         do i = ilow, iupp
            write(unit,fmt) trim(var(i)%label), var(i)%avs2, var(i)%fls2
         end do
      end if

   case (7)  ! after simulation

      if (iopt == 0) then
         write(unit,format) 'quantity', 'average', '  precision'
         write(unit,format) '--------', '-------', '  ---------'
         do i = ilow, iupp
!            write(unit,fmt) trim(var(i)%label), var(i)%avs1, var(i)%avsd              ! fixed blocklen
             write(unit,fmt) trim(var(i)%label), var(i)%av_s1(1), var(i)%av_sd_extrap  ! variable blocklen
!            write(unit,fmt) trim(var(i)%label), var(i)%avs1, var(i)%av_sd_extrap      ! "hybrid" blocklen
         end do
      else if (iopt == 1) then
         write(unit,format) 'quantity', 'average', '  precision', 'fluctuation', 'precision', 'stat eff'
         write(unit,format) '--------', '-------', '  ---------', '-----------', '---------', '--------'
         do i = ilow, iupp
            nblocklen = var(i)%nblocklen
!           write(unit,fmt) trim(var(i)%label), var(i)%avs1, var(i)%avsd, var(i)%fls1, var(i)%flsd ! fixed blocklen
            write(unit,fmt) trim(var(i)%label), var(i)%av_s1(1), var(i)%av_sd_extrap, var(i)%fl_s1(nblocklen), &
                                                var(i)%fl_sd_extrap,var(i)%av_sd_stateff           ! variable blocklen
!           write(uout,fmt) trim(var(i)%label), var(i)%avs1, var(i)%av_sd_extrap, var(i)%fl_s1(nblocklen), &
!                                                var(i)%fl_sd_extrap, var(i)%av_sd_stateff         ! "hybrid" blocklen
         end do
      end if

   end select

end subroutine ScalarWrite

! #if defined (_PAR_)
! !************************************************************************
!> \page statistics statistics.F90
!! **ScalarAllReduce**
!! *make all_reduce of scalar quantities*
! !************************************************************************


! subroutine ScalarAllReduce(iStage, ilow, iupp, var, raux)

!    use StatisticsModule
!    implicit none

!    integer(4),       intent(in) :: iStage
!    integer(4),       intent(in) :: ilow           ! lower variable
!    integer(4),       intent(in) :: iupp           ! upper variable
!    type(scalar_var), intent(inout) :: var(*)      ! scalar variable
!    real(8),          intent(inout) :: raux(*)     ! temporary variable

!    integer(4) :: i

!    select case (iStage)
!    case (6)  ! after a macrostep

!       do i = ilow, iupp
!          !TODO: par_all_reduce_ints requires an integer as the second parameter
!          call par_allreduce_ints(var(i)%nsamp2, raux, 1)
!          call par_allreduce_reals(var(i)%avs2, raux, 1)
!          call par_allreduce_reals(var(i)%fls2, raux, 1)
!       end do

!    case (7)  ! after simulation

!       do i = ilow, iupp
!          write(*,*) 'i, var(i)%nsamp2', i, var(i)%nsamp2
!          call stop('ScalarAllReduce', 'atempt to allreduce variable-blocklengh data', 6)
! !         call par_allreduce_ints(var(i)%nsamp2, raux, 1)
! !         call par_allreduce_reals(var(i)%av_s1, raux, 1)
! !         call par_allreduce_reals(var(i)%av_sd, raux, 1)
!       end do

!    end select

! end subroutine ScalarAllReduce
! #endif

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncSample**
!! *sample distribution functions*
!************************************************************************


subroutine DistFuncSample(iStage, nvar, var)

   use StatisticsModule
   use MollibModule, only: InvInt
   implicit none

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0
   integer(4),   intent(in)    :: iStage
   integer(4),   intent(in)    :: nvar     ! number of distribution functions
   type(df_var), intent(inout) :: var(*)   ! distribution functions
   integer(4) :: i, ibin
   real(8)    :: norm, norm1

   select case (iStage)
   case (2)  ! read input

      do i = 1, nvar
         var(i)%bin = (var(i)%max - var(i)%min) / var(i)%nbin          ! initiate bin
         var(i)%bini = One/var(i)%bin                                  ! initiate bini
      end do

   case (3)  ! before simulation

      do i = 1, nvar
         var(i)%nsamp1 = 0                                             ! initiate nsamp1
         var(i)%avs1(-1:var(i)%nbin) = Zero                            ! initiate avs1
         var(i)%avsd(-1:var(i)%nbin) = Zero                            ! initiate avsd
      end do

   case (4)  ! before a macrostep

      do i = 1, nvar
         var(i)%nsamp2 = 0                                             ! initiate nsamp2
         var(i)%avs2(-1:var(i)%nbin) = Zero                            ! initiate avs2
         var(i)%nsampbin(-1:var(i)%nbin) = Zero                        ! initiate nsampbin
      end do

   case (6)  ! after a macrostep

      do i = 1, nvar
         var(i)%nsamp1 = var(i)%nsamp1 + 1                             ! update nsamp1
         norm = InvInt(var(i)%nsamp2)
         do ibin = -1, var(i)%nbin
            var(i)%avs2(ibin) = var(i)%avs2(ibin)*norm                 ! form average avs2
            var(i)%avs1(ibin) = var(i)%avs1(ibin)+var(i)%avs2(ibin)    ! sum up avs1
            var(i)%avsd(ibin) = var(i)%avsd(ibin)+var(i)%avs2(ibin)**2 ! sum up avsd
         end do
      end do

   case (7)  ! after simulation

      do i = 1, nvar
         norm = InvInt(var(i)%nsamp1)
         norm1 = InvInt(var(i)%nsamp1-1)
         do ibin = -1, var(i)%nbin
            var(i)%avs1(ibin) = var(i)%avs1(ibin)*norm                 ! form average avs1
            var(i)%avsd(ibin) = sqrt(max(Zero,var(i)%avsd(ibin)*norm-var(i)%avs1(ibin)**2)*norm1) ! form sd
         end do
      end do

   end select

end subroutine DistFuncSample

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncNorm**
!! *normalize distribution functions*
!************************************************************************


!     bin * sum(avs2) = nsamp2

subroutine DistFuncNorm(ilow, iupp, var)

   use StatisticsModule
   implicit none

   real(8), parameter :: Zero = 0.0d0
   integer(4),   intent(in)    :: ilow     ! lower distribution functions
   integer(4),   intent(in)    :: iupp     ! upper distribution functions
   type(df_var), intent(inout) :: var(*)   ! distribution functions
   integer(4) :: i
   real(8)    :: fac, norm

   if ((ilow <= 0) .or. (iupp <= 0)) return

   do i = ilow, iupp
      fac = sum(var(i)%avs2(-1:var(i)%nbin))
      norm = Zero
      if (fac/= Zero) norm = real(var(i)%nsamp2)/(var(i)%bin*fac)
      var(i)%avs2(-1:var(i)%nbin) = var(i)%avs2(-1:var(i)%nbin)*norm
   end do

end subroutine DistFuncNorm

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncHead**
!! *write heading of distribution functions*
!************************************************************************


subroutine DistFuncHead(nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: nvar        ! number of distribution functions
   type(df_var), intent(in) :: var(*)      ! distribution functions
   integer(4),   intent(in) :: unit        ! output unit number
   integer(4) :: i

   write(unit,'()')
   write(unit,'(a,t7,a,t40,a,t55,a,t70,a,t85,a)') 'no', 'function', 'min', 'max', 'nbin', 'step'
   write(unit,'(a,t7,a,t40,a,t55,a,t70,a,t85,a)') '--', '--------', '---', '---', '----', '----'
   write(unit,'(i2,t7,a,t35,f10.3,t50,f10.3,t65,i9,t80,f10.3)')  &
                (i,var(i)%label,var(i)%min,var(i)%max,var(i)%nbin,var(i)%bin,i = 1,nvar)

end subroutine DistFuncHead

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncWrite**
!! *control the calls of DistFuncShow, DistFuncPlot, and DistFuncList*
!************************************************************************


subroutine DistFuncWrite(txheading, nvar, var, lout, llist, ishow, iplot, ilist)

   use StatisticsModule
   implicit none

   character(*), intent(in) :: txheading        ! heading
   integer(4),   intent(in) :: nvar             ! number of distribution functions
   type(df_var), intent(in) :: var(*)           ! distribution functions
   integer(4),   intent(in) :: lout             ! fout unit number
   integer(4),   intent(in) :: llist            ! flist unit number
   integer(4),   intent(in) :: ishow            ! controls listing on fout
   integer(4),   intent(in) :: iplot            ! controls plotint on fout
   integer(4),   intent(in) :: ilist            ! controls listing on flist

   if (ishow/= 0) call DistFuncShow(ishow, txheading, nvar, var, lout)
   if (iplot/= 0) call DistFuncPlot(nvar, var, lout)
   if (ilist/= 0) call DistFuncList(ilist, txheading, nvar, var, llist)

end subroutine DistFuncWrite

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncShow**
!! *list bin numbers, x, y, and dy of distribution functions*
!************************************************************************


subroutine DistFuncShow(il, txheading, nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: il                 ! write every il:th bin
   character(*), intent(in) :: txheading          ! heading
   integer(4),   intent(in) :: nvar               ! number of distribution functions
   type(df_var), intent(in) :: var(*)             ! distribution functions
   integer(4),   intent(in) :: unit               ! output unit number
   integer(4), save :: ncol = 4
   integer(4) :: mone, mlow, mupp, i, ibin, nbin

   mone =-1
   write(unit,'()')
   write(unit,'(a)') txheading
   write(unit,'(100a)') ('-',i = 1,len(trim(txheading)))
   do mlow = 1, nvar, ncol
      mupp = min(int(nvar),int(mlow+ncol-1))
      write(unit,'()')
      write(unit,'(a,t10,a,t38,a,t66,a,t94,a)') '   ', var(mlow:mupp)%label
      write(unit,'(a,t10,4(a,7x,a,7x,a,10x))') ' bin', ('x','y','dy',i = mlow,mupp)
      write(unit,'(i4,4(3x,a,f10.4,1x,f7.4,2x))') mone, ('below',var(i)%avs1(-1),var(i)%avsd(-1),i = mlow,mupp)
      nbin = var(mlow)%nbin
      do ibin = 0, nbin-1, il
         write(unit,'(i4,4(f8.2,f10.4,1x,f7.4,2x))') &
         ibin, (var(i)%min+(ibin+0.5)*var(i)%bin,var(i)%avs1(ibin),var(i)%avsd(ibin),i = mlow,mupp)
      end do
      write(unit,'(i4,4(3x,a,f10.4,1x,f7.4,2x))') nbin, ('above',var(i)%avs1(nbin),var(i)%avsd(nbin),i = mlow,mupp)
   end do

end subroutine DistFuncShow

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncPlot**
!! *plot distribution functions*
!************************************************************************


subroutine DistFuncPlot(nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: nvar               ! number of distribution functions
   type(df_var), intent(in) :: var(*)             ! distribution functions
   integer(4),   intent(in) :: unit               ! output unit number
   integer(4) :: i
   real(8)    :: vl, vu, dum

   do i = 1, nvar
      vl = var(i)%min
      vu = var(i)%min + real(var(i)%nbin) * var(i)%bin
      call Plot(var(i)%label, var(i)%nbin, var(i)%avs1(0), 'none', vl, vu, dum, dum, unit)
   end do

end subroutine DistFuncPlot

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncList**
!! *list x, y, and dy values of distribution functions*
!************************************************************************


subroutine DistFuncList(il, txheading, nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: il                 ! write every il:th bin
   character(*), intent(in) :: txheading          ! heading
   integer(4),   intent(in) :: nvar               ! number of distribution functions
   type(df_var), intent(in) :: var(*)             ! distribution functions
   integer(4),   intent(in) :: unit               ! output unit number
   integer(4) :: i, ibin

   write(unit,'(a)') txheading
   write(unit,'(i5)') nvar
   do i = 1, nvar
      write(unit,'(a)') var(i)%label
      write(unit,'(i5)') 1+(var(i)%nbin-1)/il
      write(unit,'(g15.5,a,g15.5,a,g15.5)') &
      (var(i)%min+(0.5+ibin)*var(i)%bin, char(9), var(i)%avs1(ibin), char(9), var(i)%avsd(ibin), ibin = 0,var(i)%nbin-1,il)
   end do

end subroutine DistFuncList

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncAverValue**
!! *calculate and write average of distribution functions*
!************************************************************************


subroutine DistFuncAverValue(nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: nvar              ! number of distribution functions
   type(df_var), intent(in) :: var(*)            ! distribution functions
   integer(4),   intent(in) :: unit              ! output unit number
   integer(4) :: i, ibin
   real(8)    :: value, aver, aver2

   write(unit,'()')
   write(unit,'(a,t7,a,t35,a,t50,a,t65,a,t80,a)') 'no', 'function', 'average', 'rms'
   write(unit,'(a,t7,a,t35,a,t50,a,t65,a,t80,a)') '--', '--------', '-------', '---'
   do i = 1, nvar
      aver = 0.0d0
      aver2 = 0.0d0
      do ibin =-1, var(i)%nbin
         value = var(i)%min+(ibin+0.5)*var(i)%bin
         aver = aver + var(i)%avs1(ibin)*value*var(i)%bin
         aver2 = aver2 + var(i)%avs1(ibin)*value**2*var(i)%bin
      end do
      write(unit,'(i2,t7,a,t30,2f15.5)') i, var(i)%label, aver, sqrt(aver2)
   end do

end subroutine DistFuncAverValue

!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncAverDist**
!! *calculate average and spread among distribution functions*
!************************************************************************


subroutine DistFuncAverDist(nvar2, ilow, iupp, var, var2, var2_spread)

   use StatisticsModule
   use MollibModule, only: InvInt
   implicit none

   integer(4),   intent(in)  :: nvar2                 ! number of distribution functions
   integer(4),   intent(in)  :: ilow(*)               ! lower distribution functions
   integer(4),   intent(in)  :: iupp(*)               ! upper distribution functions
   type(df_var), intent(in)  :: var(*)                ! underlaying distribution functions
   type(df_var), intent(inout) :: var2(*)               ! average of var from ilow to iupp
   real(8)     , intent(out) :: var2_spread(*)        ! var%avsd averaged over 0 to nbin-1
   integer(4) :: i, ibin, ncount, il, iu

   do i = 1, nvar2

      il = ilow(i)
      iu = iupp(i)

! ... distribution function of the average

      do ibin =-1, var2(i)%nbin
         var2(i)%avs1(ibin) = sum(var(il:iu)%avs1(ibin)) * InvInt(iu-il+1)
      end do

! ... distribution function of the sd

      do ibin =-1, var2(i)%nbin
         var2(i)%avsd(ibin) = sum( (var(il:iu)%avs1(ibin)-var2(i)%avs1(ibin))**2 ) * InvInt((iu-il+1)-1)
      end do
      var2(i)%avsd(-1:var2(i)%nbin) = sqrt(var2(i)%avsd(-1:var2(i)%nbin))

! ... global spread among the distribution functions

      var2_spread(i) = 0.0d0
      ncount = 0
      do ibin = 0, var2(i)%nbin-1
         if (maxval(var(il:iu)%avs1(ibin)) > 0.0d0) then                     ! exclude all Zeros
            var2_spread(i) = var2_spread(i)+var2(i)%avsd(ibin)**2
            ncount = ncount+1
         end if
      end do
      var2_spread(i) = sqrt(var2_spread(i)*InvInt(ncount-1))

   end do

end subroutine DistFuncAverDist

#if defined (_PAR_)
!************************************************************************
!> \page statistics statistics.F90
!! **DistFuncAllReduce**
!! *make all_reduce of df quantities*
!************************************************************************


subroutine DistFuncAllReduce(iStage, ilow, iupp, var, raux)

   use StatisticsModule
   implicit none

   integer(4),   intent(in) :: iStage
   integer(4),   intent(in) :: ilow           ! lower variable
   integer(4),   intent(in) :: iupp           ! upper variable
   type(df_var), intent(inout) :: var(*)      ! scalar variable
   real(8),      intent(inout) :: raux(*)     ! temporary variable

   integer(4) :: i

   select case (iStage)
   case (6)  ! after a macrostep

      do i = ilow, iupp
         call par_allreduce_reals(var(i)%avs1(-1), raux, mnbin_df+2)
         call par_allreduce_reals(var(i)%avsd(-1), raux, mnbin_df+2)
      end do

   end select

end subroutine DistFuncAllReduce
#endif

!************************************************************************
!> \page statistics statistics.F90
!! **DistFunc2DSample**
!! *sample two-dimensional distribution functions*
!************************************************************************


subroutine DistFunc2DSample(iStage, nvar, var)

   use StatisticsModule
   use MollibModule, only: InvInt
   implicit none

   real(8), parameter :: Zero = 0.0d0, One = 1.0d0
   integer(4),     intent(in)    :: iStage
   integer(4),     intent(in)    :: nvar         ! number of 2d distribution functions
   type(df2d_var), intent(inout) :: var(*)       ! 2d distribution functions
   integer(4) :: i, ibin1, ibin2
   real(8)    :: norm, norm1

   select case (iStage)
   case (2)  ! read input

! ... initialize bin and bini

      do i = 1, nvar
         var(i)%bin = (var(i)%max - var(i)%min) / var(i)%nbin
         var(i)%bini = One/var(i)%bin
      end do

   case (3)  ! before simulation

! ... initialize nsamp1, avs1, and avsd

      do i = 1, nvar
         var(i)%nsamp1 = 0
         var(i)%avs1(-1:var(i)%nbin(1),-1:var(i)%nbin(2)) = Zero
         var(i)%avsd(-1:var(i)%nbin(1),-1:var(i)%nbin(2)) = Zero
      end do

   case (4)  ! before a macrostep

! ... initialize nsamp2 and avs2

      do i = 1, nvar
         var(i)%nsamp2 = 0
         var(i)%avs2(-1:var(i)%nbin(1),-1:var(i)%nbin(2)) = Zero
      end do

   case (6)  ! after a macrostep

! ... average avs2 and sum up avs1 and avsd

      do i = 1, nvar
         var(i)%nsamp1 = var(i)%nsamp1 + 1
         norm = InvInt(var(i)%nsamp2)
         do ibin2 = -1, var(i)%nbin(2)
            do ibin1 = -1, var(i)%nbin(1)
               var(i)%avs2(ibin1,ibin2) = var(i)%avs2(ibin1,ibin2)*norm
               var(i)%avs1(ibin1,ibin2) = var(i)%avs1(ibin1,ibin2)+var(i)%avs2(ibin1,ibin2)
               var(i)%avsd(ibin1,ibin2) = var(i)%avsd(ibin1,ibin2)+var(i)%avs2(ibin1,ibin2)**2
            end do
         end do
      end do

   case (7)  ! after simulation

! ... average avs1 and calculate avsd

      do i = 1, nvar
         norm = InvInt(var(i)%nsamp1)
         norm1 = InvInt(var(i)%nsamp1-1)
         do ibin2 = -1, var(i)%nbin(2)
            do ibin1 = -1, var(i)%nbin(1)
               var(i)%avs1(ibin1,ibin2) = var(i)%avs1(ibin1,ibin2)*norm
               var(i)%avsd(ibin1,ibin2) = sqrt(max(Zero,var(i)%avsd(ibin1,ibin2)*norm-var(i)%avs1(ibin1,ibin2)**2)*norm1)
            end do
         end do
      end do

   end select

end subroutine DistFunc2DSample

!************************************************************************
!> \page statistics statistics.F90
!! **DistFunc2DNorm**
!! *normalize two-dimensinal distribution functions*
!************************************************************************


!     bin(1) * bin(2) * sum(avs2) = nsamp2

subroutine DistFunc2DNorm(ilow, iupp, var)

   use StatisticsModule
   implicit none

   real(8), parameter :: Zero = 0.0d0
   integer(4),     intent(in)    :: ilow        ! lower 2d distribution functions
   integer(4),     intent(in)    :: iupp        ! upper 2d distribution functions
   type(df2d_var), intent(inout) :: var(*)      ! 2d distribution functions
   integer(4) :: i
   real(8)    :: fac, norm

   if ((ilow <= 0) .or. (iupp <= 0)) return

   do i = ilow, iupp
      fac = sum(var(i)%avs2(-1:var(i)%nbin(1),-1:var(i)%nbin(2)))
      norm = Zero
      if (fac/= Zero) norm = real(var(i)%nsamp2)/(var(i)%bin(1)*var(i)%bin(2)*fac)
      var(i)%avs2(-1:var(i)%nbin(1),-1:var(i)%nbin(2)) = var(i)%avs2(-1:var(i)%nbin(1),-1:var(i)%nbin(2))*norm
   end do

end subroutine DistFunc2DNorm

!************************************************************************
!> \page statistics statistics.F90
!! **DistFunc2DHead**
!! *write heading of two-dimensional distribution functions*
!************************************************************************


subroutine DistFunc2DHead(nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),     intent(in) :: nvar           ! number of distribution functions
   type(df2d_var), intent(in) :: var(*)         ! distribution functions
   integer(4),     intent(in) :: unit           ! output unit number
   integer(4) :: i, j

   write(unit,'()')
   write(unit,'(a,t7,a,t35,a,t45,a,t55,a,t65,a,t75,a,t85,a,t95,a,t105,a)') &
    'no', 'function', 'min1', 'max1', 'nbin1', 'step1', 'min2', 'max2', 'nbin2', 'step2'
   write(unit,'(a,t7,a,t35,a,t45,a,t55,a,t65,a,t75,a,t85,a,t95,a,t105,a)') &
    '--', '--------', '----', '----', '-----', '-----', '----', '----', '-----', '-----'
   write(unit,'(i2,t7,a,t30,f10.3,t40,f10.3,t50,i9,t60,f10.3,t70,f10.3,t80,f10.3,t90,i9,t100,f10.3)') &
   (i,var(i)%label,(var(i)%min(j),var(i)%max(j),var(i)%nbin(j),var(i)%bin(j),j = 1,2),i = 1,nvar)

end subroutine DistFunc2DHead

!************************************************************************
!> \page statistics statistics.F90
!! **DistFunc2DShow**
!! *list bin numbers and z as well as bin numbers and \ref dz of two-dimensional distribution functions*
!************************************************************************


subroutine DistFunc2DShow(il, txheading, nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),     intent(in) :: il                 ! write every il:th bin
   character(*),   intent(in) :: txheading          ! heading
   integer(4),     intent(in) :: nvar               ! number of distribution functions
   type(df2d_var), intent(in) :: var(*)             ! 2d distribution functions
   integer(4),     intent(in) :: unit               ! output unit number
   integer(4) ::   i, ibin1, ibin2

   write(unit,'()')
   write(unit,'(a)') txheading
   write(unit,'(100a)') ('-',i = 1,len(trim(txheading)))

   do i = 1, nvar
      write(unit,'(a)') var(i)%label
      write(unit,'((42i8))') (ibin2,ibin2 = -1,var(i)%nbin(2),il)
      do ibin1 = -1, var(i)%nbin(1), il
         write(unit,'(i3,(42f8.4))') ibin1, var(i)%avs1(ibin1,-1:var(i)%nbin(2):il)
      end do
      write(unit,'()')
      write(unit,'(a)') 's.d. of '//var(i)%label
      write(unit,'((42i8))') (ibin2,ibin2 = 0,var(i)%nbin(2)-1,il)
      do ibin1 = 0, var(i)%nbin(1)-1, il
         write(unit,'(i3,(42f8.4))') ibin1, var(i)%avsd(ibin1,0:var(i)%nbin(2)-1:il)
      end do
   end do

end subroutine DistFunc2DShow

!************************************************************************
!> \page statistics statistics.F90
!! **DistFunc2DList**
!! *list bin numbers and z as well as bin numbers and \ref dz of two-dimensional distribution functions*
!************************************************************************


subroutine DistFunc2DList(il, txheading, nvar, var, unit)

   use StatisticsModule
   implicit none

   integer(4),     intent(in) :: il                 ! write every il:th bin
   character(*),   intent(in) :: txheading          ! heading
   integer(4),     intent(in) :: nvar               ! number of distribution functions
   type(df2d_var), intent(in) :: var(*)             ! 2d distribution functions
   integer(4),     intent(in) :: unit               ! output unit number
   integer(4) ::   i, ibin1, ibin2, m
   character(5) :: txvar

   if(il <= 0) return

   txvar = '  '

   write(unit,'(a)') txheading
   write(unit,'(i5)') nvar
   do i = 1, nvar
      write(unit,'(a)') var(i)%label
      if (txvar == 'xy') then                                        ! x and y variables
         write(unit,'(i5)') 1+(var(i)%nbin(1)-1)/il
         write(unit,'(8x,(42f8.4))') (var(i)%min(2)+(0.5+ibin2)*var(i)%bin(2),ibin2 = 0,var(i)%nbin(2)-1,il)
         do ibin1 = 0, var(i)%nbin(1)-1, il
            write(unit,'(42f8.4)') var(i)%min(1)+(0.5+ibin1)*var(i)%bin(1), var(i)%avs1(ibin1,0:var(i)%nbin(2)-1:il)
         end do
      else if (txvar == 'number') then                               ! ibin(1) and ibin(2)
         write(unit,'(i5)') 1+(var(i)%nbin(1)-1)/il
         write(unit,'((42i8))') (ibin2,ibin2 = 0,var(i)%nbin(2)-1,il)
         do ibin1 = 0, var(i)%nbin(1)-1, il
            write(unit,'(i3,(42f8.4))') ibin1, var(i)%avs1(ibin1,0:var(i)%nbin(2)-1:il)
         end do
      else                                                           ! nothing
         write(unit,'(2(i5,2f10.4))') (1+(var(i)%nbin(m)-1)/il, var(i)%min(m), var(i)%max(m), m=1,2)
         do ibin1 = 0, var(i)%nbin(1)-1, il
            write(unit,'(42f8.4)') var(i)%avs1(ibin1,0:var(i)%nbin(2)-1:il)
         end do
      end if
   end do

end subroutine DistFunc2DList

