!************************************************************************
!*                                                                      *
!*     SSOModule                                                        *
!*                                                                      *
!************************************************************************

! ... Module for the SSO simulation

 module SSOModule
   implicit none
   private
   public SSOSetup, DoSSOUpdate

   type :: step
      integer(8)  :: n  !number of steps
      real(8)     :: d2  !squared displacement
      real(8)     :: d4 !displacement**4
   end type step

   !define what happens if two variables of type step are added 
   interface operator(+)
      module procedure stepadd
   end interface operator(+)

   type(step), allocatable :: tots(:)     !total steps
   type(step), allocatable :: ssos(:,:)     !sso steps

   real(8),       allocatable    :: invrsso(:)

   contains

      !define what happens if two variables of type step are added 
      pure function stepadd(s1, s2) result(add)
         type(step), intent(in) :: s1, s2
         type(step)  :: add
         add%n = s1%n + s2%n
         add%d2 = s1%d2 + s2%d2
         add%d4 = s1%d4 + s2%d4
      end function stepadd


      !************************************************************************
      !*                                                                      *
      !*     SSODriver                                                        *
      !*                                                                      *
      !************************************************************************

      ! ... Driver for the SSO simulation

      subroutine SSOSetup(iStage)

         use MolModule, only: ltrace, txstart
         use MolModule, only: iReadInput, iWriteInput, iBeforeSimulation, iBeforeMacrostep, iSimulationStep, iAfterMacrostep, iAfterSimulation
         use MolModule, only: uin, uout, ucnf
         use MolModule, only: Zero, Half, One, Two, Three, Five
         use MolModule, only: npt, nstep, istep
         use MolModule, only: lbcbox, lbcrd, lbcto, lbcsph, lbcell, lbccyl, boxlen, cellside, sphrad, ellrad, cylrad, cyllen
         use MCModule, only: curdtranpt, lssopt
         implicit none

         integer(4), intent(in) :: iStage

         character(40), parameter :: txroutine ='SSODriver'

         integer(4)  :: ipt
         integer(4)  :: maxbin, upperbin, lowerbin ! bin where the local mobility is found and upper and lower boundary
         integer(4)  :: ibin, ipart

         !input variables--------------------------------------------------------------------------
         real(8), allocatable, save :: dtransso(:) ! initial displacement parameters
         integer(4), save           :: nstepzero     ! length of first part of simulation where local mobility is measured
         integer(4), save           :: nssobin
         integer(4), save           :: nstepend
         logical, save              :: ltestsso
         real(8), allocatable, save :: maxdtransso(:)! maximal allowed displacement parameters
         real(8), allocatable, save :: dtranfac(:)   !increase in displacement parameter
         !-----------------------------------------------------------------------------------------

         type  :: ssopart
            real(8)     :: fac      !increment of part length
            integer(4)  :: nextstep !step at which next part starts
            integer(4)  :: i        !current part
            integer(4)  :: n        !number of parts
         end type ssopart
         type(ssopart), save  :: part

         type  :: ssoparam
            real(8)     :: used     ! used dtran
            real(8)     :: opt      ! dtran with the highest mobility
            real(8)     :: err   ! accuracy of opt
         end type ssoparam
         type(ssoparam), save, allocatable  :: param(:,:)

         type  :: mobility
            real(8)     :: val         !value
            real(8)     :: error       !error
            real(8)     :: smooth      !smooth
         end type mobility
         type(mobility), allocatable, save :: mob(:)


         namelist /nmlSPartSSO/ dtransso, nstepzero, nssobin, nstepend, ltestsso, maxdtransso, dtranfac

         if (ltrace) call WriteTrace(2, txroutine, iStage)

         select case (iStage)
         case (iReadInput)

            ! read in nmlSpartSSO------------------------------------------------------------------
            if (.not. allocated(dtransso)) allocate(dtransso(npt))
            if (.not. allocated(maxdtransso)) allocate(maxdtransso(npt))
            if (.not. allocated(dtranfac)) allocate(dtranfac(npt))

            dtransso  = One
            nstepzero   = ceiling(sqrt(real(nstep)))
            nstepend    = max(nstepzero, int(0.1*nstep)
            nssobin     = 20
            ltestsso    = .false.
            if (lbcbox) then
               maxdtransso = minval(boxlen(1:3))
            else if (lbcrd) then
               maxdtransso = cellside
            else if (lbcto) then
               maxdtransso = cellside
            else if (lbcsph) then
               maxdtransso = Two*sphrad
            else if (lbcell) then
               maxdtransso = Two*minval(ellrad(1:3))
            else if (lbccyl) then
               maxdtransso = min(Two*cylrad,cyllen)
            end if
            dtranfac = Two

            rewind(uin)
            read(uin,nmlSPartSSO)
            !--------------------------------------------------------------------------------------

         case (iWriteInput) 

            ! check conditions---------------------------------------------------------------------
            if (nstepzero + nstepend > nstep) then
               write(uout, *) "nstepzero: ", nstepzero, "; nstepend: ", nstepend, "; nstep: ", nstep
               call Stop(txroutine, 'nstepzero + nstepend > nstep', uout)
            end if
            if (nstepend < nstepzero) then
               write(uout, *) "nstepzero: ", nstepzero, "; nstepend: ", nstepend, "; nstep: ", nstep
               call Stop(txroutine, 'nstepend < nstepzero', uout)
            end if
            dtransso = abs(dtransso) !dtransso must always be positive
            !--------------------------------------------------------------------------------------

            ! calculate part lengths---------------------------------------------------------------
            if((nstepzero .le. 0) .or. (nstepend .le. 0) .or. (nstepzero + nstepend > nstep) ) then
               call Warn(txroutine, "stepzero and nstepend are wrong. Doing only one SSO part", uout)
               part%fac = One
               part%n = 1
            else if (nstepzero == nstepend) then
               part%fac = One
               part%n = nstep/nstepend
            else
               part%fac = -real(nstepzero - nstep)/real(nstep - nstepend)
               part%n = int(log(real(nstepend/real(nstepzero)))/log(part%fac))
               if(nstepzero * (One - part%fac**(part%n + 1))/(One - part%fac) < nstep) part%n = part%n + 1
            end if
            ! -------------------------------------------------------------------------------------

            ! allcoate-----------------------------------------------------------------------------
            if(.not. allocated(curdtranpt)) allocate(curdtranpt(npt))
            if(.not. allocated(invrsso)) allocate(invrsso(npt))

            if(.not. allocated(tots)) allocate(tots(npt))
            if(.not. allocated(ssos)) allocate(ssos(nssobin, npt))

            if(.not. allocated(mob)) allocate(mob(0:nssobin))
            if(.not.allocated(param)) allocate(param(npt,part%n))
            ! -------------------------------------------------------------------------------------



         case (iBeforeSimulation)

            ! initialize values--------------------------------------------------------------------
            if (txstart == 'continue') then
               read(ucnf) curdtranpt, part%i, ssos, tots, param
            else
               curdtranpt(1:npt) = dtransso(1:npt)
               ssos=step(0, Zero, Zero)
               tots=step(0, Zero, Zero)
               part%i = 1
               part%nextstep = nstepzero
               param = ssoparam(Zero, Zero, Zero)
            end if
            if(part%i == part%n) part%nextstep = nstep
            do ipt = 1, npt
               if (lssopt(ipt)) then
                  if (curdtranpt(ipt) > 0) then
                     invrsso(ipt) = real(nssobin) * Two / curdtranpt(ipt)
                  else
                     invrsso(ipt) = Zero
                  end if
               end if
            end do
            !--------------------------------------------------------------------------------------


         case (iBeforeMacrostep)                       ! ignore macrosteps

            continue

         case (iSimulationStep)

            if(istep == part%nextstep) then

               if(ltestsso) then
                  call WriteHead(2, 'SSO - Results of current part', uout)
               end if

               do ipt = 1, npt
                  if (lssopt(ipt)) then

                     param(ipt,part%i)%used = curdtranpt(ipt) !store used dtran

                     call CalcLocMob(ipt)                      !calculate local mobility

                     ! get bin with maximum mobility--------------------------------------------------
                     maxbin = 0
                     do ibin = 1, nssobin - 1
                        if(mob(ibin)%smooth > mob(ibin + 1)%smooth .and. mob(ibin)%smooth > mob(maxbin)%smooth ) then !maximum found if the nextbin has a lower mobility and the current bin is the highest
                           maxbin = ibin
                        end if
                     end do
                     if(maxbin == 0) then  ! no local maximum found, maximum at border
                        maxbin = nssobin
                     end if
                     !--------------------------------------------------------------------------------

                     ! get upper boundary (relative to maxbin)----------------------------------------
                     upperbin = 0
                     if(maxbin == nssobin) then
                        upperbin = ceiling((dtranfac(ipt) - One)*nssobin) ! use dtranfac to estimate upper boundary if maxbin is at border
                     else
                        do ibin = maxbin + 1 , nssobin !find position where locmob is significantly different from maximum
                           if ( (mob(maxbin)%smooth - mob(maxbin)%error) .ge. (mob(ibin)%smooth + mob(ibin)%error) ) then
                              upperbin = (ibin - maxbin) + 1
                              exit
                           end if
                        end do
                        if (upperbin == 0) then        ! no significantly different position found
                           upperbin = 1 + ceiling(Two*mob(maxbin)%error*(nssobin - maxbin)/real((mob(maxbin)%smooth + mob(maxbin)%error) - (mob(nssobin)%smooth + mob(nssobin)%error)) ) ! extrapolate position
                        end if
                     end if
                     ! -------------------------------------------------------------------------------

                     ! get lower boundary (relative to maxbin)----------------------------------------
                     lowerbin = 0
                     do ibin = maxbin - 1 , 0, -1           !find lower position where locmob is significantly different from maximum
                        if (mob(maxbin)%smooth - mob(maxbin)%error .ge. mob(ibin)%smooth + mob(ibin)%error) then
                           lowerbin = maxbin - ibin
                           exit
                        end if
                     end do
                     !--------------------------------------------------------------------------------

                     !store results in param----------------------------------------------------------
                     param(ipt,part%i)%opt = Two*ssorad(maxbin,ipt)
                     param(ipt,part%i)%err = ssorad(max((upperbin + lowerbin),2),ipt)
                     !--------------------------------------------------------------------------------

                     !print tests---------------------------------------------------------------------
                     if(ltestsso) then
                        write(uout,'(a,I0,a,I0,a,I0)') 'mobility of pt ',ipt, " part ", part%i, "; current step: ", istep
                        write(uout,'(a,g15.5)') 'current dtran ', curdtranpt(ipt)
                        write(uout,'(a,I5)')  'number of bins: ', nssobin
                        write(uout,'(a, a, a, a, a, a)') "trans. rad", "mob", "msd", "error", "smoothed", "number of steps"
                        write(uout,'(i0, g15.5,a,g15.5,a,g15.5,a,g15.5,a,g15.5)') &
                        (ibin, ssorad(ibin,ipt), char(9), mob(ibin)%val, char(9), ssos(ibin, ipt)%d2, char(9), mob(ibin)%smooth, char(9), ssos(ibin, ipt)%n, ibin = 1,nssobin)
                        write(uout,'(a,I4)')  'average msd of pt ',ipt
                        write(uout,'(g15.5)') tots(ipt)%d2/real(tots(ipt)%n)
                        write(uout,'(a)')  ''
                        write(uout,'(a)')  ''
                     end if
                     !--------------------------------------------------------------------------------

                     !set next displacement parameter-------------------------------------------------
                     curdtranpt(ipt) = min(Two*ssorad((maxbin + max(upperbin,1)),ipt),maxdtransso(ipt))
                     invrsso(ipt) = real(nssobin) * Two / real(curdtranpt(ipt))
                     !--------------------------------------------------------------------------------
                  end if
               end do

               !prepare variables for next part----------------------------------------------------
               ssos = step(0, Zero, Zero)
               tots  = step(0, Zero, Zero)
               part%i = part%i + 1
               part%nextstep = part%nextstep + int(nstepzero*part%fac**(part%i - 1))
               if(part%i == part%n) then
                  part%nextstep = nstep     ! let last part end at end of simulation
               end if
               !-----------------------------------------------------------------------------------

            end if

         case (iAfterMacrostep)

               write(ucnf) curdtranpt, part%i, ssos, tots, param

         case (iAfterSimulation)

            call WriteHead(2, txroutine, uout)

            write(uout,'(a)') 'SSO - optimal displacement parameter'
            write(uout,'(a)') '------------------------------------'
            write(uout,'(a15,a15)') 'particle type' , 'optimal dtran'
            write(uout,'(a15,a15)') '---------------' , '---------------'
            write(uout,'(i15,g15.5)') (ipt, param(ipt,part%n)%opt, ipt = 1, npt)
            write(uout,'(a)') ''
            write(uout,'(a)') ''
            write(uout,'(a)')'Number of SSO parts'
            write(uout,'(i5)') part%n
            write(uout,'(a)') ''
            write(uout,'(a)') ''

            do ipt = 1, npt
               write(uout,'(a,I4)')  'displacement parameters of',ipt
               write(uout,'(a15,a,a15,a,a15,a,a20)')  'sso-part' , char(9), 'used dran' , char(9), 'optimal dtran' , char(9) , 'error on opt. dtran'
               write(uout,'(a15,a,a15,a,a15,a,a20)')  '---------------' , char(9), '---------------' , char(9), '---------------' , char(9) , '--------------------'
               write(uout,'(i15,a,g15.5,a,g15.5,a,g20.5)') &
               (ipart, char(9), param(ipt,ipart)%used, char(9),param(ipt,ipart)%opt, char(9), param(ipt,ipart)%err, ipart = 1 ,part%n)
               write(uout,'(a)') ''
               write(uout,'(a)') ''
            end do

         end select

      !........................................................................

         contains

         subroutine CalcLocMob(ipt)

            implicit none
            
            ! particle type to be calculated
            integer(4), intent(in)  :: ipt

            ! store the inverted number of steps
            real(8)              :: invntot
            ! required for Smooth subroutine:
            real(8)              :: btmp(0:nssobin)
            real(8)              :: ctmp(0:nssobin)
            real(8)              :: dtmp(0:nssobin)
            real(8)              :: x(0:nssobin)
            real(8)              :: y(0:nssobin)
            real(8)              :: dy(0:nssobin)
            real(8)              :: a(0:nssobin)

            integer(8)           :: ibin
            type(step)  :: stepbin ! steps done


            stepbin = step(0, Zero, Zero) 
            mob(0) = mobility(Zero, Zero, Zero)

            do ibin = 1, nssobin
               stepbin = stepbin + ssos(ibin,ipt)
               if(stepbin%n < 1) then
                  !assume 100% acceptance rate if no move was done
                  mob(ibin)%val = 0.15d0 * (real(ibin) / real(nssobin)  * curdtranpt(ipt))**2 ! 0.15 comes from 3/5 (from mean squared displacement in a sphere) and 1/4 (from squared diameter to squared radius
                  mob(ibin)%error = Zero
               else
                  !else: average displacement is displacement divided by number of steps; error calculated from variance
                  invntot = One/real(stepbin%n)
                  mob(ibin)%val = stepbin%d2 * invntot
                  mob(ibin)%error = sqrt((stepbin%d4 * invntot - (mob(ibin)%val)**2)*invntot)
               end if
            end do

            x = real( (/ (ibin, ibin = 0, nssobin) /) , KIND=8  )
            y = mob(0:nssobin)%val
            dy = mob(0:nssobin)%error
            !smooth using spline
            call Smooth(nssobin + 1, x, y, dy, real( (nssobin + 1) - sqrt(Two*(nssobin + 1)) , kind=8), a, btmp, ctmp, dtmp)

            mob(0:nssobin)%smooth = a


         end subroutine CalcLocMob

      !........................................................................

         real(8) function ssorad(bin,ipt)
               implicit none
               integer, intent(in)  :: bin
               integer, intent(in)  :: ipt
               ssorad = (Half * real(bin)/real(nssobin) * curdtranpt(ipt))
         end function ssorad

      !........................................................................

      end subroutine SSOSetup

      subroutine DoSSOUpdate(lacc, ipt, dr)
         implicit none
         logical, intent(in)  :: lacc      ! event of SSO-Move
         integer, intent(in)  :: ipt        ! number of moving particles
         real(8),    intent(in) :: dr(3)  ! suggested particle move

         integer(4)  :: ibin
         real(8)  :: d2
         d2=sum(dr(1:3)**2)
         ibin = ceiling(sqrt(d2)*invrsso(ipt))

         tots(ipt)%n = tots(ipt)%n + 1
         ssos(ibin, ipt)%n = ssos(ibin, ipt)%n + 1
         if ( lacc ) then
            ssos(ibin, ipt)%d2 = ssos(ibin, ipt)%d2 + d2
            tots(ipt)%d2 = tots(ipt)%d2 + d2
            ssos(ibin, ipt)%d4 = ssos(ibin, ipt)%d4 + d2**2
            tots(ipt)%d4 = tots(ipt)%d4 + d2**2
         end if

      end subroutine DoSSOUpdate

end module SSOModule

subroutine SSOUpdate(lacc, ipt, dr)
   use SSOModule, only: DoSSOUpdate
   implicit none
   logical, intent(in)  :: lacc      ! event of SSO-Move
   integer, intent(in)  :: ipt        ! number of moving particles
   real(8),    intent(in) :: dr(3)  ! suggested particle move

   call DoSSOUpdate(lacc, ipt, dr(1:3))

end subroutine

subroutine SSODriver(iStage)
   use SSOModule, only: SSOSetup
   implicit none
   integer(4), intent(in)  :: iStage

   call SSOSetup(iStage)

end subroutine
