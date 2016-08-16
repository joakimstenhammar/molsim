!************************************************************************
!*                                                                      *
!*     SSOModule                                                        *
!*                                                                      *
!************************************************************************

! ... Module for the SSO simulation

 module SSOModule
   use MCModule
!  integer(8),    allocatable    :: naccsso(:,:)  ! has to be of kind 8 to handle
!  integer(8),    allocatable    :: ntotsso(:,:)  ! long simulations
   integer(8),    allocatable    :: steptot(:)
   real(8), allocatable          :: d2tot(:)      ! long simulations
   real(8),       allocatable    :: dssostep(:,:) ! has to be of kind 8 to handle
   real(8),       allocatable    :: dssostep2(:,:)! has to be of kind 8 to handle
   integer(8),    allocatable    :: nssostep(:,:) ! long simulations
   real(8), allocatable          :: invrsso(:)
   real(8), allocatable          :: curdtranpt(:)
!  real(8), allocatable          :: dtranfac(:)
end module SSOModule

!************************************************************************
!*                                                                      *
!*     SSODriver                                                        *
!*                                                                      *
!************************************************************************

! ... Driver for the SSO simulation

subroutine SSODriver(iStage)

   use SSOModule

   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SSODriver'

   integer(4)  :: ipt
   integer     :: mbin, dbin, tmpstep, ibin, nstepend, lbin
   real(8)     :: InvInt, InvFlt
   logical, save  :: ltestsso

   real(8), allocatable, save :: locmob(:)
   real(8), allocatable, save :: locmobe(:)
   real(8), allocatable, save :: locmobs(:)
   real(8), allocatable, save :: dtranout(:,:,:)  ! output of dtrans
   real(8), allocatable, save :: dtransso(:) ! initial displacement parameters
   real(8), allocatable, save :: maxdtransso(:)! maximal allowed displacement parameters
   real(8), allocatable, save :: dtranfac(:)   !increase in displacement parameter
   integer(4), save           :: nstepzero     ! length of first part of simulation where local mobility is measured
   integer(4), save           :: nssobin
   real(8), save              :: partfac       ! increment in length of parts

   integer(8), save           :: istepnext
   integer(4), save           :: issopart
   integer(4), save           :: nssopart

   namelist /nmlSPartSSO/ dtransso, nstepzero, nssobin, nstepend, ltestsso, maxdtransso, dtranfac

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   select case (iStage)
!  case (iReadInput)   is omited, as lpspartsso is set during iStage == iWriteInput

   case (iWriteInput) ! Question: can the reading be moved to iReadInput?
!  Reply: not yet: firt the content of SSODriver has to be moved to IOMC, but let's wait

      if (.not. allocated(dtransso)) then 
         allocate(dtransso(npt))
         dtransso = 0.0E+00
      end if
      if (.not. allocated(maxdtransso)) then 
         allocate(maxdtransso(npt))
         maxdtransso = 0.0E+00
      end if
      if (.not. allocated(dtranfac)) then 
         allocate(dtranfac(npt))
         dtranfac = 0.0E+00
      end if

      dtransso  = One
      nstepzero   = ceiling(sqrt(real(nstep)))
      nstepend    = int(0.1*nstep)
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

      if (nstepzero + nstepend > nstep) call stop(txroutine, 'nstepzero + nstepend > nstep', uout)
      if(nstepend < nstepzero) call stop(txroutine, 'nstepend < nstepzero', uout)
      if(nstepzero == 0 .or. nstepend == 0) then
         partfac = 1
         nssopart = 1
      else if (nstepzero == nstepend) then
         partfac = One
         nssopart = int(nstep * InvInt(nstepend))
      else
         partfac = -real(nstepzero - nstep)*InvInt(nstep - nstepend)
         nssopart = int(log(real(nstepend*InvInt(nstepzero)))/log(partfac))
         if(nstepzero * (One - partfac**(nssopart + 1))/(One - partfac) < nstep) nssopart = nssopart + 1
      end if

      if(.not. allocated(curdtranpt)) then 
         allocate(curdtranpt(npt))
         curdtranpt = 0.0E+00
      end if
      if(.not. allocated(invrsso)) then 
         allocate(invrsso(npt))
         invrsso = 0.0E+00
      end if

      if(.not. allocated(d2tot)) then 
         allocate(d2tot(npt))
         d2tot = 0.0E+00
      end if
      if(.not. allocated(steptot)) then 
         allocate(steptot(npt))
         steptot = 0
      end if

      if(.not. allocated(locmob)) then 
         allocate(locmob(0:nssobin))
         locmob = 0.0E+00
      end if
      if(.not. allocated(locmobe)) then 
         allocate(locmobe(0:nssobin))
         locmobe = 0.0E+00
      end if
      if(.not. allocated(locmobs)) then 
         allocate(locmobs(0:nssobin))
         locmobs = 0.0E+00
      end if

      if(.not. allocated(nssostep)) then 
         allocate(nssostep(npt,nssobin))
         nssostep = 0
      end if
      if(.not. allocated(dssostep)) then 
         allocate(dssostep(npt,nssobin))
         dssostep = 0.0E+00
      end if
      if(.not. allocated(dssostep2)) then 
         allocate(dssostep2(npt,nssobin))
         dssostep2 = 0.0E+00
      end if

      if(.not.allocated(dtranout)) then 
         allocate(dtranout(npt,nssopart,3))
         dtranout = 0.0E+00
      end if
!     if(.not.allocated(dtranfac)) allocate(dtranfac(npt))

   case (iBeforeSimulation)

      if (txstart == 'continue') then
         read(ucnf) curdtranpt, dssostep, dssostep2, nssostep, istepnext, issopart, d2tot, steptot, dtranout
      else
         curdtranpt(1:npt) = dtransso(1:npt)
         dssostep = Zero
         dssostep2 = Zero
         nssostep = 0
!        naccsso = 0
!        ntotsso = 0
         istepnext = nstepzero
         issopart = 1
         dtranout = Zero
         d2tot = Zero
         steptot = 0
      end if

      if(nssopart == 1) istepnext = nstep

      do ipt = 1, npt
         invrsso(ipt) = InvFlt(curdtranpt(ipt)) * real(nssobin) * Two
      end do
!     dtranfac(1:npt) = max((nssobin + 3)*InvInt(nssobin),1.5)


   case (iBeforeMacrostep)                       ! ignore macrosteps

      continue

   case (iSimulationStep)

      if(istep == istepnext) then

         if(ltestsso) then
         end if

         do ipt = 1, npt

            if(.not. pspartsso(ipt) > Zero) cycle

            dtranout(ipt,issopart,1) = curdtranpt(ipt)

            call CalcLocMob(ipt)

            mbin = 0
            dbin = 0
            do ibin = 1, nssobin - 1
               if(locmobs(ibin) > locmobs(ibin + 1) .and. locmobs(ibin) > locmobs(mbin) ) mbin = ibin
            end do

            if(mbin == 0) then  ! no local maximum found, maximum at border
               mbin = nssobin
               dbin = ceiling((dtranfac(ipt) - One)*nssobin)
            else
               do ibin = mbin + 1 , nssobin !find position where locmob is significantly different from maximum
                  if (locmobs(mbin) - locmobe(mbin) .ge. locmobs(ibin) + locmobe(ibin)) then
                     dbin = (ibin - mbin) + 1
!                      dtranfac(ipt) = Half*(One + sqrt(Four*dtranfac(ipt) - Three)) !change dtranfac? - NO
                     exit
                  end if
               end do
               if (dbin == 0) then        ! no significantly different position found
                  dbin = ceiling(0.2*nssobin)
               end if
            end if

            do ibin = mbin - 1 , 0, -1           !find lower position where locmob is significantly different from maximum
               if (locmobs(mbin) - locmobe(mbin) .ge. locmobs(ibin) + locmobe(ibin)) then
                  lbin = mbin - ibin
                  exit
               end if
            end do

            dtranout(ipt,issopart,2) = Two*ssorad(mbin,ipt)
            dtranout(ipt,issopart,3) = ssorad(max((dbin + lbin),2),ipt)

            if(ltestsso) then
               write(ulist,'(a,a,I5,I5, I5)') '#','locmob of',ipt, issopart, istepnext
               write(ulist,'(a,I5)')  '#', nssobin
               write(ulist,'(g15.5,3(a9,g15.5))')
               write(ulist,'(g15.5,a,g15.5,a,g15.5,a,g15.5)') &
               (ssorad(ibin,ipt), char(9), locmob(ibin), char(9), locmobe(ibin), char(9), locmobs(ibin), ibin = 1,nssobin)
               write(ulist,'(a)')  ''
               write(ulist,'(a)')  ''
               write(ulist,'(a,a,I4)')  '#', 'dtranout of',ipt
               write(ulist,'(a,i5)')  '#', 1
               write(ulist,'(g15.5,a,g15.5,a,g15.5)') dtranout(ipt,issopart,1), char(9), dtranout(ipt,issopart,2), char(9), dtranout(ipt,issopart,3)
!              write(ulist,'(a)')  ''
!              write(ulist,'(a)')  ''
!              write(ulist,'(a,a,I4)') '#',  'naccrej of',ipt
!              write(ulist,'(a,i5)')  '#', nssobin
!              write(ulist,'(g15.5,a,I15,a,I15)') (ssorad(ibin,ipt), char(9), naccsso(ipt,ibin), char(9), ntotsso(ipt,ibin), ibin = 1,nssobin)
               write(ulist,'(a)')  ''
               write(ulist,'(a)')  ''
               write(ulist,'(a,a,I4)') '#',  'd2 of',ipt
               write(ulist,'(a,i5)')  '#', nssobin
               write(ulist,'(g15.5)') d2tot(ipt)*InvInt(steptot(ipt))
               write(ulist,'(a)')  ''
               write(ulist,'(a)')  ''
            end if

            curdtranpt(ipt) = min(Two*ssorad((mbin + max(dbin,1)),ipt),maxdtransso(ipt))

         end do

         do ipt = 1, npt
            invrsso(ipt) = InvFlt(curdtranpt(ipt)) * real(nssobin) * Two
         end do

!        naccsso = 0
!        ntotsso = 0

         nssostep = 0
         dssostep = Zero
         dssostep2 = Zero
         d2tot = Zero
         steptot = 0

         issopart = issopart + 1

         istepnext = istepnext + int(nstepzero*partfac**(issopart - 1))

         if(issopart == nssopart) istepnext = nstep     ! let last part end at end of simulation

      end if

   case (iAfterMacrostep)

         write(ucnf) curdtranpt, dssostep, dssostep2, nssostep, istepnext, issopart, d2tot, steptot, dtranout

   case (iAfterSimulation)

      call WriteHead(2, txroutine, uout)

      write(uout,'(a)') 'SSO - optimal displacement parameter'
      write(uout,'(a)') '------------------------------------'
      write(uout,'(a15,a15)') 'particle type' , 'optimal dtran'
      write(uout,'(a15,a15)') '---------------' , '---------------'
      write(uout,'(i15,g15.5)') (ipt, dtranout(ipt,nssopart,2), ipt = 1, npt)
      write(uout,'(a)') ''
      write(uout,'(a)') ''
      write(uout,'(a)')'Number of SSO parts'
      write(uout,'(i5)') nssopart
      write(uout,'(a)') ''
      write(uout,'(a)') ''

      do ipt = 1, npt
         write(uout,'(a,I4)')  'displacement parameters of',ipt
         write(uout,'(a15,a,a15,a,a15,a,a20)')  'sso-part' , char(9), 'used dran' , char(9), 'optimal dtran' , char(9) , 'error on opt. dtran'
         write(uout,'(a15,a,a15,a,a15,a,a20)')  '---------------' , char(9), '---------------' , char(9), '---------------' , char(9) , '--------------------'
         write(uout,'(i15,a,g15.5,a,g15.5,a,g20.5)') &
         (issopart, char(9), dtranout(ipt,issopart,1), char(9),dtranout(ipt,issopart,2), char(9), dtranout(ipt,issopart,3), issopart = 1,nssopart)
         write(uout,'(a)') ''
         write(uout,'(a)') ''
      end do

!       if(allocated(curdtranpt)) deallocate(curdtranpt)
!       if(allocated(invrsso)) deallocate(invrsso)

!       if(allocated(d2tot)) deallocate(d2tot)
!       if(allocated(steptot)) deallocate(steptot)

!       if(allocated(locmob)) deallocate(locmob)
!       if(allocated(locmobe)) deallocate(locmobe)
!       if(allocated(locmobs)) deallocate(locmobs)

!       if(allocated(nssostep)) deallocate(nssostep)
!       if(allocated(dssostep)) deallocate(dssostep)
!       if(allocated(dssostep2)) deallocate(dssostep2)

!       if(allocated(dtranout)) deallocate(dtranout)

   end select

!........................................................................

   contains

   subroutine CalcLocMob(issopt)

      use SSOModule
      implicit none

      integer(4), intent(in)  :: issopt

      real(8)              :: fac
      real(16)             :: tmpsum
      real(16)             :: tmpsume
      real(8)              :: InvInt, InvFlt
      real(8)              :: invntot
      real(8)              :: btmp(0:nssobin)
      real(8)              :: ctmp(0:nssobin)
      real(8)              :: dtmp(0:nssobin)

      integer(8)           :: ibin, ntot

      tmpsum = Zero
      tmpsume = Zero
      ntot = 0
      locmob(0) = Zero
      locmobe(0) = Zero

      do ibin = 1, nssobin
         ntot = ntot + nssostep(issopt,ibin)
         tmpsum = tmpsum + dssostep(issopt,ibin)
         tmpsume = tmpsume + dssostep2(issopt,ibin)
         if(ntot < 1) then
            locmob(ibin) = Three/Five * (ibin * InvInt(nssobin) * Half * curdtranpt(ipt))**2
            locmobe(ibin) = Zero
         else
            invntot = InvInt(ntot)
            locmob(ibin) = tmpsum * invntot
            locmobe(ibin) = sqrt((tmpsume * invntot - (locmob(ibin))**2)*invntot)
         end if
      end do
      call Smooth(nssobin + 1,real( (/ (ibin, ibin = 0, nssobin) /) , KIND=8  ) , locmob , locmobe , real( (nssobin + 1) - sqrt(Two*(nssobin + 1)) , kind=8), locmobs, btmp, ctmp, dtmp)

   end subroutine CalcLocMob

!........................................................................

real(8) function ssorad(bin,ipt)
      use SSOModule
      implicit none
      integer, intent(in)  :: bin
      integer, intent(in)  :: ipt
      real(8)  InfInt
      ssorad = (bin * InvInt(nssobin) * Half * curdtranpt(ipt))
end function ssorad

!........................................................................

end subroutine SSODriver

!************************************************************************
!*                                                                      *
!*     SSOMove                                                          *
!*                                                                      *
!************************************************************************

! ... perform one single-particle SSO trial move

subroutine SSOMove(iStage)

   use SSOModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='SSOMove'
   logical    :: lboxoverlap, lhsoverlap, lhepoverlap
   integer(4) :: iploc, dnpcl
   real(8)    :: weight, MCWeight, UmbrellaWeight, MCPmfWeight

   integer(4), save :: jptsph = 1  ! type of particle at which the grafted chains are attached
   integer(4) :: jpsph, ipsurf
   integer(4) :: ihost
   integer(4) :: ibin, ipt
   real(8) :: dx, dy, dz, norm, d2, dtr
   real(8) :: Random

   if (ltrace) call WriteTrace(3, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   imovetype = ispartsso

#if defined (_PAR_)
! ... the following is not (yet) allowed since iseed will not be syncronized
   ipnptm(1:np) = 0                                 ! necessary for later allreduce
#endif

! .............. define particle(s) to be moved ..............

! ... consider particle ipmove

   nptm          = 1
   ipnptm(1) = ipmove
   lptm(ipmove)  =.true.
   ipt = iptpn(ipmove)


! .............. calculate a trial configuration ...............

! ... get translational displacement

   dtr = abs(curdtranpt(ipt))

   do
      dx = Random(iseed)-Half
      dy = Random(iseed)-Half
      dz = Random(iseed)-Half
      d2 = dx**2 + dy**2 + dz**2
      if ( d2 < Fourth) exit
   end do

   dx = dx*dtr
   dy = dy*dtr
   dz = dz*dtr
   d2 = d2 * dtr**2

! ... get trial coordinates and store translational displacemnt

   rotm(1:3,1) = (/ ro(1,ipmove)+dx , ro(2,ipmove)+dy , ro(3,ipmove)+dz /)
   drotm(1:3,1) = (/ dx , dy, dz /)

!    if (lweakcharge) call TrialCharge(nptm, ipnptm, iptpn, .false., latweakcharge, laz, laztm)
!    if (lfixzcoord(iptmove)) call FixedZCoord
!    if (lfixxycoord(iptmove)) call FixedXYCoord
!    if (lfixchainstartspart) call FixedChainStart
!
!    if(lfixzcoord(iptmove) .or. lfixxycoord(iptmove) .or. lfixchainstartspart) d2 = sum(drotm(1:3,1)**2)

!-------------------------------------------------------------------------------------

   call CheckPartBCTM(nptm, rotm, lboxoverlap)

!    if (lpolyatom .or. lellipsoid .or. lsuperball) &
!       call GetRandomTrialOri(drot(iptmove), iseed, ori(1,1,ipmove), oritm(1,1,iploc))
!
!    if (lfixedori) then           ! lfixedori atains its value in coordinate.F90
!        call AddNeighbours
!        call UpdateOri
!    end if

   call SetTrialAtomProp
!    if (lradatbox) call CheckAtomBCTM(natm, rtm, lboxoverlap)
   if (itestmc == 2) call TestMCMove(uout)

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

   if (lboxoverlap) goto 200

! ............. evaluate energy difference ...............

   call DUTotal(lhsoverlap, lhepoverlap)
   if (lhsoverlap .or. lhepoverlap) goto 200

! ............. calculate nonenergetic weights .............

   weight = One
   dnpcl = Zero
   if (lcl1spart(iptmove)) then
      if (lvlist) call ClusterMember('new', .false., .false., radcl1, pselectcl1)       ! calculate npclnew
      if (lllist) call ClusterMemberLList('new', .false., .false., radcl1, pselectcl1)  ! calculate npclnew
      dnpcl = npclnew-npclold
      if (dnpcl /= 0) weight = weight*(One-pselectcl1(iptmove))**dnpcl
   end if

   if (lmcweight) weight = weight*MCWeight()
   if (lautumb) weight = weight*UmbrellaWeight(1)
   if (lmcpmf) weight = weight*MCPmfWeight(1)

! ............. decide new configuration .............

200 continue
   call Metropolis(lboxoverlap, lhsoverlap, lhepoverlap, weight, du%tot*beta)

! .............. update .............

   ibin = ceiling(sqrt(d2)*invrsso(ipt))
   steptot(ipt) = steptot(ipt) + 1
   nssostep(ipt,ibin) = nssostep(ipt,ibin) + 1
   if (ievent == imcaccept) then
      dssostep(ipt,ibin) = dssostep(ipt,ibin) + d2
      dssostep2(ipt,ibin) = dssostep2(ipt,ibin) + d2**2
      d2tot(ipt) = d2tot(ipt) + d2
      call MCUpdate       ! update energies and coordinates
   end if

!    if (lautumb) call UmbrellaUpdate              ! update weight function for umbrella potential
!    if (lmcpmf) call MCPmfUpdate                  ! update weight function for mc pmf
!    if (lshiftzcom(iptmove)) call ShiftZDirection

end subroutine SSOMove

