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

! relation among energy routines

!     UTotal
!       !
!       !   lcharge
!       !--------------
!       !             !
!       !             !---------- UTwoBody(A/ALList/ACellList/P)
!       !             !   lewald
!       !             !---------- UEwald
!       !             !   lrf
!       !             !---------- UIntraReac
!       !
!       !   lweakcharge
!       !--------------
!       !             !
!       !             !--------- UWeakCharge(A/P)
!       !             !  lewald
!       !             !--------- UEwald
!       !
!       !   ldipole
!       !--------------
!       !             !
!       !             !---------- UTwoBodyP
!       !             !
!       !             !---------- UDipole
!       !                            !
!       !                            !---------- UDipoleP     (field, efg, vir from q, sdm)
!       !                            !   lewald
!       !                            !---------- UDipoleEwald (field, efg, vir from q, sdm)
!       !   lpolarization
!       !--------------
!       !             !
!       !             !---------- UTwoBodyP
!       !             !
!       !             !---------- UManyBodyP
!       !                           !
!       !                           !-------- FieldStat      (pot, field from q, sdm)
!       !                           !-------- FieldStatEwald (pot, field from q, sdm)
!       !                           !
!       !                           !-------- IterIdm
!       !                           !           !
!       !                           !           !-------- FieldIdm      (field from idm)
!       !                           !           !  lewald
!       !                           !           !-------- FieldIdmEwald (field from idm)
!       !                           !
!       !                           !-------- FieldTot      (field, efg, vir from q, sdm, idm)
!       !                           ! lewald
!       !                           !-------- FieldTotEwald (field, efg, vir from q, sdm, idm)
!       !   ldipolesph
!       !--------------
!       !             !
!       !             !---------- UDipoleSph
!       !   ldieldis
!       !--------------
!       !             !
!       !             !---------- UDielDis
!       !   lchain
!       !---------- UBond
!       !
!       !   lchain
!       !---------- UAngle
!       !
!       !   lclink
!       !---------- UCrossLink
!       !
!       !   luext
!       !---------- UExternal

!************************************************************************
!> \page energy energy.F90
!! **EnergyModule**
!! *module for energy*
!************************************************************************


module EnergyModule

# ifdef F03_CBIND
   use, intrinsic :: iso_c_binding           ! for FFTW
# endif
   use MolModule
   use CardinalBSplineModule

! ... Ewald: setup and constants

# ifdef F03_CBIND
   include 'fftw3.f03'
# endif
   real(8)       :: ualpha2                  ! ualpha**2
   real(8)       :: ewaldfac1                ! One/(sqrt(Pi)*ualpha) * (Two*ualpha2)
   real(8)       :: ewaldfac2                ! One/(sqrt(Pi)*ualpha) * (Two*ualpha2)**2 / Three
   real(8)       :: ewaldfac3                ! One/(sqrt(Pi)*ualpha) * (Two*ualpha2)**3 / (Three*Five)
   real(8)       :: ewaldselffac1            ! Two*ualpha/sqrt(Pi)
   real(8)       :: ewaldselffac2            ! Four*ualpha**3/(Three*sqrt(Pi))
   integer(4)    :: nkvec                    ! number of k-vectors in one quadrant
   integer(4)    :: nkvec_word               ! number of k-vector words
   integer(4)    :: nkvec2d                  ! number of k-vectors in one quadrant (2d correction)
   real(8)       :: time_ewald               ! for profiling of ewald summation

! ... real space

   real(8), allocatable :: potstat(:)        ! electrostatic potential, from charges and static dipoles
   real(8), allocatable :: estat(:,:)        ! electrostatic field, from charges and static dipoles
   real(8), allocatable :: eidm(:,:)         ! electrostatic field, from induced dipoles
   real(8), allocatable :: etot(:,:)         ! electrostatic field, total
   real(8), allocatable :: efg(:,:)          ! electrostatic field gradient (xx,yy,zz,xy,xz,yz), total
   real(8)              :: virelec           ! virial contribution from electrostatic interactions

   real(8), allocatable :: utwobnew(:)       ! new two-body energy (temporary use)
   real(8), allocatable :: utwobold(:)       ! old two-body energy (temoprary use)

   real(8), allocatable :: uaux(:)           ! temporary use

! ... Ewald: reciprocal space: standard

   integer(4) :: naewald
   real(8),    allocatable :: kfac(:)        ! factor for k-summation
   complex(8), allocatable :: eikx(:,:)      ! exp(ikx)
   complex(8), allocatable :: eiky(:,:)      ! exp(iky)
   complex(8), allocatable :: eikz(:,:)      ! exp(ikz)
   complex(8), allocatable :: eikyzm(:)      ! exp(-iky)*exp(ikz), temporary use
   complex(8), allocatable :: eikyzp(:)      ! exp(+iky)*exp(ikz), temporary use
   complex(8), allocatable :: eikr(:,:)      ! temporary use
   complex(8), allocatable :: sumeikr(:,:)   ! sum(q*exp(sx*ikx)*exp(sy*iky)*exp(ikz)) (sx,sy) = (-1,-1),(-1,1),(1,-1),(1,1)
   complex(8), allocatable :: sumeikrd(:,:)  ! as sumeikr but also includes dipoles

   complex(8), allocatable :: eikxtm(:,:)    ! exp(ikx), trial configuration
   complex(8), allocatable :: eikytm(:,:)    ! exp(iky), trial configuration
   complex(8), allocatable :: eikztm(:,:)    ! exp(ikz), trial configuration
   complex(8), allocatable :: eikyzmtm(:)    ! exp(-iky)*exp(ikz), trial configuration, temporary use
   complex(8), allocatable :: eikyzptm(:)    ! exp(+iky)*exp(ikz), trial configuration, temporary use
   complex(8), allocatable :: eikrtm(:,:)    ! temporary use
   complex(8), allocatable :: sumeikrtm(:,:) ! sum(q*exp(sx*ikx)*exp(sx*iky)*exp(ikz), trial configuration (sx,sy) = (-1,-1),(-1,1),(1,-1),(1,1)

   complex(8), allocatable :: eikraux(:)     ! temporary use

! ... Ewald: reciprocal space: standard: 2d layer correction

   real(8), allocatable :: kfac2d(:)         ! factor for k-summation (2d correction)
   real(8), allocatable :: sinkx(:,:)        ! sin(kx)
   real(8), allocatable :: coskx(:,:)        ! cos(kx)
   real(8), allocatable :: sinky(:,:)        ! sin(ky)
   real(8), allocatable :: cosky(:,:)        ! cos(ky)
   real(8), allocatable :: sumtrig(:,:)      ! sum(q*sin(kx)*sin(ky)*sinh(kp*z)) (xyz = sss, css, scs, ccs, ssc, csc, scc, ccc)
   real(8), allocatable :: termsss(:)        ! sum(q*sin(kx)*sin(ky)*sinh(kp*z))
   real(8), allocatable :: termcss(:)        ! sum(q*cos(kx)*sin(ky)*sinh(kp*z))
   real(8), allocatable :: termscs(:)        ! sum(q*sin(kx)*cos(ky)*sinh(kp*z))
   real(8), allocatable :: termccs(:)        ! sum(q*cos(kx)*cos(ky)*sinh(kp*z))
   real(8), allocatable :: termssc(:)        ! sum(q*sin(kx)*sin(ky)*cosh(kp*z))
   real(8), allocatable :: termcsc(:)        ! sum(q*cos(kx)*sin(ky)*cosh(kp*z))
   real(8), allocatable :: termscc(:)        ! sum(q*sin(kx)*cos(ky)*cosh(kp*z))
   real(8), allocatable :: termccc(:)        ! sum(q*cos(kx)*cos(ky)*cosh(kp*z))

   real(8), allocatable :: sinkxtm(:,:)      ! sin(kx), trial configuration
   real(8), allocatable :: coskxtm(:,:)      ! cos(kx), trial configuration
   real(8), allocatable :: sinkytm(:,:)      ! sin(ky), trial configuration
   real(8), allocatable :: coskytm(:,:)      ! cos(ky), trial configuration
   real(8), allocatable :: sumtrigtm(:,:)    ! sum(q*sin(kx)*sin(ky)*sinh(kp*z)), trial configuration (xyz = sss, css, scs, ccs, ssc, csc, scc, ccc)

# ifdef F03_CBIND
! ... Ewald: reciprocal space: smooth particle mesh

   integer(4) :: meshsize(3)                 ! length of the mesh in each dimension
   real(8) :: dkdr(3)                        ! mesh units per real units
   real(8) :: d2kdr(6)                       ! dkdr squared

! ... for FFTW3 we need C compatible data types   JH July 2012

   real(C_DOUBLE), pointer :: Qmesh(:,:,:) => null()              ! generalized charge distribution
   real(C_DOUBLE), pointer :: QmeshTM(:,:,:) => null()            ! generalized charge distribution (trial move)
   complex(C_DOUBLE_COMPLEX), pointer :: FQmesh(:,:,:) => null()  ! transformed charge distribution
   type(C_PTR) :: qmesh_ptr, qmeshtm_ptr, fqmesh_ptr              ! c-pointers for the allocation of the above
   type(C_PTR) :: plan_fwd, plan_bwd                              ! pointers to store the plan for forward and backward
   logical, save :: plan_fwd_done=.false., plan_bwd_done=.false.
# ifdef _TEST_
   integer(C_INT), parameter :: plan_methode = FFTW_ESTIMATE      ! select a planning strategy
# else
   integer(C_INT), parameter :: plan_methode = FFTW_MEASURE       ! select a planning strategy
# endif
!  integer(C_INT), parameter :: plan_methode = FFTW_PATIENT       ! select a planning strategy

   real(8),   allocatable :: meshfac(:,:,:)  ! influence function in reciprocal space
   real(8),   allocatable :: energyfac(:,:,:)! reduced meshfac to compensate for doubly counted elements
   real(8),   allocatable :: virfac(:,:,:)   ! virial contribution in reciprocal space
   integer(4),allocatable :: meshmax(:,:)    ! maximal mesh position of spline
   integer(4),allocatable :: meshmaxtm(:,:)  !        "                                      (trial move)
   real(8),   allocatable :: spline(:,:,:)   ! spline value of exp(ikr) for a given r
   real(8),   allocatable :: splinetm(:,:,:) !        "                                      (trial move)
   real(8),   allocatable :: deriv(:,:,:)    ! spline value of d[exp(ikr)]/dr for a given r
   real(8),   allocatable :: derivtm(:,:,:)  !        "                                      (trial move)
   real(8),   allocatable :: deriv2(:,:,:)   ! spline value of d^2[exp(ikr)]/dr^2 for a given r
   real(8),   allocatable :: deriv2tm(:,:,:) !        "                                      (trial move)
   integer(4)             :: MeshMemSize     ! size of mesh in memory
   real(8),   allocatable :: meshaux(:)      ! meshaux for packreduce
# endif

   real(8) :: urecold                        ! reciprocal energies before mc trial move
   real(8) :: urecnew                        ! reciprocal energies after mc trial move

! ... for UExternalSphDielBoundary_q, UExternalSphDielBoundary_p, UExternalSphDielBoundary

   real(8), save :: epsi1FourPi              ! Epsi0FourPi/epsi1
   real(8), save :: epsi2FourPi              ! Epsi0FourPi/epsi2
   real(8), save :: eta                      ! epsi1 / epsi2
   real(8), save :: delta                    ! (1-eta)/(1+eta)
   real(8), save :: elunit

   integer(4), parameter :: mlmax = 50
   complex(8) QQ(0:mlmax, 0:mlmax)           ! electrostatic multipole moment
   complex(8) QQtm(0:mlmax, 0:mlmax)         ! electrostatic multipole moment, temporary use
   real(8) :: lfac(0:mlmax)                  ! factor containing relative permittivities
   real(8) :: mfac(0:mlmax)                  ! 1 for m = 0 else 2
   real(8) :: lfac1                          ! factor containing relative permittivities
   real(8) :: lfac2                          ! factor containing relative permittivities
   real(8) :: lfac3(mlmax)                   ! factor containing relative permittivities
   integer(4) :: jalowUExternal, jauppUExternal

! ... for interpolation of external potential

   integer(4), parameter :: nhext = 100
   real(8), save :: hext(0:nhext), uext(0:nhext), dhexti

!   interface
!      subroutine CardinalBSpline(order, pos, kmax, values, derivs, derivs2)
!         integer(4), intent(in)  :: order                     ! order of the cardinal B-spline
!         real(8),    intent(in)  :: pos                       ! position
!         integer(4), intent(out) :: kmax                      ! floor(pos)
!         real(8),  intent(inout) :: values(0:order-1)         ! values at M_n(pos-kmin)...M_n(pos-(kmin+n-1))
!         real(8), intent(out), optional :: derivs(0:order-1)  ! values of first derivative
!         real(8), intent(out), optional :: derivs2(0:order-1) ! values of second derivative
!      end subroutine CardinalBSpline
!   end interface

end module EnergyModule

!************************************************************************
!> \page energy energy.F90
!! **UTotal**
!! *calculate energies, forces, torques, virial, and pressure*
!************************************************************************


!     calculated properties:   u%tot
!                              u%twob()
!                              u%oneb()
!                              u%rec
!                              u%stat
!                              u%pol
!                              u%bond
!                              u%angle
!                              u%crosslink
!                              u%external
!                              force(1:3,1:na),  forceo(1:3,1:np)
!                              torque(1:3,1:na), torqueo(1:3,1:np)
!                              virial,            prsr
!

subroutine UTotal(iStage)

   use EnergyModule
   implicit none

   integer(4), intent(in) :: iStage

   character(40), parameter :: txroutine ='UTotal'

   integer(4) :: ip, ia
   real(8) :: fac, dx, dy, dz, virmolecule

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(.not.allocated(uaux)) then
      allocate(uaux(2*((nptpt+1+npt+1)+8)))  ! 2* for scratch use
      uaux = 0.0E+00
   end if
   if(.not.allocated(u%twob)) then
      allocate(u%twob(0:nptpt), du%twob(0:nptpt), u%oneb(0:npt), du%oneb(0:npt))
   end if
   if(ldipole .or. ldipolesph) then
      if(.not.allocated(potstat)) then
         allocate(potstat(na_alloc), estat(3,na_alloc), efg(6,na_alloc))
         potstat = 0.0E+00
         estat = 0.0E+00
         efg = 0.0E+00
      end if
   end if
   if(lpolarization) then
      if(.not.allocated(potstat)) then
         allocate(potstat(na_alloc), estat(3,na_alloc), eidm(3,na_alloc), etot(3,na_alloc), efg(6,na_alloc))
         potstat = 0.0E+00
         estat = 0.0E+00
         eidm = 0.0E+00
         etot = 0.0E+00
         efg = 0.0E+00
      end if
   end if

! ... initiate

   u%tot            = Zero
   u%twob           = Zero
   force (1:3,1:na) = Zero
   torque(1:3,1:na) = Zero
   diptot(1:3,1:na) = Zero
   virial           = Zero
   u%oneb(0:npt)    = Zero
   u%bond           = Zero
   u%angle          = Zero
   u%crosslink      = Zero
   u%external       = Zero

! .............. select appropiate energy routines ............

   if (lcharge) then   ! atoms possessing charges

      if (lmonoatom) then
         if (lvlist) call UTwoBodyA
         if (lllist) call UTwoBodyALList
         if (lCellList) call UTwoBodyACellList
      else
         call UTwoBodyP
      end if
      if (lewald) call UEwald
      if (lrf   ) call UIntraReac

   else if (lweakcharge) then  ! atoms possesing weak charges (titrating system)

      if (lmonoatom) then
        call UWeakChargeA
      else
        call UWeakChargeP
      end if

      if (lewald) call UEwald

if (itest == 90) then
      call writehead(3,txroutine, uout)                                     !cc
      write(uout,'(a,100l4)') ' laz',laz(1:na)                              !cc
      write(uout,'(a,10f10.5)') 'u%twob(0:nptpt)', u%twob(0:nptpt)          !cc
end if

   else if (ldipole) then  ! atoms possening charges and dipoles

      call UTwoBodyP
      call UDipole

   else if (lpolarization) then  ! atoms posessing charges, dipoles, and polarizablities

      call UTwoBodyP
      call UManyBodyP

   else if (ldipolesph) then  ! atoms possening charges and dipoles, spherical boundary condition, and image charge

      call UDipoleSph

   else if (ldieldis) then  ! atoms possesing charge in a system with dielectric discontinuities

      call UDielDis

   end if

   if (lchain) call UBond
   if (lchain) call UAngle
   if (lclink) call UCrossLink

   if (luext) call UExternal

#if defined (_PAR_)
   if (ltrace) call WriteTrace(4, trim(txroutine)//'_before Pack', iStage)
   if (ltime) call CpuAdd('start', 'comm', 2, uout)
   call PackReduceU(nptpt+1, npt+1, u%tot, u%twob, u%oneb, u%rec, u%stat, u%pol, u%bond, u%angle, u%crosslink, u%external, uaux)
   call par_allreduce_reals(force , vaux, 3*na)
   if (ldipole .or. lpolarization .or. ldipolesph) call par_allreduce_reals(torque, vaux, 3*na)
   call par_allreduce_real(virial, raux)
   if (ltime) call CpuAdd('stop', 'comm', 4, uout)
   if (ltrace) call WriteTrace(4, trim(txroutine)//'_after Pack', iStage)
#endif

! ............... calculate forces, torques, virial, and pressure ...............

   virmolecule = Zero
   forceo(1:3,1:np) = Zero
   torqueo(1:3,1:np) = Zero
   fac = (sclene/scllen)/sclfor
   do ia = 1, na
      ip = ipnan(ia)
      dx = r(1,ia)-ro(1,ip)
      dy = r(2,ia)-ro(2,ip)
      dz = r(3,ia)-ro(3,ip)
      virmolecule = virmolecule + dx*force(1,ia) + dy*force(2,ia) + dz*force(3,ia)
      force(1,ia) = force(1,ia) * fac
      force(2,ia) = force(2,ia) * fac
      force(3,ia) = force(3,ia) * fac
      forceo(1,ip)  = forceo(1,ip)  + force(1,ia)
      forceo(2,ip)  = forceo(2,ip)  + force(2,ia)
      forceo(3,ip)  = forceo(3,ip)  + force(3,ia)
      torqueo(1,ip) = torqueo(1,ip) + dy*force(3,ia)-dz*force(2,ia) + torque(1,ia)*fac
      torqueo(2,ip) = torqueo(2,ip) + dz*force(1,ia)-dx*force(3,ia) + torque(2,ia)*fac
      torqueo(3,ip) = torqueo(3,ip) + dx*force(2,ia)-dy*force(1,ia) + torque(3,ia)*fac
   end do

   virial = virial + virmolecule

! ... if (lmd) temp from previous time step is used

   prsr = (np*Boltz*temp*scltem-virial*sclene/(Three*AvNo))/(vol*sclvol)/sclpre

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine UTotal

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UTwoBodyA**
!! *calculate two-body potential energy; only monoatomic particles*
!************************************************************************


subroutine UTwoBodyA

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UTwoBodyA'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf, Getnpmyid
   real(8)    :: dx, dy, dz, r2, d, usum, fsum, virtwob

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (lmc) then
            if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 > rcut2) cycle
         if (r2 < r2atat(iptjpt)) then
            usum = 1d10                ! emulate hs overlap
         else if (r2 < r2umin(iptjpt)) then
            call StopUTwoBodyA
         else
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
               if (ibuf > nbuf) call StopIbuf('txptpt',iptjpt)
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                   d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
            fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                   d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))
         end if

         u%twob(iptjpt) = u%twob(iptjpt) + usum
         force(1,ip) = force(1,ip) + (fsum * dx)
         force(2,ip) = force(2,ip) + (fsum * dy)
         force(3,ip) = force(3,ip) + (fsum * dz)
         force(1,jp) = force(1,jp) - (fsum * dx)
         force(2,jp) = force(2,jp) - (fsum * dy)
         force(3,jp) = force(3,jp) - (fsum * dz)
         virtwob     = virtwob     - (fsum * r2)

      end do

   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot     + u%twob(0)
   virial    = virial    + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUTwoBodyA
   character(40), parameter :: txroutine ='StopUTwoBodyA'
   write(uout,'(a,i5)')  'ip', ip
   write(uout,'(a,i5)')  'jp', jp
   write(uout,'(a,3e15.5)') 'ro(ip)        = ', ro(1:3,ip)
   write(uout,'(a,3e15.5)') 'ro(jp)        = ', ro(1:3,jp)
   write(uout,'(a,3e15.5)') 'boxlen2       = ', boxlen2
   write(uout,'(a,3e15.5)') 'dpbc          = ', dpbc
   write(uout,'(a,i15)')    'uout           = ', uout
   write(uout,'(a,i15)')    'myid           = ', myid
   write(uout,'(a,i15)')    'iptjpt         = ', iptjpt
   write(uout,'(a,e15.5)')  'r2             = ', r2
   write(uout,'(a,e15.5)')  'r2umin(iptjpt) = ', r2umin(iptjpt)
   call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)
end subroutine StopUTwoBodyA

!........................................................................

end subroutine UTwoBodyA

!************************************************************************
!> \page energy energy.F90
!! **UTwoBodyALList**
!! *calculate two-body potential energy; only monoatomic particles*
!************************************************************************

!     linked list version

subroutine UTwoBodyALList

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UTwoBodyALList'
   integer(4) :: ip, iploc, ipt, jp, jpt, iptjpt, ibuf, Getnpmyid, icell, Getncellllist
   real(8)    :: dx, dy, dz, r2, d, usum, fsum, virtwob

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero
   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)

      call SetNCell(ro(1:3,ip), rcut, lcellllist)
      do icell = 1, Getncellllist()
         if (.not.lcellllist(icell)) cycle
         jp = headllist(icell)
         do
            if (jp == 0) exit
            if (jp <= ip) goto 199
            jpt = iptpn(jp)
            iptjpt = iptpt(ipt,jpt)
            dx = ro(1,ip)-ro(1,jp)
            dy = ro(2,ip)-ro(2,jp)
            dz = ro(3,ip)-ro(3,jp)
            call PBCr2(dx,dy,dz,r2)
            if (r2 > rcut2) goto 199
            if (r2 < r2atat(iptjpt)) then
               usum = 1d10                ! emulate hs overlap
            else if (r2 < r2umin(iptjpt)) then
               call StopUTwoBodyALList
            else
               ibuf = iubuflow(iptjpt)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
                  if (ibuf > nbuf) call StopIbuf('txptpt',iptjpt)
               end do
               d = r2-ubuf(ibuf)
               usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                      d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
               fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                      d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))
            end if

            u%twob(iptjpt) = u%twob(iptjpt) + usum
            force(1,ip) = force(1,ip) + (fsum * dx)
            force(2,ip) = force(2,ip) + (fsum * dy)
            force(3,ip) = force(3,ip) + (fsum * dz)
            force(1,jp) = force(1,jp) - (fsum * dx)
            force(2,jp) = force(2,jp) - (fsum * dy)
            force(3,jp) = force(3,jp) - (fsum * dz)
            virtwob     = virtwob     - (fsum * r2)

  199       continue
            jp = jpllist(jp)
         end do
      end do

   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot         + u%twob(0)
   virial    = virial        + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUTwoBodyALList
   character(40), parameter :: txroutine ='StopUTwoBodyALList'
   write(uout,'(a,i5)')  'ip', ip
   write(uout,'(a,i5)')  'jp', jp
   write(uout,'(a,3e15.5)') 'ro(ip)         = ', ro(1:3,ip)
   write(uout,'(a,3e15.5)') 'ro(jp)         = ', ro(1:3,jp)
   write(uout,'(a,3e15.5)') 'boxlen2        = ', boxlen2
   write(uout,'(a,3e15.5)') 'dpbc           = ', dpbc
   write(uout,'(a,i15)')    'uout           = ', uout
   write(uout,'(a,i15)')    'myid           = ', myid
   write(uout,'(a,i15)')    'iptjpt         = ', iptjpt
   write(uout,'(a,e15.5)')  'r2             = ', r2
   write(uout,'(a,e15.5)')  'r2umin(iptjpt) = ', r2umin(iptjpt)
   call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)
end subroutine StopUTwoBodyALList

!........................................................................

end subroutine UTwoBodyALList

!************************************************************************
!> \page energy energy.F90
!! **UTwoBodyACellList**
!! *calculate two-body potential energy*
!************************************************************************

!     only monoatomic particles
!     cell list version

subroutine UTwoBodyACellList

   use MolModule,      only: ro, iptpn, iptpt, rcut2, Zero, r2atat, r2umin
   use MolModule,      only: nptpt, np
   use MolModule,      only: nproc, myid, uout, ltime
   use MolModule,      only: iubuflow, ubuf, lmonoatom, nbuf
   use MolModule,      only: u, force, virial
   use CellListModule, only: pcellro, cell_type, cellip, ipnext

   implicit none
   character(40), parameter :: txroutine ='UTwoBodyACellList'

   integer(4) :: ip, jp, ipt, jpt, iptjpt, jploc, ibuf
   real(8)    :: dr(3), r2, d, usum, fsum, virtwob

   type(cell_type), pointer :: icell, ncell
   integer(4)               :: incell

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero

   do ip = 1, np
      ipt = iptpn(ip)
      icell => cellip(ip)%p
      do incell = 1 + myid, icell%nneighcell, nproc ! increment with nproc to have parallel execution
         ncell => icell%neighcell(incell)%p
         jp = ncell%iphead
         do jploc = 1, ncell%npart
            if (jp > ip) then
               jpt = iptpn(jp)
               iptjpt = iptpt(ipt,jpt)
               dr(1:3) = ro(1:3,ip)-ro(1:3,jp)
               call PBCr2(dr(1), dr(2), dr(3), r2)
               if (r2 < rcut2) then
                  if (r2 < r2atat(iptjpt)) then ! emulate hs overlap
                     usum = 1d10
                     fsum = 1d10
                  else if (r2 < r2umin(iptjpt)) then
                     call StopUTwoBodyACellList
                  else
                     ibuf = iubuflow(iptjpt)
                     do
                        if (r2 >= ubuf(ibuf)) exit
                        ibuf = ibuf+12
                        if (ibuf > nbuf) call StopIbuf('txptpt',iptjpt)
                     end do
                     d = r2-ubuf(ibuf)
                     usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                            d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
                     fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                            d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))
                  end if

                  u%twob(iptjpt) = u%twob(iptjpt) + usum
                  force(1:3,ip) = force(1:3,ip) + (fsum * dr(1:3))
                  force(1:3,jp) = force(1:3,jp) - (fsum * dr(1:3))
                  virtwob     = virtwob     - (fsum * r2)
               end if
            end if
            jp = ipnext(jp)
         end do
      end do
   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot         + u%twob(0)
   virial    = virial        + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   implicit none
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUTwoBodyACellList
   use MolModule, only: boxlen2, dpbc
   implicit none
   write(uout,'(a,i5)')  'ip', ip
   write(uout,'(a,i5)')  'jp', jp
   write(uout,'(a,3e15.5)') 'ro(ip)         = ', ro(1:3,ip)
   write(uout,'(a,3e15.5)') 'ro(jp)         = ', ro(1:3,jp)
   write(uout,'(a,3e15.5)') 'boxlen2        = ', boxlen2
   write(uout,'(a,3e15.5)') 'dpbc           = ', dpbc
   write(uout,'(a,i15)')    'uout           = ', uout
   write(uout,'(a,i15)')    'myid           = ', myid
   write(uout,'(a,i15)')    'iptjpt         = ', iptjpt
   write(uout,'(a,e15.5)')  'r2             = ', r2
   write(uout,'(a,e15.5)')  'r2umin(iptjpt) = ', r2umin(iptjpt)
   call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)
end subroutine StopUTwoBodyACellList

!........................................................................


end subroutine UTwoBodyACellList

!************************************************************************
!> \page energy energy.F90
!! **UTwoBodyP**
!! *calculate two-body potential energy; general particles*
!************************************************************************


subroutine UTwoBodyP

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UTwoBodyP'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf, Getnpmyid
   integer(4) :: ia, ialow, iaupp, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum, fsum, virtwob

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (lmc) then
           if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         usum = Zero
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            iat = iatan(ia)
            do ja = jalow, jaupp
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               if (r2 < r2umin(iatjat)) call StopUTwoBodyP
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
                  if (ibuf > nbuf) call StopIbuf('txatat',iatjat)
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                      d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
               fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                      d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))

               force(1,ia) = force(1,ia) + (fsum * dx)
               force(2,ia) = force(2,ia) + (fsum * dy)
               force(3,ia) = force(3,ia) + (fsum * dz)
               force(1,ja) = force(1,ja) - (fsum * dx)
               force(2,ja) = force(2,ja) - (fsum * dy)
               force(3,ja) = force(3,ja) - (fsum * dz)
               virtwob     = virtwob     - (fsum * r2)
            end do
         end do
         u%twob(iptjpt) = u%twob(iptjpt)+ usum
      end do

   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot     + u%twob(0)
   virial    = virial    + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUTwoBodyP
   character(40), parameter :: txroutine ='StopUTwoBodyP'
   write(uout,'(a,i15)')   'uout          = ', uout
   write(uout,'(a,i15)')   'myid          = ', myid
   write(*,'(a,i15)')      'myid          = ', myid
   write(uout,'(a,i15)')   'iatjat        = ', iatjat
   write(uout,'(a,i15)')   'ia            = ', ia
   write(uout,'(a,i15)')   'ja            = ', ja
   write(uout,'(a,e15.5)') 'r2            = ', r2
   write(uout,'(a,e15.5)') 'r2umin(iatjat) = ', r2umin(iatjat)
   call Stop(txroutine, 'r2 < r2umin(iatjat)', uout)
end subroutine StopUTwoBodyP

!........................................................................

end subroutine UTwoBodyP

!************************************************************************
!> \page energy energy.F90
!! **UEwald**
!! *calculate potential energy, forces, and virial from charges; k-space*
!************************************************************************


subroutine UEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UEwald'
   real(8)    :: virrec, SecondsSinceStart

   if (ltime) call CpuAdd('start', txroutine, 2, uout)
   time_ewald = SecondsSinceStart()

   if (itest == 3 .and. master) call TestUEwald(u%tot, force, virial, One/EpsiFourPi, '(real)', uout)
   if (itest == 3 .and. slave ) call TestUEwald(u%tot, force, virial, One/EpsiFourPi, '(real)_slave', uout)

! ... initiate

   u%rec          = Zero
   virrec         = Zero

! ... calculate long range term of electrostatic energy in reciprocal space

   if (txewaldrec == 'std') then
       call UEwaldRecStd
       if (lewald2dlc) call UEwaldRecStd2dlc
   else if (txewaldrec == 'spm') then
       call UEwaldRecSPM
   end if

   if (itest == 3 .and. master) call TestUEwald(u%rec, force, virrec, One, '(rec)', uout)
   if (itest == 3 .and. slave ) call TestUEwald(u%rec, force, virrec, One, '(rec)_slave', uout)

! ... calculate self interaction term of electrostatic energy in reciprocal space
   call UEwaldSelf

   if (itest == 3 .and. master) call TestUEwald(u%rec, force, virrec, One, '(rec + self)', uout)
   if (itest == 3 .and. slave ) call TestUEwald(u%rec, force, virrec, One, '(rec + self)_slave', uout)

! ... calculate surface term of electrostatic energy in reciprocal space
   if (lsurf) call  UEwaldSurf

   if (itest == 3 .and. master) call TestUEwald(u%rec, force, virrec, One, '(rec + self + sur)', uout)
   if (itest == 3 .and. slave ) call TestUEwald(u%rec, force, virrec, One, '(rec + self + sur)_slave', uout)

! ... update

   u%rec = EpsiFourPi*u%rec
   u%tot  = u%tot  + u%rec
   virial = virial + EpsiFourPi*virrec

   if (itest == 3 .and. master) call TestUEwald(u%tot, force, virial, One/EpsiFourPi, '(real + rec + self + sur)', uout)
   if (itest == 3 .and. slave ) call TestUEwald(u%tot, force, virial, One/EpsiFourPi, '(real + rec + self + sur)_slave', uout)

   time_ewald = SecondsSinceStart() - time_ewald
   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine UEwaldRecStd
   character(40), parameter :: txroutine ='UEwaldRecStd'
    integer(4) :: nx, ny, nz, kn, ia
    real(8)    :: kx, ky, kz, k2, term
    real(8)    :: kfac2, facmm, facmp, facpm, facpp

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate eikx, eiky, and eikz

   if (lmd .or. lmcall .or. lbd) call EwaldSetArray(iamyid(1), iamyid(2))
   if (lmc                     ) call EwaldSetArray(1, na)

#if defined (_PAR_)

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny)) * eikz(ia,nz)
            eikyzp(ia) =       eiky(ia,ny)  * eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle   ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn+1
            do ia = iamyid(1), iamyid(2)
               sumeikr(kn,1) = sumeikr(kn,1) + az(ia) * conjg(eikx(ia,nx)) * eikyzm(ia)
               sumeikr(kn,2) = sumeikr(kn,2) + az(ia) * conjg(eikx(ia,nx)) * eikyzp(ia)
               sumeikr(kn,3) = sumeikr(kn,3) + az(ia) *       eikx(ia,nx)  * eikyzm(ia)
               sumeikr(kn,4) = sumeikr(kn,4) + az(ia) *       eikx(ia,nx)  * eikyzp(ia)
            end do
         end do
      end do
   end do

   call par_allreduce_comps(sumeikr, eikraux, nkvec_word)

   kn = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny)) * eikz(ia,nz)
            eikyzp(ia) =       eiky(ia,ny)  * eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle   ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn+1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = az(ia) * conjg(eikx(ia,nx)) * eikyzm(ia)
               eikr(ia,2) = az(ia) * conjg(eikx(ia,nx)) * eikyzp(ia)
               eikr(ia,3) = az(ia) *       eikx(ia,nx)  * eikyzm(ia)
               eikr(ia,4) = az(ia) *       eikx(ia,nx)  * eikyzp(ia)
            end do

            term = kfac(kn)*sum (real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2)
            if (master) u%rec = u%rec + term

            kfac2 = EpsiFourPi*Two*kfac(kn)
            kx = kfac2*nx*TwoPiBoxi(1)
            ky = kfac2*ny*TwoPiBoxi(2)
            kz = kfac2*nz*TwoPiBoxi(3)
            do ia = iamyid(1), iamyid(2)
               facmm = -real(eikr(ia,1))*aimag(sumeikr(kn,1)) + aimag(eikr(ia,1))*real(sumeikr(kn,1))
               facmp = -real(eikr(ia,2))*aimag(sumeikr(kn,2)) + aimag(eikr(ia,2))*real(sumeikr(kn,2))
               facpm = -real(eikr(ia,3))*aimag(sumeikr(kn,3)) + aimag(eikr(ia,3))*real(sumeikr(kn,3))
               facpp = -real(eikr(ia,4))*aimag(sumeikr(kn,4)) + aimag(eikr(ia,4))*real(sumeikr(kn,4))
               force(1,ia) = force(1,ia) + kx*(-facmm - facmp + facpm + facpp)
               force(2,ia) = force(2,ia) + ky*(-facmm + facmp - facpm + facpp)
               force(3,ia) = force(3,ia) + kz*(+facmm + facmp + facpm + facpp)
            end do
            k2 = (nx*TwoPiBoxi(1))**2+(ny*TwoPiBoxi(2))**2+(nz*TwoPiBoxi(3))**2
            if (master) virrec = virrec - (One - k2/(Two*ualpha2))*term

         end do
      end do
   end do

#else

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny)) * eikz(ia,nz)
            eikyzp(ia) =       eiky(ia,ny)  * eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle    ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn+1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = az(ia) * conjg(eikx(ia,nx)) * eikyzm(ia)
               eikr(ia,2) = az(ia) * conjg(eikx(ia,nx)) * eikyzp(ia)
               eikr(ia,3) = az(ia) *       eikx(ia,nx)  * eikyzm(ia)
               eikr(ia,4) = az(ia) *       eikx(ia,nx)  * eikyzp(ia)
               sumeikr(kn,1) = sumeikr(kn,1) + eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + eikr(ia,4)
            end do

            term = kfac(kn)*sum (real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2)
            u%rec = u%rec + term

            kfac2 = EpsiFourPi*Two*kfac(kn)
            kx = kfac2*nx*TwoPiBoxi(1)
            ky = kfac2*ny*TwoPiBoxi(2)
            kz = kfac2*nz*TwoPiBoxi(3)
            do ia = iamyid(1), iamyid(2)
               facmm = -real(eikr(ia,1))*aimag(sumeikr(kn,1)) + aimag(eikr(ia,1))*real(sumeikr(kn,1))
               facmp = -real(eikr(ia,2))*aimag(sumeikr(kn,2)) + aimag(eikr(ia,2))*real(sumeikr(kn,2))
               facpm = -real(eikr(ia,3))*aimag(sumeikr(kn,3)) + aimag(eikr(ia,3))*real(sumeikr(kn,3))
               facpp = -real(eikr(ia,4))*aimag(sumeikr(kn,4)) + aimag(eikr(ia,4))*real(sumeikr(kn,4))
               force(1,ia) = force(1,ia) + kx*(-facmm - facmp + facpm + facpp)
               force(2,ia) = force(2,ia) + ky*(-facmm + facmp - facpm + facpp)
               force(3,ia) = force(3,ia) + kz*(+facmm + facmp + facpm + facpp)
            end do
            k2 = (nx*TwoPiBoxi(1))**2+(ny*TwoPiBoxi(2))**2+(nz*TwoPiBoxi(3))**2
            virrec = virrec - (One - k2/(Two*ualpha2))*term

         end do
      end do
   end do

#endif

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine UEwaldRecStd

!........................................................................

subroutine UEwaldRecStd2dlc

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UEwaldRecStd2dlc'
   integer(4) :: ia, nx, ny, kn
   real(8)    :: kx, ky, kp, kxx, kyy, kpp, kfac2
   real(8)    :: sinhkpz, coshkpz, term, termi

#if defined (_PAR_)
   call Stop(txroutine, 'lewald2dlc not yet implemented for _PAR_', uout)
#endif
   if (lbcrd .or. lbcto) call Stop(txroutine, 'lewald2dlc not implemented for RD and TO boundary conditions',uout)

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... calculate sinkx, coskx, sinky, and cosky

   if (lmd .or. lmcall .or. lbd) call EwaldSetArray2d(iamyid(1), iamyid(2))
   if (lmc                     ) call EwaldSetArray2d(1, na)

   kn = 0
   do ny = 0, ncut2d
      ky = ny*TwoPiBoxi(2)
      do nx = 0, ncut2d
         kx = nx*TwoPiBoxi(1)
         if (nx**2+ny**2 > ncut2d2) cycle
         if (nx == 0 .and. ny == 0) cycle
         kn = kn+1
         kp = sqrt(kx**2+ky**2)            ! k parallel to the surface

         sumtrig(kn,1) = Zero
         sumtrig(kn,2) = Zero
         sumtrig(kn,3) = Zero
         sumtrig(kn,4) = Zero
         sumtrig(kn,5) = Zero
         sumtrig(kn,6) = Zero
         sumtrig(kn,7) = Zero
         sumtrig(kn,8) = Zero
         do ia = iamyid(1), iamyid(2)
            term  = exp(kp*r(3,ia))
            termi = One/term
            sinhkpz = Half*(term - termi)
            coshkpz = Half*(term + termi)
            sumtrig(kn,1) = sumtrig(kn,1) + az(ia) * sinkx(ia,nx) * sinky(ia,ny) * sinhkpz
            sumtrig(kn,2) = sumtrig(kn,2) + az(ia) * coskx(ia,nx) * sinky(ia,ny) * sinhkpz
            sumtrig(kn,3) = sumtrig(kn,3) + az(ia) * sinkx(ia,nx) * cosky(ia,ny) * sinhkpz
            sumtrig(kn,4) = sumtrig(kn,4) + az(ia) * coskx(ia,nx) * cosky(ia,ny) * sinhkpz
            sumtrig(kn,5) = sumtrig(kn,5) + az(ia) * sinkx(ia,nx) * sinky(ia,ny) * coshkpz
            sumtrig(kn,6) = sumtrig(kn,6) + az(ia) * coskx(ia,nx) * sinky(ia,ny) * coshkpz
            sumtrig(kn,7) = sumtrig(kn,7) + az(ia) * sinkx(ia,nx) * cosky(ia,ny) * coshkpz
            sumtrig(kn,8) = sumtrig(kn,8) + az(ia) * coskx(ia,nx) * cosky(ia,ny) * coshkpz
         end do
         term = kfac2d(kn) * (-sumtrig(kn,1)**2 - sumtrig(kn,2)**2 - sumtrig(kn,3)**2 - sumtrig(kn,4)**2 &
                              +sumtrig(kn,5)**2 + sumtrig(kn,6)**2 + sumtrig(kn,7)**2 + sumtrig(kn,8)**2 )
         u%rec = u%rec + term

         kfac2 = EpsiFourPi*Two*kfac2d(kn)
         kxx = kfac2*kx
         kyy = kfac2*ky
         kpp = kfac2*kp
         do ia = iamyid(1), iamyid(2)
            term  = exp(kp*r(3,ia))
            termi = One/term
            sinhkpz = Half*(term - termi)
            coshkpz = Half*(term + termi)
            force(1,ia) = force(1,ia) - kxx*az(ia) * (-(+coskx(ia,nx)) * (+sinky(ia,ny)) * sinhkpz * sumtrig(kn,1) &
                                                      -(-sinkx(ia,nx)) * (+sinky(ia,ny)) * sinhkpz * sumtrig(kn,2) &
                                                      -(+coskx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,3) &
                                                      -(-sinkx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,4) &
                                                      +(+coskx(ia,nx)) * (+sinky(ia,ny)) * coshkpz * sumtrig(kn,5) &
                                                      +(-sinkx(ia,nx)) * (+sinky(ia,ny)) * coshkpz * sumtrig(kn,6) &
                                                      +(+coskx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,7) &
                                                      +(-sinkx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,8) )
            force(2,ia) = force(2,ia) - kyy*az(ia) * (-(+sinkx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,1) &
                                                      -(+coskx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,2) &
                                                      -(+sinkx(ia,nx)) * (-sinky(ia,ny)) * sinhkpz * sumtrig(kn,3) &
                                                      -(+coskx(ia,nx)) * (-sinky(ia,ny)) * sinhkpz * sumtrig(kn,4) &
                                                      +(+sinkx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,5) &
                                                      +(+coskx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,6) &
                                                      +(+sinkx(ia,nx)) * (-sinky(ia,ny)) * coshkpz * sumtrig(kn,7) &
                                                      +(+coskx(ia,nx)) * (-sinky(ia,ny)) * coshkpz * sumtrig(kn,8) )
            force(3,ia) = force(3,ia) - kpp*az(ia) * (-(+sinkx(ia,nx)) * (+sinky(ia,ny)) * coshkpz * sumtrig(kn,1) &
                                                      -(+coskx(ia,nx)) * (+sinky(ia,ny)) * coshkpz * sumtrig(kn,2) &
                                                      -(+sinkx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,3) &
                                                      -(+coskx(ia,nx)) * (+cosky(ia,ny)) * coshkpz * sumtrig(kn,4) &
                                                      +(+sinkx(ia,nx)) * (+sinky(ia,ny)) * sinhkpz * sumtrig(kn,5) &
                                                      +(+coskx(ia,nx)) * (+sinky(ia,ny)) * sinhkpz * sumtrig(kn,6) &
                                                      +(+sinkx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,7) &
                                                      +(+coskx(ia,nx)) * (+cosky(ia,ny)) * sinhkpz * sumtrig(kn,8) )
          end do
      end do
   end do

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

end subroutine UEwaldRecStd2dlc

!........................................................................

subroutine UEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none

   character(40), parameter :: txroutine ='UEwaldRecSPM'
   integer(4) :: ia, m, ix, iy, iz, nx, ny, nz
   real(8)    :: q, qz, qyz
   real(8)    :: splx, sply, splz
   real(8)    :: esum(3)
   real(8)    :: dsdx, dsdy, dsdz, d2sdx, d2sdy, d2sdz, Qval

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... make generalized charge mesh

   if (ltime) call CpuAdd('start', 'MakeGQMesh', 4, uout)
   Qmesh = Zero
   do ia = iamyid(1), iamyid(2)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*r(m,ia), meshmax(m,ia), spline(:,m,ia), deriv(:,m,ia), deriv2(:,m,ia))
         meshmax(m,ia) = mod(meshmax(m,ia)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         qz = az(ia)*spline(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            qyz = qz*spline(iy,2,ia)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               q = qyz*spline(ix,1,ia)
               QMesh(nx+1,ny+1,nz+1) = QMesh(nx+1,ny+1,nz+1) + q
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'MakeGQMesh', 4, uout)

#if defined (_PAR_)
   call par_allreduce_reals(QMesh, meshaux, MeshMemSize)
#endif

! ... make Fourier transformation, reciprocal space operations, and back FFT

   call SPMFFTRec(lmc, .false., .true., 'MakeFFT', 'CalcGIF', 5, u%rec, virrec)

! ... calculate potential, field, field gradient, and virial

   if (ltime) call CpuAdd('start', 'CalcProp', 4, uout)
   do ia = iamyid(1), iamyid(2)
      esum = Zero
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz =  spline(iz,3,ia)
         dsdz =  deriv(iz,3,ia)
         d2sdz = deriv2(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply =  spline(iy,2,ia)
            dsdy =   deriv(iy,2,ia)
            d2sdy = deriv2(iy,2,ia)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx =  spline(ix,1,ia)
               dsdx =   deriv(ix,1,ia)
               d2sdx = deriv2(ix,1,ia)
               Qval = QMesh(nx+1,ny+1,nz+1)

               esum(1) = esum(1) + dsdx*sply*splz*Qval
               esum(2) = esum(2) + splx*dsdy*splz*Qval
               esum(3) = esum(3) + splx*sply*dsdz*Qval
            end do
         end do
      end do
      force(1:3,ia) = force(1:3,ia) - az(ia)*esum(1:3)*dkdr(1:3)*EpsiFourPi   ! E = -dU/dr = -dU/dk * dk/dr
   end do
   if (ltime) call CpuAdd('stop', 'CalcProp', 4, uout)

   if (lmc) then
       QMesh = QMeshTM ! unmodifed particle mesh
       urecold = u%rec
   end if

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

# endif

end subroutine UEwaldRecSPM

!........................................................................

subroutine UEwaldSelf

   real(8), external :: ErfLocal
   real(8)    :: fac
   integer(4) :: ip, ipt, ia, ialow, iaupp, ja
   real(8)    :: dx, dy, dz, r2, r1, r1i, r2i, r3i, ex

! ... atomic contribution

   fac = ualpha/sqrt(Pi)
   u%rec = u%rec - fac*sum(az(iamyid(1):iamyid(2))**2)

! ... molecular contribution

   do ip = ipmyid(1), ipmyid(2)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         do ja = ia+1, iaupp
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2

            r1 = sqrt(r2)
            r1i = One/r1
            r2i = r1i**2
            ex = exp(-ualpha2*r2)
            r1i = r1i*ErfLocal(ualpha*r1)
            r3i = r2i*(r1i - ewaldfac1*ex)
            u%rec = u%rec - (az(ia)*az(ja)*r1i)

            fac = EpsiFourPi*az(ia)*az(ja)*r3i
            force(1,ia) = force(1,ia) - (fac*dx)
            force(2,ia) = force(2,ia) - (fac*dy)
            force(3,ia) = force(3,ia) - (fac*dz)
            force(1,ja) = force(1,ja) + (fac*dx)
            force(2,ja) = force(2,ja) + (fac*dy)
            force(3,ja) = force(3,ja) + (fac*dz)
         end do
       end do
    end do

end subroutine UEwaldSelf

!........................................................................

subroutine UEwaldSurf

   real(8)    :: fac, sumqrx, sumqry, sumqrz, term
   if (.not.lewald2dlc) then                    ! 3d-periodic system
      fac = TwoPi/(Three*vol)
      sumqrx = sum(az(1:na)*r(1,1:na))
      sumqry = sum(az(1:na)*r(2,1:na))
      sumqrz = sum(az(1:na)*r(3,1:na))
      term = fac*(sumqrx**2+sumqry**2+sumqrz**2)
      if (master) u%rec = u%rec + term
      force(1,iamyid(1):iamyid(2)) = force(1,iamyid(1):iamyid(2)) - EpsiFourPi*(Two*fac)*az(iamyid(1):iamyid(2))*sumqrx
      force(2,iamyid(1):iamyid(2)) = force(2,iamyid(1):iamyid(2)) - EpsiFourPi*(Two*fac)*az(iamyid(1):iamyid(2))*sumqry
      force(3,iamyid(1):iamyid(2)) = force(3,iamyid(1):iamyid(2)) - EpsiFourPi*(Two*fac)*az(iamyid(1):iamyid(2))*sumqrz
      if (master) virrec = virrec - term
   else                                        ! 2d-periodic system
      fac = TwoPi/vol
      sumqrz = sum(az(1:na)*r(3,1:na))
      term = fac*sumqrz**2
      if (master) u%rec = u%rec + term
      force(3,iamyid(1):iamyid(2)) = force(3,iamyid(1):iamyid(2)) - EpsiFourPi*(Two*fac)*az(iamyid(1):iamyid(2))*sumqrz
      if (master) virrec = virrec - term
   end if
end subroutine UEwaldSurf

!........................................................................

subroutine TestUEwald(u, force, virial, norm, label, unit)
   real(8),      intent(in) :: u
   real(8),      intent(in) :: force(3,*)
   real(8),      intent(in) :: virial
   real(8),      intent(in) :: norm
   character(*), intent(in) :: label
   integer(4),   intent(in) :: unit
   integer(4) :: ia
   call WriteHead(3, 'Test'//trim(txroutine)//'  '//label, unit)
   write(unit,'(a,2x,f12.7)') 'energy', u*norm
   write(unit,'(a,2x,f12.7)') 'virial', virial*norm
   write(unit,'(a)')  'atom                   force'
   write(unit,'((i5,2x,3f12.7))') (ia, force(1:3,ia)/EpsiFourPi, ia=1,min(6,int(na)))
end subroutine TestUEwald

!........................................................................

end subroutine UEwald

!************************************************************************
!> \page energy energy.F90
!! **EwaldSetArray**
!! *calculate eikx, eiky, and eikz arrays used for ewald summation; k-space*
!************************************************************************


subroutine EwaldSetArray(ialow, iaupp)

   use EnergyModule
   implicit none

   integer(4), intent(in) :: ialow, iaupp

   integer(4) :: ia, icut

   do ia = ialow, iaupp
      eikx(ia,0) = cmplx(One,Zero)
      eiky(ia,0) = cmplx(One,Zero)
      eikz(ia,0) = cmplx(One,Zero)
      eikx(ia,1) = cmplx(cos(TwoPiBoxi(1)*r(1,ia)),sin(TwoPiBoxi(1)*r(1,ia)))
      eiky(ia,1) = cmplx(cos(TwoPiBoxi(2)*r(2,ia)),sin(TwoPiBoxi(2)*r(2,ia)))
      eikz(ia,1) = cmplx(cos(TwoPiBoxi(3)*r(3,ia)),sin(TwoPiBoxi(3)*r(3,ia)))
   end do
   do icut = 2, ncut
      do ia = ialow, iaupp
         eikx(ia,icut) = eikx(ia,icut-1)*eikx(ia,1)
         eiky(ia,icut) = eiky(ia,icut-1)*eiky(ia,1)
         eikz(ia,icut) = eikz(ia,icut-1)*eikz(ia,1)
      end do
   end do

end subroutine EwaldSetArray

!************************************************************************
!> \page energy energy.F90
!! **EwaldSetArray2d**
!! *calculate sinkx, coskx, sinky, and cosky arrays used for ewald summation; k-space*
!************************************************************************


subroutine EwaldSetArray2d(ialow, iaupp)

   use EnergyModule
   implicit none

   integer(4), intent(in) :: ialow, iaupp
   integer(4) :: ia, icut

   do ia = ialow, iaupp
      sinkx(ia,0) = One             ! ok!!!
      coskx(ia,0) = One
      sinky(ia,0) = One             ! ok!!!
      cosky(ia,0) = One
      sinkx(ia,1) = sin(TwoPiBoxi(1)*r(1,ia))
      coskx(ia,1) = cos(TwoPiBoxi(1)*r(1,ia))
      sinky(ia,1) = sin(TwoPiBoxi(2)*r(2,ia))
      cosky(ia,1) = cos(TwoPiBoxi(2)*r(2,ia))
   end do

   do icut = 2, ncut2d
      do ia = ialow, iaupp
         sinkx(ia,icut) = sinkx(ia,icut-1)*coskx(ia,1) + coskx(ia,icut-1)*sinkx(ia,1)
         coskx(ia,icut) = coskx(ia,icut-1)*coskx(ia,1) - sinkx(ia,icut-1)*sinkx(ia,1)
         sinky(ia,icut) = sinky(ia,icut-1)*cosky(ia,1) + cosky(ia,icut-1)*sinky(ia,1)
         cosky(ia,icut) = cosky(ia,icut-1)*cosky(ia,1) - sinky(ia,icut-1)*sinky(ia,1)
      end do
   end do

end subroutine EwaldSetArray2d

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UDipole**
!! *calculate potential energy from charges and dipoles*
!************************************************************************


!     the code presuposses that charges are given in elementary units

subroutine UDipole

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDipole'
   integer(4) :: ia
   real(8)    :: forcex, forcey, forcez, torquex, torquey, torquez
   real(8)    :: ureal

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... calculate potential, field, field gradient, and virial from charges and static dipoles

   call UDipoleP

! ... temporary: save real space contribution of the electrostatic potential energy

   ureal = Zero
   do ia = iamyid(1), iamyid(2)
      ureal = ureal + az(ia)*potstat(ia) - (dip(1,ia)*estat(1,ia) + dip(2,ia)*estat(2,ia) + dip(3,ia)*estat(3,ia))
   end do
   ureal = EpsiFourPi*Half*ureal

   if (lewald) call UDipoleEwald

#if defined (_PAR_)
   if (ltime) call CpuAdd('start', 'comm', 3, uout)
   call par_allreduce_reals(potstat, vaux, na  )
   call par_allreduce_reals(estat,   vaux, 3*na)
   call par_allreduce_reals(efg    , vaux, 6*na)
   if (ltime) call CpuAdd('stop', 'comm', 3, uout)
#endif

    u%stat = Zero

    do ia = iamyid(1), iamyid(2)

! ... calculate electrostatic potential energy: u%stat = q pot - u e
!                                                                 a a

      u%stat = u%stat + az(ia)*potstat(ia) - (dip(1,ia)*estat(1,ia) + dip(2,ia)*estat(2,ia) + dip(3,ia)*estat(3,ia))

! ... calculate forces: f  = q e    + u    efg
!                        a      a      b      ab

      forcex = az(ia)*estat(1,ia) + (dip(1,ia)*efg(1,ia) + dip(2,ia)*efg(4,ia) + dip(3,ia)*efg(5,ia))
      forcey = az(ia)*estat(2,ia) + (dip(1,ia)*efg(4,ia) + dip(2,ia)*efg(2,ia) + dip(3,ia)*efg(6,ia))
      forcez = az(ia)*estat(3,ia) + (dip(1,ia)*efg(5,ia) + dip(2,ia)*efg(6,ia) + dip(3,ia)*efg(3,ia))
      force(1,ia)  = force(1,ia)  + EpsiFourPi*forcex
      force(2,ia)  = force(2,ia)  + EpsiFourPi*forcey
      force(3,ia)  = force(3,ia)  + EpsiFourPi*forcez

! ... calculate torques: t = u x e

      torquex = dip(2,ia)*estat(3,ia) - dip(3,ia)*estat(2,ia)
      torquey = dip(3,ia)*estat(1,ia) - dip(1,ia)*estat(3,ia)
      torquez = dip(1,ia)*estat(2,ia) - dip(2,ia)*estat(1,ia)
      torque(1,ia) = torque(1,ia) + EpsiFourPi*torquex
      torque(2,ia) = torque(2,ia) + EpsiFourPi*torquey
      torque(3,ia) = torque(3,ia) + EpsiFourPi*torquez

   end do

   u%stat = EpsiFourPi*Half*u%stat

! ... temporary: get reciprocal space contribution of the electrostatic potential energy

   u%rec = u%stat - ureal

! ... update

   u%tot  = u%tot  + u%stat
   virial = virial + EpsiFourPi*virelec

   if (itest == 3 .and. master) call TestUDipole(uout)

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUDipole(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a)')  'atom                     r'
   write(unit,'((i5,(2x,3f12.7)))') (ia, r(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                 force'
   write(unit,'((i5,1(2x,3f12.7)))') (ia, force(1:3,ia)/EpsiFourPi, ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                 torque'
   write(unit,'((i5,1(2x,3f12.7)))') (ia, torque(1:3,ia)/EpsiFourPi, ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'u%stat    ', u%stat/EpsiFourPi
   write(unit,'(a,2x,f12.7)') 'virelec  ', virelec
   write(unit,'(a,2x,f12.7)') 'virialtot', virial/EpsiFourPi
end subroutine TestUDipole

!........................................................................

end subroutine UDipole

!************************************************************************
!> \page energy energy.F90
!! **UDipoleP**
!! *calculate potential, field, and field gradient from charges and dipoles*
!************************************************************************


subroutine UDipoleP

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDipoleP'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, Getnpmyid
   integer(4) :: ia, ialow, iaupp, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc
   real(8)    :: r1, r2, r1i, r2i, r3i, r5i, r7i, threer5i, ex
   real(8)    :: doti, dotj, dotir5i, dotjr5i, dotir7i, dotjr7i
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
   real(8), external :: ErfLocal

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... initiate

   potstat(1:na)   = Zero
   estat(1:3,1:na) = Zero
   efg(1:6,1:na)   = Zero
   virelec         = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (lmc) then
           if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2

               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
                  r7i = r2i*(r5i + ewaldfac3*ex)
               else if (lrf) then
                  r3i = r2i*r1i
                  r5i = r2i*r3i
                  r7i = r2i*r5i
                  r1i = r1i + rffac2*r2
                  r3i = r3i - rffac           ! not to be used for field gradient
               else
                  r3i = r2i*r1i
                  r5i = r2i*r3i
                  r7i = r2i*r5i
               end if
               threer5i = Three*r5i

               doti = dip(1,ia)*dx+dip(2,ia)*dy+dip(3,ia)*dz
               dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
               dotir5i = three*r5i*doti
               dotjr5i = three*r5i*dotj
               dotir7i = (Three*Five*r7i)*doti
               dotjr7i = (Three*Five*r7i)*dotj

! ... calculate potential (potstat)

!                                  3
!     pot(r) = q / r  +  u   r  / r
!                          a   a
!                                  3
!     pot(0) = q / r  -  u   r  / r
!                          a   a

               potstat(ia) = potstat(ia) + (az(ja)*r1i + dotj*r3i)
               potstat(ja) = potstat(ja) + (az(ia)*r1i - doti*r3i)

! ... calculate field (estat)

!                   -3                         2    -5
!     e (r) = q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab
!
!                   -3                         2    -5
!     e (r) =-q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab

! ... calculate virial

!     vir = -q e  r
!               a  a

               fldx = (+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx
               fldy = (+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy
               fldz = (+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz
               estat(1,ia) = estat(1,ia) + fldx
               estat(2,ia) = estat(2,ia) + fldy
               estat(3,ia) = estat(3,ia) + fldz
               virelec = virelec - az(ia)*(fldx*dx + fldy*dy + fldz*dz)
               fldx = (-az(ia)*dx - dip(1,ia))*r3i + dotir5i*dx
               fldy = (-az(ia)*dy - dip(2,ia))*r3i + dotir5i*dy
               fldz = (-az(ia)*dz - dip(3,ia))*r3i + dotir5i*dz
               estat(1,ja) = estat(1,ja) + fldx
               estat(2,ja) = estat(2,ja) + fldy
               estat(3,ja) = estat(3,ja) + fldz

! ... calculate field gradient (stat)

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    -  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    +  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

! ... calculate virial

!     vir  = -u  e   r
!              a  ab  b

               efgxx = az(ja)*r3i + (-az(ja)*dx**2 + Two*dip(1,ja)*dx + dotj)*threer5i - dotjr7i*dx**2
               efgyy = az(ja)*r3i + (-az(ja)*dy**2 + Two*dip(2,ja)*dy + dotj)*threer5i - dotjr7i*dy**2
               efgzz = az(ja)*r3i + (-az(ja)*dz**2 + Two*dip(3,ja)*dz + dotj)*threer5i - dotjr7i*dz**2
               efgxy =              (-az(ja)*dx*dy + dip(1,ja)*dy + dip(2,ja)*dx)*threer5i - dotjr7i*dx*dy
               efgxz =              (-az(ja)*dx*dz + dip(1,ja)*dz + dip(3,ja)*dx)*threer5i - dotjr7i*dx*dz
               efgyz =              (-az(ja)*dy*dz + dip(2,ja)*dz + dip(3,ja)*dy)*threer5i - dotjr7i*dy*dz
               efg(1,ia) = efg(1,ia) + efgxx
               efg(2,ia) = efg(2,ia) + efgyy
               efg(3,ia) = efg(3,ia) + efgzz
               efg(4,ia) = efg(4,ia) + efgxy
               efg(5,ia) = efg(5,ia) + efgxz
               efg(6,ia) = efg(6,ia) + efgyz
               virelec  = virelec -(dip(1,ia)*efgxx*dx + dip(1,ia)*efgxy*dy + dip(1,ia)*efgxz*dz &
                                  + dip(2,ia)*efgxy*dx + dip(2,ia)*efgyy*dy + dip(2,ia)*efgyz*dz &
                                  + dip(3,ia)*efgxz*dx + dip(3,ia)*efgyz*dy + dip(3,ia)*efgzz*dz)
               efgxx = az(ia)*r3i + (-az(ia)*dx**2 - Two*dip(1,ia)*dx - doti)*threer5i + dotir7i*dx**2
               efgyy = az(ia)*r3i + (-az(ia)*dy**2 - Two*dip(2,ia)*dy - doti)*threer5i + dotir7i*dy**2
               efgzz = az(ia)*r3i + (-az(ia)*dz**2 - Two*dip(3,ia)*dz - doti)*threer5i + dotir7i*dz**2
               efgxy =              (-az(ia)*dx*dy - dip(1,ia)*dy - dip(2,ia)*dx)*threer5i + dotir7i*dx*dy
               efgxz =              (-az(ia)*dx*dz - dip(1,ia)*dz - dip(3,ia)*dx)*threer5i + dotir7i*dx*dz
               efgyz =              (-az(ia)*dy*dz - dip(2,ia)*dz - dip(3,ia)*dy)*threer5i + dotir7i*dy*dz
               efg(1,ja) = efg(1,ja) + efgxx
               efg(2,ja) = efg(2,ja) + efgyy
               efg(3,ja) = efg(3,ja) + efgzz
               efg(4,ja) = efg(4,ja) + efgxy
               efg(5,ja) = efg(5,ja) + efgxz
               efg(6,ja) = efg(6,ja) + efgyz

            end do
         end do
      end do

! ... add contribution from center molecule

      if (lrf) then
         do ia = ialow, iaupp
            do ja = ialow, iaupp
               dx = r(1,ja)-r(1,ia)
               dy = r(2,ja)-r(2,ia)
               dz = r(3,ja)-r(3,ia)
               potstat(ia) = potstat(ia) + rffac2*az(ja)*(dx**2+dy**2+dz**2) &
                                         + rffac*(dx*dip(1,ja)+dy*dip(2,ja)+dz*dip(3,ja))
               estat(1,ia) = estat(1,ia) + rffac*(az(ja)*dx + dip(1,ja))
               estat(2,ia) = estat(2,ia) + rffac*(az(ja)*dy + dip(2,ja))
               estat(3,ia) = estat(3,ia) + rffac*(az(ja)*dz + dip(3,ja))
            end do
         end do
      end if

   end do

   if (itest == 3 .and. master) call TestUDipoleP(uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine TestUDipoleP(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'  (real)', unit)
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'virial', virelec
end subroutine TestUDipoleP

!........................................................................

end subroutine UDipoleP

!************************************************************************
!> \page energy energy.F90
!! **UDipoleEwald**
!! *calculate potential, field, and field gradient from charges and dipoles; k-space*
!************************************************************************


subroutine UDipoleEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDipoleEwald'
   real(8) :: SecondsSinceStart

   if (ltime) call CpuAdd('start', txroutine, 3, uout)
   time_ewald = SecondsSinceStart()

   if (txewaldrec == 'std') then
       call UDipoleEwaldRecStd
   else if (txewaldrec == 'spm') then
       call UDipoleEwaldRecSPM
   end if
   if (itest == 3 .and. master) call TestUDipoleEwald('(real + rec)', uout)
   call UDipoleEwaldSelf
   if (itest == 3 .and. master) call TestUDipoleEwald('(real + rec + self)', uout)
   if (lsurf) call UDipoleEwaldSurf
   if (itest == 3 .and. master) call TestUDipoleEwald('(real + rec + self + sur)', uout)

   time_ewald = SecondsSinceStart() - time_ewald
   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine UDipoleEwaldRecStd
   character(40), parameter :: txroutine ='UDipoleEwaldRecStd'
   integer(4) :: nx, ny, nz, kn, ia
   real(8)    :: kfac2, facmm, facmp, facpm, facpp
   real(8)    :: termx, termy, termz
   real(8)    :: kx, ky, kz, kxkx, kyky, kzkz, kxky, kxkz, kykz
   complex(8) :: cmm, cmp, cpm, cpp

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... calculate eikx, eiky, and eikz

   if (lmd .or. lmcall .or. lbd) call EwaldSetArray(iamyid(1), iamyid(2))
   if (lmc                     ) call EwaldSetArray(1, na)

#if defined (_PAR_)

   sumeikr = cmplx(Zero,Zero)
   sumeikrd = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            kx = TwoPiBoxi(1)*nx
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*dip(1,ia)
               termy = ky*dip(2,ia)
               termz = kz*dip(3,ia)
               facmm = -termx - termy + termz
               facmp = -termx + termy + termz
               facpm = +termx - termy + termz
               facpp = +termx + termy + termz
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), facmm) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), facmp) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), facpm) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), facpp) * eikr(ia,4)
               sumeikrd(kn,1) = sumeikrd(kn,1) + cmplx(Zero, facmm) * eikr(ia,1)
               sumeikrd(kn,2) = sumeikrd(kn,2) + cmplx(Zero, facmp) * eikr(ia,2)
               sumeikrd(kn,3) = sumeikrd(kn,3) + cmplx(Zero, facpm) * eikr(ia,3)
               sumeikrd(kn,4) = sumeikrd(kn,4) + cmplx(Zero, facpp) * eikr(ia,4)
            end do
         end do
      end do
   end do

   call par_allreduce_comps(sumeikr, eikraux, nkvec_word)
   call par_allreduce_comps(sumeikrd, eikraux, nkvec_word)

   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         ky = TwoPiBoxi(2)*ny
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            kx = TwoPiBoxi(1)*nx
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
            end do

            kfac2 = Two*kfac(kn)
            kxkx = kfac2*kx**2
            kyky = kfac2*ky**2
            kzkz = kfac2*kz**2
            kxky = kfac2*kx*ky
            kxkz = kfac2*kx*kz
            kykz = kfac2*ky*kz
            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               potstat(ia) = potstat(ia) + kfac2*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))

               estat(1,ia) = estat(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               estat(2,ia) = estat(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               estat(3,ia) = estat(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))

               efg(1,ia) = efg(1,ia) + kxkx*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(2,ia) = efg(2,ia) + kyky*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(3,ia) = efg(3,ia) + kzkz*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(4,ia) = efg(4,ia) + kxky*(+real(cmm) - real(cmp) - real(cpm) + real(cpp))
               efg(5,ia) = efg(5,ia) + kxkz*(-real(cmm) - real(cmp) + real(cpm) + real(cpp))
               efg(6,ia) = efg(6,ia) + kykz*(-real(cmm) + real(cmp) - real(cpm) + real(cpp))
            end do

            if (master) then
               termx = sum( real(sumeikrd(kn,1:4)*real(sumeikr(kn,1:4))) + aimag(sumeikrd(kn,1:4)*aimag(sumeikr(kn,1:4))) )
               termy = sum( real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2 )
               termz = One-(kx**2+ky**2+kz**2)/(Two*ualpha2)
               virelec = virelec - kfac(kn)*(Two*termx + termy*termz)
            end if

         end do
      end do
   end do

#else

   sumeikr = cmplx(Zero,Zero)
   sumeikrd = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*dip(1,ia)
               termy = ky*dip(2,ia)
               termz = kz*dip(3,ia)
               facmm = -termx - termy + termz
               facmp = -termx + termy + termz
               facpm = +termx - termy + termz
               facpp = +termx + termy + termz
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), facmm) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), facmp) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), facpm) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), facpp) * eikr(ia,4)
               sumeikrd(kn,1) = sumeikrd(kn,1) + cmplx(Zero, facmm) * eikr(ia,1)
               sumeikrd(kn,2) = sumeikrd(kn,2) + cmplx(Zero, facmp) * eikr(ia,2)
               sumeikrd(kn,3) = sumeikrd(kn,3) + cmplx(Zero, facpm) * eikr(ia,3)
               sumeikrd(kn,4) = sumeikrd(kn,4) + cmplx(Zero, facpp) * eikr(ia,4)
            end do

            kfac2 = Two*kfac(kn)
            kxkx = kfac2*kx**2
            kyky = kfac2*ky**2
            kzkz = kfac2*kz**2
            kxky = kfac2*kx*ky
            kxkz = kfac2*kx*kz
            kykz = kfac2*ky*kz
            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               potstat(ia) = potstat(ia) + kfac2*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))

               estat(1,ia) = estat(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               estat(2,ia) = estat(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               estat(3,ia) = estat(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))

               efg(1,ia) = efg(1,ia) + kxkx*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(2,ia) = efg(2,ia) + kyky*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(3,ia) = efg(3,ia) + kzkz*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(4,ia) = efg(4,ia) + kxky*(+real(cmm) - real(cmp) - real(cpm) + real(cpp))
               efg(5,ia) = efg(5,ia) + kxkz*(-real(cmm) - real(cmp) + real(cpm) + real(cpp))
               efg(6,ia) = efg(6,ia) + kykz*(-real(cmm) + real(cmp) - real(cpm) + real(cpp))
            end do

            termx = sum( real(sumeikrd(kn,1:4)*real(sumeikr(kn,1:4))) + aimag(sumeikrd(kn,1:4)*aimag(sumeikr(kn,1:4))) )
            termy = sum( real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2 )
            termz = One-(kx**2+ky**2+kz**2)/(Two*ualpha2)
            virelec = virelec - kfac(kn)*(Two*termx + termy*termz)

         end do
      end do
   end do

#endif

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

end subroutine UDipoleEwaldRecStd

!........................................................................

subroutine UDipoleEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none

   character(40), parameter :: txroutine ='UDipoleEwaldRecSPM'
   integer(4) :: m, ia, ix, iy, iz, nx, ny, nz
   real(8)    :: q, q1, dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
   real(8)    :: splx, sply, splz
   real(8)    :: psum, esum(3), efgsum(6)
   real(8)    :: dsdx, dsdy, dsdz, d2sdx, d2sdy, d2sdz, Qval, Qdip

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... make generalized charge mesh

   if (ltime) call CpuAdd('start', 'MakeGQMesh', 5, uout)
   Qmesh = Zero
   do ia = iamyid(1), iamyid(2)
      dipx2 = dip(1,ia)*dkdr(1)
      dipy2 = dip(2,ia)*dkdr(2)
      dipz2 = dip(3,ia)*dkdr(3)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*r(m,ia), meshmax(m,ia), spline(:,m,ia), deriv(:,m,ia), deriv2(:,m,ia))
         meshmax(m,ia) = mod(meshmax(m,ia)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz = spline(iz,3,ia)
         dsdz = deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply = spline(iy,2,ia)
            q1 = az(ia)*sply*splz
            dipx1 = dipx2*sply*splz
            dipy1 = dipy2*deriv(iy,2,ia)*splz
            dipz1 = dipz2*sply*dsdz
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx = spline(ix,1,ia)
               q = q1*splx
               dipx = dipx1*deriv(ix,1,ia)
               dipy = dipy1*splx
               dipz = dipz1*splx
               QMesh(nx+1,ny+1,nz+1) = &
                    QMesh(nx+1,ny+1,nz+1) + (q + dipx + dipy + dipz)
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'MakeGQMesh', 5, uout)

#if defined (_PAR_)
   if (ltime) call CpuAdd('start', 'MeshComm', 5, uout)
   call par_allreduce_reals(QMesh, meshaux, MeshMemSize)
   if (ltime) call CpuAdd('stop', 'MeshComm', 5, uout)
#endif

    if (lmc) QMeshTM = QMesh

! ... make Fourier transformation, reciprocal space operations, and back FFT

   call SPMFFTRec(.false., lmc, .true., 'MakeFFT', 'CalcGIF', 5, urecold, virelec)

   ! ... calculate potential, field, field gradient, and virial

   if (ltime) call CpuAdd('start', 'CalcProp', 5, uout)
   do ia = iamyid(1), iamyid(2)
      dipx2 = dip(1,ia)*dkdr(1)
      dipy2 = dip(2,ia)*dkdr(2)
      dipz2 = dip(3,ia)*dkdr(3)
       psum = Zero
       esum = Zero
       efgsum = Zero
       do iz = 0,order-1
           nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
           splz =  spline(iz,3,ia)
           dsdz =  deriv(iz,3,ia)
           d2sdz = deriv2(iz,3,ia)
           do iy = 0,order-1
               ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
               sply =  spline(iy,2,ia)
               dsdy =   deriv(iy,2,ia)
               d2sdy = deriv2(iy,2,ia)
               do ix = 0,order-1
                   nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
                   splx =  spline(ix,1,ia)
                   dsdx =   deriv(ix,1,ia)
                   d2sdx = deriv2(ix,1,ia)
                   Qval = QMesh(nx+1,ny+1,nz+1)

                   psum = psum + splx*sply*splz*Qval

                   esum(1) = esum(1) + dsdx*sply*splz*Qval
                   esum(2) = esum(2) + splx*dsdy*splz*Qval
                   esum(3) = esum(3) + splx*sply*dsdz*Qval

                   efgsum(1) = efgsum(1) + d2sdx*sply*splz*Qval
                   efgsum(2) = efgsum(2) + splx*d2sdy*splz*Qval
                   efgsum(3) = efgsum(3) + splx*sply*d2sdz*Qval
                   efgsum(4) = efgsum(4) + dsdx*dsdy*splz*Qval
                   efgsum(5) = efgsum(5) + dsdx*sply*dsdz*Qval
                   efgsum(6) = efgsum(6) + splx*dsdy*dsdz*Qval

                   Qdip = dipx2*dsdx*sply*splz + dipy2*splx*dsdy*splz + dipz2*splx*sply*dsdz

                   virelec = virelec - Qdip*Qval
               end do
           end do
       end do
       potstat(ia) = potstat(ia) + psum
       estat(1:3,ia) = estat(1:3,ia) - esum(1:3)*dkdr(1:3)               ! E = -dU/dr = -dU/dk * dk/dr
       efg(1:6,ia) = efg(1:6,ia) - efgsum(1:6)*d2kdr(1:6)
   end do
   if (ltime) call CpuAdd('stop', 'CalcProp', 5, uout)

   if (lmc) QMesh = QMeshTM

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

# endif

end subroutine UDipoleEwaldRecSPM

!........................................................................

subroutine UDipoleEwaldSelf

   integer(4) :: ip, ipt, ia, ialow, iaupp, ja
   real(8)    :: dx, dy, dz, r2, r1, r1i, r2i, r3i, r5i, r7i, ex, threer5i
   real(8)    :: doti, dotj, dotir5i, dotjr5i, dotjr7i, dotir7i
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
   real(8), external :: ErfLocal

! ... atomic contribution

   potstat(iamyid(1):iamyid(2)) = potstat(iamyid(1):iamyid(2)) - ewaldselffac1*az(iamyid(1):iamyid(2))
   estat(1,iamyid(1):iamyid(2)) = estat(1,iamyid(1):iamyid(2)) + ewaldselffac2*dip(1,iamyid(1):iamyid(2))
   estat(2,iamyid(1):iamyid(2)) = estat(2,iamyid(1):iamyid(2)) + ewaldselffac2*dip(2,iamyid(1):iamyid(2))
   estat(3,iamyid(1):iamyid(2)) = estat(3,iamyid(1):iamyid(2)) + ewaldselffac2*dip(3,iamyid(1):iamyid(2))
   efg(1,iamyid(1):iamyid(2)) = efg(1,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))
   efg(2,iamyid(1):iamyid(2)) = efg(2,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))
   efg(3,iamyid(1):iamyid(2)) = efg(3,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))

! ... molcular contribution (_not_ constant)

   do ip = ipmyid(1), ipmyid(2)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         do ja = ia+1, iaupp
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2

            r1 = sqrt(r2)
            r1i = One/r1
            r2i = r1i**2
            ex = exp(-ualpha2*r2)
            r1i = r1i*ErfLocal(ualpha*r1)
            r3i = r2i*(r1i - ewaldfac1*ex)
            r5i = r2i*(r3i - ewaldfac2*ex)
            r7i = r2i*(r5i - ewaldfac3*ex)
            threer5i = Three*r5i

            doti = dip(1,ia)*dx+dip(2,ia)*dy+dip(3,ia)*dz
            dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
            dotir5i = three*r5i*doti
            dotjr5i = three*r5i*dotj
            dotir7i = (Three*Five*r7i)*doti
            dotjr7i = (Three*Five*r7i)*dotj

            potstat(ia) = potstat(ia) - (az(ja)*r1i + dotj*r3i)
            potstat(ja) = potstat(ja) - (az(ia)*r1i - doti*r3i)

            fldx = -((az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx)
            fldy = -((az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy)
            fldz = -((az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz)
            estat(1,ia) = estat(1,ia) + fldx
            estat(2,ia) = estat(2,ia) + fldy
            estat(3,ia) = estat(3,ia) + fldz
            virelec = virelec - (az(ia)*(fldx*dx + fldy*dy + fldz*dz))
            fldx = -((-az(ia)*dx - dip(1,ia))*r3i + dotir5i*dx)
            fldy = -((-az(ia)*dy - dip(2,ia))*r3i + dotir5i*dy)
            fldz = -((-az(ia)*dz - dip(3,ia))*r3i + dotir5i*dz)
            estat(1,ja) = estat(1,ja) + fldx
            estat(2,ja) = estat(2,ja) + fldy
            estat(3,ja) = estat(3,ja) + fldz

            efgxx = -( az(ja)*r3i + (-az(ja)*dx**2 + Two*dip(1,ja)*dx + dotj)*threer5i - dotjr7i*dx**2 )
            efgyy = -( az(ja)*r3i + (-az(ja)*dy**2 + Two*dip(2,ja)*dy + dotj)*threer5i - dotjr7i*dy**2 )
            efgzz = -( az(ja)*r3i + (-az(ja)*dz**2 + Two*dip(3,ja)*dz + dotj)*threer5i - dotjr7i*dz**2 )
            efgxy = -( (-az(ja)*dx*dy + dip(1,ja)*dy + dip(2,ja)*dx)*threer5i - dotjr7i*dx*dy )
            efgxz = -( (-az(ja)*dx*dz + dip(1,ja)*dz + dip(3,ja)*dx)*threer5i - dotjr7i*dx*dz )
            efgyz = -( (-az(ja)*dy*dz + dip(2,ja)*dz + dip(3,ja)*dy)*threer5i - dotjr7i*dy*dz )
            efg(1,ia) = efg(1,ia) + efgxx
            efg(2,ia) = efg(2,ia) + efgyy
            efg(3,ia) = efg(3,ia) + efgzz
            efg(4,ia) = efg(4,ia) + efgxy
            efg(5,ia) = efg(5,ia) + efgxz
            efg(6,ia) = efg(6,ia) + efgyz
            virelec = virelec - ( dip(1,ia)*efgxx*dx + dip(1,ia)*efgxy*dy + dip(1,ia)*efgxz*dz   &
                                + dip(2,ia)*efgxy*dx + dip(2,ia)*efgyy*dy + dip(2,ia)*efgyz*dz   &
                                + dip(3,ia)*efgxz*dx + dip(3,ia)*efgyz*dy + dip(3,ia)*efgzz*dz )
            efgxx = -( az(ia)*r3i + (-az(ia)*dx**2 - Two*dip(1,ia)*dx - doti)*threer5i + dotir7i*dx**2 )
            efgyy = -( az(ia)*r3i + (-az(ia)*dy**2 - Two*dip(2,ia)*dy - doti)*threer5i + dotir7i*dy**2 )
            efgzz = -( az(ia)*r3i + (-az(ia)*dz**2 - Two*dip(3,ia)*dz - doti)*threer5i + dotir7i*dz**2 )
            efgxy = -( (-az(ia)*dx*dy - dip(1,ia)*dy - dip(2,ia)*dx)*threer5i + dotir7i*dx*dy )
            efgxz = -( (-az(ia)*dx*dz - dip(1,ia)*dz - dip(3,ia)*dx)*threer5i + dotir7i*dx*dz )
            efgyz = -( (-az(ia)*dy*dz - dip(2,ia)*dz - dip(3,ia)*dy)*threer5i + dotir7i*dy*dz )
            efg(1,ja) = efg(1,ja) + efgxx
            efg(2,ja) = efg(2,ja) + efgyy
            efg(3,ja) = efg(3,ja) + efgzz
            efg(4,ja) = efg(4,ja) + efgxy
            efg(5,ja) = efg(5,ja) + efgxz
            efg(6,ja) = efg(6,ja) + efgyz
         end do
      end do
   end do
end subroutine UDipoleEwaldSelf

!........................................................................

subroutine UDipoleEwaldSurf
   real(8)    :: fac, sumqrx, sumqry, sumqrz, sumdx, sumdy, sumdz
   fac = FourPiThird/vol
   sumqrx = sum(az(1:na)*r(1,1:na))
   sumqry = sum(az(1:na)*r(2,1:na))
   sumqrz = sum(az(1:na)*r(3,1:na))
   sumdx = sum(dip(1,1:na))
   sumdy = sum(dip(2,1:na))
   sumdz = sum(dip(3,1:na))
   potstat(iamyid(1):iamyid(2)) = potstat(iamyid(1):iamyid(2)) +                &
                              fac*(r(1,iamyid(1):iamyid(2))*(sumqrx+sumdx) +    &
                                   r(2,iamyid(1):iamyid(2))*(sumqry+sumdy) +    &
                                   r(3,iamyid(1):iamyid(2))*(sumqrz+sumdz))
   estat(1,iamyid(1):iamyid(2)) = estat(1,iamyid(1):iamyid(2)) - fac*(sumqrx + sumdx)
   estat(2,iamyid(1):iamyid(2)) = estat(2,iamyid(1):iamyid(2)) - fac*(sumqry + sumdy)
   estat(3,iamyid(1):iamyid(2)) = estat(3,iamyid(1):iamyid(2)) - fac*(sumqrz + sumdz)
   if (master)  virelec  = virelec - TwoPi/(Three*vol) *                         &
                    (        (sumqrx**2 + sumqry**2 + sumqrz**2)                &
                     +   Two*(Two*(sumqrx*sumdx + sumqry*sumdy + sumqrz*sumdz)) &
                     + Three*(sumdx**2 + sumdy**2 + sumdz**2 )    )
end subroutine UDipoleEwaldSurf

!........................................................................

subroutine TestUDipoleEwald(label, unit)
   integer(4), intent(in) :: unit
   character(*), intent(in) :: label
   integer(4) :: ia
   call WriteHead(3, 'Test'//trim(txroutine)//'  '//label, unit)
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'virelec', virelec
end subroutine TestUDipoleEwald

!........................................................................

end subroutine UDipoleEwald

!************************************************************************
!> \page energy energy.F90
!! **UManyBodyP**
!! *calculate many-body potential energy; general particles*
!************************************************************************


!     includes also:
!     two-body electrostatic contribution
!     two-body and many-body reciprocal contribution when ewald is set
!     reaction field contribution from dielectrica when lrf is set

!     the code presuposses that charges are given in elementary units and polarizabilities in volume units

subroutine UManyBodyP

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UManyBodyP'
   integer(4), save :: icounter = 0    ! force the first prediction of idm to be iterative
   logical, save    :: first = .true.  ! force the first iteratative prediction of idm to invoke at least two iterations
   integer(4) :: ia
   real(8)    :: forcex, forcey, forcez, torquex, torquey, torquez

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... calculate potential and field from charges and static dipoles

   call FieldStat

   if (lewald) call FieldStatEwald

#if defined (_PAR_)
   if (ltime) call CpuAdd('start', 'comm', 3, uout)
   call par_allreduce_reals(potstat, vaux, na  )
   call par_allreduce_reals(estat, vaux, 3*na)
   if (ltime) call CpuAdd('stop', 'comm', 3, uout)
#endif

   if (itest == 3 .and. master) call TestUManyBodyP1(uout)

! ............ calculate induced dipole moments ...........

   icounter = mod(icounter,npolit)
   if (icounter == 0) then          ! iterate
      call IterIdm(first)
      if (.not.lidmconv) call Warn('IterIdm', 'iteration failed on max iteration', uout)
   else                            ! predict
      idm(1:3,1:na) = Two*idm1(1:3,1:na) - idm2(1:3,1:na)
   end if
   icounter = icounter + 1

! ... calculate total dipole moment

   diptot(1:3,1:na) = dip(1:3,1:na) + idm(1:3,1:na)

! ... calculate field, field gradient, and virial from charges and total dipoles

   call FieldTot

   if (lewald) call FieldTotEwald

#if defined (_PAR_)
   if (ltime) call CpuAdd('start', 'comm', 3, uout)
   call par_allreduce_reals(etot   , vaux, 3*na)
   call par_allreduce_reals(efg    , vaux, 6*na)
   if (ltime) call CpuAdd('stop', 'comm', 3, uout)
#endif

   u%stat = Zero
   u%pol = Zero

   do ia = iamyid(1), iamyid(2)

!                                                                               0
! ... calculate potential energy from charges and static dipoles: u%stat = q pot - u e
!                                                                                   a a

      u%stat = u%stat + az(ia)*potstat(ia) - (dip(1,ia)*estat(1,ia) + dip(2,ia)*estat(2,ia) + dip(3,ia)*estat(3,ia))

!                                                  ind 0
! ... calculate polarization energy: u%pol = -0.5 u   e
!                                                  a   a

      u%pol = u%pol - (idm(1,ia)*estat(1,ia) + idm(2,ia)*estat(2,ia) + idm(3,ia)*estat(3,ia))

!                               tot    tot
! ... calculate forces: f  = q e    + u    efg
!                        a      a      b      ab

      forcex = az(ia)*etot(1,ia) + diptot(1,ia)*efg(1,ia) + diptot(2,ia)*efg(4,ia) + diptot(3,ia)*efg(5,ia)
      forcey = az(ia)*etot(2,ia) + diptot(1,ia)*efg(4,ia) + diptot(2,ia)*efg(2,ia) + diptot(3,ia)*efg(6,ia)
      forcez = az(ia)*etot(3,ia) + diptot(1,ia)*efg(5,ia) + diptot(2,ia)*efg(6,ia) + diptot(3,ia)*efg(3,ia)
      force(1,ia)  = force(1,ia)  + EpsiFourPi*forcex
      force(2,ia)  = force(2,ia)  + EpsiFourPi*forcey
      force(3,ia)  = force(3,ia)  + EpsiFourPi*forcez

!                             tot    tot
! ... calculate torques: t = u    x e

      torquex = diptot(2,ia)*etot(3,ia) - diptot(3,ia)*etot(2,ia)
      torquey = diptot(3,ia)*etot(1,ia) - diptot(1,ia)*etot(3,ia)
      torquez = diptot(1,ia)*etot(2,ia) - diptot(2,ia)*etot(1,ia)
      torque(1,ia) = torque(1,ia) + EpsiFourPi*torquex
      torque(2,ia) = torque(2,ia) + EpsiFourPi*torquey
      torque(3,ia) = torque(3,ia) + EpsiFourPi*torquez

   end do

   u%stat = EpsiFourPi*Half*u%stat
   u%pol = EpsiFourPi*Half*u%pol

! ... update

   u%tot  = u%tot  + (u%stat + u%pol)
   virial = virial + EpsiFourPi*virelec

! ... correct the induced moments with the total fields and update idm1 and idm2

   do ia = 1, na
      idm(1,ia) = poltens(1,ia)*etot(1,ia) + poltens(4,ia)*etot(2,ia) + poltens(5,ia)*etot(3,ia)
      idm(2,ia) = poltens(4,ia)*etot(1,ia) + poltens(2,ia)*etot(2,ia) + poltens(6,ia)*etot(3,ia)
      idm(3,ia) = poltens(5,ia)*etot(1,ia) + poltens(6,ia)*etot(2,ia) + poltens(3,ia)*etot(3,ia)
      if (first) then
         idm2(1,ia) = idm(1,ia)
         idm2(2,ia) = idm(2,ia)
         idm2(3,ia) = idm(3,ia)
      else
         idm2(1,ia) = idm1(1,ia)
         idm2(2,ia) = idm1(2,ia)
         idm2(3,ia) = idm1(3,ia)
      endif
      idm1(1,ia) = idm(1,ia)
      idm1(2,ia) = idm(2,ia)
      idm1(3,ia) = idm(3,ia)
   end do

   first =.false.

! ... calculate induced dipole moment: total and particle

   call CalcIndDipMom

   if (itest == 3 .and. master) call TestUManyBodyP2(uout)

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUManyBodyP1(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1', unit)
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestUManyBodyP1

!........................................................................

subroutine TestUManyBodyP2(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a)')  'atom                     r                                   dip                                  idm'
   write(unit,'((i5,3(2x,3f12.7)))') (ia, r(1:3,ia), dip(1:3,ia), idm(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom      potstat'
   write(unit,'((i5,2x,f12.7))') (ia, potstat(ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                   estat                                 etot'
   write(unit,'((i5,2(2x,3f12.7)))') (ia, estat(1:3,ia), etot(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                 force'
   write(unit,'((i5,1(2x,3f12.7)))') (ia, force(1:3,ia)/EpsiFourPi, ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                 torque'
   write(unit,'((i5,1(2x,3f12.7)))') (ia, torque(1:3,ia)/EpsiFourPi, ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'u%stat    ', u%stat/EpsiFourPi
   write(unit,'(a,2x,f12.7)') 'u%pol     ', u%pol/EpsiFourPi
   write(unit,'(a,2x,f12.7)') 'virelec  ', virelec
   write(unit,'(a,2x,f12.7)') 'virialtot', virial/EpsiFourPi
end subroutine TestUManyBodyP2

!........................................................................

end subroutine UManyBodyP

!************************************************************************
!> \page energy energy.F90
!! **FieldStat**
!! *calculate potential and field from charges and static dipoles*
!************************************************************************


subroutine FieldStat

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldStat'
   real(8), parameter :: OneSixth = 1.0d0/6.0d0
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, Getnpmyid
   integer(4) :: ia, ialow, iaupp, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc
   real(8)    :: r1, r2, r1i, r2i, r3i, r5i, ex
   real(8)    :: doti, dotj, dotir5i, dotjr5i
   real(8)    :: s, v, d2, d1, d0
   real(8), external :: ErfLocal

   if (ltime) call CpuAdd('start', txroutine,  3, uout)

! ... initiate

   potstat(1:na)   = Zero
   estat(1:3,1:na) = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2

               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
               else if (lrf) then
                  r3i = r2i*r1i
                  r5i = r3i*r2i
                  r1i = r1i + rffac2*r2
                  r3i = r3i - rffac
               else
                  r3i = r2i*r1i
                  r5i = r3i*r2i
               end if

               if (ldamping) then
                  s  = 1.4*1.662*(poltens(1,ia)*poltens(1,ja))**(OneSixth)
                  v  = min(One,r1/s)
                  d2 = v**4
                  d1 = 4*v**3 - 3*d2
                  d0 = d2 - 2*v**3 + 2*v
                  r1i = d0*r1i
                  r3i = d1*r3i
                  r5i = d2*r5i
               end if

               doti = dip(1,ia)*dx+dip(2,ia)*dy+dip(3,ia)*dz
               dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
               dotir5i = (Three*r5i) * doti
               dotjr5i = (Three*r5i) * dotj

! ... calculate potential (potstat)

!                                  3
!     pot(r) = q / r  +  u   r  / r
!                          a   a
!                                  3
!     pot(0) = q / r  -  u   r  / r
!                          a   a

               potstat(ia) = potstat(ia) + (az(ja)*r1i + dotj*r3i)
               potstat(ja) = potstat(ja) + (az(ia)*r1i - doti*r3i)

! ... calculate field (estat)

!                   -3                         2    -5
!     e (r) = q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab
!
!                   -3                         2    -5
!     e (r) =-q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab


               estat(1,ia) = estat(1,ia) + (+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx
               estat(2,ia) = estat(2,ia) + (+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy
               estat(3,ia) = estat(3,ia) + (+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz
               estat(1,ja) = estat(1,ja) + (-az(ia)*dx - dip(1,ia))*r3i + dotir5i*dx
               estat(2,ja) = estat(2,ja) + (-az(ia)*dy - dip(2,ia))*r3i + dotir5i*dy
               estat(3,ja) = estat(3,ja) + (-az(ia)*dz - dip(3,ia))*r3i + dotir5i*dz

            end do
         end do
      end do

! ... add contributions from center molecule

      if (lrf) then
         do ia = ialow, iaupp
            do ja = ialow, iaupp
               dx = r(1,ja)-r(1,ia)
               dy = r(2,ja)-r(2,ia)
               dz = r(3,ja)-r(3,ia)
               potstat(ia) = potstat(ia) + rffac2*az(ja)*(dx**2+dy**2+dz**2) &
                                         + rffac*(dx*dip(1,ja)+dy*dip(2,ja)+dz*dip(3,ja))
               estat(1,ia) = estat(1,ia) + rffac*(az(ja)*dx + dip(1,ja))
               estat(2,ia) = estat(2,ia) + rffac*(az(ja)*dy + dip(2,ja))
               estat(3,ia) = estat(3,ia) + rffac*(az(ja)*dz + dip(3,ja))
            end do
         end do
      end if

   end do

   if (itest == 3 .and. master) call TestFieldStat(uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine TestFieldStat(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'  (real)', unit)
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestFieldStat

!........................................................................

end subroutine FieldStat

!************************************************************************
!> \page energy energy.F90
!! **FieldStatEwald**
!! *calculate potential and field from charges and static dipoles; k-space*
!************************************************************************


subroutine FieldStatEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldStatEwald'

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

   if (txewaldrec == 'std') then
       call FieldStatEwaldRecStd
   else if (txewaldrec == 'spm') then
       call FieldStatEwaldRecSPM
   end if
   if (itest == 3 .and. master) call TestFieldStatEwald('(real + rec)', uout)
   call FieldStatEwaldSelf
   if (itest == 3 .and. master) call TestFieldStatEwald('(real + rec + self)', uout)
   if (lsurf) call FieldStatEwaldSurf
   if (itest == 3 .and. master) call TestFieldStatEwald('(real + rec + self + sur)', uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine FieldStatEwaldRecStd
   character(40), parameter :: txroutine ='FieldStatEwaldRecStd'
   integer(4) :: nx, ny, nz, kn, ia
   real(8)    :: kfac2, kx, ky, kz, termx, termy, termz
   complex(8) :: cmm, cmp, cpm, cpp

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... calculate eikx, eiky, and eikz

   call EwaldSetArray(iamyid(1), iamyid(2))

#if defined (_PAR_)

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*dip(1,ia)
               termy = ky*dip(2,ia)
               termz = kz*dip(3,ia)
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), -termx - termy + termz) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), -termx + termy + termz) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), +termx - termy + termz) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), +termx + termy + termz) * eikr(ia,4)
            end do
         end do
      end do
   end do

   call par_allreduce_comps(sumeikr, eikraux, nkvec_word)

   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
            end do

            kfac2 = Two*kfac(kn)

            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               potstat(ia) = potstat(ia) + kfac2*(+real(cmm)+real(cmp)+real(cpm)+real(cpp))
               estat(1,ia) = estat(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               estat(2,ia) = estat(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               estat(3,ia) = estat(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))
            end do

         end do
      end do
   end do

#else

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*dip(1,ia)
               termy = ky*dip(2,ia)
               termz = kz*dip(3,ia)
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), -termx - termy + termz) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), -termx + termy + termz) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), +termx - termy + termz) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), +termx + termy + termz) * eikr(ia,4)
            end do

            kfac2 = Two*kfac(kn)

            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               potstat(ia) = potstat(ia) + kfac2*(+real(cmm)+real(cmp)+real(cpm)+real(cpp))
               estat(1,ia) = estat(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               estat(2,ia) = estat(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               estat(3,ia) = estat(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))
            end do

         end do
      end do
   end do

#endif

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

end subroutine FieldStatEwaldRecStd

!........................................................................

subroutine FieldStatEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none

   character(40), parameter :: txroutine ='FieldStatEwaldRecSPM'
   integer(4) :: ia, m, ix, iy, iz, nx, ny, nz
   real(8)    :: q, q1, dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
   real(8)    :: splx, sply, splz
   real(8)    :: psum, esum(3), efgsum(6)
   real(8)    :: dsdx, dsdy, dsdz, Qval
   real(8)    :: raux2

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... make generalized charge mesh

   if (ltime) call CpuAdd('start', 'StatMakeGQMesh', 5, uout)
   Qmesh = Zero
   do ia = 1, na
      dipx2 = dip(1,ia)*dkdr(1)
      dipy2 = dip(2,ia)*dkdr(2)
      dipz2 = dip(3,ia)*dkdr(3)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*r(m,ia), meshmax(m,ia), spline(:,m,ia), deriv(:,m,ia), deriv2(:,m,ia))
         meshmax(m,ia) = mod(meshmax(m,ia)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz = spline(iz,3,ia)
         dsdz = deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply = spline(iy,2,ia)
            q1 = az(ia)*sply*splz
            dipx1 = dipx2*sply*splz
            dipy1 = dipy2*deriv(iy,2,ia)*splz
            dipz1 = dipz2*sply*dsdz
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx = spline(ix,1,ia)
               q = q1*splx
               dipx = dipx1*deriv(ix,1,ia)
               dipy = dipy1*splx
               dipz = dipz1*splx
               QMesh(nx+1,ny+1,nz+1) = &
                    QMesh(nx+1,ny+1,nz+1) + (q + dipx + dipy + dipz)
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'StatMakeGQMesh', 5, uout)

! ... make Fourier transformation, reciprocal space operations, and back FFT

   call SPMFFTRec(.false., .false., .false., 'StatMakeFFT', 'StatCalcGIF', 5, raux, raux2)

! ... calculate potential, field, field gradient, and virial

   if (ltime) call CpuAdd('start', 'StatCalcProp', 5, uout)
   do ia=1,na
      psum = Zero
      esum = Zero
      efgsum = Zero
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz =  spline(iz,3,ia)
         dsdz =  deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply =  spline(iy,2,ia)
            dsdy =   deriv(iy,2,ia)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx =  spline(ix,1,ia)
               dsdx =   deriv(ix,1,ia)
               Qval = QMesh(nx+1,ny+1,nz+1)
               psum = psum + splx*sply*splz*Qval
               esum(1) = esum(1) + dsdx*sply*splz*Qval
               esum(2) = esum(2) + splx*dsdy*splz*Qval
               esum(3) = esum(3) + splx*sply*dsdz*Qval
            end do
         end do
      end do
      potstat(ia) = potstat(ia) + psum
      estat(1:3,ia) = estat(1:3,ia) - esum(1:3)*dkdr(1:3)               ! E = -dU/dr = -dU/dk * dk/dr
   end do
   if (ltime) call CpuAdd('stop', 'StatCalcProp', 5, uout)

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

# endif

end subroutine FieldStatEwaldRecSPM

!........................................................................

subroutine FieldStatEwaldSelf

   real(8), external :: ErfLocal
   integer(4) :: ip, ipt, ia, ialow, iaupp, ja
   real(8)    :: dx, dy, dz, ex, r1, r2, r1i, r2i, r3i, r5i, doti, dotj, dotir5i, dotjr5i

! ... atomic contribution

   potstat(iamyid(1):iamyid(2)) = potstat(iamyid(1):iamyid(2)) - ewaldselffac1*az(iamyid(1):iamyid(2))
   estat(1,iamyid(1):iamyid(2)) = estat(1,iamyid(1):iamyid(2)) + ewaldselffac2*dip(1,iamyid(1):iamyid(2))
   estat(2,iamyid(1):iamyid(2)) = estat(2,iamyid(1):iamyid(2)) + ewaldselffac2*dip(2,iamyid(1):iamyid(2))
   estat(3,iamyid(1):iamyid(2)) = estat(3,iamyid(1):iamyid(2)) + ewaldselffac2*dip(3,iamyid(1):iamyid(2))

! ...  molecular contribution (constant)

   do ip = ipmyid(1), ipmyid(2)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         do ja = ia+1, iaupp
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2
            r1 = sqrt(r2)
            r1i = One/r1
            r2i = r1i**2
            ex = exp(-ualpha2*r2)
            r1i = r1i*ErfLocal(ualpha*r1)
            r3i = r2i*(r1i - ewaldfac1*ex)
            r5i = r2i*(r3i - ewaldfac2*ex)

            doti = dip(1,ia)*dx+dip(2,ia)*dy+dip(3,ia)*dz
            dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
            dotir5i = (Three*r5i)*doti
            dotjr5i = (Three*r5i)*dotj

            potstat(ia) = potstat(ia) - (az(ja)*r1i + dotj*r3i)
            potstat(ja) = potstat(ja) - (az(ia)*r1i - doti*r3i)

            estat(1,ia) = estat(1,ia) - ((+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx)
            estat(2,ia) = estat(2,ia) - ((+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy)
            estat(3,ia) = estat(3,ia) - ((+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz)
            estat(1,ja) = estat(1,ja) - ((-az(ia)*dx - dip(1,ia))*r3i + dotir5i*dx)
            estat(2,ja) = estat(2,ja) - ((-az(ia)*dy - dip(2,ia))*r3i + dotir5i*dy)
            estat(3,ja) = estat(3,ja) - ((-az(ia)*dz - dip(3,ia))*r3i + dotir5i*dz)
         end do
      end do
   end do

end subroutine FieldStatEwaldSelf

!........................................................................

subroutine FieldStatEwaldSurf

   real(8)    :: fac, sumqrx, sumqry, sumqrz, sumdx, sumdy, sumdz

   fac = FourPiThird/vol
   sumqrx = sum(az(1:na)*r(1,1:na))
   sumqry = sum(az(1:na)*r(2,1:na))
   sumqrz = sum(az(1:na)*r(3,1:na))
   sumdx = sum(dip(1,1:na))
   sumdy = sum(dip(2,1:na))
   sumdz = sum(dip(3,1:na))

   potstat(iamyid(1):iamyid(2)) = potstat(iamyid(1):iamyid(2)) +                         &
                                  fac*(r(1,iamyid(1):iamyid(2))*(sumqrx+sumdx) +         &
                                       r(2,iamyid(1):iamyid(2))*(sumqry+sumdy) +         &
                                       r(3,iamyid(1):iamyid(2))*(sumqrz+sumdz))
   estat(1,iamyid(1):iamyid(2)) = estat(1,iamyid(1):iamyid(2)) - fac*(sumqrx + sumdx)
   estat(2,iamyid(1):iamyid(2)) = estat(2,iamyid(1):iamyid(2)) - fac*(sumqry + sumdy)
   estat(3,iamyid(1):iamyid(2)) = estat(3,iamyid(1):iamyid(2)) - fac*(sumqrz + sumdz)

end subroutine FieldStatEwaldSurf

!........................................................................

subroutine TestFieldStatEwald(label,unit)
   character(*), intent(in) :: label
   integer(4),   intent(in) :: unit
   integer(4) :: ia
   call WriteHead(3, 'Test'//trim(txroutine)//'  '//label, unit)
   write(unit,'(a)')  'atom      potstat                    estat'
   write(unit,'(i5,2x,f12.7,2x,3f12.7)') (ia, potstat(ia), estat(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestFieldStatEwald

!........................................................................

end subroutine FieldStatEwald

!************************************************************************
!> \page energy energy.F90
!! **IterIdm**
!! *iterate induced dipole moments self-consistently*
!************************************************************************

!     include reciprocal contribution if lewald = .true.
!     include reaction field contribution if lrf = .true.

!     induced moment =  polarizability*( field from charges
!                                     +  field from static dipole moments
!                                     +  field from induced dipole moments )

subroutine IterIdm(lforceiter)

   use EnergyModule
   implicit none

   logical, intent(in) :: lforceiter        !=.ture. => force an additonal iteration

   character(40), parameter :: txroutine ='IterIdm'
   integer(4) :: iter, ia
   real(8)    :: ex, ey, ez, idmdiff, idmsize

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... loop over induced dipole iterations

   do iter = 1, mpolit

! ... calculate field from induced dipoles

      call FieldIdm

      if (lewald) call FieldIdmEwald

#if defined (_PAR_)
      if (ltime) call CpuAdd('start', 'comm', 4, uout)
      call par_allreduce_reals(eidm, vaux, 3*na)
      if (ltime) call CpuAdd('stop', 'comm', 4, uout)
#endif

      if (itest == 3 .and. master) call TestIterIdm1(uout)

! ... construct new induced dipole moments

      do ia = 1, na
         ex = estat(1,ia) + eidm(1,ia)
         ey = estat(2,ia) + eidm(2,ia)
         ez = estat(3,ia) + eidm(3,ia)
         vaux(1,ia) = poltens(1,ia)*ex + poltens(4,ia)*ey + poltens(5,ia)*ez
         vaux(2,ia) = poltens(4,ia)*ex + poltens(2,ia)*ey + poltens(6,ia)*ez
         vaux(3,ia) = poltens(5,ia)*ex + poltens(6,ia)*ey + poltens(3,ia)*ez
      end do

! ... test loop for no convergence

      lidmconv =.true.

      do ia = 1, na

         if (.not.lapolarization(ia)) cycle

! ... get difference of the idm and the size of the old idm

         idmdiff = sqrt( (idm(1,ia)-vaux(1,ia))**2 + (idm(2,ia)-vaux(2,ia))**2 + (idm(3,ia)-vaux(3,ia))**2 )
         idmsize = sqrt( (idm(1,ia))**2            + (idm(2,ia))**2            + (idm(3,ia))**2 )

! ... check if not converged, then exit the test loop
         if(idmsize == Zero) then
            lidmconv =.false.
            exit
         end if
         if ((idmdiff/idmsize > tpolit) .or. (lforceiter.and.(iter==1)) ) then
            lidmconv =.false.
            exit
         end if

      end do

! ... update idm

      idm(1:3,1:na) = vaux(1:3,1:na)

      if (itest == 3 .and. master) call TestIterIdm2(uout)

! ... check if converged, then exit the iteration loop

      if (lidmconv) exit

   end do

   if (iter > mpolit) lidmconv =.false.

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine TestIterIdm1(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'1', unit)
   write(unit,'(a)')  'atom                    eidm'
   write(unit,'(i5,2x,3f12.7)') (ia, eidm(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestIterIdm1

!........................................................................

subroutine TestIterIdm2(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'2', unit)
   write(unit,'(a)')  'atom                     idm'
   write(unit,'(i5,2x,3f12.7)') (ia, idm(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestIterIdm2

!........................................................................

end subroutine IterIdm

!************************************************************************
!> \page energy energy.F90
!! **FieldIdm**
!! *calculate field from induced dipoles*
!************************************************************************


subroutine FieldIdm

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldIdm'
   real(8), parameter :: OneSixth = 1.0d0/6.0d0
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, Getnpmyid
   integer(4) :: ia, ialow, iaupp, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r1, r2, r1i, r2i, r3i, r5i, ex
   real(8)    :: doti, dotj
   real(8)    :: s, v, d2, d1
   real(8), external :: ErfLocal

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... initiate

   eidm(1:3,1:na) = Zero

! ... calculate field from induced dipoles

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2

               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
               else if (lrf) then
                  r3i = r2i*r1i
                  r5i = r3i*r2i
                  r3i = r3i - rffac
               else
                  r3i = r2i*r1i
                  r5i = r3i*r2i
               end if

               if (ldamping) then
                  s  = 1.4*1.662*(poltens(1,ia)*poltens(1,ja))**(OneSixth)
                  v  = min(One,r1/s)
                  d2 = v**4
                  d1 = 4*v**3 - 3*d2
                  r3i = d1*r3i
                  r5i = d2*r5i
               end if

               doti = (Three*r5i)*(idm(1,ia)*dx+idm(2,ia)*dy+idm(3,ia)*dz)
               dotj = (Three*r5i)*(idm(1,ja)*dx+idm(2,ja)*dy+idm(3,ja)*dz)

               eidm(1,ia) = eidm(1,ia) + (dotj*dx - idm(1,ja)*r3i)
               eidm(2,ia) = eidm(2,ia) + (dotj*dy - idm(2,ja)*r3i)
               eidm(3,ia) = eidm(3,ia) + (dotj*dz - idm(3,ja)*r3i)
               eidm(1,ja) = eidm(1,ja) + (doti*dx - idm(1,ia)*r3i)
               eidm(2,ja) = eidm(2,ja) + (doti*dy - idm(2,ia)*r3i)
               eidm(3,ja) = eidm(3,ja) + (doti*dz - idm(3,ia)*r3i)

            end do
         end do
      end do

! ... add contribution from center molecule

      if (lrf) then
         do ia = ialow, iaupp
            do ja = ialow, iaupp
               eidm(1,ia) = eidm(1,ia) + rffac*idm(1,ja)
               eidm(2,ia) = eidm(2,ia) + rffac*idm(2,ja)
               eidm(3,ia) = eidm(3,ia) + rffac*idm(3,ja)
            end do
         end do
      end if

   end do

   if (itest == 3 .and. master) call TestFieldIdm(uout)

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

contains

!........................................................................

subroutine TestFieldIdm(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'  (real)', unit)
   write(unit,'(a)')  'atom                    eidm'
   write(unit,'(i5,2x,3f12.7)') (ia, eidm(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestFieldIdm

!........................................................................

end subroutine FieldIdm

!************************************************************************
!> \page energy energy.F90
!! **FieldIdmEwald**
!! *calculate field from induced dipoles; k-space*
!************************************************************************

!     eikx, eiky, and eikz are presumed calculated already

subroutine FieldIdmEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldIdmEwald'

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

   if (txewaldrec == 'std') then
       call FieldIdmEwaldRecStd
   else if (txewaldrec == 'spm') then
       call FieldIdmEwaldRecSPM
   end if
   call FieldIdmEwaldSelf
   if (lsurf) call FieldIdmEwaldSurf
   if (itest == 3 .and. master) call TestFieldIdmEwald(uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine FieldIdmEwaldRecStd

   character(40), parameter :: txroutine ='FieldIdmEwaldRecStd'
   integer(4) :: nx, ny, nz, kn, ia
   real(8)    :: kfac2, kx, ky, kz
   real(8)    :: termx, termy, termz, facmm, facmp, facpm, facpp

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

#if defined (_PAR_)

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*idm(1,ia)
               termy = ky*idm(2,ia)
               termz = kz*idm(3,ia)
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(Zero, (-termx - termy + termz)) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(Zero, (-termx + termy + termz)) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(Zero, (+termx - termy + termz)) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(Zero, (+termx + termy + termz)) * eikr(ia,4)
            end do
         end do
      end do
   end do

   call par_allreduce_comps(sumeikr, eikraux, nkvec_word)

   kn = 0
   do nz = 0, ncut
      do ny = 0, ncut
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            ky = TwoPiBoxi(2)*ny
            kz = TwoPiBoxi(3)*nz
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
            end do

            kfac2 = Two*kfac(kn)

            do ia = iamyid(1), iamyid(2)
               facmm = -real(eikr(ia,1))*aimag(sumeikr(kn,1)) + aimag(eikr(ia,1))*real(sumeikr(kn,1))
               facmp = -real(eikr(ia,2))*aimag(sumeikr(kn,2)) + aimag(eikr(ia,2))*real(sumeikr(kn,2))
               facpm = -real(eikr(ia,3))*aimag(sumeikr(kn,3)) + aimag(eikr(ia,3))*real(sumeikr(kn,3))
               facpp = -real(eikr(ia,4))*aimag(sumeikr(kn,4)) + aimag(eikr(ia,4))*real(sumeikr(kn,4))
               eidm(1,ia) = eidm(1,ia) + kfac2*kx*(-facmm - facmp + facpm + facpp)
               eidm(2,ia) = eidm(2,ia) + kfac2*ky*(-facmm + facmp - facpm + facpp)
               eidm(3,ia) = eidm(3,ia) + kfac2*kz*(+facmm + facmp + facpm + facpp)
            end do

         end do
      end do
   end do

#else

   sumeikr = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            kx = TwoPiBoxi(1)*nx
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*idm(1,ia)
               termy = ky*idm(2,ia)
               termz = kz*idm(3,ia)
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(Zero, (-termx - termy + termz)) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(Zero, (-termx + termy + termz)) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(Zero, (+termx - termy + termz)) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(Zero, (+termx + termy + termz)) * eikr(ia,4)
            end do

            kfac2 = Two*kfac(kn)

            do ia = iamyid(1), iamyid(2)
               facmm = -real(eikr(ia,1))*aimag(sumeikr(kn,1)) + aimag(eikr(ia,1))*real(sumeikr(kn,1))
               facmp = -real(eikr(ia,2))*aimag(sumeikr(kn,2)) + aimag(eikr(ia,2))*real(sumeikr(kn,2))
               facpm = -real(eikr(ia,3))*aimag(sumeikr(kn,3)) + aimag(eikr(ia,3))*real(sumeikr(kn,3))
               facpp = -real(eikr(ia,4))*aimag(sumeikr(kn,4)) + aimag(eikr(ia,4))*real(sumeikr(kn,4))
               eidm(1,ia) = eidm(1,ia) + kfac2*kx*(-facmm - facmp + facpm + facpp)
               eidm(2,ia) = eidm(2,ia) + kfac2*ky*(-facmm + facmp - facpm + facpp)
               eidm(3,ia) = eidm(3,ia) + kfac2*kz*(+facmm + facmp + facpm + facpp)
            end do

         end do
      end do
   end do

#endif

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)


end subroutine FieldIdmEwaldRecStd

!........................................................................

subroutine FieldIdmEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none

   character(40), parameter :: txroutine ='FieldIdmEwaldRecSPM'
   integer(4) :: ia, m, ix, iy, iz, nx, ny, nz
   real(8)    :: dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
   real(8)    :: splx, sply, splz
   real(8)    :: psum, esum(3), efgsum(6)
   real(8)    :: dsdx, dsdy, dsdz, Qval
   real(8)    :: raux2

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... make generalized charge mesh

   if (ltime) call CpuAdd('start', 'IdmMakeGQMesh', 5, uout)
   Qmesh = Zero
   do ia = 1, na
      dipx2 = idm(1,ia)*dkdr(1)
      dipy2 = idm(2,ia)*dkdr(2)
      dipz2 = idm(3,ia)*dkdr(3)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*r(m,ia), meshmax(m,ia), spline(:,m,ia), deriv(:,m,ia), deriv2(:,m,ia))
         meshmax(m,ia) = mod(meshmax(m,ia)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz = spline(iz,3,ia)
         dsdz = deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply = spline(iy,2,ia)
            dipx1 = dipx2*sply*splz
            dipy1 = dipy2*deriv(iy,2,ia)*splz
            dipz1 = dipz2*sply*dsdz
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx = spline(ix,1,ia)
               dipx = dipx1*deriv(ix,1,ia)
               dipy = dipy1*splx
               dipz = dipz1*splx
               QMesh(nx+1,ny+1,nz+1) = &
   !                QMesh(nx+1,ny+1,nz+1) + q + dipx + dipy + dipz
                    QMesh(nx+1,ny+1,nz+1) + dipx + dipy + dipz
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'IdmMakeGQMesh', 5, uout)

! ... make Fourier transformation, reciprocal space operations, and back FFT

   call SPMFFTRec(.false., .false., .false.,'IdmMakeFFT', 'IdmCalcGIF', 5, raux, raux2)

! ... calculate potential, field, field gradient, and virial

   if (ltime) call CpuAdd('start', 'IdmCalcProp', 5, uout)
   do ia=1,na
      psum = Zero
      esum = Zero
      efgsum = Zero
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz =  spline(iz,3,ia)
         dsdz =  deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply =  spline(iy,2,ia)
            dsdy =   deriv(iy,2,ia)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx =  spline(ix,1,ia)
               dsdx =   deriv(ix,1,ia)
               Qval = QMesh(nx+1,ny+1,nz+1)
               esum(1) = esum(1) + dsdx*sply*splz*Qval
               esum(2) = esum(2) + splx*dsdy*splz*Qval
               esum(3) = esum(3) + splx*sply*dsdz*Qval
            end do
         end do
      end do
      eidm(1:3,ia) = eidm(1:3,ia) - esum(1:3)*dkdr(1:3)               ! E = -dU/dr = -dU/dk * dk/dr
   end do
   if (ltime) call CpuAdd('stop', 'IdmCalcProp', 5, uout)

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

# endif

end subroutine FieldIdmEwaldRecSPM

!........................................................................

subroutine FieldIdmEwaldSelf

   real(8), external :: ErfLocal
   integer(4) :: ip, ipt, ia, ialow, iaupp, ja
   real(8)    :: dx, dy, dz, ex, r1, r2, r1i, r2i, r3i, r5i, doti, dotj

! ... atomic contribution

   eidm(1,iamyid(1):iamyid(2)) = eidm(1,iamyid(1):iamyid(2)) + ewaldselffac2*idm(1,iamyid(1):iamyid(2))
   eidm(2,iamyid(1):iamyid(2)) = eidm(2,iamyid(1):iamyid(2)) + ewaldselffac2*idm(2,iamyid(1):iamyid(2))
   eidm(3,iamyid(1):iamyid(2)) = eidm(3,iamyid(1):iamyid(2)) + ewaldselffac2*idm(3,iamyid(1):iamyid(2))

! ... molecular contribution (_not_ constant)

   do ip = ipmyid(1), ipmyid(2)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         do ja = ia+1, iaupp
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2

            r1 = sqrt(r2)
            r1i = One/r1
            r2i = r1i**2
            ex = exp(-ualpha2*r2)
            r1i = r1i*ErfLocal(ualpha*r1)
            r3i = r2i*(r1i - ewaldfac1*ex)
            r5i = r2i*(r3i - ewaldfac2*ex)

            doti = (Three*r5i)*(idm(1,ia)*dx+idm(2,ia)*dy+idm(3,ia)*dz)
            dotj = (Three*r5i)*(idm(1,ja)*dx+idm(2,ja)*dy+idm(3,ja)*dz)

            eidm(1,ia) = eidm(1,ia) - (dotj*dx - idm(1,ja)*r3i)
            eidm(2,ia) = eidm(2,ia) - (dotj*dy - idm(2,ja)*r3i)
            eidm(3,ia) = eidm(3,ia) - (dotj*dz - idm(3,ja)*r3i)
            eidm(1,ja) = eidm(1,ja) - (doti*dx - idm(1,ia)*r3i)
            eidm(2,ja) = eidm(2,ja) - (doti*dy - idm(2,ia)*r3i)
            eidm(3,ja) = eidm(3,ja) - (doti*dz - idm(3,ia)*r3i)
         end do
      end do
   end do
end subroutine FieldIdmEwaldSelf

!........................................................................

subroutine FieldIdmEwaldSurf
   real(8)    :: fac, sumdx, sumdy, sumdz
   fac = FourPiThird/vol
   sumdx = sum(idm(1,1:na))
   sumdy = sum(idm(2,1:na))
   sumdz = sum(idm(3,1:na))
   eidm(1,iamyid(1):iamyid(2)) = eidm(1,iamyid(1):iamyid(2)) - fac*sumdx
   eidm(2,iamyid(1):iamyid(2)) = eidm(2,iamyid(1):iamyid(2)) - fac*sumdy
   eidm(3,iamyid(1):iamyid(2)) = eidm(3,iamyid(1):iamyid(2)) - fac*sumdz
end subroutine FieldIdmEwaldSurf

!........................................................................

subroutine TestFieldIdmEwald(unit)
   integer(4), intent(in) :: unit
   integer(4) :: ia
   call WriteHead(3, 'Test'//trim(txroutine)//'  (real + rec + self + sur)', unit)
   write(unit,'(a)')  'atom                   eidm'
   write(unit,'(i5,2x,3f12.7)') (ia, eidm(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestFieldIdmEwald

!........................................................................

end subroutine FieldIdmEwald

!************************************************************************
!> \page energy energy.F90
!! **FieldTot**
!! *calculate field, field gradient, and virial from charges and total dipoles*
!************************************************************************


subroutine FieldTot

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldTot'
   real(8), parameter :: OneSixth = 1.0d0/6.0d0
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, Getnpmyid
   integer(4) :: ia, ialow, iaupp, ja, jalow, jaupp
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc
   real(8)    :: r1, r2, r1i, r2i, r3i, r5i, r7i, threer5i, ex
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
   real(8)    :: doti, dotj, dotir5i, dotjr5i, dotir7i, dotjr7i
   real(8)    :: s, v, d2, d1, d0, d3
   real(8), external :: ErfLocal

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

! ... initiate

   etot(1:3,1:na) = Zero
   efg(1:6,1:na)  = Zero
   virelec        = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2

               if (lewald) then
                  ex = exp(-ualpha2*r2)
                  r1i = r1i*(One-ErfLocal(ualpha*r1))
                  r3i = r2i*(r1i + ewaldfac1*ex)
                  r5i = r2i*(r3i + ewaldfac2*ex)
                  r7i = r2i*(r5i + ewaldfac3*ex)
               else if (lrf) then
                  r3i = r2i*r1i
                  r5i = r2i*r3i
                  r7i = r2i*r5i
                  r3i = r3i - rffac           ! not to be used for field gradient
               else
                  r3i = r2i*r1i
                  r5i = r2i*r3i
                  r7i = r2i*r5i
               end if

               if (ldamping) then
                  s  = 1.4*1.662*(poltens(1,ia)*poltens(1,ja))**(OneSixth)
                  v  = min(One,r1/s)
                  d2 = v**4
                  d1 = 4*v**3 - 3*d2
                  d0 = d2 - 2*v**3 + 2*v
                  d3 = 0.2d0*d2
                  r3i = d1*r3i
                  r5i = d2*r5i
                  r7i = d3*r7i
               end if
               threer5i = Three*r5i

               doti = diptot(1,ia)*dx+diptot(2,ia)*dy+diptot(3,ia)*dz
               dotj = diptot(1,ja)*dx+diptot(2,ja)*dy+diptot(3,ja)*dz
               dotir5i = three*r5i*doti
               dotjr5i = three*r5i*dotj
               dotir7i = (Three*Five*r7i)*doti
               dotjr7i = (Three*Five*r7i)*dotj

! ... calculate field (etot)

!                   -3                         2    -5
!     e (r) = q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab
!
!                   -3                         2    -5
!     e (r) =-q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab

! ... calculate virial

!     vir = -q e  r
!               a  a

               fldx = (+az(ja)*dx - diptot(1,ja))*r3i + dotjr5i*dx
               fldy = (+az(ja)*dy - diptot(2,ja))*r3i + dotjr5i*dy
               fldz = (+az(ja)*dz - diptot(3,ja))*r3i + dotjr5i*dz
               etot(1,ia) = etot(1,ia) + fldx
               etot(2,ia) = etot(2,ia) + fldy
               etot(3,ia) = etot(3,ia) + fldz
               virelec = virelec - az(ia)*(fldx*dx + fldy*dy + fldz*dz)
               fldx = (-az(ia)*dx - diptot(1,ia))*r3i + dotir5i*dx
               fldy = (-az(ia)*dy - diptot(2,ia))*r3i + dotir5i*dy
               fldz = (-az(ia)*dz - diptot(3,ia))*r3i + dotir5i*dz
               etot(1,ja) = etot(1,ja) + fldx
               etot(2,ja) = etot(2,ja) + fldy
               etot(3,ja) = etot(3,ja) + fldz

! ... calculate field gradient (total)

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    -  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    +  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

! ... calculate virial

!     vir  = -u  e   r
!              a  ab  b

               efgxx = az(ja)*r3i + (-az(ja)*dx**2 + Two*diptot(1,ja)*dx + dotj)*threer5i - dotjr7i*dx**2
               efgyy = az(ja)*r3i + (-az(ja)*dy**2 + Two*diptot(2,ja)*dy + dotj)*threer5i - dotjr7i*dy**2
               efgzz = az(ja)*r3i + (-az(ja)*dz**2 + Two*diptot(3,ja)*dz + dotj)*threer5i - dotjr7i*dz**2
               efgxy =              (-az(ja)*dx*dy + diptot(1,ja)*dy + diptot(2,ja)*dx)*threer5i - dotjr7i*dx*dy
               efgxz =              (-az(ja)*dx*dz + diptot(1,ja)*dz + diptot(3,ja)*dx)*threer5i - dotjr7i*dx*dz
               efgyz =              (-az(ja)*dy*dz + diptot(2,ja)*dz + diptot(3,ja)*dy)*threer5i - dotjr7i*dy*dz
               efg(1,ia) = efg(1,ia) + efgxx
               efg(2,ia) = efg(2,ia) + efgyy
               efg(3,ia) = efg(3,ia) + efgzz
               efg(4,ia) = efg(4,ia) + efgxy
               efg(5,ia) = efg(5,ia) + efgxz
               efg(6,ia) = efg(6,ia) + efgyz
               virelec  = virelec -(diptot(1,ia)*efgxx*dx + diptot(1,ia)*efgxy*dy + diptot(1,ia)*efgxz*dz &
                                  + diptot(2,ia)*efgxy*dx + diptot(2,ia)*efgyy*dy + diptot(2,ia)*efgyz*dz &
                                  + diptot(3,ia)*efgxz*dx + diptot(3,ia)*efgyz*dy + diptot(3,ia)*efgzz*dz)
               efgxx = az(ia)*r3i + (-az(ia)*dx**2 - Two*diptot(1,ia)*dx - doti)*threer5i + dotir7i*dx**2
               efgyy = az(ia)*r3i + (-az(ia)*dy**2 - Two*diptot(2,ia)*dy - doti)*threer5i + dotir7i*dy**2
               efgzz = az(ia)*r3i + (-az(ia)*dz**2 - Two*diptot(3,ia)*dz - doti)*threer5i + dotir7i*dz**2
               efgxy =              (-az(ia)*dx*dy - diptot(1,ia)*dy - diptot(2,ia)*dx)*threer5i + dotir7i*dx*dy
               efgxz =              (-az(ia)*dx*dz - diptot(1,ia)*dz - diptot(3,ia)*dx)*threer5i + dotir7i*dx*dz
               efgyz =              (-az(ia)*dy*dz - diptot(2,ia)*dz - diptot(3,ia)*dy)*threer5i + dotir7i*dy*dz
               efg(1,ja) = efg(1,ja) + efgxx
               efg(2,ja) = efg(2,ja) + efgyy
               efg(3,ja) = efg(3,ja) + efgzz
               efg(4,ja) = efg(4,ja) + efgxy
               efg(5,ja) = efg(5,ja) + efgxz
               efg(6,ja) = efg(6,ja) + efgyz
            end do
         end do
      end do

! ... add contribution from center molecule

      if (lrf) then
         do ia = ialow, iaupp
            do ja = ialow, iaupp
               dx = r(1,ja)-r(1,ia)
               dy = r(2,ja)-r(2,ia)
               dz = r(3,ja)-r(3,ia)
               etot(1,ia) = etot(1,ia) + rffac*(az(ja)*dx + diptot(1,ja))
               etot(2,ia) = etot(2,ia) + rffac*(az(ja)*dy + diptot(2,ja))
               etot(3,ia) = etot(3,ia) + rffac*(az(ja)*dz + diptot(3,ja))
            end do
         end do
      end if

   end do

   if (itest == 3 .and. master) call TestFieldTot(uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine TestFieldTot(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine)//'  (real)', unit)
   write(unit,'(a)')  'atom                    etot'
   write(unit,'(i5,2x,3f12.7)') (ia, etot(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'virial', virelec
end subroutine TestFieldTot

!........................................................................

end subroutine FieldTot

!************************************************************************
!> \page energy energy.F90
!! **FieldTotEwald**
!! *calculate field, field gradient, and virial from charges and total dipoles; k-space*
!************************************************************************

!     this routine assumes that the eikx, eiky, and eikz arrays are generated

subroutine FieldTotEwald

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='FieldTotEwald'

   if (ltime) call CpuAdd('start', txroutine, 3, uout)

   if (txewaldrec == 'std') then
       call FieldTotEwaldRecStd
   else if (txewaldrec == 'spm') then
       call FieldTotEwaldRecSPM
   end if
   if (itest == 3 .and. master) call TestFieldTotEwald('(real + rec)', uout)
   call FieldTotEwaldSelf
   if (itest == 3 .and. master) call TestFieldTotEwald('(real + rec + self)', uout )
   if (lsurf) call FieldTotEwaldSurf
   if (itest == 3 .and. master) call TestFieldTotEwald('(real + rec + self + sur)', uout)

   if (ltime) call CpuAdd('stop', txroutine, 3, uout)

contains

!........................................................................

subroutine FieldTotEwaldRecStd

   character(40), parameter :: txroutine ='FieldTotEwaldRecStd'
   integer(4) :: nx, ny, nz, kn, ia
   real(8)    :: kfac2, facmm, facmp, facpm, facpp
   real(8)    :: termx, termy, termz
   real(8)    :: kx, ky, kz, kxkx, kyky, kzkz, kxky, kxkz, kykz
   complex(8) :: cmm, cmp, cpm, cpp

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

#if defined (_PAR_)

   sumeikr = cmplx(Zero,Zero)
   sumeikrd = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            kx = TwoPiBoxi(1)*nx
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*diptot(1,ia)
               termy = ky*diptot(2,ia)
               termz = kz*diptot(3,ia)
               facmm = -termx - termy + termz
               facmp = -termx + termy + termz
               facpm = +termx - termy + termz
               facpp = +termx + termy + termz
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), facmm) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), facmp) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), facpm) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), facpp) * eikr(ia,4)
               sumeikrd(kn,1) = sumeikrd(kn,1) + cmplx(Zero, facmm) * eikr(ia,1)
               sumeikrd(kn,2) = sumeikrd(kn,2) + cmplx(Zero, facmp) * eikr(ia,2)
               sumeikrd(kn,3) = sumeikrd(kn,3) + cmplx(Zero, facpm) * eikr(ia,3)
               sumeikrd(kn,4) = sumeikrd(kn,4) + cmplx(Zero, facpp) * eikr(ia,4)
            end do
         end do
      end do
   end do

   call par_allreduce_comps(sumeikr, eikraux, nkvec_word)
   call par_allreduce_comps(sumeikrd, eikraux, nkvec_word)

   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =      eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            kx = TwoPiBoxi(1)*nx
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
            end do

            kfac2 = Two*kfac(kn)
            kxkx = kfac2*kx**2
            kyky = kfac2*ky**2
            kzkz = kfac2*kz**2
            kxky = kfac2*kx*ky
            kxkz = kfac2*kx*kz
            kykz = kfac2*ky*kz
            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               etot(1,ia) = etot(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               etot(2,ia) = etot(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               etot(3,ia) = etot(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))

               efg(1,ia) = efg(1,ia) + kxkx*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(2,ia) = efg(2,ia) + kyky*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(3,ia) = efg(3,ia) + kzkz*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(4,ia) = efg(4,ia) + kxky*(+real(cmm) - real(cmp) - real(cpm) + real(cpp))
               efg(5,ia) = efg(5,ia) + kxkz*(-real(cmm) - real(cmp) + real(cpm) + real(cpp))
               efg(6,ia) = efg(6,ia) + kykz*(-real(cmm) + real(cmp) - real(cpm) + real(cpp))
            end do

            if (master) then
               termx = sum( real(sumeikrd(kn,1:4)*real(sumeikr(kn,1:4))) + aimag(sumeikrd(kn,1:4)*aimag(sumeikr(kn,1:4))) )
               termy = sum( real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2 )
               termz = One-(kx**2+ky**2+kz**2)/(Two*ualpha2)
               virelec = virelec - kfac(kn)*(Two*termx + termy*termz)
            end if

         end do
      end do
   end do

#else

   sumeikr  = cmplx(Zero,Zero)
   sumeikrd = cmplx(Zero,Zero)
   kn = 0
   do nz = 0, ncut
      kz = TwoPiBoxi(3)*nz
      do ny = 0, ncut
         ky = TwoPiBoxi(2)*ny
         if (ny**2+nz**2 > ncut2) cycle
         do ia = iamyid(1), iamyid(2)
            eikyzm(ia) = conjg(eiky(ia,ny))*eikz(ia,nz)
            eikyzp(ia) =       eiky(ia,ny) *eikz(ia,nz)
         end do
         do nx = 0, ncut
            kx = TwoPiBoxi(1)*nx
            if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle      ! only even nx+ny+nz for RD and TO bc
            if (nx**2+ny**2+nz**2 > ncut2) cycle
            if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
            kn = kn + 1
            do ia = iamyid(1), iamyid(2)
               eikr(ia,1) = conjg(eikx(ia,nx))*eikyzm(ia)
               eikr(ia,2) = conjg(eikx(ia,nx))*eikyzp(ia)
               eikr(ia,3) =       eikx(ia,nx) *eikyzm(ia)
               eikr(ia,4) =       eikx(ia,nx) *eikyzp(ia)
               termx = kx*diptot(1,ia)
               termy = ky*diptot(2,ia)
               termz = kz*diptot(3,ia)
               facmm = -termx - termy + termz
               facmp = -termx + termy + termz
               facpm = +termx - termy + termz
               facpp = +termx + termy + termz
               sumeikr(kn,1) = sumeikr(kn,1) + cmplx(az(ia), facmm) * eikr(ia,1)
               sumeikr(kn,2) = sumeikr(kn,2) + cmplx(az(ia), facmp) * eikr(ia,2)
               sumeikr(kn,3) = sumeikr(kn,3) + cmplx(az(ia), facpm) * eikr(ia,3)
               sumeikr(kn,4) = sumeikr(kn,4) + cmplx(az(ia), facpp) * eikr(ia,4)
               sumeikrd(kn,1) = sumeikrd(kn,1) + cmplx(Zero, facmm) * eikr(ia,1)
               sumeikrd(kn,2) = sumeikrd(kn,2) + cmplx(Zero, facmp) * eikr(ia,2)
               sumeikrd(kn,3) = sumeikrd(kn,3) + cmplx(Zero, facpm) * eikr(ia,3)
               sumeikrd(kn,4) = sumeikrd(kn,4) + cmplx(Zero, facpp) * eikr(ia,4)
            end do

            kfac2 = Two*kfac(kn)
            kxkx = kfac2*kx**2
            kyky = kfac2*ky**2
            kzkz = kfac2*kz**2
            kxky = kfac2*kx*ky
            kxkz = kfac2*kx*kz
            kykz = kfac2*ky*kz
            do ia = iamyid(1), iamyid(2)
               cmm = eikr(ia,1)*conjg(sumeikr(kn,1))
               cmp = eikr(ia,2)*conjg(sumeikr(kn,2))
               cpm = eikr(ia,3)*conjg(sumeikr(kn,3))
               cpp = eikr(ia,4)*conjg(sumeikr(kn,4))

               etot(1,ia) = etot(1,ia) + kfac2*kx*(-aimag(cmm) - aimag(cmp) + aimag(cpm) + aimag(cpp))
               etot(2,ia) = etot(2,ia) + kfac2*ky*(-aimag(cmm) + aimag(cmp) - aimag(cpm) + aimag(cpp))
               etot(3,ia) = etot(3,ia) + kfac2*kz*(+aimag(cmm) + aimag(cmp) + aimag(cpm) + aimag(cpp))

               efg(1,ia) = efg(1,ia) + kxkx*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(2,ia) = efg(2,ia) + kyky*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(3,ia) = efg(3,ia) + kzkz*(+real(cmm) + real(cmp) + real(cpm) + real(cpp))
               efg(4,ia) = efg(4,ia) + kxky*(+real(cmm) - real(cmp) - real(cpm) + real(cpp))
               efg(5,ia) = efg(5,ia) + kxkz*(-real(cmm) - real(cmp) + real(cpm) + real(cpp))
               efg(6,ia) = efg(6,ia) + kykz*(-real(cmm) + real(cmp) - real(cpm) + real(cpp))
            end do

            termx = sum( real(sumeikrd(kn,1:4)*real(sumeikr(kn,1:4))) + aimag(sumeikrd(kn,1:4)*aimag(sumeikr(kn,1:4))) )
            termy = sum( real(sumeikr(kn,1:4))**2 + aimag(sumeikr(kn,1:4))**2 )
            termz = One-(kx**2+ky**2+kz**2)/(Two*ualpha2)
            virelec = virelec - kfac(kn)*(Two*termx + termy*termz)

         end do
      end do
   end do

#endif

   if (ltime) call CpuAdd('stop', txroutine, 4, uout)

end subroutine FieldTotEwaldRecStd

!........................................................................

subroutine FieldTotEwaldRecSPM

# ifdef F03_CBIND

   use EnergyModule, s=>meshsize
   implicit none

   character(40), parameter :: txroutine ='FieldTotEwaldRecSPM'
   integer(4) :: ia, m, ix, iy, iz, nx, ny, nz
   real(8)    :: q, q1, dipx, dipx1, dipx2, dipy, dipy1, dipy2, dipz, dipz1, dipz2
   real(8)    :: splx, sply, splz
   real(8)    :: psum, esum(3), efgsum(6)
   real(8)    :: dsdx, dsdy, dsdz, d2sdx, d2sdy, d2sdz, Qval, Qdip

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

! ... make generalized charge mesh

   if (ltime) call CpuAdd('start', 'TotMakeGQMesh', 5, uout)
   Qmesh = Zero
   do ia = 1, na
      dipx2 = diptot(1,ia)*dkdr(1)
      dipy2 = diptot(2,ia)*dkdr(2)
      dipz2 = diptot(3,ia)*dkdr(3)
      do m = 1,3                    ! Cardinal B-spline, used for calculating exp(i*k(1:3)*r(1:3))
         call CardinalBSpline(order, dkdr(m)*r(m,ia), meshmax(m,ia), spline(:,m,ia), deriv(:,m,ia), deriv2(:,m,ia))
         meshmax(m,ia) = mod(meshmax(m,ia)+s(m),s(m))
      end do
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz = spline(iz,3,ia)
         dsdz = deriv(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply = spline(iy,2,ia)
            q1 = az(ia)*sply*splz
            dipx1 = dipx2*sply*splz
            dipy1 = dipy2*deriv(iy,2,ia)*splz
            dipz1 = dipz2*sply*dsdz
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx = spline(ix,1,ia)
               q = q1*splx
               dipx = dipx1*deriv(ix,1,ia)
               dipy = dipy1*splx
               dipz = dipz1*splx
               QMesh(nx+1,ny+1,nz+1) = QMesh(nx+1,ny+1,nz+1) + q + dipx + dipy + dipz
            end do
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', 'TotMakeGQMesh', 5, uout)


! ... make Fourier transformation, reciprocal space operations, and back FFT

   call SPMFFTRec(.false., .false., .true., 'TotMakeFFT', 'TotCalcGIF', 5, raux, virelec)

! ... calculate potential, field, field gradient, and virial

   if (ltime) call CpuAdd('start', 'TotCalcProp', 5, uout)
   do ia=1,na
      psum = Zero
      esum = Zero
      efgsum = Zero
      dipx2 = diptot(1,ia)*dkdr(1)
      dipy2 = diptot(2,ia)*dkdr(2)
      dipz2 = diptot(3,ia)*dkdr(3)
      do iz = 0,order-1
         nz = mod(s(3)+meshmax(3,ia)-iz, s(3))
         splz =  spline(iz,3,ia)
         dsdz =  deriv(iz,3,ia)
         d2sdz = deriv2(iz,3,ia)
         do iy = 0,order-1
            ny = mod(s(2)+meshmax(2,ia)-iy, s(2))
            sply =  spline(iy,2,ia)
            dsdy =   deriv(iy,2,ia)
            d2sdy = deriv2(iy,2,ia)
            do ix = 0,order-1
               nx = mod(s(1)+meshmax(1,ia)-ix, s(1))
               splx =  spline(ix,1,ia)
               dsdx =   deriv(ix,1,ia)
               d2sdx = deriv2(ix,1,ia)
               Qval = QMesh(nx+1,ny+1,nz+1)

               esum(1) = esum(1) + dsdx*sply*splz*Qval
               esum(2) = esum(2) + splx*dsdy*splz*Qval
               esum(3) = esum(3) + splx*sply*dsdz*Qval

               efgsum(1) = efgsum(1) + d2sdx*sply*splz*Qval
               efgsum(2) = efgsum(2) + splx*d2sdy*splz*Qval
               efgsum(3) = efgsum(3) + splx*sply*d2sdz*Qval
               efgsum(4) = efgsum(4) + dsdx*dsdy*splz*Qval
               efgsum(5) = efgsum(5) + dsdx*sply*dsdz*Qval
               efgsum(6) = efgsum(6) + splx*dsdy*dsdz*Qval

               Qdip = dipx2*dsdx*sply*splz + dipy2*splx*dsdy*splz + dipz2*splx*sply*dsdz
               virelec = virelec - Qdip*Qval
            end do
         end do
      end do
      etot(1:3,ia) = etot(1:3,ia) - esum(1:3)*dkdr(1:3)               ! E = -dU/dr = -dU/dk * dk/dr
      efg(1:6,ia) = efg(1:6,ia) - efgsum(1:6)*d2kdr(1:6)
   end do
   if (ltime) call CpuAdd('start', 'TotCalcProp', 5, uout)

   if (ltime) call CpuAdd('start', txroutine, 4, uout)

# endif

end subroutine FieldTotEwaldRecSPM

!........................................................................

subroutine FieldTotEwaldSelf
   real(8), external :: ErfLocal
   integer(4) :: ip, ipt, ia, ialow, iaupp, ja
   real(8)    :: dx, dy, dz, ex, r1, r2, r1i, r2i, r3i, r5i, r7i, doti, dotj, dotir5i, dotjr5i, dotir7i, dotjr7i, threer5i
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz

! ... atomic contribution

   etot(1,iamyid(1):iamyid(2)) = etot(1,iamyid(1):iamyid(2)) + ewaldselffac2*diptot(1,iamyid(1):iamyid(2))
   etot(2,iamyid(1):iamyid(2)) = etot(2,iamyid(1):iamyid(2)) + ewaldselffac2*diptot(2,iamyid(1):iamyid(2))
   etot(3,iamyid(1):iamyid(2)) = etot(3,iamyid(1):iamyid(2)) + ewaldselffac2*diptot(3,iamyid(1):iamyid(2))
   efg(1,iamyid(1):iamyid(2)) = efg(1,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))
   efg(2,iamyid(1):iamyid(2)) = efg(2,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))
   efg(3,iamyid(1):iamyid(2)) = efg(3,iamyid(1):iamyid(2)) - ewaldselffac2*az(iamyid(1):iamyid(2))

! ... molecular contribution (_not_ constant)

   do ip = ipmyid(1), ipmyid(2)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         do ja = ia+1, iaupp
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2

            r1 = sqrt(r2)
            r1i = One/r1
            r2i = r1i**2
            ex = exp(-ualpha2*r2)
            r1i = r1i*ErfLocal(ualpha*r1)
            r3i = r2i*(r1i - ewaldfac1*ex)
            r5i = r2i*(r3i - ewaldfac2*ex)
            r7i = r2i*(r5i - ewaldfac3*ex)
            threer5i = Three*r5i

            doti = diptot(1,ia)*dx+diptot(2,ia)*dy+diptot(3,ia)*dz
            dotj = diptot(1,ja)*dx+diptot(2,ja)*dy+diptot(3,ja)*dz
            dotir5i = threer5i*doti
            dotjr5i = threer5i*dotj
            dotir7i = (Three*Five*r7i)*doti
            dotjr7i = (Three*Five*r7i)*dotj

            fldx = -((az(ja)*dx - diptot(1,ja))*r3i + dotjr5i*dx)
            fldy = -((az(ja)*dy - diptot(2,ja))*r3i + dotjr5i*dy)
            fldz = -((az(ja)*dz - diptot(3,ja))*r3i + dotjr5i*dz)
            etot(1,ia) = etot(1,ia) + fldx
            etot(2,ia) = etot(2,ia) + fldy
            etot(3,ia) = etot(3,ia) + fldz
            virelec = virelec - (az(ia)*(fldx*dx + fldy*dy + fldz*dz))
            fldx = -((-az(ia)*dx - diptot(1,ia))*r3i + dotir5i*dx)
            fldy = -((-az(ia)*dy - diptot(2,ia))*r3i + dotir5i*dy)
            fldz = -((-az(ia)*dz - diptot(3,ia))*r3i + dotir5i*dz)
            etot(1,ja) = etot(1,ja) + fldx
            etot(2,ja) = etot(2,ja) + fldy
            etot(3,ja) = etot(3,ja) + fldz

            efgxx = -( az(ja)*r3i + (-az(ja)*dx**2 + Two*diptot(1,ja)*dx + dotj)*threer5i - dotjr7i*dx**2 )
            efgyy = -( az(ja)*r3i + (-az(ja)*dy**2 + Two*diptot(2,ja)*dy + dotj)*threer5i - dotjr7i*dy**2 )
            efgzz = -( az(ja)*r3i + (-az(ja)*dz**2 + Two*diptot(3,ja)*dz + dotj)*threer5i - dotjr7i*dz**2 )
            efgxy = -( (-az(ja)*dx*dy + diptot(1,ja)*dy + diptot(2,ja)*dx)*threer5i - dotjr7i*dx*dy )
            efgxz = -( (-az(ja)*dx*dz + diptot(1,ja)*dz + diptot(3,ja)*dx)*threer5i - dotjr7i*dx*dz )
            efgyz = -( (-az(ja)*dy*dz + diptot(2,ja)*dz + diptot(3,ja)*dy)*threer5i - dotjr7i*dy*dz )
            efg(1,ia) = efg(1,ia) + efgxx
            efg(2,ia) = efg(2,ia) + efgyy
            efg(3,ia) = efg(3,ia) + efgzz
            efg(4,ia) = efg(4,ia) + efgxy
            efg(5,ia) = efg(5,ia) + efgxz
            efg(6,ia) = efg(6,ia) + efgyz
            virelec = virelec - ( diptot(1,ia)*efgxx*dx + diptot(1,ia)*efgxy*dy + diptot(1,ia)*efgxz*dz   &
                                + diptot(2,ia)*efgxy*dx + diptot(2,ia)*efgyy*dy + diptot(2,ia)*efgyz*dz   &
                                + diptot(3,ia)*efgxz*dx + diptot(3,ia)*efgyz*dy + diptot(3,ia)*efgzz*dz )
            efgxx = -( az(ia)*r3i + (-az(ia)*dx**2 - Two*diptot(1,ia)*dx - doti)*threer5i + dotir7i*dx**2 )
            efgyy = -( az(ia)*r3i + (-az(ia)*dy**2 - Two*diptot(2,ia)*dy - doti)*threer5i + dotir7i*dy**2 )
            efgzz = -( az(ia)*r3i + (-az(ia)*dz**2 - Two*diptot(3,ia)*dz - doti)*threer5i + dotir7i*dz**2 )
            efgxy = -( (-az(ia)*dx*dy - diptot(1,ia)*dy - diptot(2,ia)*dx)*threer5i + dotir7i*dx*dy )
            efgxz = -( (-az(ia)*dx*dz - diptot(1,ia)*dz - diptot(3,ia)*dx)*threer5i + dotir7i*dx*dz )
            efgyz = -( (-az(ia)*dy*dz - diptot(2,ia)*dz - diptot(3,ia)*dy)*threer5i + dotir7i*dy*dz )
            efg(1,ja) = efg(1,ja) + efgxx
            efg(2,ja) = efg(2,ja) + efgyy
            efg(3,ja) = efg(3,ja) + efgzz
            efg(4,ja) = efg(4,ja) + efgxy
            efg(5,ja) = efg(5,ja) + efgxz
            efg(6,ja) = efg(6,ja) + efgyz
         end do
      end do
   end do
end subroutine FieldTotEwaldSelf

!........................................................................

subroutine FieldTotEwaldSurf
   real(8) :: fac, sumqrx, sumqry, sumqrz, sumdx, sumdy, sumdz
   fac = FourPiThird/vol
   sumqrx = sum(az(1:na)*r(1,1:na))
   sumqry = sum(az(1:na)*r(2,1:na))
   sumqrz = sum(az(1:na)*r(3,1:na))
   sumdx = sum(diptot(1,1:na))
   sumdy = sum(diptot(2,1:na))
   sumdz = sum(diptot(3,1:na))
   etot(1,iamyid(1):iamyid(2)) = etot(1,iamyid(1):iamyid(2)) - fac*(sumqrx + sumdx)
   etot(2,iamyid(1):iamyid(2)) = etot(2,iamyid(1):iamyid(2)) - fac*(sumqry + sumdy)
   etot(3,iamyid(1):iamyid(2)) = etot(3,iamyid(1):iamyid(2)) - fac*(sumqrz + sumdz)

   if (master) virelec  = virelec - TwoPi/(Three*vol) *                              &
                        (        (sumqrx**2 + sumqry**2 + sumqrz**2)                &
                         +   Two*(Two*(sumqrx*sumdx + sumqry*sumdy + sumqrz*sumdz)) &
                         + Three*(sumdx**2 + sumdy**2 + sumdz**2 )    )
end subroutine FieldTotEwaldSurf

!........................................................................

subroutine TestFieldTotEwald(label, unit)
   character(*), intent(in) :: label
   integer(4), intent(in) :: unit
   integer(4) :: ia
   call WriteHead(3, 'Test'//trim(txroutine)//'  '//label, unit)
   write(unit,'(a)')  'atom                    etot'
   write(unit,'(i5,2x,3f12.7)') (ia, etot(1:3,ia), ia=1,min(6,int(na)))
   write(unit,'(a)')  'atom                                      efg'
   write(unit,'(i5,2x,6f12.7)') (ia, efg(1:6,ia), ia=1,min(6,int(na)))
   write(unit,'(a,2x,f12.7)') 'virelec', virelec
end subroutine TestFieldTotEwald

!........................................................................

end subroutine FieldTotEwald

!************************************************************************
!> \page energy energy.F90
!! **UIntraReac**
!! *calculate intramolecular contribution to energy and force in rf*
!************************************************************************


subroutine UIntraReac

   use EnergyModule
   implicit none

   integer(4) ::  ip, iploc, ipt, ia, ialow, iaupp, ja, Getnpmyid
   real(8)    ::  dx, dy, dz, r2, uintra, fintrax, fintray, fintraz

   uintra = Zero
   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ianpn(ip)+napt(ipt)-1
      do ia = ialow, iaupp
         fintrax = Zero
         fintray = Zero
         fintraz = Zero
         do ja = ialow, iaupp
            if (ja == ia) cycle
            dx = r(1,ia)-r(1,ja)
            dy = r(2,ia)-r(2,ja)
            dz = r(3,ia)-r(3,ja)
            r2 = dx**2 + dy**2 + dz**2
            uintra = uintra + az(ia)*az(ja)*rffac2*r2
            fintrax = fintrax - az(ja)*dx
            fintray = fintray - az(ja)*dy
            fintraz = fintraz - az(ja)*dz
         end do
         force(1,ia) = force(1,ia) + EpsiFourPi*az(ia)*rffac*fintrax
         force(2,ia) = force(2,ia) + EpsiFourPi*az(ia)*rffac*fintray
         force(3,ia) = force(3,ia) + EpsiFourPi*az(ia)*rffac*fintraz
      end do
   end do

   u%tot = u%tot + Half*EpsiFourPi*uintra

end subroutine UIntraReac

#if defined (_PAR_)
!************************************************************************
!> \page energy energy.F90
!! **PackReduceU**
!! *pack the different energy contributions and perform all_reduce*
!************************************************************************


subroutine PackReduceU(n1, n2, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, uaux)

   implicit none

   integer(4), intent(in)    :: n1, n2
   real(8)   , intent(inout) :: e1, e2(n1), e3(n2), e4, e5, e6, e7, e8, e9, e10, uaux(*)
   integer(4) :: nred        ! length of packed array

   nred = n1+n2+8
   uaux(1)       = e1
   uaux(2:n1+1)  = e2(1:n1)
   uaux(n1+2:n1+n2+1) = e3(1:n2)
   uaux(n1+n2+2) = e4
   uaux(n1+n2+3) = e5
   uaux(n1+n2+4) = e6
   uaux(n1+n2+5) = e7
   uaux(n1+n2+6) = e8
   uaux(n1+n2+7) = e9
   uaux(n1+n2+8) = e10
   call par_allreduce_reals(uaux,uaux(nred+1),nred)
   e1       = uaux(1)
   e2(1:n1) = uaux(2:n1+1)
   e3(1:n2) = uaux(n1+2:n1+n2+1)
   e4       = uaux(n1+n2+2)
   e5       = uaux(n1+n2+3)
   e6       = uaux(n1+n2+4)
   e7       = uaux(n1+n2+5)
   e8       = uaux(n1+n2+6)
   e9       = uaux(n1+n2+7)
   e10      = uaux(n1+n2+8)

end subroutine PackReduceU
#endif

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UDipoleSph**
!! *calculate charge and dipole potential energy using image charge approximation*
!************************************************************************


subroutine UDipoleSph

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDipoleSph'
   integer(4) :: ip, ipt, jp, jpt, iptjpt, ibuf
   real(8)    :: dx, dy, dz
   real(8)    :: r1, r2, r1i, r2i, r3i, r5i, r7i, threer5i, d, usum, fsum, virtwob
   real(8)    :: forcex, forcey, forcez, torquex, torquey, torquez
   real(8)    :: doti, dotj, dotir5i, dotjr5i, dotir7i, dotjr7i
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
   real(8), external :: ErfLocal

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   u%stat           = Zero
   virtwob          = Zero
   potstat(1:na)    = Zero
   estat(1:3,1:na)  = Zero
   efg(1:6,1:na)    = Zero
   virelec          = Zero

   do ip = 1, np
      ipt = iptpn(ip)
      do jp = ipmyid(1), ipmyid(2)
         if (jp == ip) cycle
         if (lmc) then
           if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         r2 = dx**2+dy**2+dz**2
         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         if (r2 < r2umin(iptjpt)) call StopUDipoleSph
         ibuf = iubuflow(iptjpt)
         do
           if (r2 >= ubuf(ibuf)) exit
           ibuf = ibuf+12
           if (ibuf > nbuf) call StopIbuf('txptpt',iptjpt)
         end do
         d = r2-ubuf(ibuf)
         usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
         fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))

         force(1,ip) = force(1,ip) + (fsum * dx)
         force(2,ip) = force(2,ip) + (fsum * dy)
         force(3,ip) = force(3,ip) + (fsum * dz)
         force(1,jp) = force(1,jp) - (fsum * dx)
         force(2,jp) = force(2,jp) - (fsum * dy)
         force(3,jp) = force(3,jp) - (fsum * dz)
         virtwob     = virtwob     - (fsum * r2)

         u%twob(iptjpt) = u%twob(iptjpt)+ usum

         r3i = r2i*r1i
         r5i = r2i*r3i
         r7i = r2i*r5i
         threer5i = Three*r5i

         doti = dip(1,ip)*dx+dip(2,ip)*dy+dip(3,ip)*dz
         dotj = dip(1,jp)*dx+dip(2,jp)*dy+dip(3,jp)*dz
         dotir5i = three*r5i*doti
         dotjr5i = three*r5i*dotj
         dotir7i = (Three*Five*r7i)*doti
         dotjr7i = (Three*Five*r7i)*dotj

! ... calculate potential (potstat)

!                                  3
!     pot(r) = q / r  +  u   r  / r
!                          a   a
!                                  3
!     pot(0) = q / r  -  u   r  / r
!                          a   a

         potstat(ip) = potstat(ip) + (az(jp)*r1i + dotj*r3i)
         potstat(jp) = potstat(jp) + (az(ip)*r1i - doti*r3i)

! ... calculate field (estat)

!                   -3                         2    -5
!     e (r) = q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab
!
!                   -3                         2    -5
!     e (r) =-q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab

! ... calculate virial

!     vir = -q e  r
!               a  a

         fldx = (+az(jp)*dx - dip(1,jp))*r3i + dotjr5i*dx
         fldy = (+az(jp)*dy - dip(2,jp))*r3i + dotjr5i*dy
         fldz = (+az(jp)*dz - dip(3,jp))*r3i + dotjr5i*dz
         estat(1,ip) = estat(1,ip) + fldx
         estat(2,ip) = estat(2,ip) + fldy
         estat(3,ip) = estat(3,ip) + fldz
         virelec = virelec - az(ip)*(fldx*dx + fldy*dy + fldz*dz)
         fldx = (-az(ip)*dx - dip(1,ip))*r3i + dotir5i*dx
         fldy = (-az(ip)*dy - dip(2,ip))*r3i + dotir5i*dy
         fldz = (-az(ip)*dz - dip(3,ip))*r3i + dotir5i*dz
         estat(1,jp) = estat(1,jp) + fldx
         estat(2,jp) = estat(2,jp) + fldy
         estat(3,jp) = estat(3,jp) + fldz

! ... calculate field gradient (stat)

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    -  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    +  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

! ... calculate virial

!     vir  = -u  e   r
!              a  ab  b

         efgxx = az(jp)*r3i + (-az(jp)*dx**2 + Two*dip(1,jp)*dx + dotj)*threer5i - dotjr7i*dx**2
         efgyy = az(jp)*r3i + (-az(jp)*dy**2 + Two*dip(2,jp)*dy + dotj)*threer5i - dotjr7i*dy**2
         efgzz = az(jp)*r3i + (-az(jp)*dz**2 + Two*dip(3,jp)*dz + dotj)*threer5i - dotjr7i*dz**2
         efgxy =              (-az(jp)*dx*dy + dip(1,jp)*dy + dip(2,jp)*dx)*threer5i - dotjr7i*dx*dy
         efgxz =              (-az(jp)*dx*dz + dip(1,jp)*dz + dip(3,jp)*dx)*threer5i - dotjr7i*dx*dz
         efgyz =              (-az(jp)*dy*dz + dip(2,jp)*dz + dip(3,jp)*dy)*threer5i - dotjr7i*dy*dz
         efg(1,ip) = efg(1,ip) + efgxx
         efg(2,ip) = efg(2,ip) + efgyy
         efg(3,ip) = efg(3,ip) + efgzz
         efg(4,ip) = efg(4,ip) + efgxy
         efg(5,ip) = efg(5,ip) + efgxz
         efg(6,ip) = efg(6,ip) + efgyz
         virelec  = virelec -(dip(1,ip)*efgxx*dx + dip(1,ip)*efgxy*dy + dip(1,ip)*efgxz*dz &
                            + dip(2,ip)*efgxy*dx + dip(2,ip)*efgyy*dy + dip(2,ip)*efgyz*dz &
                            + dip(3,ip)*efgxz*dx + dip(3,ip)*efgyz*dy + dip(3,ip)*efgzz*dz)
         efgxx = az(ip)*r3i + (-az(ip)*dx**2 - Two*dip(1,ip)*dx - doti)*threer5i + dotir7i*dx**2
         efgyy = az(ip)*r3i + (-az(ip)*dy**2 - Two*dip(2,ip)*dy - doti)*threer5i + dotir7i*dy**2
         efgzz = az(ip)*r3i + (-az(ip)*dz**2 - Two*dip(3,ip)*dz - doti)*threer5i + dotir7i*dz**2
         efgxy =              (-az(ip)*dx*dy - dip(1,ip)*dy - dip(2,ip)*dx)*threer5i + dotir7i*dx*dy
         efgxz =              (-az(ip)*dx*dz - dip(1,ip)*dz - dip(3,ip)*dx)*threer5i + dotir7i*dx*dz
         efgyz =              (-az(ip)*dy*dz - dip(2,ip)*dz - dip(3,ip)*dy)*threer5i + dotir7i*dy*dz
         efg(1,jp) = efg(1,jp) + efgxx
         efg(2,jp) = efg(2,jp) + efgyy
         efg(3,jp) = efg(3,jp) + efgzz
         efg(4,jp) = efg(4,jp) + efgxy
         efg(5,jp) = efg(5,jp) + efgxz
         efg(6,jp) = efg(6,jp) + efgyz
      end do
   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot     + u%twob(0)
   virial    = virial    + virtwob

   if (laimage) then

   do ip = ipmyid(1), ipmyid(2)     ! add contribution to pot, field, and field gradient from image charges and dipoles
      do jp = 1, np
         dx = ro(1,ip)-rimg(1,jp)
         dy = ro(2,ip)-rimg(2,jp)
         dz = ro(3,ip)-rimg(3,jp)
         r2 = dx**2+dy**2+dz**2
         r1 = sqrt(r2)
         r1i = One/r1
         r2i = r1i**2
         r3i = r2i*r1i
         r5i = r2i*r3i
         r7i = r2i*r5i
         threer5i = Three*r5i

         doti = dip(1,ip)*dx+dip(2,ip)*dy+dip(3,ip)*dz
         dotj = dipimg(1,jp)*dx+dipimg(2,jp)*dy+dipimg(3,jp)*dz
         dotir5i = three*r5i*doti
         dotjr5i = three*r5i*dotj
         dotir7i = (Three*Five*r7i)*doti
         dotjr7i = (Three*Five*r7i)*dotj

         potstat(ip) = potstat(ip) + (zimg(jp)*r1i + dotj*r3i)

         fldx = (+zimg(jp)*dx - dipimg(1,jp))*r3i + dotjr5i*dx
         fldy = (+zimg(jp)*dy - dipimg(2,jp))*r3i + dotjr5i*dy
         fldz = (+zimg(jp)*dz - dipimg(3,jp))*r3i + dotjr5i*dz
         estat(1,ip) = estat(1,ip) + fldx
         estat(2,ip) = estat(2,ip) + fldy
         estat(3,ip) = estat(3,ip) + fldz
         virelec = virelec - az(ip)*(fldx*dx + fldy*dy + fldz*dz)

         efgxx = zimg(jp)*r3i + (-zimg(jp)*dx**2 + Two*dipimg(1,jp)*dx + dotj)*threer5i - dotjr7i*dx**2
         efgyy = zimg(jp)*r3i + (-zimg(jp)*dy**2 + Two*dipimg(2,jp)*dy + dotj)*threer5i - dotjr7i*dy**2
         efgzz = zimg(jp)*r3i + (-zimg(jp)*dz**2 + Two*dipimg(3,jp)*dz + dotj)*threer5i - dotjr7i*dz**2
         efgxy =                (-zimg(jp)*dx*dy + dipimg(1,jp)*dy + dipimg(2,jp)*dx)*threer5i - dotjr7i*dx*dy
         efgxz =                (-zimg(jp)*dx*dz + dipimg(1,jp)*dz + dipimg(3,jp)*dx)*threer5i - dotjr7i*dx*dz
         efgyz =                (-zimg(jp)*dy*dz + dipimg(2,jp)*dz + dipimg(3,jp)*dy)*threer5i - dotjr7i*dy*dz
         efg(1,ip) = efg(1,ip) + efgxx
         efg(2,ip) = efg(2,ip) + efgyy
         efg(3,ip) = efg(3,ip) + efgzz
         efg(4,ip) = efg(4,ip) + efgxy
         efg(5,ip) = efg(5,ip) + efgxz
         efg(6,ip) = efg(6,ip) + efgyz
         virelec  = virelec -(dip(1,ip)*efgxx*dx + dip(1,ip)*efgxy*dy + dip(1,ip)*efgxz*dz &
                            + dip(2,ip)*efgxy*dx + dip(2,ip)*efgyy*dy + dip(2,ip)*efgyz*dz &
                            + dip(3,ip)*efgxz*dx + dip(3,ip)*efgyz*dy + dip(3,ip)*efgzz*dz)

         end do
     end do

     end if

#if defined (_PAR_)
   if (ltime) call CpuAdd('start', 'comm', 3, uout)
   call par_allreduce_reals(potstat, vaux, na  )
   call par_allreduce_reals(estat,   vaux, 3*na)
   call par_allreduce_reals(efg    , vaux, 6*na)
   if (ltime) call CpuAdd('stop', 'comm', 3, uout)
#endif

    do ip = ipmyid(1), ipmyid(2)

      u%stat = u%stat + az(ip)*potstat(ip) - (dip(1,ip)*estat(1,ip) + dip(2,ip)*estat(2,ip) + dip(3,ip)*estat(3,ip))

      forcex = az(ip)*estat(1,ip) + (dip(1,ip)*efg(1,ip) + dip(2,ip)*efg(4,ip) + dip(3,ip)*efg(5,ip))
      forcey = az(ip)*estat(2,ip) + (dip(1,ip)*efg(4,ip) + dip(2,ip)*efg(2,ip) + dip(3,ip)*efg(6,ip))
      forcez = az(ip)*estat(3,ip) + (dip(1,ip)*efg(5,ip) + dip(2,ip)*efg(6,ip) + dip(3,ip)*efg(3,ip))
      force(1,ip)  = force(1,ip)  + EpsiFourPi*forcex
      force(2,ip)  = force(2,ip)  + EpsiFourPi*forcey
      force(3,ip)  = force(3,ip)  + EpsiFourPi*forcez

      torquex = dip(2,ip)*estat(3,ip) - dip(3,ip)*estat(2,ip)
      torquey = dip(3,ip)*estat(1,ip) - dip(1,ip)*estat(3,ip)
      torquez = dip(1,ip)*estat(2,ip) - dip(2,ip)*estat(1,ip)
      torque(1,ip) = torque(1,ip) + EpsiFourPi*torquex
      torque(2,ip) = torque(2,ip) + EpsiFourPi*torquey
      torque(3,ip) = torque(3,ip) + EpsiFourPi*torquez

   end do

   u%stat = EpsiFourPi*Half*u%stat

! ... update

   u%tot  = u%tot  + u%stat
   virial = virial + EpsiFourPi*virelec

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUDipoleSph
   character(40), parameter :: txroutine ='StopUDipoleSph'
   write(uout,'(a,i15)')   'uout          = ', uout
   write(uout,'(a,i15)')   'myid          = ', myid
   write(*,'(a,i15)')   'myid          = ', myid
   write(uout,'(a,i15)')   'iptjpt        = ', iptjpt
   write(uout,'(a,e15.5)') 'r2            = ', r2
   write(uout,'(a,e15.5)') 'r2umin(iptjpt) = ', r2umin(iptjpt)
   call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)
end subroutine StopUDipoleSph

!........................................................................

end subroutine UDipoleSph

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UDielDis**
!! *calculate energy for charge and in a system with dielectric discontinuities*
!************************************************************************


subroutine UDielDis

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDielDis'

   if (.not.lmc) call Stop(txroutine, '.not.lmc', uout)   ! MC is required
   if (txbc == 'xy') Call UDielDisPlane
   if (txbc == 'sph') Call UDielDisSph

end subroutine UDielDis

!************************************************************************
!> \page energy energy.F90
!! **UDielDisPlane**
!! *calculate coulomb energy in a system with dielectric discontinuity at z = 0*
!************************************************************************


subroutine UDielDisPlane

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDielDisPlane'
   logical, save :: first = .true.
   integer(4) :: ipt, jpt, ip, jp, iptjpt, iploc, jploc, Getnpmyid
   real(8) :: dx, dy, dz, r2, ri, rip

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... setup

   if (first) then
      first = .false.
      elunit = Epsi0FourPi
      if (itest == 98) elunit = one  ! set electrostatic unit; =1 for test purpose; => value of Green's function)
      epsi1FourPi = elunit/epsi1
      epsi2FourPi = elunit/epsi2
      eta = epsi1/epsi2
      delta = (one-eta)/(one+eta)
   end if

! ... calculate energy

   u%twob = Zero
   u%oneb = Zero

   do iploc = 1, Getnpmyid()         ! triangular loop (note U(q_1,q_2(image) = U(q_2,q1(image)) )
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (jp <= ip) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         ri = one/sqrt(r2)
         if ((r(3,ip) < Zero) .and. (r(3,jp) < Zero)) then  ! ion--ion and ion--image interaction
            dz = ro(3,ip)+ro(3,jp)       ! image location
            call PBCr2(dx,dy,dz,r2)
            rip = one/sqrt(r2)
            u%twob(iptjpt) = u%twob(iptjpt) + epsi1FourPi*az(ip)*az(jp)*(ri - delta*rip)
         elseif ((r(3,ip) > Zero) .and. (r(3,jp) > Zero)) then
            dz = ro(3,ip)+ro(3,jp)       ! image location
            call PBCr2(dx,dy,dz,r2)
            rip = one/sqrt(r2)
            u%twob(iptjpt) = u%twob(iptjpt) + epsi2FourPi*az(ip)*az(jp)*(ri + delta*rip)
         else
            u%twob(iptjpt) = u%twob(iptjpt) + epsi1FourPi*az(ip)*az(jp)*(ri - delta*ri)
         end if
      end do
      if (r(3,ip) < Zero) then                              ! ion--self-image interaction and Born energy
        u%oneb(ipt) = u%oneb(ipt) + (half*epsi1FourPi*az(ip)**2)*delta/(+two*r(3,ip))
        if (epsi1 < epsi2) u%oneb(ipt) = u%oneb(ipt) + (half*epsi1FourPi*az(ip)**2)*(one-eta)/radat(ipt)
      else
        u%oneb(ipt) = u%oneb(ipt) + (half*epsi2FourPi*az(ip)**2)*delta/(+two*r(3,ip))
        if (epsi1 > epsi2) u%oneb(ipt) = u%oneb(ipt) + (half*epsi2FourPi*az(ip)**2)*(one-one/eta)/radat(ipt)
      end if
   end do

   u%oneb(1:npt) = u%oneb(1:npt)/nproc
   u%twob(0) = sum(u%twob(1:nptpt))
   u%oneb(0) = sum(u%oneb(1:npt))

! ... update

   u%tot  = u%tot  + u%twob(0) +  u%oneb(0)

!  call TestUDielDisPlane(uout)

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUDielDisPlane(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,1f12.3)') 'uTotal    =',sum(u%twob(1:nptpt)) + sum(u%oneb(1:npt))
   write(unit,'(a,6f12.3)') 'uTwoBody  =',u%twob(1:nptpt)
   write(unit,'(a,3f12.3)') 'uOneBody  =',u%oneb(1:npt)
end subroutine TestUDielDisPlane

!........................................................................

end subroutine UDielDisPlane

!************************************************************************
!> \page energy energy.F90
!! **UDielDisSph**
!! *calculate coulomb energy in a system with spherical dielectric discontinuity*
!************************************************************************


subroutine UDielDisSph

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UDielDisSph'
   logical, save :: first = .true.
   integer(4) :: ipt, jpt, ip, jp, iptjpt, iploc, jploc, Getnpmyid
   real(8) :: r1, r2, r12, fac, cosa, ImageIntSph

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... setup

   if (first) then
      first = .false.
      elunit = Epsi0FourPi
      if (itest == 98) elunit = one  ! set electrostatic unit; =1 for test purpose; => value of Green's function)
      epsi1FourPi = elunit/epsi1
      epsi2FourPi = elunit/epsi2
      eta = epsi1/epsi2
      ubgfac = -nppt(1)*az(1)/boundaryrad**3 ! ubgfac: for interaction with a uniform volume charge density
   end if

! ... calculate energy

   u%twob = Zero
   u%oneb = Zero

   do iploc = 1, Getnpmyid()                           ! triangular loop (note U(q_1,q_2(image) = U(q_2,q1(image)) )
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      r1 = sqrt(r(1,ip)**2+r(2,ip)**2+r(3,ip)**2)
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (jp <= ip) cycle
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         r2 = sqrt(r(1,jp)**2+r(2,jp)**2+r(3,jp)**2)
         fac = one/(r1*r2)
         cosa = max(-one, min(fac*sum(r(1:3,ip)*r(1:3,jp)), one))
         r12 = sqrt( (r(1,ip)-r(1,jp))**2 + (r(2,ip)-r(2,jp))**2 + (r(3,ip)-r(3,jp))**2 )
                                                 ! ion--ion and ion--image interaction
         if ((r1 < boundaryrad) .and. (r2 < boundaryrad)) then
            u%twob(iptjpt) = u%twob(iptjpt) + epsi2FourPi*(az(ip)*az(jp)) * &
               (One/(r12*eta) + ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r2,cosa))
         else
            u%twob(iptjpt) = u%twob(iptjpt) + epsi2FourPi*(az(ip)*az(jp)) * &
               (One/r12 + ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r2,cosa))
         end if
      end do
                                                 ! ion--self-image interaction and Born energy
      if (r1 < boundaryrad) then
         u%oneb(ipt) = u%oneb(ipt) + (half*epsi2FourPi*az(ip)**2)*ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r1,One)
         if (epsi1 < epsi2) u%oneb(ipt) = u%oneb(ipt) + (half*epsi1FourPi*az(ip)**2)*(one-eta)/radat(ipt)
      else
         u%oneb(ipt) = u%oneb(ipt) + (half*epsi2FourPi*az(ip)**2)*ImageIntSph(lmaxdiel,boundaryrad,eta,r1,r1,One)
         if (epsi1 > epsi2) u%oneb(ipt) = u%oneb(ipt) + (half*epsi2FourPi*az(ip)**2)*(one-one/eta)/radat(ipt)
      end if

      if(lbg) then   ! interaction with the uniform volume charge density with the total charge -nppt(1)*az(1)
         if (r1 < boundaryrad) then                                                  ! inside
            u%oneb(ipt) = u%oneb(ipt) + epsi1FourPi*ubgfac*az(ip)*((half+eta)*boundaryrad**2-half*r1**2)
         else                                                                        ! outside
            u%oneb(ipt) = u%oneb(ipt) + epsi1FourPi*ubgfac*az(ip)*eta*boundaryrad**3/r1
         end if
      end if

   end do

   u%oneb(1:npt) = u%oneb(1:npt)/nproc
   u%twob(0) = sum(u%twob(1:nptpt))
   u%oneb(0) = sum(u%oneb(1:npt))

! ... update

   u%tot  = u%tot  + u%twob(0) +  u%oneb(0)

!  call TestUDielDisSph(uout)

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUDielDisSph(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,1f12.3)') 'uTotal    =',sum(u%twob(1:nptpt)) + sum(u%oneb(1:npt))
   write(unit,'(a,6f12.3)') 'uTwoBody  =',u%twob(1:nptpt)
   write(unit,'(a,3f12.3)') 'uOneBody  =',u%oneb(1:npt)
end subroutine TestUDielDisSph

!........................................................................

end subroutine UDielDisSph

!************************************************************************
!> \page energy energy.F90
!! **ImageIntSph**
!! *get the image interaction for spherical geometry*
!************************************************************************


real(8) function ImageIntSph(lmaxdiel, boundaryrad, eta, r1, r2, cosa)

   implicit none

   integer(4), intent(in) :: lmaxdiel
   real(8),    intent(in) :: boundaryrad
   real(8),    intent(in) :: eta
   real(8),    intent(in) :: r1
   real(8),    intent(in) :: r2
   real(8),    intent(in) :: cosa

   real(8), parameter :: Half = 0.5d0, One = 1.0d0, Two = 2.0d0
   character(40), parameter :: txroutine ='ImageIntSph'
   logical, save :: first = .true.
   real(8), save :: lfac1
   real(8), save :: lfac2
   real(8), save, allocatable :: lfac3(:)
   real(8), save :: lmaxdielallocated
   integer(4) :: l
   real(8)    :: fac, rdiv, t, tlow, thigh, tfac, term1, term2, term3, PL

! ... setup lfac1, lfac2, and lfac3

   if (first) then
      first = .false.
      allocate(lfac3(lmaxdiel))
      lfac3 = 0.0E+00
      lfac1 = (one-eta)/(one+eta)
      lfac2 = lfac1/(one+eta)
      do l = 1, lmaxdiel
         lfac3(l) = lfac2*eta / ((l+one) * ((one+eta)*l+one))
      end do
      lmaxdielallocated = lmaxdiel
   else
      if(lmaxdiel > lmaxdielallocated) call Stop(txroutine, 'lmaxdiel > lmaxdielallocated', 6)
   end if

! ... evaluation

   if ((r1 < boundaryrad) .and. (r2 < boundaryrad)) then        ! inside/inside
      rdiv = one/boundaryrad
      t = r1*r2*rdiv**2
      fac = sqrt(one-two*t*cosa+t**2)
      tlow = one/t
      thigh = one
      term1 = -rdiv*lfac1*thigh/(eta*fac)
   else if ((r1 > boundaryrad) .and. (r2 > boundaryrad)) then   ! ouside/ouside
      rdiv = one/boundaryrad
      t = boundaryrad**2/(r1*r2)
      fac = sqrt(one-two*t*cosa+t**2)
      tlow = one
      thigh = t
      term1 = rdiv*lfac1*thigh/fac
   else                                                         ! inside/outside
      rdiv = one/max(r1,r2)
      t = min(r1,r2)*rdiv
      fac = sqrt(one-two*t*cosa+t**2)
      tlow = one/t
      thigh = one
      term1 = rdiv*lfac1*thigh/fac
   end if

   if (t < 1.0d-10) then
      term2 = 0.0d0
   else
      if(cosa > -one+1d-10) then
         term2 = rdiv*lfac2*tlow*log((fac-t+cosa)/(one+cosa))
      else
         term2 = rdiv*lfac2*tlow*log(one-t)
      end if
   end if
   term3 = -rdiv*eta*lfac2*thigh
   tfac = rdiv*thigh
   do l = 1, lmaxdiel
      tfac = tfac*t
      term3 = term3 - lfac3(l)*tfac*PL(l,cosa)
   end do

   ImageIntSph = term1 + term2 + term3

end function ImageIntSph

!************************************************************************
!> \page energy energy.F90
!! **UTwoBodyPair**
!! *calculate two-body potential energy and force between two particles*
!************************************************************************


subroutine UTwoBodyPair(ip, jp, uuu, fforce)

   use EnergyModule
   implicit none

   integer(4), intent(in)  :: ip         ! particle 1
   integer(4), intent(in)  :: jp         ! particle 2
   real(8),    intent(out) :: uuu        ! their pair potential energy
   real(8),    intent(out) :: fforce(3)  ! force acting on particle ip

   character(40), parameter :: txroutine ='UTwoBodyPair'
   integer(4) :: ipt, jpt, iptjpt, ibuf
   integer(4) :: ia, ialow, iaupp, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum, fsum
   logical    :: EllipsoidOverlap, SuperballOverlap

   uuu = Zero
   fforce(1:3) = Zero
   ipt = iptpn(ip)
   jpt = iptpn(jp)
   iptjpt = iptpt(ipt,jpt)
   dx = ro(1,ip)-ro(1,jp)
   dy = ro(2,ip)-ro(2,jp)
   dz = ro(3,ip)-ro(3,jp)
   call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
   dx = dx-dxopbc
   dy = dy-dyopbc
   dz = dz-dzopbc
   r2 = dx**2+dy**2+dz**2
!    write(uout,'(a,3f10.5)') 'UTwoBodyPair: r2',r2
   if (r2 > rcut2) return
   if (lellipsoid) then
      if (EllipsoidOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp),radellipsoid2,aellipsoid)) then
      uuu = 1e10
      fforce = zero
      return
      end if
   end if
   if (lsuperball) then
      if (SuperballOverlap(r2,[dx,dy,dz],ori(1,1,ip),ori(1,1,jp))) then
      uuu = 1e10
      fforce = zero
      return
      end if
   end if

   if (lcharge) then   ! atoms possessing charges
      call UTwoBodyPPair
      if (lewald) call stop(txroutine, 'lewald = .true. is not appropiate', uout)
      if (lrf   ) call stop(txroutine, 'lrf = .true. is not appropiate', uout)
   else if (lweakcharge) then  ! atoms possesing weak charges (titrating system)
      call stop(txroutine, 'lwealcharge = .true. is not appropiate', uout)
   else if (ldipole) then  ! atoms possening charges and dipoles
      if(ldipole) then
         if(.not.allocated(potstat)) then
            allocate(potstat(na_alloc), estat(3,na_alloc), efg(6,na_alloc))
            potstat = 0.0E+00
            estat = 0.0E+00
            efg = 0.0E+00
         end if
      end if
      call UTwoBodyPPair
      call UDipolePair
   else if (lpolarization) then  ! atoms posessing charges, dipoles, and polarizablities
      call stop(txroutine, 'lpolarization = .true. is not appropiate', uout)
   else if (ldipolesph) then  ! atoms possening charges and dipoles, spherical boundary condition, and image charge
      call stop(txroutine, 'ldipolesph = .true. is not appropiate', uout)
   else if (ldieldis) then  ! atoms possesing charge in a system with dielectric discontinuities
      call stop(txroutine, 'ldeldis = .true. is not appropiate', uout)
   end if

contains

!........................................................................

subroutine UTwoBodyPPair

   character(40), parameter :: txroutine ='UTwoBodyPPair'
   ialow = ianpn(ip)
   iaupp = ialow+napt(ipt)-1
   jalow = ianpn(jp)
   jaupp = jalow+napt(jpt)-1
   do ia = ialow, iaupp
      iat = iatan(ia)
      do ja = jalow, jaupp
         jat = iatan(ja)
         iatjat = iatat(iat,jat)
         dx = r(1,ia)-r(1,ja)-dxopbc
         dy = r(2,ia)-r(2,ja)-dyopbc
         dz = r(3,ia)-r(3,ja)-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 < r2atat(iatjat)) then
            usum = 1d10
            fsum = zero
         else
            if (r2 < r2umin(iatjat)) call StopUTwoBodyPair
            ibuf = iubuflow(iatjat)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
                if (ibuf > nbuf) call StopIbuf('txatat',iatjat)
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                   d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
            fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                   d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))
         end if
         uuu = uuu+usum
         fforce(1) = fforce(1) + (fsum * dx)
         fforce(2) = fforce(2) + (fsum * dy)
         fforce(3) = fforce(3) + (fsum * dz)
      end do

   end do

end subroutine UTwoBodyPPair

!........................................................................

subroutine UDipolePair
   character(40), parameter :: txroutine ='UDipolePair'
   integer(4) :: ia, ialow, iaupp, jalow, jaupp
   real(8)    :: dx, dy, dz
   real(8)    :: r1, r2, r1i, r2i, r3i, r5i, r7i, threer5i, ustat, forcex, forcey, forcez
   real(8)    :: doti, dotj, dotir5i, dotjr5i, dotir7i, dotjr7i
   real(8)    :: fldx, fldy, fldz, efgxx, efgyy, efgzz, efgxy, efgxz, efgyz
   real(8), external :: ErfLocal

         potstat(1:na)   = Zero
         estat(1:3,1:na) = Zero
         ialow = ianpn(ip)
         iaupp = ialow+napt(ipt)-1
         ipt = iptpn(ip)
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         jpt = iptpn(jp)
         do ia = ialow, iaupp
            do ja = jalow, jaupp
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               r1 = sqrt(r2)
               r1i = One/r1
               r2i = r1i**2

               r3i = r2i*r1i
               r5i = r2i*r3i
               r7i = r2i*r5i
               threer5i = Three*r5i

               doti = dip(1,ia)*dx+dip(2,ia)*dy+dip(3,ia)*dz
               dotj = dip(1,ja)*dx+dip(2,ja)*dy+dip(3,ja)*dz
               dotir5i = three*r5i*doti
               dotjr5i = three*r5i*dotj
               dotir7i = (Three*Five*r7i)*doti
               dotjr7i = (Three*Five*r7i)*dotj

! ... calculate potential (potstat)

!                                  3
!     pot(r) = q / r  +  u   r  / r
!                          a   a
!                                  3
!     pot(0) = q / r  -  u   r  / r
!                          a   a

               potstat(ia) = potstat(ia) + (az(ja)*r1i + dotj*r3i)
               potstat(ja) = potstat(ja) + (az(ia)*r1i - doti*r3i)

! ... calculate field (estat)

!                   -3                         2    -5
!     e (r) = q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab
!
!                   -3                         2    -5
!     e (r) =-q r  r    +  u  ( 3 r  r  - d   r  ) r
!      a         a          b      a  b    ab

! ... calculate virial

!     vir = -q e  r
!               a  a

               fldx = (+az(ja)*dx - dip(1,ja))*r3i + dotjr5i*dx
               fldy = (+az(ja)*dy - dip(2,ja))*r3i + dotjr5i*dy
               fldz = (+az(ja)*dz - dip(3,ja))*r3i + dotjr5i*dz
               estat(1,ia) = estat(1,ia) + fldx
               estat(2,ia) = estat(2,ia) + fldy
               estat(3,ia) = estat(3,ia) + fldz
               virelec = virelec - az(ia)*(fldx*dx + fldy*dy + fldz*dz)
               fldx = (-az(ia)*dx - dip(1,ia))*r3i + dotir5i*dx
               fldy = (-az(ia)*dy - dip(2,ia))*r3i + dotir5i*dy
               fldz = (-az(ia)*dz - dip(3,ia))*r3i + dotir5i*dz
               estat(1,ja) = estat(1,ja) + fldx
               estat(2,ja) = estat(2,ja) + fldy
               estat(3,ja) = estat(3,ja) + fldz

! ... calculate field gradient (stat)

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    -  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

!                                    2   -5                         2                          -7
!     efg  (r) = - q ( 3 r  r  - d  r ) r    +  u  3 ( 5 r  r  r - r (r  d  + r  d  + r  d  ) r
!        ab               a  b    ab             g        a  b  g      a  bg   b  ag   g  ab

! ... calculate virial

!     vir  = -u  e   r
!              a  ab  b

               efgxx = az(ja)*r3i + (-az(ja)*dx**2 + Two*dip(1,ja)*dx + dotj)*threer5i - dotjr7i*dx**2
               efgyy = az(ja)*r3i + (-az(ja)*dy**2 + Two*dip(2,ja)*dy + dotj)*threer5i - dotjr7i*dy**2
               efgzz = az(ja)*r3i + (-az(ja)*dz**2 + Two*dip(3,ja)*dz + dotj)*threer5i - dotjr7i*dz**2
               efgxy =              (-az(ja)*dx*dy + dip(1,ja)*dy + dip(2,ja)*dx)*threer5i - dotjr7i*dx*dy
               efgxz =              (-az(ja)*dx*dz + dip(1,ja)*dz + dip(3,ja)*dx)*threer5i - dotjr7i*dx*dz
               efgyz =              (-az(ja)*dy*dz + dip(2,ja)*dz + dip(3,ja)*dy)*threer5i - dotjr7i*dy*dz
               efg(1,ia) = efg(1,ia) + efgxx
               efg(2,ia) = efg(2,ia) + efgyy
               efg(3,ia) = efg(3,ia) + efgzz
               efg(4,ia) = efg(4,ia) + efgxy
               efg(5,ia) = efg(5,ia) + efgxz
               efg(6,ia) = efg(6,ia) + efgyz
               virelec  = virelec -(dip(1,ia)*efgxx*dx + dip(1,ia)*efgxy*dy + dip(1,ia)*efgxz*dz &
                                  + dip(2,ia)*efgxy*dx + dip(2,ia)*efgyy*dy + dip(2,ia)*efgyz*dz &
                                  + dip(3,ia)*efgxz*dx + dip(3,ia)*efgyz*dy + dip(3,ia)*efgzz*dz)
               efgxx = az(ia)*r3i + (-az(ia)*dx**2 - Two*dip(1,ia)*dx - doti)*threer5i + dotir7i*dx**2
               efgyy = az(ia)*r3i + (-az(ia)*dy**2 - Two*dip(2,ia)*dy - doti)*threer5i + dotir7i*dy**2
               efgzz = az(ia)*r3i + (-az(ia)*dz**2 - Two*dip(3,ia)*dz - doti)*threer5i + dotir7i*dz**2
               efgxy =              (-az(ia)*dx*dy - dip(1,ia)*dy - dip(2,ia)*dx)*threer5i + dotir7i*dx*dy
               efgxz =              (-az(ia)*dx*dz - dip(1,ia)*dz - dip(3,ia)*dx)*threer5i + dotir7i*dx*dz
               efgyz =              (-az(ia)*dy*dz - dip(2,ia)*dz - dip(3,ia)*dy)*threer5i + dotir7i*dy*dz
               efg(1,ja) = efg(1,ja) + efgxx
               efg(2,ja) = efg(2,ja) + efgyy
               efg(3,ja) = efg(3,ja) + efgzz
               efg(4,ja) = efg(4,ja) + efgxy
               efg(5,ja) = efg(5,ja) + efgxz
               efg(6,ja) = efg(6,ja) + efgyz

            end do
         end do

    ustat = Zero
    do ia = ialow, ialow+napt(ipt)-1

! ... calculate electrostatic potential energy: u%stat = q pot - u e
!                                                                 a a

      ustat = ustat + az(ia)*potstat(ia) - (dip(1,ia)*estat(1,ia) + dip(2,ia)*estat(2,ia) + dip(3,ia)*estat(3,ia))

! ... calculate forces: f  = q e    + u    efg
!                        a      a      b      ab

      forcex = az(ia)*estat(1,ia) + (dip(1,ia)*efg(1,ia) + dip(2,ia)*efg(4,ia) + dip(3,ia)*efg(5,ia))
      forcey = az(ia)*estat(2,ia) + (dip(1,ia)*efg(4,ia) + dip(2,ia)*efg(2,ia) + dip(3,ia)*efg(6,ia))
      forcez = az(ia)*estat(3,ia) + (dip(1,ia)*efg(5,ia) + dip(2,ia)*efg(6,ia) + dip(3,ia)*efg(3,ia))
      fforce(1)  = fforce(1)  + EpsiFourPi*forcex
      fforce(2)  = fforce(2)  + EpsiFourPi*forcey
      fforce(3)  = fforce(3)  + EpsiFourPi*forcez

! ...  torques: t = u x e (NOT CALCULATED)

   end do

   ustat = EpsiFourPi*ustat
   uuu  = uuu  + ustat

end subroutine UDipolePair

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

!........................................................................

subroutine StopUTwoBodyPair
   character(40), parameter :: txroutine ='StopUTwoBodyPair'
   write(uout,'(a,i15)')   'iatjat        = ', iatjat
   write(uout,'(a,e15.5)') 'r2            = ', r2
   write(uout,'(a,e15.5)') 'r2umin(iatjat) = ', r2umin(iatjat)
   call Stop(txroutine, 'r2 < r2umin(iatjat)', uout)
end subroutine StopUTwoBodyPair

!........................................................................

end subroutine UTwoBodyPair

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UBond**
!! *calculate bond potential energy, forces, and virial*
!************************************************************************


subroutine UBond

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UBond'
   integer(4) :: ic, ict, iseg, ip, jp
   real(8)    :: dx, dy, dz, r1, r2, term, fac, virbond

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%bond = Zero
   virbond = Zero
   do ic = icmyid(1), icmyid(2)
      ict = ictcn(ic)
      iseg = 1
      ip = ipnsegcn(iseg,ic)
      do iseg = 2, npct(ict)
         jp = ipnsegcn(iseg,ic)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         r1 = sqrt(r2)

         term = r1-bond(ict)%eq
         fac = bond(ict)%k*term**(bond(ict)%p-1)
         u%bond = u%bond + fac*term
         fac =-bond(ict)%p*fac/r1

         force(1,ip) = force(1,ip) + (fac * dx)
         force(2,ip) = force(2,ip) + (fac * dy)
         force(3,ip) = force(3,ip) + (fac * dz)
         force(1,jp) = force(1,jp) - (fac * dx)
         force(2,jp) = force(2,jp) - (fac * dy)
         force(3,jp) = force(3,jp) - (fac * dz)
         virbond     = virbond     - (fac * r2)
         ip = jp
      end do
   end do

   if (itest == 3 .and. master) call TestUBond(uout)

! ... update

   u%tot  = u%tot  + u%bond
   virial = virial + virbond

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUBond(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,2x,es12.5)') 'u%bond ', u%bond
   write(unit,'(a,2x,es12.5)') 'virbond', virbond
   write(unit,'(a)')  'atom                   force'
   write(unit,'((i5,2x,3es12.4))') (ip, force(1:3,ip), ip=1,min(6,int(na)))
end subroutine TestUBond

!........................................................................

end subroutine UBond

!************************************************************************
!> \page energy energy.F90
!! **UAngle**
!! *calculate angle potential energy, forces, and virial*
!************************************************************************


subroutine UAngle

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UAngle'
   integer(4) :: ic, ict, iseg, ip, jm, jp
   real(8)    :: dxm, dym, dzm, dmi, dxp, dyp, dzp, dpi, cosine, theta
   real(8)    :: term, fac, forcexm, forceym, forcezm, forcexp, forceyp, forcezp

   u%angle = Zero
   if (count(angle(1:nct)%k /= Zero) == 0) return

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   do ic = icmyid(1), icmyid(2)
      ict = ictcn(ic)

      iseg = 2
      jm = ipnsegcn(iseg-1,ic)
      ip = ipnsegcn(iseg,ic)
      dxm = ro(1,ip)-ro(1,jm)
      dym = ro(2,ip)-ro(2,jm)
      dzm = ro(3,ip)-ro(3,jm)
      call PBC(dxm,dym,dzm)
      dmi = One/sqrt(dxm**2+dym**2+dzm**2)
      dxm = dxm*dmi
      dym = dym*dmi
      dzm = dzm*dmi

      do iseg = 2, npct(ict)-1
         jp = ipnsegcn(iseg+1,ic)
         dxp = ro(1,ip)-ro(1,jp)
         dyp = ro(2,ip)-ro(2,jp)
         dzp = ro(3,ip)-ro(3,jp)
         call PBC(dxp,dyp,dzp)
         dpi = One/sqrt(dxp**2+dyp**2+dzp**2)
         dxp = dxp*dpi
         dyp = dyp*dpi
         dzp = dzp*dpi

         cosine = dxm*dxp + dym*dyp + dzm*dzp
         theta = acos(max(-One,min(cosine,One)))
         term = theta-angle(ict)%eq
         fac = angle(ict)%k*term**(angle(ict)%p-1)

         u%angle = u%angle + fac*term
!        theta = min(Pi,theta)        ! to avoid division by Zero, (ok for angle%eq = Pi)
         fac = -angle(ict)%p*fac/sin(theta)
         forcexm = (fac*dmi)*(dxp-cosine*dxm)
         forceym = (fac*dmi)*(dyp-cosine*dym)
         forcezm = (fac*dmi)*(dzp-cosine*dzm)
         forcexp = (fac*dpi)*(dxm-cosine*dxp)
         forceyp = (fac*dpi)*(dym-cosine*dyp)
         forcezp = (fac*dpi)*(dzm-cosine*dzp)
         force(1,jm) = force(1,jm) + forcexm
         force(2,jm) = force(2,jm) + forceym
         force(3,jm) = force(3,jm) + forcezm
         force(1,jp) = force(1,jp) + forcexp
         force(2,jp) = force(2,jp) + forceyp
         force(3,jp) = force(3,jp) + forcezp
         force(1,ip) = force(1,ip) - (forcexm + forcexp)
         force(2,ip) = force(2,ip) - (forceym + forceyp)
         force(3,ip) = force(3,ip) - (forcezm + forcezp)

         jm = ip
         ip = jp
         dmi = dpi
         dxm = -dxp
         dym = -dyp
         dzm = -dzp
      end do
   end do

   if (itest == 3 .and. master) call TestUAngle(uout)

! ... update (no contribution to the virial)

   u%tot  = u%tot + u%angle

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine TestUAngle(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,2x,es12.5)') 'u%angle ', u%angle
   write(unit,'(a,2x,es12.5)') 'virangle', Zero
   write(unit,'(a)')  'atom                   force'
   write(unit,'((i5,2x,3es12.4))') (ip, force(1:3,ip), ip=1,min(6,int(na)))
end subroutine TestUAngle

!........................................................................

end subroutine UAngle

!************************************************************************
!> \page energy energy.F90
!! **UCrossLink**
!! *calculate crosslink potential energy, forces, and virial*
!************************************************************************


subroutine UCrossLink

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UCrossLink'
   integer(4)   :: ip, jp, icl
   real(8)      :: dx, dy, dz, r1, r2, term, fac, vircrosslink

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   u%crosslink = Zero
   vircrosslink = Zero

   do ip = ipmyid(1), ipmyid(2)
      do icl = 1, nbondcl(ip)
         jp = bondcl(icl,ip)
         if (jp < ip) cycle
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         r1 = sqrt(r2)

         term = r1-clink%eq
         fac = clink%k*term**(clink%p-1)
         u%crosslink = u%crosslink + fac*term
         fac =-clink%p*fac/r1

         force(1,ip) = force(1,ip) + (fac * dx)
         force(2,ip) = force(2,ip) + (fac * dy)
         force(3,ip) = force(3,ip) + (fac * dz)
         force(1,jp) = force(1,jp) - (fac * dx)
         force(2,jp) = force(2,jp) - (fac * dy)
         force(3,jp) = force(3,jp) - (fac * dz)
         vircrosslink = vircrosslink - (fac * r2)
      end do
   end do

   if (itest == 3 .and. master) call TestUCrossLink(uout)

! ... update

   u%tot  = u%tot  + u%crosslink
   virial = virial + vircrosslink

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine TestUCrossLink(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,2x,es12.5)') 'u%crosslink ', u%crosslink
   write(unit,'(a,2x,es12.5)') 'vircrosslink', vircrosslink
   write(unit,'(a)')  'atom                   force'
   write(unit,'((i5,2x,3es12.4))') (ip, force(1:3,ip), ip=1,min(6,int(na)))
end subroutine TestUCrossLink

!........................................................................

end subroutine UCrossLink

!**********************************************************************************************************************

!************************************************************************
!> \page energy energy.F90
!! **UExternal**
!! *calculate external potential energy*
!************************************************************************


subroutine UExternal

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UExternal'
   integer(4) :: ipt, ia, iat, ialow, iaupp, jalow, jaupp
   real(8)    :: r1, r2, fsum
   integer(4) :: imyid(1:2)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

! ... initiate

   u%external = Zero

! ... select appropriate external potential

   do ipt = 1, npt
      ialow = ianpn(ipnpt(ipt))                ! lower atom number
      iaupp = ialow+nppt(ipt)*napt(ipt)-1      ! upper atom number
      call LoadBalanceRealSpace(myid, master, nproc, ialow, iaupp, imyid, 0, 'iaptmyid', uout)  ! set imyid
      ialow = imyid(1)                         ! lower atom number of myid
      iaupp = imyid(2)                         ! upper atom number of myid
      jalow = ianpn(ipnpt(ipt))                ! lower atom number
      jaupp = jalow+nppt(ipt)*napt(ipt)-1      ! upper atom number
      if (txuext(ipt) == 'wall_z') then
         call UExternalWallZ
      else if (txuext(ipt) == 'gravitation_wall_z') then
         call UExternalGravitationWallZ
      else if (txuext(ipt) == 'superball_wall_z') then
         call UExternalSuperballWallZ
      else if (txuext(ipt) == 'laura') then
         call UExternalSuperballGravitationEstatFieldWallZ
      else if (txuext(ipt) == 'ramp_wall_z') then
         call UExternalRampWallZ
      else if (txuext(ipt) == 'sw_wall_zlow') then
         call UExternalSquareWellZlow
      else if (txuext(ipt) == 'lj_wall_z') then
         call UExternalLJWallZ
      else if (txuext(ipt) == 'lj_wall_z_ts') then
         call UExternalLJWallZts
      else if (txuext(ipt) == 'lj_wall_z_mod') then
         call UExternalLJWallZMod
      else if (txuext(ipt) == 'lj_wall_zlow') then
         call UExternalLJWallZlow
      else if (txuext(ipt) == 'lj_wall_desorb') then
         call UExternalLJWallZlowDesorb
      else if (txuext(ipt) == 'estat_field') then
         call UExternalEstatField
      else if (txuext(ipt) == 'hom_charged_walls') then
         call UExternalHomChargedWall
      else if (txuext(ipt) == 'i_soft_sphere') then
         call UExternalISoftSphere
      else if (txuext(ipt) == 'Gunnar_soft_sphere') then
         call UExternalGunnarSoftSphere
      else if (txuext(ipt) == 'out_hard_ellipsoid') then
         call UExternalOutHardEllipsoid
      else if (txuext(ipt) == 'capsid_shell') then
         call UExternalCapsidShell
      else if (txuext(ipt) == 'uniform_shell') then
         call UExternalUniformShell
      else if (txuext(ipt) == 'sphdielboundary_q') then
         call UExternalSphDielBoundary_q
      else if (txuext(ipt)(1:17) == 'sphdielboundary_p') then
         call UExternalSphDielBoundary_p(txuext(ipt)(18:18))
      else if (txuext(ipt) == 'core_shell') then
         call UExternalCoreShell
      else if (txuext(ipt) == 'insulating_sphere') then
         call UExternalInsulatingSphere
      else if (txuext(ipt) == 'hollow_sphere') then
         call UExternalHollowSphere
      else if (txuext(ipt) == 'orienting_field') then
         call UExternalOrientingField
      else if (txuext(ipt) == 'hard_cylinder') then
         call UExternalHardCylinder
      else if (txuext(ipt) == 'lekkerkerker-tuinier') then
         call UExternalLTDepletion
      end if
   end do

   if (itest == 3 .and. master) call TestUExternal(uout)

! ... update

   u%tot = u%tot + u%external

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine UExternalWallZ                       ! wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalWallZ'
   do ia = ialow, iaupp
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
   end do
end subroutine UExternalWallZ

!........................................................................

subroutine UExternalGravitationWallZ            ! graviatation force and wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalGraviationWallZ'
   do ia = ialow, iaupp
      iat = iatan(ia)
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
      u%external = u%external - gravitation_force*(r(3,ia) + wall_z_ext)
      force(3,ia) = force(3,ia) + gravitation_force
   end do
end subroutine UExternalGravitationWallZ

!........................................................................

subroutine UExternalSuperballWallZ              ! superball and nonverlapping wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalSuperballWallZ'
   real(8) :: dz, oriim(3,3)
   logical :: SuperballOverlap
   if (.not. lsuperball) call stop(txroutine, '.not.lsuperball', uout)
   do ia = ialow, iaupp
      iat = iatan(ia)
      dz = two*(r(3,ia) - sign(wall_z_ext,r(3,ia))) ! distance between atom and its image across the closest z-wall
      oriim(1:3,1:2) = ori(1:3,1:2,ipnan(ia))
      oriim(1:3,3) = -ori(1:3,3,ipnan(ia))
      if (SuperballOverlap(dz**2,[zero,zero,dz],ori(1,1,ipnan(ia)),oriim(1,1))) then
         write(*,'(a,3i8           )') 'ia, ip', ia, ipnan(ia)
         write(*,'(a,3f8.3,2x,3f8.3)') 'r, dz', r(1:3,ia), dz
         write(*,'(a,9f8.3,2x,9f8.3)') 'ori, oriim', ori(1:3,1:3,ipnan(ia)), oriim(1:3,1:3)
         call Stop(txroutine, 'superball overlap with z-wall', uout)
      end if
   end do
end subroutine UExternalSuperballWallZ

!........................................................................

subroutine UExternalSuperballGravitationEstatFieldWallZ              ! superball, gravitation, external field, and nonoverlapping wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalSuperballGravitationWallZ'
   real(8) :: dz, oriim(3,3)
   logical :: SuperballOverlap
   if (.not. lsuperball) call stop(txroutine, '.not.lsuperball', uout)
   do ia = ialow, iaupp
      iat = iatan(ia)
      dz = two*(r(3,ia) - sign(wall_z_ext,r(3,ia))) ! distance between atom and its image across the closest z-wall
      oriim(1:3,1:2) = ori(1:3,1:2,ipnan(ia))
      oriim(1:3,3) = -ori(1:3,3,ipnan(ia))
      if (SuperballOverlap(dz**2,[zero,zero,dz],ori(1,1,ipnan(ia)),oriim(1,1))) then
         write(*,'(a,3i8           )') 'ia, ip', ia, ipnan(ia)
         write(*,'(a,3f8.3,2x,3f8.3)') 'r, dz', r(1:3,ia), dz
         write(*,'(a,9f8.3,2x,9f8.3)') 'ori, oriim', ori(1:3,1:3,ipnan(ia)), oriim(1:3,1:3)
         call Stop(txroutine, 'superball overlap with z-wall', uout)
      end if
      u%external = u%external - gravitation_force*(r(3,ia) + wall_z_ext)
      force(3,ia) = force(3,ia) + gravitation_force
      u%external = u%external + az(ia)*sum(r(1:3,ia)*efield_ext(1:3))
      force(1,ia) = force(1,ia) - az(ia)*efield_ext(1)
      force(2,ia) = force(2,ia) - az(ia)*efield_ext(2)
      force(3,ia) = force(3,ia) - az(ia)*efield_ext(3)
      u%external = u%external - sum(dip(1:3,ia)*efield_ext(1:3))     ! torque should be included for dynamics
   end do
end subroutine UExternalSuperballGravitationEstatFieldWallZ

!........................................................................

subroutine UExternalRampWallZ                      ! ramp wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalRampWallZ'
   real(8) :: wall_hs, dz, width, slope, uterm, fterm
   do ia = ialow, iaupp
      iat = iatan(ia)
      wall_hs = wall_z_ext - radat(iat)            ! hs interaction at z = wall_hs
      if (abs(r(3,ia)) > wall_hs) call Stop(txroutine, 'hs overalp with z-wall', uout)
      width = (lambda_ramp_ext-One)*Two*radat(iat)
      dz = abs(r(3,ia)) - (wall_hs-width)
      if (dz > Zero .and. dz <= width) then         ! on the slope
         slope = epsilon_ramp_ext/width
         uterm = -slope*dz
         fterm = sign(One,-r(3,ia))*slope
         u%external = u%external + uterm
         force(3,ia) = force(3,ia) + fterm
      end if
   end do
end subroutine UExternalRampWallZ

!........................................................................

subroutine UExternalSquareWellZlow                 ! square-well wall at z = -wall_z_ext
   character(40), parameter :: txroutine ='UTwoBodynlSquareWellZlow'
   real(8) :: wall_hs, width, dz, uterm, fterm
   do ia = ialow, iaupp
      iat = iatan(ia)
      wall_hs = wall_z_ext - radat(iat)            ! hs interaction at z = wall_hs
      if (abs(r(3,ia)) > wall_hs) then
          write(*,*) 'istep1, istep2', istep1, istep2
          write(*,*) 'ia, r(3,ia)',ia, r(3,ia)
          call Stop(txroutine, 'hs overlap with z-wall', uout)
      end if
      width = (lambda_sw_ext-One)*Two*radat(iat)
      dz = r(3,ia) + (wall_hs-width)               ! dz < 0 => in the well
      uterm = -epsilon_sw_ext*(Half - (One/Pi)*atan(dz/alpha_sw_ext))
      fterm = epsilon_sw_ext/Pi * alpha_sw_ext / (alpha_sw_ext**2 + dz**2)
      u%external = u%external + uterm
      force(3,ia) = force(3,ia) + fterm
   end do
end subroutine UExternalSquareWellZlow

!........................................................................

subroutine UExternalLJWallZ                      ! LJ wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalLJWallZ'
   real(8) :: z1, zi1, zi3, ulj, flj
   do ia = ialow, iaupp
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
      iat = iatan(ia)
      z1 = wall_z_ext-abs(r(3,ia))           ! lj-wall at +-wall_z_ext
      zi1 = One/z1
      zi3 = zi1**3
      ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3
      flj = (Three*z3coeff_ext(iat)*zi3 + 9.0d0*z9coeff_ext(iat)*zi3**3)*zi1
      if (r(3,ia) > Zero) flj = -flj
      u%external = u%external + ulj
      force(3,ia) = force(3,ia) + flj
   end do
end subroutine UExternalLJWallZ

!........................................................................

subroutine UExternalLJWallZts                ! truncated and shifted LJ wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UTwoBodynlLJWallZts'
   real(8) :: z1, zi1, zi3, ulj, flj
   do ia = ialow, iaupp
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
      iat = iatan(ia)
      z1 = wall_z_ext-abs(r(3,ia))           ! lj-wall at +-wall_z_ext
      if (z1  <  zmin_ext(iat)) then
         zi1 = One/z1
         zi3 = zi1**3
         ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3 + delta_ext(iat)
         flj = (Three*z3coeff_ext(iat)*zi3 + 9.0d0*z9coeff_ext(iat)*zi3**3)*zi1
         if (r(3,ia) > Zero) flj = -flj
         u%external = u%external + ulj
         force(3,ia) = force(3,ia) + flj
      end If
   end do
end subroutine UExternalLJWallZts

!........................................................................

subroutine UExternalLJWallZMod             ! heterogeneous LJ wall at abs(z) = wall_z_ext
   character(40), parameter :: txroutine ='UExternalLJWallZMod'
   real(8) :: z1, zi1, zi3, ulj, flj, hx, hy, Onephxhy, utemp, ftemp(3)
   do ia = ialow, iaupp
      if (abs(r(3,ia)) > wall_z_ext) write(*,*) 'ia,r(3,ia)',ia,r(3,ia)
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
      iat = iatan(ia)
      z1 = wall_z_ext-abs(r(3,ia))           ! lj-wall at +-wall_z_ext
      zi1 = One/z1
      zi3 = zi1**3
      ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3
      flj = (Three*z3coeff_ext(iat)*zi3 + 9.0d0*z9coeff_ext(iat)*zi3**3)*zi1
      hx = cos(TwoPilxi_ext*r(1,ia))
      hy = cos(TwoPilyi_ext*r(2,ia))
      Onephxhy = One + hx*hy
      if (z1 < zmin_ext(iat) ) then
         utemp = ulj+(delta_ext(iat)*c_ext)*Onephxhy
         ftemp(1) = (delta_ext(iat)*c_ext)*sin(TwoPilxi_ext*r(1,ia))*hy*TwoPilxi_ext
         ftemp(2) = (delta_ext(iat)*c_ext)*sin(TwoPilyi_ext*r(2,ia))*hx*TwoPilyi_ext
         ftemp(3) = flj
      else
         utemp = ulj*(One-c_ext*Onephxhy)
         ftemp(1) = -ulj*c_ext*sin(TwoPilxi_ext*r(1,ia))*hy*TwoPilxi_ext
         ftemp(2) = -ulj*c_ext*sin(TwoPilyi_ext*r(2,ia))*hx*TwoPilyi_ext
         ftemp(3) = flj*(One-c_ext*Onephxhy)
      end if
      if (itest == 87) then
         write(*,'(a,3f10.5)') 'c_ext',c_ext
         write(*,'(a,3f10.5)') 'lx_ext,ly_ext',lx_ext,ly_ext
         write(*,'(a,3f10.5)') 'r(1:3,ia)',r(1:3,ia)
         write(*,'(a,3f10.5)') 'hx',hx
         write(*,'(a,3f10.5)') 'hx',hy
         write(*,'(a,3f10.5)') 'Onephxhy',Onephxhy
         write(*,'(a,3f10.5)') 'ulj, utemp',ulj, utemp
         write(*,'(a,3f10.5)') 'flj',flj
         write(*,'(a,3f10.5)') 'ftemp(1:3)',ftemp(1:3)
         stop
      end if
      u%external = u%external + utemp
      if (r(3,ia) > Zero) ftemp(3) = -ftemp(3)
      force(1:3,ia) = force(1:3,ia) + ftemp(1:3)
   end do
end subroutine UExternalLJWallZMod

!........................................................................

subroutine UExternalLJWallZlow           ! LJ + ramp wall at abs(z) = -wall_z-ext
   character(40), parameter :: txroutine ='UExternalLJWallZlow'
   real(8) :: zi1, zi3, ulj, flj, udrift, fdrift
   do ia = ialow, iaupp
      if (abs(r(3,ia)) > wall_z_ext) call Stop(txroutine, 'outside external z-wall', uout)
      iat = iatan(ia)
      zi1 = One/(r(3,ia)-(-wall_z_ext))  ! lj-wall at -wall_z_ext
      zi3 = zi1**3
      ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3
      flj = (Three*z3coeff_ext(iat)*zi3 + 9.0d0*z9coeff_ext(iat)*zi3**3)*zi1
      udrift = u0_drift + u1_drift*r(3,ia)
      fdrift = -u1_drift
      u%external = u%external + ulj + udrift
      force(3,ia) = force(3,ia) + flj + fdrift
   end do
end subroutine UExternalLJWallZlow

!........................................................................

subroutine UExternalLJWallZlowDesorb                    ! Niklas 2006-12-20 desorb
   character(40), parameter :: txroutine ='UTwoBodyA'
   real(8), parameter :: Sixth  = 1.0d0/6.0d0, TwoThree = 2.0d0/3.0d0
   real(8), allocatable, save :: z_min(:), uz_min(:), wallcut(:)
   real(8) :: zi1, zi3, ulj, flj
   if(.not.allocated(z_min)) then
      allocate(z_min(nat), uz_min(nat), wallcut(nat))
      z_min = 0.0E+00
      uz_min = 0.0E+00
      wallcut = 0.0E+00
   end if
   do iat = 1, nat
      z_min(iat) = ((-3*z9coeff_ext(iat) / z3coeff_ext(iat))**(Sixth))
      uz_min(iat) = TwoThree*((-z3coeff_ext(iat)**3) / (3*z9coeff_ext(iat)))**(Half)
      wallcut(iat) = -wall_z_ext + z_min(iat)
   end do
   do ia = ialow, iaupp
      iat = iatan(ia)
      if (abs(r(3,ia)) > wall_z_ext) then
         u%external = u%external + 1.0d10               ! mimick hs-overlap
      else if (r(3,ia) > wallcut(iat)) then              ! truncation at U_min
         cycle
      else
         zi1 = One/(r(3,ia)-(-wall_z_ext))  ! lj-wall at -wall_z_ext
         zi3 = zi1**3
         ulj = z3coeff_ext(iat)*zi3 + z9coeff_ext(iat)*zi3**3 + uz_min(iat)
         flj = (Three*z3coeff_ext(iat)*zi3 + 9.0d0*z9coeff_ext(iat)*zi3**3)*zi1
         u%external = u%external + ulj
         force(3,ia) = force(3,ia) + flj
      end if
   end do
end subroutine UExternalLJWallZlowDesorb

!........................................................................

subroutine UExternalEstatField        ! homogeneous electrical field
   do ia = ialow, iaupp
      u%external = u%external + az(ia)*sum(r(1:3,ia)*efield_ext(1:3))
      force(1,ia) = force(1,ia) - az(ia)*efield_ext(1)
      force(2,ia) = force(2,ia) - az(ia)*efield_ext(2)
      force(3,ia) = force(3,ia) - az(ia)*efield_ext(3)
      u%external = u%external - sum(dip(1:3,ia)*efield_ext(1:3))     ! torque should be included for dynamics
   end do
end subroutine UExternalEstatField

!........................................................................

! Calculation of electrostatic energy (eq. 16 of Jonsson et al, JPC, 84, 2179, 1980)
! Contribution from the interaction with the walls and long-range contribution

subroutine UExternalHomChargedWall
   real(8) :: scd, a, a2, b, z, zsb, zsb2, longrangecontr
   integer(4) :: i
   scd = surfchargeden/(ech*1.d20)
   a = boxlen2(1)
   a2 = a**2
   b = boxlen2(3)
   do ia = ialow, iaupp
      iat = iatan(ia)
      z = r(3,ia)
      do i = -1,1,2
         zsb = abs(i*b-z)
         zsb2 = zsb**2
         u%external = u%external + EpsiFourPi*zat(iat)*scd*(8.d0*a*log( (sqrt(Two*a2+zsb2)+a)/(sqrt(a2+zsb2)))- &
         Two*zsb*(asin( (a2**2-zsb2**2-Two*a2*zsb2)/(a2+zsb2)**2 )+Half*pi))
      end do
      if (llongrangecontr) then
         u%external = u%external + EpsiFourPi*zat(iat)*longrangecontr(boxlen2(1), z, scd, mninchden, zdist, chden)
      endif
   end do
end subroutine UExternalHomChargedWall

!........................................................................

subroutine UExternalISoftSphere        ! external, soft, and spherical wall
   do ia = ialow, iaupp
      r1 = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if (r1 > ruext(1,ipt)) then
         fsum = auext*(r1-ruext(1,ipt))**(nuext-1)
         u%external = u%external + fsum*(r1-ruext(1,ipt))
         fsum = nuext*fsum/r1
         force(1,ia) = force(1,ia) - fsum*r(1,ia)
         force(2,ia) = force(2,ia) - fsum*r(2,ia)
         force(3,ia) = force(3,ia) - fsum*r(3,ia)
      end if
   end do
end subroutine UExternalISoftSphere

!........................................................................

subroutine UExternalGunnarSoftSphere ! external, soft, and spherical wall (Gunnar)
   do ia = ialow, iaupp
      r1 = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      fsum = auext*(r1/ruext(1,ipt))**nuext
      u%external = u%external + fsum
      fsum = nuext*ruext(1,ipt)*fsum/r1**2
      force(1,ia) = force(1,ia) - fsum*r(1,ia)
      force(2,ia) = force(2,ia) - fsum*r(2,ia)
      force(3,ia) = force(3,ia) - fsum*r(3,ia)
   end do
end subroutine UExternalGunnarSoftSphere

!........................................................................

subroutine UExternalOutHardEllipsoid    ! hard ellipsoidal wall
   character(40), parameter :: txroutine ='UExternalOutHardEllipsoid'
   do ia = ialow, iaupp
      r2 = (r(1,ia)*ruexti(1,ipt))**2+(r(2,ia)*ruexti(2,ipt))**2+(r(3,ia)*ruexti(3,ipt))**2
      if (r2 < One) call Stop(txroutine, 'violation of hard-ellipsoid potential', uout)
   end do
end subroutine UExternalOutHardEllipsoid

!........................................................................

subroutine UExternalCapsidShell        ! hard spherical capsid shell
   character(40), parameter :: txroutine ='UExternalCapsidShell'
   do ia=ialow,iaupp
      r1=sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if (r1>=(rcap-radat(ipt)) .and. r1<=(rcap+dcap+radat(ipt))) &
         call Stop(txroutine, 'hard-core overlap with capsid', uout)
   end do
end subroutine UExternalCapsidShell

!........................................................................

subroutine UExternalUniformShell        ! hard capsid shell with a uniform surface charge density
   character(40), parameter :: txroutine ='UExternalUniformShell'
   real(8), save :: rcapchage = 2.0d0   ! radial location of capside charge from inner surface
   do ia=ialow,iaupp
      r1=sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if (r1 < (rcap-radat(ipt))) then
         u%external = u%external + zat(ipt)*zat(1)*EpsiFourPi/(rcap+rcapchage)
      else if (r1 >= (rcap-radat(ipt)) .and. r1 <= (rcap+dcap+radat(ipt))) then
         call Stop(txroutine, 'hard-core overlap with capsid', uout)
      else if (r1 > (rcap+dcap+radat(ipt))) then
         u%external = u%external + zat(ipt)*zat(1)*EpsiFourPi/r1
      end if
   end do
end subroutine UExternalUniformShell

!........................................................................

subroutine UExternalSphDielBoundary_q   ! dielectric sphere (multipole expansion)
   character(40), parameter :: txroutine ='UExternalSpheDielBoundary_q'
   real(8) :: fac, rfac, rr, theta, phi, rratio, sum, term
   complex(8) xCCLM
   logical :: first = .true.
   integer(4) :: l, m

   if (slave) return  ! generation and use of QQ are made on master

! ... setup lfac and mfac

   if (first) then
      if (lmaxdiel > mlmax) call Stop(txroutine, 'lmaxdiel > mlmax', uout)
      do l = 0, lmaxdiel
         fac = real(l)/real(l+1)
        lfac(l) = fac*(epsi2-epsi1)/(epsi2+fac*epsi1)/boundaryrad  ! from solution of Poissons equ
        lfac(l) = (EpsiFourPi/Two)*lfac(l)                         ! unit transformation and div by 2 xyz
      end do
      mfac = two
      mfac(0) = one
      first = .false.
   end if

! ... calculate the multipole moments

   QQ(0:lmaxdiel,0:lmaxdiel) = cmplx(Zero, Zero)
   ialow = ianpn(ipnpt(ipt))                ! lower atom number
   iaupp = ialow+nppt(ipt)*napt(ipt)-1      ! upper atom number
   do ia = ialow, iaupp
      rratio = boundaryrad/sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      call CarToSph('rad',r(1,ia),r(2,ia),r(3,ia),rr,theta,phi)
      rfac = One
      do l = 0, lmaxdiel
         rfac = rfac*rratio
         do m = 0, l
            QQ(m,l) = QQ(m,l) + az(ia)*rfac*xCCLM(l,m,theta,phi,0)
         end do
      end do
   end do
!      call WriteStd(mlmax,lmaxdiel,QQ,'QQ',0,6)

! ... calculate external energy

   do l = 0, lmaxdiel
      sum = zero
      do m = 0,l
         term = mfac(m)*(real(QQ(m,l))**2 + aimag(QQ(m,l))**2)
         sum = sum + term
      end do
      u%external = u%external + lfac(l)*sum
   end do

end subroutine UExternalSphDielBoundary_q

!........................................................................

subroutine UExternalSphDielBoundary_p(str)  ! dielectric sphere (pairwise interaction)
   character(1), intent(in) :: str          ! select equations to be used
                                            ! '1': short expansion
                                            ! '2': long expansion
   character(40), parameter :: txroutine ='UExternalSphDielBoundary_p'
   real(8) :: r1, r2, cosa, t, fac, tfac, sum, term1, term2, term3, PL
   logical :: first = .true.
   integer(4) :: ia, ja, l

! ... setup llfac:s

   if (first) then
      if (lmaxdiel > mlmax) call Stop(txroutine, 'lmaxdiel > mlmax', uout)
      fac = (EpsiFourPi/Two)/boundaryrad
      eta = epsi1/epsi2
      lfac1 = fac*(one-eta)/(one+eta)
      lfac2 = fac*(one-eta)/(one+eta)**2
      if (str == '1') then
         do l = 1, lmaxdiel
            lfac3(l) = -fac * (one/(one+eta))**2 / l**2 * ((one-eta)/(one+one/l+eta))
         end do
      else if (str == '2') then
         do l = 1, lmaxdiel
            lfac3(l) = fac * ((one-eta)/(one+eta))**2 * eta / l**2 * (one/(one+one/l+eta)) / (one+one/l)
         end do
      end if
      first = .false.
      if (itest == 98) then
         call WriteHead(3, txroutine, 6)
         write(*,*) 'str = ', str
         write(*,'(a,8f12.7)') 'fac',fac
         write(*,'(a,8f12.7)') 'eta',eta
         write(*,'(a,8f12.7)') 'lfac1',lfac1
         write(*,'(a,8f12.7)') 'lfac2',lfac2
         write(*,'(a,8f12.7)') 'lfac3',lfac3(1:min(8,lmaxdiel))
      end if
   end if

! ... calculate external energy

   sum = zero
   do ia = ialow, iaupp
      r1 = sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      do ja = jalow, jaupp
         r2 = sqrt(r(1,ja)**2+r(2,ja)**2+r(3,ja)**2)
         fac = one/(r1*r2)
         cosa = (r(1,ia)*r(1,ja)+r(2,ia)*r(2,ja)+r(3,ia)*r(3,ja))*fac
         cosa = max(-one, min(cosa, one))
         t = boundaryrad**2*fac
         fac = sqrt(one-two*cosa*t+t**2)
         if (str == '1') then
            term1 = lfac1*t*(one/fac - one)
            term2 = lfac2*t*log(half*(fac+one-cosa*t))
            term3 = zero
         else if (str == '2') then
            term1 = lfac1*t*(one/fac - one)
            if (ja /= ia) term2 = lfac2*log((one-cosa)/(fac+t-cosa))
            if (ja == ia) term2 = lfac2*log(one-t)
            term3 = lfac2*t
         end if
         if (itest == 98) then
            write(*,*)
            write(*,'(a,i5,2f12.5)') 'ia, r1            ',ia, r1
            write(*,'(a,i5,2f12.5)') 'ja, r2            ',ja, r2
            write(*,'(a,2f12.7)') 'l, term1, sum:        0',term1, term1
            write(*,'(a,2f12.7)') 'l, term2, sum:        0',term2, term1+term2
            write(*,'(a,2f12.7)') 'l, term3, sum:        0',term3, term1+term2+term3
         end if
         tfac = t
         do l = 1, lmaxdiel
            tfac = tfac*t
            term3 = term3 - lfac3(l)*tfac*PL(l,cosa)
            if (itest == 98) write(*,'(a,i5,2f12.7)') 'l, term3, sum:    ',l, term3, term1+term2+term3
         end do
         sum =  sum + (az(ia)*az(ja))*(term1 + term2 + term3)
      end do
   end do
   u%external = u%external + sum

end subroutine UExternalSphDielBoundary_p

!........................................................................

subroutine UExternalCoreShell     ! hard inner and outer spherical walls
   character(40), parameter :: txroutine ='UExternalCoreShell'
   do ia = ialow, iaupp
      r1=sqrt(r(1,ia)**2+r(2,ia)**2+r(3,ia)**2)
      if ((r1 < rChargeIn) .or. (r1 > rChargeOut)) then
         call Stop(txroutine, 'hard-core overlap with core/shell', uout)
      end if
   end do
end subroutine UExternalCoreShell

!........................................................................

subroutine UExternalInsulatingSphere ! penetrable uniformly charged sphere
   real(8):: sum
   sum = Zero
   do ia = ialow, iaupp
      iat = iatan(ia)
      r2=r(1,ia)**2+r(2,ia)**2+r(3,ia)**2
      r1=sqrt(r2)
      if (r1 < rInsSphere) then
         sum = sum + zat(iat)/(two*rInsSphere)*(three-r2/rInsSphere**2)
      else
         sum = sum + zat(iat)/r1
      end if
   end do
   u%external = (zInsSphere*EpsiFourPi)*sum
end subroutine UExternalInsulatingSphere

!........................................................................

subroutine UExternalHollowSphere
   real(8):: sum
   sum = Zero
   do ia = ialow, iaupp
     iat = iatan(ia)
     r2 = r(1,ia)**2+r(2,ia)**2+r(3,ia)**2
     r1 = sqrt(r2)
     if (r1 < rInSphere) then
         sum = sum - zat(iat)*threehalf*(rInSphere**2-rOutSphere**2)/(rOutSphere**3-rInSphere**3)
     else if (r1 <= rOutSphere) then
         sum = sum - zat(iat)/(rOutSphere**3-rInSphere**3) * &
            (half*r2+rInSphere**3/r1-threehalf*rInSphere**2+threehalf*(rInSphere**2-rOutSphere**2))
     else if (r1 > rOutSphere) then
         sum = sum + zat(iat)/r1
     end if
   end do
   u%external = (zInsSphere*EpsiFourPi)*sum
end subroutine UExternalHollowSphere

!........................................................................

!  1) torque should be included for dynamics
!  2) works only for a single particle type

subroutine UExternalOrientingField
   integer :: ip
      do ip=1, np
         u%external = u%external - ofstrength*(ofaxis(1)*ori(3,1,ip)+ofaxis(2)*ori(3,2,ip)+ofaxis(3)*ori(3,3,ip))
      end do
end subroutine UExternalOrientingField

!........................................................................

subroutine UExternalHardCylinder
   character(40), parameter :: txroutine ='UExternalHardCylinder'
   real(8) :: rl2
   do ia = ialow, iaupp
      r2 = r(1,ia)**2+r(2,ia)**2
      if (r2 < rCylinder**2) call Stop(txroutine, 'hard-core overlap with cylinder', uout)  ! check if particle is inside cylinder     (forbidden)
      if (zCylinder /= zero) then
        rl2 = sqrt(r2 + cyllen2**2)
        u%external = u%external + EpsiFourPi*zCylinder*zat(iatan(ia))*log((rl2+cyllen2)/(rl2-cyllen2))
      end if
   end do
end subroutine UExternalHardCylinder

!........................................................................

subroutine UExternalLTDepletion
   character(40), parameter :: txroutine ='UExternalLTDepletion'
   real(8) :: hlow, hupp, dhext, LinInter
   real(8) :: h, u0, u1, u2
   logical, save :: first = .true.
   integer(4) :: ih

   if (first) then ! generate external potential table
      first = .false.
      iat = 1
      hlow = zero
      hupp = two*rad_dep
      dhext = (hupp - hlow)/nhext
      dhexti = one/dhext
      do ih = 0, nhext
         h = hlow + ih*dhext
         call calc_LT_depletion('sw', radat(iat), rad_dep, rho_dep, factor_dep, h, u0, u1, u2)
         hext(ih) = h
         uext(ih) = u0/beta
      end do
   end if

   do ia = ialow, iaupp
      iat = iatan(ia)
      h = boxlen2(3) - abs(r(3,ia))        ! select nearest surface
      if (h < zero) call Stop(txroutine, 'outside z-wall', uout)
      ih = int(h*dhexti)
      if (ih < nhext) u%external = u%external + LinInter(hext(ih), hext(ih+1), h, uext(ih), uext(ih+1))
   end do
end subroutine UExternalLTDepletion

!........................................................................

subroutine TestUExternal(unit)
   integer(4), intent(in) :: unit
   call WriteHead(3, 'Test'//trim(txroutine), unit)
   write(unit,'(a,2x,f12.7)') 'u%external', u%external
   write(unit,'(a)')  'atom                   force'
   write(unit,'((i5,2x,3f12.7))') (ia, force(1:3,ia), ia=1,min(6,int(na)))
end subroutine TestUExternal

!........................................................................

end subroutine UExternal

!************************************************************************
!> \page energy energy.F90
!! **UreactionSphere**
!! *calculate the reaction energy of an ion and a dielectric sphere*
!************************************************************************


real(8) function UreactionSphere(fac, e1, e2, r, rad, tol,itest,unit)

   implicit none

   real(8),    intent(in) :: fac
   real(8),    intent(in) :: e1
   real(8),    intent(in) :: e2
   real(8),    intent(in) :: r
   real(8),    intent(in) :: rad
   real(8),    intent(in) :: tol
   integer(4), intent(in) :: itest
   integer(4), intent(in) :: unit

   real(8) :: term, termold, lfac
   integer(4) :: i,l

   if (itest==1) write(unit,'(a)') 'Legendre summation'
   UreactionSphere = 0.0d0
   termold = 1.d10
   l = -1
   do i = 1, 100
     l = l + 1
     lfac = real(l)/real(l+1)
     if (r < rad) then
        term = -0.50d0*fac*(e2-e1)/(e2+lfac*e1) / (e1*rad) * (r/rad)**(2*l)
     else
        term = 0.50d0*fac*lfac*(e2-e1)/(e2+lfac*e1) / (e2*r) * (rad/r)**(2*l+1)
     end if
     UreactionSphere = UreactionSphere + term
     if (itest==1) write(unit,'(a,i2,2f14.8)') 'l, term,  sum:       ',l, term, UreactionSphere
     if (abs(term-termold) < tol) exit
     termold = term
   end do

end function UreactionSphere

!************************************************************************
!> \page energy energy.F90
!! **UExternalUpdate**
!! *update arrays for external energy*
!************************************************************************


subroutine UExternalUpdate(iptmove)

   use EnergyModule
   implicit none

   integer(4), intent(in) :: iptmove

   if (txuext(iptmove) == 'sphdielboundary_q') call DUExternalSphDielBoundaryUpdate

contains

!........................................................................

subroutine  DUExternalSphDielBoundaryUpdate
   QQ(0:lmaxdiel,0:lmaxdiel) = QQtm(0:lmaxdiel,0:lmaxdiel)
end subroutine  DUExternalSphDielBoundaryUpdate

!........................................................................

end subroutine UExternalUpdate

!************************************************************************
!> \page energy energy.F90
!! **xCCLM**
!! *return the spherical harmonics*
!************************************************************************

!             norm = 0,  modified spherical harmonics c(l,m)(theta,phi)
!             nrom = 1,  spherical harmonics y(l,m)(theta,phi)

complex(8) function xCCLM(l,m,theta,phi,norm)

   implicit none

   integer(4), intent(in) :: l, m, norm
   real(8),    intent(in) :: theta, phi

   real(8), parameter :: zero = 0.0d0 , one = 1.0d0
   real(8), parameter :: pi = 3.14159265359d0 , facpi = 4.0d0*pi
   integer(4), parameter :: lmax = 100, l2max = 2*lmax
   logical, save :: first=.true.
   real(8), save :: s(0:l2max), si(0:l2max), vnorm(0:lmax)
   integer(4) :: i, mabs
   real(8) :: x, xPLM

   if (first) then
      s(0) = one
      si(0) = one
      do i = 1, l2max
         s(i) = s(i-1)*sqrt(dble(i))     ! sqrt(n!)
         si(i) = one/s(i)                ! 1/sqrt(n!)
      end do
      do i = 0, lmax
         vnorm(i) = sqrt((2*i+one)/facpi)
      end do
      first = .false.
   end if
   mabs = abs(m)
   x = cos(theta)
   xCCLM = (-1)**mabs*s(l-mabs)*si(l+mabs)*xPLM(l,mabs,x)*cdexp(dcmplx(zero,mabs*phi))
   if (m < 0) xCCLM = (-1)**mabs*conjg(xCCLM)
   if (norm == 1) xCCLM = vnorm(l)*xCCLM
end function xCCLM

!************************************************************************
!> \page energy energy.F90
!! **xPLM**
!! *return the associate legendre polynomial p(l,m) at x*
!************************************************************************

!     p(m,m)(x) = (2m-1)!!(1-x*x)*(m/2) and the recurrence relation
!     (l-m+1)*p(l+1,m)(x)=(2*l+1)*x*p(l,m)(x)-(l+m)*p(l-1,m)(x)

real(8) function xPLM(l,m,x)

   implicit none

   integer(4), intent(in) :: l, m
   real(8),    intent(in) :: x

   real(8), parameter :: one = 1.0d0 , two = 2.0d0
   integer(4) :: i, j
   real(8)    :: pmm, factor, fac, pmp1m, pjm

   pmm = one
   if (m > 0) then
      factor = one
      fac = sqrt((one-x)*(one+x))
      do i = 1,m
         pmm = pmm*factor*fac
         factor = factor+two
      end do
   end if
   if (l == m) then
      xPLM = pmm
   else
      pmp1m = (2*m+1)*x*pmm
      if (l == m+1) then
         xPLM = pmp1m
      else
         do j = m+2,l
!        do j = l,l
            pjm = (x*(2*j-1)*pmp1m-(j+m-1)*pmm)/(j-m)
            pmm = pmp1m
            pmp1m = pjm
         end do
         xPLM = pjm
      end if
   end if
end function xPLM

!************************************************************************
!> \page energy energy.F90
!! **Longrangecontr**
!! *calculate long-range contribution of the electrostatic energy*
!************************************************************************

!     (eq. 17 of Jonsson et al, JPC, 84, 2179, 1980)

real(8) function Longrangecontr(a, z, scd, n, zdist, chden)

   implicit none

   real(8),    intent(in) :: a          ! half box length i z-direction
   real(8),    intent(in) :: z          ! z-coordinate
   real(8),    intent(in) :: scd        ! surface charge density
   integer(4), intent(in) :: n          ! number of points
   real(8),    intent(in) :: zdist(*)   ! z-value
   real(8),    intent(in) :: chden(*)   ! charge density

   real(8), parameter :: pi = 3.14159265359d0
   real(8) :: a2, sum, h, g, z1, dz, dz2
   integer(4) :: i

   a2 = a**2
   sum = 0.d0

! ... contribution from ion distribution

   do i = 1, n-1
     h = abs(zdist(i+1)-zdist(i))
     g = 0.5d0*(chden(i+1)+chden(i))
     z1 = 0.5d0*(zdist(i+1)+zdist(i))
     dz = abs(z-z1)
     dz2 = dz**2
     sum = sum + h*g*(-pi*dz-8.d0*a*log( (sqrt(2.d0*a2+dz2)+a)/(sqrt(a2+dz2))) &
              +2.d0*dz*asin( (a2**2-dz2**2-2.d0*a2*dz2)/(a2+dz2)**2))
   end do

! ... contribution from walls

   do i = 1, n, n-1
      z1 = zdist(i)
      dz = abs(z-z1)
      dz2 = dz**2
      sum = sum + scd*(-pi*dz-8.d0*a*log( (sqrt(2.d0*a2+dz2)+a)/(sqrt(a2+dz2))) &
                  +2.d0*dz*asin( (a2**2-dz2**2-2.d0*a2*dz2)/(a2+dz2)**2))
   end do

   longrangecontr = sum

end function Longrangecontr

!************************************************************************
!> \page energy energy.F90
!! **SPMFFTRec**
!! *make Fourier transformation, reciprocal space operations, and back FFT*
!************************************************************************


# ifdef F03_CBIND
subroutine SPMFFTRec(lsave, linit, lenergy, txFFT, txRec, level, uloc, virloc)

   use EnergyModule, s=>meshsize
   implicit none

   logical, intent(in) :: lsave       ! .true. => save QMesh in QMeshTM
   logical, intent(in) :: linit       ! .true. => initialize reciprocal energy
   logical, intent(in) :: lenergy     ! .true. => calcuate reciprocal energy and reciprocal virial
   character(*), intent(in) :: txFFT  ! label of FFT for CpuAdd
   character(*), intent(in) :: txRec  ! label of real space operations for CpuAdd
   integer(4), intent(in) :: level    ! indentation level for CpuAdd
   real(8), intent(out) :: uloc       ! reciprocal energy
   real(8), intent(out) :: virloc     ! reciprocal virial

   integer(4) :: nx, ny, nz, nya, nza
   real(8)    :: Qsum
   complex(8)  :: v1, v2, v3, v4

   if (ltime) call CpuAdd('start', txFFT, level, uout)
   call fftw_execute_dft_r2c(plan_fwd, Qmesh, FQmesh)   ! transform Qmesh into FQmesh
   if (ltime) call CpuAdd('stop',  txFFT, level, uout)

   if (lsave) QMeshTM = QMesh
   if (linit) uloc = Zero

! ... calculate the generalized influence function in the reciprocal space

   if (ltime) call CpuAdd('start', txRec, level, uout)
   uloc = Zero             ! uloc needs to be initialized!
   do nz = 0,s(3)/2
      nza = mod(s(3)-nz,s(3))
      do ny = 0,s(2)/2
         nya = mod(s(3)-ny,s(3))
         do nx = 0,s(1)/2
            v1 = FQMesh(nx+1, ny +1, nz +1)
            v2 = FQMesh(nx+1, ny +1, nza+1)
            v3 = FQMesh(nx+1, nya+1, nz +1)
            v4 = FQMesh(nx+1, nya+1, nza+1)

            if (master .and. lenergy) then
               Qsum = real(v1)**2+aimag(v1)**2 &
                    + real(v2)**2+aimag(v2)**2 &
                    + real(v3)**2+aimag(v3)**2 &
                    + real(v4)**2+aimag(v4)**2
               uloc = uloc + Qsum*energyfac(nx,ny,nz)
               virloc = virloc - Qsum*virfac(nx,ny,nz)
            end if
            FQMesh(nx+1, ny +1, nz +1 ) =  v1 * meshfac(nx,ny,nz)
            FQMesh(nx+1, ny +1, nza+1 ) =  v2 * meshfac(nx,ny,nz)
            FQMesh(nx+1, nya+1, nz +1 ) =  v3 * meshfac(nx,ny,nz)
            FQMesh(nx+1, nya+1, nza+1 ) =  v4 * meshfac(nx,ny,nz)
         end do
      end do
   end do
   if (ltime) call CpuAdd('stop', txRec, level, uout)

! ... calculate the convoluted (Green function)*(multipole) in real space

   if (ltime) call CpuAdd('start', txFFT, level, uout)
   call fftw_execute_dft_c2r(plan_bwd, FQMesh, QMesh )
   if (ltime) call CpuAdd('stop',  txFFT, level, uout)

end subroutine SPMFFTRec
# endif

!************************************************************************
!> \page energy energy.F90
!! **EwaldSetup**
!! *setup for ewald summation*
!************************************************************************


! includes self-energy for molecules with partial charges
! surface term can optionally be included
! spherical or cubic k-space summation (controlled by ncut2)

subroutine EwaldSetup

# ifdef F03_CBIND
   use EnergyModule, s=>meshsize
# else
   use EnergyModule
# endif
   implicit none

   character(40), parameter :: txroutine ='EwaldSetup'
   real(8)    :: fac0, facx, facy, facz, fac, k2, kp, EwaldErrorUReal
   real(8)    :: m2, c, spl(0:order-1), b2(0:nmesh-1,3), expf, mf
   integer(4) :: nx, ny, nz, ipt, ialoc, i, m, k, z, ierr
   complex(8) :: a
# ifdef F03_CBIND
   integer(C_SIZE_T) :: size_fftw_alloc
# endif

! ... determination of lq2sum and q2sum (for error analysis)

   q2sum = sum(zat(iatan(1:na))**2)/np
   if (q2sum > Zero) then
      lq2sum = 0                                   ! system contains charges
   else
      q2sum = Zero
      do ipt = 1, npt
         do ialoc = 1, napt(ipt)
            q2sum = q2sum + nppt(ipt)*sum(dipa(1:3,ialoc,ipt)**2)
         end do
      end do
      q2sum = q2sum/np
      if (q2sum > Zero) then
         lq2sum = 1                                         ! system contains no charges but dipoles
      else
         lq2sum = -1                                        ! system contains no charges and no dipoles
      end if
   end if
   if (lq2sum == -1) call stop(txroutine, 'no charges or permanent dipoles', uout)

   if (txewaldrec == 'std') then

! .............. reciprocal space: standard ..............

! ... set ewald parameters ualpha, rcut, and/or ncut
!
!     iewaldopt  uewaldtol  ualphared  ualpha  rcut  ncut
!     ---------  ---------  ---------  ------  ----  ----
!        0                      x                x     x    (all parametes taken from input)
!        1           d          x                x     d    (error analysis)
!        2           x                   x       d     d       "
!        3           x                   d       x     d       "
!        4           x                   d       d     x       "
!
!   x = parameter used; d = parameter determined (ualphared and ualpha are considered to be same parameter)
!
!   error analysis: coulomb system: kolafa and perram, molecular simulation 1992, 9 351
!                   dipolar system: wang and holm jcp 2001, 6351, 115

      if (iewaldopt == 0) then       ! all parameters taken from input  (ualphared, rcut, and ncut)
         if (ualphared <= zero) call Stop(txroutine, 'ualphared <= zero', 6)
         if (rcut <= zero) call Stop(txroutine, 'rcut <= zero', 6)
         if (ncut <= 0) call Stop(txroutine, 'ncut <= 0', 6)
         ualpha  = ualphared/rcut
      else if (iewaldopt == 1) then  ! (ulphared and rcut) -> (uewaldtol and ncut)
         if (ualphared <= zero) call Stop(txroutine, 'ualphared <= zero', 6)
         if (rcut <= zero) call Stop(txroutine, 'rcut <= zero', 6)
         ualpha  = ualphared/rcut
         uewaldtol =  EwaldErrorUReal(lq2sum,q2sum,boxlenshort,ualpha,rcut,EpsiFourPi)
         call AlphaToNcut()
      else if (iewaldopt == 2) then  ! (uewaldtol and ualpha) -> (rcut and ncut)
         if (uewaldtol <= zero) call Stop(txroutine, 'uewaldtol <= zero', 6)
         if (ualpha <= zero) call Stop(txroutine, 'ualpha <= zero', 6)
         call AlphaToRcut()
         call AlphaToNcut()
         ualphared = ualpha*rcut
      else if (iewaldopt == 3) then  ! (uewaldtol and rcut) -> (alpha and ncut)
         if (uewaldtol <= zero) call Stop(txroutine, 'uewaldtol <= zero', 6)
         if (rcut <= zero) call Stop(txroutine, 'rcut <= zero', 6)
         call RcutToAlpha()
         call AlphaToNcut()
         ualphared = ualpha*rcut
      else if (iewaldopt == 4) then  ! (uewaldtol and ncut) -> (alpha and rcut)
         if (uewaldtol <= zero) call Stop(txroutine, 'uewaldtol <= zero', 6)
         if (ncut <= 0) call Stop(txroutine, 'ncut <= 0', 6)
         call NcutToAlpha()
         call AlphaToRcut()
         ualphared = ualpha*rcut
      end if
      rcut2 = rcut**2

! ... allocate memory

      naewald = na
      if (ncut > int((huge(nkvec_word)/4.0d0)**(1.0d0/3.0d0)-1.0d0)) then
         call Stop(txroutine,'integer overflow expected for nkvec_word',uout)
      else
         nkvec = (ncut+1)**3   ! for cube
      end if
      nkvec_word = 4*nkvec     ! word length for communication with MPI

      if(allocated(kfac)) deallocate(kfac)
      if(allocated(eikx)) deallocate(eikx, eiky, eikz)
      if(allocated(eikyzm)) deallocate(eikyzm, eikyzp)
      if(allocated(eikr)) deallocate(eikr)
      if(allocated(sumeikr)) deallocate(sumeikr)
      if(allocated(sumeikrd)) deallocate(sumeikrd)

      if(lmc .or. lmcall) then
          if(allocated(eikxtm)) deallocate(eikxtm, eikytm, eikztm)
          if(allocated(eikyzmtm)) deallocate(eikyzmtm, eikyzptm)
          if(allocated(eikrtm)) deallocate(eikrtm)
          if(allocated(sumeikrtm)) deallocate(sumeikrtm)
      end if

      if(allocated(eikraux)) deallocate(eikraux)

      allocate(kfac(nkvec), stat = ierr)
      kfac = 0.0E+00
      allocate(eikx(naewald,0:ncut),eiky(naewald,0:ncut), eikz(naewald,0:ncut), stat = ierr)
      eikx = cmplx(Zero,Zero)
      eiky = cmplx(Zero,Zero)
      eikz = cmplx(Zero,Zero)
      allocate(eikyzm(naewald), eikyzp(naewald), stat = ierr)
      eikyzm = cmplx(Zero,Zero)
      eikyzp = cmplx(Zero,Zero)
      allocate(eikr(naewald,4), stat = ierr)
      eikr = cmplx(Zero,Zero)
      allocate(sumeikr(nkvec,4), stat = ierr)
      sumeikr = cmplx(Zero,Zero)
      allocate(sumeikrd(nkvec,4), stat = ierr)
      sumeikrd = cmplx(Zero,Zero)

      if(lmc .or. lmcall) then
          allocate(eikxtm(naewald,0:ncut),eikytm(naewald,0:ncut), eikztm(naewald,0:ncut), stat = ierr)
          eikxtm = cmplx(Zero,Zero)
          eikytm = cmplx(Zero,Zero)
          eikztm = cmplx(Zero,Zero)
          allocate(eikyzmtm(naewald), eikyzptm(naewald), stat = ierr)
          eikyzmtm = cmplx(Zero,Zero)
          eikyzptm = cmplx(Zero,Zero)
          allocate(eikrtm(naewald,4), stat = ierr)
          eikrtm = cmplx(Zero,Zero)
          allocate(sumeikrtm(nkvec,4), stat = ierr)
          sumeikrtm = cmplx(Zero,Zero)
      end if

      allocate(eikraux(nkvec_word), stat = ierr)
      eikraux = cmplx(Zero,Zero)

      if (ierr /= 0) call WriteIOStat(txroutine, 'memory allocation failed', ierr, 2, 6)

! ... determine ncut2

      if (ncutregion == 'sphere') then
         ncut2 = ncut**2
      else if (ncutregion == 'cube') then
         ncut2 = 3*ncut**2
      else
         call Stop(txroutine, 'error in ncutregion', uout)
      end if

! ... determine nkvec and kfac for the reciprocal space

      nkvec = 0
      fac0 = TwoPi/vol
      do nz = 0, ncut
         facz = Two
         if (nz == 0) facz = One
         do ny = 0, ncut
            facy = One
            if (ny == 0) facy = Half
            do nx = 0, ncut
               if ((lbcrd .or. lbcto) .and. (mod((nx+ny+nz),2) /= 0)) cycle ! only even nx+ny+nz for RD and TO bc
               if (nx**2+ny**2+nz**2 > ncut2) cycle
               if (nx == 0 .and. ny == 0 .and. nz == 0) cycle
               nkvec = nkvec+1
               facx = One
               if (nx == 0) facx = Half
               k2 = (nx*TwoPiBoxi(1))**2+(ny*TwoPiBoxi(2))**2+(nz*TwoPiBoxi(3))**2
               kfac(nkvec) = fac0*facz*facy*facx*exp(-k2/(Four*ualpha**2))/k2
            end do
         end do
      end do

      if (lewald2dlc) then

         ncut2d = ncut            ! later increase ncut2d
         nkvec2d = (ncut2d+1)**2  ! for cube

         if(allocated(kfac2d)) deallocate(kfac2d)
         if(allocated(sinkx)) deallocate(sinkx, coskx, sinky, cosky)
         if(allocated(sinkxtm)) deallocate(sinkxtm, coskxtm, sinkytm, coskytm)
         if(allocated(sumtrigtm)) deallocate(sumtrigtm)
         allocate(kfac2d(nkvec2d), stat = ierr)
         kfac2d = 0.0E+00
         allocate(sinkx(naewald,0:ncut2d), coskx(naewald,0:ncut2d), sinky(naewald,0:ncut2d), cosky(naewald,0:ncut2d), stat = ierr)
         sinkx = 0.0E+00
         coskx = 0.0E+00
         sinky = 0.0E+00
         cosky = 0.0E+00
         allocate(sumtrig(nkvec2d,9), stat = ierr)
         sumtrig = 0.0E+00
         allocate(termsss(nkvec2d), termcss(nkvec2d), termscs(nkvec2d), termccs(nkvec2d), &
                  termssc(nkvec2d), termcsc(nkvec2d), termscc(nkvec2d), termccc(nkvec2d), stat = ierr)
         termsss = 0.0E+00
         termcss = 0.0E+00
         termscs = 0.0E+00
         termccs = 0.0E+00
         termssc = 0.0E+00
         termcsc = 0.0E+00
         termscc = 0.0E+00
         termccc = 0.0E+00
         allocate(sinkxtm(naewald,0:ncut2d),&
            coskxtm(naewald,0:ncut2d),&
            sinkytm(naewald,0:ncut2d),&
            coskytm(naewald,0:ncut2d),&
            stat = ierr)
         sinkxtm = 0.0E+00
         coskxtm = 0.0E+00
         sinkytm = 0.0E+00
         coskytm = 0.0E+00
         allocate(sumtrigtm(nkvec2d,9), stat = ierr)
         sumtrigtm = 0.0E+00
         if (ierr /= 0) call WriteIOStat(txroutine, 'memory allocation failed', ierr, 2, 6)

! ... determine ncut2d

         if (ncutregion == 'sphere') then
            ncut2d2 = ncut2d**2
         else if (ncutregion == 'cube') then
            ncut2d2 = 2*ncut2d**2
         else
            call Stop(txroutine, 'error in ncutregion, 2d-correction', uout)
         end if

! ... determine nkvec2d and kfac2d

         nkvec2d = 0
!        fac0 = TwoPi*EpsiFourPi/(boxlen(1)*boxlen(2))    ! TO BE REMOVED
         fac0 = TwoPi/(boxlen(1)*boxlen(2))               ! problems with holm's equation
         do ny = 0, ncut2d
            do nx = 0, ncut2d
               if (nx**2+ny**2 > ncut2d2) cycle
               if (nx == 0 .and. ny == 0) cycle
               nkvec2d = nkvec2d+1
               fac = One
               if (nx/= 0 .and. ny/= 0) fac = Four*fac
               kp = sqrt((nx*TwoPiBoxi(1))**2+(ny*TwoPiBoxi(2))**2)
               kfac2d(nkvec2d) =-fac0*fac/(kp*(exp(kp*boxlen(3))-One))
            end do
         end do

      end if

   else if (txewaldrec == 'spm') then

# ifndef F03_CBIND
      call stop(txroutine, 'SPME is selected but the fftw3 package was not linked into the executable', uout)
# else

! .............. reciprocal space: smooth particle mesh ..............

! ... set ewald parameters ualpha, rcut, and/or nmesh and order
!
!     iewaldopt  uewaldtol  ualphared  ualpha  rcut  nmesh  order
!     ---------  ---------  ---------  ------  ----  -----  -----
!        0                      x                x     x      x     (all parametes taken from input)
!        3           x          d                x     x      x     (error analysis)
!
!   x = parameter used; d = parameter determined (ualphared and ualpha are considered to be same parameter)
!
!   error analysis: coulomb system: kolafa and perram, molecular simulation 1992, 9 351
!                   dipolar system: wang and holm jcp 2001, 6351, 115

      if (iewaldopt == 0) then       ! all parameters taken from input  (ualphared, rcut, nmesh, and order)
         if (ualphared <= Zero) call Stop(txroutine, 'ualphared <= 0', uout)
         if (rcut <= zero) call Stop(txroutine, 'rcut <= zero', 6)
         if (nmesh <= 0) call Stop(txroutine,'nmesh <= 0',uout)
         if (order <= 0) call Stop(txroutine,'order <= 0',uout)
         ualpha  = ualphared/rcut
      else if (iewaldopt == 3) then  ! (uewaldtol and rcut) -> (alpha)
         if (uewaldtol <= zero) call Stop(txroutine, 'uewaldtol <= zero', 6)
         if (rcut <= zero) call Stop(txroutine, 'rcut <= zero', 6)
         call RcutToAlpha()
         if (nmesh <= 0) call Stop(txroutine,'nmesh <= 0',uout)
         if (order <= 0) call Stop(txroutine,'order <= 0',uout)
         ualphared = ualpha*rcut
      else
         call Stop(txroutine, 'iewaldopt out of order', uout)
      end if
      if (mod(nmesh,2) /= 0 ) call Stop(txroutine,'nmesh must be even', uout)
      s(1:3) = nmesh ! FIXME non-cubic mesh

! ... allocate memory

      if (associated(QMesh)) then
         call fftw_free(qmesh_ptr)
         Qmesh => null()
      endif
      size_fftw_alloc = s(1) * s(2) * int ( s(3), C_SIZE_T)
      qmesh_ptr = fftw_alloc_real( size_fftw_alloc )
      call c_f_pointer(qmesh_ptr, Qmesh, [ s(1), s(2), s(3) ])

      if (associated(FQmesh)) then
         call fftw_free(fqmesh_ptr)
         FQmesh => null()
      endif
      fqmesh_ptr = fftw_alloc_complex( int( (s(1)/2+1)*s(2)*s(3), C_SIZE_T ))
      call c_f_pointer(fqmesh_ptr, FQmesh, [ s(1)/2+1, s(2), s(3) ])
                                                                           ! Memory
      if (allocated(meshfac)) deallocate(meshfac,energyfac,virfac)
      allocate(meshfac(0:s(1)/2,0:s(2)/2,0:s(3)/2))                        ! 8 M / 8
      meshfac = 0.0E+00
      allocate(energyfac(0:s(1)/2,0:s(2)/2,0:s(3)/2))                      ! 8 M / 8
      energyfac = 0.0E+00
      allocate(virfac(0:s(1)/2,0:s(2)/2,0:s(3)/2))                         ! 8 M / 8
      virfac = 0.0E+00
      if (lmc) then
         if (associated(QMeshTM)) then
             call fftw_free(QmeshTM_PTR)
             QmeshTM => null()
         endif
      QmeshTM_ptr = fftw_alloc_real( int( s(1)*s(2)*s(3), C_SIZE_T ) )
      call c_f_pointer(QmeshTM_ptr, QmeshTM, [ s(1), s(2), s(3) ])
      end if

! ... creating FFTW plans, we must reverse

      if (plan_fwd_done) call fftw_destroy_plan(plan_fwd)
      if (plan_bwd_done) call fftw_destroy_plan(plan_bwd)
      plan_fwd = fftw_plan_dft_r2c_3d ( s(3), s(2), s(1), QMesh, FQMesh, plan_methode )
      plan_bwd = fftw_plan_dft_c2r_3d ( s(3), s(2), s(1), FQMesh, QMesh, plan_methode )
      plan_fwd_done = .true.
      plan_bwd_done = .true.

      MeshMemSize = s(1) * s(2) * s(3)

#if defined (_PAR_)
      if (allocated(meshaux)) deallocate(meshaux)
      allocate(meshaux(MeshMemSize))
      meshaux = 0.0E+00
#endif

      if (allocated(spline)) deallocate(spline, deriv, deriv2, meshmax)
      allocate(spline(0:order-1,3,na_alloc))                               ! 8 * N * O**3
      spline = 0.0E+00
      allocate(deriv(0:order-1,3,na_alloc))                                ! 8 * N * O**3
      deriv = 0.0E+00
      allocate(deriv2(0:order-1,3,na_alloc))                               ! 8 * N * O**3
      deriv2 = 0.0E+00
      allocate(meshmax(3,na_alloc))                                        ! 8 * 3 * N
      meshmax = 0
      if (lmc) then
         natm = na ! FIXME
         if (allocated(splinetm)) deallocate(splinetm, derivtm, deriv2tm, meshmaxtm)
         allocate(splinetm(0:order-1,3,natm))                              ! 8 * N * O**3
         splinetm = 0.0E+00
         allocate(derivtm(0:order-1,3,natm))                               ! 8 * N * O**3
         derivtm = 0.0E+00
         allocate(deriv2tm(0:order-1,3,natm))                              ! 8 * N * O**3
         deriv2tm = 0.0E+00
         allocate(meshmaxtm(3,natm))                                       ! 8 * 3 * N
         meshmaxtm = 0
      end if
!                                                                        --------------------
!                                                                        28 M + 48 NO**3 + 48 N
!
!                                                     assume M = N and O = 6 => 10N kB

!                                                     CPU time : Make Mesh Order(NO**3)
!                                                                Make FFT  Order(M ln M)
!                                                                Calc GIF  Order(M)
!                                                                Calc Prop Order(NO**3)

! ... asign mesh units per length units

      dkdr(1:3) = s(1:3)/boxlen(1:3)
      d2kdr(1:6) = [dkdr(1)**2, dkdr(2)**2, dkdr(3)**2, dkdr(1)*dkdr(2), dkdr(1)*dkdr(3), dkdr(2)*dkdr(3)]

! ... calculate b(m)**2, eq (4.4) of JCP 1995, 103, 8577

      call CardinalBSpline(order, Zero, z, spl)
      do i = 1,3
         do m = 0,s(i)-1
            if ((mod(order,2) == 1) .and. (2*m==s(i))) cycle
            a = Zero
            do k = 0,order-2
               a=a+spl(k+1)*exp(cmplx(0,TwoPi*m*k/s(i)))
            end do
            b2(m,i) = One/(real(a)**2+aimag(a)**2)
         end do
      end do

! ... calculate B*C, eqs (3.9) and (4.8) of JCP 1995, 103, 8577

      fac0 = FourPi/vol
      expf = -One/(Four*ualpha**2)
      do nz = 0,s(3)/2
         do ny = 0,s(2)/2
            do nx = 0,s(1)/2
               if ((nx == 0) .and. (ny == 0) .and. (nz == 0)) then
                  meshfac(nx,ny,nz) = Zero
                  energyfac(nx,ny,nz) = Zero
                  virfac(nx,ny,nz) = Zero
               else
                  m2 = (nx*TwoPiBoxi(1))**2+(ny*TwoPiBoxi(2))**2+(nz*TwoPiBoxi(3))**2
                  c = exp(m2*expf)/m2
                  mf = fac0*c*b2(nx,1)*b2(ny,2)*b2(nz,3)
                  meshfac(nx,ny,nz) = mf
                  if (nx == mod(s(1)-nx,s(1))) mf=mf/2   ! compensate for double counts in CalcGIF section
                  if (ny == mod(s(2)-ny,s(2))) mf=mf/2   ! compensate for double counts in CalcGIF section
                  if (nz == mod(s(3)-nz,s(3))) mf=mf/2   ! compensate for double counts in CalcGIF section
                  energyfac(nx,ny,nz) = mf
                  virfac(nx,ny,nz) = (One + 2*m2*expf)*mf
               end if
            end do
         end do
      end do
# endif
   end if

! ... check upper limit of rcut

   if (lbcbox .and. (rcut > minval(boxlen2(1:3)))) then
      write(uout,'(a,1g15.5)') 'rcut = ', rcut
      write(uout,'(a,3g15.5)') 'boxlen2 = ', boxlen2(1:3)
      call Warn(txroutine, 'lewald .and. (rcut > minval(boxlen2(1:3)))', uout)
   else if (lbcrd .and. (rcut > sqrt(Two/Three)*cellside)) then
      write(uout,'(a,1g15.5)') 'rcut = ', rcut
      write(uout,'(a,3g15.5)') 'sqrt(2/3)*cellside', sqrt(Two/Three)*cellside
      call Warn(txroutine, 'lewald .and. (rcut > sqrt(2/3)*cellside)', uout)
   else if (lbcto .and. (rcut > sqrt(Three/Two)*cellside)) then
      write(uout,'(a,1g15.5)') 'rcut = ', rcut
      write(uout,'(a,3g15.5)') 'sqrt(3/2)*cellside', sqrt(Three/Two)*cellside
      call Warn(txroutine, 'lewald .and. (rcut > sqrt(3/2)*cellside)', uout)
   end if

! ... set other numerical constants for Ewald summation

   ualpha2 = ualpha**2               ! ualpha is now available for all options of txewaldrec
   fac = One /(sqrt(Pi)*ualpha)
   ewaldfac1 = fac*Two*ualpha2
   ewaldfac2 = fac*(Two*ualpha2)**2/Three
   ewaldfac3 = fac*(Two*ualpha2)**3/(Three*Five)
   ewaldselffac1 = (Two*ualpha)/sqrt(Pi)
   ewaldselffac2 = (Four*ualpha**3)/(Three*sqrt(Pi))

contains

!........................................................................

subroutine AlphaToRcut()
   character(40), parameter :: txroutine ='AlphaToRcut'
   integer(4) :: i
   real(8) :: y, xold, xnew
   if (lq2sum == 0) y = (uewaldtol/(q2sum*EpsiFourPi))*sqrt(Two*boxlenshort**3)*ualpha**2
   if (lq2sum == 1) y = (sqrt(15.0d0)/4.0d0)*(uewaldtol/(q2sum*EpsiFourPi))*sqrt(boxlenshort**3)/ualpha**2
   xold = One
   do i = 1, 40
      if (lq2sum == 0) xnew = (One/ualpha)*sqrt(-log(y*xold**1.5))
      if (lq2sum == 1) xnew = (One/ualpha)*sqrt(-log(y/xold**0.5))
      if (abs(xold-xnew)/xnew < 1.0d-6) exit
      xold = xnew
   end do
   if (i > 40) call Stop(txroutine,'i > 40',6)
   rcut = xnew
end subroutine AlphaToRcut

!........................................................................

subroutine AlphaToNcut()
   character(40), parameter :: txroutine ='AlphaToNcut'
   integer(4) :: i
   real(8) :: y, xold, xnew
   if (lq2sum == 0) y = (uewaldtol/(q2sum*EpsiFourPi))*(Pi**2/(ualpha**2*boxlenshort))
   if (lq2sum == 1) y = (sqrt(15.0d0)/4.0d0)*(uewaldtol/(q2sum*EpsiFourPi))*boxlenshort/ualpha**2
   xold = One
   do i = 1, 40
      if (lq2sum == 0) xnew = (ualpha*boxlenshort/Pi)*sqrt(-log(y*xold/(One+One/(ualpha*boxlenshort*sqrt(xold)))))
      if (lq2sum == 1) xnew = (ualpha*boxlenshort/Pi)*sqrt(-log(y/xold/(One+One/(ualpha*boxlenshort*sqrt(xold)))))
      if (abs(xold-xnew)/xnew < 1.0d-6) exit
      xold = xnew
   end do
   if (i > 40) call Stop(txroutine,'i > 40',6)
   ncut = int(xnew + One)
end subroutine AlphaToNcut

!........................................................................

subroutine RcutToAlpha()
   character(40), parameter :: txroutine ='RCutToAlpha'
   integer(4) :: i
   real(8) :: y, xold, xnew
   if (lq2sum == 0) y = (uewaldtol/(q2sum*EpsiFourPi))*sqrt(Two*(boxlenshort*rcut)**3)
   if (lq2sum == 1) y = (3.0d0/4.0d0)*(uewaldtol/(q2sum*EpsiFourPi))*sqrt((boxlenshort**3/rcut))
   xold = 1.0d-1
   do i = 1, 40
      if (lq2sum == 0) xnew = (One/rcut)*sqrt(-log(y*xold**2))
      if (lq2sum == 1) xnew = (One/rcut)*sqrt(-log(y/xold**2))
      if (abs(xold-xnew)/xnew < 1.0d-6) exit
      xold = xnew
   end do
   if (i > 40) call Stop(txroutine,'i > 40',6)
   ualpha = xnew
end subroutine RcutToAlpha

!........................................................................

subroutine NcutToAlpha()
   character(40), parameter :: txroutine ='NcutToAlpha'
   integer(4) :: i
   real(8) :: y, xold, xnew
   if (lq2sum == 0) y = (uewaldtol/(q2sum*EpsiFourPi))*(Pi**2*ncut)/boxlenshort
   if (lq2sum == 1) y = (3.0d0/4.0d0)*(uewaldtol/(q2sum*EpsiFourPi))*boxlenshort/ncut
   xold = One
   do i = 1, 40
      if (lq2sum == 0) xnew = (Pi*ncut/boxlenshort)/sqrt(-log(y/xold**2/(1.0+1.0/(xold*boxlenshort*sqrt(real(ncut))))))
      if (lq2sum == 1) xnew = (Pi*ncut/boxlenshort)/sqrt(-log(y/xold**2/(1.0+1.0/(xold*boxlenshort*sqrt(real(ncut))))))
      if (abs(xold-xnew)/xnew < 1.0d-6) exit
      xold = xnew
   end do
   if (i > 40) call Stop(txroutine,'i > 40',6)
   ualpha = xnew
end subroutine NcutToAlpha

end subroutine EwaldSetup

!************************************************************************
!> \page energy energy.F90
!! **Getnkvec**
!! *return the value of nkvec*
!************************************************************************


function Getnkvec()
   use EnergyModule
   implicit none
   integer(4) :: Getnkvec
   Getnkvec = nkvec
end function Getnkvec

!************************************************************************
!> \page energy energy.F90
!! **Getnkvec2d**
!! *return the value of nkvec2d*
!************************************************************************


function Getnkvec2d()
   use EnergyModule
   implicit none
   integer(4) :: Getnkvec2d
   Getnkvec2d = nkvec2d
end function Getnkvec2d

!************************************************************************
!> \page energy energy.F90
!! **Gettime_ewald**
!! *return the value of ncut2d*
!************************************************************************


function Gettime_ewald()
   use EnergyModule
   implicit none
   real(8) :: Gettime_ewald
   Gettime_ewald = time_ewald
end function Gettime_ewald

! ........................... Under Development ...............................

!************************************************************************
!> \page energy energy.F90
!! **UWeakChargeA**
!! *calculate two-body potential energy; only monoatomic particles*
!************************************************************************

!     weak charges, limited to charged hard spheres

subroutine UWeakChargeA

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UWeakChargeA'
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf, Getnpmyid
   real(8)    :: dx, dy, dz, r2, d, usum, fsum, virtwob

   if (.not.lmonoatom) call Stop(txroutine, '.not.lmonoatom', uout)

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      if (.not.laz(ip)) cycle  ! ip uncharged
      do jploc = 1, nneighpn(iploc)
         jp = jpnlist(jploc,iploc)
         if (.not.laz(jp)) cycle  ! jp uncharged
         if (lmc) then
            if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBCr2(dx,dy,dz,r2)
         if (r2 > rcut2) cycle
         if (r2 < r2umin(iptjpt)) call StopUWeakChargeA
            ibuf = iubuflow(iptjpt)
            do
               if (r2 >= ubuf(ibuf)) exit
               ibuf = ibuf+12
               if (ibuf > nbuf) call StopIbuf('txptpt',iptjpt)
            end do
            d = r2-ubuf(ibuf)
            usum = ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                   d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
            fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                   d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))

         u%twob(iptjpt) = u%twob(iptjpt) + usum
         force(1,ip) = force(1,ip) + (fsum * dx)
         force(2,ip) = force(2,ip) + (fsum * dy)
         force(3,ip) = force(3,ip) + (fsum * dz)
         force(1,jp) = force(1,jp) - (fsum * dx)
         force(2,jp) = force(2,jp) - (fsum * dy)
         force(3,jp) = force(3,jp) - (fsum * dz)
         virtwob     = virtwob     - (fsum * r2)

      end do

   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot     + u%twob(0)
   virial    = virial    + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUWeakChargeA
   character(40), parameter :: txroutine ='StopUWeakChargeA'
   write(uout,'(a,i5)')  'ip', ip
   write(uout,'(a,i5)')  'jp', jp
   write(uout,'(a,3e15.5)') 'ro(ip)        = ', ro(1:3,ip)
   write(uout,'(a,3e15.5)') 'ro(jp)        = ', ro(1:3,jp)
   write(uout,'(a,3e15.5)') 'boxlen2       = ', boxlen2
   write(uout,'(a,3e15.5)') 'dpbc          = ', dpbc
   write(uout,'(a,i15)')    'uout           = ', uout
   write(uout,'(a,i15)')    'myid           = ', myid
   write(uout,'(a,i15)')    'iptjpt         = ', iptjpt
   write(uout,'(a,e15.5)')  'r2             = ', r2
   write(uout,'(a,e15.5)')  'r2umin(iptjpt) = ', r2umin(iptjpt)
   call Stop(txroutine, 'r2 < r2umin(iptjpt)', uout)
end subroutine StopUWeakChargeA

!........................................................................

end subroutine UWeakChargeA

!************************************************************************
!> \page energy energy.F90
!! **UWeakChargeP**
!! *calculate two-body potential energy; general particles*
!************************************************************************

!     weak charges, limited to charged hard spheres

subroutine UWeakChargeP

   use EnergyModule
   implicit none

   character(40), parameter :: txroutine ='UWeakChargeP'
   logical, parameter :: lintrapartint = .true.                   ! enable intraparticle interaction
   integer(4) :: ip, iploc, ipt, jp, jploc, jpt, iptjpt, ibuf, Getnpmyid
   integer(4) :: ia, ialow, iaupp, iat, ja, jalow, jaupp, jat, iatjat
   real(8)    :: dx, dy, dz, dxopbc, dyopbc, dzopbc, r2, d, usum, fsum, virtwob

   if (ltime) call CpuAdd('start', txroutine, 2, uout)

   u%twob(0:nptpt) = Zero
   virtwob         = Zero

   do iploc = 1, Getnpmyid()
      ip = ipnploc(iploc)
      ipt = iptpn(ip)
      ialow = ianpn(ip)
      iaupp = ialow+napt(ipt)-1
      do jploc = 0, nneighpn(iploc)
         if (jploc == 0) then
            if (.not.lintrapartint) cycle
            jp = ip
         else
            jp = jpnlist(jploc,iploc)
         end if
         if (lmc) then
           if (jp < ip) cycle
         end if
         jpt = iptpn(jp)
         iptjpt = iptpt(ipt,jpt)
         dx = ro(1,ip)-ro(1,jp)
         dy = ro(2,ip)-ro(2,jp)
         dz = ro(3,ip)-ro(3,jp)
         call PBC2(dx,dy,dz,dxopbc,dyopbc,dzopbc)
         dx = dx-dxopbc
         dy = dy-dyopbc
         dz = dz-dzopbc
         r2 = dx**2+dy**2+dz**2
         if (r2 > rcut2) cycle

         usum = Zero
         jalow = ianpn(jp)
         jaupp = jalow+napt(jpt)-1
         do ia = ialow, iaupp
            if (.not.laz(ia)) cycle  ! ia uncharged
            iat = iatan(ia)
            do ja = jalow, jaupp
               if (ja <= ia) cycle
               if (.not.laz(ja)) cycle  ! ja uncharged
               jat = iatan(ja)
               iatjat = iatat(iat,jat)
               dx = r(1,ia)-r(1,ja)-dxopbc
               dy = r(2,ia)-r(2,ja)-dyopbc
               dz = r(3,ia)-r(3,ja)-dzopbc
               r2 = dx**2+dy**2+dz**2
               if (r2 < r2umin(iatjat)) call StopUWeakChargeP
               ibuf = iubuflow(iatjat)
               do
                  if (r2 >= ubuf(ibuf)) exit
                  ibuf = ibuf+12
                  if (ibuf > nbuf) call StopIbuf('txatat',iatjat)
               end do
               d = r2-ubuf(ibuf)
               usum = usum+ubuf(ibuf+1)+d*(ubuf(ibuf+2)+d*(ubuf(ibuf+3)+ &
                      d*(ubuf(ibuf+4)+d*(ubuf(ibuf+5)+d*ubuf(ibuf+6)))))
               fsum = ubuf(ibuf+7)+d*(ubuf(ibuf+8)+d*(ubuf(ibuf+9)+ &
                      d*(ubuf(ibuf+10)+d*ubuf(ibuf+11))))

               force(1,ia) = force(1,ia) + (fsum * dx)
               force(2,ia) = force(2,ia) + (fsum * dy)
               force(3,ia) = force(3,ia) + (fsum * dz)
               force(1,ja) = force(1,ja) - (fsum * dx)
               force(2,ja) = force(2,ja) - (fsum * dy)
               force(3,ja) = force(3,ja) - (fsum * dz)
               virtwob     = virtwob     - (fsum * r2)
    !    if (itest == 90) write(uout,'(a,2i5,2f12.5)') 'ia, ja, r2, usum', ia, ja, r2, usum   !cc
            end do
         end do
         u%twob(iptjpt) = u%twob(iptjpt) + usum
      end do

   end do

   u%twob(0) = sum(u%twob(1:nptpt))

   u%tot     = u%tot     + u%twob(0)
   virial    = virial    + virtwob

   if (ltime) call CpuAdd('stop', txroutine, 2, uout)

contains

!........................................................................

subroutine StopIbuf(txstring,i)
   character(*), intent(in) :: txstring
   integer(4),   intent(in) :: i
   write(uout,*)
   write(uout,'(a,i5)') txstring, i
   call Stop(txroutine, 'ibuf > nbuf', uout)
end subroutine StopIbuf

subroutine StopUWeakChargeP
   character(40), parameter :: txroutine ='StopUWeakChargeP'
   write(uout,'(a,i15)')   'uout          = ', uout
   write(uout,'(a,i15)')   'myid          = ', myid
   write(*,'(a,i15)')      'myid          = ', myid
   write(uout,'(a,i15)')   'iatjat        = ', iatjat
   write(uout,'(a,i15)')   'ia            = ', ia
   write(uout,'(a,i15)')   'ja            = ', ja
   write(uout,'(a,e15.5)') 'r2            = ', r2
   write(uout,'(a,e15.5)') 'r2umin(iatjat) = ', r2umin(iatjat)
   call Stop(txroutine, 'r2 < r2umin(iatjat)', uout)
end subroutine StopUWeakChargeP

!........................................................................

end subroutine UWeakChargeP

