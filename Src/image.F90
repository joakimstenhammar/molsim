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
!*     ImageDriver                                                      *
!*                                                                      *
!************************************************************************

! ... driver of image preparation routines

subroutine ImageDriver(iStage)

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='ImageDriver'
   integer(4), intent(in) :: iStage
   logical,    save :: lvrml, lvtf, limageuser
   integer(4), save :: iimage
   integer(4), save :: nsamp1, nsamp2

   namelist /nmlImage/ iimage, lvrml, lvtf, limageuser

   if (slave) return

   if (ltrace) call WriteTrace(1, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 0, uout)

   select case (iStage)
   case (iReadInput)

      iimage     = 1
      lvrml      = .false.
      lvtf       = .false.
      limageuser = .false.

      rewind(uin)
      read(uin,nmlImage)

! ... check conditions

! ... iimage may not be chosen to be 0

      if (iimage <= 0) call Stop(txroutine, 'iimage may not be chosen to be smaller or equal 0', uout)

! ... ImageVTF requires grouping

      if (lvtf .and. .not.lgroup) call Stop(txroutine, 'lvtf is selected, but no group division', uout)

      call ImageDriverSub(iimage)

   case (iWriteInput)

      call ImageDriverSub(iimage)

   case (iBeforeSimulation)

      nsamp1  = 0
      if (lsim .and. master .and. txstart == 'continue') read(ucnf) nsamp1

      call ImageDriverSub(iimage)

   case (iBeforeMacrostep)

      nsamp2  = 0
      call ImageDriverSub(iimage)

   case (iSimulationStep)

      if (mod(istep2,iimage) == 0) then
         nsamp2 = nsamp2+1
         call ImageDriverSub(iimage)
      end if

   case (iAfterMacrostep)

      nsamp1  = nsamp1+1
      if (lsim .and. master) write(ucnf) nsamp1
      call ImageDriverSub(iimage)

   case (iAfterSimulation)

      if (master) then
         call WriteHead(2, 'image preparation: general', uout)
         write(uout,'(a,t35,i10)') 'image interval                 = ', iimage
         write(uout,'(a,t35,i10)') 'no of time steps/passes used   = ', nsamp2*nsamp1
         write(uout,'()')
         write(uout,'(a)') 'image preparation routines used'
         write(uout,'(a)') '-------------------------------'
         if (lvrml     ) write(uout,'(a)') '   ImageVRML  '
         if (lvtf      ) write(uout,'(a)') '   ImageVTF   '
         if (limageuser) write(uout,'(a)') '   ImageUser  '
      end if
      call ImageDriverSub(iimage)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 0, uout)

contains

!........................................................................

subroutine ImageDriverSub(iimage)
   integer(4), intent(in) :: iimage
   if (lvrml      .and. master) call ImageVRML(iStage, iimage)
   if (lvtf       .and. master) call ImageVTF(iStage,iimage)
   if (limageuser .and. master) call ImageUser(iStage)
end subroutine ImageDriverSub

!........................................................................

end subroutine ImageDriver

!************************************************************************
!*                                                                      *
!*     ImageVRML                                                        *
!*                                                                      *
!************************************************************************

! ... generate input files for vrml 97 viewers

subroutine ImageVRML(iStage, iimage)

   use MolModule
   implicit none

   integer(4), intent(in) :: iStage
   integer(4), intent(in) :: iimage

   character(40), parameter :: txroutine ='ImageVRML'
   character(80), parameter :: txheading ='preparation of vrml file'
   character(9),  save :: txfile
   character(12), save :: txwhen
   real(8), allocatable, save :: atsize(:), rgbcolor(:,:)
   real(8),       save :: blmax, bondr
   character(20), save :: tximage(3)
   character(8),  save :: txwrap(2) =['C&_FILE=', 'C&_END  ' ]
   integer(4),    save :: iframe
   logical,       save :: lgr
   character(20) :: string
   integer(4)    :: iat, m, m2, nfac

   namelist /nmlVRML/ txfile, txwhen, atsize, rgbcolor, blmax, bondr, lgr, tximage

   select case (iStage)
   case (iReadInput)

      if (.not.allocated(atsize)) then
         allocate(atsize(nat), rgbcolor(3,nat))
      end if

! ... set default values

      txfile = 'merged'
      txwhen = 'after_run'
      atsize(1:nat) = radat(1:nat)
      tximage = ['frame','     ','     ']
      do iat = 1, 3
         if(iat > nat) exit
         rgbcolor(1:3,iat) = [ Zero,  Zero, Zero ]
         rgbcolor(iat,iat) = One
      end do
      do iat = 4, 6
         if(iat > nat) exit
         rgbcolor(1:3,iat) = [ One,  One, One]
         rgbcolor(iat-3,iat) = Zero
      end do
      do iat = 7, nat
        rgbcolor(1:3,iat) = [ (One/(iat+m), m = 1,3) ]
      end do

      blmax = 1.5d0
      bondr = 0.3d0
      lgr   =.false.

! ... read input data

      rewind(uin)
      read(uin,nmlVRML)

      call LowerCase(txwhen)

      if (tximage(2)(5:9) == 'plane') tximage(2)(5:11) = 'surface'   ! backward compability

   case (iWriteInput)

! ... open FWRL

      if (master) call FileOpen(uwrl, fwrl, 'form/noread')

   case (iBeforeSimulation)

      iframe = -1
      if (lsim .and. master) then
         if (txwhen == 'after_run') then
            continue                                                              ! do nothing
         else if ((txwhen == 'after_macro') .or. (txwhen == 'after_iimage')) then
            if ((txstart == 'setconf') .or. (txstart == 'zero')) then             ! write start frame
               call ImageVRMLSub
            else if(txstart == 'continue') then                                   ! forward ucnf
               read(ucnf) iframe
               if (txwhen == 'after_macro') nfac = 1
               if (txwhen == 'after_iimage') nfac = nstep2/iimage
               do m = 0, nfac*(nstep1beg-1)                                       ! number of frames to be read
                  do m2 = 1, 10000
                     read(uwrl,'(a)') string
                     if (string(1:6) == txwrap(2)) exit
                  end do
                  if (m2 > 10000) call stop(txroutine, 'error in advancing '//fwrl, uout)
               end do
               if (m > nfac*(nstep1beg-1)+1) call stop(txroutine, 'error in advancing '//fwrl, uout)
            end if
         end if
      end if

      if (lana .and. (txwhen == 'after_iimage')) then
         call IOCnf('open')
         call IOCnf('read')
         call IOCnf('close')
         call SetAtomPos(1, np, .false.)
         call ImageVRMLSub
      end if

      if (lsim .and. master .and. txstart == 'continue') then
      end if


   case (iSimulationStep)

      if (txwhen == 'after_iimage') then
         call ImageVRMLSub
      end if

   case (iAfterMacrostep)

      if (txwhen == 'after_macro') then
         call ImageVRMLSub
      end if

      if (lsim .and. master) write(ucnf) iframe

   case (iAfterSimulation)

      if (txwhen == 'after_run') then
         call ImageVRMLSub
      end if

      call WriteHead(2, txheading, uout)
      write(uout,'(a,a   )') 'vrml files are                   ', txfile
      write(uout,'(a,a   )') 'making vrml files                ', txwhen
      write(uout,'(a,2a  )') 'options (tximage)                ', tximage
      write(uout,'(a,f8.2)') 'max bond length                = ', blmax
      write(uout,'(a,f8.2)') 'bond radius                    = ', bondr
      write(uout,'()')
      write(uout,'(a,t15,a,t30,a,t50,a)') 'atom type no', 'atom type', 'size (radius)', 'rgbcolor'
      write(uout,'(a,t15,a,t30,a,t50,a)') '------------', '---------', '-------------', '--------'
      write(uout,'(i3,t15,a,t30,f8.3,t45,3f6.2)') (iat,txat(iat),atsize(iat),rgbcolor(1:3,iat),iat = 1,nat)
      write(uout,'()')
      write(uout,'(a,i4)')   'number of images made          = ', iframe+1

      if (allocated(atsize)) deallocate(atsize, rgbcolor)

      close(uwrl)

   end select

contains

! ... call of VRMLSUB and writing wrapper for mixed.f90 (if txfile == 'merged')

subroutine ImageVRMLSub

   character(40), parameter :: txroutine ='ImageVRMLSub'
   character(10) :: str
   integer(4)    :: unit

   iframe = iframe+1
   str = '      '
   write(str,'(i10)') iframe                              ! get frame number
   if (txfile == 'merged') then                           ! a single vrml file
      if (txwhen == 'after_run') then
         call VRMLSub(tximage, 'image after run', atsize, rgbcolor, blmax, bondr, lgr, uwrl)
      else if (txwhen == 'after_macro') then
         write(uwrl,'(a)') txwrap(1)//'macrostep.'//trim(adjustl(str))//'.wrl'
         call VRMLSub(tximage, 'image after macrostep'//str, atsize, rgbcolor, blmax, bondr, lgr, uwrl)
         write(uwrl,'(a)') txwrap(2)
      else if (txwhen == 'after_iimage') then
         write(uwrl,'(a)') txwrap(1)//'frame.'//trim(adjustl(str))//'.wrl'
         call VRMLSub(tximage, 'image after frame'//str, atsize, rgbcolor, blmax, bondr, lgr, uwrl)
         write(uwrl,'(a)') txwrap(2)
      else
         call Stop(txroutine, 'unsupported value of txwhen', uout)
      end if
   else if (txfile == 'separated') then                   ! separated vrml files
      unit = 55
      if (txwhen == 'after_run') then
         open(unit, file='run.wrl')
         call VRMLSub(tximage, 'image after run', atsize, rgbcolor, blmax, bondr, lgr, unit)
      else if (txwhen == 'after_macro') then
         open(unit, file='macrostep.'//trim(adjustl(str))//'.wrl')
         call VRMLSub(tximage, 'image after macrostep'//str, atsize, rgbcolor, blmax, bondr, lgr, unit)
      else if (txwhen == 'after_iimage') then
         open(unit, file='image.'//trim(adjustl(str))//'.wrl')
         call VRMLSub(tximage, 'image after frame'//str, atsize, rgbcolor, blmax, bondr, lgr, unit)
      else
         call Stop(txroutine, 'unsupported value of txwhen', uout)
      end if
      close(55)
   else
      call Stop(txroutine, 'unsupported value of txfile', uout)
   end if

end subroutine ImageVRMLSub

end subroutine ImageVRML

!************************************************************************
!*                                                                      *
!*     VRMLSub                                                          *
!*                                                                      *
!************************************************************************

! ... writes a vrml-97 file

subroutine VRMLSub(tximage, txlabel, atsize, rgbcolor, blmax, bondr, lgr, unit)

   use MolModule
   implicit none

   character(*), intent(in) :: tximage(*)       ! tximage(1) == ''              : nothing;
                                                ! tximage(1) == 'frame'         : draw a frame around the system
                                                ! tximage(2) == ''              : nothing;
                                                ! tximage(2) == 'one_surface'   : draw one surface at -boxlen(3)/2;
                                                ! tximage(2) == 'two_surface'   : draw surfaces at -boxlen(3)/2 and boxlen(3)/2
                                                ! tximage(3) == ''              : nothing;
                                                ! tximage(3) == 'undopbc'       : periodic boundray conditions are not applied to chains
   character(*), intent(in) :: txlabel          ! label
   real(8),      intent(in) :: atsize(*)        ! size of atom type
   real(8),      intent(in) :: rgbcolor(3,*)    ! rgb color components of each atom type
   real(8),      intent(in) :: blmax            ! maximal bond length
   real(8),      intent(in) :: bondr            ! bond radius
   logical,      intent(in) :: lgr              ! if true then coloring is following the group assignment
   integer(4),   intent(in) :: unit             ! output unit

   integer(4),    save :: mnbond, nbond              ! number of bonds
   integer(4), allocatable, save :: bondlist(:,:)    ! pair of atoms joind by a bond

   character(40), parameter :: txroutine ='VRMLSub'
   character(20) :: txcolor
   integer(4), allocatable :: icount(:)
   integer(4) :: ia, iat, ja, ib, icolor, iangle, nangle, i, ip, ipt, icorner, t(3)
   real(8) :: xdir, ydir, zdir, xnorm, ynorm, znorm, xc, yc, zc, height, arg, dangle, rrr(3), mat(4,4), dir
   real(8)     ,  parameter :: cornerref(3,8) = reshape( &
           [-One, One, One,   One, One, One,   One, One,-One,  -One, One,-One, &
            -One,-One, One,   One,-One, One,   One,-One,-One,  -One,-One,-One ] , [3,8] )
   real(8) :: corner(1:3,1:8)
   real(8), save :: rgbcolor_dipole(1:3,1:2) = reshape([One, One, Zero,  Zero, One, One],[3,2])

! ... initializing label

   write(unit,'(a)') '#VRML V2.0 utf8'
   write(unit,'(a)') '#-------------------------------------------'
   write(unit,'(a)') '#   start: '//txlabel
   write(unit,'(a)') '#-------------------------------------------'
   write(unit,'(a)')

! ... set view point

   write(unit,'(a)') '# ... viewpoint'
   if (lbcbox) then
      write(unit,'(a,3f12.3,a)') 'Viewpoint {description "View 1" position ', Two*boxlen(1), Zero, Zero
      write(unit,'(a)')         ' orientation 0 1 0 1.57079633 }'    ! viewing in the +x direction
   else if (lbcsph) then
      write(unit,'(a,3f12.3,a)') 'Viewpoint {description "View 1" position ', Zero, Zero, Three*sphrad,'}'
   else if (lbccyl) then
      write(unit,'(a,3f12.3,a)') 'Viewpoint {description "View 1" position ', Zero, Zero, 1.5*cyllen,'}'
   end if
   write(unit,'(a)')

! ... draw frame and surfaces

   if (tximage(1) == 'frame') then
      if (lbcbox) then
         call DrawBox
         if (tximage(2)(5:11) == 'surface') call DrawSurface
      else if (lbcrd) then
         call DrawRhombicDodecahedron
      else if (lbcto) then
         call DrawTruncatedOctahedron
      else if (lbcsph) then
         call DrawSphere
      else if (lbccyl) then
         call DrawCylinder(cylrad,cyllen)
         if (luext) then
            if (txuext(1) == 'hard_cylinder') call DrawCylinder(rCylinder, cyllen)
         end if
      end if
   end if

! ... undo periodic boundary conditions for chains and hierarchical structures

   vaux(1:3,1:na) = r(1:3,1:na)
   if (tximage(3) == 'undopbc') then
      call UndoPBC(vaux)
   end if

! ... make bonds

   mnbond = 2*na                                          ! maximal number of bonds that can be made
   if(.not.allocated(bondlist)) then
      allocate(bondlist(2,mnbond))
      bondlist = 0
   end if
   call MakeBondList(vaux, mnbond, blmax, nbond, bondlist)

! ... draw objects

   if (lellipsoid) then
      call DrawElliposid
   else if (lsuperball) then
      call DrawSuperball
   else
      if (txuser(1:6) == 'jasper') then
         call DrawAtom_jasper
      else
         call DrawAtom
      end if
   end if

! ... draw bonds

   if (nbond > 1) then
      call DrawBond
   end if

! ... finalizing label

   write(unit,'(a)') '#-------------------------------------------'
   write(unit,'(a)') '#   end: '//txlabel
   write(unit,'(a)') '#-------------------------------------------'

contains

!........................................................................

subroutine DrawBox
   corner(1,1:8) = boxlen2(1)*cornerref(1,1:8)
   corner(2,1:8) = boxlen2(2)*cornerref(2,1:8)
   corner(3,1:8) = boxlen2(3)*cornerref(3,1:8)
   write(unit,'(a)')          '# ... frame: box'
   write(unit,'(a)')          'Shape {appearance Appearance { material Material {'
   write(unit,'(a)')          ' diffuseColor 1.0 1.0 1.0 transparency 0.8 }}'
   write(unit,'(a,3f12.3,a)') ' geometry Box { size', boxlen(1), boxlen(2), boxlen(3), '}}'
   write(unit,'(a)')          'Shape {appearance Appearance {material Material {'
   write(unit,'(a)')          ' emissiveColor 0.8 0.8 0.8 }}'
   write(unit,'(a)')          ' geometry IndexedLineSet {coord Coordinate {'
   write(unit,'(a)')          ' point ['
   write(unit,'(3f12.3)')     (corner(1:3,icorner), icorner = 1, 8)
   write(unit,'(a)')          ' ]}'
   write(unit,'(a)')          ' coordIndex ['
   write(unit,'(a)')          ' 0, 1, 2, 3, 0, -1'
   write(unit,'(a)')          ' 4, 5, 6, 7, 4, -1'
   write(unit,'(a)')          ' 0, 4, -1'
   write(unit,'(a)')          ' 1, 5, -1'
   write(unit,'(a)')          ' 2, 6, -1'
   write(unit,'(a)')          ' 3, 7, -1'
   write(unit,'(a)')          ' ]}}'
   write(unit,'(a)')
end subroutine DrawBox

!........................................................................

subroutine DrawRhombicDodecahedron
   real(8) :: s, d
   d = sqrt(Two/Three)*cellside
   s = d/SqTwo
   write(unit,'(a)')        '# ... frame: rhombic dodecahedron'
   write(unit,'(a)')        'Shape {appearance Appearance {material Material {'
   write(unit,'(a)')        ' emissiveColor 0.8 0.8 0.8 }}'
   write(unit,'(a)')        ' geometry IndexedLineSet {coord Coordinate {'
   write(unit,'(a)')        ' point ['
   write(unit,'(3f12.3)')   Zero,  Zero,  SqTwo*d   !0
   write(unit,'(3f12.3)')   Zero,  Zero, -SqTwo*d   !1
   write(unit,'(3f12.3)')   Zero,    -d,        s   !2
   write(unit,'(3f12.3)')      d,    -d,     Zero   !3
   write(unit,'(3f12.3)')   Zero,    -d,       -s   !4
   write(unit,'(3f12.3)')     -d,    -d,     Zero   !5
   write(unit,'(3f12.3)')   Zero,     d,        s   !6
   write(unit,'(3f12.3)')      d,     d,     Zero   !7
   write(unit,'(3f12.3)')   Zero,     d,       -s   !8
   write(unit,'(3f12.3)')     -d,     d,     Zero   !9
   write(unit,'(3f12.3)')     -d,  Zero,        s   !10
   write(unit,'(3f12.3)')     -d,  Zero,       -s   !11
   write(unit,'(3f12.3)')      d,  Zero,        s   !12
   write(unit,'(3f12.3)')      d,  Zero,       -s   !13
   write(unit,'(a)')          ' ]}'
   write(unit,'(a)')          ' coordIndex ['
   write(unit,'(a)')          ' 2, 3, 4, 5, 2, -1'
   write(unit,'(a)')          ' 6, 7, 8, 9, 6, -1'
   write(unit,'(a)')          ' 10, 9, 11, 5, 10, -1'
   write(unit,'(a)')          ' 12, 7, 13, 3, 12, -1'
   write(unit,'(a)')          ' 0, 2, -1'
   write(unit,'(a)')          ' 0, 6, -1'
   write(unit,'(a)')          ' 0, 10, -1'
   write(unit,'(a)')          ' 0, 12, -1'
   write(unit,'(a)')          ' 1, 4, -1'
   write(unit,'(a)')          ' 1, 8, -1'
   write(unit,'(a)')          ' 1, 11, -1'
   write(unit,'(a)')          ' 1, 13, -1'
   write(unit,'(a)')          ' ]}}'
   write(unit,'(a)')
end subroutine DrawRhombicDodecahedron

!........................................................................

subroutine DrawTruncatedOctahedron
   real(8) :: s, d
   s = cellside/SqTwo
   d = SqTwo*cellside
   write(unit,'(a)')        '# ... frame: truncated octahedron'
   write(unit,'(a)')        'Shape {appearance Appearance {material Material {'
   write(unit,'(a)')        ' emissiveColor 0.8 0.8 0.8 }}'
   write(unit,'(a)')        ' geometry IndexedLineSet {coord Coordinate {'
   write(unit,'(a)')        ' point ['
   write(unit,'(3f12.3)')      s,  Zero,     d   !0
   write(unit,'(3f12.3)')   Zero,     s,     d   !1
   write(unit,'(3f12.3)')     -s,  Zero,     d   !2
   write(unit,'(3f12.3)')   Zero,    -s,     d   !3
   write(unit,'(3f12.3)')      s,  Zero,    -d   !4
   write(unit,'(3f12.3)')   Zero,     s,    -d   !5
   write(unit,'(3f12.3)')     -s,  Zero,    -d   !6
   write(unit,'(3f12.3)')   Zero,    -s,    -d   !7
   write(unit,'(3f12.3)')      d,     s,  Zero   !8
   write(unit,'(3f12.3)')      d,  Zero,     s   !9
   write(unit,'(3f12.3)')      d,    -s,  Zero   !10
   write(unit,'(3f12.3)')      d,  Zero,    -s   !11
   write(unit,'(3f12.3)')     -d,     s,  Zero   !12
   write(unit,'(3f12.3)')     -d,  Zero,     s   !13
   write(unit,'(3f12.3)')     -d,    -s,  Zero   !14
   write(unit,'(3f12.3)')     -d,  Zero,    -s   !15
   write(unit,'(3f12.3)')      s,     d,  Zero   !16
   write(unit,'(3f12.3)')   Zero,     d,     s   !17
   write(unit,'(3f12.3)')     -s,     d,  Zero   !18
   write(unit,'(3f12.3)')   Zero,     d,    -s   !19
   write(unit,'(3f12.3)')      s,    -d,  Zero   !20
   write(unit,'(3f12.3)')   Zero,    -d,     s   !21
   write(unit,'(3f12.3)')     -s,    -d,  Zero   !22
   write(unit,'(3f12.3)')   Zero,    -d,    -s   !23
   write(unit,'(a)')          ' ]}'
   write(unit,'(a)')          ' coordIndex ['
   write(unit,'(a)')          ' 0, 1, 2, 3, 0, -1'
   write(unit,'(a)')          ' 4, 5, 6, 7, 4, -1'
   write(unit,'(a)')          ' 8, 9, 10, 11, 8, -1'
   write(unit,'(a)')          ' 12, 13, 14, 15, 12, -1'
   write(unit,'(a)')          ' 16, 17, 18, 19, 16, -1'
   write(unit,'(a)')          ' 20, 21, 22, 23, 20, -1'
   write(unit,'(a)')          ' 0, 9, -1'
   write(unit,'(a)')          ' 1, 17, -1'
   write(unit,'(a)')          ' 2, 13, -1'
   write(unit,'(a)')          ' 3, 21, -1'
   write(unit,'(a)')          ' 4, 11, -1'
   write(unit,'(a)')          ' 5, 19, -1'
   write(unit,'(a)')          ' 6, 15, -1'
   write(unit,'(a)')          ' 7, 23, -1'
   write(unit,'(a)')          ' 8, 16, -1'
   write(unit,'(a)')          ' 10, 20, -1'
   write(unit,'(a)')          ' 12, 18, -1'
   write(unit,'(a)')          ' 14, 22, -1'
   write(unit,'(a)')          ' ]}}'
   write(unit,'(a)')
end subroutine DrawTruncatedOctahedron

!........................................................................

subroutine DrawSphere
   write(unit,'(a)')          '# ... frame: sphere'
   write(unit,'(a)')          'Shape {appearance Appearance { material Material {'
   write(unit,'(a)')          ' diffuseColor 1.0 1.0 1.0 transparency 0.8 }}'
   write(unit,'(a,1f12.3,a)') ' geometry Sphere { radius', sphrad, '}}'
   do i = 1, 3
      write(unit,'(a)')       'Shape {appearance Appearance {material Material {'
      write(unit,'(a)')       ' emissiveColor 0.8 0.8 0.8 }}'
      write(unit,'(a)')       ' geometry IndexedLineSet {coord Coordinate {'
      write(unit,'(a)')       ' point ['
      nangle = 128
      dangle = Two*pi/nangle
      do iangle = 1, nangle
         rrr(1+mod(i+0,3)) = sphrad*cos(iangle*dangle)
         rrr(1+mod(i+1,3)) = sphrad*sin(iangle*dangle)
         rrr(1+mod(i+2,3)) = Zero
         write(unit,'(3f12.3)') rrr(1:3)
      end do
      write(unit,'(a)')        ' ]}'
      write(unit,'(a)')        ' coordIndex ['
      write(unit,'(200(i3,a))') (iangle-1, ',', iangle = 1, nangle), (iangle-1, ',', iangle = 1, 0, -1)
      write(unit,'(a)')        ' ]}}'
   end do
   write(unit,'(a)')
end subroutine DrawSphere

!........................................................................

subroutine DrawCylinder(rad,len)
   real(8), intent(in) :: rad
   real(8), intent(in) :: len
   write(unit,'(a)')          '# ... frame: cylinder'
   write(unit,'(a,4(f10.2))') 'Transform  { rotation', One, Zero, Zero, Pi/Two
   write(unit,'(a,3(f10.2))') '          translation', Zero, Zero, Zero
   write(unit,'(a)')          ' children [ Shape{ appearance Appearance { material Material {'
   write(unit,'(a)')          ' diffuseColor 1.0 1.0 1.0 transparency 0.8 }}'
   write(unit,'(2(a,f10.2))') ' geometry  Cylinder { radius', rad, ' height', len
   write(unit,'(a)')          ' side TRUE   top FALSE  bottom FALSE } } ] }'
   do i = 1, 2
      write(unit,'(a)')       'Shape {appearance Appearance {material Material {'
      write(unit,'(a)')       ' emissiveColor 0.8 0.8 0.8 }}'
      write(unit,'(a)')       ' geometry IndexedLineSet {coord Coordinate {'
      write(unit,'(a)')       ' point ['
      nangle = 128
      dangle = Two*pi/nangle
      do iangle = 1, nangle
         rrr(1) = rad*cos(iangle*dangle)
         rrr(2) = rad*sin(iangle*dangle)
         rrr(3) = (-1)**i*Half*len
         write(unit,'(3f12.3)') rrr(1:3)
      end do
      write(unit,'(a)')        ' ]}'
      write(unit,'(a)')        ' coordIndex ['
      write(unit,'(200(i3,a))') (iangle-1, ',', iangle = 1, nangle), (iangle-1, ',', iangle = 1, 0, -1)
      write(unit,'(a)')        ' ]}}'
   end do
   write(unit,'(a)')          'Shape {appearance Appearance {material Material {'
   write(unit,'(a)')        ' emissiveColor 0.8 0.8 0.8 }}'
   write(unit,'(a)')          ' geometry IndexedLineSet {coord Coordinate {'
   write(unit,'(a)')          ' point ['
   write(unit,'(3f12.3)')     -rad, Zero, -Half*len
   write(unit,'(3f12.3)')     -rad, Zero, +Half*len
   write(unit,'(3f12.3)')     +rad, Zero, -Half*len
   write(unit,'(3f12.3)')     +rad, Zero, +Half*len
   write(unit,'(3f12.3)')     Zero, -rad, -Half*len
   write(unit,'(3f12.3)')     Zero, -rad, +Half*len
   write(unit,'(3f12.3)')     Zero, +rad, -Half*len
   write(unit,'(3f12.3)')     Zero, +rad, +Half*len
   write(unit,'(a)')          ' ]}'
   write(unit,'(a)')          ' coordIndex ['
   write(unit,'(a)')          '  0, 1, -1'
   write(unit,'(a)')          '  2, 3, -1'
   write(unit,'(a)')          '  4, 5, -1'
   write(unit,'(a)')          '  6, 7, -1'
   write(unit,'(a)')          ' ]}}'
   write(unit,'(a)')
end subroutine DrawCylinder

!........................................................................

subroutine DrawSurface
   real(8) :: thickness, zadjust
   thickness = 0.01d0*boxlen(3)
   zadjust = Zero
   write(unit,'(a)')          '# ... frame: surface'
   write(unit,'(a,3f12.3)')   'Transform { translation ', Zero, Zero, -(boxlen2(3)-zadjust)
   write(unit,'(a)')          'children [ DEF surface Shape{appearance Appearance { material Material {'
   write(unit,'(a)')          'diffuseColor 0.2 0.4 0.7 transparency 0.0 emissiveColor 0.2 0.4 0.7 shininess 1 }}'
   write(unit,'(a,3f12.3,a)') 'geometry Box { size', boxlen(1), boxlen(2), thickness, '}}]}'
   write(unit,'(a,3f12.3,a)') 'Transform { translation ', Zero, Zero, -(boxlen2(3)-zadjust), ' children [ USE surface ]}'
   if (tximage(2)(1:11) == 'one_surface') then
   else if (tximage(2)(1:11) == 'two_surface') then
      write(unit,'(a,3f12.3,a)') 'Transform { translation ', Zero, Zero, (boxlen2(3)-zadjust), ' children [ USE surface ]}'
   end if
   write(unit,'(a)')
end subroutine DrawSurface

!........................................................................

subroutine DrawElliposid
   write(unit,'(a)') '# ... object: elliposid'
   do ip = 1, np
      ipt = iptpn(ip)
      icolor = ipt
      txcolor = char(96+icolor)
      xnorm = -ori(2,3,ip)
      ynorm = ori(1,3,ip)
      znorm = Zero
      dangle = acos(ori(3,3,ip))
      if (ip == ipnpt(ipt)) then         ! first ellipsoid of type ipt
         write(unit,'(a,3f12.3)')   'Transform { translation ', ro(1:3,ip)
         write(unit,'(a,3f12.3)')   ' scale ', One, One, aellipsoid
         write(unit,'(a,4f10.4)')   ' rotation ', xnorm, ynorm, znorm, dangle
         write(unit,'(a,a,a)')      ' children [ DEF ', trim(txcolor), ' Shape{ appearance Appearance { material'
         write(unit,'(a,3f12.3,a)') ' Material {diffuseColor ', rgbcolor(1:3,icolor), '}}'
         write(unit,'(a,f6.2,a)')   ' geometry Sphere { radius', radellipsoid, '} } ] }'
      end if                             ! fix to obtain the first particle of a new color
      write(unit,'(a,3f12.3)') 'Transform { translation ', ro(1:3,ip)
      write(unit,'(a,3f12.3)')   ' scale ', One, One, aellipsoid
      write(unit,'(a,4f10.4)')   ' rotation ', xnorm, ynorm, znorm, dangle
      write(unit,'(3a)') ' children [ USE ', trim(txcolor), ']}'
   end do
   write(unit,'(a)')
end subroutine DrawElliposid

!........................................................................

subroutine DrawSuperball
   type(TriMesh) :: Mesh         ! triangle mesh with DOP-tree
   real(8) :: dist, dipref(3), dipdir(3)
   integer(4) :: m
   character(9) :: sbmesh, dipend

   write(unit,'(a)') '# ... object: superball'
   do ip = 1, np
      write(sbmesh,'(i4)') ip
      sbmesh=trim('sbmesh'//adjustl(sbmesh))
      Mesh =  superBallMesh
      Mesh%c(1,:) = ori(1,1,ip)*superBallMesh%c(1,:) + ori(1,2,ip)*superBallMesh%c(2,:) + ori(1,3,ip)*superBallMesh%c(3,:)
      Mesh%c(2,:) = ori(2,1,ip)*superBallMesh%c(1,:) + ori(2,2,ip)*superBallMesh%c(2,:) + ori(2,3,ip)*superBallMesh%c(3,:)
      Mesh%c(3,:) = ori(3,1,ip)*superBallMesh%c(1,:) + ori(3,2,ip)*superBallMesh%c(2,:) + ori(3,3,ip)*superBallMesh%c(3,:)
      Mesh%c(1,:) = ro(1,ip) + Mesh%c(1,:)
      Mesh%c(2,:) = ro(2,ip) + Mesh%c(2,:)
      Mesh%c(3,:) = ro(3,ip) + Mesh%c(3,:)
      write(unit,'(a)') 'Switch { choice [ Shape {'
      write(unit,'(a)') 'geometry DEF '//sbmesh//' IndexedFaceSet {'
      write(unit,'(a)') 'coord Coordinate { point ['
      do i = 1,size(Mesh%c,2)
          write(unit,'(3f12.3)') Mesh%c(:,i)
      end do
      write(unit,'(a)') ']}'
      write(unit,'(a)') ' coordIndex ['
      do i = 1,size(superBallMesh%t,2)
          t = superBallMesh%t(:,i)
          dir = product(sum(superBallMesh%c(:,t),2))
          if(dir < 0) t([1,2]) = t([2,1])
          write(unit,'(3i5,3x,a)') (t-1), '-1'
      end do
      write(unit,'(a)') ']'
      write(unit,'(a)') '}'
      write(unit,'(a)') '}]}'
      ipt = iptpn(ip)
      mat(1:3,1:3) = ori(:,:,ip)
      mat(4,1:3) = ro(1:3,ip)
      mat(:,4) = [0d0, 0d0, 0d0, 1d0]
      write(unit,*) 'Shape { '
      write(unit,*) 'appearance Appearance { material'
      write(unit,'(a,3f12.3,a)') ' Material {diffuseColor ', rgbcolor(1:3,ipt), '}}'
      write(unit,*) 'geometry USE '//sbmesh//' }'
   end do

   if (ldipole) then
      do m = 1, 2                             ! draw the two dipole ends
         do ip = 1, np
            ipt = iptpn(ip)
            dist = sum(dipa(1:3,1,ipt)**qsuperball2)**(-Half*qsuperballi)
            dipref(1:3) = radsuperball*dist*dipa(1:3,1,ipt)      ! intersection dipole direction supeball surface, ref frame
            dipdir(1:3) = ori(1:3,1,ip)*dipref(1) + ori(1:3,2,ip)*dipref(2) + ori(1:3,3,ip)*dipref(3) ! rotate
            dipdir(1:3) = ro(1:3,ip) -(-1)**m*dipdir(1:3) ! translate and scale with radsuperball
            write(dipend,'(i1)') m
            dipend=trim('dipend'//adjustl(dipend))
            if (ip == 1) then    ! fix to obtain the first particle of a new color
               write(unit,'(a,3f12.3)')   'Transform { translation ', dipdir(1:3)
               write(unit,'(a,a,a)')      ' children [ DEF ', dipend, ' Shape{appearance Appearance { material'
               write(unit,'(a,3f12.3,a)') ' Material {diffuseColor ', rgbcolor_dipole(1:3,m), '}}'
               write(unit,'(a,f6.2,a)')   ' geometry Sphere { radius', 0.1d0, '} } ] }'
            end if
            write(unit,'(a,3f12.3,3a)') 'Transform { translation ', dipdir(1:3), ' children [ USE ', dipend, ']}'
         end do
      end do
   end if

   write(unit,'(a)')

end subroutine DrawSuperball

!........................................................................

subroutine DrawAtom_jasper
      real(8) :: rgbloc(3)
      write(unit,'(a)') '# ... object: atom'
      do ia = 1, na
         iat = iatan(ia)
         rgbloc(1:3) = 1.0d0 - drotm(1,ia)
         write(txcolor,'(i5)') ia
         txcolor ='col_'//trim(adjustl(txcolor))
         write(*,*) 'DrawAtom_jasper, ia, txcolor, op, rgbloc', ia, txcolor, drotm(1,ia), rgbloc
         write(unit,'(a,3f12.3)')   'Transform { translation ', vaux(1:3,ia)
         write(unit,'(a,a,a)')      ' children [ DEF ', trim(txcolor), ' Shape{ appearance Appearance { material'
         write(unit,'(a,3f12.3,a)') ' Material {diffuseColor ', rgbloc(1:3), '}}'
         write(unit,'(a,f6.2,a)')   ' geometry Sphere { radius', atsize(iat), '} } ] }'
         write(unit,'(a,3f12.3,3a)') 'Transform { translation ', vaux(1:3,ia), ' children [ USE ', trim(txcolor), ']}'
      end do
      write(unit,'(a)')
end subroutine DrawAtom_jasper

!........................................................................

subroutine DrawAtom
      if(.not.allocated(icount)) then
         allocate(icount(nat))
      end if
      icount = 0
      write(unit,'(a)') '# ... object: atom'
      do ia = 1, na
         iat = iatan(ia)
         if (lgr) then                     ! each group defines a new color
            if (igrpn(ia,1) > 0) then
               icolor = igrpn(ia,1)
               txcolor = char(96+icolor)
             else
               cycle
            end if
         else                              ! each atom type defines a new color
            icolor = iat
            txcolor = txat(iat)
         end if
         icount(icolor) = icount(icolor) + 1
         if (icount(icolor) == 1) then     ! first atom of a color
            write(unit,'(a,3f12.3)')   'Transform { translation ', vaux(1:3,ia)
            write(unit,'(a,a,a)')      ' children [ DEF ', trim(txcolor), ' Shape{ appearance Appearance { material'
            write(unit,'(a,3f12.3,a)') ' Material {diffuseColor ', rgbcolor(1:3,icolor), '}}'
            write(unit,'(a,f6.2,a)')   ' geometry Sphere { radius', atsize(iat), '} } ] }'
         end if                            ! fix to obtain the first particle of a new color
         write(unit,'(a,3f12.3,3a)') 'Transform { translation ', vaux(1:3,ia), ' children [ USE ', trim(txcolor), ']}'
      end do
      write(unit,'(a)')
      deallocate(icount)
end subroutine DrawAtom

!........................................................................

subroutine DrawBond
   write(unit,'(a)') '# ... bond'
   do ib = 1, nbond
      ia = bondlist(1,ib)
      ja = bondlist(2,ib)
      if (lgr) then
         if (igrpn(ia,1) == 0) cycle
      end if
      xdir = vaux(1,ja) - vaux(1,ia)
      ydir = vaux(2,ja) - vaux(2,ia)
      zdir = vaux(3,ja) - vaux(3,ia)
      xc = Half*(vaux(1,ja) + vaux(1,ia))
      yc = Half*(vaux(2,ja) + vaux(2,ia))
      zc = Half*(vaux(3,ja) + vaux(3,ia))
      xnorm = zdir
      ynorm = Zero
      znorm =-xdir
      height = dsqrt(xdir**2+ydir**2+zdir**2)
      arg = ydir/height
      write(unit,'(a,4(f10.3))')      'Transform  { rotation', xnorm, ynorm, znorm, dacos(arg)
      write(unit,'(a,3(f10.3))')      '          translation', xc, yc, zc
      write(unit,'(a)')               'children [ Shape{ appearance Appearance { material Material { } }'
      write(unit,'(a,f10.3,a,f10.3)') 'geometry Cylinder { radius ', bondr, ' height', height
      write(unit,'(a)')               ' side TRUE   top FALSE  bottom FALSE } } ] }'
   end do
   write(unit,'(a)')
end subroutine DrawBond

!........................................................................

end subroutine VRMLSub

!************************************************************************
!*                                                                      *
!*     ImageVTF                                                         *
!*                                                                      *
!************************************************************************

! ... generate vtf file(s) and tcl-script for VMD
!
! ... VMD (Visual Molecular Dynamics) is (after registration) free to use, please see:
! ... http://www.ks.uiuc.edu/Research/vmd/
!
! ... Any published work which utilizes VMD shall include the following reference:
! ... "Humphrey, W., Dalke, A. and Schulten, K., `VMD -Visual Molecular
! ...  Dynamics', J. Molecular Graphics, 1996, vol. 14, pp. 33-38."

subroutine ImageVTF(iStage,iimage)

   use MolModule
   implicit none

   integer(4),     intent(in) :: iStage
   integer(4),     intent(in) :: iimage

   character(40),   parameter :: txroutine ='ImageVTF'
   character(80),   parameter :: txheading ='preparation of vtf file'

   character(5),         save :: txfile
   character(12),        save :: txwhen
   character(10),        save :: tximage(3)
   real(8), allocatable, save :: atsize(:), rgbcolor(:,:)
   real(8),              save :: blmax, bondr, bondres, sphres
   logical,              save :: lframezero
   logical,              save :: lgr

   integer(4),           save :: iframe
   integer(4),           save :: nframe              ! number of frames in the whole simulation

   logical,              save :: lsplitvtf = .false. ! flag for splitting the vtf file in multiple files

   character(16),   parameter :: txwrap(3) = ['# Start of image','timestep ordered','# End Image     ' ]

   integer(4),      parameter :: itypegr = 1 ! use reference group division (= 1 for ref, = 2 for field)
                                             ! ... the idea is to make itypegr an input parameter

   character(1),    parameter :: vmdname(36) = [ '0','1','2','3','4','5','6','7','8'&
                                                ,'9','A','B','C','D','E','F','G','H'&
                                                ,'I','J','K','L','M','N','O','P','Q'&
                                                ,'R','S','T','U','V','W','X','Y','Z' ]
   integer(4)                 :: igr, m
   integer(4)                 :: iat

   namelist /nmlVTF/ txfile, txwhen, tximage, atsize, rgbcolor, blmax, bondr, bondres, sphres, lframezero, lgr

   if (ltrace) call WriteTrace(2, txroutine, iStage)

   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   select case (iStage)
   case (iReadInput)

! ... allocations

      if (.not.allocated(atsize)) then
         allocate(atsize(ngr(itypegr)), rgbcolor(3,ngr(itypegr)))
      end if

! ... usually the default values of the namelist would be set here and afterwards read from the input file
! ... in order to operate with group related properties this has been shifted to iStage = iWriteInput

   case (iWriteInput)

! ... set default values

      txfile                 = 'merge'     ! alternatively choose "txfile = 'split'", for example for dynamic grouping
      txwhen                 = 'after_run' ! alternatively choose "txwhen = 'after_macro'|'after_iimage'"
      tximage                = ['frame     ','          ','          '] ! define here which kind of options shall be applied
      atsize(1:ngr(itypegr)) = radat(iatgr(1:nat,itypegr))
      do igr = 1, 3
         if (igr > ngr(itypegr)) exit
         rgbcolor(1:3,igr)   = [ (Zero, m = 1,3) ]
         rgbcolor(igr,igr)   = One
      end do
      do igr = 4, 6
         if (igr > ngr(itypegr)) exit
         rgbcolor(1:3,igr)   = [ (One, m = 1,3) ]
         rgbcolor(igr-3,igr) = Zero
      end do
      do igr = 7, ngr(itypegr)
         rgbcolor(1:3,igr)   = [ (One/(igr+m), m = 1,3) ]
      end do
      blmax                  = Zero    ! maximum bond length
      bondr                  = 0.3d0   ! bond
      bondres                = 12.0    ! number of prisms of which drawn bonds are set up of
      sphres                 = 12.0    ! number of triangles of which drawn spheres are set up of
      lframezero             = .true.  ! set to .false. to exclude the frame containing the initial configuration
      lgr                    = .true.  ! set to .false. to force coloring according to atom types

! ... read input data

      rewind(uin)
      read(uin,nmlVTF)

! ... number of frames to be made

      select case (txwhen)
      case ('after_run')
         nframe = 1
      case ('after_macro')
         nframe = nstep1
      case ('after_iimage')
         nframe = nstep/iimage
      case default
         call Stop(txroutine,'unsupported value of txwhen',uout)
      end select

! ... initialization of first frame

      iframe = -1 ! first frame iframe = 0: iframe is incremented before the first frame is stored

! ... if vtf file shall be split, prepare some variables

      if (txfile == 'merge') then
         continue ! fvtf = project.vtf (see inititalization in molsim.F90)
      else if (txfile == 'split') then
         call UpdateVTFFileName(iframe,nframe) ! Prepare format string [with (!) optional argument nframe]
         lsplitvtf = .true.
      else
         call Stop(txroutine,'unsupported value of txfile',uout)
      end if

! ... open, write and close tcl-script

      if (master .and. (txstart == 'setconf' .or. txstart == 'zero' .or. txstart == 'readfin')) then
         call FileOpen(utcl,ftcl,'form/noread')
         call WriteTCLScript(iStage,rgbcolor,bondr,bondres,sphres,tximage,vmdname,lgr,utcl) ! TODO
         close(utcl)
      end if

   case (iBeforeSimulation)

      if (master) then
         if (txstart == 'setconf' .or. txstart == 'zero' .or. txstart == 'readfin') then
            if (lsplitvtf) then
               if (lframezero) call ImageVTFSub
            else
               call FileOpen(uvtf, fvtf, 'form/noread')
               call WriteVTFHeader(atsize,blmax,vmdname,uvtf)
               if (lframezero) call ImageVTFSub
            end if
         else if (txstart == 'continue') then
            read(ucnf) iframe
            call FileOpen(uvtf, fvtf, 'form/read')
         end if
      end if

   case (iSimulationStep)

      if (txwhen == 'after_iimage') call ImageVTFSub

   case (iAfterMacrostep)

      if (txwhen == 'after_macro') call ImageVTFSub

      if (lsim .and. master) write(ucnf) iframe

   case (iAfterSimulation)

      if (txwhen == 'after_run') call ImageVTFSub

      call WriteHead(2, txheading, uout)
      write(uout,'(a,a   )') 'generating vtf file                 ', txwhen
      write(uout,'(a,4a  )') 'options (tximage)                   ', tximage
      write(uout,'(a,f8.2)') 'maximum bond length               = ', blmax
      write(uout,'(a,f8.2)') 'bond radius                       = ', bondr
      write(uout,'()')
      write(uout,'(a,t15,a,t30,a,t50,a)') 'atom type no', 'atom type', 'size (radius)', 'rgbcolor'
      write(uout,'(a,t15,a,t30,a,t50,a)') '------------', '---------', '-------------', '--------'
      write(uout,'(i3,t15,a,t30,f8.3,t45,3f6.2)') (iat,txat(iat),atsize(iat),rgbcolor(1:3,iat),iat = 1,nat)
      write(uout,'()')
      write(uout,'(a,i4)')   'number of images made             = ', iframe+1

      if (allocated(atsize)) deallocate(atsize,rgbcolor)

      if (master) close(uvtf)

      ! call system ('sed -i "s?trajectory.vtf?$(ls -1|grep vtf)?g" *.tcl')  ! change from trajectory.vtf to $Job.vtf in *.tcl files
!     call execute_command_line ('sed -i "s?trajectory.vtf?$(ls -1|grep vtf)?g" *.tcl', wait=.false.)

   end select

   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

contains

! ........................................................................

subroutine ImageVTFSub

   character(40), parameter :: txroutine = 'ImageVTFSub'
   character(10) :: str

   iframe = iframe+1

   if (lsplitvtf) then
      call UpdateVTFFileName(iframe,nframe)
      call FileOpen(uvtf, fvtf, 'form/noread')
      call WriteVTFHeader(atsize,blmax,vmdname,uvtf)
   end if

   str = '      '
   write(str, '(i10)') iframe

   if (txwhen == 'after_run') then
      write (uvtf,'(a)') txwrap(2)
      call WriteVTFCoordinates(tximage, uvtf)
   else if (txwhen == 'after_macro') then
      write (uvtf,'(a)') txwrap(1)//trim(adjustl(str))
      write (uvtf,'(a)') txwrap(2)
      call WriteVTFCoordinates(tximage, uvtf)
      write (uvtf,'(a)') txwrap(3)
   else if (txwhen == 'after_iimage') then
      write (uvtf,'(a)') txwrap(1)//trim(adjustl(str))
      write (uvtf,'(a)') txwrap(2)
      call WriteVTFCoordinates(tximage, uvtf)
      write (uvtf,'(a)') txwrap(3)
   else
      call stop(txroutine, 'unsupported value of txwhen', uout)
   end if

   if (lsplitvtf) close (uvtf)

end subroutine ImageVTFSub

! .......................................................................

end subroutine ImageVTF

!************************************************************************
!*                                                                      *
!*     UpdateVTFFileName                                                *
!*                                                                      *
!************************************************************************

! ... update vtf file name with respect to a given frame iframe

subroutine UpdateVTFFileName(iframe,nframe)

   use MolModule
   implicit none

   integer(4), intent(in) :: iframe
   integer(4), intent(in) :: nframe

   character(6),     save :: framefmt

   character(1)           :: ndigit
   character(10)          :: ifrfmt

   logical,          save :: first = .true.

   if (first) then
      write(ndigit,'(i1)') int(floor(log10(real(nframe))))+1
      framefmt = '(i'//ndigit//'.'//ndigit//')'
      first = .false.
   else
      write(ifrfmt,framefmt) iframe
      fvtf = adjustl(trim(project))//'.'//adjustl(trim(ifrfmt))//'.vtf'
   end if

end subroutine UpdateVTFFileName

!************************************************************************
!*                                                                      *
!*     WriteVTFHeader                                                   *
!*                                                                      *
!************************************************************************

! ... write header of the vtf file

subroutine WriteVTFHeader(atsize, blmax, vmdname, unit)

   use MolModule
   implicit none

   real(8),      intent(in)      :: atsize(*)   ! size of atom type
   real(8),      intent(in)      :: blmax       ! maximal bond length
   character(1), intent(in)      :: vmdname(*)  ! label to identify different particle types
   integer(4)                    :: unit        ! output unit

   character(40),      parameter :: txroutine = 'WriteVTFHeader'
   integer(4)                    :: ibond, ia
   integer(4),              save :: nbond, mnbond   ! actual and maximal number of bonds
   integer(4), allocatable, save :: bondlist(:,:)   ! pair of atoms joined by a bond
   character(10)                 :: str = '          '

! ... undo periodic boundary conditions for chains and hierarchical structures

   vaux(1:3,1:na) = r(1:3,1:na)
   call UndoPBC(vaux)

! ... declare atoms

   write(unit,'(a5,i5,a8,f6.3,a6,a11,a6,a)') &
      ('atom ',ia-1,' radius ',atsize(iatan(ia)),' type ',txat(iatan(ia)),' name ',vmdname(iatan(ia)), ia = 1, na)
   write(unit,'(/)')

! ... determine connectivity of atoms

   mnbond = 2*na
   if (.not. allocated(bondlist)) then
      allocate(bondlist(2,mnbond))
      bondlist = 0
   end if
   call MakeBondList(vaux,mnbond,blmax,nbond,bondlist)

! ... declare bonds

   do ibond = 1, nbond
      write(str,'(i10)') bondlist(2,ibond)-1
      write(unit,'(a5,i5,a)') 'bond ',bondlist(1,ibond)-1,':'//trim(adjustl(str))
   end do
   write(unit,'(/)')

end subroutine WriteVTFHeader

!************************************************************************
!*                                                                      *
!*     WriteVTFCoordinates                                              *
!*                                                                      *
!************************************************************************

! ... writes current atom coordinates to vtf-file

subroutine WriteVTFCoordinates(tximage, unit)

   use MolModule
   use CoordinateModule  ! needed ?

   implicit none

   character(20),           intent(in) :: tximage(4)       !
   integer(4),              intent(in) :: unit             ! output unit

   character(40),            parameter :: txroutine = 'WriteVTFCoordinates'

   real(8),          allocatable, save :: ro_vtf(:,:)
   real(8)                             :: rcom(3)
   integer(4)                          :: ip

   if (.not. allocated(ro_vtf)) then
      allocate(ro_vtf(1:3,1:na))
      ro_vtf = 0.0E+00
   end if

! ... write coordinates in ro_vtf

   ro_vtf(1:3,1:na) = r(1:3,1:na)

! ... if intended undo periodic boundary conditions

   if (tximage(3) == 'undopbc') call UndoPBC(ro_vtf)

! ... if centering of the atoms is required

   rcom(1:3) = Zero

   if (tximage(4) == 'center') then
      rcom(1)   = sum(ro_vtf(1,1:na))/real(na)
      rcom(2)   = sum(ro_vtf(2,1:na))/real(na)
      rcom(3)   = sum(ro_vtf(3,1:na))/real(na)
      ro_vtf(1,1:na) = ro_vtf(1,1:na)-rcom(1)
      ro_vtf(2,1:na) = ro_vtf(2,1:na)-rcom(2)
      ro_vtf(3,1:na) = ro_vtf(3,1:na)-rcom(3)
   else if (tximage(4) == 'xycenter') then
      rcom(1)   = sum(ro_vtf(1,1:na))/real(na)
      rcom(2)   = sum(ro_vtf(2,1:na))/real(na)
      ro_vtf(1,1:na) = ro_vtf(1,1:na)-rcom(1)
      ro_vtf(2,1:na) = ro_vtf(2,1:na)-rcom(2)
   ! else if (tximage(4) == 'centernw' ) then               ! center network
   !    if (tximage(3) /= 'undopbc') call UndoPBC(ro_vtf)
   !    if (lnomassnw) then
   !       rcom(1) = sum(ro_vtf(1,1:np) - ro_vtf(1,ipref), MASK=lptinnw(iptpn(1:np)))
   !       rcom(2) = sum(ro_vtf(2,1:np) - ro_vtf(2,ipref), MASK=lptinnw(iptpn(1:np)))
   !       rcom(3) = sum(ro_vtf(3,1:np) - ro_vtf(3,ipref), MASK=lptinnw(iptpn(1:np)))
   !    else
   !       rcom(1) = sum((ro_vtf(1,1:np) - ro_vtf(1,ipref)) * masspt(iptpn(1:np)), MASK=lptinnw(iptpn(1:np)))
   !       rcom(2) = sum((ro_vtf(2,1:np) - ro_vtf(2,ipref)) * masspt(iptpn(1:np)), MASK=lptinnw(iptpn(1:np)))
   !       rcom(3) = sum((ro_vtf(3,1:np) - ro_vtf(3,ipref)) * masspt(iptpn(1:np)), MASK=lptinnw(iptpn(1:np)))
   !    end if
   !    rcom(1:3) = rcom(1:3)*invmass + ro_vtf(1:3,ipref)
   !    do ip = 1, np
   !       ro_vtf(1:3,ip) = ro_vtf(1:3,ip) - rcom(1:3)
   !       call PBC(ro_vtf(1,ip),ro_vtf(2,ip),ro_vtf(3,ip))
   !    end do
   end if

! ... write coordinate block into vtf-file

   write(unit,'(3f20.5)') (ro_vtf(1:3,ip), ip = 1, np)

end subroutine WriteVTFCoordinates

!************************************************************************
!*                                                                      *
!*     WriteTCLScript                                                   *
!*                                                                      *
!************************************************************************

! ... TCL: "tool command language"
!
! ... write VTF-accompanying TCL-script to be executed in VMD to adjust colors, bond radius,
! ... bond and atom resolution, and insert objects such as frames, planes

subroutine WriteTCLScript(iStage,rgbcolor,bondr,bondres,sphres,tximage,vmdname,lgr,unit)

   use MolModule
   implicit none

   integer(4),    intent(in) :: iStage
   real(8),       intent(inout) :: rgbcolor(3,*)
   real(8),       intent(in) :: bondr
   real(8),       intent(in) :: bondres
   real(8),       intent(in) :: sphres
   character(10), intent(in) :: tximage(4)
   character(1) , intent(in) :: vmdname(*)
   logical,       intent(in) :: lgr
   integer(4),    intent(in) :: unit

   character(40),  parameter :: txroutine = 'WriteTCLScript'
   integer(4)                :: iat, icolor, icube, ird
   real(8)      ,  parameter :: cornerref(3,8) = reshape( &
           [ One, One, One,   -One, One, One,  -One, One,-One,  One, One,-One, &
             One,-One, One,   -One,-One, One,  -One,-One,-One,  One,-One,-One ] , [3,8] )
   real(8) :: corner(1:3,1:14)

   select case (iStage)
   case (iWriteInput)

      write(unit,'(a)') 'set mode load'

      write(unit,'(a)') 'if { $mode == "load" } {'

! ... define names used to apply coloring scheme

      write(unit,'(a)') '   set init [mol new atoms 1]'
      write(unit,'(a)') '   set sel  [atomselect $init all]'
      write(unit,'(a17,a1)') ('   $sel set name ', vmdname(iat), iat = 1, nat)
      write(unit,'(a)') '   $sel delete'
      write(unit,'(a)') '   mol delete $init'

! ... rgbcolors should be set between 0 and 1

      if(any(rgbcolor(1:3,1:nat) > One)) rgbcolor(1:3,1:nat) = rgbcolor(1:3,1:nat)/255.0

 ! ... write script header

      write(unit,'(a)') '   mol delete all                  '   ! begin a fresh session
      write(unit,'(a)') '   display update off              '   ! update of display after modifications

! ... load scene

      write(unit,'(a)') '   mol load vtf trajectory.vtf     '   ! load frames
      write(unit,'(a)') '   color Display Background white  '   ! background default color is black -> turn to white
      write(unit,'(a)') '   color Axes Labels black         '   ! axes labels default color is white -> turn to black
      write(unit,'(a)') '   display depthcue off            '   ! disable depth cueing
      write(unit,'(a)') '   set molID [molinfo top]         '   ! get ID of molecule

! ... adjust drawing style, bond radius, bond resolution, sphere resolution

      write(unit,'(a)')             '   mol delrep 0 $molID        '
      write(unit,'(a30,f6.1)')      '   mol representation VDW 1.0 ', sphres          ! 1.0 is a radius scaling factor - the atom radius is declared in WriteVTFHeader
      write(unit,'(a)')             '   mol addrep $molID          '
      write(unit,'(a28,f5.2,f6.1)') '   mol representation Bonds   ', bondr, bondres
      write(unit,'(a)')             '   mol addrep $molID          '

! ... adjust colors

      do iat = 1, nat
         if(lgr) then
            if(igrpn(ipnpt(iptat(iat)),1) > Zero) then
               icolor = igrpn(ipnpt(iptat(iat)),1)       ! colouring follows group assignment
            end if
         else
            icolor = iat                                 ! colouring according to atom types
         end if
         write(unit,'(a20,i4,3f6.3)') '   color change rgb ',1025-iat, rgbcolor(1:3,icolor) ! overwrite colors (count-down from color no. 1024 in order to prevent the modification of already used colo
         write(unit,'(a14,a2,i5)')    '   color Name ', vmdname(iat), 1025-icolor
      end do

! ... insert graphical objects as provided by tximage

! ... draw frame according to the geometry of the simulation cell

      if (tximage(1) == 'frame') then
         write(unit,'(a)') '   mol load graphics frame      '
         write(unit,'(a)') '   set frameID [molinfo top]'      ! store the ID of the frame graphics
         write(unit,'(a)') '   mol top $molID               '  ! molecules can be active or inactive - the actual molecule should be the active one
         write(unit,'(a)') '   graphics $frameID color black'
         if (lbcbox) then
            call DrawBoxTCL
         else if (lbcrd) then
            call DrawRhombicDodecahedronTCL
  !      else if (lbcto) then
  !         call DrawTruncatedOctahedronTCL
  !      else if (lbcsph) then
  !         call DrawSphereTCL
  !      else if (lbccyl) then
  !         call DrawCylinderTCL
         else
            call Warn(txroutine,'Drawing of simulation cell currently only for lbcbox or lbcrd',uout)
         end if
      end if

! ... show result

      write(unit,'(a)') '   display update on    ' ! show modifications
      write(unit,'(a)') '}'

   case (iSimulationStep)

      continue

   case (iAfterSimulation)

      continue

   end select

contains

! .......................................................................

subroutine DrawBoxTCL

   corner(1,1:8) = boxlen2(1)*cornerref(1,1:8)
   corner(2,1:8) = boxlen2(2)*cornerref(2,1:8)
   corner(3,1:8) = boxlen2(3)*cornerref(3,1:8)

   do icube = 0, 3
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(icube,4)+1),' } { ', corner(1:3,mod(icube+1,4)+1), '} width 2'
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(icube,4)+5),' } { ', corner(1:3,mod(icube+1,4)+5), '} width 2'
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(icube,4)+1),' } { ', corner(1:3,mod(icube,4)+5)  , '} width 2'
   end do

end subroutine DrawBoxTCL

! .......................................................................

subroutine DrawRhombicDodecahedronTCL

   real(8) :: s, d

   d = sqrt(Two/Three)*cellside
   s = d/SqTwo

   corner(1:3,1)  = [    d,  Zero,       -s ]
   corner(1:3,2)  = [ Zero,  Zero, -SqTwo*d ]
   corner(1:3,3)  = [   -d,  Zero,       -s ]
   corner(1:3,4)  = [   -d,     d,     Zero ]
   corner(1:3,5)  = [ Zero,     d,        s ]
   corner(1:3,6)  = [    d,     d,     Zero ]
   corner(1:3,7)  = [    d,    -d,     Zero ]
   corner(1:3,8)  = [ Zero,    -d,       -s ]
   corner(1:3,9)  = [   -d,    -d,     Zero ]
   corner(1:3,10) = [   -d,  Zero,        s ]
   corner(1:3,11) = [ Zero,  Zero,  SqTwo*d ]
   corner(1:3,12) = [    d,  Zero,        s ]
   corner(1:3,13) = [ Zero,    -d,        s ]
   corner(1:3,14) = [ Zero,     d,       -s ]

   do ird = 0, 5
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(ird,6)+1), ' } { ', corner(1:3,mod(ird+1,6)+1),  '} width 2'
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(ird,6)+7), ' } { ', corner(1:3,mod(ird+1,6)+7),  '} width 2'
      write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(ird,6)+1), ' } { ', corner(1:3,mod(ird,6)+7),    '} width 2'
      if      ( mod(ird,2) == 1 ) then
         write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(ird,6)+1), ' } { ', corner(1:3,14),  '} width 2'
      else if ( mod(ird,2) == 0 ) then
         write(unit,'(a28,3f8.1,a6,3f8.1,a12)') '   graphics $frameID line { ', corner(1:3,mod(ird,6)+7), ' } { ', corner(1:3,13),  '} width 2'
      end if
   end do

end subroutine DrawRhombicDodecahedronTCL

! ......................................................................

end subroutine WriteTCLScript

module UndoPBCModule

   use MolModule, only: np, ro, iptpn
   use MolModule, only: lhierarchical, lclink, bondnn, nbondcl, bondcl
   use MolModule, only: na, ianpn, napt, r
   implicit none
   public   :: ipatcenter
   public   :: UndoclPBC
   private  :: ipclose

   real(8)  , allocatable, public  :: rotmp(:,:)
   logical  , allocatable, public  :: lundo(:)

   logical, private, allocatable :: loclundoip(:)


   contains

   function ipatcenter(ip) result(ipcenter)
      implicit none
      integer, intent(in)  :: ip
      integer  :: ipcenter
      real  :: r2

      if (.not. allocated(loclundoip)) then
         allocate(loclundoip(np))
      end if
      loclundoip = .false.

      call ipclose(ip, ipcenter,r2)

      deallocate(loclundoip)

   end function ipatcenter

   recursive subroutine ipclose(ip, ipcenter, r2min)
      implicit none
      integer, intent(in)  :: ip
      integer, intent(out)  :: ipcenter
      real,    intent(out)  :: r2min
      integer  :: jpcenter
      integer  :: jp, icl, ib
      real  :: r2

      if (loclundoip(ip) == .true.) then
         ipcenter = 0
         r2min = huge(r2min)
      else

         loclundoip(ip) = .true.
         ipcenter = ip
         r2min = (ro(1,ip)**2 + ro(2,ip)**2 + ro(3,ip)**2)

         do ib = 1, 2      ! check bonded partners
            jp = bondnn(ib,ip)
            if(jp /= 0 ) then
               call ipclose(jp, jpcenter, r2)
               if( (jpcenter > 0) .and. (r2 < r2min))  then
                  ipcenter = jpcenter
                  r2min = r2
               end if
            end if
         end do

         if(lhierarchical .or. lclink) then
            do icl = 1, nbondcl(ip) ! check crosslinked partners
               jp = bondcl(icl,ip)
               call ipclose(jp, jpcenter, r2)
               if( (jpcenter > 0) .and. (r2 < r2min))  then
                  ipcenter = jpcenter
                  r2min = r2
               end if
            end do
         end if

      end if

   end subroutine ipclose
!
!  UndoClPBC undoes all periodic boundary conditions along the bonds and crosslinked particles
!
   recursive subroutine UndoClPBC(rref, ip)
      implicit none
      real(8), intent(in)     :: rref(3)     ! reference point for the undo
      integer, intent(in)     :: ip     ! particle to be undone for the undo

      real(8)  ::  dr(3)
      integer  :: ia, ipt, ib, jp, icl
      real(8)  :: rip(3)     ! reference point for the undo

      if(lundo(ip)) then
         return
      else
         dr(1:3) = ro(1:3,ip) - rref(1:3)
         call PBC(dr(1), dr(2), dr(3))
         rip(1:3) = rref(1:3) + dr(1:3)

         do ia = ianpn(ip), ianpn(ip) + napt(iptpn(ip)) - 1
            rotmp(1:3, ia) = r(1:3,ia) - ro(1:3,ip) + rip(1:3)
         end do

         lundo(ip) = .true.

         do ib = 1, 2      ! check bonded partners
            jp = bondnn(ib,ip)
            if(jp /= 0 ) then
               call UndoClPBC(rip(1:3), jp)
            end if
         end do

         if(lhierarchical .or. lclink) then
            do icl = 1, nbondcl(ip) ! check crosslinked partners
               jp = bondcl(icl,ip)
               call UndoClPBC(rip(1:3), jp)
            end do
         end if

      end if

   end subroutine UndoClPBC

end module UndoPBCModule

!************************************************************************
!*                                                                      *
!*     UndoPBC                                                          *
!*                                                                      *
!************************************************************************

! ... undo periodic boundary conditions

subroutine UndoPBC(vhelp)

   use UndoPBCModule
   implicit none
   real(8),    intent(out)  :: vhelp(1:3,*)       ! undone atom position
   integer :: ip, ipmin, jp


   if(.not. allocated(lundo)) then
      allocate(lundo(np), rotmp(3,na))
   end if

   lundo = .false.
   rotmp = 0.0

   do ip = 1, np
      if(.not. lundo(ip)) then
!          jp = ip
         ipmin = ipatcenter(ip)
         call UndoClPBC(ro(1:3,ipmin),ipmin)
      end if
   end do

   vhelp(1:3,1:na) = rotmp(1:3, 1:na)

   deallocate(lundo, rotmp)

end subroutine UndoPBC

!************************************************************************
!*                                                                      *
!*     MakeBondList                                                     *
!*                                                                      *
!************************************************************************

! ... make bond list

subroutine MakeBondList(vhelp, mnbond, blmax, nbond, bondlist)

   use MolModule
   implicit none

   character(40), parameter :: txroutine ='MakeBondList'
   real(8),    intent(in)  :: vhelp(1:3,*)     ! atom position
   integer(4), intent(in)  :: mnbond           ! maximum number of bonds
   real(8),    intent(in)  :: blmax            ! maximum bond length
   integer(4), intent(out) :: nbond            ! number of bonds
   integer(4), intent(out) :: bondlist(2,*)    ! pair of atoms jointed by a bond
   integer(4) :: i, ic, ict, ip, jp, ia, ja, ialow, iaupp, iseg
   real(8)    :: dx, dy, dz, r2

   nbond = 0

! ... make bonds between atoms residing in the same particles and within separation blmax

   if (blmax > Zero) then
      do ip = 1, np
         ialow = ianpn(ip)
         iaupp = ianpn(ip)+napt(iptpn(ip))-1
         do ia = ialow, iaupp
            do ja = ia+1, iaupp
               r2 = (vhelp(1,ia)-vhelp(1,ja))**2+(vhelp(2,ia)-vhelp(2,ja))**2+(vhelp(3,ia)-vhelp(3,ja))**2
               if (r2 < blmax**2) then
                  nbond = nbond+1
                  if (nbond > mnbond) call Stop(txroutine, 'nbond > mnbond', 6)
                  bondlist(1,nbond) = ia
                  bondlist(2,nbond) = ja
               end if
            end do
         end do
      end do
   end if

! ... make bonds between connected particles in chains but exclude those passing the edge of the box
!     works for lpolyatom

   if (lchain) then
      do ic = 1, nc
         ict = ictcn(ic)
         do iseg = 1, npct(ict)-1
            ip = ipnsegcn(iseg,ic)
            jp = ipnsegcn(iseg+1,ic)
            dx = vhelp(1,ip)-vhelp(1,jp)
            dy = vhelp(2,ip)-vhelp(2,jp)
            dz = vhelp(3,ip)-vhelp(3,jp)
            if (txbc == 'xyz' .and. (abs(dx) > boxlen2(1) .or. abs(dy) > boxlen2(2) .or. abs(dz) > boxlen2(3))) cycle
            if (txbc == 'xy' .and. (abs(dx) > boxlen2(1) .or. abs(dy) > boxlen2(2))) cycle
            if (txbc == 'z' .and. abs(dz) > boxlen2(3)) cycle
            nbond = nbond+1
            if (nbond > mnbond) call Stop(txroutine, 'nbond > mnbond', 6)
            bondlist(1,nbond) = ianpn(ip)
            bondlist(2,nbond) = ianpn(jp)
         end do
      end do
   end if

! ... make crosslinks between particles but exclude those passing the edge of the box
!     works for lpolyatom

   if (lclink) then
      do ip = 1, np
         do i = 1, nbondcl(ip)
            jp = bondcl(i,ip)
            if (jp < ip) cycle
            dx = vhelp(1,ip)-vhelp(1,jp)
            dy = vhelp(2,ip)-vhelp(2,jp)
            dz = vhelp(3,ip)-vhelp(3,jp)
            if (txbc == 'xyz' .and. (abs(dx) > boxlen2(1) .or. abs(dy) > boxlen2(2) .or. abs(dz) > boxlen2(3))) cycle
            if (txbc == 'xy' .and. (abs(dx) > boxlen2(1) .or. abs(dy) > boxlen2(2))) cycle
            if (txbc == 'z' .and. abs(dz) > boxlen2(3)) cycle
            nbond = nbond+1
            if (nbond > mnbond) call Stop(txroutine, 'nbond > mnbond', 6)
            bondlist(1,nbond) = ianpn(ip)
            bondlist(2,nbond) = ianpn(jp)
         end do
      end do
   end if

end subroutine MakeBondList
