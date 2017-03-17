module CellListModule

implicit none
private
public UpdateCellip, InitCellList, SetCellList, CellListAver
public cell, pcellro, cell_type, cell_pointer_array, nneighcell

integer(4), parameter  :: nneighcell = 27
type cell_pointer_array
   type(cell_type), pointer  :: p => null()
end type cell_pointer_array

type cell_type
   integer(4)  :: id       ! fore easy recognition
   integer(4)  :: npart
   type(cell_pointer_array) :: neighcell(nneighcell)
   integer(4), allocatable :: ip(:)
end type cell_type

type(cell_type), target, allocatable  :: cell(:,:,:)    !cells
type(cell_pointer_array), allocatable  :: cellip(:) !cell of each particle
integer(4)  :: ncell(3)                !number of cells in x y z in each octant
real(8)     :: ircell                  !inverse edge length of each cell

contains

subroutine InitCellList(rcut, iStage)

   use MolModule, only: ltrace, ltime, uout
   use MolModule, only: lbcbox, boxlen2
   use MolModule, only: np

   real(8), intent(in)  :: rcut
   integer(4), intent(in)  :: iStage
   character(40), parameter :: txroutine ='InitCellList'
   integer(4)  :: ix, iy, iz, jx, jy, jz, neigh(3), ineigh, id

   if (ltrace) call WriteTrace(1, txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(rcut .le. 0.0d0) then
      call Stop(txroutine, 'rcut .le. 0.0', uout)
   end if
   if(.not. lbcbox) then
      call Stop(txroutine, 'celllist needs cubic box (lbcbox should be true)', uout)
   end if
   ncell(1:3) = ceiling(boxlen2(1:3)/rcut)
   ircell = 1.0d0/rcut

   if(allocated(cell)) then
      deallocate(cell)
   end if
   if(allocated(cellip)) then
      deallocate(cellip)
   end if
   allocate(cell(-ncell(1):(ncell(1)-1),-ncell(2):(ncell(2)-1),-ncell(3):(ncell(3)-1)))
   allocate(cellip(1:np))
   ! cell structure when ncell(1:3) = 2
   !        +-----+-----+-----+-----+
   !       /.., 1/     /     /     /|
   !      +-----+-----+-----+-----+ |
   !     /.., 0/     /     /     /| +
   !    +-----+-----+-----+-----+ |/|
   !   /..,-1/     /     /     /| + |
   !  +-----+-----+-----+-----+ |/| +
   ! /..,-2/     /     /     /| + |/|
   !+-----+-----+-----+-----+ |/| + |
   !|-2, 1|-1, 1| 0, 1| 1, 1| + |/| +
   !| ,-2 | ,-2 | ,-2 | ,-2 |/| + |/|
   !+-----+-----+-----+-----+ |/| + |
   !|-2, 0|-1, 0| 0, 0| 1, 0| + |/| +
   !| ,-2 | ,-2 | ,-2 | ,-2 |/| + |/
   !+-----+-----+-----+-----+ |/| +
   !|-2,-1|-1,-1| 0,-1| 1,-1| + |/
   !| ,-2 | ,-2 | ,-2 | ,-2 |/| +
   !+-----+-----+-----+-----+ |/
   !|-2,-2|-1,-2| 0,-2| 1,-2| +
   !| ,-2 | ,-2 | ,-2 | -2  |/
   !+-----+-----+-----+-----+

   !loop over all cells
   id = 0
   do ix = -ncell(1), ncell(1) -1
      do iy = -ncell(2), ncell(2) -1
         do iz = -ncell(3), ncell(3) -1

            id = id + 1
            cell(ix,iy,iz)%id = id
            !allocate cell
            allocate(cell(ix,iy,iz)%ip(np))
            cell(ix,iy,iz)%npart = 0

            !make pointers to neighbouring cells
            ineigh = 0
            do jx = -1, 1
               do jy = -1, 1
                  do jz = -1, 1
                     neigh(1:3) = (/ ix + jx, iy + jy , iz + jz /)
                     where (neigh .ge. ncell)
                        neigh = -ncell
                     elsewhere (neigh < -ncell)
                        neigh = ncell - 1
                     end where
                     ineigh = ineigh + 1
                     cell(ix,iy,iz)%neighcell(ineigh)%p => cell(neigh(1), neigh(2), neigh(3))
                  end do
               end do
            end do

         end do
      end do
   end do



end subroutine InitCellList

pure function pcellro(ro) result(icell)
   implicit none
   real(8), intent(in)  :: ro(3)
   type(cell_type), pointer :: icell
   integer(4)  :: i(3)

   i = min(ncell-1,floor(ro*ircell))
   icell => cell(i(1), i(2), i(3))
end function pcellro

subroutine AddIpToCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), target, intent(inout) :: icell

   icell%npart = icell%npart + 1
   icell%ip(icell%npart) = ip
   cellip(ip)%p => icell

end subroutine AddIpToCell

subroutine RmIpFromCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), intent(inout) :: icell
   integer(4)  :: npart

   npart = icell%npart
   icell%ip(1:(npart -1)) = pack(icell%ip(1:npart), (icell%ip(1:npart) .ne. ip))
   icell%npart = npart - 1
   cellip(ip)%p => null()

end subroutine RmIpFromCell

subroutine UpdateCellip(ip)
   use MolModule, only: ro
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), pointer :: cellold
   type(cell_type), pointer :: cellnew

   cellold => cellip(ip)%p
   cellnew => pcellro(ro(1:3,ip))

   if(cellold%id .ne. cellnew%id) then !do not update cell if cell has not changed
      call RmIpFromCell(ip, cellold)
      call AddIpToCell(ip, cellnew)
   end if
end subroutine UpdateCellip

subroutine SetCellList

   use MolModule, only: np, ro
   implicit none

   integer(4)  :: ip
   type(cell_type), pointer :: celltmp

   cell(:,:,:)%npart = 0
   do ip = 1, np
      celltmp => pcellro(ro(1:3,ip))
      call AddiptoCell(ip, celltmp)
   end do

end subroutine SetCellList

subroutine CellListAver(iStage)

   use MolModule, only: uout
   use MolModule, only: np
   use MolModule, only: master, nstep1
   use MolModule, only: iWriteInput, iAfterMacrostep, iAfterSimulation
   implicit none
   integer(4), intent(in) :: iStage
   character(80), parameter :: txheading ='cell list statistics'
   real(8), save    :: nppp(4)         ! number of particles per processor
   real(8), save    :: nneigh(4)       ! number of neigbours per particle and processor
   type(cell_type), pointer   :: tmpcell
   integer(4)  :: ip, icell
   real(8)     :: nneighip


   if (master) then
      select case (iStage)
      case (iWriteInput)

         nppp(1:2)    = 0.0d0
         nppp(3)    =+huge(1.0d0)
         nppp(4)    =-huge(1.0d0)
         nneigh(1:2)  = 0.0d0
         nneigh(3)  =+huge(1.0d0)
         nneigh(4)  =-huge(1.0d0)

      case (iAfterMacrostep)

         nppp(1) = nppp(1) + sum(cell(:,:,:)%npart)
         nppp(2) = nppp(2) + sum(cell(:,:,:)%npart**2)
         nppp(3) = min(nppp(3), real(minval(cell(:,:,:)%npart)))
         nppp(4) = max(nppp(4), real(maxval(cell(:,:,:)%npart)))

         do ip=1, np
            nneighip = -1.0d0 !do not count the particle itself as a neighbour
            do icell = 1, nneighcell
               tmpcell => cellip(ip)%p%neighcell(icell)%p
               nneighip = nneighip + tmpcell%npart
            end do
            nneigh(1) = nneigh(1) + nneighip
            nneigh(2) = nneigh(2) + nneighip**2
            nneigh(3) = min(nneigh(3), nneighip)
            nneigh(4) = max(nneigh(4), nneighip)
         end do

      case (iAfterSimulation)

         nppp(1)   = nppp(1)/(nstep1*8*product(ncell(1:3)))
         nppp(2)   = sqrt(nppp(2)/(nstep1*8*product(ncell(1:3)))-nppp(1)**2)
         nneigh(1) = nneigh(1)/(nstep1*np)
         nneigh(2) = sqrt(nneigh(2)/(nstep1*np)-nneigh(1)**2)

         call WriteHead(2, txheading, uout)
         write(uout,'()')
         write(uout,'(t12,a,t34,a)') &
         'no of part per cell', 'no of neighbours per part'
         write(uout,'(t12,a,t34,a)') &
         '-------------------', '----------------------------------'
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'average', nppp(1), nneigh(1)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'one sd ', nppp(2), nneigh(2)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'minimum', nppp(3), nneigh(3)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'maximum', nppp(4), nneigh(4)

      end select
   end if

   end subroutine CellListAver
end module CellListModule
