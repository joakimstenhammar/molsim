module CellList

implicit none
private
public UpdateCellip

type cell_pointer
   type(cell), pointer  :: p => null()
end type cell_pointer

type cell_type
   integer(4)  :: id       ! fore easy recognition
   integer(4)  :: npart
   type(cell_pointer), pointer  :: neighcell(27)
   integer(4), allocatable :: ip(:)
end type cell_type

type(cell_pointer)  :: cellip(:) !cell of each particle
integer(4)  :: ncell(3)                !number of cells in x y z in each octant
real(8)     :: rcell                   !edge length of each cell
real(8)     :: ircell                  !inverse edge length of each cell
type(cell_type), target, allocatable  :: cell(:,:,:)    !cells

contains

subroutine InitCellList(rcut, iStage)

   use MolModule, only: ltrace, ltime, uout
   use MolModule, only: lbcbox, boxlen2

   real(8), intent(in)  :: rcut
   integer(4), intent(in)  :: iStage
   character(40), parameter :: txroutine ='InitCellList'
   integer(4)  :: ix, iy, iz, jx, jy, jy, neigh(3), ineigh, id

   if (ltrace) call WriteTrace(1, txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(rcut .le. 0.0d0) then
      call Stop(txroutine, 'rcut .le. 0.0', uout)
   end if
   ncell(1:3) = ceiling(boxlen2(1:3)/rcut)
   rcell = rcut
   ircell = 1.0d0/rcell

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
            cel(ix,iy,iz)%id = id
            !allocate cell
            allocate(cell(ix,iz,iy)%ip(np))
            cell(ix,iy,iz)%npart = 0

            !make pointers to neighbouring cells
            ineigh = 0
            do jx = -1, 1
               do jy = -1, 1
                  do jz = -1, 1
                     neigh(1:3) = (/ ix + jx, iy + jy , iz + jz \)
                     where (neigh .ge. ncell)
                        neigh = -ncell
                     elsewhere (neigh < -ncell)
                        neigh = ncell - 1
                     end where
                     ineigh = ineigh + 1
                     cell(ix,iz,iy)%neighcell(ineigh)%p => cell(neigh(1), neigh(2), neigh(3))
                  end do
               end do
            end do

         end do
      end do
   end do



end subroutine InitCellList

pure, function pcellro(ro), result(icell)
   real(8), intent(in)  :: ro(3)
   type(cell_type), pointer :: icell
   integer(4)  :: i(3)

   i = min(ncell,floor(ro*ircell))
   icell => cell(i(1), i(2), i(3))
end function pcellro

subroutine addiptocell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), intent(inout) :: icell

   icell%npart = icell%npart + 1
   icell%ip(call%npart) = ip
   icellip(ip)%p => icell

end subroutine addiptocell

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
   implicit none
   use MolModule, only: ro
   integer(4), intent(in)  :: ip
   type(cell_type), pointer :: cellold
   type(cell_type), pointer :: cellnew

   cellold = cellip(ip)%p
   cellnew = pcellro(ro(1:3,ip))

   if(cellold%id .ne. cellnew%id) then !do not update cell if cell has not changed
      call RmIpFromCell(ip, cellold)
      call AddIpToCell(ip, cellnew)
   end if

subroutine SetCellList

   use MolModule, only: np
   implicit none

   integer(4)  :: ip

   cell(:)%npart = 0
   do ip = 1, np
      call AddiptoCell(ip, pcellro(ip)%p)
   end do

end subroutine SetCellList



!needed:
      !call SetVList
      !call VListAver(iStage)
         !if (lllist) call UTwoBodyALList
            !call SetLList(rcut+drnlist)
            !call LListAver(iStage)
!       !             !              !----- DUTwoBody(A/ALList/P)New
end module CellList
