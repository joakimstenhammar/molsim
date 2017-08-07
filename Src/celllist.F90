!************************************************************************
!*                                                                      *
!*     CellList                                                         *
!*                                                                      *
!************************************************************************

! ... a new improved implementation of the linked list

! In each cell we store the particles in a linked list. The first in that cell is the head of the cell (given by iphead). Each
! particle is linked to its next and previous particle in the list by ipnext and ipprev. The last particle in the list is the tail.
! ipprev of the head is 0 and ipnext of the tail is 0. Therefore we cann loop through all particles of the linked list, by starting
! with handling iphead, and proceeding by handling the next particle of the particle just considered. When we reach particle 0, we
! have considered all particles in the list.
! See Allen, M. P. & Tildesley, D. J., Computer Simulation of Liquids, Oxford University Press, USA, 1987, 149.

module CellListModule

implicit none
private
public UpdateCellIp, InitCellList, SetCellList, CellListAver, TestCellList
public pcellro, cell_type, cellip, ipnext

integer(4)                               :: maxneighcell ! maximum number of neighbouring cells

type cell_pointer_array
   type(cell_type), pointer              :: p => null()
end type cell_pointer_array

type cell_type
   integer(4)                            :: id           ! for easy recognition
   integer(4)                            :: npart
   integer(4)                            :: nneighcell
   type(cell_pointer_array), allocatable :: neighcell(:)
   integer(4)                            :: iphead
end type cell_type

integer(4), allocatable                  :: ipnext(:)
integer(4), allocatable                  :: ipprev(:)

type(cell_type), target, allocatable     :: cell(:,:,:)             ! cells
type(cell_pointer_array), allocatable    :: cellip(:)               ! cell of each particle
integer(4)                               :: ncell(3) = 0            ! number of cells in x y z in each octant
real(8)                                  :: cellSize(3) = 0.0       ! inverse edge length of each cell
real(8)                                  :: cellSizei(3) = 0.0      ! inverse edge length of each cell

contains

subroutine InitCellList(rcell, iStage)

   use MolModule, only: ltrace, ltime, uout
   use MolModule, only: lbcbox, boxlen, lPBC, dpbc
   use MolModule, only: np
   use MolModule, only: rcut, rcut2

   real(8),    intent(in)                :: rcell
   integer(4), intent(in)                :: iStage
   character(40), parameter              :: txroutine ='InitCellList'
   integer(4),               allocatable :: directions(:,:) ,directionindex(:), tmpidneigh(:)
   type(cell_pointer_array), allocatable :: icellid(:)
   type(cell_type), pointer              :: icell
   integer(4)                            :: idir
   integer(4)                            :: ix, iy, iz, neigh(3), ineigh, id, ixdir, iydir, izdir
   integer(4)                            :: minneighcell
   real(8)                               :: r2, dr(3)

   if (ltrace) call WriteTrace(1, txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(rcell .le. 0.0d0) then
      call Stop(txroutine, 'rcell .le. 0.0', uout)
   end if
   if(.not. lbcbox) then
      call Stop(txroutine, 'CellList needs cubic box (lbcbox should be true)', uout)
   end if

   ncell(1:3) = max((/1, 1, 1/), floor(boxlen(1:3)/rcell)) !floor to underestimate the number of cells
   ! therefore the cellSize is >= rcell
   ! underestimation as when rcell = rcut one wants to have the cells larger than rcut
   ! but at least one cell in each direction is needed

   cellSize(1:3) = boxlen(1:3)/ncell(1:3) ! cellSize is the the cell size (larger than rcell)
   cellSizei(1:3) = 1.0d0/cellSize(1:3)   ! cellSizei is the inverse of the cell size (smaller than 1/rcell)

   call AllocateCellStructure(ncell, np)

   ! cell structure when ncell(1:3) = 4
   !        +-----+-----+-----+-----+
   !       /.., 3/     /     /     /|
   !      +-----+-----+-----+-----+ |
   !     /.., 2/     /     /     /| +
   !    +-----+-----+-----+-----+ |/|
   !   /.., 1/     /     /     /| + |
   !  +-----+-----+-----+-----+ |/| +
   ! /.., 0/     /     /     /| + |/|
   !+-----+-----+-----+-----+ |/| + |
   !| 0, 3| 1, 3| 2, 3| 3, 3| + |/| +
   !| , 0 | , 0 | , 0 | , 0 |/| + |/|
   !+-----+-----+-----+-----+ |/| + |
   !| 0, 2| 1, 2| 2, 2| 3, 2| + |/| +
   !| , 0 | , 0 | , 0 | , 0 |/| + |/
   !+-----+-----+-----+-----+ |/| +
   !| 0, 1| 1, 1| 2, 1| 3, 1| + |/
   !| , 0 | , 0 | , 0 | , 0 |/| +
   !+-----+-----+-----+-----+ |/
   !| 0, 0| 1, 0| 2, 0| 3, 0| +
   !| , 0 | , 0 | , 0 |  0  |/
   !+-----+-----+-----+-----+

   !allocate cells and set the id of the cells
   allocate(icellid(product(ncell(1:3))))
   id = 0
   do ix = 0, ncell(1) - 1
      do iy = 0, ncell(2) - 1
         do iz =  0, ncell(3) - 1
            id = id + 1
            icellid(id)%p        => cell(ix,iy,iz)
            icellid(id)%p%id     = id ! sets the id
            icellid(id)%p%iphead = 0
            icellid(id)%p%npart  = 0        ! initialize cell particles
         end do
      end do
   end do

   ! get the directions

   ! the maximum possible number neighbouring cells in any direction is ceiling(rcut/cellsize)
   ! as the cells have neighbours in both positive and negative direction we have multiply the number of cells by two
   ! in addition the central cell is also part if the neighbouring cells (+ 1)
   ! therefore we have 2*ceiling(rcut/cellsize) + 1 cells in each direction
   ! (when the cellsize is larger than rcut (drnlist .ge. 0.0), we have 3 cells in each direction)
   ! the total number of neighbouring cells is the product the number of cells in each direction
   maxneighcell = product(2*ceiling(rcut*cellSizei(1:3))+1)

   allocate(directions(3,maxneighcell))
   allocate(directionindex(maxneighcell))
   allocate(tmpidneigh(maxneighcell))

   ! loop over all possible neighbouring positions
   idir = 0
   do ixdir = -ceiling(rcut*cellSizei(1)), ceiling(rcut*cellSizei(1))
      do iydir = -ceiling(rcut*cellSizei(2)), ceiling(rcut*cellSizei(2))
         do izdir = -ceiling(rcut*cellSizei(3)), ceiling(rcut*cellSizei(3))

            dr(1:3) = max((/0, 0, 0/),abs((/ixdir, iydir, izdir/))-1)*cellSize(1:3) !distance to closest part of cell

            if(lPBC) then
               call PBCr2(dr(1), dr(2), dr(3),r2)
            else
               r2 = sum(dr(1:3)**2)
            end if

            if(r2 > rcut2) then ! distance of the cells so large, that it is not needed
               cycle
            end if

            idir = idir + 1
            ! store direction as one which is a neighbouring cell
            directions(1:3,idir) = (/ixdir, iydir, izdir/)

         end do
      end do
   end do

   maxneighcell = idir

   ! sort the directions to have the neighbours with the smallest distance at the beginning
   ! therefore the hard core overlaps occur as early as possible
   call HeapSortIndex(maxneighcell, real(sum(abs(directions(1:3,1:maxneighcell)),dim=1),kind=8), directionindex(1:maxneighcell))

   minneighcell = huge(minneighcell)

   !set the neighbours
   do ix = 0, ncell(1) - 1
      do iy = 0, ncell(2) - 1
         do iz =  0, ncell(3) - 1
            icell => cell(ix,iy,iz)

            !get cell neighbors
            ineigh = 0
            icell%nneighcell = 0
            tmpidneigh = 0
            do idir = 1, maxneighcell
               neigh(1:3) = (/ ix, iy, iz/) + directions(1:3, directionindex(idir))
               if(lPBC) then
                  ! apply periodic boundary conditions in directions where dpbc equals the boxlen,
                  ! to apply only the right periodic boundary conditions
                  where (((neigh >= ncell) .or. (neigh < 0)) .and. (dpbc == boxlen))
                     neigh = modulo(neigh,ncell)
                  end where
               end if
               if(any(neigh(1:3) < 0) .or. any(neigh(1:3) >= ncell)) then !neighbour is out of bounds
                  cycle
               end if

               if(all(tmpidneigh(1:ineigh) /= cell(neigh(1), neigh(2), neigh(3))%id)) then
                  ! the neighbouring cell is not already a neighbour
                  ineigh = ineigh + 1
                  tmpidneigh(ineigh) = cell(neigh(1), neigh(2), neigh(3))%id
               end if
            end do
            icell%nneighcell = ineigh

            minneighcell = min(minneighcell, ineigh)

            !allocate memory for the current cell
            if(allocated(icell%neighcell)) then
               deallocate(icell%neighcell)
            end if
            allocate(icell%neighcell(icell%nneighcell))

            ! now assign neighbours
            do ineigh = 1, icell%nneighcell
               icell%neighcell(ineigh)%p => icellid(tmpidneigh(ineigh))%p
            end do
         end do
      end do
   end do

   if(minneighcell < nproc) then
      call Warn(txroutine, 'the number of neighbouring cells is smaller than the number of processors, parallel execution
      will be very inefficient', uout)
   end if

   deallocate(directions, directionindex, tmpidneigh, icellid)
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine InitCellList

subroutine AllocateCellStructure(ncell, np)
   implicit none
   integer(4), intent(in) :: ncell(3)
   integer(4), intent(in) :: np

   if(allocated(cell)) then
      deallocate(cell)
   end if
   if(allocated(cellip)) then
      deallocate(cellip)
   end if
   if(allocated(ipnext)) then
      deallocate(ipnext)
   end if
   if(allocated(ipprev)) then
      deallocate(ipprev)
   end if
   allocate(cell(0:(ncell(1)-1),0:(ncell(2)-1),0:(ncell(3)-1)))
   allocate(cellip(1:np))
   allocate(ipnext(1:np))
   allocate(ipprev(1:np))
   ipnext = 0
   ipprev = 0
end subroutine

function pcellro(ro) result(icell)
   use MolModule, only: boxlen2
   implicit none
   real(8), intent(in)  :: ro(3)
   type(cell_type), pointer :: icell
   integer(4)  :: i(3)

   i = floor((ro+ boxlen2)*cellSizei)
   icell => cell(i(1), i(2), i(3))
end function pcellro

subroutine AddIpToCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), target, intent(inout) :: icell
   integer(4)  :: jp

   icell%npart = icell%npart + 1
   jp = icell%iphead  ! current head
   icell%iphead = ip  ! make ip the head
   if(jp .ne. 0) then ! when a head is present
      ipprev(jp) = ip ! move the head one down
      ipnext(ip) = jp ! make old head particle next to ip
   else               ! new particle is alone in cell
      ipnext(ip) = 0
   end if
   ipprev(ip) = 0     ! there is no particle before the head

   cellip(ip)%p => icell !associate particle with cell

end subroutine AddIpToCell

subroutine RmIpFromCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), intent(inout) :: icell
   integer(4)  :: nextp, prevp

   icell%npart = icell%npart - 1
   nextp = ipnext(ip)
   prevp = ipprev(ip)
   if (nextp .eq. 0) then ! ip is at the tail
      if( icell%iphead .eq. ip ) then ! ip is at head
         ! no particles are left
         icell%iphead = 0
      else
         ! make the previous particle the tail
         ipnext(prevp) = 0
      end if
   else ! ip is not at the tail
      if( icell%iphead .eq. ip ) then ! ip is at head
         !make next particle the head
         icell%iphead = nextp
         ipprev(nextp) = 0
      else ! ip is neither at the tail nor at the head
         ! connect next and previous particle
         ipprev(nextp) = prevp
         ipnext(prevp) = nextp
      end if
   end if

   !reset ip
   ipprev(ip) = 0
   ipnext(ip) = 0
   cellip(ip)%p => null()

end subroutine RmIpFromCell

subroutine UpdateCellIp(ip)
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

end subroutine UpdateCellIp

subroutine SetCellList

   use MolModule, only: np, ro
   use MolModule, only: ltime, uout
   implicit none

   character(40), parameter :: txroutine ='SetCellList'
   integer(4)               :: ip
   type(cell_type), pointer :: celltmp

   if (ltime) call CpuAdd('start', txroutine, 1, uout)
   cell(:,:,:)%npart = 0
   do ip = 1, np
      celltmp => pcellro(ro(1:3,ip))
      call AddIpToCell(ip, celltmp)
   end do
   if (ltime) call CpuAdd('stop', txroutine, 1, uout)

end subroutine SetCellList

subroutine TestCellList(output)
   use MolModule, only: ro, np
   use MolModule, only: rcut2
   implicit none
   integer(4), intent(in) :: output
   character(80), parameter :: txheading ='cell list testing'
   character(40), parameter :: txroutine = 'TestCellList'
   integer(4)  :: ix, iy, iz, ineigh
   type(cell_type), pointer   :: icell, locncell
   integer(4)  :: ip, jp, incell, jploc, jpneigh
   real(8)  :: dr(3), r2
   logical  :: lipjpneighbour

   call WriteHead(2, txheading, output)
   write(output,'()')
   write(output,'(a,i0)') "Number of cells ", size(cell)
   write(output,'()')
   write(output,'(tr2,a49)') repeat('-',49)
   write(output,'(3(tr2,a15))') 'cell id', 'ith neighbour', 'id of neighbour'
   write(output,'(3(tr2,a15))') repeat('-',15), repeat('-', 15), repeat('-',15)
   do ix = lbound(cell,dim=1), ubound(cell,dim=1)
      do iy = lbound(cell,dim=2), ubound(cell,dim=2)
         do iz = lbound(cell,dim=3), ubound(cell,dim=3)
            icell => cell(ix,iy,iz)
            do ineigh = 1, icell%nneighcell
               write(output,'(3(tr2,i15))') icell%id, ineigh, icell%neighcell(ineigh)%p%id
            end do
         end do
      end do
   end do
   write(output,'(tr2,a49)') repeat('-',49)
   write(output,'()')
   write(output,'(tr2,a49)') repeat('-',49)
   do ip = 1, np
      do jp = 1, np
         dr = ro(1:3,ip)-ro(1:3,jp)
         call PBCr2(dr(1), dr(2), dr(3), r2)
         if(r2 .le. rcut2) then !check if ip and jp are neighbours
            lipjpneighbour = .false.
            icell => cellip(ip)%p
            neighbourcells: do incell = 1, icell%nneighcell
               locncell => icell%neighcell(incell)%p
               jpneigh = locncell%iphead
               do jploc = 1, locncell%npart
                  if(jpneigh .eq. jp) then
                     lipjpneighbour = .true.
                     exit neighbourcells
                  end if
                  jpneigh = ipnext(jpneigh)
               end do
            end do neighbourcells
            if(.not. lipjpneighbour) then
               write(output, *) "Error ip ", ip, " and jp ", jp , "are not in neighbouring cells!"
               write(output, *) "cell(ip)%id: ", cellip(ip)%p%id, "cell(jp)%id): ",cellip(jp)%p%id
               write(output, *) "ro(ip): ",ro(1:3,ip), " ro(jp) ",ro(1:3,jp)
               call Stop(txroutine, 'found two particles which should be neighbours but which are not', output)
            else
               continue
            end if
         end if
      end do
   end do
   write(output,'(tr2,a49)') repeat('-',49)
   write(output,'()')
   write(output,'(tr2,a49)') repeat('-',49)
   write(output,'(4(tr2,a15))') 'particle id', 'cell id', 'prev particle id', 'next particle id'
   write(output,'(4(tr2,a15))') repeat('-',15), repeat('-', 15), repeat('-',15), repeat('-',15)
   do ip = 1, np
      write(output,'(4(tr2,i15))') ip, cellip(ip)%p%id, ipprev(ip), ipnext(ip)
   end do
   write(output,'(tr2,a49)') repeat('-',49)

end subroutine

subroutine CellListAver(iStage)

   use MolModule, only: uout
   use MolModule, only: np
   use MolModule, only: master, nstep1
   use MolModule, only: iWriteInput, iAfterMacrostep, iAfterSimulation
   implicit none
   integer(4), intent(in)   :: iStage
   character(80), parameter :: txheading ='cell list statistics'

   ! variables to average the cell statistics
   ! 1: average ; 2: standard deviation ; 3: minimum ; 4: maximum
   real(8), save            :: npPerCell(4)     ! number of particles
   real(8), save            :: nNeighPerPart(4) ! number of neigbours per particle

   type(cell_type), pointer :: tmpcell
   integer(4)               :: ip, icell
   real(8)                  :: nNeighPerPartIp


   if (master) then
      select case (iStage)
      case (iWriteInput)

         npPerCell(1:2)     = 0.0d0
         npPerCell(3)       = +huge(1.0d0)
         npPerCell(4)       = -huge(1.0d0)
         nNeighPerPart(1:2) = 0.0d0
         nNeighPerPart(3)   = +huge(1.0d0)
         nNeighPerPart(4)   = -huge(1.0d0)

      case (iAfterMacrostep)

         npPerCell(1) = npPerCell(1) + sum(cell(:,:,:)%npart)
         npPerCell(2) = npPerCell(2) + sum(cell(:,:,:)%npart**2)
         npPerCell(3) = min(npPerCell(3), real(minval(cell(:,:,:)%npart)))
         npPerCell(4) = max(npPerCell(4), real(maxval(cell(:,:,:)%npart)))

         do ip = 1, np
            nNeighPerPartIp = -1.0d0 !do not count the particle itself as a neighbour
            do icell = 1, cellip(ip)%p%nNeighPerPartcell
               tmpcell => cellip(ip)%p%neighcell(icell)%p
               nNeighPerPartIp = nNeighPerPartIp + tmpcell%npart
            end do
            nNeighPerPart(1) = nNeighPerPart(1) + nNeighPerPartIp
            nNeighPerPart(2) = nNeighPerPart(2) + nNeighPerPartIp**2
            nNeighPerPart(3) = min(nNeighPerPart(3), nNeighPerPartIp)
            nNeighPerPart(4) = max(nNeighPerPart(4), nNeighPerPartIp)
         end do

      case (iAfterSimulation)

         npPerCell(1)     = npPerCell(1)/(nstep1*size(cell))
         npPerCell(2)     = sqrt(npPerCell(2)/(nstep1*size(cell))-npPerCell(1)**2)
         nNeighPerPart(1) = nNeighPerPart(1)/(nstep1*np)
         nNeighPerPart(2) = sqrt(nNeighPerPart(2)/(nstep1*np)-nNeighPerPart(1)**2)

         call WriteHead(2, txheading, uout)
         write(uout,'()')
         write(uout,'(a,i0)') "Number of cells ", size(cell)
         write(uout,'()')
         write(uout,'(a,3f14.4)') "Size of cells in x,y,z", 1.0d0/cellSizei(1:3)
         write(uout,'()')
         write(uout,'(t12,a,t34,a)') &
         'no of part per cell', 'no of neighbours per part'
         write(uout,'(t12,a,t34,a)') &
         '-------------------', '----------------------------------'
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'average', npPerCell(1), nNeighPerPart(1)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'one sd ', npPerCell(2), nNeighPerPart(2)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'minimum', npPerCell(3), nNeighPerPart(3)
         write(uout,'(a,t10,f15.1,10x,f15.1,15x,2f15.1)') 'maximum', npPerCell(4), nNeighPerPart(4)

      end select
   end if

   end subroutine CellListAver
end module CellListModule
