module CellListModule

implicit none
private
public UpdateCellip, InitCellList, SetCellList, CellListAver, TestCellList
public pcellro, cell_type, cell_pointer_array, cellip, ipnext

integer(4), parameter  :: maxneighcell = 27
type cell_pointer_array
   type(cell_type), pointer  :: p => null()
end type cell_pointer_array

type cell_type
   integer(4)  :: id       ! fore easy recognition
   integer(4)  :: npart
   integer(4)  :: nneighcell
   type(cell_pointer_array) :: neighcell(maxneighcell)
   integer(4)  :: iphead
end type cell_type

integer(4), allocatable :: ipnext(:)
integer(4), allocatable :: ipprev(:)

type(cell_type), target, allocatable  :: cell(:,:,:)    !cells
type(cell_pointer_array), allocatable  :: cellip(:) !cell of each particle
integer(4)  :: ncell(3)                !number of cells in x y z in each octant
real(8)     :: ircell(3)                  !inverse edge length of each cell

contains

subroutine InitCellList(rcut, iStage)

   use MolModule, only: ltrace, ltime, uout
   use MolModule, only: lbcbox, boxlen
   use MolModule, only: np

   real(8), intent(in)  :: rcut
   integer(4), intent(in)  :: iStage
   character(40), parameter :: txroutine ='InitCellList'
   integer(4)  :: ix, iy, iz, neigh(3), ineigh, id, jneigh
   logical  :: alreadyneighbour
   integer(4) :: directions(3,maxneighcell), directionindex(maxneighcell), idir !

   if (ltrace) call WriteTrace(1, txroutine, iStage)
   if (ltime) call CpuAdd('start', txroutine, 1, uout)

   if(rcut .le. 0.0d0) then
      call Stop(txroutine, 'rcut .le. 0.0', uout)
   end if
   if(.not. lbcbox) then
      call Stop(txroutine, 'celllist needs cubic box (lbcbox should be true)', uout)
   end if
   ncell(1:3) = max(1,floor(boxlen(1:3)/rcut)) !floor to underestimate the number of cells, therefore the cellsize is >= rcut
   ! but at least one cell in each direction is needed

   !boxlen2/ncell is the cell-size (larger than rcut)
   ircell(1:3) = ncell(1:3)/boxlen(1:3) !ircell is the inverse of the cell-size (smaller than 1/rcut)

   if(allocated(cell)) then
      deallocate(cell)
   end if
   if(allocated(cellip)) then
      deallocate(cellip)
   end if
   allocate(cell(0:(ncell(1)-1),0:(ncell(2)-1),0:(ncell(3)-1)))
   allocate(cellip(1:np))
   allocate(ipnext(1:np))
   allocate(ipprev(1:np))
   ipnext = 0
   ipprev = 0
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
   id = 0
   do ix = lbound(cell,dim=1), ubound(cell,dim=1)
      do iy = lbound(cell,dim=2), ubound(cell,dim=2)
         do iz = lbound(cell,dim=3), ubound(cell,dim=3)
            id = id + 1
            cell(ix,iy,iz)%id = id ! sets the id
            !allocate(cell(ix,iy,iz)%ip(np)) ! allocate cell particles
            cell(ix,iy,iz)%iphead = 0
            cell(ix,iy,iz)%npart = 0        ! initialize cell particles
         end do
      end do
   end do

   !generate all possible directions of neighbouring cells
   idir = 0
   do ix = -1, 1
      do iy = -1, 1
         do iz = -1, 1
            idir = idir + 1
            directions(1:3,idir) = (/ix, iy, iz/)
         end do
      end do
   end do

   ! sort the directions to have the neighbours with the smallest distance at the beginning
   ! therefore the hard core overlaps occur as early as possible
   call HeapSortIndex(maxneighcell, real(sum(abs(directions),dim=1),kind=8), directionindex)

   !set the neighbours
   do ix = lbound(cell,dim=1), ubound(cell,dim=1)
      do iy = lbound(cell,dim=2), ubound(cell,dim=2)
         do iz = lbound(cell,dim=3), ubound(cell,dim=3)
            !make pointers to neighbouring cells
            ineigh = 0
            do idir = 1, maxneighcell
               neigh(1:3) = (/ ix , iy , iz /) + directions(1:3,directionindex(idir))
               where (neigh > ubound(cell))
                  neigh = lbound(cell)
               elsewhere (neigh < lbound(cell))
                  neigh = ubound(cell)
               end where

               alreadyneighbour = .false.
               do jneigh = 1, ineigh
                  if (cell(ix,iy,iz)%neighcell(jneigh)%p%id .eq. cell(neigh(1), neigh(2), neigh(3))%id) then !neighbour is already set
                     alreadyneighbour = .true.
                     exit
                  end if
               end do

               if(.not. alreadyneighbour) then
                  ineigh = ineigh + 1
                  cell(ix,iy,iz)%nneighcell = ineigh
                  cell(ix,iy,iz)%neighcell(ineigh)%p => cell(neigh(1), neigh(2), neigh(3))
               end if
            end do
         end do
      end do
   end do

end subroutine InitCellList

pure function pcellro(ro) result(icell)
   use MolModule, only: boxlen2
   implicit none
   real(8), intent(in)  :: ro(3)
   type(cell_type), pointer :: icell
   integer(4)  :: i(3)

   i = floor((ro+ boxlen2)*ircell)
   icell => cell(i(1), i(2), i(3))
end function pcellro

subroutine AddIpToCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), target, intent(inout) :: icell
   integer(4)  :: jp

   icell%npart = icell%npart + 1
   jp = icell%iphead
   ipprev(jp) = ip
   icell%iphead = ip
   ipnext(ip) = jp
   ipprev(ip) = 0
   cellip(ip)%p => icell

end subroutine AddIpToCell

subroutine RmIpFromCell(ip, icell)
   implicit none
   integer(4), intent(in)  :: ip
   type(cell_type), intent(inout) :: icell
   integer(4)  :: nextp, prevp

   icell%npart = icell%npart - 1
   nextp = ipnext(ip)
   prevp = ipprev(ip)
   if( prevp .eq. 0) then ! ip is at head
      !make next particle to head
      icell%iphead = nextp
      ipprev(nextp) = 0
   else
      ! connext next and previous particle
      ipprev(nextp) = prevp
      ipnext(prevp) = nextp
   end if
   ipprev(ip) = 0
   ipnext(ip) = 0
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

subroutine TestCellList(output)
   use MolModule, only: ro, np
   use MolModule, only: rcut2
   implicit none
   integer(4), intent(in) :: output
   character(80), parameter :: txheading ='cell list testing'
   character(40), parameter :: txroutine ='InitCellList'
   integer(4)  :: ix, iy, iz, ineigh
   type(cell_type), pointer   :: icell, tncell
   integer(4)  :: ip, jp, incell, jploc, jpneigh
   real(8)  :: dr(3), r2
   logical  :: lipjpneighbour

   call WriteHead(2, txheading, output)
   write(output,'()')
   write(output,'(a,i0)') "Number of cells ", size(cell)
   write(output,'()')
   write(output,'(tr2,a49)') &
   repeat('-',49)
   write(output,'(tr2,a15,tr2,a15,tr2,a15)') &
   'cell id', 'i neighbour', 'id of neighbour'
   write(output,'(tr2,a15,tr2,a15,tr2,a15)') &
   repeat('-',15), repeat('-', 15), repeat('-',15)
   do ix = lbound(cell,dim=1), ubound(cell,dim=1)
      do iy = lbound(cell,dim=2), ubound(cell,dim=2)
         do iz = lbound(cell,dim=3), ubound(cell,dim=3)
            icell => cell(ix,iy,iz)
            do ineigh = 1, icell%nneighcell
               write(output,'(tr2,i15,tr2,i15,tr2,i15)') icell%id, ineigh, icell%neighcell(ineigh)%p%id
            end do
         end do
      end do
   end do
   write(output,'(tr2,a49)') &
   repeat('-',49)
   write(output,'()')
   write(output,'(tr2,a49)') &
   repeat('-',49)
   do ip = 1, np
      do jp = 1, np
         dr = ro(1:3,ip)-ro(1:3,jp)
         call PBCr2(dr(1), dr(2), dr(3), r2)
         if(r2 .le. rcut2) then !check if ip and jp are neighbours
            lipjpneighbour = .false.
            icell => cellip(ip)%p
            do incell = 1, icell%nneighcell
               tncell => icell%neighcell(incell)%p
               jpneigh = tncell%iphead
               do jploc = 1, tncell%npart
                  if(jpneigh .eq. jp)) then
                     lipjpneighbour = .true.
                     exit
                  end if
                  jpneigh = ipnext(jpneigh)
               end do
               if(lipjpneighbour) then
                  exit
               end if
            end do
            if(.not. lipjpneighbour) then
               write(output, *) "Error ip ", ip, " and jp ", jp , "are not in neighbouring cells! cell(ip)%id ", cellip(ip)%p%id, " cell(jp)%id) ",cellip(jp)%p%id
               write(output, *) "ro(ip): ",ro(1:3,ip), " ro(jp) ",ro(1:3,jp)
                  call Stop(txroutine, 'found two particles which should be neighbours but which are not', output)
            else
               write(output, *) "ip ", ip, " and jp ", jp , "are neighbours"
            end if
         end if
      end do
   end do
   write(output,'(tr2,a49)') &
   repeat('-',49)
   write(output,'()')
   write(output,'(tr2,a49)') &
   repeat('-',49)
   write(output,'(tr2,a15,tr2,a15,tr2,a15,tr2,a15)') &
   'particle id', 'cell id', 'prev particle id', 'next particle id'
   write(output,'(tr2,a15,tr2,a15,tr2,a15,tr2,a15)') &
   repeat('-',15), repeat('-', 15), repeat('-',15), repeat('-',15)
   do ip = 1, np
      write(output,'(tr2,i15,tr2,i15,tr2,i15,tr2,i15)') ip, cellip(ip)%p%id, ipprev(ip), ipnext(ip)
   end do
   write(output,'(tr2,a49)') &
   repeat('-',49)

end subroutine


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
            do icell = 1, cellip(ip)%p%nneighcell
               tmpcell => cellip(ip)%p%neighcell(icell)%p
               nneighip = nneighip + tmpcell%npart
            end do
            nneigh(1) = nneigh(1) + nneighip
            nneigh(2) = nneigh(2) + nneighip**2
            nneigh(3) = min(nneigh(3), nneighip)
            nneigh(4) = max(nneigh(4), nneighip)
         end do

      case (iAfterSimulation)

         nppp(1)   = nppp(1)/(nstep1*size(cell))
         nppp(2)   = sqrt(nppp(2)/(nstep1*size(cell))-nppp(1)**2)
         nneigh(1) = nneigh(1)/(nstep1*np)
         nneigh(2) = sqrt(nneigh(2)/(nstep1*np)-nneigh(1)**2)

         call WriteHead(2, txheading, uout)
         write(uout,'()')
         write(uout,'(a,i0)') "Number of cells ", size(cell)
         write(uout,'()')
         write(uout,'(a,3f14.4)') "Size of Cells in x,y,z", 1.0d0/ircell(1:3)
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
