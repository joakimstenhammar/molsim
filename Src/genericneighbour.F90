module GenericNeighbourModule

   use CellListModule, only: cell_type
   implicit none
   private
   public SetupNeighs, jpneighloc

   !saved variables for the cell list
   type(cell_type), pointer   :: icell
   integer(4)  :: incell, jpold
   integer(4)  :: ipcenter

   contains

   subroutine SetupNeighs(ip, nneigh, posip)
      use MolModule, only: lclist
      use MolModule, only: nneighpn
      use MolModule, only: myid, nproc
      use CellListModule, only: pcellro, cellip
      implicit none
      integer(4), intent(in)  :: ip
      integer(4), intent(out) :: nneigh
      real(8), optional, intent(in)  :: posip(3)
      integer(4) :: incelltmp

      ipcenter = ip
      if(lclist) then
         if(present(posip)) then
            icell => pcellro(posip(1:3))
         else
            icell => cellip(ip)%p
         end if
         nneigh = 0
         do incelltmp = 1 + myid, icell%nneighcell, nproc
            !(increment with nproc to have parallel execution)
            nneigh = nneigh + icell%neighcell(incelltmp)%p%npart
         end do
         nneigh = nneigh - 1 !exclude particle ip
      else
         nneigh = nneighpn(ip)
      end if

   end subroutine SetupNeighs

   function jpneighloc(jploc) result(jp)

      use MolModule, only: lclist
      use MolModule, only: jpnlist
      use MolModule, only: myid, nproc
      use CellListModule, only: ipnext
      implicit none
      integer(4), intent(in)  :: jploc
      integer(4) :: jp

      if(lclist) then
         !either choose the first head, or the next particle
         if(jploc == 1) then
            incell = 1 + myid
            jp = icell%neighcell(incell)%p%iphead
         else
            jp = ipnext(jpold)
         end if
         ! skip particle ipcenter
         if(jp == ipcenter) then
            jp = ipnext(jp)
         end if

         ! if at tail
         if(jp == 0) then
            !change to next cell
            do while (incell .le. icell%nneighcell - nproc)
               incell = incell + nproc
               jp = icell%neighcell(incell)%p%iphead
               !if a new head is found
               if (jp .ne. 0) then
                  exit
               end if
            end do
         end if
         jpold = jp
      else
         jp = jpnlist(jploc, ipcenter)
      end if
   end function

end module GenericNeighbourModule
