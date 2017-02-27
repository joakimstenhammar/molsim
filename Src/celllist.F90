module CellList

   implicit none
   private
   public

!needed:
         if (lllist) call UTwoBodyALList
            call SetLList(rcut+drnlist)
            call LListAver(iStage)
!       !             !              !----- DUTwoBody(A/ALList/P)New
end module CellList
