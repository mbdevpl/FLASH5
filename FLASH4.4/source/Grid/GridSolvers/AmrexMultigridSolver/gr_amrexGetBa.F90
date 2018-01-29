!!****if* source/Grid/GridSolvers/AmrexMultigridSolver/gr_amrexGetBa
!!
!!  NAME 
!!
!! gr_amrexGetBa
!!
!!  SYNOPSIS
!!
!!  call gr_amrexGetBa()
!!
!!
!!  DESCRIPTION 
!! This routine returns the amrex_boxarray associated with the multifab unk
!! 
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  Solver settings used ?? DESCRIPTION??
!!
!!***

subroutine gr_amrexGetBa(ba)
      use amrex_boxarray_module,     ONLY : amrex_boxarray
      use gr_physicalMultifabs, ONLY : unk
      implicit none
      type(amrex_boxarray), intent(inout) :: ba(:)
      ba=unk%ba
      return
end subroutine gr_amrexGetBa
