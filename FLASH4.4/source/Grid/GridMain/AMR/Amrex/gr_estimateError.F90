!!****if* source/Grid/GridMain/AMR/Amrex/gr_estimateError
!!
!! NAME
!!  gr_estimateError
!!
!! SYNOPSIS
!!
!!  gr_estimateError(real(INOUT) :: error(MAXBLOCKS),
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!  For each block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.  Update the
!!  corresponding value in the error array to be the maximum of the incoming
!!  value and the value calculated here.
!!
!!  ARGUMENTS 
!!
!!    error - indexed by block IDs.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

subroutine gr_estimateError(error, iref, refine_filter)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN)    :: iref
  real,    intent(IN)    :: refine_filter
  real,    intent(INOUT) :: error(MAXBLOCKS)
 
  call Driver_abortFlash("[gr_estimateError] Not implemented with AMReX")
end subroutine gr_estimateError

