!!****if* source/Grid/GridMain/AMR/gr_estimateError
!!
!! NAME
!!  gr_estimateBlkError
!!
!! SYNOPSIS
!!
!!  call gr_estimateBlkError(real(INOUT) :: error,
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!
!!  DESCRIPTION
!!
!!  For one block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.
!!
!!  ARGUMENTS
!!
!!    error - the computed error, a scalar value for the current variable
!!            (given by iref) and current block.
!!
!!    blockDesc - describes the block.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically
!!
!!  NOTES
!!
!!    In the case of the PARAMESH Grid implementation, this routine is
!!    called from gr_estimateError.
!!
!!  SEE ALSO
!!
!!    Grid_markRefineDerefine
!!    Grid_markRefineDerefineCallback
!!
!!***

subroutine gr_estimateBlkError(error, tileDesc, iref, refine_filter)
  use flash_tile, ONLY : flash_tile_t

  implicit none

  integer, intent(IN) :: iref
  type(flash_tile_t),intent(IN) :: tileDesc
  real, intent(IN) ::  refine_filter
  real,intent(INOUT) :: error

end subroutine gr_estimateBlkError

