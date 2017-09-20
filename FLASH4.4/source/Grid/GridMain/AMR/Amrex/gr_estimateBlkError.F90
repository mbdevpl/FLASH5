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

!!REORDER(4): solnData

subroutine gr_estimateBlkError(error, blockDesc, iref, refine_filter)
  use block_metadata, ONLY : block_metadata_t
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"  
#ifdef INDEXREORDER
  integer,parameter::IX=1,IY=2,IZ=3
#else
  integer,parameter::IX=2,IY=3,IZ=4
#endif  
  integer, intent(IN) :: iref
  type(block_metadata_t),intent(IN) :: blockDesc
  real, intent(IN) ::  refine_filter
  real,intent(INOUT) :: error

  call Driver_abortFlash("[gr_estimateBlkError] Not yet implemented for AMReX") 
end subroutine gr_estimateBlkError

