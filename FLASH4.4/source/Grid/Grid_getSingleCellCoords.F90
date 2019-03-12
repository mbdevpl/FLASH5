!!****f* source/Grid/Grid_getSingleCellCoords
!!
!! NAME
!!  Grid_getSingleCellCoords
!!
!! SYNOPSIS
!!
!!  call Grid_getSingleCellCoords(integer(in)       ::  ind(MDIM), 
!!                           integer(in)       ::  level,
!!                           integer(in)       ::  edge,
!!                           real(OUT)         ::  coords(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Returns the coordinates of a single cell of a given block.
!!
!!
!! ARGUMENTS
!!
!!  ind - array holding the indices of the cell whose coordinates to return.
!!        In this implementation, global indices (for the given level) are assumed.
!! 
!!  level - refinement level (1 based)
!!
!!  edge - indicates if user wants the left, center or right edge of the block
!!         options are LEFT_EDGE, RIGHT_EDGE, or CENTER.  These constants are
!!         defined in constants.h
!!
!!    
!!  coords : returned coordinates of the specificed cell
!!
!!***

#include "constants.h"

subroutine Grid_getSingleCellCoords(ind, level, edge, coords)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)  :: ind(1:MDIM)
  integer, intent(in)  :: level
  integer, intent(in)  :: edge
  real,    intent(out) :: coords(1:MDIM)

  coords(:) = 0.0

  call Driver_abortFlash("[Grid_getSingleCellCoords] DEPRECATED")
end subroutine Grid_getSingleCellCoords
