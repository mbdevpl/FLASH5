!!****f* source/Grid/Grid_getSingleCellCoords_lev
!!
!! NAME
!!  Grid_getSingleCellCoords_lev
!!
!! SYNOPSIS
!!
!!  call Grid_getSingleCellCoords_lev(integer(in)       ::  ind(MDIM), 
!!                           integer(in)       ::  level,
!!                           integer(in)       ::  edge,
!!                           real(OUT)         ::  coords(MDIM))
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
!! NOTES
!!
!!  The specific routine Grid_getSingleCellCoords_lev is also available under the
!!  generic name Grid_getSingleCellCoords.
!!
!!***

subroutine Grid_getSingleCellCoords_lev(ind, level,edge, coords)

  implicit none

#include "constants.h"

  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: level, edge
  real, dimension(MDIM), intent(out) :: coords

  return
end subroutine Grid_getSingleCellCoords_lev
