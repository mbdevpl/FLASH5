!!****if* source/Grid/GridMain/AMR/paramesh/AmrexTransition/Grid_getSingleCellCoords_lev
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
#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
#define DEBUG_GRID

#include "constants.h"
#include "Flash.h"

subroutine Grid_getSingleCellCoords_lev(ind, level,edge, coords)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_getDeltas
  use Grid_data, ONLY : gr_globalDomain

  use Grid_data, ONLY : maxRefine => gr_maxRefine

  implicit none



  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: level, edge
  real, dimension(MDIM), intent(out) :: coords

  integer :: axis
  integer :: stride
  real    :: first
  real    :: fineDeltas(MDIM)

#ifdef DEBUG_GRID
!!$  print*,' inside Grid_getSingleCellCoords_lev', ind, level, edge
#ifdef FLASH_GRID_UG
  if((level < 1)) then
     print*,"Grid_getSingleCellCoords_lev :invalid level "
     call Driver_abortFlash("Grid_getSingleCellCoords_lev :invalid level ")
  end if
#endif
  if((edge/=LEFT_EDGE).and.(edge/=CENTER).and.(edge/=RIGHT_EDGE))&
       call Driver_abortFlash('Grid_getSingleCellCoods : invalid edge')

!  print*, 'leaving the DEBUG_GRID statement'
#endif

  ! DEVNOTE: We wouldn't need the following if we had gr_minCellSizes initialized... - KW
  call Grid_getDeltas(maxRefine, fineDeltas)

  stride = 2**(maxRefine - level)


  first = 0
  if (edge==CENTER) then
     first = stride * 0.5
  else if (edge==RIGHT_EDGE) then
     first = stride
  end if

  coords(:) = 0.0
  do axis = 1,NDIM
     coords(axis) = gr_globalDomain(LOW,axis) + (first+(ind(axis)-1)*stride) * fineDeltas(axis)
  end do

end subroutine Grid_getSingleCellCoords_lev









