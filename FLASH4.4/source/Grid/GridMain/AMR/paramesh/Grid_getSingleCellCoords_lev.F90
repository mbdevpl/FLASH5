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
  use Grid_data, ONLY : gr_delta, gr_globalDomain

#ifdef FLASH_GRID_PARAMESH
  use tree,      ONLY : maxRefine => lrefine_max
#else
  use Grid_data, ONLY : maxRefine => gr_maxRefine
#endif

  implicit none



  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: level, edge
  real, dimension(MDIM), intent(out) :: coords

  integer :: axis
  integer :: stride
  real    :: first

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


!!$  do axis = 1,MDIM
!!$     do i = 1,gr_nBlockX-1
!!$        coords = ( gr_globalDomain(LOW,axis)*(gr_nBlockX-i)  + gr_globalDomain(HIGH,axis) * i ) / real(gr_nBlockX)
!!$     end do
!!$  end do

  stride = 2**(maxRefine - level)


  first = 0
  if (edge==CENTER) then
     first = stride * 0.5
  else if (edge==RIGHT_EDGE) then
     first = stride
  end if

  coords(:) = 0.0
  do axis = 1,NDIM
     coords(axis) = gr_globalDomain(LOW,axis) + (first+(ind(axis)-1)*stride) * gr_delta(axis,maxRefine)
  end do

end subroutine Grid_getSingleCellCoords_lev









