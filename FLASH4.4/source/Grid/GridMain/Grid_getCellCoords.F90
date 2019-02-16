!!****if* source/Grid/GridMain/Grid_getCellCoords
!!
!! NAME
!!  Grid_getCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getCellCoords(integer(IN)  :: axis,
!!                      block_metadata_t(IN) :: block
!!                      integer(IN):: edge, 
!!                      logical(IN):: guardcell, 
!!                      real(OUT)  :: coordinates(size),
!!                      integer(IN):: size)
!!  
!!  
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block.
!!    Coordinates are retrieved one axis at a time, 
!!    meaning you can get the i, j, _or_ k coordinates with one call.  
!!    If you want all the coordinates, all axes, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!    The code carries coordinates at cell centers as well as faces.
!!    It is possible to get coordinates for CENTER, only LEFT_EDGE,
!!    only RIGHT_EDGE or for all FACES along a dimension.
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   block - derived type containing metadata of block of interest
!!
!!   edge - integer value with one of four values, 
!!          LEFT_EDGE, RIGHT_EDGE, CENTER or FACES
!!          The edge argument specifies what side of the zone to get, 
!!          the CENTER point, the LEFT_EDGE  or the RIGHT_EDGE of the zone.
!!          FACES gets the left and right face of each cell, but since 
!!          two adjacent cells have a common face, there are only N+1
!!          unique values if N is the number of cells.
!!
!!   guardcell - logical input. If true coordinates for guardcells are returned
!!          along with the interior cells, if false, only the interior coordinates 
!!          are returned.
!!
!!          
!!   coordinates - The array holding the data returning the coordinate values
!!                 coordinates must be at least as big as "size" (see below)
!!           
!!   size - integer specifying the size of the coordinates array.
!!          if edge = CENTER/LEFT_EDGE/RIGHT_EDGE then
!!                If guardcell true then size =  interior cells + 2*guardcells
!!                otherwise size = number of interior cells
!!          If edge=FACES 
!!                If guardcell true then size =  interior cells + 2*guardcells+1
!!                otherwise size = number of interior cells+1
!!
!!               
!!  EXAMPLE 
!!
!!  1. Getting cell centered values
!!
!!   #include "constants.h"
!!   #include "Flash.h"
!!
!!      
!!      integer :: coordSize
!!      integer :: xCoord(coordSize) !sized to be number of coords returned
!!      
!!      
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)
!!          call Grid_getCellCoords(IAXIS, cid, stride, CENTER, .true., xCoord, coordSize) 
!!
!!     end do    
!!
!!  2. Getting face values
!! 
!!   #include "constants.h"
!!   #include "Flash.h"
!!
!!      
!!      integer :: coordSize
!!      integer :: xCoord(coordSize) !sized to be number of coords returned
!!      
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)+1
!!          call Grid_getCellCoords(IAXIS, cid, stride, FACES, .true., xCoord, coordSize) 
!!
!!     end do    
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the FORTRAN module Grid_data. Variables that are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!***

#ifdef DEBUG
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_getCellCoords(axis, edge, level, lo, hi, coordinates)
  use Grid_interface,   ONLY : Grid_getDeltas
  use Grid_data,        ONLY : gr_globalDomain
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)  :: level
  integer, intent(in)  :: axis
  integer, intent(in)  :: edge
  integer, intent(in)  :: lo(1:MDIM)
  integer, intent(in)  :: hi(1:MDIM)
  real,    intent(out) :: coordinates(:)

  real    :: shift
  integer :: nElements
  real    :: deltas(1:MDIM)
  integer :: i

#ifdef DEBUG_GRID
  if (      (edge /= LEFT_EDGE) .AND. (edge /= RIGHT_EDGE) &
      .AND. (edge /= FACES)     .AND. (edge /= CENTER)) then
     call Driver_abortFlash("[Grid_getCellCoords] invalid edge specification, must be LEFT_EDGE &
          RIGHT_EDGE, or CENTER")
  end if

  if((axis /= IAXIS) .AND. (axis /= JAXIS) .AND. (axis /= KAXIS)) then
     call Driver_abortFlash("[Grid_getCellCoords] invalid axis, must be IAXIS, JAXIS or KAXIS ")
  end if
#endif
  
  ! Calculate number of cells
  nElements = hi(axis) - lo(axis) + 1
  if      (edge == FACES)  then
      shift = 2.0
      ! Get number of faces
      nElements = nElements + 1
  else if (edge == LEFT_EDGE) then
      shift = 2.0
  else if (edge == CENTER) then
      shift = 1.5
  else if (edge == RIGHT_EDGE) then
      shift = 1.0
  else
      call Driver_abortFlash('[Grid_getCellCoords] Invalid edge')
  end if

  if (SIZE(coordinates) < nElements) then
      call Driver_abortFlash("[Grid_getCellCoords] coordinates is too small")
  end if

  call Grid_getDeltas(level, deltas)

  associate (x0 => gr_globalDomain(LOW, axis), &
             dx => deltas(axis))
      do i = 1, nElements
          coordinates(i) = x0 + (lo(axis) + i - shift) * dx
      end do
  end associate
end subroutine Grid_getCellCoords

