!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getCellCoords
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
!!
!!    Coordinates are retrieved one axis at a time, 
!!    meaning you can get the i, j, _or_ k coordinates with one call.  
!!    If you want all the coordinates, all axes, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!    The code carries coordinates at cell centers as well as faces.
!!    It is possible to get coordinates for CENTER, only LEFT_EDGE,
!!    only RIGHT_EDGE or for all FACES along a dimension.
!!
!!    Note that radial coordinates obtained for guardcells beyond the r=0
!!    boundary will be negative.
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
!!      integer :: coordSize
!!      integer :: xCoord(coordSize) !sized to be number of coords returned
!!      
!!          !holds the number of cells returned in idir
!!          coordSize = blkLimitsGC(HIGH, IAXIS)+1
!!          call Grid_getCellCoords(IAXIS, cid, stride, FACES, .true., xCoord, coordSize) 
!!
!!     end do    
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getCellCoords(axis, block, edge, guardcell, coordinates, size)
  use amrex_amrcore_module,  ONLY : amrex_geom
  use amrex_geometry_module, ONLY : amrex_problo

  use Driver_interface,      ONLY : Driver_abortFlash
  use Grid_interface,        ONLY : Grid_getGeometry
  use block_metadata,        ONLY : block_metadata_t

#include "constants.h"
#include "Flash.h"

  implicit none

  integer,                intent(in)  :: axis
  type(block_metadata_t), intent(in)  :: block
  integer,                intent(in)  :: edge
  logical,                intent(in)  :: guardcell
  integer,                intent(in)  :: size
  real,                   intent(out) :: coordinates(size)

  real    :: shift
  real    :: x
  integer :: i

  integer :: geometry

#ifdef DEBUG_GRID
  print*,' get coordinates', axis, edge, guardcell, size
  ! DEV: TODO Implement for AMReX
!  if((blockID<1).or.(blockID>MAXBLOCKS)) then
!     call Driver_abortFlash("Grid_getCellCoords :invalid blockID ")
!  end if

  if((axis/=IAXIS) .and. (axis/=JAXIS) .and. (axis/=KAXIS)) then
     call Driver_abortFlash("Get Coords : invalid axis, must be IAXIS, JAXIS or KAXIS ")
  end if
#endif
  
  associate (x0   => amrex_problo(axis), &
             x_lo => block%limits(LOW, axis), &
             dx   => amrex_geom(block%level - 1)%dx(axis))
    ! x_lo is 1-based cell-index of lower-left cell in block 

    if      (edge == LEFT_EDGE) then
      shift = 0.0
    else if (edge == CENTER) then
      shift = 0.5
    else if (edge == RIGHT_EDGE) then
      shift = 1.0
    else if (edge == FACES) then
      shift = 0.0
    else
      call Driver_abortFlash('[Grid_getCellCoods] invalid edge')
    end if

    ! First index should be given in 1-based global indices adjusted according
    ! to edge.  However, make 0-based to simplify calculation in loop.
    x = x_lo + shift - 1.0
    if (guardcell) then
      x = x - NGUARD
    end if

    coordinates(:) = 0.0
    do i = 1, size
        coordinates(i) = x0 + x*dx
        x = x + 1.0
    end do
  end associate

end subroutine Grid_getCellCoords

