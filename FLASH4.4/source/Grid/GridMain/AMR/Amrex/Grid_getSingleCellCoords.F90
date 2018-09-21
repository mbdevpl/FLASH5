!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getSingleCellCoords
!!
!! NAME
!!  Grid_getSingleCellCoords
!!
!! SYNOPSIS
!!
!!  Grid_getSingleCellCoords(integer(in)       ::  ind(MDIM), 
!!                           integer(in)       ::  blockid,
!!                           integer(in)       ::  edge,
!!                           integer(in)       ::  beginCount,
!!                           real(OUT)         ::  coords(MDIM))
!!  
!! DESCRIPTION 
!!
!!  Returns the coordinates of a single cell of a given block.
!!
!!  Note that radial coordinates obtained for guardcells beyond the r=0
!!  boundary will be negative.
!!
!! ARGUMENTS
!!
!!  ind - array holding the indices of the cell whose coordinates to return
!! 
!!  blockid - local block number
!!
!!  edge - indicates if user wants the left, center or right edge of the block
!!         options are LEFT_EDGE, RIGHT_EDGE, or CENTER.  These constants are
!!         defined in constants.h
!!
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!
!!    
!!  coords : returned coordinates of the specificed cell
!!
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_getSingleCellCoords(ind, blockId,edge, beginCount,coords)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer,dimension(MDIM), intent(in) :: ind
  integer, intent(in) :: blockId, edge
  integer, intent(in) :: beginCount
  real, dimension(MDIM), intent(out) :: coords

  coords(:) = 0.0
  CALL Driver_abortFlash("[Grid_getSingleCellCoords] AMReX does *not* deal in block IDs")
end subroutine Grid_getSingleCellCoords

subroutine Grid_getSingleCellCoords_Itor(ind, block, edge, beginCount, coords)
  use amrex_amrcore_module,  ONLY : amrex_geom
  use amrex_geometry_module, ONLY : amrex_problo

  use Driver_interface,      ONLY : Driver_abortFlash
  use Grid_interface,        ONLY : Grid_getGeometry
  use block_metadata,        ONLY : block_metadata_t

  implicit none
  
  type(block_metadata_t), intent(in) :: block
  integer, intent(in)  :: ind(MDIM)
  integer, intent(in)  :: edge
  integer, intent(in)  :: beginCount
  real,    intent(out) :: coords(MDIM)

  real    :: ind_t(1:MDIM)
  real    :: shift
  integer :: axis
  integer :: geometry

#ifdef DEBUG_GRID
  print*,' inside Grid_getSingleCellCoords', ind, edge, beginCount
  if(beginCount == EXTERIOR) then
     if((ind(IAXIS)<GRID_ILO_GC).or.(ind(IAXIS)>GRID_IHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<GRID_JLO_GC).or.(ind(JAXIS)>GRID_JHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<GRID_KLO_GC).or.(ind(KAXIS)>GRID_KHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else if(beginCount == INTERIOR) then
     if((ind(IAXIS)<1).or.(ind(IAXIS)>NXB))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<1).or.(ind(JAXIS)>NYB))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<1).or.(ind(KAXIS)>NZB))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else
     coords(:) = 0.0
     call Driver_abortFlash("Grid_getSingleCellCoords, incorrect value for beginCount")
  end if
#endif

  associate (x0       => amrex_problo, &
             x_blk_lo => block%limits(LOW, :), &
             dx       => amrex_geom(block%level - 1)%dx)
    ! x_blk_lo is 1-based cell-index of lower-left cell in block 

    if      (edge == LEFT_EDGE) then
      shift = 0.0
    else if (edge == CENTER) then
      shift = 0.5
    else if (edge == RIGHT_EDGE) then
      shift = 1.0
    else
      coords(:) = 0.0
      call Driver_abortFlash('[Grid_getSingleCellCoods] invalid edge')
    end if

    ! Translate indices to 1-based global indices adjusted according to edge
    ind_t = (ind + x_blk_lo - 1) + shift
    if (beginCount == EXTERIOR) then
      ind_t = ind_t - NGUARD
    end if

    coords(:) = 0.0
    coords(1:NDIM) = x0(1:NDIM) + (ind_t(1:NDIM) - 1) * dx(1:NDIM)
  end associate

end subroutine Grid_getSingleCellCoords_Itor

