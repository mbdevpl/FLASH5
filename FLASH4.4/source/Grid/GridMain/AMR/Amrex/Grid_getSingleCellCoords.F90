!!****if* source/Grid/GridMain/Chombo/Grid_getSingleCellCoords
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

  CALL Driver_abortFlash("[Grid_getSingleCellCoords] AMReX does *not* deal in block IDs")
end subroutine Grid_getSingleCellCoords

subroutine Grid_getSingleCellCoords_Itor(ind, block, edge, beginCount, coords)
  use amrex_amrcore_module,  ONLY : amrex_geom
  use amrex_geometry_module, ONLY : amrex_problo

  use block_metadata,        ONLY : block_metadata_t

  implicit none
  
  type(block_metadata_t), intent(in) :: block
  integer, intent(in)  :: ind(MDIM)
  integer, intent(in)  :: edge
  integer, intent(in)  :: beginCount
  real,    intent(out) :: coords(MDIM)

  real    :: ind_t(1:MDIM) = 0.0d0
  real    :: shift = 0.0d0
  integer :: axis = 1

#ifdef DEBUG_GRID
  print*,' inside Grid_getSingleCellCoords', ind, edge, beginCount
  if ((edge /= LEFT_EDGE) .and. (edge /= CENTER) .and. (edge /= RIGHT_EDGE)) then
    call Driver_abortFlash('Grid_getSingleCellCoods : invalid edge')
  end if
  if(beginCount == EXTERIOR) then
     if((ind(IAXIS)<GRID_ILO_GC).or.(ind(IAXIS)>GRID_IHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<GRID_JLO_GC).or.(ind(JAXIS)>GRID_JHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<GRID_KLO_GC).or.(ind(KAXIS)>GRID_KHI_GC))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else if(beginCount == INTERIOR) then
     if((ind(IAXIS)<GRID_ILO).or.(ind(IAXIS)>GRID_IHI))&
          call Driver_abortFlash('GetSingleCellCoords : I index out of blkLimits')
     if((ind(JAXIS)<GRID_JLO).or.(ind(JAXIS)>GRID_JHI))&
          call Driver_abortFlash('GetSingleCellCoords : J index out of blkLimits')
     if((ind(KAXIS)<GRID_KLO).or.(ind(KAXIS)>GRID_KHI))&
          call Driver_abortFlash('GetSingleCellCoords : K index out of blkLimits')
  else
     call Driver_abortFlash("Grid_getSingleCellCoords, incorrect value for beginCount")
  end if
#endif

  associate (x0       => amrex_problo, &
             x_blk_lo => block%limits(LOW, :), &
             dx       => amrex_geom(block%level - 1)%dx)
    ! x_blk_lo is 1-based cell-index of lower-left cell in block 
    
    shift = 0.0d0
    if (edge == CENTER) then
      shift = 0.5d0
    else if (edge == RIGHT_EDGE) then
      shift = 1.0d0
    end if

    ! Translate indices to 1-based global indices adjusted according to edge
    ind_t = (ind + x_blk_lo - 1) + shift
    if (beginCount == EXTERIOR) then
      ind_t = ind_t - NGUARD
    end if

    coords(:) = 0.0d0
    coords(1:NDIM) = x0(1:NDIM) + (ind_t(1:NDIM) - 1) * dx(1:NDIM)
  end associate
 
end subroutine Grid_getSingleCellCoords_Itor

