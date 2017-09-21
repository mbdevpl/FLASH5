!!****if* source/Grid/GridMain/paramesh/AmrexTransition/Grid_getBlkBoundBox_desc
!!
!! NAME
!!  Grid_getBlkBoundBox_desc
!!
!! SYNOPSIS
!! 
!! 
!!  call Grid_getBlkBoundBox_desc(integer(IN)  :: blockDesc,
!!                         real(OUT) :: boundBox(2, MDIM))
!!  
!! DESCRIPTION 
!!
!!  Gets the physical domain bounding box of the block identified 
!!  by blockDesc.  For each dimension the left (lower or forward) 
!!  physical coordinate of the block edge and the right (upper or back) 
!!  physical coordinate of the block edge is returned.  See arguments
!!  below for more detail.
!!
!! ARGUMENTS
!!
!!  blockDesc -local block number
!!
!!  boundBox - returned array holding the boundBox coordinates in
!!             each dimension
!!
!!            for readability, in constants.h we define IAXIS = 1, JAXIS = 2, KAXIS = 3
!!
!!            boundBox(1,IAXIS) = left edge coordinate of block in x direction
!!            boundBox(2,IAXIS) = right edge coordinate of block in x direction
!!            boundBox(1,JAXIS) = top edge coordinate of block in y direction
!!            boundBox(2,JAXIS) = bottom edge coordinate of block in y direction
!!            boundBox(1,KAXIS) = front edge coordinate of block in z direction
!!            boundBox(2,KAXIS) = back edge coordinate of block in z direction
!!
!! EXAMPLE
!!  
!!   In 2 dimensions, if physical coordinates are ...
!!    
!!     ________________(0.5 1.0)
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |                |
!!    |_______________ |
!!  (-0.5, 0.0)
!!
!!
!!
!!     boundBox(1, IAXIS) = -0.5
!!     boundBox(2, IAXIS) = 0.5
!!     boundBox(1, JAXIS) = 0.0
!!     boundBox(2, JAXIS) = 1.0
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!     boundBox(1, KAXIS) = 1 !returned as 1 because only 2 dims
!!
!!
!!***


subroutine Grid_getBlkBoundBox_desc(blockDesc, boundBox)

  use tree, ONLY : bnd_box
  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata, ONLY : block_metadata_t
  use Grid_interface, ONLY : Grid_getSingleCellCoords_lev
  implicit none

#include "constants.h"


  type(block_metadata_t), intent(in) :: blockDesc
  real,dimension(2,MDIM),intent(out) :: boundBox


  integer, dimension(LOW:HIGH,MDIM) :: lim
  integer :: blockId,level

  blockid = blockDesc%id

  if(blockID > 0) then
     boundBox = bnd_box(:,:,blockId)
  else
!!$     print *, "blockId = ", blockID
!!$     call Driver_abortFlash("Error: Grid_getBlkBoundBox_desc, blockId out of bounds")

     lim=blockDesc%limits
     level=blockDesc%level
     call Grid_getSingleCellCoords_lev(lim(LOW,:), level,LEFT_EDGE, boundBox(LOW,:))
     call Grid_getSingleCellCoords_lev(lim(HIGH,:), level,RIGHT_EDGE, boundBox(HIGH,:))
  end if

  return
end subroutine Grid_getBlkBoundBox_desc
