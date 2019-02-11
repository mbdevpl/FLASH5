
!!****if* source/Grid/GridMain/paramesh/Grid_getBlkCenterCoords
!!
!! NAME
!!  Grid_getBlkCenterCoords
!!
!! SYNOPSIS
!!
!!  Grid_getBlkCenterCoords(type(block_metadata_t)(IN) :: block,
!!                          real(OUT)   :: blockCenter(MDIM))
!!
!!
!! DESCRIPTION 
!!   Gets the coordinates of the center of the block identified by
!!   blockId.  Returns the coordinates in an array blockCenter
!!
!!
!! ARGUMENTS
!!  block - block metadata
!!  blockCenter - returned array of size MDIM holding the blockCenter coords
!!
!! Example
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
!!  then the values returned in blockCenter are 
!!  blockCenter(IAXIS) = 0.0
!!  blockCenter(JAXIS) = 0.5
!!  blockCenter(KAXIS) = 0.0 since the dimension is not included  
!!
!!***


subroutine Grid_getBlkCenterCoords(block, blockCenter)

  use tree, ONLY : coord
  use block_metadata, ONLY : block_metadata_t
  implicit none

#include "constants.h"

  type(block_metadata_t),intent(in) :: block
  real,dimension(MDIM),intent(out) :: blockCenter
  blockCenter=coord(:,block%id)

end subroutine Grid_getBlkCenterCoords
