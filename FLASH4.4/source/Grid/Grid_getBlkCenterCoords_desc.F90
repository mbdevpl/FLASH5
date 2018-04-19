!!****f* source/Grid/Grid_getBlkCenterCoords
!!
!! NAME
!!  Grid_getBlkCenterCoords_desc
!!
!! SYNOPSIS
!!  Grid_getBlkCenterCoords_desc(integer(IN) :: blockDesc
!!                          real(OUT)   :: blockCenter(MDIM))
!!  
!! DESCRIPTION 
!!   Gets the coordinates of the center of the block identified by
!!   blockDesc.  Returns the coordinates in an array blockCenter
!!
!! ARGUMENTS
!!  blockDesc - block_metadata_t of the block. for UG always 1
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

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkCenterCoords_desc(blockDesc, blockCenter)
  use block_metadata,        ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN)  :: blockDesc
  real,    intent(OUT) :: blockCenter(MDIM)
  blockCenter=0.0
end subroutine Grid_getBlkCenterCoords_desc
