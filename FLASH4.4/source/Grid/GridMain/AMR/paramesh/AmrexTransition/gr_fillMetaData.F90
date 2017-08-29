#include "constants.h"

subroutine gr_fillMetaData(blockID, blockDesc)
  use Grid_interface
  use block_metadata, ONLY : block_metadata_t
  implicit none
  integer,intent(IN) :: blockID
  type(block_metadata_t),intent(OUT) :: blockDesc

  blockDesc%id = blockID
  blockDesc%grid_index = -1 !DEV: ??
  call Grid_getBlkRefineLevel(blockID, blockDesc%level)
  call Grid_getBlkIndexLimits(blockID, blockDesc%limits, blockDesc%limitsGC)

  ! localLimits should end up being the same as limits, and same for GC
  blockDesc%localLimits(LOW, :)   = blockDesc%limits(LOW, :)   - blockDesc%limitsGC(LOW, :) + 1
  blockDesc%localLimits(HIGH, :)  = blockDesc%limits(HIGH, :)  - blockDesc%limitsGC(LOW, :) + 1
  blockDesc%localLimitsGC(LOW, :) = 1
  blockDesc%localLimitsGC(HIGH, :)= blockDesc%limitsGC(HIGH, :)- blockDesc%limitsGC(LOW, :) + 1

  call Grid_getBlkPtr        (blockID, blockDesc%fp,     CENTER)

  call Grid_getBlkCornerID   (blockID, blockDesc%cid,    blockDesc%stride)
end subroutine gr_fillMetaData
