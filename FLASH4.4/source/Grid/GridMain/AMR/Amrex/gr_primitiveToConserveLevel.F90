#include "constants.h"

subroutine gr_primitiveToConserveLevel(level)
  use gr_interface,      ONLY : gr_getBlkIterator, gr_releaseBlkIterator
  use gr_amrexInterface, ONLY : gr_primitiveToConserve
  use gr_iterator,       ONLY : gr_iterator_t
  use block_metadata,    ONLY : block_metadata_t

  implicit none

  integer, intent(IN) :: level

  type(gr_iterator_t)    :: itor
  type(block_metadata_t) :: blockDesc

  call gr_getBlkIterator(itor, nodetype=ALL_BLKS, level=level+1)
  do while (itor%is_valid())
    call itor%blkMetaData(blockDesc)
    call gr_primitiveToConserve(blockDesc)
    call itor%next()
  end do
  call gr_releaseBlkIterator(itor)

end subroutine gr_primitiveToConserveLevel

