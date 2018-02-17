!!****if* source/Grid/GridMain/AMR/paramesh/gr_getBlkIterator
!!
!! NAME
!!  gr_getBlkIterator
!!
!! SYNOPSIS
!!  gr_getBlkIterator(gr_iterator_t(OUT)    :: itor,
!!                    integer(IN), optional :: nodetype,
!!                    integer(IN), optional :: level,
!!                    logical(IN), optional :: tiling)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of all blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  gr_releaseBlkIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block iterator
!!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).  If
!!             no block class is given, then the iterator defaults to iterating
!!             over all blocks.
!!  level    - iterate only over all blocks/tiles located at this level of
!!             refinement.
!!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
!!             walk across all associated blocks on a tile-by-tile basis *if*
!!             the implementation supports this feature.  If a value is not
!!             given, is FALSE, or the implementation does not support tiling,
!!             the iterator will iterate on a block-by-block basis.
!!
!! NOTES
!!  gr_releaseBlkIterator.F90
!!  constants.h
!!
!!***

#include "constants.h"

subroutine gr_getBlkIterator(itor, nodetype, level, tiling)
  use gr_iterator, ONLY : gr_iterator_t, build_iterator

  implicit none

  type(gr_iterator_t), intent(OUT)          :: itor
  integer,             intent(IN), optional :: nodetype
  integer,             intent(IN), optional :: level
  logical,             intent(IN), optional :: tiling

  integer :: ntype

  if (present(nodetype)) then
    ntype = nodetype
  else
    ntype = ALL_BLKS
  end if

  call build_iterator(itor, ntype, level, tiling)
end subroutine gr_getBlkIterator

