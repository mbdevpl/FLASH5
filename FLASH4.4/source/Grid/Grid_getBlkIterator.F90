!!****if* source/Grid/Grid_getBlkIterator
!!
!! NAME
!!  Grid_getBlkIterator
!!
!! SYNOPSIS
!!  Grid_getBlkIterator(block_iterator_t(OUT) :: itor,
!!                      integer(IN)           :: nodetype,
!!                      integer(IN), optional :: level,
!!                      logical(IN), optional :: tiling)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  Grid_releaseBlkIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block iterator
!!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
!!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
!!             iterate only over blocks/tiles located at this level of
!!             refinement.
!!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
!!             walk across all associated blocks on a tile-by-tile basis *if*
!!             the implementation supports this feature.  If a value is not
!!             given, is FALSE, or the implementation does not support tiling,
!!             the iterator will iterate on a block-by-block basis.
!!
!! NOTES
!!  Grid_releaseBlkIterator.F90
!!  constants.h
!!
!!***

subroutine Grid_getBlkIterator(itor, nodetype, level, tiling)
  use block_iterator, ONLY : block_iterator_t

  implicit none

  type(block_iterator_t), intent(OUT)          :: itor
  integer,                intent(IN)           :: nodetype
  integer,                intent(IN), optional :: level
  logical,                intent(IN), optional :: tiling

  return
end subroutine Grid_getBlkIterator

