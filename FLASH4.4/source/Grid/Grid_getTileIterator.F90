!!****if* source/Grid/Grid_getTileIterator
!!
!! NAME
!!  Grid_getTileIterator
!!
!! SYNOPSIS
!!  Grid_getTileIterator(leaf_iterator_t(OUT)  :: itor,
!!                       integer(IN)           :: nodetype,
!!                       integer(IN), optional :: level,
!!                       logical(IN), optional :: tiling)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching leaf block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  Grid_releaseTileIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block/tile iterator
!!  nodetype - the class of blocks to iterate over.  The options for AMReX are
!!                - ALL_BLKS - all blocks
!!                - LEAF     - only leaf blocks
!!  level    - iterate only over leaf blocks/tiles located at this level of
!!             refinement.
!!             A level value of UNSPEC_LEVEL is equivalent to omitting
!!             this optional argument.
!!  tiling   - an optional optimization hint.  If TRUE, then the iterator will
!!             walk across all associated blocks on a tile-by-tile basis *if*
!!             the implementation supports this feature.  If a value is not
!!             given, is FALSE, or the implementation does not support tiling,
!!             the iterator will iterate on a block-by-block basis.
!!
!! NOTES
!!  Grid_releaseTileIterator.F90
!!
!!***

subroutine Grid_getTileIterator(itor, nodetype, level, tiling, tileSize)
  use flash_iterator, ONLY : flash_iterator_t

  implicit none

  type(flash_iterator_t), intent(OUT)          :: itor
  integer,                intent(IN)           :: nodetype
  integer,                intent(IN), optional :: level
  logical,                intent(IN), optional :: tiling
  integer,                intent(IN), optional :: tileSize(1:MDIM)

  return
end subroutine Grid_getTileIterator

