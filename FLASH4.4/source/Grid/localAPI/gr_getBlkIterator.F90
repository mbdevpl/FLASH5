!!****if* source/Grid/localAPI/gr_getBlkIterator
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
!!  Construct an iterator for walking across a specific subset of blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  gr_releaseBlkIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block iterator
!!  nodetype - the class of blocks to iterate over.  The options are
!!                - ALL_BLKS             - all blocks
!!                - LEAF                 - only leaf blocks
!!                - PARENT_BLK           - only those blocks with at least one
!!                                         child block that is a leaf
!!                - ACTIVE_BLKS          - only blocks that are leaves or 
!!                                         parents
!!                - ANCESTOR             - all blocks that are neither a leaf
!!                                         nor a parent
!!                - REFINEMENT           - all blocks on a given refinement
!!                                         level
!!                - IBDRY_BLKS           - only active blocks with at least one
!!                                         face on the I-axis boundaries
!!                - JBDRY_BLKS           - only active blocks with at least one 
!!                                         face on the J-axis boundaries
!!                - KBDRY_BLKS           - only active blocks with at least one
!!                                         face on the K-axis boundaries
!!                - ANY_BDRY_BLKS        - only those active blocks with at
!!                                         least one face on a boundary
!!                - TRAVERSED            - (paramesh)
!!                - TRAVERSED_AND_ACTIVE - (paramesh)
!!                - INREGION             - (paramesh) 
!!             The default value is ALL_BLKS.  Please see the documentation
!!             for a specific Grid implementation to see which options are 
!!             implemented.
!!  level    - iterate only over all blocks/tiles of the correct nodetype
!!             that are located at this level of refinement.  Note that the
!!             level value must be given with respect to FLASH's 1-based
!!             level index scheme.  If no level value is given, then
!!             iteration is not restricted to any level.
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

subroutine gr_getBlkIterator(itor, nodetype, level, tiling)
  use gr_iterator, ONLY : gr_iterator_t, build_iterator

  implicit none

  type(gr_iterator_t), intent(OUT)          :: itor
  integer,             intent(IN), optional :: nodetype
  integer,             intent(IN), optional :: level
  logical,             intent(IN), optional :: tiling

  return
end subroutine gr_getBlkIterator

