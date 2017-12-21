!!****if* source/Grid/GridMain/AMR/Grid_getBlkIterator
!!
!! NAME
!!  Grid_getBlkIterator
!!
!! SYNOPSIS
!!  Grid_getBlkIterator(block_iterator_t(OUT) :: itor
!!                      integer(IN)           :: nodetype,
!!                      level(IN), optional   :: level)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of blocks
!!  within the current mesh octree structure.  The iterator is already
!!  set to the first matching block.
!!
!!  Once finished, the iterator should be destroyed with
!!  Grid_releaseBlkIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block iterator
!!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS)
!!  level    - if nodetype is LEAF, PARENT, ANCESTOR, or REFINEMENT, then 
!!             iterate only over blocks located at this level of 
!!             octree structure
!!
!! NOTES
!!  Grid_releaseBlkIterator.F90
!!  constants.h
!!
!!***

subroutine Grid_getBlkIterator(itor, nodetype, level)
  use block_iterator, ONLY : block_iterator_t

  implicit none

  type(block_iterator_t), intent(OUT)          :: itor
  integer,                intent(IN)           :: nodetype
  integer,                intent(IN), optional :: level

  itor = block_iterator_t(nodetype, level)
end subroutine Grid_getBlkIterator

