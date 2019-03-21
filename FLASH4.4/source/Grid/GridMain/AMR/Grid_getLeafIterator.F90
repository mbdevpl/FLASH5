!!****if* source/Grid/GridMain/AMR/Grid_getLeafIterator
!!
!! NAME
!!  Grid_getLeafIterator
!!
!! SYNOPSIS
!!  Grid_getLeafIterator(leaf_iterator_t(OUT)  :: itor,
!!                       integer(IN), optional :: level,
!!                       logical(IN), optional :: tiling)
!!  
!! DESCRIPTION 
!!  Construct an iterator for walking across a specific subset of leaf blocks or
!!  tiles within the current octree structure.  The iterator is already
!!  set to the first matching leaf block/tile.
!!
!!  Once finished, the iterator should be destroyed with
!!  Grid_releaseLeafIterator
!!
!! ARGUMENTS 
!!  itor     - the requested block iterator
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
!!  Grid_releaseLeafIterator.F90
!!  constants.h
!!
!!***

subroutine Grid_getLeafIterator(itor, level, tiling)
  use Driver_interface, ONLY : Driver_abortFlash
  use leaf_iterator,    ONLY : leaf_iterator_t, build_iterator

  implicit none

  type(leaf_iterator_t), intent(OUT)          :: itor
  integer,               intent(IN), optional :: level
  logical,               intent(IN), optional :: tiling

!  call Driver_abortFlash("[Grid_getLeafIterator] This iterator is DEPRECATED.  Use Grid_iterator_t.")

  call build_iterator(itor, level, .FALSE.)
end subroutine Grid_getLeafIterator

