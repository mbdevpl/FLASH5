!!****if* source/Grid/GridMain/AMR/Grid_releaseLeafIterator
!!
!! NAME
!!  Grid_releaseLeafIterator
!!
!! SYNOPSIS
!!  Grid_releaseLeafIterator(leaf_iterator_t(INOUT) :: itor)
!!  
!! DESCRIPTION 
!!  Destroy given leaf block iterator.
!!
!! ARGUMENTS 
!!  itor - the block iterator to destroy.
!!
!! SEE ALSO
!!  Grid_getLeafIterator
!!
!!***

subroutine Grid_releaseLeafIterator(itor)
  use leaf_iterator, ONLY : leaf_iterator_t, destroy_iterator

  implicit none

  type(leaf_iterator_t), intent(INOUT) :: itor

  call destroy_iterator(itor)
end subroutine Grid_releaseLeafIterator

