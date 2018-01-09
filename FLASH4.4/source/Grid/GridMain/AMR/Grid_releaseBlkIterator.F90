!!****if* source/Grid/GridMain/AMR/Grid_releaseBlkIterator
!!
!! NAME
!!  Grid_releaseBlkIterator
!!
!! SYNOPSIS
!!  Grid_releaseBlkIterator(block_iterator_t(INOUT) :: itor)
!!  
!! DESCRIPTION 
!!  Destroy given block iterator.
!!
!! ARGUMENTS 
!!  itor - the block iterator to destroy.
!!
!! SEE ALSO
!!  Grid_getBlkIterator
!!
!!***

subroutine Grid_releaseBlkIterator(itor)
  use block_iterator, ONLY : block_iterator_t, destroy_iterator

  implicit none

  type(block_iterator_t), intent(INOUT) :: itor

  call destroy_iterator(itor)
end subroutine Grid_releaseBlkIterator

