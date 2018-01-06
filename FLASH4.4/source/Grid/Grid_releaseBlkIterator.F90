!!****if* source/Grid/Grid_releaseBlkIterator
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
  use block_iterator, ONLY : block_iterator_t

  implicit none

  type(block_iterator_t), intent(INOUT) :: itor

  return
end subroutine Grid_releaseBlkIterator

