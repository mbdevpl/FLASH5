!!****if* source/Grid/GridMain/AMR/Amrex/gr_releaseBlkIterator
!!
!! NAME
!!  gr_releaseBlkIterator
!!
!! SYNOPSIS
!!  gr_releaseBlkIterator(gr_iterator_t(INOUT) :: itor)
!!  
!! DESCRIPTION 
!!  Destroy given block iterator.
!!
!! ARGUMENTS 
!!  itor - the block iterator to destroy.
!!
!! SEE ALSO
!!  gr_getBlkIterator
!!
!!***

subroutine gr_releaseBlkIterator(itor)
  use gr_iterator, ONLY : gr_iterator_t, destroy_iterator

  implicit none

  type(gr_iterator_t), intent(INOUT) :: itor

  call destroy_iterator(itor)
end subroutine gr_releaseBlkIterator

