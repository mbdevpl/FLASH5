!!****if* source/Grid/GridMain/AMR/gr_releaseBlkIterator
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
!!  gr_getBlkIterator.F90
!!
!!***

subroutine gr_releaseBlkIterator(itor)
  use gr_iterator, ONLY : gr_iterator_t, destroy_iterator

  implicit none

  type(gr_iterator_t), intent(INOUT) :: itor

  ! Destroy explicitly ONLY for compilers that do not implement destructors
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
  call destroy_iterator(itor)
#endif
end subroutine gr_releaseBlkIterator

