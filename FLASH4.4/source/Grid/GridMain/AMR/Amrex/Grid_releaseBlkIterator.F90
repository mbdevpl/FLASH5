!!****if* source/Grid/GridMain/AMR/Amrex/Grid_releaseBlkIterator
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
!! NOTES
!!  Grid_getBlkIterator.F90
!!  constants.h
!!
!!***

subroutine Grid_releaseBlkIterator(itor)
  use block_iterator, ONLY : block_iterator_t

  implicit none

  type(block_iterator_t), intent(INOUT) :: itor

  ! Destroy explicitly ONLY for compilers that do not implement destructors
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
  call itor%destroy_iterator()
#endif
end subroutine Grid_releaseBlkIterator

