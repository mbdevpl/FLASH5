!!****if* source/Grid/GridMain/paramesh/Grid_getLocalNumBlks
!!
!! NAME
!!  Grid_getLocalNumBlks
!!
!! SYNOPSIS
!!
!!  Grid_getLocalNumBlks(integer(OUT) :: numBlocks)
!!
!! DESCRIPTION
!!  Get the number of local blocks on a processor
!!
!! ARGUMENTS
!!  numBlocks : The number of blocks currently in use on myProcessor
!!
!! NOTES
!!  There should be a better (more efficient) way to do this!
!!  The current version iterates through a loop and counts interations.
!!  There should be (and maybe already is) a better way to get the
!!  information directly from AMReX.
!!
!!  The blocks counted include LEAF as well as covered blocks.
!!  This is contrary to the convention for block iterators, where
!!  only LEAF blocks are exposed in the public Grid_getLeafIterator
!!  interface. But it is consistent with the functionality of the
!!  PARAMESH version, and used in this way in the friendly IO unit.
!!
!!  An alternative version Grid_getLocalNumLeafBlks is coded below,
!!  but currently (2018-02-26) not yet made public via Grid_interface.
!!***


subroutine Grid_getLocalNumBlks(numBlocks)

  use gr_interface, ONLY : gr_getBlkIterator, gr_releaseBlkIterator
  use gr_iterator, ONLY : gr_iterator_t

  implicit none

  integer,intent(out) :: numBlocks

  type(gr_iterator_t)  :: itor
  integer :: lb

  lb = 0

  call gr_getBlkIterator(itor)
  do while (itor%is_valid())
     lb = lb + 1
     call itor%next()
  enddo
  call gr_releaseBlkIterator(itor)

  numBlocks = lb

  return
end subroutine Grid_getLocalNumBlks



subroutine Grid_getLocalNumLeafBlks(numBlocks)

  use Grid_interface, ONLY : Grid_getLeafIterator, Grid_releaseLeafIterator
  use leaf_iterator, ONLY : leaf_iterator_t

  implicit none

  integer,intent(out) :: numBlocks

  type(leaf_iterator_t)  :: itor
  integer :: lb

  lb = 0

  call Grid_getLeafIterator(itor)
  do while (itor%is_valid())
     lb = lb + 1
     call itor%next()
  enddo
  call Grid_releaseLeafIterator(itor)

  numBlocks = lb

  return
end subroutine Grid_getLocalNumLeafBlks
