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
!!  This is consistent with the functionality of the PARAMESH version, and used
!!  in this way in the friendly IO unit.
!!
!!  An alternative version Grid_getLocalNumLeafBlks is coded below,
!!  but currently (2018-02-26) not yet made public via Grid_interface.
!!***

#include "constants.h"

subroutine Grid_getLocalNumBlks(numBlocks)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use flash_iterator, ONLY : flash_iterator_t

  implicit none

  integer, intent(OUT) :: numBlocks

  type(flash_iterator_t) :: itor

  numBlocks = 0

  call Grid_getTileIterator(itor, ALL_BLKS, tiling=.FALSE.)
  do while (itor%isValid())
     numBlocks = numBlocks + 1
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_getLocalNumBlks

subroutine Grid_getLocalNumLeafBlks(numBlocks)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use flash_iterator, ONLY : flash_iterator_t

  implicit none

  integer, intent(OUT) :: numBlocks

  type(flash_iterator_t) :: itor

  numBlocks = 0

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while (itor%isValid())
     numBlocks = numBlocks + 1
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_getLocalNumLeafBlks

