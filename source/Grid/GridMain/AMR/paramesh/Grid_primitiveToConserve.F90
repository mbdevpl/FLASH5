!!****if* source/Grid/GridMain/Grid_primitiveToConserve
!!
!! NAME
!!
!!  Grid_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  call Grid_primitiveToConserve(integer(in) :: blkList(count),
!!                           integer(in) :: count,
!!                           logical(in) :: force)
!!
!!
!! DESCRIPTION
!!
!!  Calls gr_primitiveToConserve
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   force - whether to force conversion
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_primitiveToConserve(blkList,count,force)
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator

  use Grid_data, ONLY: gr_convertToConsvdForMeshCalls

  implicit none
  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList 
  logical,intent(IN) :: force
  logical :: tempSwap

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  ! DEV: TODO This needs to be rethought or modernized to work with iterators
  tempSwap = (force .and. (.not.gr_convertToConsvdForMeshCalls))
  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .true.
  end if

  call Grid_getTileIterator(itor, ACTIVE_BLKS, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call gr_primitiveToConserve(tileDesc)

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  if (tempSwap) then
     gr_convertToConsvdForMeshCalls = .false.
  end if
end subroutine Grid_primitiveToConserve
