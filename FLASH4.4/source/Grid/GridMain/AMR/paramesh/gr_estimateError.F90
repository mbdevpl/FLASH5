!!****if* source/Grid/GridMain/AMR/gr_estimateError
!!
!! NAME
!!  gr_estimateError
!!
!! SYNOPSIS
!!
!!  gr_estimateError(real(INOUT) :: error(MAXBLOCKS),
!!                   integer(IN) :: iref,
!!                   real(IN)    :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!  For each block, estimate the error associated with the given variable to
!!  help determine if the block needs refinement or derefinement.  Update the
!!  corresponding value in the error array to be the maximum of the incoming
!!  value and the value calculated here.
!!
!!  ARGUMENTS 
!!
!!    error - indexed by block IDs.
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

#include "constants.h"  

subroutine gr_estimateError(error, iref, refine_filter)
  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use gr_interface,   ONLY : gr_estimateBlkError
  use Grid_iterator,  ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t

  implicit none

  integer, intent(IN)    :: iref
  real,    intent(IN)    :: refine_filter
  real,    intent(INOUT) :: error(MAXBLOCKS)
  
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

!==============================================================================

! A non-directional guardcell fill for CENTER (and also EOS calls for
! all block cells, including guardcells, if any refinement variables
! refine_var_# require this to be current) must have been performed
! when this routine is invoked. Moreover, there must not be any
! intervening calls that would modify the solution data in unk (at
! least, for the variables to be used for refinement criteria).
! Finally, this should be true for both LEAF and PARENT blocks
! (node types 1 and 2).
! Normally the caller of this routine, Grid_markRefineDerefine, takes care
! of all that.
!
! If this routine must be used in a situation where the conditions above
! are not true, the simplest (but probably inefficient) way of adapting
! this code to that situation would be uncommenting the following line:
!!$  call Grid_fillGuardCells(CENTER_FACES,ALLDIR)


! We are using more cell layers, including guardcells, from unk.

     
  !==============================================================================

  call Grid_getTileIterator(itor, ACTIVE_BLKS, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call gr_estimateBlkError(error(tileDesc%id), tileDesc, iref, refine_filter)

     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)
end subroutine gr_estimateError

