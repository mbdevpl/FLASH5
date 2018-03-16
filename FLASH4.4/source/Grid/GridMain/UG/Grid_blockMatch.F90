!!****if* source/Grid/GridMain/UG/Grid_blockMatch
!!
!! NAME
!!
!!  Grid_blockMatch
!!
!! SYNOPSIS
!!
!!  ans = Grid_blockMatch(integer(in)  :: blkid,
!!                        integer(in)  :: ntype,
!!                        integer(in),OPTIONAL  :: refinementlevel)
!!
!! DESCRIPTION
!!
!!  Test whether a block matches a criterion.
!!
!! ARGUMENTS
!!
!!   blkid : block ID, should be 1 for UG
!!
!!   ntype : block type, or type of requested match.
!!
!!              For the UG Grid implementation, valid values are :
!!
!!              ALL_BLKS    all local blocks on a processor
!!
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces
!!                          on a physical boundary.
!!              ACTIVE_BLKS all currently active blocks (In the paramesh
!!                          context this means parent and leaf blocks.)
!!
!!              values that have meaning for paramesh are :
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!              the type of node in the oct-tree managing the blocks.
!!              REFINEMENT the refinement level
!!
!!              All of these constants are defined in constants.h
!!
!!   refinementlevel : level used as additional criterion.
!!
!!                     Only used if ntype = REFINEMENT.
!!
!! RETURN TYPE
!!
!!  LOGICAL
!!
!! NOTES
!!
!!  The only non-stub implementation provided is in the UG Grid implementation.
!!  For the paramesh Grid implementation, see gr_blockMatch instead.
!!
!! SEE ALSO
!!
!!   Grid_getListOfBlocks
!!
!!***

function Grid_blockMatch(blkID,ntype,refinementLevel) result(match)

  use Grid_data, ONLY : gr_blockType, lnblocks
  use Grid_data, ONLY : gr_axisMe, gr_axisNumProcs

  implicit none

#include "constants.h"
  logical               :: match
  integer,intent(in   ) :: blkID
  integer,intent(in   ) :: ntype
  integer,intent(in   ),OPTIONAL :: refinementLevel

  integer :: i
  logical :: isBnd

  i = blkID
  if (i > lnblocks) then
     match = .FALSE.
     return                   !block number is not valid, return immediately!
  end if

  select case (ntype)
  case (LEAF)
     match = .TRUE.
  case (PARENT_BLK,ANCESTOR)
     match = .FALSE.
  case(ALL_BLKS,ACTIVE_BLKS)
     match = .TRUE.
  case (IBDRY_BLKS)
     match = ((gr_axisMe(IAXIS)==0).or.(gr_axisMe(IAXIS)==gr_axisNumProcs(IAXIS)-1))
  case (JBDRY_BLKS)
     match = ((gr_axisMe(JAXIS)==0).or.(gr_axisMe(JAXIS)==gr_axisNumProcs(JAXIS)-1))
  case(KBDRY_BLKS)
     match = ((gr_axisMe(KAXIS)==0).or.(gr_axisMe(KAXIS)==gr_axisNumProcs(KAXIS)-1))
  case(ANY_BDRY_BLKS)
     match = ANY((gr_axisMe(IAXIS:IAXIS+NDIM-1)==0).or.(gr_axisMe(IAXIS:IAXIS+NDIM-1)==gr_axisNumProcs(IAXIS:IAXIS+NDIM-1)-1))
  case(TRAVERSED)
     match = (gr_blockType==TRAVERSED) !DEV: What does this mean?
  case(TRAVERSED_AND_ACTIVE)
     match = (gr_blockType==TRAVERSED_AND_ACTIVE) !DEV: What does this mean?
  case(REFINEMENT)
     if (present(refinementLevel)) then
        match = (refinementLevel .LE. 1)
     else
        call Driver_abortFlash("[Grid_getListofBlocks] with ntype REFINEMENT optional argument refinementlevel must be present")
     end if
  case default
     match = .FALSE.
     call Driver_abortFlash("[Grid_blockMatch] ntype argument not recognized")
  end select

end function Grid_blockMatch
