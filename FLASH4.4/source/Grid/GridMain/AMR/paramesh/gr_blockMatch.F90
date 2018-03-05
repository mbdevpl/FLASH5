!!****if* source/Grid/GridMain/paramesh/gr_blockMatch
!!
!! NAME
!!
!!  Grid_blockMatch
!!
!! SYNOPSIS
!!
!!  ans = Grid_blockMatch(integer  :: blkid,
!!                        integer  :: ntype,
!!                        integer,OPTIONAL  :: refinementlevel)
!!
!! DESCRIPTION
!!
!!  Test whether a block matches a criterion.
!!
!! ARGUMENTS
!!
!!   blkid : block ID
!!
!!   ntype : block type, or type of requested match.
!!
!!              For the paramesh implementation, valid values are :
!!
!!              ALL_BLKS    all local blocks on a processor
!!
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces
!!              on a physical boundary.
!!              ACTIVE_BLKS all currently active blocks, in paramesh
!!              context that means parent and leaf blocks
!!
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!              the type of node in the Oct-tree managing the blocks.
!!              REFINEMENT refinement level, optional argument must be present.
!!
!!              All of these constants are defined in constants.h
!!
!!   refinementlevel : level used as additional criterion.
!!
!!                     This is used as additional criterion if present.
!!
!! NOTES
!!
!!   Only implemented for the PARAMESH AMR Grid implementation
!!
!! RETURN TYPE
!!
!!  LOGICAL
!!
!! SEE ALSO
!!
!!   Grid_getListOfBlocks
!!
!!***

function gr_blockMatch(blkID,ntype,refinementLevel) result(match)

  use tree, ONLY : nodetype,lnblocks,neigh,lrefine
  use gr_specificData, ONLY : gr_oneBlock

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
  case (LEAF,PARENT_BLK,ANCESTOR)
     match = (nodetype(i)==ntype)
  case(ALL_BLKS)
     match = .TRUE.
  case (IBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,1,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,2,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case (JBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,3,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,4,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case(KBDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        isBnd = (neigh(1,5,i) .LE. PARAMESH_PHYSICAL_BOUNDARY)
        isBnd = isBnd.or.(neigh(1,6,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
        match = isBnd
     end if
  case(ANY_BDRY_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
     if (match) then
        match = ANY(neigh(1,1:6,i).LE.PARAMESH_PHYSICAL_BOUNDARY)
     end if
  case(ACTIVE_BLKS)
     match = ((nodetype(i)==LEAF).or.(nodetype(i)==PARENT_BLK))
  case(TRAVERSED)
     match = (gr_oneBlock(i)%blockType==TRAVERSED) !DEV: What does this mean?
  case(TRAVERSED_AND_ACTIVE)
     match = (gr_oneBlock(i)%blockType==TRAVERSED_AND_ACTIVE) !DEV: What does this mean?
  case(REFINEMENT)
     if (present(refinementLevel)) then
        match = .TRUE.          ! test for matching level is done below
     else
        call Driver_abortFlash("[Grid_getListofBlocks] with ntype REFINEMENT optional argument refinementlevel must be present")
     end if
  case default
     match = .FALSE.
     call Driver_abortFlash("[Grid_blockMatch] ntype argument not recognized")
  end select

  if (match .AND. present(refinementLevel)) then
     if (refinementLevel > 0) then
        match = (lrefine(i) == refinementLevel)
     end if
  end if

end function gr_blockMatch
