!!****f* source/Grid/Grid_blockMatch
!!
!! NAME
!!
!!  Grid_blockMatch
!!
!! SYNOPSIS
!!
!!  match = Grid_blockMatch(integer(in)  :: blkid,
!!                          integer(in)  :: ntype,
!!                          integer(in),OPTIONAL :: refinementlevel)
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
!!              For all Grid implementations, valid values are :
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
!!              values that have meaning only for paramesh are :
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!              the type of node in the Oct-tree managing the blocks.
!!              REFINEMENT the refinement level
!!
!!              All of these constants are defined in constants.h
!!
!!   refinementlevel : level used as additional criterion.
!!
!!                     Only used if ntype = REFINEMENT.
!!
!! RETURN VALUE
!!
!!  LOGICAL match  - indicates whether the block matches the criterion
!!
!! SEE ALSO
!!
!!   Grid_getListOfBlocks
!!
!!***

function Grid_blockMatch(blkID,ntype,refinementLevel) result(match)

  implicit none

  logical               :: match
  integer,intent(in   ) :: blkID
  integer,intent(in   ) :: ntype
  integer,intent(in   ),OPTIONAL :: refinementLevel

  match = 0

end function Grid_blockMatch
