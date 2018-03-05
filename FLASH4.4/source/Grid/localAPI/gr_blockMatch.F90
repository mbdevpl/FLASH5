!!****if* source/Grid/localAPI/gr_blockMatch
!!
!! NAME
!!
!!  Grid_blockMatch
!!
!! SYNOPSIS
!!
!!  ans = gr_blockMatch(integer  :: blkid,
!!                      integer  :: ntype,
!!                      integer,OPTIONAL  :: refinementlevel)
!!
!! DESCRIPTION
!!
!!  Test whether a block matches a criterion.  Paramesh only.
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

function gr_blockMatch(blkID, ntype, refinementLevel) result(match)
  implicit none

  integer, intent(IN)           :: blkID
  integer, intent(IN)           :: ntype
  integer, intent(IN), OPTIONAL :: refinementLevel
  logical                       :: match

  match = .FALSE.
end function gr_blockMatch

