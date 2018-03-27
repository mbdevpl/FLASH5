!!****if* source/Grid/localAPI/gr_blockMatch
!!
!! NAME
!!
!!  Grid_blockMatch
!!
!! SYNOPSIS
!!
!!  ans = gr_blockMatch(integer(IN)  :: blkid,
!!                      integer(IN)  :: ntype,
!!                      integer(IN),OPTIONAL  :: refinementlevel)
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
!!              For the paramesh implementation, valid values are :
!!
!!              ALL_BLKS    all local blocks on a processor
!!
!!              IBDRY_BLKS  blocks that are on physical boundary along IAXIS
!!              JBDRY_BLKS  blocks that are on physical boundary along JAXIS
!!              KBDRY_BLKS  blocks that are on physical boundary along KAXIS
!!              ANY_BDRY_BLKS  blocks that have any of their faces
!!                          on a physical boundary.
!!              ACTIVE_BLKS all currently active blocks; in the paramesh
!!                          context this means parent and leaf blocks.
!!
!!              values that have meaning only for paramesh are :
!!              LEAF, PARENT_BLK or ANCESTOR  representing
!!                          the type of node in the oct-tree managing the blocks.
!!              REFINEMENT refinement level, the optional refinementlevel
!!                                           argument must be present and valid.
!!
!!              All of these constants are defined in constants.h
!!
!!   refinementlevel : If valid, refinement level to be used as match criterion.
!!                     Values greater that 0 are considered valid.
!!
!!                     If ntype is REFINEMENT, the caller must provide a present
!!                     and valid refinementlevel, and all blocks at that level
!!                     are considered matching.
!!                     For other choices of ntype, the refinementlevel is used
!!                     as additional criterion if present and valid.
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

