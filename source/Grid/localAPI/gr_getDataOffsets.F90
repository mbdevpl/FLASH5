!!****if* source/Grid/GridMain/gr_getDataOffsets
!!
!! NAME
!!  gr_getDataOffsets
!!
!! SYNOPSIS
!!
!!  gr_getDataOffsets(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: startingPos(MDIM),
!!                    integer(IN) :: length(MDIM),
!!                    integer(IN) :: beginCount, 
!!                    integer(OUT):: begOffset(MDIM),
!!                    logical(OUT):: getBlkPtr )
!!  
!! DESCRIPTION 
!!  
!!  This routine determines the offset of the data withing the block from
!!  where to fetch. These offsets are used by data get/put routines to 
!!  correctly determine the starting index in all dimensions
!!
!! ARGUMENTS 
!!
!!  blockID   : my block number
!!
!!  gridDataStruct : integer value specifying the type of data desired.
!!
!!  startingPos(MDIM):
!!           specifies the starting position in each dimension of 
!!           the plane of data being fetched.
!!   
!!           startingPos(1) = i
!!           startingPos(2) = j
!!           startingPos(3) = k
!!
!!           If a problem is 2 dimensions startingPos(3) is irrelevant and
!!           ignored.  If a problem is only 1 dimension this routine doesn't
!!           make any sense and an error is returned.
!!
!! length : the length of the data in each dimension
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See examples
!!               in get/put data routines for more details.
!!               if GLOBALIDX1 is specified... !DEV : incomplete
!!
!! begOffset - the calculated offset values
!!
!! getBlkPtr - indicates to the data get/put routines whether blk pointer needs to 
!!          fetched. It is always false in permanent GC mode, and in UG. Depending
!!          upon the starting index, the block may or may not need to be fetched in
!!          the non permanent mode when using PARAMESH 4.
!!
!!***
subroutine gr_getDataOffsets(block,gridDataStruct,startingPos,length,beginCount,begOffset,getIntPtr)
  use block_metadata, ONLY : block_metadata_t

#include "constants.h"
#include "Flash.h"
  implicit none
  type(block_metadata_t), intent(IN)  :: block
  integer, intent(IN)  :: gridDataStruct
  integer, intent(IN)  :: beginCount
  integer,dimension(MDIM),intent(IN) :: startingPos
  integer,dimension(MDIM),intent(IN) :: length
  integer,dimension(MDIM),intent(OUT) :: begOffset
  logical,intent(OUT) :: getIntPtr

  begOffset(:) = 0
  getIntPtr = .FALSE.
end subroutine gr_getDataOffsets
