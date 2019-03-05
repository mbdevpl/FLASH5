!!****if* source/Grid/localAPI/gr_getInteriorBlkPtr
!!
!! NAME
!!  gr_getInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_getInteriorBlkPtr(integer(IN)            :: blockID,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to a single block of simulation data from the
!!  paramesh data structures when there are no permanent guard cells
!!
!!  When using Paramesh 4 in NO_PERMANENT_GUARDCELLS mode, it is important to
!!  release the block pointer for a block before getting it for another block.
!!  For example if pointer to block 1 is not yet released and the user
!!  tries to get a pointer to block 2, the routine will abort.
!!
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!
!!
!!
!! NOTES
!!
!!  gr_getInteriorBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

subroutine gr_getInteriorBlkPtr_blk(blockDesc, dataPtr, gridDataStruct)
  use block_metadata, ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(IN)         :: blockDesc
  real,                              pointer :: dataPtr(:,:,:,:)
  integer,                intent(in)         :: gridDataStruct
  
  return
end subroutine gr_getInteriorBlkPtr_blk

subroutine gr_getInteriorBlkPtr(tileDesc, dataPtr, gridDataStruct)
  use flash_tile, ONLY : flash_tile_t

  implicit none

  type(flash_tile_t), intent(IN)         :: tileDesc
  real,                          pointer :: dataPtr(:,:,:,:)
  integer,            intent(IN)         :: gridDataStruct
  
  return
end subroutine gr_getInteriorBlkPtr

