!!****if* source/Grid/GridMain/Grid_getBlkPtr
!!
!! NAME
!!  Grid_getBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPtr(type(block_metadta_t)(IN) :: block,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!  Gets a pointer to a single block of simulation data from the
!!  specified Grid data structure. The block includes guard cells.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
!!
!!  When using Paramesh 4 in NO_PERMANENT_GUARDCELLS mode, it is important to
!!  release the block pointer for a block before getting it for another block.
!!  For example if pointer to block 1 is not yet released and the user
!!  tries to get a pointer to block 2, the routine will abort.
!!
!! ARGUMENTS 
!!
!!   block - block metadata
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell centered variables (default)
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!
!!
!!
!! NOTES
!!
!!  Grid_getBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseBlkPtr when you are finished with it!
!!
!!***

!!REORDER(5): unk, facevar[xyz], scratch_ctr, scratch_facevar[xyz]
!!REORDER(4): dataPtr
!!FOR FUTURE: Add REORDER for unk, facevar[xyz]1, etc.?

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkPtr(block,dataPtr, gridDataStruct)
  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata, ONLY : block_metadata_t

  implicit none
  
  type(block_metadata_t), intent(in) :: block
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, optional,intent(in) :: gridDataStruct
  
  call Driver_abortFlash("[Grid_getBlkPtr] DEPRECATED: Use tile object")
end subroutine Grid_getBlkPtr

