!!****if* source/Grid/GridMain/Chombo/gr_releaseInteriorBlkPtr
!!
!! NAME
!!  gr_releaseInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_releaseInteriorBlkPtr(integer(in) :: blockID
!!                           real(pointer) :: dataPtr(:,:,:,:)
!!                           integer(in) :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
!!  
!! ARGUMENTS 
!!
!!  blockID - The block being pointed to.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - The type of block being pointed to.
!!
!!
!! SEE ALSO
!!  gr_getInteriorBlkPtr
!!***

#include "constants.h"
#include "Flash.h"

! DEVNOTE: Need REORDER directive here?
subroutine gr_releaseInteriorBlkPtr(block,dataPtr,gridDataStruct)
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                               Grid_getNumVars
  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t),intent(in) :: block
  real, pointer :: dataPtr(:,:,:,:)
  integer, intent(in) :: gridDataStruct

  call Driver_abortFlash("[gr_releaseInteriorBlkPtr]: Not implemented yet!")
end subroutine gr_releaseInteriorBlkPtr
