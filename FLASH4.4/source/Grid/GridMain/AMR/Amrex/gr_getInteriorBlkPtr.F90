!!****if* source/Grid/GridMain/Chombo/gr_getInteriorBlkPtr
!!
!! NAME
!!  gr_getInteriorBlkPtr
!!
!! SYNOPSIS
!!
!!  gr_getInteriorBlkPtr(integer(IN)      :: blockID,
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

#include "constants.h"
#include "Flash.h"

! DEVNOTE: Need REORDER directive here?
subroutine gr_getInteriorBlkPtr(block, dataPtr, gridDataStruct)
  use Grid_interface,   ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                               Grid_getNumVars
  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata,   ONLY : block_metadata_t
 
  implicit none

  type(block_metadata_t), intent(in) :: block
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer, intent(in) :: gridDataStruct

  call Driver_abortFlash("[gr_getInteriorBlkPtr]: Not implemented yet!")
end subroutine gr_getInteriorBlkPtr

