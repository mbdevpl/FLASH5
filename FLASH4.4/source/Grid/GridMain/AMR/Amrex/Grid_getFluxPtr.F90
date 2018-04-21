!!****if* source/Grid/GridMain/AMR/Amrex/Grid_getFluxPtr
!!
!! NAME
!!  Grid_getFluxPtr
!!
!! SYNOPSIS
!!
!!  Grid_getFluxPtr(block_metadata_t(IN)   :: block,
!!                  real(pointer)        :: fluxPtrX(:,:,:,:),
!!                  real(pointer)        :: fluxPtrY(:,:,:,:),
!!                  real(pointer)        :: fluxPtrZ(:,:,:,:))
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
!!  block : derived type containing metadata for block whose data we need to
!!          access
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
!!  Grid_getFluxPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call Grid_releaseFluxPtr when you are finished with it!
!!
!!  This implementations stops with an error message. That is because calling
!!  Grid_getFluxPtr with a blockID as first argument is a leftover from FLASH4
!!  Grid implementations that is not supported any more when the Grid implementation
!!  is purely based on AMReX!
!!***

subroutine Grid_getFluxPtr(blockID, fluxPtrX, fluxPtrY, fluxPtrZ)
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none
  
  integer, intent(IN)            :: blockID
  real, pointer                      :: fluxPtrX(:,:,:,:)
  real, pointer                      :: fluxPtrY(:,:,:,:)
  real, pointer                      :: fluxPtrZ(:,:,:,:)
  
  call Driver_abortFlash("[Grid_getFluxPtr] AMReX does *not* deal in block IDs")
end subroutine Grid_getFluxPtr

