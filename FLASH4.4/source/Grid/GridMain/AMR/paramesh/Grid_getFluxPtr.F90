!!****f* source/Grid/GridMain/AMR/Amrex/Grid_getFluxPtr
!!
!! NAME
!!  Grid_getFluxPtr
!!
!! SYNOPSIS
!!
!!  Grid_getFluxPtr(block_metadata_t(IN) :: blockDesc,
!!                  real(pointer)        :: fluxPtrX(:,:,:,:),
!!                  real(pointer)        :: fluxPtrY(:,:,:,:),
!!                  real(pointer)        :: fluxPtrZ(:,:,:,:))
!!  
!! DESCRIPTION 
!!  
!!  Obtain pointers to the X, Y, and Z flux data for a single block.
!!  The data structures will only pertain to flux in the block interior and
!!  on the block faces.
!!
!!  Note that the data stored in these structures is flux and *not* flux 
!!  densities.
!! 
!! ARGUMENTS 
!!
!!  blockDesc - the block metatdata object associated with the block
!!              whose flux data is to be obtained
!!  fluxPtrX - Pointer to the the Flux data along the IAXIS
!!  fluxPtrY - Pointer to the the Flux data along the JAXIS
!!  fluxPtrZ - Pointer to the the Flux data along the KAXIS
!!
!! NOTES
!!
!!  Once finished with the pointer, call Grid_releaseFluxPtr to release the
!!  pointer.
!!
!! SEE ALSO
!!  Grid_releaseFluxPtr
!!
!!***

#include "constants.h"

subroutine Grid_getFluxPtr(blockDesc, fluxPtrX, fluxPtrY, fluxPtrZ)
  use block_metadata,   ONLY : block_metadata_t
  use gr_specificData, ONLY : gr_flxx, gr_flxy, gr_flxz
  implicit none
  

  type(block_metadata_t), intent(IN) :: blockDesc
  real, pointer                      :: fluxPtrX(:,:,:,:)
  real, pointer                      :: fluxPtrY(:,:,:,:)
  real, pointer                      :: fluxPtrZ(:,:,:,:)

  integer :: blockID
  blockID=blockDesc%id

  fluxPtrX => gr_flxx(:,:,:,:,blockID)
  fluxPtrY => gr_flxy(:,:,:,:,blockID)
  fluxPtrZ => gr_flxz(:,:,:,:,blockID) 
end subroutine Grid_getFluxPtr

