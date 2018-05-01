!!****f* source/Grid/GridMain/AMR/Amrex/Grid_getFluxPtr
!!
!! NAME
!!  Grid_getFluxPtr
!!
!! SYNOPSIS
!!
!!  call Grid_getFluxPtr(block_metadata_t(IN) :: blockDesc,
!!                       real(pointer)        :: fluxPtrX(:,:,:,:),
!!                       real(pointer)        :: fluxPtrY(:,:,:,:),
!!                       real(pointer)        :: fluxPtrZ(:,:,:,:))
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
!!  The indexing depends on the INDEXREORDER preprocessor symbol, which needs
!!  to be passed by the compiler on the command line if ncessary.
!!
!! SEE ALSO
!!  Grid_releaseFluxPtr
!!
!!***

#include "constants.h"

subroutine Grid_getFluxPtr(blockDesc, fluxPtrX, fluxPtrY, fluxPtrZ)
  use block_metadata,   ONLY : block_metadata_t
  use gr_specificData, ONLY : gr_flxx, gr_flxy, gr_flxz
  use gr_specificData, ONLY : gr_loFl
  implicit none
  

  type(block_metadata_t), intent(IN) :: blockDesc
  real, pointer                      :: fluxPtrX(:,:,:,:)
  real, pointer                      :: fluxPtrY(:,:,:,:)
  real, pointer                      :: fluxPtrZ(:,:,:,:)

  integer :: blockID
  blockID=blockDesc%id

  associate (lo => gr_loFl)
#ifdef INDEXREORDER
    fluxPtrX(lo(1):, lo(2):, lo(3):, 1:) => gr_flxx(:,:,:,:,blockID)
    fluxPtrY(lo(1):, lo(2):, lo(3):, 1:) => gr_flxy(:,:,:,:,blockID)
    fluxPtrZ(lo(1):, lo(2):, lo(3):, 1:) => gr_flxz(:,:,:,:,blockID)
#else
    fluxPtrX(1:, lo(1):, lo(2):, lo(3):) => gr_flxx(:,:,:,:,blockID)
    fluxPtrY(1:, lo(1):, lo(2):, lo(3):) => gr_flxy(:,:,:,:,blockID)
    fluxPtrZ(1:, lo(1):, lo(2):, lo(3):) => gr_flxz(:,:,:,:,blockID) 
#endif
  end associate

end subroutine Grid_getFluxPtr

