!!****f* source/Grid/Grid_getFluxPtr_desc
!!
!! NAME
!!  Grid_getFluxPtr_desc
!!
!! SYNOPSIS
!!
!!  Grid_getFluxPtr_desc(block_metadata_t(IN) :: blockDesc,
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

subroutine Grid_getFluxPtr_desc(blockDesc, fluxPtrX, fluxPtrY, fluxPtrZ)
    use block_metadata,   ONLY : block_metadata_t

    implicit none

    type(block_metadata_t), intent(IN) :: blockDesc
    real, pointer                      :: fluxPtrX(:,:,:,:)
    real, pointer                      :: fluxPtrY(:,:,:,:)
    real, pointer                      :: fluxPtrZ(:,:,:,:)

end subroutine Grid_getFluxPtr_desc

