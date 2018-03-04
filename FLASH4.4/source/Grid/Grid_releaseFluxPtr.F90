!!****f* source/Grid/Grid_releaseFluxPtr
!!
!! NAME
!!  Grid_releaseFluxPtr
!!
!! SYNOPSIS
!!
!!  Grid_releaseFluxPtr(block_metadata_t(IN) :: blockDesc,
!!                      real(pointer)        :: fluxPtrX(:,:,:,:),
!!                      real(pointer)        :: fluxPtrY(:,:,:,:),
!!                      real(pointer)        :: fluxPtrZ(:,:,:,:))
!!  
!! DESCRIPTION 
!!  Releases pointers to single blocks of flux data.
!!  
!! ARGUMENTS 
!!  blockDesc - the block metatdata object associated with the block
!!              whose flux data is to be obtained
!!  fluxPtrX - Pointer to the the Flux data along the IAXIS
!!  fluxPtrY - Pointer to the the Flux data along the JAXIS
!!  fluxPtrZ - Pointer to the the Flux data along the KAXIS
!!
!! NOTES
!!
!!  Some implementations actually do more than just releasing the pointer.
!!
!! SEE ALSO
!!  Grid_getFluxPtr
!!
!!***

subroutine Grid_releaseFluxPtr(blockDesc, fluxPtrX, fluxPtrY, fluxPtrZ)
    use block_metadata, ONLY : block_metadata_t

    implicit none

    type(block_metadata_t), intent(IN)    :: blockDesc
    real, pointer,          intent(INOUT) :: fluxPtrX(:,:,:,:)
    real, pointer,          intent(INOUT) :: fluxPtrY(:,:,:,:)
    real, pointer,          intent(INOUT) :: fluxPtrZ(:,:,:,:)

    nullify(fluxPtrX)
    nullify(fluxPtrY)
    nullify(fluxPtrZ)
end subroutine Grid_releaseFluxPtr

