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
!!  fluxPtrX - Pointer to the the Flux data along the IAXIS.
!!  fluxPtrY - Pointer to the the Flux data along the JAXIS
!!  fluxPtrZ - Pointer to the the Flux data along the KAXIS
!!
!! NOTES
!!
!!  The given pointers must not be associated and if a pointer is not needed
!!  due to the dimensionality of the problem, it will be set to NULL.
!!
!!  Once finished with the pointer, call Grid_releaseFluxPtr to release the
!!  pointer.
!!
!! SEE ALSO
!!  Grid_releaseFluxPtr
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getFluxPtr(blockDesc, fluxPtrX, fluxPtrY, fluxPtrZ)
    use Driver_interface,     ONLY : Driver_abortFlash
    use block_metadata,       ONLY : block_metadata_t
    use gr_physicalMultifabs, ONLY : fluxes

    implicit none

    type(block_metadata_t), intent(IN) :: blockDesc
    real, pointer                      :: fluxPtrX(:,:,:,:)
    real, pointer                      :: fluxPtrY(:,:,:,:)
    real, pointer                      :: fluxPtrZ(:,:,:,:)

    ! Avoid possible memory leaks
    if (     associated(fluxPtrX) &
        .OR. associated(fluxPtrY) & 
        .OR. associated(fluxPtrZ)) then
        call Driver_abortFlash("[Grid_getFluxPtr] All pointers must be NULL")
    end if
    nullify(fluxPtrX)
    nullify(fluxPtrY)
    nullify(fluxPtrZ)

    ! Note that limits(LOW, :) is the index of the lower-leftmost cell in the
    ! given block.  We are assigning this same index to the face just to the
    ! left of it (for X) or below it (for Y).
    !
    ! For an 8x8 cell block with the lower-leftmost cell indexed as (1, 1), 
    ! the lower X faces would be indexed from (1,1) to (9,1).
    associate (lo   => blockDesc%limits(LOW, :), &
               ilev => blockDesc%level - 1, &
               igrd => blockDesc%grid_index)
        fluxPtrX(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, IAXIS)%dataptr(igrd)
#if NDIM > 1
        fluxPtrY(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, JAXIS)%dataptr(igrd)
#if NDIM == 3
        fluxPtrZ(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, KAXIS)%dataptr(igrd)
#endif
#endif
    end associate
end subroutine Grid_getFluxPtr

