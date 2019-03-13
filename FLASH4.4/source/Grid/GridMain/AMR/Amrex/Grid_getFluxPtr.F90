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

    real, pointer                      :: auxPtr(:,:,:,:)

    ! Avoid possible memory leaks
    if (     associated(fluxPtrX) &
        .OR. associated(fluxPtrY) & 
        .OR. associated(fluxPtrZ)) then
        call Driver_abortFlash("[Grid_getFluxPtr] All pointers must be NULL")
    end if
    nullify(fluxPtrX)
    nullify(fluxPtrY)
    nullify(fluxPtrZ)

    ! THE WHOLE BLOCK CASE:
    !
    ! Note that limits(LOW, :) is the index of the lower-leftmost cell in the
    ! given block.  We are assigning this same index to the face just to the
    ! left of it (for X) or below it (for Y).
    !
    ! For an 8x8 cell block with the lower-leftmost cell indexed as (1, 1), 
    ! the lower X faces would be indexed from (1,1) to (9,1).

    ! WHEN THE DESCRIPTOR DESCRIBES A TILE (smaller than a block):
    !
    ! Fluxes are currently stored without guard cells,
    ! in face-centered multifabs.
    ! We use the limits (not limitsGC) component of the
    ! descriptor to determine the region to point to.
    ! Note that in the case of tiling (when the descriptor
    ! describes not a full block but a tile), the region
    ! pointed to will in general NOT encompass all the
    ! interior cells of a block, but only the interior
    ! cells included in the tile (similar to the AMReX
    ! method tilebox()). ALSO NOTE that in the latter
    ! case, the region pointed to, and thus in Fortran
    ! 2008 terminology the pointer, will NOT be
    ! CONTIGUOUS.
    ! FURTHER NOTE that false sharing of some cell faces
    ! between neighboring tiles can occur, this is is
    ! problem that must be addressed before tiling can
    ! be used for computing fluxes.
    associate (lo   => blockDesc%limits(LOW, :), &
               ilev => blockDesc%level - 1, &
               igrd => blockDesc%grid_index)
!!$      print*,' ========== ilev=',ilev,', igrd=',igrd,' =========='
!!$      print*,'desc%       lo     :',lo
!!$      print*,'desc%limits(HIGH,:):',blockDesc%limits(HIGH,:)
        auxPtr                               => fluxes(ilev, IAXIS)%dataptr(igrd)
!!$      print*,'lbound(auxPtr):',lbound(auxPtr)
!!$      print*,'ubound(auxPtr):',ubound(auxPtr)
        fluxPtrX(lo(1):, lo(2):, lo(3):, 1:) => auxPtr(lo(1)-1:, lo(2)-K2D:, lo(3)-K3D:, :)
!!$      print*,'lbound(fluxPtrX):',lbound(fluxPtrX)
!!$      print*,'ubound(fluxPtrX):',ubound(fluxPtrX)
#if NDIM > 1
        auxPtr                               => fluxes(ilev, JAXIS)%dataptr(igrd)
        fluxPtrY(lo(1):, lo(2):, lo(3):, 1:) => auxPtr(lo(1)-1:, lo(2)-1:, lo(3)-K3D:, :)
#if NDIM == 3
        auxPtr                               => fluxes(ilev, KAXIS)%dataptr(igrd)
        fluxPtrZ(lo(1):, lo(2):, lo(3):, 1:) => auxPtr(lo(1)-1:, lo(2)-1:, lo(3)-1:, :)
#endif
#endif
    end associate
end subroutine Grid_getFluxPtr

