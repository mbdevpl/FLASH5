#include "Flash.h"

subroutine Grid_zeroFluxData()
  use amrex_amrcore_module, ONLY : amrex_get_finest_level

  use Grid_interface,       ONLY : Grid_getFluxPtr, Grid_releaseFluxPtr
  use gr_interface,         ONLY : gr_getBlkIterator, &
                                   gr_releaseBlkIterator
  use gr_iterator,          ONLY : gr_iterator_t
  use block_metadata,       ONLY : block_metadata_t

  implicit none

  type(gr_iterator_t)    :: itor
  type(block_metadata_t) :: block

  real, pointer :: fluxDataX(:, :, : , :) => null()
  real, pointer :: fluxDataY(:, :, : , :) => null()
  real, pointer :: fluxDataZ(:, :, : , :) => null()

  integer :: level
  integer :: finest_level

  if(NFLUXES < 1)   RETURN

  finest_level = amrex_get_finest_level()
  do level=1, finest_level 
    call gr_getBlkIterator(itor, level=level, tiling=.FALSE.)
    do while (itor%is_valid())
      call itor%blkMetaData(block)
         
      call Grid_getFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)
      if (associated(fluxDataX)) then
          fluxDataX(:,:,:,:) = 0.0
      end if
      if (associated(fluxDataY)) then
          fluxDataY(:,:,:,:) = 0.0
      end if
      if (associated(fluxDataZ)) then
          fluxDataZ(:,:,:,:) = 0.0
      end if
      call Grid_releaseFluxPtr(block, fluxDataX, fluxDataY, fluxDataZ)

      call itor%next()
    end do
    call gr_releaseBlkIterator(itor)
  end do
end subroutine Grid_zeroFluxData

