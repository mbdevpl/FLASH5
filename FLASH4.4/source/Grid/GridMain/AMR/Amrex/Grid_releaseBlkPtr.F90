#include "constants.h"
#include "FortranLangFeatures.fh"

subroutine Grid_releaseBlkPtr(blockID, blkPtr, gridDataStruct)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(in)              :: blockID
  real,    intent(inout), pointer  :: blkPtr(:, :, :, :)
  integer, intent(in),    optional :: gridDataStruct

  call Driver_abortFlash("[Grid_releaseBlkPtr] AMReX does *not* deal in block IDs")
end subroutine Grid_releaseBlkPtr

subroutine Grid_releaseBlkPtr_Itor(block, blkPtr, gridDataStruct)

!#include "Flash.h"
!#ifdef Grid_releaseBlkPtr
!! disabling drift macro expansion, see: drift
!#undef Grid_releaseBlkPtr
!#endif

  use amrex_fort_module, ONLY : wp => amrex_real

  use block_metadata, ONLY : block_metadata_t
#if DRIFT_ENABLE
  use Driver_interface, only: Driver_driftBlock, &
                              Driver_abortFlash
  use Driver_data, only: dr_driftSrcFile, dr_driftSrcLine
#endif

  implicit none

  ! DEV: How to match data types for blkPtr with FLASH?
  ! DEV: FIXME Need to use POINTER_INTENT_INOUT
  type(block_metadata_t), intent(in)              :: block
  real(wp),               intent(inout), pointer  :: blkPtr(:, :, :, :)
  integer,                intent(in),    optional :: gridDataStruct

  integer :: gds
 
  if (present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

#ifdef FL_NON_PERMANENT_GUARDCELLS
  ! DEV: TODO Will we use NONPERMANENT GUARDCELLS with AMReX?
  call Driver_abortFlash("[Grid_releaseBlkPtr_desc] NON-PERMANENT GCs not implemented yet")
#else
  if (      (gds /= CENTER) &
      .AND. (gds /= FACEX) .AND. (gds /= FACEY) .AND. (gds /= FACEZ)  &
      .AND. (gds /= SCRATCH_CTR)                                    ) then
     call Driver_abortFlash("[Grid_releaseBlkPtr_desc] gridDataStruct not implemented yet")
  end if
#endif

#if DRIFT_ENABLE
  ! DEV: TODO If this is to stay here, we need to convert Driver_driftBlock
  ! to take a block instead of blockID.
  if(dr_driftSrcLine >= 0) then
    call Driver_driftBlock(dr_driftSrcFile, dr_driftSrcLine, blockId, &
      blkPtr(:,&
        lbound(blkPtr,2)+K1D*NGUARD : ubound(blkPtr,2)-K1D*NGUARD,&
        lbound(blkPtr,3)+K2D*NGUARD : ubound(blkPtr,3)-K2D*NGUARD,&
        lbound(blkPtr,4)+K3D*NGUARD : ubound(blkPtr,4)-K3D*NGUARD),&
      gds)
  end if
  dr_driftSrcLine = -1
#endif

  ! always destroy the pointer, because the other users will have their 
  ! own
  nullify(blkPtr)

end subroutine Grid_releaseBlkPtr_Itor

