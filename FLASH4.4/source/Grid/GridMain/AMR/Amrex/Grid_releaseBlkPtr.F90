!!****if* source/Grid/GridMain/AMR/Amrex/Grid_releaseBlkPtr
!!
!! NAME
!!  Grid_releaseBlkPtr
!!
!! SYNOPSIS
!!
!!  Grid_releaseBlkPtr(integer(IN)   :: blockId,
!!                     real(pointer) :: dataPtr(:,:,:,:),
!!                     integer(IN),optional :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  Releases a pointer to a block.
!!  
!! ARGUMENTS 
!!
!!  blockId - ID of the block, should be the same ID was used in the
!!            corresponding Grid_getBlkPtr call.
!!  dataPtr - Pointer to be released.
!!  gridDataStruct - an optional argument that designates the type of grid data
!!                   structure to handle (i.e. facevar, unknown, scratch...)
!!
!!
!! NOTES
!!
!!  This implementation actually does more than just releasing the pointer.
!!
!! SEE ALSO
!!  Grid_getBlkPtr
!!***

subroutine Grid_releaseBlkPtr(blockId, dataPtr, gridDataStruct)
  implicit none

  integer, intent(in)              :: blockId
  real,    intent(inout), pointer  :: dataPtr(:,:,:,:)
  integer, intent(in),    optional :: gridDataStruct

  ! TODO: Remove this version once blockIDs are removed from interface.
  write(*,*) "AMReX does *not* deal in blockIDs"
  stop
end subroutine Grid_releaseBlkPtr

subroutine Grid_releaseBlkPtr_Itor(block, dataPtr, gridDataStruct)

#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling drift macro expansion, see: drift
#undef Grid_releaseBlkPtr
#endif

#include "constants.h"

  use block_metadata, ONLY : block_metadata_t
#if DRIFT_ENABLE
  use Driver_interface, only: Driver_driftBlock
  use Driver_data, only: dr_driftSrcFile, dr_driftSrcLine
#endif

  implicit none

  type(block_metadata_t), intent(in)              :: block
  real,                   intent(inout), pointer  :: dataPtr(:, :, :, :)
  integer,                intent(in),    optional :: gridDataStruct

  integer :: gds

  if(present(gridDataStruct)) then
     gds = gridDataStruct
  else
     gds = CENTER
  end if

#if DRIFT_ENABLE
  ! TODO: If this is to stay here, we need to convert Driver_driftBlock
  ! to take a block instead of blockID.
  if(dr_driftSrcLine >= 0) then
    call Driver_driftBlock(dr_driftSrcFile, dr_driftSrcLine, blockId, &
      dataPtr(:,&
        lbound(dataPtr,2)+K1D*NGUARD : ubound(dataPtr,2)-K1D*NGUARD,&
        lbound(dataPtr,3)+K2D*NGUARD : ubound(dataPtr,3)-K2D*NGUARD,&
        lbound(dataPtr,4)+K3D*NGUARD : ubound(dataPtr,4)-K3D*NGUARD),&
      gds)
  end if
  dr_driftSrcLine = -1
#endif

  ! always destroy the pointer, because the other users will have their 
  ! own
  nullify(dataPtr)

end subroutine Grid_releaseBlkPtr_Itor

