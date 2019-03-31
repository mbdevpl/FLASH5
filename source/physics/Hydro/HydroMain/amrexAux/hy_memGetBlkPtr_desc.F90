!!****if* source/physics/Hydro/HydroMain/hy_memGetBlkPtr
!!
!! NAME
!!  hy_memGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memGetBlkPtr(Grid_tile(IN)          :: blockDesc,
!!                       real(pointer)(:,:,:,:) :: dataPtr,
!!                       integer(IN),optional   :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Gets a pointer to a single block of allocated scratch data that conforms with
!!  the specified Grid data structure.
!!
!!  The block data may include zero, one, or several layers of guard cells,
!!  dependent on how hy_memAllocScratch was called for the gridDataStruct.
!!  If the optional argument "gridDataStruct" is not specified,
!!  it returns a block from cell centered data structure.
!!
!! ARGUMENTS
!!
!!  blockDesc : describes the local block
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure.
!!                   The options are defined in constants.h and they are :
!!                   CENTER cell-centered scratch space similar to UNK
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space for facex variables
!!                   SCRATCH_FACEY scratch space for facey variables
!!                   SCRATCH_FACEZ scratch space for facez variables
!!
!!
!!
!! NOTES
!!
!!  hy_memGetBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface.
!!
!!  If you call hy_memGetBlkPtr, you should also  call hy_memReleaseBlkPtr when you
!!  are done using the pointer.
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine hy_memGetBlkPtr_desc(tileDesc,dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  use Grid_tile, ONLY : Grid_tile_t
  implicit none

  type(Grid_tile_t), intent(IN) :: tileDesc
  real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
  integer, optional,intent(in) :: gridDataStruct

  call tileDesc%getDataPtr(dataPtr, gridDataStruct)
end subroutine hy_memGetBlkPtr_desc


! Note: there is currently no Amrex-based implementation of the
!       following interface:
!!$subroutine hy_memGetBlk5Ptr_desc(blockDesc,data5Ptr, gridDataStruct)
!!$  use Grid_tile,   ONLY : Grid_tile_t
!!$  implicit none
!!$  type(Grid_tile_t), intent(IN) :: blockDesc
!!$  real, dimension(:,:,:,:,:), pointer :: data5Ptr
!!$  integer, optional,intent(in) :: gridDataStruct
!!$end subroutine hy_memGetBlk5Ptr_desc
