!!****if* source/physics/Hydro/HydroMain/hy_memGetBlkPtr
!!
!! NAME
!!  hy_memGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memGetBlkPtr(block_metadata(IN)  :: blockDesc,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
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

!!REORDER(5):hy_memArrayScratch, hy_memArrayScratch_ctr, hy_memArrayScratch_facevar[xyz]

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine hy_memGetBlkPtr_desc(blockDesc,dataPtr, gridDataStruct)

#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"

  use hy_memInterface, ONLY : hy_memGetBlkPtr
  use block_metadata,   ONLY : block_metadata_t
  implicit none

  type(block_metadata_t), intent(IN) :: blockDesc
  real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
  integer, optional,intent(in) :: gridDataStruct

  real, dimension(:,:,:,:), pointer :: medPtr
  integer,dimension(MDIM+1) :: lo
  integer :: blockID
#ifdef INDEXREORDER
  integer, parameter :: iX = 1
#else
  integer, parameter :: iX = 2
#endif


  blockID = blockDesc%id

  call hy_memGetBlkPtr(blockID,medPtr,gridDataStruct)

  lo = lbound(medPtr)
  lo(iX:ix+MDIM-1) = lo(iX:ix+MDIM-1) + blockDesc%limitsGC(LOW,:) &
                                       -blockDesc%localLimitsGC(LOW,:)

  dataPtr(lo(1):,lo(2):,lo(3):,lo(4):) => medPtr

  return
end subroutine hy_memGetBlkPtr_desc


subroutine hy_memGetBlk5Ptr_desc(blockDesc,data5Ptr, gridDataStruct)

  use hy_memInterface, ONLY : hy_memGetBlkPtr
  use block_metadata,   ONLY : block_metadata_t
  implicit none

  type(block_metadata_t), intent(IN) :: blockDesc

  real, dimension(:,:,:,:,:), pointer :: data5Ptr
  integer, optional,intent(in) :: gridDataStruct

  integer :: blockID

  blockID = blockDesc%id

  call hy_memGetBlkPtr(blockID,data5Ptr,gridDataStruct)

  !! DEV: Should probably do some index-shifting as in hy_memGetBlkPtr_desc above.

end subroutine hy_memGetBlk5Ptr_desc
