!!****ih* source/physics/Hydro/localAPI/hy_memInterface
!!
!! NAME
!!  hy_memInterface
!!
!! SYNOPSIS
!!  use hy_memInterface
!!
!! DESCRIPTION
!!  This is a local interface for memory management routines.
!!
!!***


Module hy_memInterface

  use hy_memScratchData, ONLY: hy_memScratchInitialize, hy_memAllocScratch, hy_memDeallocScratch

  implicit none

#include "FortranLangFeatures.fh"

  interface hy_memGetBlkPtr
     subroutine hy_memGetBlkPtr(blockID,dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlkPtr
     subroutine hy_memGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlk5Ptr
     subroutine hy_memGetBlkPtr_desc(tileDesc,dataPtr, gridDataStruct)
       use flash_tile, ONLY : flash_tile_t 
       implicit none
       type(flash_tile_t), intent(IN) :: tileDesc
       real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlkPtr_desc
  end interface

  interface hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlkPtr(blockId, dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlk5Ptr
     subroutine hy_memReleaseBlkPtr_desc(tileDesc, dataPtr, gridDataStruct)
       use flash_tile, ONLY : flash_tile_t 
       implicit none
       type(flash_tile_t), intent(IN) :: tileDesc
       real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlkPtr_desc
     subroutine hy_memReleaseBlk5Ptr_desc(blockDesc, data5Ptr, gridDataStruct)
       use block_metadata,   ONLY : block_metadata_t
       implicit none
       type(block_metadata_t), intent(IN) :: blockDesc
       real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlk5Ptr_desc
  end interface

End Module hy_memInterface
