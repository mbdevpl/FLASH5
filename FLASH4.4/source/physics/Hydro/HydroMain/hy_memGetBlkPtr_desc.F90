!!****if* source/physics/Hydro/HydroMain/hy_memGetBlkPtr
!!
!! NAME
!!  hy_memGetBlkPtr
!!
!! SYNOPSIS
!!
!!  call hy_memGetBlkPtr(integer(IN)            :: blockID,
!!                 real(pointer)(:,:,:,:) :: dataPtr,
!!                 integer(IN),optional   :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Gets a pointer to a single block of simulation data from the
!!  specified Grid data structure. The block includes guard cells.
!!  If the optional argument "gridDataStructure" is not specified,
!!  it returns a block from cell centered data structure.
!!
!!  When using Paramesh 4 in NO_PERMANENT_GUARDCELLS mode, it is important to
!!  release the block pointer for a block before getting it for another block.
!!  For example if pointer to block 1 is not yet released and the user
!!  tries to get a pointer to block 2, the routine will abort.
!!
!! ARGUMENTS
!!
!!  blockID : the local blockid
!!
!!  dataPtr : Pointer to the data block
!!
!!  gridDataStruct : optional integer value specifying data structure.
!!                   The options are defined in constants.h and they are :
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!
!!
!!
!! NOTES
!!
!!  hy_memGetBlkPtr is an accessor function that passes a pointer
!!  as an argument and requires an explicit interface for most compilers.
!!
!!  Don't forget to call hy_memReleaseBlkPtr when you are finished with it!
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


contains
  subroutine AssoMed(pp, mm, leafNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:,:)
    integer,intent(in) :: leafNo
    call AssoFin(pp,mm(:,:,:,:,leafNo),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb1:,lb2:,lb3:,lb4:)
    pp => dd
  end subroutine AssoFin

end subroutine hy_memGetBlkPtr_desc


subroutine hy_memGetBlk5Ptr_desc(blockDesc,data5Ptr, gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use hy_memInterface, ONLY : hy_memGetBlkPtr
  use block_metadata,   ONLY : block_metadata_t
  implicit none

  type(block_metadata_t), intent(IN) :: blockDesc

  real, dimension(:,:,:,:,:), pointer :: data5Ptr
  integer, optional,intent(in) :: gridDataStruct

  integer :: blockID

  blockID = blockDesc%id

  call hy_memGetBlkPtr(blockID,data5Ptr,gridDataStruct)

contains
  subroutine AssoMed(pp, mm, leafNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:,:,:)
    integer,intent(in) :: leafNo
    call AssoFin(pp,mm(:,:,:,:,:,leafNo),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb1:,lb2:,lb3:,lb4:,:)
    pp => dd
  end subroutine AssoFin

end subroutine hy_memGetBlk5Ptr_desc
