!!****if* source/Grid/GridMain/Grid_makeVector
!!
!! NAME
!!  Grid_makeVector
!!
!! SYNOPSIS
!!
!!  Grid_makeVector(integer(IN)             :: maxSize
!!                  real,dimension(maxSize,numVec) (OUT) :: newVec,
!!                  integer(INOUT)          :: numVec,
!!                  OPTIONAL,integer(IN)    :: gridDataStruct)
!!  
!! DESCRIPTION 
!!  
!!
!! ARGUMENTS 
!!
!!  maxSize - 
!!            Used to declare the shape of newVec.
!!  numVec  - number of 
!!  newVec  - the repackaged data is returned here, ready to call Eos on.
!!  gridDataStruct - Should be CENTER if present. FLASH does not require
!!                   support for marshalling other datastructs (FACEX, etc.)
!!                   for Eos calls.
!!
!! NOTES
!!
!!  This won;t work it is just to show.
!!***

!!REORDER(5): unk, facevar[xyz], scratch_ctr, scratch_facevar[xyz]
!!REORDER(4): dataPtr
!!FOR FUTURE: Add REORDER for unk, facevar[xyz]1, etc.?

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_makeVector(maxSize,newVec,numVec,gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_getData
  use Grid_data, ONLY :  gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  implicit none

  integer, intent(in) :: maxSize
  integer,intent(INOUT) :: numVec
  real, dimension(maxSize,numVec) :: newVec
  integer, optional,intent(in) :: gridDataStruct

  integer :: OneBlkSize, nblks,i,j,ptr,blkID
  integer,dimension(LOW:HIGH,MDIM) :: range
  real,dimension(:,:,:,:),pointer :: dataPtr

  range(LOW,IAXIS)=gr_ilo
  range(HIGH,IAXIS)=gr_ihi
  range(LOW,JAXIS)=gr_jlo
  range(HIGH,JAXIS)=gr_jhi
  range(LOW,KAXIS)=gr_klo
  range(HIGH,KAXIS)=gr_khi

  ! Number of cells avalable:
  !    M  =  gr_blkCount    *     oneBlksize
  ! Number of cells processed in the loop below:
  !    L  =  numVec(IN) * nblks * oneBlksize
  !    L >=  gr_blkCount        * oneBlksize   since  numVec * nblks >= gr_blkCount
  !    L >=               nblks * maxSize      since  numVec * oneBlksize >= maxSize

  oneBlkSize = NXB*NYB*NZB      ! number of (interior) cells in each
                                ! solution data block
  numVEc = (maxSize+oneBlkSize-1)/oneBlkSize ! number of solution data blocks from
                                             ! FLASH needed to provide one "vector's" worth of data 
  nblks=(gr_blkCount+numVec-1)/numVec        ! total number of FLASH solution data blocks,
                                             ! rounded up to ???
  blkID=1
  do j=1,numVec
     ptr=1
     do i=1,nblks
        if(blkID.le.gr_blkCount) then
           dataPtr => unk(:,:,:,:,blockid)
           call Eos_getData(range,maxLen,ptr,solnData,gridDataStruct,newVec(ptr,j))
           ptr=ptr+oneBlkSize
           blkID=blkID+1
        end if
     end do
  end do
  return
end subroutine Grid_makeVector








