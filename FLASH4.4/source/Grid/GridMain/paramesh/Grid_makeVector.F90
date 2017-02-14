!!****if* source/Grid/GridMain/Grid_makeVector
!!
!! NAME
!!  Grid_makeVector
!!
!! SYNOPSIS
!!
!!  Grid_makeVector(integer(IN)               :: vecLen,
!!                  real,dimension(:,:) (OUT) :: newVec,
!!                  integer  (INOUT)          :: numVec,
!!                  integer(IN)               :: gridDataStruct)
!!  
!! DESCRIPTION 
!!
!!  This routine converts solution data organized as blocks into a collection of vectors
!!  The length of the vector is an input "vecLen", which can be smaller or bigger in size than 
!!  the data contained in one block. The value will typically be dictated by the constraints of
!!  the target hardware. The value "numVec" is the number of vectors that will be generated from 
!!  flattening of N blocks. Its value is N * oneBlockSize / vecLen. The newly generated vectors are
!!  stored in newVec, which is at the moment a 2D array. Part of the exercise is to determine if it 
!!  should be a 2D or 3D array and what should be the data layout for different variables. Should
!!  variable be the leading dimension or should space be the leading dimension.
!!
!! ARGUMENTS 
!!
!!   vecLen     :  the length of vector into which solution data needs to be converted
!!   newVect    : storage for newly generated vectors
!!   numVec     : number of vectors to be generated from all that data contained in all blocks
!!   gridDataStruct : whether cell centered, face centered etc. (may be deprecated later)
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

subroutine Grid_makeVector(vecLen,newVec,numVec,gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_getData
  use Grid_data, ONLY :  gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  implicit none

  integer, intent(in) :: vecLen
  integer,intent(INOUT) :: numVec
  real, dimension(vecLen,numVec),intent(OUT) :: newVec
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

  oneBlkSize = NXB*NYB*NZB
  numVec = (vecLen+oneBlkSize-1)/oneBlkSize
  nblks=(gr_blkCount+numVec-1)/numVec
  blkID=1
  do j=1,numVec
     ptr=1
     do i=1,nblks
        if(blkID.le.gr_blkCount) then
           dataPtr => unk(:,:,:,:,blockid)
           call Eos_getData(range,vecLen,ptr,dataPtr,gridDataStruct,newVec(ptr,j))
           nullify(dataPtr)
!!           call Eos_getData(range,oneBlkSize,ptr,solnData,gridDataStruct,newVec(ptr,j))
           ptr=ptr+oneBlkSize
           blkID=blkID+1
        end if
     end do
  end do
  return
end subroutine Grid_makeVector








