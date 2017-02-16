!!****if* source/Grid/GridMain/Grid_makeVector
!!
!! NAME
!!  Grid_makeVector
!!
!! SYNOPSIS
!!
!!  call Grid_makeVector(integer(IN)        :: vecLen,
!!                       integer(IN)        :: numVars,
!!                  real(OUT),dimension(vecLen,numVars,numVec) :: newVec,
!!                  integer(INOUT)          :: numVec,
!!                  OPTIONAL,integer(OUT)   :: vecLastFree,
!!                  OPTIONAL,integer(IN)    :: gridDataStruct)
!!  
!! DESCRIPTION 
!!
!!  This routine converts solution data organized as blocks into a collection of vectors.
!!
!!  The length of the vector is an input "vecLen", which can be smaller or bigger in size than 
!!  the data contained in one block. The value will typically be dictated by the constraints of
!!  the target hardware. The value "numVec" is the number of vectors that will be generated from 
!!  flattening of N blocks. Its value is N * oneBlockSize / vecLen. The newly generated vectors are
!!  stored in newVec, which is in this version a 3D array. Part of the exercise is to determine if it 
!!  should be a 2D or 3D array and what should be the data layout for different variables. Should
!!  variable be the leading dimension or should space be the leading dimension.
!!
!! ARGUMENTS 
!!
!!   vecLen     :  the length of vector into which solution data needs to be converted
!!   numVars    :  number of variables per cell that are stored in newVect.
!!                 For use with the Eos unit as in this example, this should be = EOS_NUM.
!!   newVect    : storage for newly generated vectors.
!!                The repackaged (and Eos_getData-preprocessed) data is returned here,
!!                ready to call Eos on.
!!   numVec     : number of vectors (of numVars variables each) to be generated
!!                from all that data contained in all blocks.
!!                On return, numVec may be modified (lowered) so it counts only
!!                vectors actually used.
!!   vecLastFree : If this optional argument is present, it will on return contain a
!!                 count of free (unused) positions of the last vector.
!!                 That is, the valid vectors returned consist of
!!                   newVec(1:vecLen            ,:,ivec) , ivec=1,numVec-1 ,
!!                   newVec(1:vecLen-vecLastFree,:,numVec);
!!                 where numVec is the count of vectors returned (possibly lowered).
!!                 If the argument is not present, no indication of incomplete vectors
!!                 is provided to the caller.
!!                 If returning all blockdata as vectors would overflow the available space
!!                 in newVect, then vecLastFree (if present) will be set to -1 to indicate
!!                 this condition.
!!   gridDataStruct : whether cell centered, face centered etc. (may be deprecated later).
!!                    Should be CENTER if present. FLASH does not require
!!                    support for marshalling other datastructs (FACEX, etc.)
!!                    for Eos calls.
!!
!! NOTES
!!
!!  This won;t work it is just to show.
!!***

!!REORDER(5): unk
!!REORDER(4): dataPtr

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_makeVector(vecLen,numVars,newVec,numVec,vecLastFree,gridDataStruct)

#include "constants.h"
#include "Flash.h"

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
  use tree,         ONLY : gr_blkCount => lnblocks ! there is no FLASH-owned reliable gr_blkCount...
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY : Eos_getData
  use Grid_data, ONLY :  gr_meshMe
  use Grid_data, ONLY :  gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi

  implicit none

  integer, intent(in) :: vecLen
  integer, intent(in) :: numVars
  integer,intent(INOUT) :: numVec
  real, dimension(vecLen,numVars,numVec),intent(OUT) :: newVec
  integer, optional,intent(OUT):: vecLastFree
  integer, optional,intent(in) :: gridDataStruct

  integer :: oneBlkSize, i,ptr,blkID
  integer :: numVecIn
  integer,dimension(LOW:HIGH,MDIM) :: range
  real,dimension(:,:,:,:),pointer :: dataPtr

  integer :: vecFree, ivec, subBlkSize
  integer,dimension(LOW:HIGH,MDIM) :: r0, sr0, eosRange
  integer,dimension(         MDIM) :: blkSize,sublen,stride,npos
  logical :: incomplete, doit
  integer :: axis, carryAxis
  integer,parameter :: hiAxis = IAXIS+NDIM-1

#define ASSERT(condition) if (condition) call Driver_abortFlash("Grid_makeVector: ASSERT failed!")
  range(LOW,IAXIS)=gr_ilo
  range(HIGH,IAXIS)=gr_ihi
  range(LOW,JAXIS)=gr_jlo
  range(HIGH,JAXIS)=gr_jhi
  range(LOW,KAXIS)=gr_klo
  range(HIGH,KAXIS)=gr_khi
  blkSize (IAXIS:KAXIS) = range(HIGH,:)-range(LOW,:)+1
  r0(LOW,:)  = 0                !zero-based range
  r0(HIGH,:)  = blkSize(:)      !Python-like convention for end of range
  sr0     (:,:)         = r0   (:,:) !zero-based subrange
  numVecIn = numVec

  do i=IAXIS,KAXIS
     stride(i) = PRODUCT(blkSize(IAXIS:i-1))
  end do

  ! Number of cells avalable:
  !    M  =  gr_blkCount    *     oneBlksize

  oneBlkSize = NXB*NYB*NZB      ! number of (interior) cells in each
                                ! solution data block
!!$  blksPerVec = (vecLen+oneBlkSize-1)/oneBlkSize ! number of solution data blocks from
!!$                                            ! FLASH needed to provide one vector's length of data
!!$  vecsPerBlk = (vecLen+oneBlkSize-1)/vecLen

  numVec = (vecLen+gr_blkCount*oneBlkSize-1)/vecLen ! number of vectors needed.
                                                    ! That is N * oneBlockSize / vecLen (rounded up)
  ! Thus numVec*oneBlkSize :                 ! number of cells per vector,
                                             ! rounded up to a multiple of oneBlkSize
  ASSERT (numVecIn .GE. numVec)              ! ouput array must have enough space!
  ASSERT (vecLen == numVec * oneBlkSize)     ! maxSize must be a multiple of oneBlkSize !


  ptr     = 1
  vecFree = vecLen

  npos(:) = 0
  blkID=1
  if (blkID .LE. gr_blkCount) then
#ifdef HAVE_PTRBNDREMAP
     dataPtr(1:,0:,0:,0:) => unk(PROP_VARS_BEGIN:,range(LOW,IAXIS):,range(LOW,JAXIS):,range(LOW,KAXIS):,blkID)
##else
     dataPtr => unk(PROP_VARS_BEGIN:,:,:,:,blkID)
#endif
  else
     nullify(dataPtr)
  end if

  iv:do ivec=1,numVec
     ptr     = 1
     vecFree = vecLen
     vf:do while (vecFree > 0)

        pre:do axis=IAXIS,hiAxis-1
           incomplete = ( npos(axis) > 0 .OR. vecFree < blkSize(axis)*stride(axis) )
           if (incomplete) then
              sr0(:,   :   ) = r0(:,:)
              sr0(LOW, axis) = npos(axis)
              sr0(HIGH,axis) = min(npos(axis)+vecFree/stride(axis),blkSize(axis))
              sublen(  :   ) = blkSize(:)
              sublen(  axis) = sr0(HIGH,axis) - sr0(LOW,axis)
              if (sublen(axis) == 0) EXIT pre
              sublen(axis+1:hiAxis) = 1

#ifdef HAVE_PTRBNDREMAP
              eosRange(LOW, :) = sr0(LOW, :)
              eosRange(HIGH,:) = sr0(HIGH,:)-1
#else
              eosRange(LOW, :) = range(LOW, :) + sr0(LOW, :)
              eosRange(HIGH,:) = range(LOW, :) + sr0(HIGH,:)-1
#endif
              subBlkSize = PRODUCT(sublen)
              call Eos_getData(eosRange,vecLen,dataPtr,gridDataStruct,newVec(ptr,1,ivec))
              ptr=ptr+subBlkSize
              vecFree = vecFree - subBlkSize
              npos(axis) = npos(axis) + sublen(axis)
              carryAxis = axis
              do carryAxis=axis,hiAxis-1
                 if (npos(carryAxis)==blkSize(carryAxis)) then !carry
                    npos(carryAxis+1) = npos(carryAxis+1) + 1
                    npos(carryAxis) = 0
                 end if
              end do
              if (npos(hiAxis)==blkSize(hiAxis)) then !this block is done
                 blkID=blkID+1
                 if (blkID .LE. gr_blkCount) then
#ifdef HAVE_PTRBNDREMAP
                    dataPtr(1:,0:,0:,0:) => &
                         unk(PROP_VARS_BEGIN:,range(LOW,IAXIS):,range(LOW,JAXIS):,range(LOW,KAXIS):,blkID)
##else
                    dataPtr => unk(PROP_VARS_BEGIN:,:,:,:,blkID)
#endif
                 else
                    nullify(dataPtr)
                    EXIT iv
                 end if
                 npos(hiAxis) = 0
              end if

              if (vecFree == 0) EXIT vf
           end if

        end do pre


        hiAx:do while (vecFree>0)
           incomplete = ( npos(hiAxis) > 0 .OR. vecFree < oneBlkSize )
           if (incomplete) then
              sr0(:,   :   ) = r0(:,:)
              sr0(LOW, hiAxis) = npos(hiAxis)
              sr0(HIGH,hiAxis) = min(npos(hiAxis)+vecFree/stride(axis),blkSize(hiAxis))
              sublen(  :   ) = blkSize(:)
              if (sublen(axis) == 0) EXIT hiAx
              sublen(  hiAxis) = sr0(HIGH,hiAxis) - sr0(LOW,hiAxis)

#ifdef HAVE_PTRBNDREMAP
              eosRange(LOW, :) = sr0(LOW, :)
              eosRange(HIGH,:) = sr0(HIGH,:)-1
#else
              eosRange(:  , :) = range(:, :)
#endif
              subBlkSize = PRODUCT(sublen)
              call Eos_getData(eosRange,vecLen,dataPtr,gridDataStruct,newVec(ptr,1,ivec))
              ptr=ptr+subBlkSize
              vecFree = vecFree - subBlkSize
           else
              sr0(:,   :   ) = r0(:,:)
              sublen(  :   ) = blkSize(:)

#ifdef HAVE_PTRBNDREMAP
              eosRange(LOW, :) = sr0(LOW, :)
              eosRange(HIGH,:) = sr0(HIGH,:)-1
#else
              eosRange(LOW, :) = range(LOW, :) + sr0(LOW, :)
              eosRange(HIGH,:) = range(LOW, :) + sr0(HIGH,:)-1
#endif
              subBlkSize = oneBlkSize
              call Eos_getData(eosRange,vecLen,dataPtr,gridDataStruct,newVec(ptr,1,ivec))
              ptr=ptr+subBlkSize
              vecFree = vecFree - subBlkSize
           end if
           npos(hiAxis) = npos(hiAxis) + sublen(hiAxis)
           if (npos(hiAxis)==blkSize(hiAxis)) then !this block is done
              blkID=blkID+1
              if (blkID .LE. gr_blkCount) then
#ifdef HAVE_PTRBNDREMAP
                 dataPtr(1:,0:,0:,0:) => &
                      unk(PROP_VARS_BEGIN:,range(LOW,IAXIS):,range(LOW,JAXIS):,range(LOW,KAXIS):,blkID)
##else
                 dataPtr => unk(PROP_VARS_BEGIN:,:,:,:,blkID)
#endif
              else
                 nullify(dataPtr)
                 EXIT iv
              end if
              npos(hiAxis) = 0
           end if

           if (vecFree == 0) EXIT vf
        end do hiAx

        post:do axis=hiAxis-1,IAXIS,-1
           incomplete = ( vecFree < blkSize(axis)*stride(axis) )
           ASSERT( incomplete )
           doit = ( vecFree .GE. stride(axis) ) 
           if (doit) then
              sr0(:,   :   ) = r0(:,:)
              sr0(LOW, axis) = npos(axis)
              sr0(HIGH,axis) = min(npos(axis)+vecFree/stride(axis),blkSize(axis))
              sublen(  :   ) = blkSize(:)
              sublen(  axis) = sr0(HIGH,axis) - sr0(LOW,axis)
              sublen(axis+1:hiAxis) = 1

#ifdef HAVE_PTRBNDREMAP
              eosRange(LOW, :) = sr0(LOW, :)
              eosRange(HIGH,:) = sr0(HIGH,:)-1
#else
              eosRange(LOW, :) = range(LOW, :) + sr0(LOW, :)
              eosRange(HIGH,:) = range(LOW, :) + sr0(HIGH,:)-1
#endif
              subBlkSize = PRODUCT(sublen)
              call Eos_getData(eosRange,vecLen,dataPtr,gridDataStruct,newVec(ptr,1,ivec))
              ptr=ptr+subBlkSize
              vecFree = vecFree - subBlkSize
           end if

        end do post
     end do vf
  end do iv

  if (present(vecLastFree)) then
     if (associated(dataPtr)) then !vecFree SHOULD always be 0 in this case.
        vecLastFree = -1
     else
        vecLastFree = vecFree
     end if
  else if (associated(dataPtr)) then !vecFree SHOULD always be 0 in this case.
99   format('PE',I6,': newVec(',I4,',:,',I4,') not large enough for',I6,' blocks, vecFree=',I5)
     print 99,gr_meshMe,vecLen,numVecIn,gr_blkCount,vecFree
     call Driver_abortFlash("Grid_makeVector: newVect not large enough for block data")
  end if

  return
end subroutine Grid_makeVector








