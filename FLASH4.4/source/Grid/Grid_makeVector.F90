!!****if* source/Grid/Grid_makeVector
!!
!! NAME
!!  Grid_makeVector
!!
!! SYNOPSIS
!!
!!  call Grid_makeVector(integer(IN)        :: vecLen,
!!                       integer(IN)        :: numVars,
!!                  real(INOUT),dimension(vecLen,numVars,numVec) :: newVec,
!!                  integer(INOUT)          :: numVec,
!!                  OPTIONAL,integer(OUT)   :: vecLastFree,
!!                  OPTIONAL,integer(IN)    :: copyDirection,
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
!!   numVars    :  number of variables per cell that are stored in newVec.
!!                 For use with the Eos unit as in this example, this should be = EOS_NUM.
!!   newVec     : storage for newly generated vectors.
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
!!                 in newVec, then vecLastFree (if present) will be set to -1 to indicate
!!                 this condition.
!!   copyDirection  : Direction of data transfer.
!!                    Use GRID_COPYDIR_TO_VECT to tranfer data from Grid blocks to vectors (default);
!!                    use GRID_COPYDIR_FROM_VECT to tranfer data from vectors back to Grid blocks.
!!   gridDataStruct : whether cell centered, face centered etc. (may be deprecated later).
!!                    Should be CENTER if present. FLASH does not require
!!                    support for marshalling other datastructs (FACEX, etc.)
!!                    for Eos calls.
!!
!! NOTES
!!
!!  This is for demonstration, it may or may not work.
!!
!!  GRID_COPYDIR_TO_VECT and GRID_COPYDIR_FROM_VECT are defined in Grid_interface.F90 .
!!
!! SEE ALSO
!!
!!  Eos_getData
!!  Eos_putData
!!  Eos
!!  Grid_interface
!!***

!!REORDER(5): unk
!!REORDER(4): dataPtr

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_makeVector(vecLen,numVars,newVec,numVec,vecLastFree,copyDirection,gridDataStruct)

#include "constants.h"
#include "Flash.h"


  implicit none

  integer, intent(in) :: vecLen
  integer, intent(in) :: numVars
  integer,intent(INOUT) :: numVec
  real, dimension(vecLen,numVars,numVec),intent(INOUT) :: newVec
  integer, optional,intent(OUT):: vecLastFree
  integer, optional,intent(in) :: copyDirection
  integer, optional,intent(in) :: gridDataStruct
end subroutine Grid_makeVector








