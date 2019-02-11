!!****if* source/Grid/GridMain/Grid_addToVar
!!
!! NAME
!!
!!  Grid_addToVar
!!
!! SYNOPSIS
!!
!!  call Grid_addToVar(integer(in) :: srcVar,
!!                     integer(in) :: destVar,
!!                     real(in) :: multFactor,
!!                     logical(in) :: reset)
!!
!! DESCRIPTION
!!   Compute solnData(srcVar,:,:,:)*multFactor and save in solnData(destVar,:,:,:).
!!
!!   If reset is true, the destination variable is first zeroed;
!!   otherwise the product is added to the existing values of destVar.
!!
!!   The operation is applied to interior cells of all LEAF blocks.
!!
!! ARGUMENTS
!!
!!
!!   srcVar : the state variables to be used in the RHS of the expression
!!
!!   destVar : the state variables to be used in the LHS of the expression
!!
!!   multFactor : multiplication factor
!!
!!   reset : indicates whether the destination variable should be zeroed first
!!
!! NOTES
!!
!!   srcVar == destVar is allowed and behaves as expected iff reset is .FALSE.
!!
!!   For a copy call Grid_addToVar(srcVar, destVar, 1.0, .true.)
!!
!!***

subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)
  use Grid_interface, ONLY : Grid_getBlkPtr,  Grid_releaseBlkPtr,&
       Grid_getLeafIterator
  use block_metadata, ONLY : block_metadata_t
  use leaf_iterator, ONLY : leaf_iterator_t
#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  implicit none
  integer, intent(in) :: srcVar, destVar
  real,    intent(in) :: multFactor
  logical, intent(in) :: reset
  integer :: blockList(MAXBLOCKS), blockCount
  integer :: blkLimits(LOW:HIGH,MDIM), blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blk, k, j, i, n
  real, dimension(:,:,:,:), pointer :: solnData
  type(block_metadata_t) :: block
  type(leaf_iterator_t) :: itor
  
  call Grid_getLeafIterator(itor)
  do while(itor%is_valid())
     call itor%blkMetaData(block)
     
     call Grid_getBlkPtr(block, solnData)
     if (reset) then 
        solnData(destVar,:,:,:) = 0.0
     end if
     solnData(destVar,:,:,:) = solnData(destVar,:,:,:) + &
          multFactor*solnData(srcVar,:,:,:)
     call Grid_releaseBlkPtr(blockList(blk), solnData)
  end do
end subroutine Grid_addToVar
