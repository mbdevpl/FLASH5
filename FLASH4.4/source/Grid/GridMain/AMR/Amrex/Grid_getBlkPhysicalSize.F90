!!****f* source/Grid/GridMain/AMR/Amrex/Grid_getBlkPhysicalSize
!!
!! NAME
!!  Grid_getBlkPhysicalSize
!!
!! SYNOPSIS
!!  Grid_getBlkPhysicalSize(integer(IN)  :: blockId,
!!                          real(OUT)    :: blockSize(MDIM))
!!  
!! DESCRIPTION 
!!  Gets the physical size of the block along each dimension by calculating
!!  the deltas of each direction times the number of zones.  It is important
!!  to calculate the size of a block consistently because rounding errors
!!  can occur in some cases.  Use this function to get the size of the 
!!  block rather than calculating it in a new way.  
!!
!! ARGUMENTS
!!  blockId -local block number
!!  blockSize - returned array of size MDIM holding the size of 
!!              each dimension of the block
!!
!! EXAMPLE
!!
!!    If there is 2-d block with 8 cells along IAXIS and 4 cells along
!!    JAXIS, and a single cell size (or delta) along each dimension is 0.1, 
!!    then the part of the physical domain occupied by the block is 
!!    0.8 x 0.4 in size, as in the block shown below
!!    
!!     - - - - - - - - (1.0,0.8)
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!    | | | | | | | | |
!!     - - - - - - - - 
!!   (0.2,0.4)
!!  
!!
!! NOTES
!!
!!  The physical block sizes returned by this routine are
!!  Calculated using deltas.
!!
!!  This interface uses constant dx times number cells in the blocks
!!  Will work when the all cells in the block are uniform size
!!  An amrex internal routine returning physical size should replace the 
!!  body of this subroutine, if such routine exists in amrex.
!!  
!!
!!***

subroutine Grid_getBlkPhysicalSize(block, blockSize)

#include "Flash.h"
#include "constants.h"
  use amrex_amrcore_module,  ONLY : amrex_geom
  use amrex_geometry_module, ONLY : amrex_problo
  use block_metadata,        ONLY : block_metadata_t

  implicit none
  type(block_metadata_t),intent(in) :: block
  real,dimension(MDIM),intent(out) :: blockSize

  associate(dx   => amrex_geom(block%level - 1)%dx, &
            lo   => block%limits(LOW,  :), &
            hi   => block%limits(HIGH, :))
    ! lo is 1-based cell-index of lower-left cell in block 
    ! hi is 1-based cell-index of upper-right cell in block
  blockSize(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) * dx(1:NDIM)
  end associate
end subroutine Grid_getBlkPhysicalSize

subroutine Grid_getBlkPhysicalSize_blkId(blockId, blockSize)

#include "Flash.h"
#include "constants.h"
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none
   integer, intent(in) :: blockId
   real,dimension(MDIM),intent(out) :: blockSize

   call Driver_abortFlash('Grid_getBlkPhysicalSize_blkId: just a stub for this Grid configuration, do not call!')

end subroutine Grid_getBlkPhysicalSize_blkId
