!!****if* source/Grid/GridMain/UG/Grid_getBlkPhysicalSize
!!
!! NAME
!!  Grid_getBlkPhysicalSize
!!
!! SYNOPSIS
!!
!!  Grid_getBlkPhysicalSize(type(block_metadta_t)(IN) :: block,
!!                          real(OUT) :: blockSize(MDIM))
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
!!  This interface is general enough to work even
!!  when the block sizes (in terms of number of cells) are not uniform.
!!  
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkPhysicalSize(block, blockSize)
  use Driver_interface, ONLY : Driver_abortFlash
  use block_metadata,   ONLY : block_metadata_t

  implicit none

  type(block_metadata_t), intent(in) :: block
  real,dimension(MDIM),intent(out) :: blockSize

  call Driver_abortFlash("[Grid_getBlkPhysicalSize] DEPRECATED: Use tile object")
end subroutine Grid_getBlkPhysicalSize

