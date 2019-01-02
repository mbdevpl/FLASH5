!!****if* source/Grid/GridMain/Chombo/Grid_getBlkBC_desc
!!
!! NAME
!!  Grid_getBlkBC_desc
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkBC(block_metadata_t(IN)  :: blockDesc,
!!                integer(OUT) :: faces(2,MDIM),
!!       optional,integer(OUT) :: onBoundary(2,MDIM))
!!                    
!! DESCRIPTION 
!!  Returns the boundary condition for each face of the block.
!!  
!!  this function finds out if a block face is
!!  on the physical boundary, and if so, returns the boundary condition.
!!  The boundary conditions are defined in the header file
!!  constants.h, ie. OUTFLOW, REFLECTING, PERIODIC.  NOT_BOUNDARY, also
!!  defined in constants.h is returned if a block face is not on a 
!!  physical boundary.
!!   
!! ARGUMENTS 
!!
!!  blockDesc - descriptor for local block
!!  faces   - array returned holding boundary conditions
!!                
!!            the first index of the array can take on values LOW or
!!            HIGH, and the second index can be IAXIS, JAXIS or KAXIS
!! onBoundary - array returned with boundary conditions including periodic
!!
!! NOTES
!!
!!   The #define constants LOW, HIGH, IAXIS, JAXIS and KAXIS
!!   are defined in constants.h and are
!!   meant to ease the readability of the code.       
!!   instead of faces(2,3) = PERIODIC the code reads
!!   faces(HIGH, KAXIS) = PERIODIC
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine Grid_getBlkBC(blockDesc, faces, onBoundary)
  use Grid_data, ONLY : gr_globalDomain, gr_domainBC
  use Grid_interface, ONLY : Grid_getBlkBoundBox
  use block_metadata, ONLY : block_metadata_t
  implicit none
  type(block_metadata_t), intent(in) :: blockDesc
  integer, dimension(2,MDIM),intent(out):: faces
  integer, optional, dimension(LOW:HIGH,MDIM),intent(out):: onBoundary
end subroutine Grid_getBlkBC
