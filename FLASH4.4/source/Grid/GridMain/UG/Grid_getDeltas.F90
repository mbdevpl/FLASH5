!!****if* source/Grid/GridMain/UG/Grid_getDeltas
!!
!! NAME
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  call Grid_getDeltas(integer(IN):: level, 
!!                      real(OUT)::del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the dx/dy/dz for a given blockId on the Grid
!!  dx is the size of one zone in the x direction of a block
!!
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!  
!!***
#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getDeltas(level,del)
  use Grid_data, ONLY : gr_delta
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"

  integer, intent(IN) :: level
  real, dimension(MDIM),intent(OUT) :: del

  if(level /= 1) call Driver_abortFlash("any value of level other than 1 is not valied")
  del = gr_delta(:,1)
  return
end subroutine Grid_getDeltas

