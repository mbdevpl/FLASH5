!!****if* source/Grid/GridBoundaryConditions/gr_bcApplyToAllBlks
!!
!! NAME
!!  gr_bcApplyToAllBlks
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToAllBlks(integer(IN) :: axis,
!!                      logical(IN) :: isWork)
!!  
!! DESCRIPTION 
!!
!!  This routine is a wrapper around gr_bcApplyToOneFace, and is used by UG and PM2.
!!  It calls gr_bcApplyToOneFace one each for lowerface and upperface, and repeats
!!  the process for all blocks in the grid.
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    isWork         - is always false for UG. In PM2, if true, indicates that
!!                     the boundary conditions should be applied to work array
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

#include "constants.h"
#include "Flash.h"
  
subroutine gr_bcApplyToAllBlks(axis,isWork)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
  integer, intent(in) :: axis
  logical, intent(in) :: isWork

  ! DEV: TODO Implement this if needed
  call Driver_abortFlash("[gr_bcApplyToAllBlks] not implemented yet")

end subroutine gr_bcApplyToAllBlks

