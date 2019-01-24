!!****f* source/Simulation/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!!
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution for a given block.
!!
!!  This is simulation-dependent, there is no general implementation.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!! SIDE EFFECTS
!!
!!  The analytical solution is computed and stored in an appropriate slot
!!  (or slots) in the solution vector, UNK.
!!
!! NOTES
!! 
!!  By default this is just a stub that does not do anything.
!!
!!***

subroutine Simulation_computeAnalytical(solnData, tileDesc, tcurr)

  use flash_tile, ONLY : flash_tile_t
  implicit none

#include "FortranLangFeatures.fh"
  
  real,dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  type(flash_tile_t), intent(in) :: tileDesc
  real,    intent(in) :: tcurr

end subroutine Simulation_computeAnalytical
