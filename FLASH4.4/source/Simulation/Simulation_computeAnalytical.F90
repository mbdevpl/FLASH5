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

subroutine Simulation_computeAnalytical(solnData, blockDesc, tcurr)

  use block_metadata,   ONLY : block_metadata_t
  implicit none

#include "FortranLangFeatures.fh"
  
  real,dimension(:,:,:,:),POINTER_INTENT_IN :: solnData
  type(block_metadata_t), intent(in) :: blockDesc
  real,    intent(in) :: tcurr

end subroutine Simulation_computeAnalytical
