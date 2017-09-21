!!****if* source/Grid/GridMain/Chombo/AMR/gr_primitiveToConserve
!!
!! NAME
!!
!!  gr_primitiveToConserve
!!
!!
!! SYNOPSIS
!!
!!  gr_primitiveToConserve(integer(in) :: blkList(count),
!!                         integer(in) :: count)
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity) to the corresponding conservative form
!!  (i.e., momentum) if gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!! NOTES
!!
!!  The variables that are converted are the named cell-centered
!!  solution variables marked to be of type PER_MASS explicitly in a
!!  Config file.  Additionally, abundances and mass scalars are
!!  considered to be of type PER_MASS.
!!
!!  For proper functioning, DENS_VAR must not be marked as PER_MASS!
!!
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_conserveToPrimitive
!!
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are now conserved.  No
!!  such mechanism exists in the code yet.
!!
!!  This routine accesses the global variable storage 
!!  array solnData directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_primitiveToConserve(blkList,count)
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none

  integer,intent(IN) :: count
  integer,dimension(count),intent(IN) :: blkList

  call Driver_abortFlash("[gr_primitiveToConservative]: Not implemented yet!")
end subroutine gr_primitiveToConserve
