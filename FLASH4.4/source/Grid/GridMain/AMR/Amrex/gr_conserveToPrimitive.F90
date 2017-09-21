!!****if* source/Grid/GridMain/Chombo/AMR/gr_conserveToPrimitive
!!
!! NAME
!!
!!  gr_conserveToPrimitive
!!
!!
!! SYNOPSIS
!!
!!  gr_conserveToPrimitive(integer(in) :: blkList(count),
!!                         integer(in) :: count,
!!                         logical(in) :: allCells)
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  convert variables which are normally represented in PER_MASS
!!  form (e.g., velocity) from the corresponding conservative form
!!  (i.e., momentum) back to the normal PER_MASS form  if 
!!  gr_convertToConsvdForMeshCalls is TRUE.
!!  Do nothing if gr_convertToConsvdForMeshCalls is FALSE.
!!
!!  Additionally,
!!   - energies (ENER_VAR and EINT_VAR) are forced to be .ge. gr_smalle,
!!   - the density (DENS_VAR) is forced to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   allCells - act on all cells, including guardcells, if .TRUE.,
!!              otherwise only modify interior cells.
!!
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
!!
!! SEE ALSO
!!
!!  Simulation_getVarnameType
!!  gr_primitiveToConserve 
!!
!!
!! BUGS
!!
!!  This routine does not set the variable attributes to 
!!  indicate that the variables are no longer conserved.  No
!!  such mechanism exists in the code yet.
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_conserveToPrimitive(blkList,count,allCells)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  logical, intent(IN):: allCells
  
  call Driver_abortFlash("[gr_conservativeToPrimitive]: Not implemented yet!")
end subroutine gr_conserveToPrimitive
