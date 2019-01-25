!!****f* source/Grid/Grid_renormAbundance
!!
!! NAME
!!
!!  Grid_renormAbundance
!!
!!
!! SYNOPSIS
!!
!!  Grid_renormAbundance(integer(IN) :: blockID,
!!                     integer(IN) :: blkLimits(2,MDIM),
!!                     real,pointer :: solnData(:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Renormalize the abundances in a given block so they sum to 1.
!!  This should be used before calling the EOS.  Each abundance is
!!  restricted to fall between smallx and 1.
!!
!!  Also check the abundance conservation and report if it is worse
!!  than abundErr.
!!
!!  This routine is called automatically by Hydro and MHD in FLASH 
!!  if the irenorm runtime parameter is set to 1.
!!
!!  Only abundances/fluids which contribute to the EOS are included.
!!
!!
!! ARGUMENTS
!!
!!  blockDesc -   the block number to renormalize
!!  blkLimits - the index limits for internal zones of the block to renormalize
!!  solnData -  Pointer to the block to be renormalized
!!
!!
!! PARAMETERS
!!
!!  smallx -    the cutoff value for the composition mass fraction
!!
!!
!! SEE ALSO 
!!
!!  Grid_limitAbundance
!!
!!
!!***

#include "constants.h"

subroutine Grid_renormAbundance(tileDesc, solnData)
  use flash_tile,   ONLY : flash_tile_t

  implicit none

  type(flash_tile_t), intent(IN)         :: tileDesc
  real,                          pointer :: solnData(:,:,:,:)

  return
end subroutine Grid_renormAbundance

