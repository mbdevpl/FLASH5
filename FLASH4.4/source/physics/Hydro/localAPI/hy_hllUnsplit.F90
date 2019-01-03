!!****if* source/physics/Hydro/localAPI/hy_hllUnsplit
!!
!! NAME
!!  hy_hllUnsplit
!!
!! SYNOPSIS
!!  call hy_hllUnsplit( integer (IN) :: blockCount,
!!                      integer (IN) :: blockList(blockCount),
!!                      real    (IN) :: dt,
!!                      real    (IN) :: dtOld  )
!!
!! DESCRIPTION
!!  Performs Hydro update in a directionally unsplit fashion over a set
!!  of blocks.
!!  dt gives the timestep over which to advance, and timeEndAdv gives the
!!  simulation time at the end of the update.
!! 
!!  This routine performs a guardcell fill and for each block: 
!!   - applies an Eos call to the guard cells (!DEV: only if/where necessary?)
!!   - computes fluxes
!!   - update all the cell values from the fluxes.
!!   - and finally, we apply an Eos call to the block (interiors).
!!
!! ARGUMENTS
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!  dtOld      -  old timestep (not used here)
!!
!! NOTES
!!  This is a simple demo version. Some of the numerous limitiations, compared
!!  to the more serious Hydro implementations available in FLASH:
!!
!!  * Flux correction is not implemented.
!!  * No reconstruction (and thus no limiting) of variables; only HLL as "Riemann solver".
!!  * No support for advecting mass fractions / abundances or other mass scalar
!!    variables.
!!  * No support for non-Cartesian geometries. If this works for any of them,
!!    it should be considered an accident.
!!  * No support for MHD, or anything related to magnetic fields.
!!  * No support for gravity or other body forces or source terms.
!!  * No support for flux-based diffusive terms.
!!  * No artificial viscosity term.
!!  * No support for eintSwitch .NE. 0.0.
!!
!! HISTORY
!!  June  2013  - created KW, outer structure derived from hy_uhd_unsplit.F90 (Dongwook)
!!
!!***

! Note: the following arrays need to be spelled exactly like this in the code below,
!       preserving case.
!!REORDER(4): Uin, Uout, face[XYZ], auxC

#include "constants.h"

Subroutine hy_hllUnsplit ( tileLimits, Uin, plo, Uout, del, dt )
  implicit none

  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: UOUT(plo(1):,plo(2):,plo(3):,plo(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt

  Uout(:, :, :, :) = 0.0
End Subroutine hy_hllUnsplit

