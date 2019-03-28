!!****if* source/Grid/Grid_addCoarseToFluxRegister
!!
!! NAME
!!  Grid_addCoarseToFluxRegister
!!
!! SYNOPSIS
!!  call Grid_addCoarseToFluxRegister(integer(IN) :: coarse_level,
!!                          optional, logical(IN) :: isDensity(:),
!!                          optional, real(IN)    :: coefficient)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  In normal
!!  use, client code could add flux data from both levels into the flux register
!!  for use with adjusting flux data on the coarse level.
!!
!!  This routine allows client code to request that the Grid unit add coarse data
!!  from the Grid unit's flux data structures to the contents of the associated
!!  flux registers.  This routine is clearly intended for use with AMR.  Note
!!  that the flux registers may choose to only store flux data that exists at 
!!  fine/coarse boundaries.
!!
!!  All data stored in the Grid unit's flux data structures as flux densities
!!  will automatically be transformed to flux before applying to the flux
!!  register.
!!
!!  Additionally, a multiple scale factor may be applied to all flux data before
!!  passing the data to the flux register.
!!
!!  It is assumed that before calling this routine, the client code has already
!!  written flux data to Grid's data structures using the Grid_getFluxPtr
!!  interface.
!!
!! ARGUMENTS
!!  coarse_level - the 1-based level index (1 is the coarsest level) indicating
!!                 which level's data should be added to the flux register as
!!                 coarse data.
!!  isDensity - a mask that identifies which physical flux quantities are
!!              actually stored in the Grid unit's flux data structures as
!!              flux densities.  If no mask is given, it is assumed that data
!!              is stored as flux.
!!  coefficient - a scaling parameter to apply to all flux data before applying
!!                the data to the flux register.
!!
!! SEE ALSO
!!   Grid_getFluxPtr/Grid_releaseFluxPtr
!!   Grid_zeroFluxRegister
!!   Grid_addCoarseToFluxRegister
!!   Grid_overwriteFluxes
!!
!!***

subroutine Grid_addCoarseToFluxRegister(coarse_level, isDensity, coefficient)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN)           :: coarse_level
  logical, intent(IN), optional :: isDensity(:)
  real,    intent(IN), optional :: coefficient

  call Driver_abortFlash("[Grid_addCoarseToFluxRegister] Prototype stub.  Do NOT use!")
end subroutine Grid_addCoarseToFluxRegister

