!!****f* source/Grid/Grid_overwriteFluxes
!!
!! NAME
!!  Grid_overwriteFluxes
!!
!! SYNOPSIS
!!  call Grid_overwriteFluxes(integer(IN) :: axis,
!!                            integer(IN) :: level,
!!                            logical(IN) :: isDensity(:))
!!
!! DESCRIPTION 
!!  Each flux register managed by the Grid unit is associated with a fine and a
!!  coarse level.  In normal use, client code will add flux data from either
!!  level to the flux register.  At a later point, the flux register data that 
!!  exists at fine/coarse boundaries will be applied in an application-specific 
!!  way to the flux data also managed the Grid unit.
!!
!!  This routine allows for overwriting coarse flux data at fine-coarse
!!  boundaries with the data in the flux register.  In particular, the applied
!!  data is the sum of the coarse flux with the fine flux averaged down to the
!!  coarse level.
!!
!! ARGUMENTS 
!!  axis - apply the flux register data only at fine/coarse boundaries 
!!         perpendicular to the direction indicated by this value.  Valid
!!         values are IAXIS, JAXIS, KAXIS, or ALLDIR.  For the latter value,
!!         the operation will be applied at all fine/coarse boundaries.
!!  level - the 1-based level index (1 is the coarsest level) that indicates
!!          the refinement level of the fluxes to be altered.
!!  isDensity - a mask that identifies which physical flux quantities are
!!              actually stored in the Grid unit's flux data structures as
!!              flux densities.  If no mask is given, it is assumed that data
!!              is stored as flux.
!!
!! SEE ALSO
!!  Grid_getFluxPtr/Grid_releaseFluxPtr
!!  Grid_zeroFluxRegister
!!  Grid_addFineToFluxRegister
!!  Grid_addCoarseToFluxRegister
!!
!!***

subroutine Grid_overwriteFluxes(axis, level, isDensity)
    implicit none

    use Driver_interface, ONLY : Driver_abortFlash

    integer, intent(IN)           :: axis
    integer, intent(IN)           :: level
    logical, intent(IN), optional :: isDensity(:)

    call Driver_abortFlash("[Grid_overwriteFluxes] Prototype stub.  Do NOT use!")
end subroutine Grid_overwriteFluxes

