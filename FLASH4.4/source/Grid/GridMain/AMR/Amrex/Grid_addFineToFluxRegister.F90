!!****if* source/Grid/Grid_addFineToFluxRegister
!!
!! NAME
!!  Grid_addFineToFluxRegister
!!
!! SYNOPSIS
!!  call Grid_addFineToFluxRegister(integer(IN) :: level,
!!                        optional, logical(IN) :: isDensity(:),
!!                        optional, real(IN)    :: coefficient)
!!
!! DESCRIPTION 
!!  Each flux register is associated with a fine and a coarse level.  In normal
!!  use, client code could add flux data from both levels into the flux register
!!  for use with adjusting flux data on the coarse level.
!!
!!  This routine allows client code to request that the Grid unit add fine data
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
!!  level - the 1-based level index (1 is the coarsest level) indicating which
!!          level's data should be added to the flux register as fine data.
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

#include "Flash.h"
#include "constants.h"

subroutine Grid_addFineToFluxRegister(level, isDensity, coefficient)
    use amrex_fort_module,    ONLY : wp => amrex_real
#if DEBUG_GRID
    use amrex_amrcore_module, ONLY : amrex_ref_ratio
#endif

    use Driver_interface,     ONLY : Driver_abortFlash
    use Grid_interface,       ONLY : Grid_getGeometry
    use gr_physicalMultifabs, ONLY : flux_registers, &
                                     fluxes

    implicit none
    
    integer, intent(IN)           :: level
    logical, intent(IN), optional :: isDensity(:)
    real,    intent(IN), optional :: coefficient

    integer  :: geometry
    real(wp) :: coef

    if (present(coefficient)) then
        coef = coefficient
    else
        coef = 1.0_wp
    end if

    if (present(isDensity)) then
        call Driver_abortFlash("[Grid_addFineToFluxRegister] isDensity not implemented")
    end if

    call Grid_getGeometry(geometry)

    select case (geometry)
    case (CARTESIAN)
      ! The scaling factor=1/r^(NDIM-1) used here assumes that the refinement
      ! ratio, r, between levels is always 2
#if DEBUG_GRID
        if (amrex_ref_ratio(level-1) /= 2) then
          call Driver_abortFlash("[Grid_addFineToFluxRegister] refinement ratio not 2")
        end if
#endif

#if   NDIM == 2
        coef = coef * 0.5_wp
#elif NDIM == 3
        coef = coef * 0.25_wp
#endif

        call flux_registers(level-1)%fineadd(fluxes(level-1, 1:NDIM), coef)
    case default
        call Driver_abortFlash("[Grid_addFineToFluxRegister] Only works with Cartesian")
    end select
end subroutine Grid_addFineToFluxRegister

