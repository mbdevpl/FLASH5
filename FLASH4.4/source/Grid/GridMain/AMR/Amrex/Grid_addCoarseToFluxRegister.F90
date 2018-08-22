!!****if* source/Grid/GridMain/AMR/Amrex/Grid_addCoarseToFluxRegister
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

#include "Flash.h"
#include "constants.h"

subroutine Grid_addCoarseToFluxRegister(coarse_level, isDensity, coefficient)
    use amrex_fort_module,    ONLY : wp => amrex_real
    use amrex_amrcore_module, ONLY : amrex_get_finest_level, &
                                     amrex_ref_ratio

    use Driver_interface,     ONLY : Driver_abortFlash
    use Grid_interface,       ONLY : Grid_getGeometry
    use gr_physicalMultifabs, ONLY : flux_registers, &
                                     fluxes

    implicit none

    integer, intent(IN)           :: coarse_level
    logical, intent(IN), optional :: isDensity(:)
    real,    intent(IN), optional :: coefficient

    integer  :: coarse
    integer  :: fine
    integer  :: geometry
    real(wp) :: coef

    if (NFLUXES < 1) then
        RETURN
    end if
 
    ! FLASH uses 1-based level index / AMReX uses 0-based index
    coarse = coarse_level - 1
    fine   = coarse_level

    ! The finest refinement level is never the coarse level of a flux register
    if ((coarse < 0) .OR. (coarse >= amrex_get_finest_level())) then
        call Driver_abortFlash("[Grid_addCoarseToFluxRegister] Invalid level")
    end if

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
      if (amrex_ref_ratio(coarse) /= 2) then
        call Driver_abortFlash("[Grid_addFineToFluxRegister] refinement ratio not 2")
      end if

#if   NDIM == 2
        coef = coef * 0.5_wp
#elif NDIM == 3
        coef = coef * 0.25_wp
#endif

        ! Flux registers index is 0-based index of fine level
        call flux_registers(fine)%crseadd(fluxes(coarse, 1:NDIM), coef)
    case default
        call Driver_abortFlash("[Grid_addFineToFluxRegister] Only works with Cartesian")
    end select
end subroutine Grid_addCoarseToFluxRegister

