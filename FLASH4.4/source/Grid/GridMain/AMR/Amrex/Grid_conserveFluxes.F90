!!****f* source/Grid/GridMain/AMR/Amrex/Grid_conserveFluxes
!!
!! NAME
!!  Grid_conserveFluxes
!!
!! SYNOPSIS
!!  call Grid_conserveFluxes(integer(IN) :: axis,
!!                           integer(IN) :: level)
!!  
!! DESCRIPTION 
!!  When FLASH is run with AMR, it is possible that some leaf blocks
!!  will have a neighboring leaf block that is refined at the next coarsest
!!  level.  To maintain conservation, the flux entering into the
!!  fine block at this shared boundary must equal the flux leaving the
!!  coarse block at the boundary.
!!
!!  To enforce this requirement, this routine overwrites in the coarse
!!  block the flux at the shared boundary with the averaged flux from the
!!  fine block, which is sensible as the fine flux data should be at least as
!!  accurate as the flux in the coarse level.
!! 
!!  It is assumed that before calling this routine, the code has already
!!  loaded the corrected fluxes for the fine level into the flux 
!!  registers using Grid_putFluxData and that the uncorrected flux for the
!!  coarse level has been stored in Grid with the Grid_getFluxPtr interface.
!!  This routine will only overwrite the flux data at fine/coarse boundaries.
!!
!! ARGUMENTS 
!!  axis - the only acceptable value for AMReX is ALLDIR.
!!  level - the 1-based level index of the coarse blocks
!!
!! SEE ALSO
!!  Grid_getFluxPtr/Grid_releaseFluxPtr
!!  Grid_putFluxData
!!
!!***

#include "constants.h"

subroutine Grid_conserveFluxes(axis, level)
    use amrex_fort_module,    ONLY : wp => amrex_real
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    use Driver_interface,     ONLY : Driver_abortFlash
    use Grid_interface,       ONLY : Grid_getGeometry
    use gr_physicalMultifabs, ONLY : fluxes, &
                                     flux_registers

    implicit none

    integer, intent(IN) :: axis
    integer, intent(IN) :: level

    integer :: geometry

    if (axis /= ALLDIR) then
        call Driver_abortFlash("[Grid_conserveFluxes] AMReX requires axis==ALLDIR")
    end if

    ! No need to conserve on the finest level in existence or any
    ! level index corresponding to a finer mesh
    !
    ! AMReX level index is 0-based
    if (level-1 >= amrex_get_finest_level())     RETURN

    call Grid_getGeometry(geometry)

    select case (geometry)
    case (CARTESIAN)
        ! The AMReX flux registers are dealing with fluxes and 
        ! *not* flux densities.  In Grid_putFluxData, flux densities at the fine
        ! level were scaled to fluxes with the assumption that the cell lengths
        ! at the coarse level are one.  Therefore, reconversion to flux densities
        ! is automatic here.
        call flux_registers(level)%overwrite(fluxes(level-1, :), 1.0_wp)
    case default
        ! DEV: TODO This routine should take an isFluxDensity array as an argument
        ! so that the routine can determine which need to be scaled and which do not
        call Driver_abortFlash("[Grid_conserveFluxes] Only works with Cartesian")
    end select
end subroutine Grid_conserveFluxes

